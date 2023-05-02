from __future__ import annotations

import re
from abc import ABC, abstractmethod

import ibis


# TODO: invert
class AbstractFilterExpr(ABC):
    """Base class for all filter expressions. Defines logical operators and interface."""

    def __and__(self, other):
        return AndFilterExpr(self, other)

    def __or__(self, other):
        return OrFilterExpr(self, other)

    def __invert__(self):
        return NotFilterExpr(self)

    @abstractmethod
    def convert(self) -> ibis.expr.deferred.Deferred:
        """Convert genomic-features filter expression to ibis deferred expression."""
        pass

    @abstractmethod
    def required_tables(self) -> set[str]:
        """Tables required by this filter."""
        pass

    @abstractmethod
    def columns(self) -> set[str]:
        """Columns required by this filter."""
        pass


class EmptyFilter(AbstractFilterExpr):
    """An empty filter that returns all rows."""

    def __repr__(self) -> str:
        return "EmptyFilter()"

    def convert(self) -> None:
        return None

    def required_tables(self) -> set[str]:
        return set()

    def columns(self) -> set[str]:
        return set()


class AbstractFilterOperatorExpr(AbstractFilterExpr):
    def __init__(self, left: AbstractFilterExpr, right: AbstractFilterExpr):
        self.left = left
        self.right = right

    def required_tables(self) -> set[str]:
        return self.left.required_tables() & self.right.required_tables()


    def columns(self) -> set[str]:
        return self.left.columns() | self.right.columns()


class AndFilterExpr(AbstractFilterOperatorExpr):
    """Logical and of two filters."""

    def __repr__(self) -> str:
        return f"({self.left} & {self.right})"

    def convert(self) -> ibis.Expr:
        return self.left.convert() & self.right.convert()


class NotFilterExpr(AbstractFilterExpr):
    """A filter that inverts the result of another filter."""

    def __init__(self, expr: AbstractFilterExpr):
        self.expr = expr

    def __repr__(self) -> str:
        return f"~{repr(self.expr)}"

    def convert(self) -> ibis.expr.deferred.Deferred:
        return ~self.expr.convert()

    def columns(self) -> set[str]:
        return self.expr.columns()

    def required_tables(self) -> set[str]:
        return self.expr.required_tables()


class OrFilterExpr(AbstractFilterOperatorExpr):
    """Logical or of two filters."""

    def __repr__(self) -> str:
        return f"({self.left} | {self.right})"

    def convert(self) -> ibis.Expr:
        return self.left.convert() | self.right.convert()


class AbstractFilterEqualityExpr(AbstractFilterExpr):
    def __init__(self, value: str | list[str]):
        self.value = value

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({self.value})"

    def convert(self) -> ibis.expr.deferred.Deferred:
        if isinstance(self.value, str):
            return ibis.deferred[list(self.columns())[0]] == self.value
        else:
            return ibis.deferred[list(self.columns())[0]].isin(self.value)


#  ------------------------ OVERLAP TYPE: any ------------------------ #
# Range:       |==============================|
# Annotation: |---|     |-------|           |-----| |---|
# Selected:    ***       *******             *****

#  ------------------------ OVERLAP TYPE: within --------------------- #
# Range:       |==============================|
# Annotation: |---|     |-------|           |-----| |---|
# Selected:              *******


class AbstractFilterRangeExpr(AbstractFilterExpr):
    def __init__(self, value: str, type: str = "any"):
        self.value = value
        self.type = type

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({self.value})"

    def convert(self) -> ibis.expr.deferred.Deferred:
        match = re.match(r"^(\w+):(\d+)-(\d+)$", self.value)
        if match is None:
            raise ValueError(
                "Invalid range format. Valid format is '{seq_name}:{start}-{end}'"
            )
        range_seq_name, range_start, range_end = match.groups()
        start_column, end_column, seq_name_column = self._range_columns
        if self.type == "any":
            return (ibis.deferred[seq_name_column] == range_seq_name) & (
                (
                    (ibis.deferred[end_column] >= int(range_start))
                    & (ibis.deferred[end_column] <= int(range_end))
                )
                | (
                    (ibis.deferred[start_column] >= int(range_start))
                    & (ibis.deferred[start_column] <= int(range_end))
                )
            )
        elif self.type == "within":
            return (ibis.deferred[seq_name_column] == range_seq_name) & (
                (ibis.deferred[end_column] <= int(range_end))
                & (ibis.deferred[start_column] >= int(range_start))
            )
        else:
            raise ValueError(
                "Invalid overlap type. Valid options are 'any' and 'within'"
            )


class GeneIDFilter(AbstractFilterEqualityExpr):
    """Filter by gene_id column."""

    def columns(self) -> set[str]:
        return {"gene_id"}

    def required_tables(self) -> set[str]:
        # TODO: Joining on gene_id is not necessary for transcript queries
        return {"gene"}


class GeneBioTypeFilter(AbstractFilterEqualityExpr):
    """Filter by gene_biotype column."""

    def columns(self) -> set[str]:
        return {"gene_biotype"}

    def required_tables(self) -> set[str]:
        return {"gene"}


class GeneRangesFilter(AbstractFilterRangeExpr):
    """
    Filter features within a genomic range

    Parameters
    ----------
    value : str
        Genomic range in the format "seq_name:start-end"
    type : str
        String indicating how overlaps are to be filters.
        Options are 'any' and 'within'. Default is 'any'
    """

    @property
    def _range_columns(self) -> list[str]:
        return ["gene_seq_start", "gene_seq_end", "seq_name"]

    def columns(self) -> set[str]:
        return set(self._range_columns)

    def required_tables(self) -> set[str]:
        return {"gene"}


# class GeneIDFilter(AbstractFilterExpr):
#     def __init__(self, gene_id: str | list[str]):
#         self.gene_id = gene_id

#     def convert(self) -> ibis.expr.deferred.Deferred:
#         if isinstance(self.gene_id, str):
#             return ibis.deferred.gene_id == self.gene_id
#         else:
#             return ibis.deferred.gene_id.isin(self.gene_id)

# def required_tables(self) -> set[str]:
#     # TODO: Joining on gene_id is not necessary for transcript queries
#     return {"gene"}


# class GeneBioTypeFilter(AbstractFilterExpr):
#     def __init__(self, gene_biotype: str | list[str]):
#         self.gene_biotype = gene_biotype

#     def convert(self) -> ibis.expr.deferred.Deferred:
#         if isinstance(self.gene_biotype, str):
#             return ibis.deferred.gene_biotype == self.gene_biotype
#         else:
#             return ibis.deferred.gene_biotype.isin(self.gene_biotype)
