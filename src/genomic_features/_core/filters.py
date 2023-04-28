from __future__ import annotations

from abc import ABC, abstractmethod

import ibis


# TODO: invert
class AbstractFilterExpr(ABC):
    def __and__(self, other):
        return AndFilterExpr(self, other)

    def __or__(self, other):
        return OrFilterExpr(self, other)

    @abstractmethod
    def convert(self) -> ibis.expr.deferred.Deferred:
        pass

    @abstractmethod
    def required_tables(self) -> set[str]:
        pass

    @abstractmethod
    def columns(self) -> set[str]:
        pass


class EmptyFilter(AbstractFilterExpr):
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
    def __repr__(self) -> str:
        return f"({self.left} & {self.right})"

    def convert(self) -> ibis.Expr:
        return self.left.convert() & self.right.convert()


class OrFilterExpr(AbstractFilterOperatorExpr):
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
        try:
            range_seq_name, range_coords = self.value.split(":")
            range_start, range_end = range_coords.split("-")
        except ValueError:
            raise ValueError(
                "Invalid range format. Valid format is 'seq_name:start-end'"
            )
        start_column, end_column, seq_name_column = list(self.columns())
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
    def columns(self) -> set[str]:
        return {"gene_id"}

    def required_tables(self) -> set[str]:
        # TODO: Joining on gene_id is not necessary for transcript queries
        return {"gene"}


class GeneBioTypeFilter(AbstractFilterEqualityExpr):
    def columns(self) -> set[str]:
        return {"gene_biotype"}

    def required_tables(self) -> set[str]:
        return {"gene"}


class RangesFilter(AbstractFilterRangeExpr):
    """
    Filter features within a genomic range

    Parameters
    ----------
    value : str
        Genomic range in the format "seq_name:start-end"
    type : str
        String indicating how overlaps are to be filteres.
        Options are 'any' and 'within'. Default is 'any'
    """

    def columns(self) -> set[str]:
        return ["gene_seq_start", "gene_seq_end", "seq_name"]

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
