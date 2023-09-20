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
    def columns(self) -> set[str]:
        """Columns required by this filter."""
        pass


class EmptyFilter(AbstractFilterExpr):
    """An empty filter that returns all rows."""

    def __repr__(self) -> str:
        return "EmptyFilter()"

    def convert(self) -> None:
        return None

    def columns(self) -> set[str]:
        return set()


class AbstractFilterOperatorExpr(AbstractFilterExpr):
    def __init__(self, left: AbstractFilterExpr, right: AbstractFilterExpr):
        self.left = left
        self.right = right

    def columns(self) -> set[str]:
        return self.left.columns() | self.right.columns()


class AndFilterExpr(AbstractFilterOperatorExpr):
    """Logical and of two filters."""

    def __repr__(self) -> str:
        return f"({self.left} & {self.right})"

    def convert(self) -> ibis.expr.deferred.Deferred:
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


class OrFilterExpr(AbstractFilterOperatorExpr):
    """Logical or of two filters."""

    def __repr__(self) -> str:
        return f"({self.left} | {self.right})"

    def convert(self) -> ibis.expr.deferred.Deferred:
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


class AbstractFilterRangeExpr(AbstractFilterExpr, ABC):
    def __init__(self, value: str, type: str = "any"):
        self.value = value
        self.type = type

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({self.value})"

    @property
    @abstractmethod
    def _range_columns(self) -> list[str]:
        pass

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
    """Filter by gene_id."""

    def columns(self) -> set[str]:
        return {"gene_id"}


class GeneBioTypeFilter(AbstractFilterEqualityExpr):
    """Filter by gene_biotype."""

    def columns(self) -> set[str]:
        return {"gene_biotype"}


class TxIDFilter(AbstractFilterEqualityExpr):
    """Filter by tx_id column."""

    def columns(self) -> set[str]:
        return {"tx_id"}


class TxBioTypeFilter(AbstractFilterEqualityExpr):
    """Filter by tx_biotype column."""

    def columns(self) -> set[str]:
        return {"tx_biotype"}


class ExonIDFilter(AbstractFilterEqualityExpr):
    """Filter by exon_id column."""

    def columns(self) -> set[str]:
        return {"exon_id"}


class GeneNameFilter(AbstractFilterEqualityExpr):
    """Filter by gene_name."""

    def columns(self) -> set[str]:
        return {"gene_name"}


class GeneRangeFilter(AbstractFilterRangeExpr):
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


class SeqNameFilter(AbstractFilterEqualityExpr):
    """Filter by seq_name (e.g. chromosome).

    Usage
    -----

    >>> ensdb.genes(filter=gf.filters.SeqNameFilter("MT"))
    """

    def __init__(self, value: str | int | list):
        if isinstance(value, int):
            value = str(value)
        elif isinstance(value, str):
            pass
        else:
            orig_value = value
            value = [str(v) for v in orig_value]

        self.value = value

    def columns(self) -> set[str]:
        return {"seq_name"}


class CanonicalTxFilter(AbstractFilterExpr):
    """Filter for canonical transcripts.

    Usage
    -----

    >>> ensdb.transcripts(filter=gf.filters.CanonicalTxFilter())
    >>> ensdb.exons(
    ...     columns=["tx_id", "exon_id", "seq_name", "exon_seq_start", "exon_seq_end"],
    ...     filter=gf.filters.CanonicalTxFilter()
    ... )
    """

    def __init__(self):
        pass

    def __repr__(self) -> str:
        return "CanonicalTxFilter()"

    def columns(self) -> set[str]:
        return {"tx_is_canonical"}

    def convert(self) -> ibis.expr.deferred.Deferred:
        return ibis.deferred["tx_is_canonical"] == 1


class UniProtIDFilter(AbstractFilterEqualityExpr):
    """Filter by UniProt ID.

    Usage
    -----

    >>> ensdb.genes(filter=gf.filters.UniProtIDFilter("P12345"))
    """

    def columns(self) -> set[str]:
        return {"uniprot_id"}


class UniProtDBFilter(AbstractFilterEqualityExpr):
    """Filter by UniProt database.

    Usage
    -----

    >>> ensdb.genes(filter=gf.filters.UniProtDBFilter("SWISSPROT"))
    """

    def columns(self) -> set[str]:
        return {"uniprot_db"}


class UniProtMappingTypeFilter(AbstractFilterEqualityExpr):
    """Filter by UniProt mapping type.

    Usage
    -----

    >>> ensdb.genes(filter=gf.filters.UniProtMappingTypeFilter("DIRECT"))
    """

    def columns(self) -> set[str]:
        return {"uniprot_mapping_type"}
