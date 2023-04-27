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


class EmptyFilter(AbstractFilterExpr):
    def convert(self) -> None:
        return None

    def required_tables(self) -> set[str]:
        return set()

    @property
    def column(self) -> str:
        return None


class AbstractFilterOperatorExpr(AbstractFilterExpr):
    def __init__(self, left: AbstractFilterExpr, right: AbstractFilterExpr):
        self.left = left
        self.right = right

    def required_tables(self) -> set[str]:
        return self.left.required_tables() & self.right.required_tables()


class AndFilterExpr(AbstractFilterOperatorExpr):
    def convert(self) -> ibis.Expr:
        return self.left.convert() & self.right.convert()


class OrFilterExpr(AbstractFilterOperatorExpr):
    def convert(self) -> ibis.Expr:
        return self.left.convert() | self.right.convert()


class AbstractFilterEqualityExpr(AbstractFilterExpr):
    def __init__(self, value: str | list[str]):
        self.value = value

    @property
    @abstractmethod
    def column(self) -> str:
        pass

    def convert(self) -> ibis.expr.deferred.Deferred:
        if isinstance(self.value, str):
            return ibis.deferred[self.column] == self.value
        else:
            return ibis.deferred[self.column].isin(self.value)


class GeneIDFilter(AbstractFilterEqualityExpr):
    @property
    def column(self) -> str:
        return "gene_id"

    def required_tables(self) -> set[str]:
        # TODO: Joining on gene_id is not necessary for transcript queries
        return {"gene"}


class GeneBioTypeFilter(AbstractFilterEqualityExpr):
    @property
    def column(self) -> str:
        return "gene_biotype"

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
