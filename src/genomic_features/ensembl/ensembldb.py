from __future__ import annotations

import ibis
import pandas as pd

from genomic_features._core import filters
from genomic_features._core.cache import retrieve_annotation


BIOC_ANNOTATION_HUB_URL = (
    "https://bioconductorhubs.blob.core.windows.net/annotationhub/"
)
ENSEMBL_URL_TEMPLATE = (
    BIOC_ANNOTATION_HUB_URL + "AHEnsDbs/v{version}/EnsDb.{species}.v{version}.sqlite"
)


def annotation(species: str, version: str | int):
    """Get an annotation database for a species and version.

    Parameters
    ----------
    species
        The species name. E.g. Hsapiens for human, Mmusculus for mouse.
    version
        The ensembl release number.

    Returns
    -------
    EnsemblDB
        The annotation database.
    """
    return EnsemblDB(
        ibis.sqlite.connect(
            retrieve_annotation(
                ENSEMBL_URL_TEMPLATE.format(species=species, version=version)
            )
        )
    )


class EnsemblDB:
    """Ensembl annotation database."""

    def __init__(self, connection: ibis.BaseBackend):
        self.db = connection

    def genes(
        self,
        cols: list = None,
        filter: filters.AbstractFilterExpr = filters.EmptyFilter(),
    ) -> pd.DataFrame:
        """Get the genes table."""
        if cols is None:
            # TODO: check why R adds entrezid
            cols = self.listColumns("gene")  # get all columns
        filter.required_tables()

        cols += [filter.column]  # add columns from filter
        cols = self.cleanColumns(cols)  # clean columns
        # TODO: handle joins
        query = self.db.table("gene").filter(filter.convert()).select(cols)

        return query.execute()

    def chromosomes(self):
        """Get chromosome information."""
        return self.db.table("chromosome").execute()

    def listColumns(self, table=None) -> list:
        """
        List all columns available in the genomic features table.
        """

        if table is None:
            table = self.db.list_tables()  # list of table names
        else:
            table = [table]  # list of tables names (only one)
        columns = []
        for t in table:
            columns += self.db.table(t).columns
        return columns

    def cleanColumns(self, columns) -> list:
        """
        Clean a list of columns to make sure they are valid.
        """
        valid_columns = self.listColumns()
        columns = [c for c in columns if c in valid_columns]
        unvalid_columns = [c for c in columns if c not in valid_columns]

        if len(unvalid_columns) > 0:
            print("The following columns are not valid and will be ignored:")
            print(unvalid_columns)

        if len(columns) == 0:
            raise ValueError("No valid columns were found.")
        return columns

    def tablesForColumns(self, cols) -> list:
        """
        Return a list of tables that contain the specified columns.
        """
        columns = self.cleanColumns(cols)
        tables = []
        for t in self.db.list_tables():
            if set(columns).issubset(self.db.table(t).columns):
                tables.append(t)
        return tables
