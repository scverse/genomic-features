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
        join_type: str = "inner",
    ) -> pd.DataFrame:
        """Get the genes table."""
        table = "gene"
        if cols is None:
            # TODO: check why R adds entrezid
            cols = self.list_columns(table)  # get all columns

        return_cols = cols
        cols = self.clean_columns(cols)
        filter.required_tables()
        cols = list(set(cols) | filter.columns())  # add columns from filter

        query = self.build_query(table, cols, filter, join_type)

        result = query.execute()

        # order columns
        return_cols = return_cols + [col for col in cols if col not in return_cols]
        result = result[return_cols]
        return result

    def chromosomes(self):
        """Get chromosome information."""
        return self.db.table("chromosome").execute()

    def build_query(self, table, cols, filter, join_type="inner"):
        """Build a query for the genomic features table."""
        # check if join is required
        if len(self.tables_for_columns(cols)) > 1:
            query = self.join_query(cols, table, join_type)
        else:
            query = self.db.table(table)
        # add filter
        query = query.filter(filter.convert()).select(cols)
        return query

    def join_query(self, cols, start_with, join_type="inner"):
        """Join tables and return a query."""
        # check tables for join
        tables = self.tables_for_columns(cols, start_with=start_with)
        # TODO: check what ordering by degree does

        if join_type == "inner":
            query = self.inner_join_query(tables)
        elif join_type == "left":
            query = self.left_join_query(tables)
        else:
            raise ValueError(f"Invalid join type: {join_type}")

        return query

    def inner_join_query(self, tables):
        """Join tables in inner join and return a query."""
        # build query
        query = self.db.table(tables[0])
        for t in tables[1:]:
            # determine common columns
            # TODO: Determine if we want to use all common columns
            common_cols = list(set(query.columns) & set(self.db.table(t).columns))
            query = query.join(self.db.table(t), predicates=common_cols, how="inner")
        return query

    def left_join_query(self, tables):
        """Join tables in left join and return a query."""
        # build query
        query = self.db.table(tables[0])
        for t in tables[1:]:
            # determine common columns
            common_cols = list(set(query.columns) & set(self.db.table(t).columns))
            query = query.join(
                self.db.table(t),
                predicates=common_cols,
                how="left",
                suffixes=("", "_y"),
            )
            query = query.drop(
                *[f"{c}_y" for c in common_cols]
            )  # drop duplicate columns
        return query

    def list_tables(self) -> list:
        """List all tables available in the genomic features database."""
        return self.db.list_tables()

    def tables_by_degree(self, tab: list[str] = None) -> list:
        """Order tables available in the genomic features database."""
        if tab is None:
            tab = self.list_tables()  # list of table names
        # check that all tables are in the database and print warning
        if not set(tab).issubset(set(self.list_tables())):
            print(
                f"Warning: tables {set(tab) - set(self.list_tables())} are not in the database."
            )
            tab = list(set(tab) & set(self.list_tables()))  # remove tables not in db

        # order tables

        table_order = {
            "gene": 1,
            "tx": 2,
            "tx2exon": 3,
            "exon": 4,
            "chromosome": 5,
            "protein": 6,
            "uniprot": 7,
            "protein_domain": 8,
            "entrezgene": 9,
            "metadata": 99,
        }

        return sorted(tab, key=lambda x: table_order[x])

    def list_columns(self, table=None) -> list:
        """List all columns available in the genomic features table."""
        if table is None:
            table = self.db.list_tables()  # list of table names
        else:
            table = [table]  # list of tables names (only one)
        columns = []
        for t in table:
            columns += self.db.table(t).columns
        return columns

    def clean_columns(self, columns) -> list:
        """Clean a list of columns to make sure they are valid."""
        valid_columns = self.list_columns()
        columns = [c for c in columns if c in valid_columns]
        unvalid_columns = [c for c in columns if c not in valid_columns]

        if len(unvalid_columns) > 0:
            print("The following columns are not valid and will be ignored:")
            print(unvalid_columns)

        if len(columns) == 0:
            raise ValueError("No valid columns were found.")
        return columns

    def tables_for_columns(self, cols: list, start_with: str = None) -> list:
        """
        Return a list of tables that contain the specified columns.

        Parameters
        ----------
        cols
            Columns that we're looking for
        start_with
            Primary table being queried
        """
        cols = self.clean_columns(cols)
        table_list = self.tables_by_degree()  # list of table names

        # remove start_with from table_list and add it to the beginning of the list
        if start_with is not None:
            # check if start_with is a valid table
            if start_with not in table_list:
                raise ValueError(f"Invalid table: {start_with}")
            # remove start_with from table_list and add it to the beginning of the list
            table_list.remove(start_with)
            table_list = [start_with] + table_list

        tables = []
        for t in table_list:
            # check if all columns are in one table
            if set(cols).issubset(self.db.table(t).columns):
                tables.append(t)
                return tables
            else:
                # check if a single column is in the table
                for c in cols.copy():
                    if c in self.db.table(t).columns:
                        if t not in tables:
                            tables.append(t)
                        cols.remove(c)  # remove column from list
        return tables
