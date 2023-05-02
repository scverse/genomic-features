from __future__ import annotations

import warnings
from functools import cached_property
from pathlib import Path

import ibis
import requests
from ibis import _
from pandas import DataFrame, Timestamp
from requests.exceptions import HTTPError

from genomic_features import filters
from genomic_features._core import filters as _filters
from genomic_features._core.cache import retrieve_annotation

PKG_CACHE_DIR = "genomic-features"

BIOC_ANNOTATION_HUB_URL = (
    "https://bioconductorhubs.blob.core.windows.net/annotationhub/"
)
ENSEMBL_URL_TEMPLATE = (
    BIOC_ANNOTATION_HUB_URL + "AHEnsDbs/v{version}/EnsDb.{species}.v{version}.sqlite"
)
ANNOTATION_HUB_URL = (
    "https://annotationhub.bioconductor.org/metadata/annotationhub.sqlite3"
)
TIMESTAMP_URL = "https://annotationhub.bioconductor.org/metadata/database_timestamp"


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
    try:
        ensdb = EnsemblDB(
            ibis.sqlite.connect(
                retrieve_annotation(
                    ENSEMBL_URL_TEMPLATE.format(species=species, version=version)
                )
            )
        )
    except HTTPError as err:
        if err.response.status_code == 404:
            raise ValueError(
                f"No Ensembl database found for {species} v{version}. Check available versions with `genomic_features.ensembl.list_versions`."
            ) from err
        else:
            raise HTTPError from err
    return ensdb


def list_ensdb_annotations(species: None | str | list[str] = None) -> DataFrame:
    """List available Ensembl gene annotations.

    Parameters
    ----------
    species
        Show gene annotations for subset of species E.g. Hsapiens for human, Mmusculus for mouse (optional)

    Returns
    -------
    DataFrame
        A table of available species and annotation versions in EnsDb.
    """
    # Get latest AnnotationHub timestamp
    db_path = Path(retrieve_annotation(ANNOTATION_HUB_URL))
    timestamp = requests.get(TIMESTAMP_URL).text
    ahdb = ibis.sqlite.connect(db_path)
    latest_ts = Timestamp(timestamp).replace(tzinfo=None)
    cached_ts = ahdb.table("timestamp").execute()["timestamp"][0]
    if latest_ts != cached_ts:
        db_path.unlink()
        ahdb = ibis.sqlite.connect(retrieve_annotation(ANNOTATION_HUB_URL))

    version_table = ahdb.table("rdatapaths").filter(_.rdataclass == "EnsDb").execute()
    version_table["Species"] = (
        version_table["rdatapath"]
        .str.split("/", expand=True)[2]
        .str.split(".", expand=True)[1]
    )
    if species is not None:
        if isinstance(species, str):
            version_table = version_table[version_table["Species"] == species]
        else:
            version_table = version_table[version_table["Species"].isin(species)]
        # check that species exist
        if version_table.shape[0] == 0:
            raise ValueError(
                f"No Ensembl database found for {species}. Check species name."
            )

    version_table["Ensembl_version"] = version_table["rdatapath"].str.split(
        "/", expand=True
    )[1]
    version_table["Ensembl_version"] = (
        version_table["Ensembl_version"].str.replace("v", "").astype(int)
    )
    return version_table[["Species", "Ensembl_version"]].sort_values(
        ["Species", "Ensembl_version"]
    )


class EnsemblDB:
    """Ensembl annotation database."""

    def __init__(self, connection: ibis.BaseBackend):
        self.db = connection

    @cached_property
    def metadata(self) -> dict:
        """Metadata for the database (e.g. who built it, when, with what resource)."""
        metadata_tbl = self.db.table("metadata").execute()
        return dict(zip(metadata_tbl["name"], metadata_tbl["value"]))

    def __repr__(self) -> str:
        d = self.metadata
        return f"EnsemblDB(organism='{d['Organism']}', ensembl_release='{d['ensembl_version']}')"

    def genes(
        self,
        cols: list = None,
        filter: _filters.AbstractFilterExpr = filters.EmptyFilter(),
        join_type: str = "inner",
    ) -> DataFrame:
        """Get the genes table."""
        table = "gene"
        if cols is None:
            # TODO: check why R adds entrezid
            cols = self.list_columns(table)  # get all columns

        return_cols = cols
        cols = list(set(cols + ["gene_id"]))  # genes always needs gene_id
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
        tables = self.get_required_tables(self.tables_for_columns(cols))
        if len(tables) > 1:
            query = self.join_query(tables, start_with=table, join_type=join_type)
        else:
            query = self.db.table(table)
        # add filter
        query = query.filter(filter.convert()).select(cols)
        return query

    def join_query(self, tables, start_with, join_type="inner"):
        """Join tables and return a query."""
        # check for intermediate tables

        tables.remove(start_with)
        tables = [start_with] + tables
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
            missing_tables = ", ".join(set(tab) - set(self.list_tables()))
            warnings.warn(
                f"The following tables are not in the database: {missing_tables}.",
                UserWarning,
                stacklevel=2,
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

    def get_required_tables(self, tab):
        # If we have exon and any other table, we need definitely tx2exon
        if "exon" in tab and len(tab) > 1:
            tab = list(set(tab + ["tx2exon"]))

        # If we have chromosome and any other table, we'll need gene
        if "chromosome" in tab and len(tab) > 1:
            tab = list(set(tab + ["gene"]))

        # If we have exon and we have gene, we'll need also tx
        if ("exon" in tab or "tx2exon" in tab) and "gene" in tab:
            tab = list(set(tab + ["tx"]))

        # Resolve the proteins: need tx to map between proteome and genome
        if any(t in ["uniprot", "protein_domain", "protein"] for t in tab) and any(
            t in ["exon", "tx2exon", "gene", "chromosome", "entrezgene"] for t in tab
        ):
            tab = list(set(tab + ["tx"]))

        # Need protein.
        if any(t in ["uniprot", "protein_domain"] for t in tab) and any(
            t in ["exon", "tx2exon", "tx", "gene", "chromosome", "entrezgene"]
            for t in tab
        ):
            tab = list(set(tab + ["protein"]))

        # entrezgene is only linked via gene
        if "entrezgene" in tab and len(tab) > 1:
            tab = list(set(tab + ["gene"]))

        return self.tables_by_degree(tab)

    def list_columns(self, tables=None) -> list:
        """List all columns available in the genomic features table."""
        if tables is None:
            tables = self.db.list_tables()  # list of table names
        elif isinstance(tables, str):
            tables = [tables]  # list of tables names (only one)
        columns = []
        columns = [c for t in tables for c in self.db.table(t).columns]
        return columns

    def clean_columns(self, columns) -> list:
        """Clean a list of columns to make sure they are valid."""
        if isinstance(columns, str):
            columns = [columns]

        valid_columns = set(self.list_columns())
        cols = list(filter(lambda c: c in valid_columns, columns))
        invalid_columns = set(columns) - valid_columns
        if invalid_columns:
            raise ValueError(
                f"The following columns are not found in any database: {invalid_columns}"
            )
        if not cols:
            raise ValueError("No valid columns were found.")
        return cols

    def tables_for_columns(self, cols: list, start_with: str = None) -> list:
        """
        Return a list of tables that contain the specified columns.

        Parameters
        ----------
        cols
            Columns that we're looking for.
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
