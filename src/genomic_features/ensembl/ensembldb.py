from __future__ import annotations

import warnings
from collections.abc import Sequence
from functools import cached_property
from itertools import product
from pathlib import Path
from typing import Final, Literal

import ibis
import ibis.expr.types as ir
import ibis.selectors as s
import requests
from ibis import deferred
from ibis.expr.types import Table as IbisTable
from pandas import DataFrame, Timestamp
from requests.exceptions import HTTPError

from genomic_features import filters
from genomic_features._core.cache import retrieve_annotation
from genomic_features._core.filters import AbstractFilterExpr

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


def annotation(
    species: str, version: str | int, backend: Literal["duckdb", "sqlite"] = "sqlite"
) -> EnsemblDB:
    """Get an annotation database for a species and version.

    Parameters
    ----------
    species
        The species name. E.g. Hsapiens for human, Mmusculus for mouse.
    version
        The ensembl release number.
    backend
        The backend to use for the database. Either "sqlite" or "duckdb".

    Returns
    -------
    The annotation database.


    Usage
    -----
    >>> gf.ensembl.annotation("Hsapiens", "108")
    """
    try:
        sqlite_file_path = retrieve_annotation(
            ENSEMBL_URL_TEMPLATE.format(species=species, version=version)
        )

        if backend == "sqlite":
            # Connect to SQLite database
            conn = ibis.sqlite.connect(sqlite_file_path)
            ensdb = EnsemblDB(conn)
        elif backend == "duckdb":
            # Connect to DuckDB through Ibis
            conn = ibis.duckdb.connect(":memory:", extensions=["sqlite"])
            conn.attach_sqlite(sqlite_file_path)
            ensdb = EnsemblDB(conn)
        else:
            raise ValueError(f"Invalid backend: {backend}")

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
        Show gene annotations for subset of species E.g. Hsapiens for human, Mmusculus
        for mouse (optional)

    Returns
    -------
    A table of available species and annotation versions in EnsDb.


    Usage
    -----
    >>> gf.ensembl.list_ensdb_annotations("Mmusculus")
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

    version_table = (
        ahdb.table("rdatapaths").filter(deferred.rdataclass == "EnsDb").execute()
    )
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
        """Metadata for the database as a dict (e.g. who built it, when, with what resource)."""
        metadata_tbl = self.db.table("metadata").execute()
        return dict(zip(metadata_tbl["name"], metadata_tbl["value"]))

    def __repr__(self) -> str:
        d = self.metadata
        return f"EnsemblDB(organism='{d['Organism']}', ensembl_release='{d['ensembl_version']}')"

    def genes(
        self,
        cols: list[str] | None = None,
        filter: AbstractFilterExpr = filters.EmptyFilter(),
        join_type: Literal["inner", "left"] = "inner",
        order_by: Sequence[str] | str | None = None,
    ) -> DataFrame:
        """Get gene annotations.

        Parameters
        ----------
        cols
            Which columns to retrieve from the database. Can be from other tables.
            Returns all gene columns if None.
        filters
            Filters to apply to the query.
        join_type
            How to perform joins during the query (if cols or filters requires them).
        order_by
            Columns to order the results by.


        Usage
        -----
        >>> ensdb.genes(cols=["gene_id", "gene_name", "tx_id"])
        """
        table: Final = "gene"
        if cols is None:
            # TODO: check why R adds entrezid
            cols = self.list_columns(table)  # get all columns

        cols = cols.copy()
        if "gene_id" not in cols:  # genes always needs gene_id
            cols.append("gene_id")

        query = self._build_query(table, cols, filter, join_type, order_by)
        return self._execute_query(query)

    def transcripts(
        self,
        cols: list[str] | None = None,
        filter: AbstractFilterExpr = filters.EmptyFilter(),
        join_type: Literal["inner", "left"] = "inner",
        order_by: Sequence[str] | str | None = None,
    ) -> DataFrame:
        """Get transcript annotations.

        Parameters
        ----------
        cols
            Which columns to retrieve from the database. Can be from other tables.
            Returns all transcript columns if None.
        filter
            Filters to apply to the query.
        join_type
            How to perform joins during the query (if cols or filters requires them).
        order_by
            Columns to order the results by.


        Usage
        -----
        >>> ensdb.transcripts(cols=["tx_id", "tx_name", "gene_id"])
        """
        table: Final = "tx"
        if cols is None:
            cols = self.list_columns(table)  # get all columns

        cols = cols.copy()

        # Require primary key in output
        if "tx_id" not in cols:
            cols.append("tx_id")
        # seq_name is required for genomic range operations
        if ("tx_seq_start" in cols or "tx_seq_end" in cols) and "seq_name" not in cols:
            cols.append("seq_name")

        query = self._build_query(table, cols, filter, join_type, order_by)
        return self._execute_query(query)

    def exons(
        self,
        cols: list[str] | None = None,
        filter: AbstractFilterExpr = filters.EmptyFilter(),
        join_type: Literal["inner", "left"] = "inner",
        order_by: Sequence[str] | str | None = None,
    ) -> DataFrame:
        """Get exons table.

        Parameters
        ----------
        cols
            Which columns to retrieve from the database. Can be from other tables.
            Returns all exon columns if None.
        filter
            Filter to apply to the query.
        join_type
            Type of join to use for the query.


        Usage
        -----
        >>> ensdb.exons()
        """
        table: Final = "exon"
        if cols is None:
            cols = self.list_columns(table)  # get all columns

        cols = cols.copy()
        # Require primary key in output
        if "exon_id" not in cols:
            cols.append("exon_id")
        # seq_name is required for genomic range operations
        if (
            "exon_seq_start" in cols or "exon_seq_end" in cols
        ) and "seq_name" not in cols:
            cols.append("seq_name")

        query = self._build_query(table, cols, filter, join_type, order_by)
        return self._execute_query(query)

    def _execute_query(self, query: IbisTable) -> DataFrame:
        """Run a query and return the results."""
        # TODO: Allow more options for returning results
        return query.distinct().execute()

    def chromosomes(self) -> DataFrame:
        """Get chromosome information (seq_name, length, etc.).

        Usage
        -----
        >>> ensdb.chromosomes()
        """
        return self.db.table("chromosome").execute()

    def _build_query(
        self,
        table: Literal["gene", "tx", "exon"],
        cols: list[str],
        filter: AbstractFilterExpr,
        join_type: Literal["inner", "left"] = "inner",
        order_by: str
        | ir.Column
        | s.Selector
        | Sequence[str]
        | Sequence[ir.Column]
        | Sequence[s.Selector]
        | None = None,
    ) -> IbisTable:
        """Build a query for the genomic features table."""
        # Finalize cols
        self._clean_columns(cols)
        for col in filter.columns():
            if col not in cols:
                cols.append(col)

        # check if join is required
        tables = self._get_required_tables(self._tables_for_columns(cols))

        # Basically just to make sure exons stay in the query
        if table not in tables:
            tables.append(table)

        if len(tables) > 1:
            query = self._join_query(tables, start_with=table, join_type=join_type)
        else:
            query = self.db.table(table)

        # add filter
        query = query.filter(filter.convert())
        query = query.select(cols)

        if order_by is not None:
            # Custom ordering is provided
            query = query.order_by(order_by)
        else:
            # Default ordering
            order_by = []
            if "seq_name" in cols:
                order_by = ["seq_name"]
            if "gene_seq_start" in cols:
                order_by.extend(["gene_seq_start"])
            if "tx_seq_start" in cols:
                order_by.extend(["tx_seq_start"])
            if "exon_seq_start" in cols:
                order_by.extend(["exon_seq_start"])

            order_by.extend([c for c in cols if "id" in c])
            query = query.order_by(order_by)

        return query

    def _join_query(
        self,
        tables: list[str],
        start_with: str,
        join_type: Literal["inner", "left"] = "inner",
    ) -> IbisTable:
        """Join tables and return a query."""
        # check for intermediate tables
        JOIN_TABLE = [
            (("gene", "tx"), "gene_id"),
            (("gene", "chromosome"), "seq_name"),
            (("tx", "tx2exon"), "tx_id"),
            (("tx2exon", "exon"), "exon_id"),
            (("tx", "protein"), "tx_id"),
            (("gene", "entrezgene"), "gene_id"),
            (("protein", "protein_domain"), "protein_id"),
            (("protein", "uniprot"), "protein_id"),
            (("uniprot", "protein_domain"), "protein_id"),
        ]
        tables = tables.copy()
        tables.remove(start_with)
        db = self.db
        current_tables = [start_with]
        query = db.table(start_with)

        while len(tables) > 0:
            for (table_names, key), t1_name, t2_name in product(  # noqa: B007
                JOIN_TABLE, current_tables, tables
            ):
                if t1_name in table_names and t2_name in table_names:
                    break
            else:
                raise ValueError(
                    f"Failed to find match for tables: {current_tables} and {tables}"
                )

            current_tables.append(t2_name)
            tables.remove(t2_name)

            t2 = db.table(t2_name)
            if join_type == "inner":
                query = query.join(t2, predicates=[key], how="inner")
            elif join_type == "left":
                query = query.join(
                    t2,
                    predicates=[key],
                    how="left",
                    rname="{name}_y",
                    # suffixes=("", "_y"),
                )
                query = query.drop(f"{key}_y")  # drop duplicate columns
            else:
                raise ValueError(f"Invalid join type: {join_type}")

        return query

    def list_tables(self) -> list:
        """List all tables available in the genomic features database."""
        return self.db.list_tables()

    def _tables_by_degree(self, tab: list[str] = None) -> list:
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

    def _get_required_tables(self, tab) -> list:
        """Given tables, get all intermediate tables required to execute the query."""
        # If we have exon and any other table, we need definitely tx2exon
        if "exon" in tab and len(tab) > 1:
            tab = list(set(tab + ["tx2exon"]))

        # If we have chromosome and any other table, we'll need gene
        if "seq_name" in tab and len(tab) > 1:
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

        return self._tables_by_degree(tab)

    def list_columns(self, tables: str | list[str] | None = None) -> list[str]:
        """List all columns available in the genomic features table."""
        if tables is None:
            tables = self.db.list_tables()  # list of table names
        elif isinstance(tables, str):
            tables = [tables]  # list of tables names (only one)
        columns = [c for t in tables for c in self.db.table(t).columns]
        return columns

    def _clean_columns(self, columns: list[str]) -> list[str]:
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

    def _tables_for_columns(self, cols: list, start_with: str | None = None) -> list:
        """
        Return a list of tables that contain the specified columns.

        Parameters
        ----------
        cols
            Columns that we're looking for.
        """
        cols = self._clean_columns(cols)
        table_list = self._tables_by_degree()  # list of table names

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
