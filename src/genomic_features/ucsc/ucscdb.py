from __future__ import annotations

import warnings
from functools import cached_property
from itertools import product
import os
from pathlib import Path
from typing import Final, Literal

import ibis
import requests
from ibis import deferred
from ibis.expr.types import Table as IbisTable
from pandas import DataFrame, Timestamp
from requests.exceptions import HTTPError

from genomic_features import filters
from genomic_features._core import filters as _filters
from genomic_features._core.cache import retrieve_annotation

PKG_CACHE_DIR = "genomic-features"

BIOC_ANNOTATION_HUB_URL = (
    "https://bioconductorhubs.blob.core.windows.net/annotationhub/"
)
ANNOTATION_HUB_URL = (
    "https://annotationhub.bioconductor.org/metadata/annotationhub.sqlite3"
)
TIMESTAMP_URL = "https://annotationhub.bioconductor.org/metadata/database_timestamp"


def annotation(species: str, bioc_version: str, assembly: str,
               ucsc_table: str) -> UCSCDB:
    try:
        ucscdb = UCSCDB(
            ibis.sqlite.connect(
                retrieve_annotation(os.path.join(
                    BIOC_ANNOTATION_HUB_URL,
                    f"ucsc/standard/{bioc_version}/TxDb.{species}.UCSC.{assembly}.{ucsc_table}.sqlite"
                ))
            )
        )
    except HTTPError as err:
        if err.response.status_code == 404:
            raise ValueError(
                f"No ucsc TxDb database found for {species} {bioc_version} {assembly} {ucsc_table}. Check available versions with `genomic_features.ucsc.list_ucscdb_annotation`."
            ) from err
        else:
            raise HTTPError from err
    return ucscdb


def list_ucscdb_annotations(species: None | str | list[str] = None) -> DataFrame:
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
    _COL_ORDERS = ['species', 'assembly', 'ucsc_table', 'bioc_version']
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
        ahdb.table("rdatapaths").filter(deferred.rdataclass == "TxDb").execute()
    )
    version_table = version_table[version_table['rdatapath'].map(lambda x: x.split('/')[0] == 'ucsc')]

    version_table["bioc_version"] = (
        version_table["rdatapath"]
        .str.split("/", expand=True)[2]
    )
    version_table["species"] = (
        version_table["rdatapath"]
        .str.split("/", expand=True)[3]
        .str.split(".", expand=True)[1]
    )
    version_table["assembly"] = (
        version_table["rdatapath"]
        .str.split("/", expand=True)[3]
        .str.split(".", expand=True)[3]
    )
    version_table["ucsc_table"] = (
        version_table["rdatapath"]
        .str.split("/", expand=True)[3]
        .str.split(".", expand=True)[4]
    )
    # `Athaliana` do not follow the normal name formatting, drop them.
    version_table = version_table[version_table['ucsc_table'] != 'sqlite']

    if species is not None:
        if isinstance(species, str):
            version_table = version_table[version_table["species"] == species]
        else:
            version_table = version_table[version_table["species"].isin(species)]
        # check that species exist
        if version_table.shape[0] == 0:
            raise ValueError(
                f'No ucsc database found for {species}. Must be in {" ".join(df["species"].unique())}.'
            )

    return version_table[_COL_ORDERS].sort_values(_COL_ORDERS)


class UCSCDB:
    """UCSC annotation database."""

    def __init__(self, connection: ibis.BaseBackend):
        self.db = connection

    @cached_property
    def metadata(self) -> dict:
        metadata_tbl = self.db.table("metadata").execute()
        return dict(zip(metadata_tbl["name"], metadata_tbl["value"]))

    def __repr__(self) -> str:
        d = self.metadata
        return f"UCSCDB(organism='{d['Organism']}', ucsc_track='{d['UCSC Track']}', genome='{d['Genome']}', ucsc_table='{d['UCSC Table']}')"

    # TODO(gamazeps): should we add some info on that ? UCSC just has tx_id
    def genes(
        self,
        cols: list[str] | None = None,
        filter: _filters.AbstractFilterExpr = filters.EmptyFilter(),
        join_type: Literal["inner", "left"] = "inner",
    ) -> DataFrame:
        table: Final = "gene"
        if cols is None:
            cols = self.list_columns(table)  # get all columns

        cols = cols.copy()
        if "gene_id" not in cols:  # genes always needs gene_id
            cols.append("gene_id")

        query = self._build_query(table, cols, filter, join_type)
        return self._execute_query(query)

    def transcripts(
        self,
        cols: list[str] | None = None,
        filter: _filters.AbstractFilterExpr = filters.EmptyFilter(),
        join_type: Literal["inner", "left"] = "inner",
    ) -> DataFrame:
        table: Final = "transcript"
        if cols is None:
            cols = self.list_columns(table)  # get all columns

        cols = cols.copy()
        # Require primary key in output
        if "_tx_id" not in cols:
            cols.append("_tx_id")
        # seq_name is required for genomic range operations
        if ("tx_start" in cols or "tx_end" in cols) and "tx_chrome" not in cols:
            cols.append("tx_chrom")

        query = self._build_query(table, cols, filter, join_type)
        return self._execute_query(query)

    def exons(
        self,
        cols: list[str] | None = None,
        filter: _filters.AbstractFilterExpr = filters.EmptyFilter(),
        join_type: Literal["inner", "left"] = "inner",
    ) -> DataFrame:
        table: Final = "exon"
        if cols is None:
            cols = self.list_columns(table)  # get all columns

        cols = cols.copy()
        # Require primary key in output
        if "_exon_id" not in cols:
            cols.append("_exon_id")
        # seq_name is required for genomic range operations
        if (
            "exon_start" in cols or "exon_end" in cols
        ) and "exon_chrom" not in cols:
            cols.append("exon_chrom")

        query = self._build_query(table, cols, filter, join_type)
        return self._execute_query(query)

    def _execute_query(self, query: IbisTable) -> DataFrame:
        # TODO: Allow more options for returning results
        return query.distinct().execute()

    def chrominfo(self) -> DataFrame:
        return self.db.table("chrominfo").execute()

    def list_tables(self) -> list:
        return self.db.list_tables()

    def _tables_by_degree(self, tab: list[str] = None) -> list:
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
            "transcript": 1,
            "cds": 2,
            "gene": 3,
            "splicing": 4,
            "exon": 5,
            "chrominfo": 6,
            "metadata": 99,
        }

        return sorted(tab, key=lambda x: table_order[x])

    def list_columns(self, tables: str | list[str] | None = None) -> list[str]:
        if tables is None:
            tables = self.db.list_tables()  # list of table names
        elif isinstance(tables, str):
            tables = [tables]  # list of tables names (only one)
        columns = [c for t in tables for c in self.db.table(t).columns]
        return columns

    def _clean_columns(self, columns: list[str]) -> list[str]:
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

    def _build_query(
        self,
        table: Literal["gene", "tx", "exon"],
        cols: list[str],
        filter: _filters.AbstractFilterExpr,
        join_type: Literal["inner", "left"] = "inner",
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
        query = query.filter(filter.convert()).select(cols)
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

