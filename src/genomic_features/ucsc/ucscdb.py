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

_TX_TABLE = 'transcript'
_EXONS_TABLE = 'exon'
_GENES_TABLE = 'gene'

_PRETTY_NAMES = {
    '_tx_id': 'tx_id',
    'tx_chrom': 'chrom',
    'tx_strand': 'strand',
    'tx_start': 'start',
    'tx_end': 'end',
    '_exon_id': 'exon_id',
    'exon_chrom': 'chrom',
    'exon_strand': 'strand',
    'exon_start': 'start',
    'exon_end': 'end',
}

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

    def chrominfo(self) -> DataFrame:
        return self.db.table("chrominfo").execute()

    def list_tables(self) -> list:
        return self.db.list_tables()

    def transcripts(
        self,
        #cols: list[str] | None = None,
        #filter: _filters.AbstractFilterExpr = filters.EmptyFilter(),
    ) -> DataFrame:
        tx = self.db.table(_TX_TABLE).execute()
        tx = tx.rename(columns=_PRETTY_NAMES)
        tx = tx.drop('tx_type', axis=1) # always None
        return tx

    def exons(
        self,
        #cols: list[str] | None = None,
        #filter: _filters.AbstractFilterExpr = filters.EmptyFilter(),
    ) -> DataFrame:
        exons = self.db.table(_EXONS_TABLE).execute()
        exons = exons.rename(columns=_PRETTY_NAMES)
        exons = exons.drop('exon_name', axis=1) # always None
        return exons

    def genes(
        self,
        #cols: list[str] | None = None,
        #filter: _filters.AbstractFilterExpr = filters.EmptyFilter(),
    ) -> DataFrame:
        genes = self.db.table(_GENES_TABLE).execute()
        return genes

    def _execute_query(self, query: IbisTable) -> DataFrame:
        # TODO: Allow more options for returning results
        return query.distinct().execute()

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
