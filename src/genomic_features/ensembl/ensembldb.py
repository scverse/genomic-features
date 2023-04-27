from __future__ import annotations

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
    db_path = retrieve_annotation(ANNOTATION_HUB_URL)
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

    def genes(
        self, filter: _filters.AbstractFilterExpr = filters.EmptyFilter()
    ) -> DataFrame:
        """Get the genes table."""
        filter.required_tables()
        # TODO: handle joins
        query = self.db.table("gene").filter(filter.convert())

        return query.execute()

    def chromosomes(self):
        """Get chromosome information."""
        return self.db.table("chromosome").execute()
