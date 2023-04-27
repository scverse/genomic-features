from __future__ import annotations

import glob
import os

import ibis
import pooch
import requests
from ibis import _
from pandas import DataFrame, Timestamp
from requests.exceptions import HTTPError

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


def list_versions(species: str) -> DataFrame:
    """List available Ensembl versions for a species.

    Parameters
    ----------
    species
        The species name. E.g. Hsapiens for human, Mmusculus for mouse.

    Returns
    -------
    DataFrame
        A table of available versions in EnsDb for species of interest.
    """
    ahdb = ibis.sqlite.connect(retrieve_annotation(ANNOTATION_HUB_URL))

    # Get latest AnnotationHub timestamp
    timestamp = requests.get(TIMESTAMP_URL).text
    ahdb = ibis.sqlite.connect(retrieve_annotation(ANNOTATION_HUB_URL))
    latest_ts = Timestamp(timestamp).replace(tzinfo=None)
    cached_ts = ahdb.table("timestamp").execute()["timestamp"][0]
    if latest_ts != cached_ts:
        cached_db = glob.glob(
            os.path.join(pooch.os_cache("genomic-features"), "*annotationhub.sqlite3")
        )[0]
        os.remove(cached_db)
        ahdb = ibis.sqlite.connect(retrieve_annotation(ANNOTATION_HUB_URL))

    version_table = ahdb.table("rdatapaths").filter(_.rdataclass == "EnsDb").execute()
    version_table = version_table[
        version_table["rdatapath"].str.contains(f"EnsDb.{species}.v")
    ]
    # check that species exists
    if version_table.shape[0] == 0:
        raise ValueError(
            f"No Ensembl database found for {species}. Check species name."
        )
    else:
        version_table["Ensembl_version"] = version_table["rdatapath"].str.split(
            "/", expand=True
        )[1]
        version_table["Species"] = species
        return version_table[["Species", "Ensembl_version"]]


class EnsemblDB:
    """Ensembl annotation database."""

    def __init__(self, connection: ibis.BaseBackend):
        self.db = connection

    def genes(self):
        """Get the genes table."""
        return self.db.table("gene").execute()

    def chromosomes(self):
        """Get chromosome information."""
        return self.db.table("chromosome").execute()
