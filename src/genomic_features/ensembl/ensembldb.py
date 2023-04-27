from __future__ import annotations

import ibis
from ibis import _
from pandas import DataFrame

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
    ANNOTATION_HUB_URL = (
        "https://annotationhub.bioconductor.org/metadata/annotationhub.sqlite3"
    )
    ahdb = ibis.sqlite.connect(retrieve_annotation(ANNOTATION_HUB_URL))
    version_table = ahdb.table("rdatapaths").filter(_.rdataclass == "EnsDb").execute()
    version_table = version_table[version_table["rdatapath"].str.contains(species)]
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
