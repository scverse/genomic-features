from __future__ import annotations

import ibis

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

    def genes(self):
        """Get the genes table."""
        return self.db.table("gene").execute()

    def chromosomes(self):
        """Get chromosome information."""
        return self.db.table("chromosome").execute()
