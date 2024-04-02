import pandas as pd
import pytest

import genomic_features as gf

SPECIES = [
    "Scerevisiae",
    "Mmusculus_nzohlltj",
    "Ggallus",
    "Hsapiens",
    "Mmusculus",
    "Celegans",
]


def test_ensdb_versions():
    species_versions = gf.ensembl.list_ensdb_annotations("Hsapiens")
    assert isinstance(species_versions, pd.DataFrame)
    species_versions = gf.ensembl.list_ensdb_annotations("Mmusculus")
    assert isinstance(species_versions, pd.DataFrame)
    species_versions = gf.ensembl.list_ensdb_annotations("Rnorvegicus")
    assert isinstance(species_versions, pd.DataFrame)


def test_missing_species():
    # test error is raised when species is not found
    with pytest.raises(ValueError):
        gf.ensembl.list_ensdb_annotations("Homo sapiens")


def test_unique_versions():
    # Check that different mouse strains are not considered the same species
    species_versions = gf.ensembl.list_ensdb_annotations("Mmusculus")
    assert species_versions.Ensembl_version.is_unique


def test_all_species():
    species_versions = gf.ensembl.list_ensdb_annotations()
    assert pd.Series(SPECIES).isin(species_versions["Species"]).all()
