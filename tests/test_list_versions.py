import pandas as pd
import pytest

import genomic_features as gf


def test_ensdb_versions():
    species_versions = gf.ensembl.list_versions("Hsapiens")
    assert isinstance(species_versions, pd.DataFrame)
    species_versions = gf.ensembl.list_versions("Mmusculus")
    assert isinstance(species_versions, pd.DataFrame)
    species_versions = gf.ensembl.list_versions("Rnorvegicus")
    assert isinstance(species_versions, pd.DataFrame)

def test_missing_species():
    # test error is raised when species is not found
    with pytest.raises(ValueError):
        gf.ensembl.list_versions("Homo sapiens")

def test_unique_versions():
    # Check that different mouse strains are not considered the same species
    species_versions = gf.ensembl.list_versions("Mmusculus")
    assert species_versions.Ensembl_version.is_unique