import pandas as pd
import pytest

import genomic_features as gf


def test_package_has_version():
    assert gf.__version__ is not None


def test_genes():
    genes = gf.ensembl.annotation("Hsapiens", 108).genes()
    assert isinstance(genes, pd.DataFrame)


def test_missing_version():
    with pytest.raises(ValueError):
        gf.ensembl.annotation("Hsapiens", 86)


def test_repr():
    result = repr(gf.ensembl.annotation("Hsapiens", 108))
    expected = "EnsemblDB(organism='Homo sapiens', ensembl_release='108')"

    assert result == expected


def test_invalid_join():
    with pytest.raises(ValueError, match=r"Invalid join type: flarb"):
        gf.ensembl.annotation("Hsapiens", 108).genes(cols=["tx_id"], join_type="flarb")
