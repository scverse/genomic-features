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
