import pandas as pd

import genomic_features as gf


def test_package_has_version():
    assert gf.__version__ is not None


def test_genes():
    genes = gf.ensembl.annotation("Hsapiens", 108).genes()

    assert isinstance(genes, pd.DataFrame)
