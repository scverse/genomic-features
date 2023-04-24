import pandas as pd

import genomic_annotations as ga


def test_package_has_version():
    assert ga.__version__ is not None


def test_genes():
    genes = ga.ensembl.annotation("Hsapiens", 108).genes()

    assert isinstance(genes, pd.DataFrame)
