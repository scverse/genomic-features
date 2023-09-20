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
        gf.ensembl.annotation("Hsapiens", 108).genes(
            columns=["tx_id"], join_type="flarb"
        )


def test_exons():
    ensdb = gf.ensembl.annotation("Hsapiens", 108)
    exons = ensdb.exons()

    pd.testing.assert_index_equal(
        exons.columns,
        pd.Index(["exon_id", "exon_seq_start", "exon_seq_end", "seq_name"]),
    )
    assert exons.shape == (888642, 4)  # Number from using R package on same DB

    exons_id = ensdb.exons(["exon_id"])

    pd.testing.assert_index_equal(exons_id.columns, pd.Index(["exon_id"]))
    assert exons_id.shape[0] == exons.shape[0]
