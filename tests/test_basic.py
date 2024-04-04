import pandas as pd
import pytest

import genomic_features as gf


def test_package_has_version():
    assert gf.__version__ is not None


@pytest.mark.parametrize("backend", ["sqlite", "duckdb"])
def test_genes(backend):
    genes = gf.ensembl.annotation("Hsapiens", 108, backend=backend).genes()
    assert isinstance(genes, pd.DataFrame)

    # Test sort order
    genes_resorted = genes.sort_values(
        ["seq_name", "gene_seq_start", "gene_id"]
    ).reset_index(drop=True)
    pd.testing.assert_frame_equal(genes, genes_resorted)


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


@pytest.mark.parametrize("backend", ["sqlite", "duckdb"])
def test_exons(backend):
    ensdb = gf.ensembl.annotation("Hsapiens", 108, backend=backend)
    exons = ensdb.exons()

    pd.testing.assert_index_equal(
        exons.columns,
        pd.Index(["exon_id", "exon_seq_start", "exon_seq_end", "seq_name"]),
    )
    assert exons.shape == (888642, 4)  # Number from using R package on same DB

    exons_id = ensdb.exons(["exon_id"])

    pd.testing.assert_index_equal(exons_id.columns, pd.Index(["exon_id"]))
    assert exons_id.shape[0] == exons.shape[0]

    # Test sort order
    exons_resorted = exons.sort_values(
        ["seq_name", "exon_seq_start", "exon_id"]
    ).reset_index(drop=True)
    pd.testing.assert_frame_equal(exons, exons_resorted)
