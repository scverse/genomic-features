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

def test_promoters():
    promoters = gf.ensembl.annotation("Hsapiens", 108).promoters()
    assert isinstance(promoters, pd.DataFrame)
    promoters = gf.ensembl.annotation("Hsapiens", 108).promoters(upstream=100, downstream=100)
    assert ((promoters.promoter_seq_end - promoters.promoter_seq_start) == 200).all()
    promoters = gf.ensembl.annotation("Hsapiens", 108).promoters(upstream=1000, downstream=100)
    assert ((promoters.promoter_seq_end - promoters.promoter_seq_start) == 1100).all()
    # Test strandedness
    promoters = gf.ensembl.annotation("Hsapiens", 108).promoters(upstream=1000, downstream=100)
    assert (promoters[promoters.seq_strand == -1].promoter_seq_start == promoters[promoters.seq_strand == -1].gene_seq_end - 100).all()
    assert (promoters[promoters.seq_strand == 1].promoter_seq_start == promoters[promoters.seq_strand == 1].gene_seq_start - 1000).all()