import pytest

import genomic_features as gf
from genomic_features import filters


@pytest.fixture(scope="module")
def hsapiens108():
    return gf.ensembl.annotation("Hsapiens", 108)


@pytest.mark.parametrize(
    "filt",
    [
        filters.GeneIDFilter("ENSG00000000003"),
        filters.GeneIDFilter("ENSG00000000460"),
        filters.GeneIDFilter("LRG_997"),
        filters.GeneBioTypeFilter("protein_coding"),
        filters.GeneBioTypeFilter("TR_C_gene"),
    ],
)
def test_equality_filter_single(hsapiens108, filt):
    result = hsapiens108.genes(filter=filt)[list(filt.columns())[0]]
    assert set(result) == {filt.value}


@pytest.mark.parametrize(
    "filt",
    [
        filters.GeneIDFilter(["ENSG00000000003", "ENSG00000093183"]),
        filters.GeneBioTypeFilter(["TR_J_gene", "TR_V_gene"]),
    ],
)
def test_equality_filter_list(hsapiens108, filt):
    result = hsapiens108.genes(filter=filt)[list(filt.columns())[0]]
    assert set(result) == set(filt.value)


# These are not working quite as expected:
# https://github.com/ibis-project/ibis/issues/6096
def test_and_filter(hsapiens108):
    assert (
        hsapiens108.genes(
            filter=(
                filters.GeneBioTypeFilter("protein_coding")
                & filters.GeneBioTypeFilter("TR_C_gene")
            )
        ).shape[0]
        == 0
    )
    assert (
        hsapiens108.genes(
            filter=(
                filters.GeneBioTypeFilter("protein_coding")
                & filters.GeneIDFilter(
                    ["LRG_997", "ENSG00000000460", "ENSG00000000003"]
                )
            )
        ).shape[0]
        == 2
    )


def test_or_filter(hsapiens108):
    assert (
        hsapiens108.genes(
            filter=(
                filters.GeneBioTypeFilter("protein_coding")
                | filters.GeneBioTypeFilter("TR_C_gene")
            )
        ).shape[0]
        == hsapiens108.genes()["gene_biotype"]
        .isin(["protein_coding", "TR_C_gene"])
        .sum()
    )
    assert (
        hsapiens108.genes(
            filter=(
                filters.GeneIDFilter("LRG_997")
                | filters.GeneIDFilter(
                    ["LRG_997", "ENSG00000000460", "ENSG00000000003"]
                )
            )
        ).shape[0]
        == 3
    )


def test_range_filter(hsapiens108):
    any_overlap_filter = hsapiens108.genes(
        filter=filters.RangesFilter("1:77000000-78000000")
    )
    within_overlap_filter = hsapiens108.genes(
        filter=filters.RangesFilter("1:77000000-78000000", type="within")
    )
    assert all(within_overlap_filter.seq_name == "1") & all(
        any_overlap_filter.seq_name == "1"
    )
    assert any_overlap_filter.shape[0] > within_overlap_filter.shape[0]
    assert (all(within_overlap_filter.gene_seq_start > 77000000)) & (
        all(within_overlap_filter.gene_seq_end < 78000000)
    )
    assert (all(any_overlap_filter.gene_seq_end >= 77000000)) & (
        all(any_overlap_filter.gene_seq_start <= 78000000)
    )
