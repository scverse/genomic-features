import pandas as pd
import pytest

import genomic_features as gf


@pytest.fixture(scope="module")
def hsapiens108():
    return gf.ensembl.annotation("Hsapiens", 108)


def test_tables_by_degree(hsapiens108):
    result = hsapiens108._tables_by_degree()
    assert result == [
        "gene",
        "tx",
        "tx2exon",
        "exon",
        "chromosome",
        "protein",
        "uniprot",
        "protein_domain",
        "entrezgene",
        "metadata",
    ]
    result = hsapiens108._tables_by_degree(tab=["protein", "exon"])
    assert result == ["exon", "protein"]
    result = hsapiens108._tables_by_degree(tab=["protein", "invalid_table"])
    assert result == ["protein"]


def test_list_columns(hsapiens108):
    result = hsapiens108.list_columns("gene")
    assert result == list(hsapiens108.db.table("gene").columns)


def test_clean_columns(hsapiens108):
    result = hsapiens108._clean_columns("gene_id")
    assert result == ["gene_id"]
    result = hsapiens108._clean_columns(["gene_id", "gene_name"])
    assert result == ["gene_id", "gene_name"]
    with pytest.raises(ValueError):
        hsapiens108._clean_columns(["gene_id", "invalid_column"])
    with pytest.raises(ValueError):
        hsapiens108._clean_columns([])


def test_tables_for_columns(hsapiens108):
    result = hsapiens108._tables_for_columns(["gene_id"])
    assert result == ["gene"]


def test_required_tables(hsapiens108):
    result = hsapiens108._get_required_tables(["gene", "tx"])
    assert result == ["gene", "tx"]
    # case where we need intermediate tables
    result = hsapiens108._get_required_tables(["gene", "protein"])
    assert result == ["gene", "tx", "protein"]


# Test simple subsetting to columns in one table gene
def test_simple_subsetting(hsapiens108):
    result = hsapiens108.genes(cols=["gene_id", "gene_name"])
    assert result.shape == (70616, 2)
    assert result.columns.tolist() == ["gene_id", "gene_name"]


# Test subsetting to columns in multiple tables
def test_multiple_table_subsetting(hsapiens108):
    # table genes and transcripts
    result = hsapiens108.genes(
        cols=["gene_id", "gene_name", "tx_id"],
        join_type="inner",
    )
    assert result.shape == (275721, 3)
    assert list(result.columns) == ["gene_id", "gene_name", "tx_id"]

    # table genes and transcripts with filter
    result = hsapiens108.genes(
        cols=["gene_id", "gene_name", "tx_id"],
        join_type="inner",
        filter=gf.filters.GeneBioTypeFilter(["protein_coding"]),
    )
    assert result.shape == (185904, 4)
    assert list(result.columns) == ["gene_id", "gene_name", "tx_id", "gene_biotype"]

    # table genes, transcripts and exons and filter
    result = hsapiens108.genes(
        cols=["gene_id", "gene_name", "exon_id"],
        join_type="inner",
        filter=gf.filters.GeneIDFilter(["ENSG00000139618"]),
    )

    assert result.shape == (78, 3)
    assert list(result.columns) == ["gene_id", "gene_name", "exon_id"]

    # test left join
    # table genes and transcripts
    result = hsapiens108.genes(
        cols=["gene_id", "gene_name", "protein_id"],
        join_type="left",
        filter=gf.filters.GeneBioTypeFilter(["protein_coding"]),
    )
    assert result.shape == (135834, 4)
    assert list(result.columns) == [
        "gene_id",
        "gene_name",
        "protein_id",
        "gene_biotype",
    ]
    assert (
        result.loc[result.gene_biotype != "protein_coding", "protein_id"].isna().all()
    )


def test_chromosome_columns(hsapiens108):
    # https://github.com/scverse/genomic-features/pull/44/files#r1196331705
    result = hsapiens108.genes(cols=["gene_id", "seq_name", "seq_length"])
    assert result.shape[0] == hsapiens108.db.table("gene").count().execute()

    chroms = hsapiens108.chromosomes()
    expected_lengths = (
        chroms.set_index("seq_name")["seq_length"]
        .loc[result["seq_name"]]
        .reset_index(drop=True)
    )
    pd.testing.assert_series_equal(result["seq_length"], expected_lengths)
