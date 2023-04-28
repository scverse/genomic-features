import pytest

import genomic_features as gf


@pytest.fixture(scope="module")
def hsapiens108():
    return gf.ensembl.annotation("Hsapiens", 108)


def test_tables_by_degree(hsapiens108):
    result = hsapiens108.tables_by_degree()
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
    result = hsapiens108.tables_by_degree(tab=["protein", "exon"])
    assert result == ["exon", "protein"]
    result = hsapiens108.tables_by_degree(tab=["protein", "invalid_table"])
    assert result == ["protein"]


def test_list_columns(hsapiens108):
    result = hsapiens108.list_columns("gene")
    assert len(result) == len(hsapiens108.genes().columns)
    assert hsapiens108.genes().columns.isin(result).all()


def test_clean_columns(hsapiens108):
    result = hsapiens108.clean_columns(["gene_id", "gene_name"])
    assert result == ["gene_id", "gene_name"]
    result = hsapiens108.clean_columns(["gene_id", "invalid_column"])
    assert result == ["gene_id"]
    with pytest.raises(ValueError):
        hsapiens108.clean_columns(["invalid_column"])
    with pytest.raises(ValueError):
        hsapiens108.clean_columns([])


def test_tables_for_columns(hsapiens108):
    result = hsapiens108.tables_for_columns(["gene_id", "gene_name"], start_with="gene")
    assert result == ["gene"]
    result = hsapiens108.tables_for_columns(["gene_id", "tx_id"])
    assert result == ["gene", "tx"]
    result = hsapiens108.tables_for_columns(["gene_name", "invalid_column"])
    assert result == ["gene"]
    with pytest.raises(ValueError):
        hsapiens108.tables_for_columns(["invalid_column"])
    with pytest.raises(ValueError):
        hsapiens108.tables_for_columns([])


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
    # TODO: check why the number changes
    result = hsapiens108.genes(
        cols=["gene_id", "gene_name", "tx_id"],
        join_type="inner",
        filter=gf.filters.GeneBioTypeFilter(["protein_coding"]),
    )
    assert result.shape == (188553, 4)
    assert list(result.columns) == ["gene_id", "gene_name", "tx_id", "gene_biotype"]

    # table genes, transcripts and exons and filter

    result = hsapiens108.genes(
        cols=["gene_id", "gene_name", "tx_id", "exon_id"],
        join_type="inner",
        filter=gf.filters.GeneIDFilter(["ENSG00000139618"]),
    )

    assert result.shape == (186, 4)
    assert list(result.columns) == ["gene_id", "gene_name", "tx_id", "exon_id"]

    # test left join
    # table genes and transcripts
    result = hsapiens108.genes(
        cols=["gene_id", "gene_name", "protein_id"],
        join_type="left",
    )
    assert result.shape == (275721, 3)
    assert list(result.columns) == ["gene_id", "gene_name", "protein_id"]
