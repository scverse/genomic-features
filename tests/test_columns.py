import pytest

import genomic_features as gf


@pytest.fixture(scope="module")
def hsapiens108():
    return gf.ensembl.annotation("Hsapiens", 108)


def test_columns_subset(hsapiens108):
    result = hsapiens108.genes(cols=["gene_id", "gene_name"])
    assert set(result.columns) == {"gene_id", "gene_name"}


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
    result = hsapiens108.tables_for_columns(["gene_id", "transcript_id"])
    assert result == ["entrezgene"]
    result = hsapiens108.tables_for_columns(["gene_name", "invalid_column"])
    assert result == ["gene"]
    with pytest.raises(ValueError):
        hsapiens108.tables_for_columns(["invalid_column"])
    with pytest.raises(ValueError):
        hsapiens108.tables_for_columns([])
