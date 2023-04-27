import pytest

import genomic_features as gf


@pytest.fixture(scope="module")
def hsapiens108():
    return gf.ensembl.annotation("Hsapiens", 108)


def test_columns_subset(hsapiens108):
    result = hsapiens108.genes(cols=["gene_id", "gene_name"])
    assert set(result.columns) == {"gene_id", "gene_name"}


def test_listColumns(hsapiens108):
    result = hsapiens108.listColumns("gene")
    assert len(result) == len(hsapiens108.genes().columns)
    assert hsapiens108.genes().columns.isin(result).all()


def test_cleanColumns(hsapiens108):
    result = hsapiens108.cleanColumns(["gene_id", "gene_name"])
    assert result == ["gene_id", "gene_name"]
    result = hsapiens108.cleanColumns(["gene_id", "invalid_column"])
    assert result == ["gene_id"]
    with pytest.raises(ValueError):
        hsapiens108.cleanColumns(["invalid_column"])
    with pytest.raises(ValueError):
        hsapiens108.cleanColumns([])


def test_tablesForColumns(hsapiens108):
    result = hsapiens108.tablesForColumns(["gene_id", "gene_name"])
    assert result == ["gene"]
    result = hsapiens108.tablesForColumns(["gene_id", "transcript_id"])
    assert result == ["entrezgene", "gene", "tx"]
    result = hsapiens108.tablesForColumns(["gene_name", "invalid_column"])
    assert result == ["gene"]
    with pytest.raises(ValueError):
        hsapiens108.tablesForColumns(["invalid_column"])
    with pytest.raises(ValueError):
        hsapiens108.tablesForColumns([])
