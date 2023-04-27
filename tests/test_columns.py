import pytest

import genomic_features as gf
from genomic_features._core import columns


@pytest.fixture(scope="module")
def hsapiens108():
    return gf.ensembl.annotation("Hsapiens", 108)


def test_listColumns(hsapiens108):
    result = columns.listColumns(hsapiens108, "gene")
    assert len(result) == len(hsapiens108.genes().columns)
    assert hsapiens108.genes().columns.isin(result).all()


def test_cleanColumns(hsapiens108):
    result = columns.cleanColumns(hsapiens108, ["gene_id", "gene_name"])
    assert result == ["gene_id", "gene_name"]
    result = columns.cleanColumns(hsapiens108, ["gene_id", "invalid_column"])
    assert result == ["gene_id"]
    with pytest.raises(ValueError):
        columns.cleanColumns(hsapiens108, ["invalid_column"])
    with pytest.raises(ValueError):
        columns.cleanColumns(hsapiens108, [])


def test_tablesForColumns(hsapiens108):
    result = columns.tablesForColumns(hsapiens108, ["gene_id", "gene_name"])
    assert result == ["gene"]
    result = columns.tablesForColumns(hsapiens108, ["gene_id", "transcript_id"])
    assert result == ["entrezgene", "gene", "tx"]
    result = columns.tablesForColumns(hsapiens108, ["gene_name", "invalid_column"])
    assert result == ["gene"]
    with pytest.raises(ValueError):
        columns.tablesForColumns(hsapiens108, ["invalid_column"])
    with pytest.raises(ValueError):
        columns.tablesForColumns(hsapiens108, [])
