import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

import genomic_features as gf

@pytest.fixture(scope="module")
def genes():
    return gf.ensembl.annotation("Hsapiens", 108).genes()

@pytest.fixture(scope="module")
def adata():
    return AnnData(
        X=np.array([[1, 2, 3], [4, 5, 6]]),
        obs=pd.DataFrame(index=["cell1", "cell2"]),
        var=pd.DataFrame(index=["ENSG00000278862", "ENSG00000212694", "ENSG00000113721"]),
    )

@pytest.fixture(scope="module")
def annotated_adata_var(adata, genes):
    return gf.annotate_anndata(adata.var, genes)

def test_output(annotated_adata_var):
    assert isinstance(annotated_adata_var, pd.DataFrame)

def test_size(annotated_adata_var, adata):
    assert annotated_adata_var.shape[0] == adata.n_vars

def test_var_names(annotated_adata_var, adata):
    assert "var_names" not in annotated_adata_var
    assert annotated_adata_var.index.equals(adata.var.index)
    assert annotated_adata_var.index.name == adata.var.index.name


def test_anndata_changes(adata, genes):
    gf.annotate_anndata(adata.var, genes)
    assert "gene_name" not in adata.var
    adata.var["gene_ids"] = adata.var_names.copy()
    gf.annotate_anndata(adata.var, genes, on="gene_ids")
    assert "gene_name" not in adata.var


def test_missing_genes(genes):
    # Test missing genes get NAs
    adata_custom_genes = AnnData(
        X=np.array([[1, 2, 3, 5], [4, 5, 6, 8]]),
        obs=pd.DataFrame(index=["cell1", "cell2"]),
        var=pd.DataFrame(
            index=["ENSG00000278862", "ENSG00000212694", "ENSG00000113721", "GENE1"]
        ),
    )
    with pytest.warns(UserWarning):
        annotated_adata_var = gf.annotate_anndata(adata_custom_genes.var, genes)
        assert annotated_adata_var.loc["GENE1"][genes.columns].isna().all()

def test_unique_ids(genes, adata):
    adata_duplicate_genes = AnnData(
        X=np.array([[1, 2, 3, 5], [4, 5, 6, 8]]),
        obs=pd.DataFrame(index=["cell1", "cell2"]),
        var=pd.DataFrame(
            index=[
                "ENSG00000278862",
                "ENSG00000212694",
                "ENSG00000113721",
                "ENSG00000212694",
            ]
        ),
    )
    with pytest.raises(ValueError):
        gf.annotate_anndata(adata_duplicate_genes.var, genes)
    with pytest.raises(ValueError):
        gf.annotate_anndata(adata.var, genes, id_column="gene_name")


def test_clashing_columns(adata, genes):
    # Test that common columns are dropped
    adata.var["gene_biotype"] = "foo"
    adata.var["seq_name"] = "chr1"
    with pytest.warns(UserWarning, match=r"gene_biotype.+seq_name"):
        annotated_var = gf.annotate_anndata(adata.var, genes)
    assert 'lncRNA' in annotated_var["gene_biotype"].values
    assert '12' in annotated_var["seq_name"].values
    
