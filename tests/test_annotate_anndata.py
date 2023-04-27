import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

import genomic_features as gf

genes = gf.ensembl.annotation("Hsapiens", 108).genes()

# Make dummy anndata object
adata = AnnData(
    X=np.array([[1, 2, 3], [4, 5, 6]]),
    obs=pd.DataFrame(index=["cell1", "cell2"]),
    var=pd.DataFrame(index=["ENSG00000278862", "ENSG00000212694", "ENSG00000113721"]),
)


def test_output():
    adata_var = gf.annotate_anndata(adata.var, genes)
    assert isinstance(adata_var, pd.DataFrame)


def test_size():
    annotated_adata_var = gf.annotate_anndata(adata.var, genes)
    assert annotated_adata_var.shape[0] == adata.n_vars


def test_var_names():
    annotated_adata_var = gf.annotate_anndata(adata.var, genes)
    assert "var_names" not in annotated_adata_var
    assert annotated_adata_var.index.equals(adata.var.index)
    assert annotated_adata_var.index.name == adata.var.index.name


def test_anndata_changes():
    gf.annotate_anndata(adata.var, genes)
    assert "gene_name" not in adata.var
    adata.var["gene_ids"] = adata.var_names.copy()
    gf.annotate_anndata(adata.var, genes, on="gene_ids")
    assert "gene_name" not in adata.var


# Test missing genes get NAs
adata_custom_genes = AnnData(
    X=np.array([[1, 2, 3, 5], [4, 5, 6, 8]]),
    obs=pd.DataFrame(index=["cell1", "cell2"]),
    var=pd.DataFrame(
        index=["ENSG00000278862", "ENSG00000212694", "ENSG00000113721", "GENE1"]
    ),
)


def test_missing_genes():
    with pytest.warns(UserWarning):
        annotated_adata_var = gf.annotate_anndata(adata_custom_genes.var, genes)
        assert annotated_adata_var.loc["GENE1"][genes.columns].isna().all()


# Test duplicate IDs
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


def test_unique_ids():
    with pytest.raises(AssertionError):
        gf.annotate_anndata(adata_duplicate_genes.var, genes)
    with pytest.raises(AssertionError):
        gf.annotate_anndata(adata.var, genes, id_column="gene_name")
