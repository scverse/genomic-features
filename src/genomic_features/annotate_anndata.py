import warnings

import pandas as pd


def annotate_anndata(
    adata_var: pd.DataFrame,
    annotation_df: pd.DataFrame,
    on: str = None,
    id_column: str = None,
) -> pd.DataFrame:
    """
    Annotate variables from AnnData object with a DataFrame of annotations.

    Parameters
    ----------
    adata_var
        The AnnData.var DataFrame.
    annotation_df
        The DataFrame of annotations (e.g. output of `EnsemblDB.genes()`).
    on
        The column name in `adata_var` to join on (useful if geneIDs are not stored in `adata.var_names`).
        If `None`, tables are joined based on index of `adata_var`.
    id_column
        The column name in `annotation_df` to join on. If `None`, the column name is inferred from the
        `annotation_df` DataFrame (i.e. the column with unique values ending in `_id`).

    Returns
    -------
    annotated_var
        The annotated AnnData.var DataFrame.
    """
    # Pick column with IDs in annotation table
    if id_column is None:
        id_column = [
            i
            for i in annotation_df.columns
            if annotation_df[i].is_unique and i.endswith("_id")
        ]
        if len(id_column) == 0:
            raise ValueError(
                "No unique ID column found in annotation_df - specify ID column with `id_column` parameter"
            )
        if len(id_column) > 1:
            raise ValueError(
                "Multiple unique ID columns found in annotation_df - specify ID column with `id_column` parameter"
            )
        else:
            id_column = id_column[0]
    else:
        if not annotation_df[id_column].is_unique:
            raise ValueError(f"Column {id_column} does not contain unique IDs")

    adata_var_merge = adata_var.copy()

    # Check for common columns between adata_var and annotation_df
    # (these will be overwritten by the merge)
    common_cols = annotation_df.columns[
        annotation_df.columns.isin(adata_var_merge.columns)
    ].tolist()
    if len(common_cols) > 0:
        warnings.warn(
            f"Columns {common_cols} are present in both adata_var and annotation_df - these will be overwritten",
            stacklevel=2,
        )
        adata_var_merge.drop(common_cols, axis=1, inplace=True)

    # Pick IDs in adata_var table
    if on is None:
        assert (
            adata_var_merge.index.is_unique
        ), "var_names do not contain unique IDs - specify ID column with `on` parameter"
        index_name = adata_var_merge.index.name
        adata_var_merge = adata_var_merge.reset_index(names="var_names")
    else:
        assert (
            on in adata_var_merge.columns
        ), "Column specified by `on` does not exist in adata_var_merge"
        assert adata_var_merge[
            on
        ].is_unique, "Column specified by `on` does not contain unique IDs"
        adata_var_merge["var_names"] = adata_var_merge[on].copy()

    # Merge
    annotated_var = (
        adata_var_merge.reset_index()
        .merge(annotation_df, how="left", right_on=[id_column], left_on=["var_names"])
        .set_index("index")
    )

    # Retain order and var_names
    if on is None:
        annotated_var = annotated_var.set_index("var_names")
        annotated_var = annotated_var.loc[adata_var_merge["var_names"]]
        annotated_var.index.name = index_name
    else:
        annotated_var = annotated_var.loc[adata_var_merge.index]
        annotated_var = annotated_var.drop("var_names", axis=1)

    # Check for missing genes
    missing_vars = annotated_var.index[
        annotated_var[annotation_df.drop([id_column] + common_cols, axis=1).columns]
        .isna()
        .all(axis=1)
    ]
    if len(missing_vars) > 0:
        warnings.warn(
            f"Missing annotations for {len(missing_vars)}/{annotated_var.shape[0]} vars",
            stacklevel=2,
        )

    return annotated_var
