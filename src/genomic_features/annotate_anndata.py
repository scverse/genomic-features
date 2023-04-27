import warnings

from pandas import DataFrame


def annotate_anndata(
    adata_var: DataFrame,
    annotation_df: DataFrame,
    on: str = None,
    id_column: str = None,
) -> DataFrame:
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
                "No unique ID column found in annotation_df - specify ID column with `on` parameter"
            )
    else:
        assert annotation_df[
            id_column
        ].is_unique, "Column specified by `id_column` does not contain unique IDs"

    # Pick IDs in adata_var table
    if on is None:
        assert (
            adata_var.index.is_unique
        ), "var_names do not contain unique IDs - specify ID column with `on` parameter"
        index_name = adata_var.index.name
        adata_var = adata_var.reset_index(names="var_names")
    else:
        assert (
            on in adata_var.columns
        ), "Column specified by `on` does not exist in adata_var"
        assert adata_var[
            on
        ].is_unique, "Column specified by `on` does not contain unique IDs"
        adata_var["var_names"] = adata_var[on].copy()

    annotated_var = adata_var.merge(
        annotation_df, how="left", right_on=id_column, left_on="var_names"
    )

    # Retain order and var_names
    if on is None:
        annotated_var = annotated_var.set_index("var_names")
        annotated_var = annotated_var.loc[adata_var["var_names"]]
        annotated_var.index.name = index_name
    else:
        annotated_var = annotated_var.loc[adata_var.index]

    # Check for missing genes
    missing_vars = annotated_var.index[
        annotated_var[annotation_df.columns].isna().all(axis=1)
    ]
    if len(missing_vars) > 0:
        warnings.warn(
            f"Missing annotations for vars {missing_vars.tolist()}", stacklevel=2
        )

    return annotated_var
