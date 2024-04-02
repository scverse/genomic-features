# API

## Ensembl

```{eval-rst}
.. module:: genomic_features.ensembl
.. currentmodule:: genomic_features

.. autosummary::
    :toctree: generated

    ensembl.annotation
    ensembl.EnsemblDB
    ensembl.list_ensdb_annotations
```

## Filters

```{eval-rst}
.. module:: genomic_features.filters
.. currentmodule:: genomic_features

.. autosummary::
    :toctree: generated

    filters.GeneIDFilter
    filters.GeneBioTypeFilter
    filters.GeneNameFilter
    filters.SeqNameFilter
    filters.GeneRangesFilter
    filters.TxIDFilter
    filters.TxBioTypeFilter
    filters.ExonIDFilter
    filters.CanonicalTxFilter
    filters.UniProtIDFilter
    filters.UniProtDBFilter
    filters.UniProtMappingTypeFilter
    filters.EmptyFilter
```

### Filter base classes (don't use these directly)

```{eval-rst}
.. module:: genomic_features._core.filters
.. currentmodule:: genomic_features

.. autosummary::
    :toctree: generated

    _core.filters.AbstractFilterExpr
    _core.filters.AbstractFilterOperatorExpr
    _core.filters.AndFilterExpr
    _core.filters.OrFilterExpr
    _core.filters.NotFilterExpr
    _core.filters.AbstractFilterEqualityExpr
```
