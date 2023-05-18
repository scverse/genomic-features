from importlib.metadata import version

from .annotate_anndata import annotate_anndata
from . import ensembl, filters

__all__ = ["ensembl"]

__version__ = version("genomic-features")
