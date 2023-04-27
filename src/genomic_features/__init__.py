from importlib.metadata import version

from . import ensembl
from .annotate_anndata import annotate_anndata

__all__ = ["ensembl"]

__version__ = version("genomic-features")
