from importlib.metadata import version

from . import ensembl, filters
from .annotate import annotate

__all__ = ["ensembl"]

__version__ = version("genomic-features")
