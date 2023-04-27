from importlib.metadata import version

from . import ensembl, filters

__all__ = ["ensembl"]

__version__ = version("genomic-features")
