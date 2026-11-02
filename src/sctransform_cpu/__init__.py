"""CPU reference implementation of the validated SCT v2 paths."""
from importlib.metadata import PackageNotFoundError, version

from .core import SCTConfig, SCTResult, sctransform
from .sampling import SelectionManifest, make_selection_manifest

try:
    __version__ = version("sctransform-cpu")
except PackageNotFoundError:  # direct source-tree imports before installation
    __version__ = "0+unknown"

__all__ = [
    "SCTConfig", "SCTResult", "SelectionManifest", "make_selection_manifest",
    "__version__", "sctransform",
]
