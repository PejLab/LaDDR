"""LaDDR package."""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("laddr")
except PackageNotFoundError:  # pragma: no cover - used only in editable/source trees
    __version__ = "unknown"
