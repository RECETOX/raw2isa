"""Compatibility and infrastructure helpers."""

from importlib.resources import files as resource_files

__all__ = ["resource_files", "tqdm"]

try:
    import tqdm
except ImportError:
    tqdm = None
