"""Make ``import ligandparam`` work from this checkout (flat package dir)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path


def package_root() -> Path:
    """This ligandparam clone."""
    return Path(__file__).resolve().parents[1]


def ensure_ligandparam() -> None:
    """Bind the local tree before ``pip install -e .``."""
    if "ligandparam" in sys.modules:
        return
    pkg = package_root()
    init = pkg / "__init__.py"
    if not init.is_file():
        raise RuntimeError(f"ligandparam package init not found at {init}")
    spec = importlib.util.spec_from_file_location(
        "ligandparam",
        init,
        submodule_search_locations=[str(pkg)],
    )
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load ligandparam from {pkg}")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["ligandparam"] = mod
    spec.loader.exec_module(mod)
