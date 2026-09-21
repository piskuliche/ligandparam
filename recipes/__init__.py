"""Public ligand parameterization recipes.

Exports are loaded lazily so importing ``ligandparam.recipes.DihedOptions``
(or ``registry``) does not pull every recipe and its stage graph.
"""

from __future__ import annotations

from typing import Any

_EXPORTS = {
    "LazyLigand": ".LazyLigand",
    "LazierLigand": ".LazierLigand",
    "FreeLigand": ".FreeLigand",
    "DPLigand": ".DpLazyLigand",
    "DPFreeLigand": ".DpFreeLigand",
    "SQMLigand": ".OptLigand",
    "available_recipes": ".Registry",
    "get_recipe": ".Registry",
}

__all__ = list(_EXPORTS)


def __getattr__(name: str) -> Any:
    mod = _EXPORTS.get(name)
    if mod is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    import importlib

    m = importlib.import_module(mod, __name__)
    return getattr(m, name)


def __dir__() -> list[str]:
    return sorted(set(__all__) | set(globals()))
