"""Multi-residue / multi-orientation RESP helpers (ParmEd + Amber MDIN)."""

from __future__ import annotations

from importlib import import_module
from typing import Any

_MODULES = {
    "EndState": ".EndState",
    "Functions": ".Functions",
    "IntermolEquiv": ".IntermolEquiv",
    "MdinUtils": ".MdinUtils",
    "ParmHelper": ".ParmHelper",
    "ResidueResp": ".ResidueResp",
    "RespFunctions": ".RespFunctions",
}
_ALIASES = {
    "endstate": "EndState",
    "functions": "Functions",
    "intermolequiv": "IntermolEquiv",
    "mdinutils": "MdinUtils",
    "parmhelper": "ParmHelper",
    "residueresp": "ResidueResp",
    "respfunctions": "RespFunctions",
}

__all__ = list(_MODULES) + list(_ALIASES)


def __getattr__(name: str) -> Any:
    canon = name if name in _MODULES else _ALIASES.get(name)
    if canon is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    module = import_module(_MODULES[canon], __name__)
    globals()[canon] = module
    if name != canon:
        globals()[name] = module
    return module


def __dir__() -> list[str]:
    return sorted(set(__all__) | set(globals()))
