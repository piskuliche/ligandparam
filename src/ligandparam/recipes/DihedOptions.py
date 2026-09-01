"""Opt-in kwargs for dihedral correction (owned by ALPS, not ligandparam).

Recipes still accept ``dihed_correct=...`` so existing callers do not
break. ligandparam does **not** import ffpopt or scission; ALPS runs
twist / fragment after parameterization.
"""

from __future__ import annotations

import logging
from typing import Any, MutableMapping


def coerce_fragment_config(value: Any):
    """Store a fragment-config mapping or pass-through object.

    ligandparam does not import scission. A dict is copied; other values
    (including a live ``FragmentConfig``) are stored as-is for ALPS.
    """
    if value is None:
        return None
    if isinstance(value, dict):
        return dict(value)
    return value


def pop_dihed_options(kwargs: MutableMapping[str, Any]) -> dict[str, Any]:
    """Pull dihedral-correction options out of recipe kwargs.

    These flags are recorded on the recipe for ALPS. ligandparam does not
    run ffpopt when they are set.
    """
    return {
        "dihed_correct": bool(kwargs.pop("dihed_correct", False)),
        "dihed_model": kwargs.pop("dihed_model", "qdpi2"),
        "dihed_maxiter": int(kwargs.pop("dihed_maxiter", 2)),
        "dihed_nprim": int(kwargs.pop("dihed_nprim", 3)),
        "dihed_delta": int(kwargs.pop("dihed_delta", 10)),
        "dihed_geometric_opt": bool(kwargs.pop("dihed_geometric_opt", True)),
        "dihed_skip_existing": bool(kwargs.pop("dihed_skip_existing", True)),
        "dihed_out_frcmod": kwargs.pop("dihed_out_frcmod", None),
        "dihed_out_dir": kwargs.pop("dihed_out_dir", None),
        "dihed_rotatable_bond_smarts": kwargs.pop("dihed_rotatable_bond_smarts", None),
        "dihed_fragment_config": coerce_fragment_config(
            kwargs.pop("dihed_fragment_config", None)
        ),
    }


def apply_dihed_options(obj: Any, kwargs: MutableMapping[str, Any]) -> None:
    """Set dihedral-correction attributes on a recipe from kwargs."""
    opts = pop_dihed_options(kwargs)
    for key, value in opts.items():
        setattr(obj, key, value)


def append_dihed_twist_stage(
    stages: list,
    *,
    recipe: Any,
    mol2,
    lib,
    frcmod,
) -> None:
    """No-op: ffpopt twist is an ALPS step, not a ligandparam stage.

    ``stages`` / path arguments are accepted so existing recipe ``setup()``
    calls stay valid.
    """
    _ = (stages, mol2, lib, frcmod)
    if not getattr(recipe, "dihed_correct", False):
        return
    msg = (
        "dihed_correct is handled by ALPS, not ligandparam. "
        "Skipping in-recipe twist. After lig-getparam, run lig-dihed-correct."
    )
    log = getattr(recipe, "logger", None)
    if log is not None:
        log.warning(msg)
    else:
        logging.getLogger("ligandparam").warning(msg)
