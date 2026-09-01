"""Recipe name -> constructor registry for CLI / drivers."""

from __future__ import annotations

from typing import Any, Callable

RecipeFactory = Callable[..., Any]

_REGISTRY: dict[str, str] = {
    "lazyligand": "ligandparam.recipes.LazyLigand:LazyLigand",
    "lazierligand": "ligandparam.recipes.LazierLigand:LazierLigand",
    "freeligand": "ligandparam.recipes.FreeLigand:FreeLigand",
    "dplazyligand": "ligandparam.recipes.DpLazyLigand:DPLigand",
    "dpfreeligand": "ligandparam.recipes.DpFreeLigand:DPFreeLigand",
    "sqmligand": "ligandparam.recipes.OptLigand:SQMLigand",
}


def available_recipes() -> list[str]:
    return sorted(_REGISTRY)


def get_recipe(name: str, **kwargs):
    """Instantiate a recipe by CLI name (case-insensitive)."""
    key = str(name).strip().lower()
    if key not in _REGISTRY:
        raise ValueError(
            f"Unknown recipe name: {name}. Available recipes: {', '.join(available_recipes())}"
        )
    path, cls_name = _REGISTRY[key].split(":")
    import importlib

    mod = importlib.import_module(path)
    cls = getattr(mod, cls_name)
    return cls(**kwargs)
