"""Parametrization drivers and shared recipe option helpers.

Classes
-------
Parametrization
    Base class that owns input paths, logging, and stage lists.
Recipe
    Thin subclass used by concrete ligand workflows.

Functions
---------
fresh_recipe_defaults
    Build a new defaults mapping with fresh mutable values.
apply_option_defaults
    Assign kwargs or defaults onto a recipe instance safely.
configure_gaussian_recipe
    Shared ``__init__`` configuration for Gaussian-based recipes.
"""

import logging
from pathlib import Path
from typing import Any, Iterable, Mapping, MutableMapping, Optional, Union
from typing_extensions import override

from ligandparam.Driver import Driver
from ligandparam.Log import get_logger, log_success_quote, set_file_logger, set_stream_logger


def fresh_recipe_defaults() -> dict[str, Any]:
    """Return recipe option defaults with newly created mutable values.

    Call this (or rely on :func:`apply_option_defaults`) whenever defaults are
    needed so ``theory`` / ``leaprc`` are never shared across instances.
    """
    return {
        "theory": {"low": "HF/6-31G*", "high": "PBE1PBE/6-31G*"},
        "leaprc": ["leaprc.gaff2"],
        "force_gaussian_rerun": False,
        "nproc": 1,
        "mem": 1,
    }


def _copy_if_mutable(value: Any) -> Any:
    """Shallow-copy dict/list values; leave immutables unchanged."""
    if isinstance(value, dict):
        return dict(value)
    if isinstance(value, list):
        return list(value)
    return value


def apply_option_defaults(
    obj: Any,
    kwargs: MutableMapping[str, Any],
    option_names: Iterable[str],
    *,
    defaults: Optional[Mapping[str, Any]] = None,
) -> None:
    """Set ``obj`` attributes from ``kwargs``, else from fresh defaults.

    User-supplied dict/list values are shallow-copied so later mutation on the
    instance (e.g. ``add_leaprc``) does not alter the caller's objects.
    """
    if defaults is None:
        defaults = fresh_recipe_defaults()
    for opt in option_names:
        if opt in kwargs:
            setattr(obj, opt, _copy_if_mutable(kwargs.pop(opt)))
        else:
            setattr(obj, opt, _copy_if_mutable(defaults[opt]))


def configure_gaussian_recipe(
    obj: Any,
    kwargs: MutableMapping[str, Any],
    *,
    with_orientation: bool = False,
    with_dihed: bool = False,
) -> None:
    """Apply shared Gaussian-recipe option parsing onto ``obj``.

    Pops ``logger`` (stages receive the recipe logger explicitly), requires
    ``net_charge``, applies theory/leaprc/Gaussian-rerun/nproc/mem defaults,
    and stores optional Gaussian path overrides. Optionally configures
    multi-RESP orientation protocol and/or dihedral-correction options.
    Remaining kwargs are assigned to ``obj.kwargs``.
    """
    kwargs.pop("logger", None)

    try:
        obj.net_charge = kwargs.pop("net_charge")
    except KeyError as exc:
        raise KeyError("Missing net_charge") from exc

    apply_option_defaults(
        obj,
        kwargs,
        ("theory", "leaprc", "force_gaussian_rerun", "nproc", "mem"),
    )

    for opt in ("gaussian_root", "gauss_exedir", "gaussian_binary", "gaussian_scratch"):
        setattr(obj, opt, kwargs.pop(opt, None))

    if with_orientation:
        from ligandparam.io.Orientations import DEFAULT_ORIENTATION_PROTOCOL

        obj.orientation_protocol = kwargs.pop(
            "orientation_protocol", DEFAULT_ORIENTATION_PROTOCOL
        )
        if obj.orientation_protocol not in ("so3_n28", "legacy_euler"):
            raise ValueError(
                "orientation_protocol must be 'so3_n28' or 'legacy_euler', "
                f"got {obj.orientation_protocol!r}"
            )
        if obj.orientation_protocol == "so3_n28":
            for key in ("alpha", "beta", "gamma"):
                kwargs.pop(key, None)

    if with_dihed:
        from ligandparam.recipes.DihedOptions import apply_dihed_options

        apply_dihed_options(obj, kwargs)

    obj.kwargs = kwargs


class Parametrization(Driver):
    """Base ligand parameterization workflow.

    Subclasses (or callers) populate ``self.stages`` and then call
    :meth:`~ligandparam.Driver.Driver.execute`. See :class:`Recipe` for the
    thin alias used by the built-in recipe modules.
    """

    @override
    def __init__(self, in_filename: Union[Path, str], cwd: Union[Path, str], *args, **kwargs):
        """
        Parameters
        ----------
        in_filename : path-like
            Input ligand structure.
        cwd : path-like
            Working directory for intermediate and output files.
        label : str, optional
            Ligand label; defaults to the stem of ``in_filename``.
        leaprc : list of str, optional
            Leaprc files; default ``["leaprc.gaff2"]``.
        logger : {"file", "stream"} or logging.Logger, optional
            Logging destination.

        Raises
        ------
        ValueError
            If ``logger`` is an unrecognized string or type.
        """
        self.in_filename = Path(in_filename).resolve()
        self.label = kwargs.get("label", self.in_filename.stem)
        self.cwd = Path(cwd)
        self.stages = []
        if "leaprc" in kwargs:
            self.leaprc = list(kwargs["leaprc"])
        else:
            self.leaprc = ["leaprc.gaff2"]
        try:
            logger = kwargs.pop("logger")
            if isinstance(logger, str):
                if logger == "file":
                    self.logger = set_file_logger(self.cwd / f"{self.label}.log")
                elif logger == "stream":
                    self.logger = set_stream_logger()
                else:
                    raise ValueError("Invalid input string for logger. Must be either 'file' or 'stream'.")
            elif isinstance(logger, logging.Logger):
                self.logger = logger
            else:
                raise ValueError("logger must be a string or a logging.Logger instance.")
        except KeyError:
            self.logger = get_logger()

    def add_leaprc(self, leaprc) -> None:
        """
        Add a leaprc file to the list of leaprc files.

        Parameters
        ----------
        leaprc : str
            The name of the leaprc file to add.
        """
        self.leaprc.append(leaprc)


class Recipe(Parametrization):
    """Convenience alias for a ligand-parametrization workflow.

    Subclasses typically implement :meth:`setup` to populate ``self.stages``.
    :meth:`execute` logs start/done around the base stage pipeline.
    """

    @override
    def execute(
        self, dry_run=False, nproc: Optional[int] = None, mem: Optional[int] = None
    ) -> Any:
        """Run all stages defined by :meth:`setup`.

        Parameters
        ----------
        dry_run : bool, optional
            If True, log planned commands without running external programs.
        nproc : int, optional
            Override the recipe processor count for this run.
        mem : int, optional
            Override the recipe memory allocation in GB for this run.
        """
        name = type(self).__name__
        self.logger.info(f"Starting the {name} recipe at {self.cwd}")
        super().execute(dry_run=dry_run, nproc=nproc, mem=mem)
        self.logger.info(f"Done with the {name} recipe")
        log_success_quote(self.logger)
