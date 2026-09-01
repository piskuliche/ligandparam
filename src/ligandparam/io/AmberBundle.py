"""Shared Amber ligand file-bundle helpers for ligandparam CLIs and stages.

After ``lig-getparam``, recipes typically write ``{stem}.mol2``, ``{stem}.lib``,
and ``{stem}.frcmod`` under ``{cwd}/{data_cwd}/{resname}/``. Dihedral correction
and scission CLIs both need to resolve that layout; keep the logic here.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


_ANCILLARY_MOL2_MARKERS = (
    ".initial.",
    ".centered.",
    ".resp.",
)


def _is_recipe_output_mol2(path: Path) -> bool:
    """Return True when ``path`` looks like a final recipe mol2 (not intermediate)."""
    name = path.name
    if any(marker in name for marker in _ANCILLARY_MOL2_MARKERS):
        return False
    if name.startswith("final_"):
        return False
    return True


@dataclass(frozen=True)
class AmberLigandBundle:
    """Parent Amber ligand triplet plus its working directory.

    Attributes
    ----------
    mol2, lib, frcmod
        Absolute paths to the charged structure, Leap library, and frcmod.
    work_dir
        Directory containing those files (typically ``{data_cwd}/{resname}``).
    """

    mol2: Path
    lib: Path
    frcmod: Path
    work_dir: Path

    @property
    def stem(self) -> str:
        """File stem shared by the mol2 / lib / frcmod trio."""
        return self.mol2.stem

    def as_input_paths(self, ligand_name: str | None = None) -> dict[str, Path | str | None]:
        """Return the Amber triplet as a mapping (ALPS builds scission InputBundle)."""
        return {
            "mol2_path": self.mol2,
            "lib_path": self.lib,
            "frcmod_path": self.frcmod,
            "ligand_name": ligand_name,
        }

    def to_scission_input(self, ligand_name: str | None = None) -> dict[str, Path | str | None]:
        """Alias of :meth:`as_input_paths` (does not import scission)."""
        return self.as_input_paths(ligand_name=ligand_name)


def resolve_getparam_bundle(
    *,
    cwd: Path | None = None,
    data_cwd: Path | str | None = None,
    resname: str | None = None,
    label: str | None = None,
    mol2: Path | str | None = None,
    lib: Path | str | None = None,
    frcmod: Path | str | None = None,
) -> AmberLigandBundle:
    """Resolve a ligandparam Amber triplet from explicit paths or getparam layout.

    Parameters
    ----------
    cwd
        Base directory for ``data_cwd`` (default: process CWD).
    data_cwd, resname
        Same ``-d`` / ``-r`` values used with ``lig-getparam``. Required unless
        all of ``mol2``, ``lib``, and ``frcmod`` are provided.
    label
        Recipe file stem (e.g. ``chaps`` from ``chaps.mol2``). Defaults to
        ``resname``, then to a unique non-ancillary ``*.mol2`` in ``work_dir``.
    mol2, lib, frcmod
        Explicit paths. If all three are set, ``work_dir`` is the mol2 parent
        and ``data_cwd`` / ``resname`` are ignored.

    Returns
    -------
    AmberLigandBundle

    Raises
    ------
    ValueError
        If neither an explicit triplet nor ``data_cwd``+``resname`` is given.
    FileNotFoundError
        If the working directory or expected output files are missing.
    """
    if mol2 is not None and lib is not None and frcmod is not None:
        mol2_p = Path(mol2).resolve()
        lib_p = Path(lib).resolve()
        frcmod_p = Path(frcmod).resolve()
        return AmberLigandBundle(
            mol2=mol2_p,
            lib=lib_p,
            frcmod=frcmod_p,
            work_dir=mol2_p.parent,
        )

    if data_cwd is None or resname is None:
        raise ValueError(
            "Provide either (mol2, lib, frcmod) or (data_cwd and resname)."
        )

    base = Path.cwd() if cwd is None else Path(cwd)
    work_dir = (base / Path(data_cwd) / resname).resolve()
    if not work_dir.is_dir():
        raise FileNotFoundError(f"Working directory does not exist: {work_dir}")

    stem = label or resname
    cand_mol2 = work_dir / f"{stem}.mol2"
    cand_lib = work_dir / f"{stem}.lib"
    cand_frcmod = work_dir / f"{stem}.frcmod"

    if not cand_mol2.is_file() and label is None:
        mol2s = sorted(p for p in work_dir.glob("*.mol2") if _is_recipe_output_mol2(p))
        if len(mol2s) == 1:
            cand_mol2 = mol2s[0]
            stem = cand_mol2.stem
            cand_lib = work_dir / f"{stem}.lib"
            cand_frcmod = work_dir / f"{stem}.frcmod"

    missing = [p for p in (cand_mol2, cand_lib, cand_frcmod) if not p.is_file()]
    if missing:
        raise FileNotFoundError(
            "Could not find ligandparam outputs in "
            f"{work_dir}. Missing: {', '.join(p.name for p in missing)}. "
            "Pass --label (input stem used by the recipe) or explicit "
            "mol2/lib/frcmod paths."
        )

    return AmberLigandBundle(
        mol2=cand_mol2.resolve(),
        lib=cand_lib.resolve(),
        frcmod=cand_frcmod.resolve(),
        work_dir=work_dir,
    )
