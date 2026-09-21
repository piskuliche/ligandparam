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
    ".minimized.",
    ".user_input.",
    ".antechamber_in.",
    ".sanitized.",
    ".gau_geom.",
)


def _is_recipe_output_mol2(path: Path) -> bool:
    """Return True when ``path`` looks like a final recipe mol2 (not intermediate)."""
    name = path.name
    if any(marker in name for marker in _ANCILLARY_MOL2_MARKERS):
        return False
    if name.startswith("final_"):
        return False
    return True


def _find_file_ci(directory: Path, filename: str) -> Path | None:
    """Return ``directory/filename``, matching the name case-insensitively."""
    wanted = filename.lower()
    matches = [
        p for p in directory.iterdir() if p.is_file() and p.name.lower() == wanted
    ]
    if not matches:
        return None
    exact = [p for p in matches if p.name == filename]
    return (exact or matches)[0]


def _triplet_for_stem(work_dir: Path, stem: str) -> tuple[Path, Path, Path] | None:
    """Return mol2/lib/frcmod for ``stem`` if all three files exist."""
    mol2 = _find_file_ci(work_dir, f"{stem}.mol2")
    lib = _find_file_ci(work_dir, f"{stem}.lib")
    frcmod = _find_file_ci(work_dir, f"{stem}.frcmod")
    if mol2 is not None and lib is not None and frcmod is not None:
        return mol2, lib, frcmod
    return None


def _complete_triplets(work_dir: Path) -> list[tuple[Path, Path, Path]]:
    """Amber mol2/lib/frcmod triples in ``work_dir`` (final recipe mol2 only)."""
    found: list[tuple[Path, Path, Path]] = []
    seen: set[str] = set()
    for path in work_dir.iterdir():
        if not path.is_file() or path.suffix.lower() != ".mol2":
            continue
        if not _is_recipe_output_mol2(path):
            continue
        trio = _triplet_for_stem(work_dir, path.stem)
        if trio is None:
            continue
        key = trio[0].name.lower()
        if key in seen:
            continue
        seen.add(key)
        found.append(trio)
    return sorted(found, key=lambda t: t[0].name.lower())


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
        Recipe file stem (e.g. ``SDS`` from ``SDS.mol2``). Case-insensitive.
        Defaults to ``resname``, then to a unique Amber triplet in ``work_dir``.
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

    stems: list[str] = []
    for stem in (label, resname):
        if stem and stem not in stems:
            stems.append(stem)

    trio = None
    for stem in stems:
        trio = _triplet_for_stem(work_dir, stem)
        if trio is not None:
            break
    if trio is None:
        available = _complete_triplets(work_dir)
        if len(available) == 1:
            trio = available[0]

    if trio is None:
        available = _complete_triplets(work_dir)
        found = (
            ", ".join(f"{m.name} + {lb.name} + {fr.name}" for m, lb, fr in available)
            if available
            else "none"
        )
        looked = ", ".join(f"{s}.mol2/{s}.lib/{s}.frcmod" for s in stems) or "n/a"
        raise FileNotFoundError(
            "Could not find ligandparam outputs in "
            f"{work_dir}. Looked for {looked} (case-insensitive). "
            f"Found Amber triplets: {found}. "
            "Pass --label matching the recipe input stem "
            "(SDS.mol2 -> SDS, not sds) or explicit mol2/lib/frcmod paths."
        )

    cand_mol2, cand_lib, cand_frcmod = trio
    return AmberLigandBundle(
        mol2=cand_mol2.resolve(),
        lib=cand_lib.resolve(),
        frcmod=cand_frcmod.resolve(),
        work_dir=work_dir,
    )
