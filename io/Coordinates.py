import warnings
from collections import defaultdict
from dataclasses import dataclass
from typing import Optional, Union
import shutil
from pathlib import Path

import numpy as np

import MDAnalysis as mda
from MDAnalysis.topology.guessers import guess_atom_element, guess_masses

ATOMIC_NUMBERS = {
    "H": 1,
    "He": 2,
    "Li": 3,
    "Be": 4,
    "B": 5,
    "C": 6,
    "N": 7,
    "O": 8,
    "F": 9,
    "Ne": 10,
    "Na": 11,
    "Mg": 12,
    "Al": 13,
    "Si": 14,
    "P": 15,
    "S": 16,
    "Cl": 17,
    "K": 19,
    "Ca": 20,
    "Fe": 26,
    "Zn": 30,
    "Br": 35,
    "I": 53,
}


class Coordinates:
    """Thin MDAnalysis wrapper for reading and transforming structure coordinates."""

    def __init__(self, filename: Union[Path, str], filetype: str = 'pdb'):
        """
        Load a structure and sanitize masses for center-of-mass operations.

        Parameters
        ----------
        filename : Union[Path, str]
            Path to the structure file to read.
        filetype : str, optional
            File type hint (default: ``'pdb'``).
        """
        self.filename = Path(filename)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.u = mda.Universe(filename)
        self.original_coords = np.array(self.get_coordinates(), dtype=float, copy=True)

        # If the mol2 comes from antechaamber, then the atom names are weird and both rdkit and mda will have trouble
        if np.any(np.isclose(self.u.atoms.masses, 0, atol=0.1)):
            self.u.guess_TopologyAttrs(to_guess=['elements'], force_guess=['masses'])
        # We tried to get correct masses but may have failed in the process. Lack of masses will fail
        # MDAnalysis's center_of_mass(), so just set them to 1.0, since the exact values are not important
        self.u.atoms.masses[np.isclose(self.u.atoms.masses, 0, atol=0.1)] = 1.0

        return

    def get_coordinates(self):
        """Return the current atomic coordinates.

        Returns
        -------
        np.ndarray
            Coordinates of the atoms in the structure.
        """
        return self.u.atoms.positions

    def get_elements(self):
        """Return atomic element symbols.

        Returns
        -------
        list
            Element symbols for each atom.
        """
        try:
            return [atom.element for atom in self.u.atoms]
        except (AttributeError, TypeError, ValueError) as exc:
            import warnings
            warnings.warn(
                f"Could not read atom.element ({type(exc).__name__}: {exc}); "
                "guessing elements from atom names",
                stacklevel=2,
            )
            return self._get_elements_from_topology()

    def _get_elements_from_topology(self):
        """Guess element symbols from atom names in the topology.

        Returns
        -------
        list
            Guessed element symbols for each atom.
        """
        from MDAnalysis.topology.guessers import guess_types
        elements = guess_types(self.u.atoms.names)
        return elements

    def update_coordinates(self, coords, original=False):
        """Replace the current atomic coordinates.

        Parameters
        ----------
        coords : np.ndarray
            New coordinates with the same shape as the current positions.
        original : bool, optional
            If True, also update the stored original coordinates used by
            :meth:`rotate`.
        """
        assert np.shape(coords) == np.shape(self.get_coordinates()), "Coordinate dimensions do not match"
        self.u.atoms.positions = coords
        if original:
            self.original_coords = coords
        return

    def rotate(self, alpha=0.0, beta=0.0, gamma=0.0):
        """Rotate coordinates about the center of mass using Euler angles.

        Rotations are applied in order alpha (x), beta (y), gamma (z), matching
        the previous MDAnalysis ``rotateby`` sequence. Angles are in degrees.

        Parameters
        ----------
        alpha : float
            Rotation about the x-axis (degrees).
        beta : float
            Rotation about the y-axis (degrees).
        gamma : float
            Rotation about the z-axis (degrees).

        Returns
        -------
        np.ndarray
            Rotated coordinates with shape ``(n_atoms, 3)``.
        """
        coords = np.asarray(self.original_coords, dtype=float)
        # COM from original geometry (masses already sanitized in __init__)
        masses = self.u.atoms.masses
        com = np.average(coords, axis=0, weights=masses)

        a, b, g = np.deg2rad([alpha, beta, gamma])
        ca, sa = np.cos(a), np.sin(a)
        cb, sb = np.cos(b), np.sin(b)
        cg, sg = np.cos(g), np.sin(g)

        # Intrinsic/extrinsic composition matching sequential Rx, Ry, Rz on positions
        rx = np.array([[1.0, 0.0, 0.0],
                       [0.0, ca, -sa],
                       [0.0, sa, ca]])
        ry = np.array([[cb, 0.0, sb],
                       [0.0, 1.0, 0.0],
                       [-sb, 0.0, cb]])
        rz = np.array([[cg, -sg, 0.0],
                       [sg, cg, 0.0],
                       [0.0, 0.0, 1.0]])
        rotation = rz @ ry @ rx

        rotated = (coords - com) @ rotation.T + com
        self.u.atoms.positions = rotated
        return rotated

    def rotate_matrix(self, rotation: np.ndarray) -> np.ndarray:
        """Rotate the original coordinates using an explicit rotation matrix.

        The rotation is applied about the mass-weighted center of mass. This
        is the path used by quaternion SO(3) orientation protocols.

        Parameters
        ----------
        rotation : np.ndarray, shape (3, 3)
            Proper orthogonal rotation matrix.

        Returns
        -------
        np.ndarray
            Rotated coordinates with shape ``(n_atoms, 3)``.

        Raises
        ------
        ValueError
            If ``rotation`` is not a valid proper rotation matrix.
        """
        rotation = np.asarray(rotation, dtype=float)
        if rotation.shape != (3, 3):
            raise ValueError(f"Expected rotation shape (3, 3), got {rotation.shape}")
        if not np.allclose(rotation.T @ rotation, np.eye(3), atol=1e-8):
            raise ValueError("Rotation matrix must be orthogonal")
        if not np.isclose(np.linalg.det(rotation), 1.0, atol=1e-8):
            raise ValueError("Rotation matrix must have determinant +1")

        coords = np.asarray(self.original_coords, dtype=float)
        com = np.average(coords, axis=0, weights=self.u.atoms.masses)
        rotated = (coords - com) @ rotation.T + com
        self.u.atoms.positions = rotated
        return rotated


def SimpleXYZ(file_obj, coordinates):
    """Write coordinates to a simple XYZ trajectory frame.

    Parameters
    ----------
    file_obj : file object
        Open file handle to write to.
    coordinates : np.ndarray
        Coordinates to write.
    """
    file_obj.write(f"{len(coordinates)}\n")
    file_obj.write("Generated by ligand_param\n")
    for i, coord in enumerate(coordinates):
        file_obj.write(f"{i + 1} {coord[0]} {coord[1]} {coord[2]}\n")
    return


class Mol2Writer:
    """Write an MDAnalysis Universe selection to a mol2 file."""

    def __init__(self, u, filename=None, selection="all"):
        """
        Parameters
        ----------
        u : MDAnalysis.Universe
            Universe to write.
        filename : str, optional
            Output mol2 path.
        selection : str, optional
            Atom selection string (default: ``'all'``).
        """
        self.u = u
        self.filename = Path(filename)
        self.selection = selection
        return

    def _write(self):
        """Write the selected atoms to mol2 via MDAnalysis."""
        ag = self.u.select_atoms(self.selection)
        ag.write(self.filename)

    def _remove_blank_lines(self):
        """Remove blank lines from the written mol2 file.

        Raises
        ------
        FileNotFoundError
            If the output file does not exist.
        """
        if Path(self.filename).exists():
            # Read the file and filter out blank lines
            with open(self.filename, 'r') as file:
                lines = file.readlines()
                non_blank_lines = [line for line in lines if line.strip()]

            # Write the non-blank lines back to the file
            with open(self.filename, 'w') as file:
                file.writelines(non_blank_lines)
        else:
            raise FileNotFoundError(f"File {self.filename} not found.")

    def write(self):
        """Write the mol2 file and strip blank lines that confuse antechamber."""
        self._write()
        self._remove_blank_lines()
        return


def _split_pdb_element_charge(token: str) -> tuple[str, str]:
    """Split Open Babel-style ``O1-`` / ``N+1`` into PDB element + charge.

    Parameters
    ----------
    token : str
        Trailing PDB element/charge field (columns 77-80 or a overflow token).

    Returns
    -------
    element : str
        One- or two-letter element symbol (``O``, ``Cl``, …).
    charge : str
        PDB charge field (``1-``, ``2+``, …) or empty.
    """
    import re

    t = (token or "").strip()
    if not t:
        return "", ""
    m = re.fullmatch(
        r"([A-Za-z]{1,2})(?:(\d+)([+-])|([+-])(\d+)|([+-]))?",
        t,
    )
    if not m:
        letters = "".join(c for c in t if c.isalpha())[:2]
        if not letters:
            return "", ""
        elem = letters[0].upper() + letters[1:].lower()
        return elem, ""
    elem = m.group(1)
    elem = elem[0].upper() + elem[1:].lower()
    if m.group(2) and m.group(3):
        charge = f"{m.group(2)}{m.group(3)}"
    elif m.group(4) and m.group(5):
        charge = f"{m.group(5)}{m.group(4)}"
    elif m.group(6):
        charge = f"1{m.group(6)}"
    else:
        charge = ""
    return elem, charge


def _rewrite_pdb_atom_line(line: str) -> str:
    """Rewrite ATOM/HETATM so cols 77-78 are element and 79-80 are charge."""
    nl = "\n" if line.endswith("\n") else ""
    raw = line.rstrip("\n")
    rec = raw[:6].strip()
    if rec not in ("ATOM", "HETATM"):
        return line
    if len(raw) < 80:
        raw = raw.ljust(80)
    elem, charge = _split_pdb_element_charge(raw[76:80])
    if not elem:
        # Open Babel often parks ``O1-`` after the B-factor instead of cols 77-80.
        elem, charge = _split_pdb_element_charge(raw[66:].strip())
    if not elem:
        atom_name = raw[12:16].strip()
        letters = "".join(c for c in atom_name if c.isalpha())[:2].upper()
        special = {"CL": "Cl", "BR": "Br", "NA": "Na", "MG": "Mg", "FE": "Fe", "ZN": "Zn"}
        if letters in special:
            elem = special[letters]
        elif letters:
            elem = letters[0]
    if not elem:
        elem = "C"
    return raw[:76] + f"{elem:>2}{charge:<2}" + nl


def _rewrite_pdb_conect_line(line: str) -> list[str]:
    """Keep CONECT, but drop duplicate partners (Open Babel S=O as ``15 15``)."""
    parts = line.split()
    if len(parts) < 2:
        return [line if line.endswith("\n") else line + "\n"]
    serial = int(parts[1])
    unique: list[int] = []
    for p in parts[2:]:
        try:
            idx = int(p)
        except ValueError:
            continue
        if idx not in unique:
            unique.append(idx)
    if not unique:
        return [f"CONECT{serial:5d}\n"]
    out = []
    for i in range(0, len(unique), 4):
        rec = f"CONECT{serial:5d}"
        for idx in unique[i : i + 4]:
            rec += f"{idx:5d}"
        out.append(rec + "\n")
    return out


def count_structure_atoms(path: Union[Path, str]) -> int:
    """Count atoms in a PDB (ATOM/HETATM) or mol2 (``@<TRIPOS>ATOM``) file."""
    p = Path(path)
    suffix = p.suffix.lower()
    if suffix == ".pdb":
        n = 0
        with open(p, encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line.startswith(("ATOM", "HETATM")):
                    n += 1
        return n
    if suffix == ".mol2":
        n = 0
        in_atom = False
        with open(p, encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line.startswith("@<TRIPOS>"):
                    in_atom = line.startswith("@<TRIPOS>ATOM")
                    continue
                if in_atom and line.strip():
                    parts = line.split()
                    if not parts:
                        continue
                    try:
                        int(parts[0])
                    except ValueError:
                        continue
                    n += 1
        return n
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return int(len(mda.Universe(str(p)).atoms))


def _element_from_mol2_name_type(name: str, atype: str) -> str:
    """Map mol2 atom name / GAFF type to a Gaussian element symbol."""
    for token in (atype, name):
        letters = "".join(c for c in (token or "") if c.isalpha())
        if not letters:
            continue
        two = letters[:2].capitalize()
        if two in {"Cl", "Br", "Na", "Mg", "Fe", "Zn", "Si"}:
            return two
        return letters[0].upper()
    return "C"


def parse_structure_atoms(path: Union[Path, str]) -> list[tuple[str, np.ndarray]]:
    """Read element symbols and coordinates from PDB or mol2 (no MDA).

    MDAnalysis element guessing can invent a hydrogen on anionic oxygen when
    building the Gaussian ``.com``. This parser only uses file records.
    """
    p = Path(path)
    suffix = p.suffix.lower()
    atoms: list[tuple[str, np.ndarray]] = []
    if suffix == ".pdb":
        with open(p, encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if not line.startswith(("ATOM", "HETATM")):
                    continue
                rewritten = _rewrite_pdb_atom_line(line)
                elem = rewritten[76:78].strip() or "C"
                xyz = np.array(
                    [float(rewritten[30:38]), float(rewritten[38:46]), float(rewritten[46:54])],
                    dtype=float,
                )
                atoms.append((elem, xyz))
        return atoms
    if suffix == ".mol2":
        in_atom = False
        with open(p, encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line.startswith("@<TRIPOS>"):
                    in_atom = line.startswith("@<TRIPOS>ATOM")
                    continue
                if not in_atom or not line.strip():
                    continue
                parts = line.split()
                if len(parts) < 6:
                    continue
                try:
                    int(parts[0])
                except ValueError:
                    continue
                name, xs, ys, zs, atype = parts[1], parts[2], parts[3], parts[4], parts[5]
                elem = _element_from_mol2_name_type(name, atype)
                atoms.append((elem, np.array([float(xs), float(ys), float(zs)], dtype=float)))
        return atoms
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        u = mda.Universe(str(p))
        elems = [str(a.element) for a in u.atoms]
        return list(zip(elems, np.asarray(u.atoms.positions, dtype=float)))


def match_current_to_reference(
    reference: list[tuple[str, np.ndarray]],
    current: list[tuple[str, np.ndarray]],
    max_dist: float = 0.75,
) -> tuple[list[str], np.ndarray]:
    """Keep current atoms that match a reference, dropping extras (added H).

    Aligns on the heavy-atom centroid so a spurious OH hydrogen does not
    shift the frame. Returned coordinates are the current (e.g. centered)
    positions in reference order.
    """
    if len(current) < len(reference):
        raise RuntimeError(
            f"Current structure has fewer atoms ({len(current)}) than the "
            f"reference ({len(reference)})."
        )
    ref_e = [str(e) for e, _ in reference]
    ref_x = np.asarray([x for _, x in reference], dtype=float)
    cur_e = [str(e) for e, _ in current]
    cur_x = np.asarray([x for _, x in current], dtype=float)

    def _heavy_centroid(elems, xyz):
        heavy = np.array(
            [row for el, row in zip(elems, xyz) if str(el).upper() != "H"],
            dtype=float,
        )
        if heavy.size == 0:
            return xyz.mean(axis=0)
        return heavy.mean(axis=0)

    ref_a = ref_x - _heavy_centroid(ref_e, ref_x)
    cur_a = cur_x - _heavy_centroid(cur_e, cur_x)
    used: set[int] = set()
    keep: list[int] = []
    for e, x in zip(ref_e, ref_a):
        best_j, best_d = None, 1e9
        for j, (ej, y) in enumerate(zip(cur_e, cur_a)):
            if j in used or str(ej).upper() != str(e).upper():
                continue
            d = float(np.linalg.norm(x - y))
            if d < best_d:
                best_d, best_j = d, j
        if best_j is None or best_d > max_dist:
            raise RuntimeError(
                f"Could not match reference {e} atom onto the current "
                f"structure (best distance {best_d:.3f} A). "
                "The topology may have changed, not just an extra hydrogen."
            )
        used.add(best_j)
        keep.append(best_j)
    elems = [cur_e[j] for j in keep]
    coords = cur_x[keep]
    return elems, coords


def normalize_element(sym: str) -> str:
    """Return a canonical element symbol (``H``, ``Cl``, …)."""
    s = "".join(c for c in str(sym) if c.isalpha())
    if not s:
        return "C"
    if len(s) == 1:
        return s.upper()
    two = s[0].upper() + s[1].lower()
    if two in ATOMIC_NUMBERS:
        return two
    return s[0].upper()


def electron_count(elements, charge: int = 0) -> int:
    """Valence electron count: sum of atomic numbers minus ``charge``."""
    total = 0
    for el in elements:
        key = normalize_element(el)
        if key not in ATOMIC_NUMBERS:
            raise ValueError(f"Unknown element symbol {el!r}")
        total += ATOMIC_NUMBERS[key]
    return int(total) - int(charge)


def closed_shell_ok(elements, charge: int = 0, multiplicity: int = 1) -> bool:
    """True if ``charge`` / ``multiplicity`` match the element list."""
    n_e = electron_count(elements, charge)
    unpaired = int(multiplicity) - 1
    return unpaired >= 0 and n_e >= unpaired and (n_e - unpaired) % 2 == 0


def anion_extra_hydrogen_indices(
    atoms: list[tuple[str, np.ndarray]],
    net_charge: int = 0,
    multiplicity: int = 1,
) -> list[int]:
    """Indices of OH hydrogens on sulfate/carboxylate that make a closed-shell anion impossible.

    Only drops hydrogens when the current electron count is already illegal
    (e.g. SDS ``C12H25SO4-`` plus a sulfate proton → 147 electrons, singlet).
    Neutral acids such as ``ROSO3H`` are left alone.
    """
    if not atoms:
        return []
    elems = [normalize_element(e) for e, _ in atoms]
    if closed_shell_ok(elems, net_charge, multiplicity):
        return []
    xyz = np.asarray([x for _, x in atoms], dtype=float)
    h_idx = [i for i, e in enumerate(elems) if e == "H"]
    candidates: list[tuple[int, float, int]] = []
    for i in h_idx:
        d = np.linalg.norm(xyz - xyz[i], axis=1)
        d[i] = np.inf
        j = int(np.argmin(d))
        if elems[j] != "O" or d[j] > 1.25:
            continue
        d_o = np.linalg.norm(xyz - xyz[j], axis=1)
        d_o[j] = np.inf
        d_o[i] = np.inf
        k = int(np.argmin(d_o))
        if elems[k] == "S" and d_o[k] < 1.90:
            candidates.append((0, float(d[j]), i))
            continue
        if elems[k] == "C" and d_o[k] < 1.60:
            n_o = sum(
                1
                for t, e in enumerate(elems)
                if e == "O" and float(np.linalg.norm(xyz[t] - xyz[k])) < 1.55
            )
            if n_o >= 2:
                candidates.append((1, float(d[j]), i))
    candidates.sort()
    dropped: list[int] = []
    keep = set(range(len(atoms)))
    for _, _, i in candidates:
        trial = [idx for idx in sorted(keep) if idx != i]
        trial_el = [elems[idx] for idx in trial]
        if closed_shell_ok(trial_el, net_charge, multiplicity):
            dropped.append(i)
            keep.remove(i)
            if closed_shell_ok([elems[idx] for idx in sorted(keep)], net_charge, multiplicity):
                break
    return dropped


def drop_anion_extra_hydrogens(
    atoms: list[tuple[str, np.ndarray]],
    net_charge: int = 0,
    multiplicity: int = 1,
) -> list[tuple[str, np.ndarray]]:
    """Return ``atoms`` without illegal anion extra hydrogens."""
    drop = set(anion_extra_hydrogen_indices(atoms, net_charge, multiplicity))
    if not drop:
        return atoms
    return [atom for i, atom in enumerate(atoms) if i not in drop]


def _gaussian_reference_paths(cwd: Union[Path, str], label: str) -> list[Path]:
    """Candidate original-geometry files, user input first."""
    base = Path(cwd)
    return [
        base / f"{label}.user_input.mol2",
        base / f"{label}.user_input.pdb",
        base / f"{label}.sanitized.pdb",
        base / f"{label}.antechamber_in.mol2",
        base / f"{label}.antechamber_in.pdb",
        base / f"{label}.initial.mol2",
    ]


def atoms_for_gaussian(
    in_path: Union[Path, str],
    cwd: Union[Path, str],
    net_charge: int = 0,
    multiplicity: int = 1,
    logger=None,
) -> tuple[list[str], np.ndarray]:
    """Elements and coords for a Gaussian ``.com``, dropping extra anion H.

    Matching ``initial.mol2`` is not enough when that file is already
    protonated. This also compares to the original user mol2/PDB and, if
    the electron count is still illegal, removes a sulfate/carboxylate OH
    hydrogen.
    """
    in_p = Path(in_path)
    current = parse_structure_atoms(in_p)
    label = in_p.name.split(".")[0]
    charge = int(round(float(net_charge)))
    mult = int(multiplicity)
    for ref_path in _gaussian_reference_paths(cwd, label):
        if not ref_path.is_file():
            continue
        if ref_path.resolve() == in_p.resolve():
            continue
        reference = parse_structure_atoms(ref_path)
        if 0 < len(reference) < len(current):
            if logger is not None:
                logger.warning(
                    "Gaussian geometry has %s atoms but %s has %s; "
                    "dropping unmatched atoms (usually an extra sulfate H)",
                    len(current),
                    ref_path.name,
                    len(reference),
                )
            elements, coords = match_current_to_reference(reference, current)
            current = list(zip(elements, coords))
            break
    dropped = anion_extra_hydrogen_indices(current, charge, mult)
    if dropped:
        if logger is not None:
            logger.warning(
                "Dropping %s extra anion hydrogen(s) so charge %s "
                "multiplicity %s is a closed shell",
                len(dropped),
                charge,
                mult,
            )
        current = drop_anion_extra_hydrogens(current, charge, mult)
    elements = [normalize_element(e) for e, _ in current]
    coords = np.asarray([x for _, x in current], dtype=float)
    if not closed_shell_ok(elements, charge, mult):
        n_e = electron_count(elements, charge)
        raise RuntimeError(
            f"Gaussian geometry has {len(elements)} atoms and {n_e} electrons "
            f"with charge {charge} multiplicity {mult}, which is impossible. "
            "For SDS-like anions this is usually an extra sulfate hydrogen "
            "(SDS is 42 atoms, C12H25SO4-)."
        )
    if logger is not None:
        logger.info("Gaussian atom count for %s: %s", in_p.name, len(elements))
    return elements, coords


def write_simple_mol2(
    path: Union[Path, str],
    elements,
    coords,
    resname: str = "LIG",
) -> Path:
    """Write a bond-free mol2 from element symbols and coordinates."""
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    coords = np.asarray(coords, dtype=float)
    res = (resname or "LIG")[:3].upper()
    lines = [
        "@<TRIPOS>MOLECULE\n",
        f"{res}\n",
        f" {len(elements)} 0 1 0 0\n",
        "SMALL\n",
        "USER_CHARGES\n",
        "\n",
        "@<TRIPOS>ATOM\n",
    ]
    for i, (el, xyz) in enumerate(zip(elements, coords), start=1):
        elem = normalize_element(el)
        lines.append(
            f"{i:7d} {elem:<8s} {xyz[0]:10.4f} {xyz[1]:10.4f} {xyz[2]:10.4f} "
            f"{elem:<7s} 1 {res:<8s} {0.0:10.4f}\n"
        )
    p.write_text("".join(lines), encoding="utf-8")
    return p


@dataclass
class _Mol2Atom:
    atom_id: int
    name: str
    x: float
    y: float
    z: float
    atype: str
    subst_id: str = "1"
    subst_name: str = "LIG"
    charge: str = "0.0000"


@dataclass
class _Mol2Bond:
    bond_id: int
    a: int
    b: int
    order: str


def _parse_mol2_records(path: Union[Path, str]) -> tuple[list[_Mol2Atom], list[_Mol2Bond], str]:
    """Parse ATOM/BOND records and the molecule name."""
    p = Path(path)
    atoms: list[_Mol2Atom] = []
    bonds: list[_Mol2Bond] = []
    name = p.stem
    section = ""
    mol_lines: list[str] = []
    with open(p, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("@<TRIPOS>"):
                section = line.strip()[9:]
                continue
            if section == "MOLECULE":
                mol_lines.append(line.rstrip("\n"))
                continue
            if not line.strip():
                continue
            parts = line.split()
            if section == "ATOM":
                if len(parts) < 6:
                    continue
                try:
                    atom_id = int(parts[0])
                except ValueError:
                    continue
                atoms.append(
                    _Mol2Atom(
                        atom_id=atom_id,
                        name=parts[1],
                        x=float(parts[2]),
                        y=float(parts[3]),
                        z=float(parts[4]),
                        atype=parts[5],
                        subst_id=parts[6] if len(parts) > 6 else "1",
                        subst_name=parts[7] if len(parts) > 7 else "LIG",
                        charge=parts[8] if len(parts) > 8 else "0.0000",
                    )
                )
            elif section == "BOND":
                if len(parts) < 4:
                    continue
                try:
                    bonds.append(
                        _Mol2Bond(
                            bond_id=int(parts[0]),
                            a=int(parts[1]),
                            b=int(parts[2]),
                            order=parts[3],
                        )
                    )
                except ValueError:
                    continue
    if mol_lines:
        name = mol_lines[0].strip() or name
    return atoms, bonds, name


def _write_mol2_records(
    path: Union[Path, str],
    name: str,
    atoms: list[_Mol2Atom],
    bonds: list[_Mol2Bond],
) -> Path:
    """Write a Sybyl mol2 from parsed ATOM/BOND records."""
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "@<TRIPOS>MOLECULE\n",
        f"{name}\n",
        f" {len(atoms)} {len(bonds)} 1 0 0\n",
        "SMALL\n",
        "USER_CHARGES\n",
        "\n",
        "@<TRIPOS>ATOM\n",
    ]
    id_map = {}
    for i, atom in enumerate(atoms, start=1):
        id_map[atom.atom_id] = i
        lines.append(
            f"{i:7d} {atom.name:<8s} {atom.x:10.4f} {atom.y:10.4f} "
            f"{atom.z:10.4f} {atom.atype:<7s} {str(atom.subst_id):>3s} "
            f"{atom.subst_name:<8s} {atom.charge:>10s}\n"
        )
    lines.append("@<TRIPOS>BOND\n")
    bi = 0
    for bond in bonds:
        a = id_map.get(bond.a)
        b = id_map.get(bond.b)
        if a is None or b is None:
            continue
        bi += 1
        lines.append(f"{bi:6d} {a:4d} {b:4d} {bond.order}\n")
    p.write_text("".join(lines), encoding="utf-8")
    return p


def mol2_bond_count(path: Union[Path, str]) -> int:
    """Number of ``@<TRIPOS>BOND`` records in a mol2 file."""
    _, bonds, _ = _parse_mol2_records(path)
    return len(bonds)


def _retype_terminal_sulfate_oxygens(atoms: list[_Mol2Atom], bonds: list[_Mol2Bond]) -> None:
    """Set unprotonated terminal sulfate oxygens to ``O.co2`` (not alcohol)."""
    by_id = {a.atom_id: a for a in atoms}
    adj: dict[int, list[int]] = defaultdict(list)
    for bond in bonds:
        adj[bond.a].append(bond.b)
        adj[bond.b].append(bond.a)
    for sulfur in atoms:
        if normalize_element(_element_from_mol2_name_type(sulfur.name, sulfur.atype)) != "S":
            continue
        oxygens = []
        for nid in adj[sulfur.atom_id]:
            neigh = by_id.get(nid)
            if neigh is None:
                continue
            if normalize_element(_element_from_mol2_name_type(neigh.name, neigh.atype)) == "O":
                oxygens.append(neigh)
        if len(oxygens) != 4:
            continue
        for oxygen in oxygens:
            partners = [by_id[n] for n in adj[oxygen.atom_id] if n in by_id]
            heavies = [
                p
                for p in partners
                if normalize_element(_element_from_mol2_name_type(p.name, p.atype)) != "H"
            ]
            hydrogens = [
                p
                for p in partners
                if normalize_element(_element_from_mol2_name_type(p.name, p.atype)) == "H"
            ]
            if len(heavies) == 1 and not hydrogens:
                oxygen.atype = "O.co2"


def sanitize_mol2_ligand(
    src: Union[Path, str],
    dst: Union[Path, str],
    net_charge: int = 0,
    multiplicity: int = 1,
) -> Path:
    """Drop illegal anion extra H and retype terminal sulfate oxygens.

    Open Babel mol2 files often leave the SDS sulfate oxygen as ``O.3``
    (alcohol) or already include the extra proton. Antechamber then writes
    a 43-atom mol2 and Gaussian dies with an odd electron count.
    """
    atoms, bonds, name = _parse_mol2_records(src)
    xyz_atoms = [
        (
            normalize_element(_element_from_mol2_name_type(a.name, a.atype)),
            np.array([a.x, a.y, a.z], dtype=float),
        )
        for a in atoms
    ]
    drop = set(anion_extra_hydrogen_indices(xyz_atoms, int(net_charge), int(multiplicity)))
    if drop:
        drop_ids = {atoms[i].atom_id for i in drop}
        atoms = [a for i, a in enumerate(atoms) if i not in drop]
        bonds = [b for b in bonds if b.a not in drop_ids and b.b not in drop_ids]
    _retype_terminal_sulfate_oxygens(atoms, bonds)
    return _write_mol2_records(dst, name, atoms, bonds)


def _filter_pdb_lines(lines: list[str], keep_atom_indices: list[int]) -> list[str]:
    """Keep selected ATOM/HETATM records and remap CONECT serials."""
    atom_lines = [ln for ln in lines if ln[:6].strip() in ("ATOM", "HETATM")]
    keep_set = set(keep_atom_indices)
    old_serials = []
    for ln in atom_lines:
        try:
            old_serials.append(int(ln[6:11]))
        except ValueError:
            old_serials.append(len(old_serials) + 1)
    serial_map = {}
    new_atom_lines = []
    new_i = 0
    for old_i, ln in enumerate(atom_lines):
        if old_i not in keep_set:
            continue
        new_i += 1
        serial_map[old_serials[old_i]] = new_i
        raw = ln.rstrip("\n")
        if len(raw) < 80:
            raw = raw.ljust(80)
        new_atom_lines.append(f"{raw[:6]}{new_i:5d}{raw[11:]}" + ("\n" if ln.endswith("\n") else ""))
    out: list[str] = []
    atom_i = -1
    keep_set = set(keep_atom_indices)
    atom_iter = iter(new_atom_lines)
    for ln in lines:
        key = ln[:6].strip()
        if key in ("ATOM", "HETATM"):
            atom_i += 1
            if atom_i in keep_set:
                out.append(next(atom_iter))
            continue
        elif key == "CONECT":
            parts = ln.split()
            if len(parts) < 2:
                continue
            try:
                serial = int(parts[1])
            except ValueError:
                continue
            if serial not in serial_map:
                continue
            rec = f"CONECT{serial_map[serial]:5d}"
            n_partners = 0
            for token in parts[2:]:
                try:
                    partner = int(token)
                except ValueError:
                    continue
                if partner in serial_map:
                    rec += f"{serial_map[partner]:5d}"
                    n_partners += 1
            if n_partners:
                out.append(rec + "\n")
        else:
            out.append(ln if ln.endswith("\n") else ln + "\n")
    return out


def sanitize_pdb_ligand(
    src: Union[Path, str],
    dst: Union[Path, str],
    net_charge: int = 0,
    multiplicity: int = 1,
) -> Path:
    """Write an Amber-safe PDB, dropping illegal extra anion hydrogens.

    Open Babel often writes the anionic sulfate oxygen as element ``O1-``
    and lists S=O twice in CONECT. Antechamber then treats that oxygen as
    an alcohol and appends a hydrogen (42-atom SDS -> 43-atom Gaussian).

    This keeps remaining ATOM/HETATM records, rewrites element/charge into
    columns 77-80, uniquifies CONECT partners, and removes a sulfate OH
    hydrogen when ``net_charge`` makes the electron count illegal.
    """
    src_p = Path(src)
    dst_p = Path(dst)
    dst_p.parent.mkdir(parents=True, exist_ok=True)
    out_lines: list[str] = []
    with open(src_p, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            key = line[:6].strip()
            if key in ("ATOM", "HETATM"):
                out_lines.append(_rewrite_pdb_atom_line(line))
            elif key == "CONECT":
                out_lines.extend(_rewrite_pdb_conect_line(line))
            else:
                out_lines.append(line if line.endswith("\n") else line + "\n")
    parsed: list[tuple[str, np.ndarray]] = []
    for line in out_lines:
        if line[:6].strip() not in ("ATOM", "HETATM"):
            continue
        elem = line[76:78].strip() or "C"
        xyz = np.array(
            [float(line[30:38]), float(line[38:46]), float(line[46:54])],
            dtype=float,
        )
        parsed.append((elem, xyz))
    drop = anion_extra_hydrogen_indices(parsed, int(net_charge), int(multiplicity))
    if drop:
        keep = [i for i in range(len(parsed)) if i not in set(drop)]
        out_lines = _filter_pdb_lines(out_lines, keep)
    with open(dst_p, "w", encoding="utf-8") as fh:
        fh.writelines(out_lines)
    return dst_p


def Remove_PDB_CONECT(filename: Union[Path, str], backup: bool = False):
    """Deprecated: previously stripped CONECT, which protonated anions.

    Now sanitizes in place (keeps unique CONECT, fixes ``O1-`` elements).
    Prefer :func:`sanitize_pdb_ligand`, which does not modify the source.
    """
    fn = Path(filename)
    if backup:
        shutil.copyfile(fn, fn.parent / f"input_{fn.name}")
    tmp = fn.with_name(fn.stem + ".sanitized_tmp.pdb")
    sanitize_pdb_ligand(fn, tmp)
    tmp.replace(fn)
    return
