from pathlib import Path
from typing import Optional, Union

from ligandparam.stages.AbstractStage import AbstractStage

try:
    from openff.interchange import Interchange
    from openff.toolkit import ForceField, Molecule
    from openff.units import unit
except ImportError:
    Interchange = None
    ForceField = None
    Molecule = None
    unit = None


class StageSageCreate(AbstractStage):
    """ Converts the final resulting mol2 to the SAGE forcefield.

    Parameters
    ----------
    stage_name: str
        The name of the stage.
    main_input: Union[Path, str]
        The main input file (mol2 format).
    cwd: Union[Path, str]
        The current working directory.
    out_parm: str
        The output parm file (parm7 format).
    use_mol2_charges: bool, optional
        If True, use charges from the MOL2 atom table and pass them directly
        into OpenFF instead of triggering OpenFF charge assignment.
    resname: str, optional
        Residue name override. If not provided, inferred from the MOL2 atom
        records (`subst_name` column).

    Attributes
    ----------
    in_mol2: Path
        The input mol2 file.
    out_parm: Path
        The output parm file.
    out_rst7: Path
        The output rst7 file.
    ff_name: str
        The name of the force field.
    
    Notes
    -----
    This stage is responsible for converting the final resulting mol2 file to the SAGE forcefield format.

    """
    def __init__(self, stage_name: str, main_input: Union[Path, str], cwd: Union[Path, str], *args, **kwargs) -> None:
        if Interchange is None or ForceField is None or Molecule is None or unit is None:
            raise ImportError(
                "StageSageCreate requires the optional 'openff-toolkit' dependency, which is not installed. "
                "Install it with `pip install ligandparam[sage]` (or `pip install openff-toolkit`)."
            )

        super().__init__(stage_name, main_input, cwd, *args, **kwargs)
        self.in_mol2 = Path(main_input)
        self.out_parm = Path(kwargs["out_parm"])
        self.out_rst7 = Path(kwargs["out_parm"].replace(".parm7", ".rst7"))
        if ".parm7" not in self.out_parm.name:
            raise ValueError("Output parameter file must have .parm7 extension")
        self.ff_name = kwargs.get("ff_name", "openff-2.2.0.offxml")
        self.use_mol2_charges = kwargs.get("use_mol2_charges", False)
        self.resname = kwargs.get("resname", None)
        self.add_required(self.in_mol2)

    @staticmethod
    def _read_mol2_atom_records(mol2_file: Path):
        """Read atom names, optional charges, and residue names from MOL2 ATOM section."""
        atom_names = []
        atom_charges = []
        atom_resnames = []
        in_atom_section = False

        with Path(mol2_file).open("r", encoding="utf-8") as handle:
            for raw_line in handle:
                line = raw_line.strip()

                if not line:
                    continue

                if line.startswith("@<TRIPOS>"):
                    section = line.upper()
                    if section == "@<TRIPOS>ATOM":
                        in_atom_section = True
                        continue

                    if in_atom_section:
                        break
                    continue

                if not in_atom_section:
                    continue

                columns = raw_line.split()
                if len(columns) < 2:
                    raise ValueError(f"Malformed MOL2 atom line in {mol2_file}: {raw_line.rstrip()}")
                atom_names.append(columns[1])
                if len(columns) >= 8:
                    atom_resnames.append(columns[7])
                else:
                    atom_resnames.append(None)
                if len(columns) >= 9:
                    try:
                        atom_charges.append(float(columns[8]))
                    except ValueError as exc:
                        raise ValueError(
                            f"Malformed MOL2 atom charge in {mol2_file}: {raw_line.rstrip()}"
                        ) from exc
                else:
                    atom_charges.append(None)

        if not atom_names:
            raise ValueError(f"No atom names found in MOL2 ATOM section: {mol2_file}")

        return atom_names, atom_charges, atom_resnames

    @staticmethod
    def _infer_resname(atom_resnames, fallback: str = "LIG") -> str:
        """Infer a single residue name from MOL2 atom records."""
        valid = []
        for name in atom_resnames:
            if name is None:
                continue
            cleaned = str(name).strip()
            if not cleaned or cleaned == "****":
                continue
            valid.append(cleaned)

        if not valid:
            return fallback

        return valid[0]

    @staticmethod
    def _apply_resname_to_parmed(parm, resname: Optional[str]) -> None:
        """Apply a single residue name to all residues in a parmed object."""
        if resname is None:
            return

        cleaned = str(resname).strip()
        if not cleaned:
            return

        for residue in parm.residues:
            residue.name = cleaned

        if hasattr(parm, "parm_data") and "RESIDUE_LABEL" in parm.parm_data:
            parm.parm_data["RESIDUE_LABEL"] = [cleaned for _ in parm.parm_data["RESIDUE_LABEL"]]
    
    def _run(self, dry_run=False, nproc: Optional[int]=None, mem: Optional[int]=None) -> None:
        from rdkit import Chem

        atom_names, atom_charges, atom_resnames = self._read_mol2_atom_records(self.in_mol2)
        mol = Chem.MolFromMol2File(str(self.in_mol2), removeHs=False)
        if mol is None:
            raise ValueError(f"RDKit failed to read MOL2 file: {self.in_mol2}")

        if mol.GetNumAtoms() != len(atom_names):
            raise ValueError(
                f"Input MOL2 atom-count mismatch for {self.in_mol2}: "
                f"RDKit sees {mol.GetNumAtoms()} atoms, MOL2 ATOM block has {len(atom_names)} names."
            )

        # Ensure RDKit has per-atom names, then propagate into OpenFF atoms.
        for atom_index, atom_name in enumerate(atom_names):
            mol.GetAtomWithIdx(atom_index).SetProp("_Name", atom_name)

        molecule = Molecule.from_rdkit(mol, allow_undefined_stereo=True, hydrogens_are_explicit=True)

        if molecule.n_atoms != len(atom_names):
            raise ValueError(
                f"OpenFF atom-count mismatch for {self.in_mol2}: "
                f"OpenFF sees {molecule.n_atoms} atoms, MOL2 ATOM block has {len(atom_names)} names."
            )

        for off_atom, atom_name in zip(molecule.atoms, atom_names):
            off_atom.name = atom_name

        if self.resname is None:
            inferred = self._infer_resname(atom_resnames)
            unique_resnames = sorted({r for r in atom_resnames if r is not None})
            if len(unique_resnames) > 1:
                self.logger.warning(
                    f"Multiple MOL2 residue names found in {self.in_mol2}: {unique_resnames}. Using '{inferred}'."
                )
            self.resname = inferred
        else:
            self.resname = str(self.resname).strip()

        if self.resname:
            molecule.name = self.resname

        charge_from_molecules = None
        if self.use_mol2_charges:
            if len(atom_charges) != molecule.n_atoms:
                raise ValueError(
                    f"Input MOL2 charge-count mismatch for {self.in_mol2}: "
                    f"OpenFF sees {molecule.n_atoms} atoms, MOL2 ATOM block has {len(atom_charges)} charges."
                )
            if any(charge is None for charge in atom_charges):
                raise ValueError(
                    f"use_mol2_charges=True requires charge values in the MOL2 ATOM section for all atoms: {self.in_mol2}"
                )
            molecule.partial_charges = unit.Quantity(atom_charges, unit.elementary_charge)
            charge_from_molecules = [molecule]

        topology = molecule.to_topology()
        ff = ForceField(self.ff_name)
        if charge_from_molecules is None:
            interchange = Interchange.from_smirnoff(
                force_field = ff,
                topology = topology
            )
        else:
            interchange = Interchange.from_smirnoff(
                force_field=ff,
                topology=topology,
                charge_from_molecules=charge_from_molecules
            )
        interchange.to_prmtop(f"{self.out_parm}")
        interchange.to_inpcrd(f"{self.out_rst7}")

        # OpenFF writes a generic residue label; rewrite to the input MOL2 resname.
        if self.resname:
            import parmed

            amber = parmed.load_file(str(self.out_parm), xyz=str(self.out_rst7))
            self._apply_resname_to_parmed(amber, self.resname)
            amber.save(str(self.out_parm), overwrite=True)
            amber.save(str(self.out_rst7), overwrite=True)

        return 
    

class StageSageToAmber(AbstractStage):
    """Convert SAGE-generated AMBER files into OFF/lib and frcmod with unique atom types.

    Parameters
    ----------
    stage_name : str
        Name of the stage.
    main_input : Union[Path, str]
        Input `parm7` path.
    cwd : Union[Path, str]
        Working directory.
    in_rst7 : str, optional
        Input `rst7` path. Defaults to `main_input` with `.rst7` suffix.
    out_lib : str
        Output OFF/lib file path.
    out_frcmod : str
        Output frcmod file path.
    type_prefix : str, optional
        Prefix for generated atom types. Defaults to empty.
    type_width : int, optional
        Total AMBER atom-type width to generate. Defaults to ``2`` for tleap/frcmod
        compatibility.
    reference_leaprc : tuple/list, optional
        leaprc files used to collect existing atom types to avoid collisions.
        Relative names/paths are resolved under ``$AMBER_HOME/dat/leap/cmd``.
        Defaults to ``("leaprc.protein.ff19SB", "leaprc.gaff2")``.
    strict_reference_check : bool, optional
        If True, fail when reference leaprc files cannot be loaded.
    resname : str, optional
        Residue name override to apply before writing OFF/lib and frcmod.
    """

    def __init__(self, stage_name: str, main_input: Union[Path, str], cwd: Union[Path, str], *args, **kwargs) -> None:
        super().__init__(stage_name, main_input, cwd, *args, **kwargs)
        self.in_parm = Path(main_input)
        self.in_rst7 = Path(kwargs.get("in_rst7", str(self.in_parm).replace(".parm7", ".rst7")))
        self.out_lib = Path(kwargs["out_lib"])
        self.out_frcmod = Path(kwargs["out_frcmod"])
        self.type_prefix = kwargs.get("type_prefix", "")
        self.type_width = int(kwargs.get("type_width", 2))
        self.reference_leaprc = tuple(kwargs.get("reference_leaprc", ("leaprc.protein.ff19SB", "leaprc.gaff2")))
        self.strict_reference_check = kwargs.get("strict_reference_check", True)
        self.resname = kwargs.get("resname", None)

        if ".parm7" not in self.in_parm.name:
            raise ValueError("Input parameter file must have .parm7 extension")
        if ".rst7" not in self.in_rst7.name:
            raise ValueError("Input coordinate file must have .rst7 extension")
        if ".lib" not in self.out_lib.name:
            raise ValueError("Output library file must have .lib extension")
        if ".frcmod" not in self.out_frcmod.name:
            raise ValueError("Output parameter modification file must have .frcmod extension")
        if self.type_width < 1 or self.type_width > 2:
            raise ValueError("`type_width` must be 1 or 2 for AMBER tleap/frcmod compatibility.")

        self.add_required(self.in_parm)
        self.add_required(self.in_rst7)

    @staticmethod
    def _base36_encode(value: int, width: int) -> str:
        if value < 0:
            raise ValueError("Base36 encoding only supports non-negative integers.")
        alphabet = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        if value == 0:
            encoded = "0"
        else:
            chars = []
            current = value
            while current:
                current, remainder = divmod(current, 36)
                chars.append(alphabet[remainder])
            encoded = "".join(reversed(chars))
        if len(encoded) > width:
            raise ValueError(f"Value {value} does not fit in base36 width {width}.")
        return encoded.rjust(width, "0")

    def _resolve_reference_leaprcs(self):
        """Resolve leaprc references to concrete files under ``$AMBER_HOME``."""
        from contextlib import nullcontext
        import os

        amber_home = os.environ.get("AMBER_HOME")
        if not amber_home:
            raise RuntimeError(
                "AMBER_HOME is not set. Please set AMBER_HOME to your AmberTools/Amber installation root."
            )

        leap_cmd_dir = Path(amber_home) / "dat" / "leap" / "cmd"
        if not leap_cmd_dir.is_dir():
            raise RuntimeError(
                f"AMBER_HOME appears invalid: could not find leap command directory at {leap_cmd_dir}"
            )

        resolved = []
        for leaprc in self.reference_leaprc:
            candidate = str(leaprc).strip()
            if not candidate:
                continue

            candidate_path = Path(candidate)
            if candidate_path.is_absolute():
                resolved_path = candidate_path
            else:
                resolved_path = leap_cmd_dir / candidate_path

            if not resolved_path.is_file():
                raise FileNotFoundError(
                    f"Reference leaprc file not found: {resolved_path} (from '{candidate}')"
                )
            resolved.append(str(resolved_path))

        return resolved, nullcontext()

    @staticmethod
    def _parse_add_atom_types_from_leaprc(leaprc_file: Union[str, Path]) -> set:
        """Parse atom-type names from a leaprc ``addAtomTypes { ... }`` block."""
        import re

        atom_types = set()
        in_block = False
        brace_depth = 0

        with Path(leaprc_file).open("r", encoding="utf-8", errors="ignore") as handle:
            for raw_line in handle:
                # Remove comments and surrounding whitespace.
                line = raw_line.split("#", 1)[0].strip()
                if not line:
                    continue

                if not in_block:
                    if line.lower().startswith("addatomtypes"):
                        in_block = True
                        brace_depth = line.count("{") - line.count("}")
                    else:
                        continue
                else:
                    brace_depth += line.count("{") - line.count("}")

                # Matches lines like: { "CA" "C" "sp2" }
                for match in re.finditer(r'\{\s*"([^"]+)"\s*"[^"]*"\s*"[^"]*"\s*\}', line):
                    atom_types.add(match.group(1).strip().upper())

                if in_block and brace_depth <= 0:
                    in_block = False

        return atom_types

    def _collect_reserved_atom_types(self) -> set:
        """Collect atom types already defined in the reference ff19SB/gaff2 leaprc files."""
        from parmed.amber import AmberParameterSet

        reserved = set()
        loaded = []
        failed = []
        resolved_leaprcs, stack = self._resolve_reference_leaprcs()

        with stack:
            for leaprc in resolved_leaprcs:
                parse_error = None
                parsed_atom_types = set()

                # Prefer lightweight parsing for concrete leaprc files to avoid
                # executing loadamberparams/loadOff commands.
                if Path(leaprc).is_file():
                    try:
                        parsed_atom_types = self._parse_add_atom_types_from_leaprc(leaprc)
                    except Exception as exc:
                        parse_error = str(exc)

                if parsed_atom_types:
                    loaded.append(f"{leaprc} [parsed]")
                    reserved.update(parsed_atom_types)
                    continue

                try:
                    params = AmberParameterSet.from_leaprc(leaprc)
                    loaded.append(leaprc)
                    reserved.update(str(atom_type).strip().upper() for atom_type in params.atom_types.keys())
                except Exception as exc:
                    if parse_error is not None:
                        failed.append((leaprc, f"addAtomTypes parse failed ({parse_error}); from_leaprc failed ({exc})"))
                    elif Path(leaprc).is_file():
                        failed.append((leaprc, f"No addAtomTypes entries parsed; from_leaprc failed ({exc})"))
                    else:
                        failed.append((leaprc, str(exc)))

        for leaprc, err in failed:
            self.logger.warning(f"Could not load reference leaprc '{leaprc}': {err}")

        if self.strict_reference_check and failed:
            failed_names = ", ".join(name for name, _ in failed)
            raise RuntimeError(
                f"Failed to load required reference leaprc files ({failed_names}) for atom-type collision checks."
            )

        if loaded:
            self.logger.info(f"Loaded reference atom types from: {', '.join(loaded)}")
        else:
            self.logger.warning(
                "No reference leaprc files were loaded; using generated atom types to avoid collisions."
            )

        return reserved

    def _generate_unique_atom_types(self, n_atoms: int, reserved_types: set) -> list:
        """Generate n unique AMBER atom types, excluding reserved types."""
        prefix = str(self.type_prefix).strip().upper()
        if len(prefix) > self.type_width - 1:
            raise ValueError(
                f"`type_prefix` must be at most {self.type_width - 1} characters for type_width={self.type_width}."
            )
        if prefix and not prefix.isalnum():
            raise ValueError("`type_prefix` must be alphanumeric.")

        width = self.type_width - len(prefix)
        if width < 1:
            raise ValueError("`type_prefix` leaves no room for unique type suffixes.")

        max_unique = 36 ** width
        if n_atoms > max_unique:
            raise ValueError(
                f"Cannot generate {n_atoms} unique atom types with prefix '{prefix}', type_width={self.type_width}. "
                f"Maximum possible is {max_unique}."
            )

        blocked = {str(atom_type).strip().upper() for atom_type in reserved_types}
        generated = []
        candidate_index = 0

        while len(generated) < n_atoms:
            if candidate_index >= max_unique:
                raise RuntimeError("Exhausted atom-type namespace before assigning all atoms.")
            candidate = f"{prefix}{self._base36_encode(candidate_index, width)}"
            candidate_index += 1
            if candidate in blocked:
                continue
            generated.append(candidate)
            blocked.add(candidate)

        return generated

    @staticmethod
    def _assign_atom_types_unique(parm, unique_types: list) -> None:
        """Assign one atom type per atom while preserving existing per-atom LJ/mass values."""
        import copy

        if len(parm.atoms) != len(unique_types):
            raise ValueError(
                f"Atom-count mismatch while assigning unique types: {len(parm.atoms)} atoms vs {len(unique_types)} types."
            )

        for atom, new_type in zip(parm.atoms, unique_types):
            atom.type = new_type
            if getattr(atom, "atom_type", None) is not None:
                atom.atom_type = copy.copy(atom.atom_type)
                atom.atom_type.name = new_type

    @staticmethod
    def _apply_resname(parm, resname: Optional[str]) -> None:
        """Apply a single residue name to all residues in the parmed object."""
        if resname is None:
            return

        cleaned = str(resname).strip()
        if not cleaned:
            return

        for residue in parm.residues:
            residue.name = cleaned

        if hasattr(parm, "parm_data") and "RESIDUE_LABEL" in parm.parm_data:
            parm.parm_data["RESIDUE_LABEL"] = [cleaned for _ in parm.parm_data["RESIDUE_LABEL"]]

    def _run(self, dry_run=False, nproc: Optional[int]=None, mem: Optional[int]=None) -> None:
        import parmed
        from parmed.tools import actions

        if dry_run:
            self.logger.info(
                f"[dry-run] Would load {self.in_parm}/{self.in_rst7}, assign unique atom types, "
                f"and write {self.out_lib} + {self.out_frcmod}"
            )
            return

        parm = parmed.load_file(str(self.in_parm), xyz=str(self.in_rst7))

        # Ensure coordinates are attached from rst7 for a complete OFF/lib write.
        if hasattr(parm, "load_rst7"):
            parm.load_rst7(str(self.in_rst7))

        self._apply_resname(parm, self.resname)

        reserved_types = self._collect_reserved_atom_types()
        unique_types = self._generate_unique_atom_types(len(parm.atoms), reserved_types)
        self._assign_atom_types_unique(parm, unique_types)

        if hasattr(parm, "remake_parm"):
            parm.remake_parm()

        actions.writeOFF(parm, str(self.out_lib)).execute()
        actions.writeFrcmod(parm, str(self.out_frcmod)).execute()
        return

