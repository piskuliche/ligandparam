import re
import shutil
import warnings
from typing import Optional, Union, Any
from pathlib import Path
import shutil as cp

import MDAnalysis as mda

from ligandparam.stages.abstractstage import AbstractStage
from ligandparam.interfaces import Antechamber
from ligandparam.io.coordinates import Mol2Writer


class StageUpdate(AbstractStage):
    """This class updates either (or both) the atom types and names in a mol2 file to match another mol2 file.

    Parameters
    ----------

    """

    def __init__(self, stage_name: str, main_input: Union[Path, str], cwd: Union[Path, str], *args, **kwargs) -> None:
        super().__init__(stage_name, main_input, cwd, *args, **kwargs)
        self.in_mol2 = Path(main_input)
        self.out_mol2 = Path(kwargs["out_mol2"])

        self.source_mol2 = Path(kwargs["source_mol2"])
        self.net_charge = kwargs.get("net_charge", 0.0)
        self.update_names = kwargs.get("update_names", False)
        self.update_types = kwargs.get("update_types", False)
        self.update_resname = kwargs.get("update_resname", False)
        self.update_charges = kwargs.get("update_charges", False)
        self.atom_type = kwargs.get("atom_type", "gaff2")
        if "molname" in kwargs:
            self.additional_args = {"rn": kwargs["molname"]}
        else:
            self.additional_args = {}

        self.add_required(Path(self.in_mol2))
        self.add_required(Path(self.source_mol2))

    def _append_stage(self, stage: "AbstractStage") -> "AbstractStage":
        return stage

    def execute(self, dry_run=False, nproc: Optional[int] = None, mem: Optional[int] = None) -> Any:
        super()._setup_execution(dry_run=dry_run, nproc=nproc, mem=mem)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            if not (self.update_names or self.update_types or self.update_charges):
                self.logger.debug("No updates requested. Exiting.")
                return
            if self.update_names and self.update_types:
                self.logger.debug("Both updates requested. This will update both atom names and types.")
            elif self.update_names:
                self.logger.debug("Only updating atom names.")
            elif self.update_types:
                self.logger.debug("Only updating atom types.")

            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                dest_u = mda.Universe(self.in_mol2, format="mol2")
                source_u = mda.Universe(self.source_mol2)
            if self.update_resname:
                dest_u.residues.resnames = source_u.residues.resnames
            for orig_atom, new_atom in zip(source_u.atoms, dest_u.atoms):
                if self.update_types:
                    self.logger.debug(
                        f"Atom {orig_atom.name} with type {orig_atom.type} and will be updated to {new_atom.type}"
                    )
                    new_atom.type = orig_atom.type
                if self.update_names:
                    self.logger.debug(f"Atom {new_atom.name} will named {orig_atom.name}")
                    new_atom.name = orig_atom.name
                if self.update_charges:
                    self.logger.debug(
                        f"Atom {new_atom.name}'s charge ({orig_atom.charge}) will be updated to {orig_atom.name}"
                    )
                    new_atom.charge = orig_atom.charge

            types_mol2 = Path(self.cwd, f"{self.out_mol2.stem}.types.mol2")
            if not dry_run:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    Mol2Writer(dest_u, filename=types_mol2).write()

            ante = Antechamber(cwd=self.cwd, logger=self.logger, nproc=self.nproc)
            ante.call(
                i=types_mol2,
                fi="mol2",
                o=self.out_mol2,
                fo="mol2",
                pf="y",
                at=self.atom_type,
                an="no",
                nc=self.net_charge,
                dry_run=dry_run,
                **self.additional_args,
            )

    def _clean(self):
        raise NotImplementedError


class StageMatchAtomNames(AbstractStage):
    """This class updates the atom names in a mol2 file to match those of another source structure file.
    It does this through text editing.

    Parameters
    ----------

    """

    def __init__(self, stage_name: str, main_input: Union[Path, str], cwd: Union[Path, str], *args, **kwargs) -> None:
        super().__init__(stage_name, main_input, cwd, *args, **kwargs)
        self.in_mol2 = Path(main_input)
        self.out_mol2 = Path(kwargs["out_mol2"])
        self.source_mol = Path(kwargs["source_mol"])
        self.add_required(Path(self.in_mol2))
        self.add_required(Path(self.source_mol))

    def _append_stage(self, stage: "AbstractStage") -> "AbstractStage":
        return stage

    def execute(self, dry_run=False, nproc: Optional[int] = None, mem: Optional[int] = None) -> Any:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            source_u = mda.Universe(self.source_mol)

        source_names = [atom.name for atom in source_u.atoms]

        with open(self.in_mol2, "r") as f:
            lines = f.readlines()


        with open(self.out_mol2, "w") as f:
            ite = iter(lines)
            while True:
                for line in ite:
                    f.write(line)
                    if line.startswith("@<TRIPOS>ATOM"):
                        break
                for line in ite:
                    match = re.search(r"^\s*\S+\s+(\S+)", line)
                    if match:
                        try:
                            name_idx = match.start(1)
                            new_name = source_names.pop(0)
                            line = line[:name_idx] + new_name + line[name_idx + len(new_name):]
                        except IndexError:
                            self.logger.warning(
                                f"Source structure ({self.source_mol}) has fewer atoms than input mol2 file ({self.in_mol2})")
                    f.write(line)
                break


    def _clean(self):
        raise NotImplementedError
