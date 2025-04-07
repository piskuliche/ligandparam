import io
import warnings
from collections import Counter
from pathlib import Path
from typing import Optional, Union, Any, Sequence

import numpy as np
from typing_extensions import override

from rdkit import Chem

from ligandparam.stages import AbstractStage
from ligandparam.utils import stderr_redirector


class SDFToPDB(AbstractStage):

    @override
    def __init__(self, stage_name: str, main_input: Union[Path, str], cwd: Union[Path, str], *args, **kwargs) -> None:
        super().__init__(stage_name, main_input, cwd, *args, **kwargs)
        self.in_sdf = Path(main_input)
        try:
            self.out_pdb = kwargs["out_pdb"]
        except KeyError:
            self.out_pdb = None
        try:
            self.out_mol2 = kwargs["out_mol2"]
        except KeyError:
            self.out_mol2 = None

        if self.out_mol2 is None and self.out_pdb is None:
            err_msg = f"Must provide either out_pdb or out_mol2"
            self.logger.error(err_msg)
            raise ValueError(err_msg)
        self.resname = kwargs.get("resname", "LIG")
        self.removeHs = kwargs.get("removeHs", False)
        self.add_conect = kwargs.get("add_conect", True)
        self.mol_idx = kwargs.get("mol_idx", 0)

    def execute(self, dry_run=False, nproc: Optional[int] = None, mem: Optional[int] = None) -> Any:
        # First, create the molecule
        try:
            mols = Chem.SDMolSupplier(str(self.in_sdf), removeHs=False)
        except Exception as e:
            err_msg = f"Failed to generate an rdkit molecule from input SDF {self.in_sdf}. Got exception: {e}"
            self.logger.error(err_msg)
            raise RuntimeError(err_msg)

        mol = mols[self.mol_idx]
        # Set metadata and write away
        mol.SetProp("_Name", self.resname)
        mi = Chem.AtomPDBResidueInfo()
        mi.SetResidueName(self.resname)
        mi.SetResidueNumber(1)
        mi.SetOccupancy(0.0)
        mi.SetTempFactor(0.0)
        [a.SetMonomerInfo(mi) for a in mol.GetAtoms()]
        flavor = 0 if self.add_conect else 2
        if self.out_pdb is not None:
            self.write_pdb(mol, flavor=flavor)

    def write_pdb(self, mol: Chem.Mol, flavor: int = 0):
        self.logger.info(f"Writing {self.in_sdf} to {self.out_pdb}")
        try:
            Chem.MolToPDBFile(mol, str(self.out_pdb), flavor=flavor)
        except Exception as e:
            self.logger.error(
                f"Failed to write to  {self.out_pdb}. Got exception: {e}")

    def _append_stage(self, stage: "AbstractStage") -> "AbstractStage":
        raise NotImplementedError

    def _clean(self):
        raise NotImplementedError



# noinspection DuplicatedCode
class SDFToPDBBatch(AbstractStage):
    
    @override
    def __init__(self, stage_name: str, main_input: Union[Path, str], cwd: Union[Path, str], *args, **kwargs) -> None:
        super().__init__(stage_name, main_input, cwd, *args, **kwargs)
        self.in_sdf = Path(main_input)

        self.removeHs = kwargs.get("removeHs", False)
        self.add_conect = kwargs.get("add_conect", True)

        self.out_pdb_template = kwargs.get("out_pdb_template", None)
        self.out_pdbs = kwargs.get("out_pdbs", None)
        self.out_pdb_read_field = kwargs.get("out_pdb_read_field", "_Name")

        self.resnames = kwargs.get("resnames", None)
        self.resname = kwargs.get("resname", None)
        self.resname_read_field = kwargs.get("resname_read_field", "_Name")

    def execute(self, dry_run=False, nproc: Optional[int] = None, mem: Optional[int] = None) -> Any:
        # First, create the molecule
        try:
            mols = Chem.SDMolSupplier(str(self.in_sdf), removeHs=False)
        except Exception as e:
            err_msg = f"Failed to generate an rdkit molecule from input SDF {self.in_sdf} Got exception: {e}"
            self.logger.error(err_msg)
            raise RuntimeError(err_msg)

        # Set up names and paths
        if self.resnames is None:
            if self.resname is None:
                self.resnames = [mol.GetProp(self.resname_read_field)[:3] for mol in mols]
            elif self.resname:
                self.resnames = [self.resname for _ in mols]
        if self.out_pdbs is None:
            if self.out_pdb_template is None:
                filenames = [f'{mol.GetProp(self.out_pdb_read_field)}.pdb' for mol in mols]
                counts = Counter(filenames)
                if np.all(np.array(list(counts.values())) == 1):
                    self.out_pdbs = [self.cwd / fn for fn in filenames]
                else:
                    err_msg = f"Multiple molecules with the same name in {self.in_sdf} Please provide `out_pdbs` or `out_pdb_template`."
                    self.logger.error(err_msg)
                    raise ValueError(err_msg)
            else:
                out_dir = self.out_pdb_template.parent
                label = self.out_pdb_template.stem
                self.out_pdbs = [out_dir / f"{label}_{i}.pdb" for i in range(0, len(self.resnames))]

        if len(self.resnames) != len(self.out_pdbs) or len(self.resnames) != len(mols):
            err_msg = f"Lengths of `out_pdbs`, `resnames`, and mols don't match: {len(self.out_pdbs)}, {len(self.resnames)}, and {len(mols)}"
            self.logger.error(err_msg)
            raise ValueError(err_msg)

        # Write each mol to a different PDB
        flavor = 0 if self.add_conect else 2
        for mol, pdb, resname in zip(mols, self.out_pdbs, self.resnames):
            # Set metadata and write away
            mol.SetProp("_Name", resname)
            mi = Chem.AtomPDBResidueInfo()
            mi.SetResidueName(resname)
            mi.SetResidueNumber(1)
            mi.SetOccupancy(0.0)
            mi.SetTempFactor(0.0)
            [a.SetMonomerInfo(mi) for a in mol.GetAtoms()]
            self.logger.info(f"Writing {self.in_sdf} to {pdb}")

            try:
                Chem.MolToPDBFile(mol, str(pdb), flavor=flavor)
            except Exception as e:
                self.logger.error(
                    f"Failed to write to  {pdb}. Got exception: {e}")

    def _append_stage(self, stage: "AbstractStage") -> "AbstractStage":
        raise NotImplementedError

    def _clean(self):
        raise NotImplementedError



############
# Known issues:
#   - openbabel adds a 1 to the residue name and there's nothing to do about that. If we really need these classes,
#       (SDFToPDBBatch and SDFToPDB), we'll have to edit each output mol2 and remove the resname
#   - openbabel doesn't error out if the output path is not valid.
############

# class SDFToMol2(AbstractStage):
#     from openbabel import pybel
#     @override
#     def __init__(self, stage_name: str, main_input: Union[Path, str], cwd: Union[Path, str], *args, **kwargs) -> None:
#         super().__init__(stage_name, main_input, cwd, *args, **kwargs)
#         self.in_sdf = Path(main_input)
#         self.out_mol2 = kwargs["out_mol2"]

#         self.resname = kwargs.get("resname", "LIG")
#         self.add_hydrogens = kwargs.get("add_hydrogens", False)
#         self.overwrite = kwargs.get("overwrite", True)
#         self.mol_idx = kwargs.get("mol_idx", 0)

#     def execute(self, dry_run=False, nproc: Optional[int] = None, mem: Optional[int] = None) -> Any:
#         # First, create the molecule
#         try:
#             mols = pybel.readfile("sdf", str(self.in_sdf))
#         except Exception as e:
#             err_msg = f"Failed to generate an rdkit molecule from input SDF {self.in_sdf}. Got exception: {e}"
#             self.logger.error(err_msg)
#             raise RuntimeError(err_msg)

#         [next(mols) for _ in range(self.mol_idx)]
#         mol = next(mols)
#         mol.title = self.resname
#         for atm in mol.atoms:
#             atm.residue.OBResidue.SetName(self.resname)
#         if self.add_hydrogens:
#             mol.OBMol.AddHydrogens()
#         self.logger.debug(f"Writing {self.in_sdf} to {self.out_mol2}")
#         try:
#             mol.write("mol2", str(self.out_mol2), overwrite=self.overwrite)
#         except Exception as e:
#             self.logger.error(
#                 f"Failed to write to  {self.out_mol2}. Got exception: {e}")

#     def _append_stage(self, stage: "AbstractStage") -> "AbstractStage":
#         raise NotImplementedError

#     def _clean(self):
#         raise NotImplementedError


# # noinspection DuplicatedCode
# class SDFToMol2Batch(AbstractStage):
#     from openbabel import pybel
#     @override
#     def __init__(self, stage_name: str, main_input: Union[Path, str], cwd: Union[Path, str], *args, **kwargs) -> None:
#         super().__init__(stage_name, main_input, cwd, *args, **kwargs)
#         self.in_sdf = Path(main_input)

#         self.add_hydrogens = kwargs.get("add_hydrogens", False)
#         self.overwrite = kwargs.get("overwrite", True)

#         self.out_mol2_template = kwargs.get("out_mol2_template", None)
#         self.out_mol2s = kwargs.get("out_mol2s", None)
#         self.out_mol2_read_field = kwargs.get("out_mol2_read_field", "ID")

#         self.resnames = kwargs.get("resnames", None)
#         self.resname = kwargs.get("resname", None)
#         self.resname_read_field = kwargs.get("resname_read_field", "ID")

#     def execute(self, dry_run=False, nproc: Optional[int] = None, mem: Optional[int] = None) -> Any:
#         # First, create the molecule
#         try:
#             mols = pybel.readfile("sdf", str(self.in_sdf))
#         except Exception as e:
#             err_msg = f"Failed to generate an rdkit molecule from input SDF {self.in_sdf} Got exception: {e}"
#             self.logger.error(err_msg)
#             raise RuntimeError(err_msg)

#         f = io.BytesIO()
#         with stderr_redirector(f):
#             mols = list(mols)
#         self.setup_resnames_and_paths(mols)

#         # Write each mol to a different mol2
#         for mol, out_mol2, resname in zip(mols, self.out_mol2s, self.resnames):
#             if self.add_hydrogens:
#                 mol.OBMol.AddHydrogens()

#             mol.title = resname
#             for atm in mol.atoms:
#                 atm.residue.OBResidue.SetName(resname)

#             self.logger.debug(f"Writing {self.in_sdf} to {out_mol2}")
#             try:
#                 mol.write("mol2", str(out_mol2), overwrite=self.overwrite)
#             except Exception as e:
#                 self.logger.error(
#                     f"Failed to write to  {out_mol2}. Got exception: {e}")


#     def setup_resnames_and_paths(self, mols: Sequence) -> None:
#         if self.resnames is None:
#             if self.resname is None:
#                 self.resnames = [mol.data[self.resname_read_field][:3] for mol in mols]
#             elif self.resname:
#                 self.resnames = [self.resname for _ in mols]
#         if self.out_mol2s is None:
#             if self.out_mol2_template is None:
#                 filenames = [f'{mol.data[self.out_mol2_read_field]}.mol2' for mol in mols]
#                 counts = Counter(filenames)
#                 if np.all(np.array(list(counts.values())) == 1):
#                     self.out_mol2s = [self.cwd / fn for fn in filenames]
#                 else:
#                     err_msg = f"Multiple molecules with the same name in {self.in_sdf} Please provide `out_mol2s` or `out_mol2_template`."
#                     self.logger.error(err_msg)
#                     raise ValueError(err_msg)
#             else:
#                 out_dir = self.out_mol2_template.parent
#                 label = self.out_mol2_template.stem
#                 self.out_mol2s = [out_dir / f"{label}_{i}.mol2" for i in range(0, len(self.resnames))]

#         if len(self.resnames) != len(self.out_mol2s) or len(self.resnames) != len(mols):
#             err_msg = f"Lengths of `out_mol2s`, `resnames`, and mols don't match: {len(self.out_mol2s)}, {len(self.resnames)}, and {len(mols)}"
#             self.logger.error(err_msg)
#             raise ValueError(err_msg)


#     def _append_stage(self, stage: "AbstractStage") -> "AbstractStage":
#         raise NotImplementedError


#     def _clean(self):
#         raise NotImplementedError
