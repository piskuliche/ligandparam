import os
from typing import Optional,  Union, Any
import logging
import warnings

import MDAnalysis as mda
from MDAnalysis.topology.guessers import guess_atom_element

from pathlib import Path
import shutil as sh
from ase.io import read
from ase.optimize import BFGS

from ligandparam.stages.abstractstage import AbstractStage
from ligandparam.io.coordinates import Coordinates, SimpleXYZ, Mol2Writer
from ligandparam.io.gaussianIO import GaussianWriter, GaussianInput, GaussianReader
from ligandparam.interfaces import Gaussian, Antechamber
from ligandparam.log import get_logger


class DPMinimize(AbstractStage):
    """
    This class uses DeepMD to minimize the ligand structure.

    Parameters
    ----------

    Returns
    -------

    """
    def __init__(self, stage_name: str, main_input: Union[Path, str], cwd: Union[Path, str], *args, **kwargs) -> None:
        super().__init__(stage_name, main_input, cwd, *args, **kwargs)
        self.in_mol2 = Path(main_input)
        self.out_xyz = Path(kwargs["out_xyz"])
        self.out_mol2 = Path(kwargs["out_mol2"])

        self.model = kwargs.get("model", "deepmd_model.pb")
        self.ftol = kwargs.get("ftol", 0.05)
        self.steps = kwargs.get("steps", 1000)
        if not getattr(self, "coord_object", None):
            self.coord_object = Coordinates(self.in_mol2, filetype="mol2")
    
    def _append_stage(self, stage: "AbstractStage") -> "AbstractStage":
        """Appends the stage.

        Args:
            stage (AbstractStage): _description_

        Returns:
            AbstractStage: _description_
        """
        return stage

    def execute(self, dry_run=False, nproc: Optional[int] = None, mem: Optional[int] = None) -> Any:
        """Execute the Gaussian calculations.

        Parameters
        ----------
        dry_run : bool, optional
            If True, the stage will not be executed, but the function will print the commands that would
        model : str, optional
            The model to be used for the calculations.

        Returns
        -------
        None

        """
        if dry_run:
            print(f"Dry run: would execute with model {self.model}")
            return
        print("Starting execute")
        elements = self.coord_object.u.atoms.elements
        with open("temp.xyz", 'w') as f:
            f.write(f"{len(self.coord_object.u.atoms)}\n\n")
            for atom in self.coord_object.u.atoms:
                f.write(f"{elements[atom.index]} {atom.position[0]} {atom.position[1]} {atom.position[2]}\n")
        calculator = self._choose_calculator()
        try:
            atoms = read("temp.xyz", format='xyz')
        except Exception as e:
            print(f"Error reading input XYZ file: {e}")
            return
        atoms.calc = calculator
        optimizer = BFGS(atoms)
        optimizer.run(fmax=self.ftol, steps=self.steps)
        atoms.write(self.out_xyz, format='xyz')
        self.replace_mol2_coords(self.in_mol2, self.out_xyz, self.out_mol2)
        print(f"Minimized coordinates written to {self.out_xyz} and {self.out_mol2}")

        return
    
    def _choose_calculator(self):
        """Choose the calculator based on the model.
        
        Raises
        ------
        ValueError
            If the model type is unknown.
        ImportError
            If DeepMD or MACE is not installed.
        """
        try:
            if '.pb' in self.model:
                from deepmd.calculator import DP
                return DP(model=self.model)
            elif '.model' in self.model:
                from mace.calculators import MACECalculator
                return MACECalculator(model_paths=self.model)
            else:
                raise ValueError(f"Unknown model type: {self.model}. Expected a .pb file.")
        except ImportError as e:
            raise ImportError("Please install DeepMD or MACE to use this stage.", e) from e
        
    @staticmethod
    def replace_mol2_coords(mol2_in, xyz_in, mol2_out):
        """Replace coordinates in a MOL2 file with those from an XYZ file.

        Parameters
        ----------
        mol2_in : str
            Input MOL2 file path.
        xyz_in : str
            Input XYZ file path containing minimized coordinates.
        mol2_out : str
            Output MOL2 file path where coordinates will be replaced.
        """
        # Read minimized coordinates from XYZ
        with open(xyz_in) as f:
            lines = f.readlines()
            # Skip first two lines (XYZ header)
            xyz_coords = [line.split()[1:4] for line in lines[2:] if line.strip()]
        
        # Read MOL2 and replace coordinates
        with open(mol2_in) as f:
            mol2_lines = f.readlines()
        
        out_lines = []
        in_atom_section = False
        atom_idx = 0
        for line in mol2_lines:
            if line.startswith("@<TRIPOS>ATOM"):
                in_atom_section = True
                out_lines.append(line)
                continue
            if line.startswith("@<TRIPOS>") and in_atom_section:
                in_atom_section = False
            if in_atom_section and line.strip() and not line.startswith("@<TRIPOS>ATOM"):
                parts = line.split()
                if atom_idx < len(xyz_coords):
                    parts[2:5] = xyz_coords[atom_idx]
                    atom_idx += 1
                    parts[2] = "{:.4f}".format(float(parts[2]))
                    parts[3] = "{:.4f}".format(float(parts[3]))
                    parts[4] = "{:.4f}".format(float(parts[4]))
                out_lines.append("{:<7} {:<8} {:>10} {:>10} {:>10} {:<6} {:>3} {:<8} {:>10}\n".format(*parts[:9]))
            else:
                out_lines.append(line)
        
        with open(mol2_out, "w") as f:
            f.writelines(out_lines)

    
    def _clean(self):
        """Clean the files generated during the stage."""
        raise NotImplementedError("clean method not implemented")
