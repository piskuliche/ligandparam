import os
from typing import Optional,  Union, Any
import logging
import warnings

import MDAnalysis as mda
from MDAnalysis.topology.guessers import guess_atom_element

from pathlib import Path
import shutil as sh
from ase.io import read
from ase.optimize import BFGS, FIRE
from ase.calculators.calculator import Calculator, all_changes

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
        self.charge = kwargs.get("charge", 0)

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
        if not getattr(self, "coord_object", None):
            self.coord_object = Coordinates(self.in_mol2, filetype="mol2")
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
        optimizer = BFGS(atoms, maxstep=0.1)
        optimizer.run(fmax=0.005, steps=self.steps)
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
                mlp_calc = None
                mlp = DPModel(self.model)
                return QDpi2Calculator(mlp, self.charge)
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


class DPModel(object):
    """ This class is a wrapper for the DeepMD model that uses xtb + deepmd to calculate the energy and forces of a molecule.
    
    This script was originally written by Timothy Giese at Rutgers University.
    
    Parameters
    ----------
    fname : str
        The path to the DeepMD model file (usually a .pb file).
    
    Attributes
    ----------
    dp : DeepPot
        The DeepMD potential object.
    cell : None
        The cell parameters, if any.
    rcut : float
        The cutoff radius for the potential.
    ntypes : int
        The number of atom types in the potential.
    tmap : list
        The type map, mapping element symbols to indices.
    
    Methods
    -------
    GetTypeIdxFromSymbol(ele):
        Returns the index of the atom type for a given element symbol.
    GetTypeIdxs(eles):
        Returns a list of indices for the atom types corresponding to a list of element symbols.
    CalcEne(eles, crds):
        Calculates the energy and forces for a given list of element symbols and their coordinates.
    
    """
    def __init__(self,fname):
        
        from deepmd.infer import DeepPot
        #from deepmd.env import reset_default_tf_session_config

        #try:
        self.dp = DeepPot(fname)
        #except:
        #    reset_default_tf_session_config(True)
        #    self.dp = DeepPot(fname)

        self.cell   = None
        self.rcut   = self.dp.get_rcut()
        self.ntypes = self.dp.get_ntypes()
        self.tmap   = self.dp.get_type_map()

    def GetTypeIdxFromSymbol(self,ele):
        """
        Returns the index of the atom type for a given element symbol.
        
        Parameters
        ----------
        ele : str
            The element symbol (e.g., 'H', 'C', 'O').
        Returns
        -------
        int or None
            The index of the atom type if found, otherwise None.
        
        """
        idx = None
        if ele in self.tmap:
            idx = self.tmap.index(ele)
        return idx

    def GetTypeIdxs(self,eles):
        """ 
        Returns a list of indices for the atom types corresponding to a list of element symbols.
        Parameters
        ----------
        eles : list of str
            A list of element symbols (e.g., ['H', 'C', 'O']).
        
        Returns
        -------
        list of int
            A list of indices corresponding to the atom types for the given element symbols.
        
        """
        return [ self.GetTypeIdxFromSymbol(ele) for ele in eles ]

    def CalcEne(self,eles,crds):
        """ Calculates the energy and forces for a given list of element symbols and their coordinates.
        
        Parameters
        ----------
        eles : list of str
            A list of element symbols (e.g., ['H', 'C', 'O']
        crds : np.ndarray
            A numpy array of shape (n, 3) containing the coordinates of the atoms,
            where n is the number of atoms. 
            
        Returns
        -------
        tuple
            A tuple containing:
            - energy (float): The calculated energy in eV.
            - forces (np.ndarray): A numpy array of shape (n, 3) containing
              the calculated forces on the atoms in eV/angstrom.

        Raises
        ------
        None

        """

        import numpy as np


        energy_convert = 1
        force_convert = 1
        
        coord = np.array(crds).reshape([1, -1])
        atype = self.GetTypeIdxs(eles)
        e, f, v = self.dp.eval(coord, self.cell, atype)
        f = f[0] * force_convert
        e = e[0][0] * energy_convert
        
        return e,f
    


class QDpi2Calculator(Calculator):
    """ A calculator that uses DeepMD and xtb to calculate the energy and forces of a molecule.
    
    This class is a wrapper for the DeepMD model that uses xtb + deepmd to calculate the energy and forces of a molecule.
    It combines the results from a DeepMD model and an xtb calculation.
    
    Parameters
    ----------
    dpmodel : DPModel
        The DeepMD model object.
    charge : int
        The charge of the system.
    **kwargs : dict
        Additional keyword arguments to be passed to the Calculator base class.
    
    Attributes
    ----------
    implemented_properties : list
        A list of properties that this calculator can compute, including 'energy', 'forces', and 'free_energy'.
    nolabel : bool
        A flag indicating that this calculator does not label atoms.   
    
    """

    implemented_properties = ['energy','forces','free_energy']
    nolabel=True
    
    def __init__(self,dpmodel,charge,**kwargs):
        from tblite.ase import TBLite
        self.dpmodel = dpmodel
        self.charge = charge
        #self.xtbcalc = XTBCalculator(charge=charge,method="GFN2-xTB")
        self.xtbcalc = TBLite(method="GFN2-xTB",charge=self.charge)
        Calculator.__init__(self,**kwargs)
        
    def calculate(self,
                  atoms=None,
                  properties=None,
                  system_changes=all_changes):
        """ 
        Calculate the energy and forces of the given atoms using DeepMD and xtb.
        
        Parameters
        ----------
        atoms : ase.Atoms
            The atoms object containing the molecular structure.
        properties : list, optional
            A list of properties to calculate. If None, all implemented properties are calculated.
        system_changes : list, optional
            A list of system changes that may affect the calculation. Defaults to all_changes.
        
        Raises
        ------
        ImportError
            If the required libraries (DeepMD, TBLite) are not installed.
        ValueError
            If the DeepMD model is not provided or is invalid.
            
        """
        import numpy as np
        import ase
        if properties is None:
            properties = self.implemented_properties
        Calculator.calculate(self, atoms, properties, system_changes)
        natoms = len(self.atoms)
        positions = self.atoms.positions
        forces = np.zeros((natoms, 3))
        energy = 0

        eles = self.atoms.get_chemical_symbols()
        crds = self.atoms.get_positions()
        
        e,f = self.dpmodel.CalcEne(eles,crds)
        
        atlist = "".join( ["%s1"%(ele) for ele in eles ] )
        atoms = ase.Atoms(atlist,positions=crds)
        atoms.calc =  self.xtbcalc
        e2 = atoms.get_potential_energy()
        f2 = atoms.get_forces()

        energy = e + e2
        forces = f + f2

        self.results['energy'] = energy
        self.results['free_energy'] = energy
        self.results['forces'] = forces