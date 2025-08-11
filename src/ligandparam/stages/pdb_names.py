from pathlib import Path
from typing import Optional, Union, Any

from ligandparam.stages import AbstractStage
from rdkit import Chem
from rdkit.Chem import rdFMCS
from rdkit.Chem.AllChem import ETKDGv3, EmbedMolecule, AlignMol
from typing_extensions import override

from ligandparam.stages import set_atom_pdb_info


class PDB_Name_Fixer(AbstractStage):
    """
    Stage for fixing PDB atom names and related metadata in ligand files.

    Parameters
    ----------
    stage_name : str
        Name of the stage.
    main_input : Union[Path, str]
        Path to the input PDB file.
    cwd : Union[Path, str]
        Current working directory.
    out_pdb : str
        Path to the output PDB file (from kwargs).
    resname : str, optional
        Residue name to use (default is 'LIG').
    reduce : bool, optional
        Whether to reduce hydrogens (default is True).
    add_conect : bool, optional
        Whether to add CONECT records (default is True).
    random_seed : int, optional
        Random seed for coordinate generation (default is None).
    reference_pdb : str, optional
        Path to reference PDB for atom name normalization (from kwargs).
    align : bool, optional
        Whether to align molecule to reference (default is False).
    """

    @override
    def __init__(self, stage_name: str, main_input: Union[Path, str], cwd: Union[Path, str], *args, **kwargs) -> None:
        """
        Initialize the PDB_Name_Fixer stage.

        Parameters
        ----------
        stage_name : str
            Name of the stage.
        main_input : Union[Path, str]
            Path to the input PDB file.
        cwd : Union[Path, str]
            Current working directory.
        *args
            Additional positional arguments.
        **kwargs
            Additional keyword arguments. Must include 'out_pdb'.
        """
        super().__init__(stage_name, main_input, cwd, *args, **kwargs)
        self.in_pdb = main_input
        self.out_pdb = Path(kwargs["out_pdb"])
        self.resname = kwargs.get("resname", "LIG")
        self.reduce = kwargs.get("reduce", True)
        self.add_conect = kwargs.get("add_conect", True)
        self.random_seed = kwargs.get("random_seed", None)

        try:
            self.reference_pdb = Path(kwargs["reference_pdb"]).resolve()
            self.add_required(self.reference_pdb)
            self.normalize_atom_names = True
            self.align = kwargs.get("align", False)
        except KeyError:
            self.normalize_atom_names = False
            self.align = False

    def execute(self, dry_run=False, nproc: Optional[int] = None, mem: Optional[int] = None) -> Any:
        """
        Execute the stage: fix atom names, normalize to reference, and write output PDB.

        Parameters
        ----------
        dry_run : bool, optional
            If True, do not perform actual execution (default is False).
        nproc : int, optional
            Number of processors to use (default is None).
        mem : int, optional
            Memory to use in MB (default is None).

        Returns
        -------
        Any
            Output of the execution (None).
        """
        super()._setup_execution(dry_run=dry_run, nproc=nproc, mem=mem)
        # First, create the molecule
        try:
            mol = Chem.MolFromPDBFile(str(self.in_pdb), removeHs=False)
        except Exception as e:
            err_msg = f"Failed to generate an rdkit molecule from input SMILES {self.in_pdb}. Got exception: {e}"
            self.logger.error(err_msg)
            raise RuntimeError(err_msg)

        
        #if self.reduce:
        #    mol = Chem.rdmolops.AddHs(mol)
        # All the atoms have their coordinates set to zero. Come up with some values
        #params = ETKDGv3()
        #if self.random_seed:
        #    params.randomSeed = self.random_seed
        #EmbedMolecule(mol, params)

        # Set metadata
        mol = set_atom_pdb_info(mol, self.resname)

        # Normalize the molecule to match the reference PDB
        if self.normalize_atom_names:
            mol = self.normalize_to_reference(mol, self.reference_pdb, self.align)

        flavor = 0 if self.add_conect else 2
        self.logger.info(f"Writing {self.in_pdb} to {self.out_pdb}")

        try:
            Chem.MolToPDBFile(mol, str(self.out_pdb), flavor=flavor)
        except Exception as e:
            self.logger.error(
                f"Failed to write to  {self.out_pdb}. Got exception: {e}")

    def _append_stage(self, stage: "AbstractStage") -> "AbstractStage":
        """
        Not implemented. Appends a stage to the workflow.

        Parameters
        ----------
        stage : AbstractStage
            Stage to append.

        Returns
        -------
        AbstractStage
            The appended stage.
        """
        raise NotImplementedError

    def _clean(self):
        """
        Not implemented. Cleans up after stage execution.
        """
        raise NotImplementedError

    def normalize_to_reference(self, mol: Chem.Mol, reference_pdb: Path, align: bool = False) -> Chem.Mol:
        """
        Normalize atom names in a molecule to match a reference PDB file.

        Parameters
        ----------
        mol : Chem.Mol
            Molecule to normalize.
        reference_pdb : Path
            Path to reference PDB file.
        align : bool, optional
            Whether to align molecule coordinates to reference (default is False).

        Returns
        -------
        Chem.Mol
            Normalized molecule.
        """
        # Normalize the atom names to match the reference PDB
        ref_mol = Chem.MolFromPDBFile(str(reference_pdb), removeHs=False)
        if not ref_mol:
            raise ValueError(f"Failed to read reference PDB file {reference_pdb}")
        if len([at for at in ref_mol.GetAtoms() if at.GetAtomicNum() == 1]) == 0:
            self.logger.warn(
                f"Reference '{reference_pdb}' does not contain any hydrogen atoms. It's not a good reference PDB.")

        mcs_mol = self.get_mcs_mol(ref_mol, mol)
        mol_match = mol.GetSubstructMatch(mcs_mol)
        ref_match = ref_mol.GetSubstructMatch(mcs_mol)
        assert len(mol_match) == len(ref_match), \
            f"Mismatch in number of common atoms. {len(mol_match)} vs {len(ref_match)}. This is likely a bug."

        # Get the mapping of the atoms in the target molecule to the reference molecule
        atom_map = list(zip(mol_match, ref_match))
        dict_atom_map = dict(zip(mol_match, ref_match))
        ref_atoms = list(ref_mol.GetAtoms())

        # Atom names will come, either from the reference molecule or from the available names per element.
        # We do this to avoid repeating atom names in the output molecule.
        available_names_per_element = self.get_available_names_per_element(ref_mol, ref_match, mol)
        for a in mol.GetAtoms():
            pdb_info = a.GetPDBResidueInfo()
            if a.GetIdx() in dict_atom_map:
                name = ref_atoms[dict_atom_map[a.GetIdx()]].GetPDBResidueInfo().GetName()
            else:
                name = available_names_per_element[a.GetAtomicNum()].pop(0)
            pdb_info.SetName(self.pad_atom_name(name))
            a.SetMonomerInfo(pdb_info)
        if align:
            AlignMol(mol, ref_mol, atomMap=atom_map)
        return mol

    @staticmethod
    def pad_atom_name(name) -> str:
        """
        Pad atom name to 4 characters for PDB format.

        Parameters
        ----------
        name : str
            Atom name.

        Returns
        -------
        str
            Padded atom name.
        """
        name = f" {name}" if len(name) < 4 else name
        return name.ljust(4)

    
    def get_available_names_per_element(self, ref_mol: Chem.Mol, ref_match, mol: Chem.Mol) -> dict[int, list[str]]:
        """
        Get available atom names per element from reference and target molecules.

        Parameters
        ----------
        ref_mol : Chem.Mol
            Reference molecule.
        ref_match : list
            Atom indices matched in reference molecule.
        mol : Chem.Mol
            Target molecule.

        Returns
        -------
        dict[int, list[str]]
            Mapping from atomic number to available atom names.
        """
        ref_atoms = list(ref_mol.GetAtoms())
        mol_atoms = list(mol.GetAtoms())
        natoms = len(mol_atoms)
        available_names_per_element = {}

        for atm in mol_atoms:
            element_number, name, number, element = self.get_element_name_and_number(atm)
            if element_number not in available_names_per_element:
                available_names_per_element[element_number] = [ f"{element}{i}" for i in range(1, natoms+1)]

        for idx in ref_match:
            element_number, name, number, element = self.get_element_name_and_number(ref_atoms[idx])
            if element_number not in available_names_per_element:
                available_names_per_element[element_number] = [ f"{element}{i}" for i in range(1, natoms+1)]
            if number not in available_names_per_element[element_number]:
                available_names_per_element[element_number].append(f"{element}{number}")
            try:
                available_names_per_element[element_number].remove(name)
            except ValueError as e:
                print(name)
                print(available_names_per_element)
                print("element_number", element_number, "name", name, "number", number, "element", element)
                self.logger.warn(f"Name '{name}' not found in available names for element {element_number}.")
                raise e

        return available_names_per_element

    @staticmethod
    def get_element_name_and_number(atom) -> tuple[int, str, int, str]:
        """
        Get atomic number, atom name, atom number, and element symbol from an atom.

        Parameters
        ----------
        atom : Chem.Atom
            RDKit atom object.

        Returns
        -------
        tuple
            (atomic number, atom name, atom number, element symbol)
        """
        element_number = atom.GetAtomicNum()
        name = atom.GetPDBResidueInfo().GetName().strip()
        if name == '':
            name = Chem.GetPeriodicTable().GetElementSymbol(element_number) + "0"
            number = 0
            element = Chem.GetPeriodicTable().GetElementSymbol(element_number)
        else:
            number = int(''.join(char for char in name if char.isdigit()))
            element = ''.join(char for char in name if not char.isdigit())
        return element_number, name, number, element
    
    @staticmethod
    def get_mcs_mol(ref_mol: Chem.Mol, mol: Chem.Mol) -> Chem.Mol:
        """
        Get the Maximum Common Substructure (MCS) molecule between reference and target.

        Parameters
        ----------
        ref_mol : Chem.Mol
            Reference molecule.
        mol : Chem.Mol
            Target molecule.

        Returns
        -------
        Chem.Mol
            MCS molecule.
        """
        mcs = rdFMCS.FindMCS([ref_mol, mol])
        common_mol = Chem.rdmolfiles.MolFromSmarts(mcs.smartsString)
        return common_mol
