import io
import json
from pathlib import Path
from typing import Optional,  Union, Any

import requests
from ligandparam.stages import AbstractStage
from ligandparam.utils import stderr_redirector
from rdkit import Chem
from rdkit.Chem import Mol
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.AllChem import AssignBondOrdersFromTemplate
from typing_extensions import override


class LigHFix(AbstractStage):
    @override
    def __init__(self, stage_name: str, main_input: Union[Path, str], cwd: Union[Path, str], *args, **kwargs) -> None:
        super().__init__(stage_name, main_input, cwd, *args, **kwargs)
        self.lig_id = main_input
        self.in_pdb = Path(kwargs["in_pdb"])
        self.out_pdb = Path(kwargs["out_pdb"])
        self.resname = kwargs.get("resname", "LIG")
        self.reduce = kwargs.get("reduce", True)
        self.add_conect = kwargs.get("add_conect", True)
        self.random_seed = kwargs.get("random_seed", None)

    def execute(self, dry_run=False, nproc: Optional[int]=None, mem: Optional[int]=None) -> Any:
        ligand_info = self.get_rcsb_small_molecule_info(ligand_id=self.lig_id)
        try:
            descriptors = ligand_info["rcsb_chem_comp_descriptor"]
        except KeyError:
            err_msg = f"Failed to get 'rcsb_chem_comp_descriptor' for ligand {self.lig_id}"
            self.logger.error(err_msg)
            raise ValueError(err_msg)
        try:
            inchi = descriptors["in_ch_i"]
        except KeyError:
            err_msg = f"Failed to get 'inchi' from ligand {self.lig_id} descriptors: {descriptors}"
            self.logger.error(err_msg)
            raise ValueError(err_msg)
        assert isinstance(inchi, str), f"Bad InChI: {inchi}"

        try:
            target_mol = Chem.MolFromInchi(inchi)
        except Exception as e:
            err_msg = f"RDKit could not create a molecule from InChi: {inchi}. Got exception: {e}"
            self.logger.error(err_msg)
            raise ValueError(err_msg)

        template_mol = Chem.MolFromPDBFile(str(self.in_pdb), removeHs=True)
        new_mol = self.assign_coordinates_and_names(template_mol, target_mol)
        self.write_mol(new_mol)

    def set_metadata(self, new_mol: Mol, template_mol: Mol, match: list[str]) -> Mol:
        """
        Sets atom names, resnames, etc. for the new molecule based on the template molecule.
        Hydrogens atom names are not copied from the template molecule.

        Args:
            new_mol: The RDKit molecule object to which metadata will be assigned.
            template_mol: The RDKit molecule object from which metadata will be derived.
            match: A list of atom indices that map the template molecule to the new molecule.

        Returns:
            A new RDKit molecule object with the assigned metadata.
        """
        new_mol.SetProp("_Name", self.resname)

        mi = Chem.AtomPDBResidueInfo()
        mi.SetResidueName(self.resname)
        mi.SetResidueNumber(1)
        mi.SetOccupancy(1.0)
        mi.SetTempFactor(0.0)
        mi.SetIsHeteroAtom(True)

        names = [atom.GetPDBResidueInfo().GetName() for atom in template_mol.GetAtoms()]
        for i, atm in enumerate(new_mol.GetAtoms()):
            # mi.SetSerialNumber(i)
            if i < len(match):
                mi.SetName(names[match[i]])
                atm.SetMonomerInfo(mi)
            else:
                mi.SetName("")

                atm.SetMonomerInfo(mi)

        return new_mol

    def write_mol(self, mol: Mol) -> None:
        flavor = 0 if self.add_conect else 2
        self.logger.info(f"Writing {self.in_pdb} to {self.out_pdb}")

        try:
            Chem.MolToPDBFile(mol, self.out_pdb, flavor=flavor)
        except Exception as e:
            self.logger.error(f"Failed to write to  {self.out_pdb}. Got exception: {e}")

    def assign_coordinates_and_names(self, template_mol, target_mol) -> Mol:
        """
        Assigns coordinates from a template molecule to a target molecule using RDKit.

        Args:
            template_mol: The RDKit molecule object with the desired coordinates.
            target_mol: The RDKit molecule object to which coordinates will be assigned.

        Returns:
            A new RDKit molecule object with the assigned coordinates, or None if the
            alignment fails or if either input molecule is None.
        """
        if template_mol is None or target_mol is None:
            err_msg = "Nor `template_mol` or `target_mol` can be None."
            self.logger.error(err_msg)
            raise ValueError(err_msg)

        # Ensure both molecules have the same number of atoms and similar structure
        if template_mol.GetNumAtoms() != target_mol.GetNumAtoms():
            err_msg = f"Error: template ({template_mol.GetNumAtoms()}) and target "
            f"({target_mol.GetNumAtoms()}) molecules must have the same number of atoms."
            self.logger.error(err_msg)
            raise ValueError(err_msg)

        # Use GetSubstructureMatch to get the atom mapping.  Crucial for correct alignment.
        match = template_mol.GetSubstructMatch(target_mol)
        if not match:
            # So I don't get bitten by tautomers
            enumerator = rdMolStandardize.TautomerEnumerator()
            target_mol = enumerator.Canonicalize(target_mol)
            # If input PDB doesn't come with CONECT records, this might help. User may get a warning saying:
            # WARNING: More than one matching pattern found - picking one
            # So we capture it and ignore it.
            f = io.BytesIO()
            with stderr_redirector(f):
                template_mol = AssignBondOrdersFromTemplate(refmol=target_mol, mol=template_mol)

            template_mol = enumerator.Canonicalize(template_mol)

            match = template_mol.GetSubstructMatch(target_mol)
            if not match:
                err_msg = "Error: could not find a substructure match between target and template molecules."
                self.logger.error(err_msg)
                raise ValueError(err_msg)

        # Create a copy of the target molecule *with* the new coordinates.
        new_mol = Chem.Mol(target_mol)

        new_mol.RemoveAllConformers()  # Important: Remove existing conformers.
        new_conf = Chem.Conformer(target_mol.GetNumAtoms())

        for i in range(template_mol.GetNumAtoms()):
            pos = template_mol.GetConformer().GetAtomPosition(match[i])
            new_conf.SetAtomPosition(i, pos)
        new_mol.AddConformer(new_conf)

        hmol = Chem.AddHs(new_mol, addCoords=True)
        hmol = self.set_metadata(hmol, template_mol, match)

        return hmol

    def get_rcsb_small_molecule_info(
            self, ligand_id: str, base_url: str = "https://data.rcsb.org/rest/v1/core/chemcomp/"
    ):
        """
        Queries the RCSB REST API for detailed information about a small molecule (ligand).

        Args:
            ligand_id: The 3-letter RCSB ligand ID (e.g., "NOV").
            base_url: The base URL for the RCSB REST API (default: "https://data.rcsb.org/rest/v1/core/chemcomp/").

        Returns:
            A dictionary containing the ligand information, or None if the query fails.
            Raises a ValueError if `ligand_id` is bad or if it can't find it.
        """
        if not isinstance(ligand_id, str) or len(ligand_id) != 3:
            raise ValueError(f"Bad ligand id ({ligand_id}). Must be a 3-character string.")

        ligand_id = ligand_id.upper()  # Ensure uppercase for consistency
        url = f"{base_url}{ligand_id}"
        try:
            response = requests.get(url)
            response.raise_for_status()  # Raise an exception for bad status codes (4xx or 5xx)
            return response.json()
        except requests.exceptions.RequestException as e:
            err_msg = f"Bad request ({url}). Got exception: {e}"
        except json.JSONDecodeError as e:
            err_msg = f"Bad request ({url}). Error decoding JSON: {e}"
        except Exception as e:
            err_msg = f"Bad request ({url}). Unexpected error occurred: {e}"
        self.logger.error(err_msg)
        raise ValueError(err_msg)

    def _append_stage(self, stage: "AbstractStage") -> "AbstractStage":
        raise NotImplementedError

    def _clean(self):
        raise NotImplementedError

    @staticmethod
    def draw(mol, filepath: Path):
        # Helper method for debugging purposes
        from rdkit.Chem import Draw
        assert filepath.suffix == ".png"
        dm = Draw.PrepareMolForDrawing(mol)
        d2d = Draw.MolDraw2DCairo(450, 400)
        d2d.DrawMolecule(dm)
        d2d.FinishDrawing()
        png = d2d.GetDrawingText()
        with open(filepath, 'wb') as outf:
            outf.write(png)
