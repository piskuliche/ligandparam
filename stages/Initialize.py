
"""
StageInitialize module
----------------------
This module provides the StageInitialize class for initializing a ligand from a PDB file and generating a mol2 file.
"""

from typing import Optional,  Union, Any

from pathlib import Path

from ligandparam.stages.AbstractStage import AbstractStage
from ligandparam.Interfaces import Antechamber
from ligandparam.io.Coordinates import (
    count_structure_atoms,
    mol2_bond_count,
    sanitize_mol2_ligand,
    sanitize_pdb_ligand,
)
import shutil


class StageInitialize(AbstractStage):
    """
    Initialize the ligand from a PDB file and generate a mol2 file.

    Parameters
    ----------
    stage_name : str
        The name of the stage.
    main_input : Union[Path, str]
        Path to the input PDB file.
    cwd : Union[Path, str]
        Current working directory.
    out_mol2 : str
        Path to the output mol2 file.
    net_charge : float, optional
        Net charge for the molecule (default: 0.0).
    assign_charges : bool, optional
        If True (default), run AM1-BCC / ``-c bcc`` (SQM). Gaussian RESP
        recipes set this False so Initialize only assigns GAFF types and
        still writes ``-nc`` / ``-m`` onto the mol2.
    atom_type : str, optional
        Atom type (default: 'gaff2').
    charge_model : str, optional
        Charge model to use ('bcc' or 'abcg2', default: 'bcc').
    sqm : bool, optional
        Whether to run secondary SQM calculation (default: False).
    molname : str, optional
        Molecule name for additional arguments.
    ek : any, optional
        Additional argument for Antechamber.

    Attributes
    ----------
    in_pdb : Path
        Path to the input PDB file.
    out_mol2 : Path
        Path to the output mol2 file.
    net_charge : float
        Net charge for the molecule.
    atom_type : str
        Atom type.
    charge_model : str
        Charge model to use.
    secondary : bool
        Whether to run secondary SQM calculation.
    additional_args : dict
        Additional arguments for Antechamber.
    """
    def __init__(self, stage_name: str, main_input: Union[Path, str], cwd: Union[Path, str], *args, **kwargs) -> None:
        """
        Initialize the StageInitialize class.

        Parameters
        ----------
        stage_name : str
            The name of the stage.
        main_input : Union[Path, str]
            Path to the input PDB file.
        cwd : Union[Path, str]
            Current working directory.
        *args
            Additional positional arguments.
        **kwargs
            Additional keyword arguments.
        """
        super().__init__(stage_name, main_input, cwd, *args, **kwargs)
        self.in_pdb = Path(main_input)
        self.add_required(self.in_pdb)
        self.out_mol2 = Path(kwargs["out_mol2"])

        self.net_charge = int(round(float(kwargs.get("net_charge", 0.0))))
        self.atom_type = kwargs.get("atom_type", "gaff2")
        self.charge_model = kwargs.get("charge_model", "bcc")
        self.multiplicity = int(kwargs.get("multiplicity", 1))
        # Gaussian RESP recipes only need GAFF types here. ``-c bcc`` launches
        # SQM and is what dies when -nc is dropped for anions.
        self.assign_charges = bool(kwargs.get("assign_charges", True))
        self.secondary = kwargs.get("sqm", False)
        if self.assign_charges and self.charge_model not in ("bcc", "abcg2"):
            raise ValueError(f"Unknown charge model '{self.charge_model}'. Must be 'bcc' or 'abcg2'")
        if "molname" in kwargs:
            self.additional_args = {"rn": kwargs["molname"]}
        else:
            self.additional_args = {}
        if "ek" in kwargs:
            self.additional_args["ek"] = kwargs["ek"]

    def _run(self, dry_run=False, nproc: Optional[int] = None, mem: Optional[int] = None) -> Any:
        """
        Execute the initialization stage to generate a mol2 file from a PDB file.

        Parameters
        ----------
        dry_run : bool, optional
            If True, the stage will not be executed, but the function will print the commands that would be run.
        nproc : int, optional
            Number of processors to use.
        mem : int, optional
            Amount of memory to use (in GB).

        Returns
        -------
        None

        Raises
        ------
        ValueError
            If the input file type is not supported.
        """
        detect_type = self.in_pdb.suffix.lower()
        if detect_type not in [".pdb", ".mol2"]:
            raise ValueError(f"Unsupported input file type: {detect_type}. Expected .pdb or .mol2.")
        user_copy = self.cwd / f"{self.in_pdb.stem}.user_input{detect_type}"
        if self.in_pdb.resolve() != user_copy.resolve():
            shutil.copy2(self.in_pdb, user_copy)
        n_user = count_structure_atoms(self.in_pdb)
        if detect_type == ".mol2":
            ftype = "mol2"
            ante_in = self.cwd / f"{self.in_pdb.stem}.antechamber_in.mol2"
            sanitize_mol2_ligand(
                self.in_pdb,
                ante_in,
                net_charge=self.net_charge,
                multiplicity=self.multiplicity,
            )
            n_in = count_structure_atoms(ante_in)
            if n_in != n_user:
                self.logger.warning(
                    "Initialize: dropped %s extra atom(s) from mol2 %s "
                    "(%s -> %s) so charge %s is a closed shell",
                    n_user - n_in,
                    self.in_pdb.name,
                    n_user,
                    n_in,
                    self.net_charge,
                )
        else:
            ftype = "pdb"
            ante_in = self.cwd / f"{self.in_pdb.stem}.sanitized.pdb"
            sanitize_pdb_ligand(
                self.in_pdb,
                ante_in,
                net_charge=self.net_charge,
                multiplicity=self.multiplicity,
            )
            n_in = count_structure_atoms(ante_in)
            self.logger.info(
                "Initialize: wrote Amber-safe PDB %s (%s atoms); "
                "kept CONECT and rewrote element/charge columns",
                ante_in,
                n_in,
            )
        ante = Antechamber(cwd=self.cwd, logger=self.logger, nproc=self.nproc)
        ante_kw = dict(
            i=ante_in,
            fi=ftype,
            o=self.out_mol2,
            fo="mol2",
            nc=self.net_charge,
            m=self.multiplicity,
            pf="y",
            at=self.atom_type,
            an="no",
            dry_run=dry_run,
            **self.additional_args,
        )
        if detect_type == ".mol2" and mol2_bond_count(ante_in) > 0:
            # Keep existing mol2 bonds so antechamber does not rebuild
            # valences and append a sulfate hydrogen.
            ante_kw["j"] = 5
        if self.assign_charges:
            ante_kw["c"] = self.charge_model
        else:
            self.logger.info(
                "Initialize: skipping antechamber -c %s (assign_charges=False); "
                "passing -nc %s -m %s for the mol2 header only",
                self.charge_model,
                self.net_charge,
                self.multiplicity,
            )
        ante.call(**ante_kw)
        if not dry_run:
            n_out = count_structure_atoms(self.out_mol2)
            if n_out != n_in:
                sanitize_mol2_ligand(
                    self.out_mol2,
                    self.out_mol2,
                    net_charge=self.net_charge,
                    multiplicity=self.multiplicity,
                )
                n_out = count_structure_atoms(self.out_mol2)
            if n_out != n_in:
                raise RuntimeError(
                    f"Initialize changed the atom count ({n_in} -> {n_out}) for "
                    f"{self.in_pdb}. Extra hydrogens usually mean the input "
                    "PDB/mol2 used an alcohol oxygen on a sulfate/carboxylate "
                    "anion. Use the sanitized file in the recipe directory and "
                    "keep explicit hydrogens plus net_charge."
                )
        if self.secondary and self.assign_charges:
            second_ante = Antechamber(cwd=self.cwd, logger=self.logger, nproc=self.nproc)
            second_ante.call(
                i="sqm.pdb",
                fi="pdb",
                o=self.out_mol2,
                fo="mol2",
                c=self.charge_model,
                nc=self.net_charge,
                m=self.multiplicity,
                pf="y",
                at=self.atom_type,
                an="no",
                dry_run=dry_run,
                **self.additional_args,
            )
            if not dry_run:
                n_out = count_structure_atoms(self.out_mol2)
                if n_out != n_in:
                    raise RuntimeError(
                        f"SQM follow-up changed the atom count ({n_in} -> {n_out}) "
                        f"for {self.in_pdb}."
                    )

