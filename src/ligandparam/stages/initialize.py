
"""
StageInitialize module
----------------------
This module provides the StageInitialize class for initializing a ligand from a PDB file and generating a mol2 file.
"""

from typing import Optional,  Union, Any

from pathlib import Path

from ligandparam.stages.abstractstage import AbstractStage
from ligandparam.interfaces import Antechamber
from ligandparam.io.coordinates import Remove_PDB_CONECT
from ligandparam.log import get_logger
from rdkit import Chem
from rdkit.Chem import AllChem


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

        self.net_charge = kwargs.get("net_charge", 0.0)
        self.atom_type = kwargs.get("atom_type", "gaff2")
        self.charge_model = kwargs.get("charge_model", "bcc")
        self.secondary = kwargs.get("sqm", False)
        if self.charge_model not in ("bcc", "abcg2"):
            raise ValueError(f"Unknown charge model '{self.charge_model}'. Must be 'bcc' or 'abcg2'")
        if "molname" in kwargs:
            self.additional_args = {"rn": kwargs["molname"]}
        else:
            self.additional_args = {}
        if "ek" in kwargs:
            self.additional_args["ek"] = kwargs["ek"]

    def _append_stage(self, stage: "AbstractStage") -> "AbstractStage":
        """
        Append a stage to the current stage.

        Parameters
        ----------
        stage : AbstractStage
            The stage to append.

        Returns
        -------
        AbstractStage
            The appended stage.
        """
        return stage

    def execute(self, dry_run=False, nproc: Optional[int]=None, mem: Optional[int]=None) -> Any:
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
        super()._setup_execution(dry_run=dry_run, nproc=nproc, mem=mem)
        Remove_PDB_CONECT(self.in_pdb)
        ante = Antechamber(cwd=self.cwd, logger=self.logger, nproc=self.nproc)
        detect_type = self.in_pdb.suffix.lower()
        if detect_type not in [".pdb", ".mol2"]:
            raise ValueError(f"Unsupported input file type: {detect_type}. Expected .pdb or .mol2.")
        if detect_type == ".mol2":
            ftype= "mol2"
        else:
            ftype = "pdb"
        ante.call(
            i=self.in_pdb,
            fi=ftype,
            o=self.out_mol2,
            fo="mol2",
            c=self.charge_model,
            nc=self.net_charge,
            pf="y",
            at=self.atom_type,
            an="no",
            dry_run=dry_run,
            **self.additional_args,
        )
        if self.secondary:
            second_ante = Antechamber(cwd=self.cwd, logger=self.logger, nproc=self.nproc)
            second_ante.call(
                i="sqm.pdb",
                fi="pdb",
                o=self.out_mol2,
                fo="mol2",
                c=self.charge_model,
                nc=self.net_charge,
                pf="y",
                at=self.atom_type,
                an="no",
                dry_run=dry_run,
                **self.additional_args,
            )

    def _clean(self):
        """
        Clean the files generated during the stage.
        """
        raise NotImplementedError("clean method not implemented")


"""
class StageSmilestoPDB(AbstractStage):
     This class is used to initialize from smiles to pdb.
    
    def __init__(self, name,=None) -> None:
        pass
    
    def _append_stage(self, stage: "AbstractStage") -> "AbstractStage":
        pass
    
    def _execute(self, dry_run=False, nproc=1, mem=1):
        pass
    
    def _clean(self):
        pass
    
    """
