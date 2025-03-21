from typing import Optional,  Union, Any

from pathlib import Path

from ligandparam.stages.abstractstage import AbstractStage
from ligandparam.interfaces import Antechamber
from ligandparam.io.coordinates import Remove_PDB_CONECT
from ligandparam.log import get_logger


class StageInitialize(AbstractStage):
    def __init__(self, stage_name: str, main_input: Union[Path, str], cwd: Union[Path, str], *args, **kwargs) -> None:
        super().__init__(stage_name, main_input, cwd, *args, **kwargs)
        self.in_pdb = Path(main_input)
        self.add_required(self.in_pdb)
        self.out_mol2 = Path(kwargs["out_mol2"])

        self.net_charge = kwargs.get("net_charge", 0.0)
        self.atom_type = kwargs.get("atom_type", "gaff2")
        self.charge_model = kwargs.get("charge_model", "bcc")
        if self.charge_model not in ("bcc", "abcg2"):
            raise ValueError(f"Unknown charge model '{self.charge_model}'. Must be 'bcc' or 'abcg2'")
        if "molname" in kwargs:
            self.additional_args = {"rn": kwargs["molname"]}
        else:
            self.additional_args = {}

    def _append_stage(self, stage: "AbstractStage") -> "AbstractStage":
        """Appends the stage."""
        return stage

    def execute(self, dry_run=False, nproc: Optional[int]=None, mem: Optional[int]=None) -> Any:
        """Execute the Gaussian calculations.

        Parameters
        ----------
        dry_run : bool, optional
            If True, the stage will not be executed, but the function will print the commands that would

        Returns
        -------
        None
        """
        super()._setup_execution(dry_run=dry_run, nproc=nproc, mem=mem)
        Remove_PDB_CONECT(self.in_pdb)
        ante = Antechamber(cwd=self.cwd, logger=self.logger, nproc=self.nproc)
        ante.call(
            i=self.in_pdb,
            fi="pdb",
            o=self.out_mol2,
            fo="mol2",
            c=self.charge_model,
            nc=self.net_charge,
            pf="y",
            at=self.atom_type,
            dry_run=dry_run,
            **self.additional_args,
        )

    def _clean(self):
        """Clean the files generated during the stage."""
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
