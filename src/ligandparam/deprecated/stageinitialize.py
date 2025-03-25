from typing import Optional,  Any
from ligandparam.abstractstage import AbstractStage
from ligandparam.interfaces import Antechamber


class StageInitialize(AbstractStage):
    """This class is used to initialize from pdb to mol2 file using Antechamber.

    Parameters
    ----------
    name : str
        Name of the stage.
    base_cls : object
        Object of the base class.

    """

    def __init__(self, name, base_cls=None) -> None:
        self.name = name
        self.base_cls = base_cls

        return

    def _append_stage(self, stage: "AbstractStage") -> "AbstractStage":
        return stage

    def execute(self, dry_run=False, nproc=1, mem=1) -> Any:
        ante = Antechamber()
        ante.call(
            i=self.base_cls.name + ".pdb",
            fi="pdb",
            o=self.base_cls.name + ".antechamber.mol2",
            fo="mol2",
            c="bcc",
            nc=self.base_cls.net_charge,
            pf="y",
            at=self.base_cls.atom_type,
            an="no",
            dry_run=dry_run,
        )

