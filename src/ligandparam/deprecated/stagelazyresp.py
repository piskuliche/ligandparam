import MDAnalysis as mda
import numpy as np

from ligandparam.abstractstage import AbstractStage
from ligandparam.interfaces import Antechamber
from ligandparam.log import get_logger


class StageLazyResp(AbstractStage):
    def __init__(self, name, base_cls=None) -> None:
        self.name = name
        self.base_cls = base_cls
        return

    def _append_stage(self, stage: "AbstractStage") -> "AbstractStage":
        return stage

    def execute(self, dry_run=False, nproc=1, mem=1) -> Any:
        self.logger.info(f"Executing {self.name} with netcharge={self.base_cls.net_charge}")
        ante = Antechamber()
        ante.call(
            i=f"gaussianCalcs/{self.base_cls.name}.log",
            fi="gout",
            o=self.base_cls.name + ".resp.mol2",
            fo="mol2",
            gv=0,
            c="resp",
            nc=self.base_cls.net_charge,
            at="gaff2",
            an="no",
            dry_run=dry_run,
        )
        return
