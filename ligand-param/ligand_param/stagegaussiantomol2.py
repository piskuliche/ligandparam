from abstractstage import AbstractStage
from antechamber_interface import Antechamber

class StageGaussiantoMol2(AbstractStage):
    def __init__(self, name, base_cls=None, dry_run = None) -> None:
        self.name = name
        self.base_cls = base_cls
        self.dry_run = dry_run

    def _append_stage(self, stage: "AbstractStage") -> "AbstractStage":
        return stage

    def execute(self, dry_run=False):

        if self.dry_run is not None:
            dry_run = self

        ante = Antechamber()
        ante.call(i=self.base_cls.base_name+'.g09', fi='g09',
                  o=self.base_cls.base_name+'.mol2', fo='mol2',
                  pf='y', at=self.base_cls.atom_type,
                  run=(not dry_run))

        # Assign the charges

        return

