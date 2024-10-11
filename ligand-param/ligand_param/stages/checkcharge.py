from ligand_param.abns

class CheckCharge(AbstractStage):
    """ This is an abstract class for all the stages. """
    def __init__(self, name, **kwargs) -> None:
        self.name = name
        pass
    
    def _append_stage(self, stage: "AbstractStage") -> "AbstractStage":
        return stage

    def _execute(self, dry_run=False):
        return

    def _clean(self):
        raise NotImplementedError("clean method not implemented")

