from ligand_param.stages.abstractstage import AbstractStage

class TestStage(AbstractStage):
    """ This is an abstract class for all the stages. """
    def __init__(self, name, **kwargs) -> None:
        self.name = name
        return
    
    def _append_stage(self, stage: "AbstractStage") -> "AbstractStage":
        return stage

    def execute(self, dry_run=False):
        print("This worked!")
        return
