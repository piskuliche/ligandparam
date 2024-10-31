from ligandparam.stages.abstractstage import AbstractStage

class TestStage(AbstractStage):
    """ This is an abstract class for all the stages. """
    def __init__(self, name, **kwargs) -> None:
        self.name = name
        self._parse_inputoptions(kwargs)
        return
    
    def _append_stage(self, stage: "AbstractStage") -> "AbstractStage":
        return stage

    def _execute(self, dry_run=False):
        print("This worked!")
        return

    def _clean(self):
        raise NotImplementedError("clean method not implemented")
