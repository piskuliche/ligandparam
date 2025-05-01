from typing import Optional,  Any

from ligandparam.stages.abstractstage import AbstractStage
from ligandparam.log import get_logger


class TestStage(AbstractStage):
    """ This is an abstract class for all the stages. """
    def __init__(self, name, **kwargs) -> None:
        self.name = name
        self._parse_inputoptions(kwargs)
        return
    
    def _append_stage(self, stage: "AbstractStage") -> "AbstractStage":
        return stage

    def execute(self, dry_run=False, nproc: Optional[int]=None, mem: Optional[int]=None) -> Any:
        self.logger.info("This worked!")
        return

    def _clean(self):
        raise NotImplementedError("clean method not implemented")
