from abc import ABCMeta, abstractmethod

class AbstractStage(metaclass=ABCMeta):
    """ This is an abstract class for all the stages. """
    @abstractmethod
    def __init__(self, name, **kwargs) -> None:
        pass
    
    @abstractmethod
    def _append_stage(self, stage: "AbstractStage") -> "AbstractStage":
        pass

    @abstractmethod
    def execute(self, dry_run=False):
        pass

    def append_stage(self, stage: "AbstractStage") -> "AbstractStage":
        return self._append_stage(stage)
    
    def _execute(self, dry_run=False) -> None:
        print("************************************")
        print(f"Executing {self.name}")
        print("************************************")
        self.execute(dry_run=dry_run)
