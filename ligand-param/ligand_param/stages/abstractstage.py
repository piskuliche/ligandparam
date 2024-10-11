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
    def _execute(self, dry_run=False):
        pass

    @abstractmethod
    def _clean(self):
        pass

    def append_stage(self, stage: "AbstractStage") -> "AbstractStage":
        return self._append_stage(stage)
    
    def execute(self, dry_run=False) -> None:
        print("************************************")
        print(f"Executing {self.name}")
        print("************************************")
        self._execute(dry_run=dry_run)

    def clean(self) -> None:
        print("************************************")
        print(f"Cleaning {self.name}")
        print("************************************")
        self._clean()
        return
