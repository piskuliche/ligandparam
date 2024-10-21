from abc import ABCMeta, abstractmethod
from pathlib import Path

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
        self.starting_files = self.list_files_in_directory(".")
        self._execute(dry_run=dry_run)
        self.ending_files = self.list_files_in_directory(".")

    def clean(self) -> None:
        print("************************************")
        print(f"Cleaning {self.name}")
        print("************************************")
        self._clean()
        return

    def _add_output(self, file):
        """ Add an output file to the stage. 

        This is used by different stages to keep track of the output files that
        it generates. This allows for easy cleaning of the files.
        
        Parameters
        ----------
        file : str
            The file to add to the output list.
        """
        if self.outputs is None:
            self.outputs = []
        self.outputs.append(file)
        return
    
    def list_files_in_directory(self, directory):
        """ List all the files in a directory. 
        
        Parameters
        ----------
        directory : str
            The directory to list the files from.
            
        """
        return [f.name for f in Path(directory).iterdir() if f.is_file()]
    
