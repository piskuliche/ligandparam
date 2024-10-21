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
        starting_files = self.list_files_in_directory(".")
        self._check_required()
        self._execute(dry_run=dry_run)
        ending_files = self.list_files_in_directory(".")
        self.new_files = [f for f in ending_files if f not in starting_files]
        print("Files generated:")
        for fnames in self.new_files:
            print(f"------> {fnames}")
        return

    def clean(self) -> None:
        print("************************************")
        print(f"Cleaning {self.name}")
        print("************************************")
        self._clean()
        return
    
    def list_files_in_directory(self, directory):
        """ List all the files in a directory. 
        
        Parameters
        ----------
        directory : str
            The directory to list the files from.
            
        """
        return [f.name for f in Path(directory).iterdir() if f.is_file()]

    def add_required(self, filename):
        """ Add a required file to the stage. 
        
        Parameters
        ----------
        filename : str
            The file to add to the required list.
        """
        if not hasattr(self, 'required'):
            self.required = []
        self.required.append(filename)
        return
    
    def _check_required(self):
        """ Check if the required files are present. """
        current_fnames = self.list_files_in_directory(".")
        if not hasattr(self, 'required'):
            return 
        for fname in self.required:
            if fname not in current_fnames:
                raise FileNotFoundError(f"ERROR: File {fname} not found.")
        return 

