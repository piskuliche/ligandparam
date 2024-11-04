from abc import ABCMeta, abstractmethod
from pathlib import Path
from ligandparam.io.coordinates import Coordinates

class AbstractStage(metaclass=ABCMeta):
    """ This is an abstract class for all the stages. """

    default_options = {
        "base_name": None,
        "pdb_filename": None,
        "nproc": 6,
        "mem": "8GB",
        "net_charge": 0.0,
        "theory": {"low": "HF/6-31G*", "high": "PBE1PBE/6-31G*"},
        "atom_type": "gaff2",
        "leaprc": ["leaprc.gaff2"],
        "target_pdb": None,
        "force_gaussian_rerun": False
    }

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
        print("\nFiles generated:")
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

    def _add_required(self, filename):
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
            
        if not hasattr(self, 'required'):
            return
        
        for fname in self.required:
            if not Path(fname).exists():
                raise FileNotFoundError(f"ERROR: File {fname} not found.")
        return 
    
    def _add_output(self, output):
        """ Add the output to the stage. 
        
        Parameters
        ----------
        output : str
            The output file to add to the stage.
        """
        if not hasattr(self, 'outputs'):
            self.outputs = []
        self.outputs.append(output)
        return
    
    def _parse_inputoptions(self, inputoptions=None, **kwargs):
        """ Parse the input options. 
        
        Parameters
        ----------
        inputoptions : dict
            A dictionary of input options.
        **kwargs: dict
            A dictionary of input options
        
        """
        for key, value in self.default_options.items():
            setattr(self, key, value)
        if inputoptions is not None:
            for key in inputoptions:
                if key not in self.default_options:
                    raise KeyError(f"ERROR: {key} is not a valid input option.")
            for key, value in inputoptions.items():
                setattr(self, key, value)
        for key, value in kwargs.items():
            if key not in self.default_options:
                raise KeyError(f"ERROR: {key} is not a valid input option.")
            setattr(self, key, value)
        self._generate_implied()
        self._check_self()
        return
    
    def _generate_implied(self):
        """ Generate the implied options. 
        
        This function generates the implied options, such as the base_name from the pdb_filename.
        
        """
        if self.base_name is None and hasattr(self, 'pdb_filename'):
            self.base_name = Path(self.pdb_filename).stem
        if self.pdb_filename is None and hasattr(self, 'base_name'):
            self.pdb_filename = f"{self.base_name}.pdb"

        self.header =  [f'%NPROC={self.nproc}', f'%MEM={self.mem}']

        try:
            self.coord_object = Coordinates(self.pdb_filename, filetype='pdb')
        except FileExistsError:
            raise FileExistsError(f"ERROR: File {self.pdb_filename} does not exist.")
        return
    
    def _check_self(self):
        pass

    def print_docs(self):
        """ Print the documentation for the stage. """
        try:
            doclines = self.execute.__doc__.split('\n')
            for line in doclines:
                if "Parameters" in line:
                    break
                print(line)
        except:
            print("No documentation available.")
        return
