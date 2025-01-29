from abc import ABCMeta, abstractmethod
from pathlib import Path
from typing import Union

from ligandparam.io.coordinates import Coordinates
import logging
from ligandparam.log import get_logger




class AbstractStage(metaclass=ABCMeta):
    """ This is an abstract class for all the stages. """

    default_options = {
        "stage_name": None,
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

    def __init__(self, stage_name: str, name: Union[Path, str], cwd: Union[Path, str], *args, **kwargs) -> None:
        self.name = Path(name)
        try:
            self.coord_object = Coordinates(self.name, filetype='pdb')
        except Exception as e:
            raise ValueError(f"ERROR: Invalid input file: {self.name}") from e

        self.stage_name = stage_name
        self.cwd = Path(cwd)
        self.nproc = getattr(kwargs, 'nproc', 12)
        self.mem = getattr(kwargs, 'mem', 48)
        self.header = [f'%NPROC={self.nproc}', f'%MEM={self.mem}GB']
        self.required = []
        self.logger = kwargs.get('logger', get_logger())

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
        self.logger.info(f"Executing {self.stage_name}")
        starting_files = self.list_files_in_directory(".")
        self._check_required()
        self._execute(dry_run=dry_run)
        ending_files = self.list_files_in_directory(".")
        self.new_files = [f for f in ending_files if f not in starting_files]
        # TODO: Write code to ctually assert that the files are there and raise an error if they are not.
        # self.logger.info("\nFiles generated:")
        # for fnames in self.new_files:
        #     self.logger.info(f"------> {fnames}")
        return

    def clean(self) -> None:
        self.logger.info(f"Cleaning {self.stage_name}")
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

    def add_required(self, filename: Union[Path, str]):
        """ Add a required file to the stage. 
        
        Parameters
        ----------
        filename : str
            The file to add to the required list.
        """
        if filename:
            self.required.append(Path(filename))
        return

    def _check_required(self):
        """ Check if the required files are present. """
        for fname in self.required:
            if not Path(fname).exists():
                raise FileNotFoundError(f"ERROR: File {fname} not found.")
        return

    def _add_outputs(self, outputs):
        """ Add the outputs to the stage. 
        
        Parameters
        ----------
        outputs : str
            The output file to add to the stage.
        """
        if not hasattr(self, 'outputs'):
            self.outputs = []
        self.outputs.append(outputs)
        return

    # def _parse_inputoptions(self, inputoptions=None, **kwargs):
    #     """ Parse the input options.
    #
    #     Parameters
    #     ----------
    #     inputoptions : dict
    #         A dictionary of input options.
    #     **kwargs: dict
    #         A dictionary of input options
    #
    #     """
    #     for key, value in self.default_options.items():
    #         setattr(self, key, value)
    #     if inputoptions is not None:
    #         for key in inputoptions:
    #             if key not in self.default_options:
    #                 raise KeyError(f"ERROR: {key} is not a valid input option.")
    #         for key, value in inputoptions.items():
    #             setattr(self, key, value)
    #     for key, value in kwargs.items():
    #         if key not in self.default_options:
    #             raise KeyError(f"ERROR: {key} is not a valid input option.")
    #         setattr(self, key, value)
    #     self._generate_implied()
    #     self._check_self()
    #     return

    def _generate_implied(self):
        """ Generate the implied options. 
        
        This function generates the implied options, such as the name from the pdb_filename.
        
        """

        return

    def _check_self(self):
        pass

    # Quite hacky, but it works.
    def __str__(self) -> str:
        # return str(type(self))
        return str(type(self)).split("'")[1].split(".")[-1]
