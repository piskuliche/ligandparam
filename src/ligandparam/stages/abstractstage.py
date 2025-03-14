from abc import ABCMeta, abstractmethod
from pathlib import Path
from typing import Union, Any

from ligandparam.io.coordinates import Coordinates
from ligandparam.log import get_logger
import warnings


class AbstractStage(metaclass=ABCMeta):
    """This is an abstract class for all the stages."""

    def __init__(self, stage_name: str, main_input: Union[Path, str], cwd: Union[Path, str], *args, **kwargs) -> None:
        # TODO Fix: we assume that all stages deal with an input file, but don't read it yet. Make in_filename a kwarg.
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                self.coord_object = Coordinates(main_input, filetype="pdb")
        except Exception:
            # TODO: Fix Not a pdb, no problem. This shouldn't be in this class.
            pass

        resname = kwargs.get("resname", "LIG")
        if resname and len(resname) > 3:
            raise ValueError(f"Bad input resname: {kwargs['resname']}")

        self.cwd = Path(cwd)
        if not self.cwd.parent.is_dir():
            raise ValueError(f"Bad input `cwd` working dir: {self.cwd}")

        self.stage_name = stage_name
        self.required = []
        self.logger = kwargs.get("logger", get_logger())
        self.nproc = kwargs.get("nproc", 1)
        self.mem = kwargs.get("mem", 512)
        self.dry_run = kwargs.get("dry_run", False)

    @abstractmethod
    def _append_stage(self, stage: "AbstractStage") -> "AbstractStage":
        pass

    @abstractmethod
    def _clean(self):
        pass

    def append_stage(self, stage: "AbstractStage") -> "AbstractStage":
        return self._append_stage(stage)

    def execute(self, dry_run=False, nproc=1, mem=512) -> Any:
        self.logger.info(f"Executing {self.stage_name}")
        starting_files = self.list_files_in_directory(self.cwd)
        self._check_required()
        self.execute(self, dry_run=False, nproc=1, mem=512)
        ending_files = self.list_files_in_directory(self.cwd)
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
        """List all the files in a directory.

        Parameters
        ----------
        directory : str
            The directory to list the files from.

        """
        return [f.name for f in Path(directory).iterdir() if f.is_file()]

    def add_required(self, filename: Union[Path, str]):
        """Add a required file to the stage.

        Parameters
        ----------
        filename : str
            The file to add to the required list.
        """
        if filename:
            self.required.append(Path(filename))
        return

    def _check_required(self):
        """Check if the required files are present."""
        for fname in self.required:
            if not Path(fname).exists():
                raise FileNotFoundError(f"ERROR: File {fname} not found.")
        return

    def _add_outputs(self, outputs):
        """Add the outputs to the stage.

        Parameters
        ----------
        outputs : str
            The output file to add to the stage.
        """
        if not hasattr(self, "outputs"):
            self.outputs = []
        self.outputs.append(outputs)
        return

    def _generate_implied(self):
        """Generate the implied options.

        This function generates the implied options, such as the name from the pdb_filename.

        """

        return

    def _check_self(self):
        pass

    # Quite hacky, but it works.
    def __str__(self) -> str:
        # return str(type(self))
        return str(type(self)).split("'")[1].split(".")[-1]
