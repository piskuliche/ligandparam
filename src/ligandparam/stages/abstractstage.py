from abc import ABCMeta, abstractmethod
from pathlib import Path
from typing import Optional, Union, Any

from ligandparam.io.coordinates import Coordinates
from ligandparam.log import get_logger
import warnings


class AbstractStage(metaclass=ABCMeta):
    """
    Abstract base class for all stages.

    Parameters
    ----------
    stage_name : str
        The name of the stage.
    main_input : Union[Path, str]
        The main input file for the stage.
    cwd : Union[Path, str]
        The current working directory for the stage.

    Attributes
    ----------
    stage_name : str
        The name of the stage.
    main_input : Path
        The main input file for the stage.
    cwd : Path
        The current working directory for the stage.
    required : list
        A list of required files for the stage.
    logger : logging.Logger
        The logger instance for the stage.
    nproc : int
        The number of processors to use.
    mem : int
        The amount of memory to use (in GB).
    dry_run : bool
        Whether to execute the stage in dry-run mode.
    """

    def __init__(self, stage_name: str, main_input: Union[Path, str], cwd: Union[Path, str], *args, **kwargs) -> None:
        """
        Initialize the AbstractStage.

        Parameters
        ----------
        stage_name : str
            The name of the stage.
        main_input : Union[Path, str]
            The main input file for the stage.
        cwd : Union[Path, str]
            The current working directory for the stage.
        """
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                self.coord_object = Coordinates(main_input, filetype="pdb")
        except Exception:
            pass

        resname = kwargs.get("resname", "LIG")
        if resname and len(resname) > 3:
            raise ValueError(f"Bad input resname: {kwargs['resname']}")

        self.cwd = Path(cwd)
        if not self.cwd.parent.is_dir():
            raise ValueError(f"Bad input `cwd` working dir: {self.cwd}")

        self.main_input = Path(main_input).resolve()
        self.stage_name = stage_name
        self.required = []
        self.logger = kwargs.get("logger", get_logger())
        self.nproc = kwargs.get("nproc", 1)
        self.mem = kwargs.get("mem", 1)
        self.dry_run = kwargs.get("dry_run", False)

    @abstractmethod
    def _append_stage(self, stage: "AbstractStage") -> "AbstractStage":
        """
        Append a stage to the current stage.

        Parameters
        ----------
        stage : AbstractStage
            The stage to append.

        Returns
        -------
        AbstractStage
            The appended stage.
        """
        pass

    @abstractmethod
    def _clean(self):
        """
        Clean up the stage.

        This method performs any necessary cleanup after the stage has been executed.
        """
        pass

    def append_stage(self, stage: "AbstractStage") -> "AbstractStage":
        """
        Append a stage to the current stage.

        Parameters
        ----------
        stage : AbstractStage
            The stage to append.

        Returns
        -------
        AbstractStage
            The appended stage.
        """
        return self._append_stage(stage)

    def _setup_execution(self, dry_run=False, nproc: Optional[int]=None, mem: Optional[int]=None) -> None:
        """
        Set up the execution environment for the stage.

        Parameters
        ----------
        dry_run : bool, optional
            Whether to execute in dry-run mode (default is False).
        nproc : int, optional
            The number of processors to use (default is None).
        mem : int, optional
            The amount of memory to use (in GB, default is None).
        """
        self.nproc = self.nproc if nproc is None else nproc
        self.mem = self.mem if mem is None else mem
        self._check_required()

    def execute(self, dry_run=False, nproc: Optional[int] = None, mem: Optional[int] = None) -> Any:
        """
        Execute the stage.

        Parameters
        ----------
        dry_run : bool, optional
            Whether to execute in dry-run mode (default is False).
        nproc : int, optional
            The number of processors to use (default is None).
        mem : int, optional
            The amount of memory to use (in GB, default is None).

        Returns
        -------
        Any
            The result of the execution.
        """
        self.logger.info(f"Executing {self.stage_name}")
        starting_files = self.list_files_in_directory(self.cwd)
        self._check_required()

        self._setup_execution(dry_run=dry_run, nproc=nproc, mem=mem)
        self.execute(self, nproc=self.nproc, mem=self.mem)
        ending_files = self.list_files_in_directory(self.cwd)
        self.new_files = [f for f in ending_files if f not in starting_files]
        return

    def clean(self) -> None:
        """
        Clean up the stage.

        This method performs any necessary cleanup after the stage has been executed.
        """
        self.logger.info(f"Cleaning {self.stage_name}")
        self._clean()
        return

    def list_files_in_directory(self, directory):
        """
        List all the files in a directory.

        Parameters
        ----------
        directory : str
            The directory to list the files from.

        Returns
        -------
        list of str
            A list of file names in the directory.
        """
        return [f.name for f in Path(directory).iterdir() if f.is_file()]

    def add_required(self, filename: Union[Path, str]):
        """
        Add a required file to the stage.

        Parameters
        ----------
        filename : str
            The file to add to the required list.
        """
        if filename:
            self.required.append(Path(filename))
        return

    def _check_required(self):
        """
        Check if the required files are present.

        Raises
        ------
        FileNotFoundError
            If any of the required files are not found.
        """
        for fname in self.required:
            if not Path(fname).exists():
                raise FileNotFoundError(f"ERROR: File {fname} not found.")
        return

    def _add_outputs(self, outputs):
        """
        Add the outputs to the stage.

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
        """
        Generate the implied options.

        This function generates the implied options, such as the name from the pdb_filename.
        """
        return

    def _check_self(self):
        """
        Perform self-checks for the stage.

        This method ensures that the stage is properly configured and ready for execution.
        """
        pass

    def __str__(self) -> str:
        """
        Return a string representation of the stage.

        Returns
        -------
        str
            The string representation of the stage.
        """
        return str(type(self)).split("'")[1].split(".")[-1]
