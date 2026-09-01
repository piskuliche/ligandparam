from abc import abstractmethod
from typing import Optional, Union, Any
from pathlib import Path


class Driver:
    """Base class for ordered stage pipelines.

    Subclasses populate :attr:`stages` and call :meth:`execute` to run them
    in order. Recipes such as :class:`~ligandparam.Parametrization.Recipe`
    inherit from this class.

    Attributes
    ----------
    in_filename : Path
        Primary input structure path.
    cwd : Path
        Working directory for intermediate and output files.
    stages : list
        Ordered stage objects to execute.
    """

    @abstractmethod
    def __init__(self, in_filename: Union[Path, str], cwd: Union[Path, str], *args, **kwargs):
        """Initialize the driver. Subclasses must set ``stages`` (often empty)."""
        pass

    def add_stage(self, stage):
        """Append a stage to the pipeline and list the current stages.

        Parameters
        ----------
        stage : AbstractStage
            Stage instance to add. Stages should already be ordered for the
            intended workflow.
        """
        self.stages.append(stage)
        self.list_stages()

    def execute(self, dry_run=False, nproc: Optional[int] = None, mem: Optional[int] = None) -> Any:
        """Run each stage in order.

        Parameters
        ----------
        dry_run : bool, optional
            If True, stages should log planned work without running external tools.
        nproc : int, optional
            Processor count override for stages that support it.
        mem : int, optional
            Memory override in GB for stages that support it.

        Raises
        ------
        RuntimeError
            If a stage raises an exception during execution.
        """
        for stage in self.stages:
            try:
                stage.execute(dry_run=dry_run, nproc=nproc, mem=mem)
            except Exception as e:
                raise RuntimeError(f"Error in stage {stage.stage_name}: {e}") from e

    def clean(self):
        """Clean stage outputs in reverse order.

        Stages without a ``clean`` implementation are skipped. Other failures
        are re-raised after logging.
        """
        for stage in reversed(self.stages):
            try:
                stage.clean()
            except NotImplementedError:
                print(f"Clean method not implemented for stage {stage.stage_name}")
                print("Skipping...")
                continue
            except Exception as e:
                print(f"Error in stage {stage.stage_name}: {e}")
                print("Exiting")
                raise

    def list_stages(self):
        """Print the current stage list to stdout."""
        print("List of Stages to Run")
        for stage in self.stages:
            print(f"-->{stage.stage_name} ({stage})")

    def remove_stage(self, stage_name):
        """Remove the first stage whose ``stage_name`` matches.

        Parameters
        ----------
        stage_name : str
            Name of the stage to remove.
        """
        for stage in self.stages:
            if stage.stage_name == stage_name:
                self.stages.remove(stage)
                print(f"Stage {stage_name} removed.")
                self.list_stages()
                return
        print(f"Stage {stage_name} not found in list of stages.")

    def insert_stage(self, newstage, stage_name, print_info=False):
        """Insert ``newstage`` immediately before the named stage.

        Parameters
        ----------
        newstage : AbstractStage
            Stage to insert.
        stage_name : str
            Existing stage name to insert before.
        print_info : bool, optional
            If True, print the updated stage list.

        Raises
        ------
        ValueError
            If ``stage_name`` is not found.
        """
        for stage in self.stages:
            if stage.stage_name == stage_name:
                idx = self.stages.index(stage)
                self.stages.insert(idx, newstage)
                if print_info:
                    print(f"Stage {newstage.stage_name} inserted before {stage_name}")
                    self.list_stages()
                return
        raise ValueError(f"Stage {stage_name} not found in list of stages.")
