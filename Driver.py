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

    def _unique_stage_ids(self) -> list[str]:
        """Stable board row ids; suffix duplicates (``Initialize#2``)."""
        seen: dict[str, int] = {}
        ids: list[str] = []
        for stage in self.stages:
            name = str(getattr(stage, "stage_name", type(stage).__name__))
            n = seen.get(name, 0) + 1
            seen[name] = n
            ids.append(name if n == 1 else f"{name}#{n}")
        return ids

    def _start_recipe_board(self):
        """Register stages and start a live stdout board, or return Nones."""
        cwd = getattr(self, "cwd", None)
        if cwd is None:
            return None, None, []
        import sys

        from ligandparam.runtime.ProgressBoard import JobBoardWatcher, JobProgressStore

        cwd = Path(cwd)
        logger = getattr(self, "logger", None)
        title = f"ligandparam recipe - live status ({type(self).__name__})"
        try:
            store = JobProgressStore(
                cwd / ".recipe_progress.json",
                collection_key="stages",
                id_header="Stage",
                title=title,
                empty_hint="no stages registered yet",
                detail_hint_label="Recipe working directory",
            )
            ids = self._unique_stage_ids()
            for stage, job_id in zip(self.stages, ids):
                store.register(
                    job_id,
                    status="queued",
                    stage="queued",
                    detail=type(stage).__name__,
                )
            watcher = JobBoardWatcher(
                store,
                board_path=cwd / "RECIPE_STATUS.txt",
                logger=logger,
                interval_sec=5.0,
                heartbeat_sec=20.0,
                stream=getattr(self, "progress_stream", None) or sys.__stdout__,
                log_root_hint=str(cwd),
                thread_name="recipe-progress-board",
            )
            watcher.start()
        except Exception as exc:
            log = getattr(self, "logger", None)
            if log is not None:
                try:
                    log.warning("Could not start recipe progress board: %s", exc)
                except Exception:
                    pass
            return None, None, []
        return store, watcher, ids

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
        store, watcher, stage_ids = self._start_recipe_board()
        try:
            for index, stage in enumerate(self.stages):
                job_id = stage_ids[index] if index < len(stage_ids) else None
                if store is not None and job_id is not None:
                    store.update(
                        job_id,
                        status="running",
                        stage="running",
                        detail=type(stage).__name__,
                    )
                try:
                    stage.execute(dry_run=dry_run, nproc=nproc, mem=mem)
                except Exception as e:
                    if store is not None and job_id is not None:
                        store.update(
                            job_id,
                            status="failed",
                            stage="failed",
                            error=str(e)[:200],
                            detail=type(stage).__name__,
                        )
                    raise RuntimeError(f"Error in stage {stage.stage_name}: {e}") from e
                if store is not None and job_id is not None:
                    store.update(
                        job_id,
                        status="done",
                        stage="finished",
                        detail=type(stage).__name__,
                    )
        finally:
            if watcher is not None:
                watcher.stop()

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
