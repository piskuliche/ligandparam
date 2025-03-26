import warnings
from typing import Optional,  Union, Any

import numpy as np
import MDAnalysis as mda

from pathlib import Path

from ligandparam.stages.abstractstage import AbstractStage
from ligandparam.interfaces import Antechamber
from ligandparam.io.coordinates import Mol2Writer
from ligandparam.log import get_logger


class StageDisplaceMol(AbstractStage):
    """Displaces a molecule based on a vector or centers it at the origin.

    Attributes:
        in_molecule (Path): Path to the input molecule file.
        out_molecule (Path): Path where the displaced/centered molecule will be written.
        displacement_vtor (Optional[np.ndarray]): The vector used for displacement.
                                                 Only set if 'vector' is provided in kwargs.
        center (bool): If True, the molecule will be centered at the origin.
                       If False, displacement is done using `displacement_vtor`.
    """

    def __init__(self, stage_name: str, main_input: Union[Path, str], cwd: Union[Path, str], *args, **kwargs) -> None:
        super().__init__(stage_name, main_input, cwd, *args, **kwargs)
        self.in_molecule = Path(main_input)
        self.out_molecule = Path(kwargs["out_mol"])

        if "vector" in kwargs:
            self.displacement_vtor = kwargs["vector"]
            if not isinstance(self.displacement_vtor, np.ndarray):
                raise TypeError("vector must be a numpy array")
            self.center = False
        else:
            self.center = True
        self.add_required(Path(self.in_molecule))

    def _append_stage(self, stage: "AbstractStage") -> "AbstractStage":
        return stage

    def execute(self, dry_run=False, nproc: Optional[int]=None, mem: Optional[int]=None) -> np.ndarray:

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            u = mda.Universe(self.in_molecule)
            if self.center:
                self.displacement_vtor = -u.atoms.center_of_mass()
            u.atoms.translate(self.displacement_vtor)
            u.atoms.write(self.out_molecule)
        return self.displacement_vtor


    def _clean(self):
        raise NotImplementedError("clean method not implemented")
