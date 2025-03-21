import MDAnalysis as mda

from ligandparam.abstractstage import AbstractStage
from ligandparam.coordinates import Coordinates
from ligandparam.gaussianIO import GaussianWriter, GaussianInput

from ligandparam.log import get_logger



class StageGaussianRotation(AbstractStage):
    def __init__(self, name, alpha = [0.0], beta = [0.0], gamma = [0.0], base_cls=None) -> None:
        self.name = name
        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma

        if base_cls.coord_object is None:
            raise ValueError(f"Error (Stage {self.name}): Coordinate object not set")

        if base_cls.name is None:
            raise ValueError(f"Error (Stage {self.name}): Base name not set")

        if base_cls.header is None:
            raise ValueError(f"Error (Stage {self.name}): Header not set")

        self.base_cls = base_cls
        
        return
    
    def _append_stage(self, stage: "AbstractStage") -> "AbstractStage":
        return stage

    def execute(self, dry_run=False, nproc=1, mem=1) -> Any:
        self.logger.info(f"Executing {self.stage_name} with alpha={self.alpha}, beta={self.beta}, and gamma={self.gamma}")

        run_apply = print

        for a in self.alpha:
            for b in self.beta:
                for g in self.gamma:
                    #TODO: add elements and header, and make sure they are consistent between steps. Probably initialized with class
                    newgau = GaussianWriter('gaussianCalcs/'+self.base_cls.name+f'_rot_{a}_{b}.com')
                    
                    newgau.add_block(GaussianInput(command=f"#P {self.base_cls.theory['low']} OPT(CalcFC)",
                                        initial_coordinates = self.base_cls.coord_object.rotate(alpha=a, beta=b),
                                        elements = self.base_cls.coord_object.get_elements(),
                                        header=self.base_cls.header))
                    newgau.write(dry_run=dry_run)
                    run_apply(newgau.get_run_command())
        
        return
