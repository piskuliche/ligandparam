import os

import MDAnalysis as mda

from ligandparam.abstractstage import AbstractStage
from ligandparam.coordinates import Coordinates
from ligandparam.gaussianIO import GaussianWriter, GaussianInput
from ligandparam.interfaces import Gaussian


class StageGaussian(AbstractStage):
    def __init__(self, name, base_cls=None) -> None:
        self.name = name
        self.base_cls = base_cls
        return

    def _append_stage(self, stage: "AbstractStage") -> "AbstractStage":
        return stage

    def execute(self, dry_run=False, nproc=1, mem=1) -> Any:
        stageheader = self.base_cls.header
        stageheader.append(f"%chk={self.base_cls.name}.antechamber.chk")
        gau = GaussianWriter(f'gaussianCalcs/{self.base_cls.name}.com')
        gau.add_block(GaussianInput(command=f"#P {self.base_cls.theory['low']} OPT(CalcFC)",
                                    initial_coordinates=self.base_cls.coord_object.get_coordinates(),
                                    elements=self.base_cls.coord_object.get_elements(),
                                    header=stageheader))
        gau.add_block(GaussianInput(command=f"#P {self.base_cls.theory['high']} OPT(CalcFC) GEOM(ALLCheck) Guess(Read)",
                                    header=stageheader))
        gau.add_block(GaussianInput(
            command=f"#P {self.base_cls.theory['low']} GEOM(AllCheck) Guess(Read) NoSymm Pop=mk IOp(6/33=2) GFInput GFPrint",
            header=stageheader))

        has_run = gau.write(dry_run=dry_run)

        if not has_run:
            gau_run = Gaussian()
            gau_run.call(inp_pipe=self.base_cls.name + '.com',
                         out_pipe=self.base_cls.name + '.log',
                         dry_run=dry_run)
        os.chdir('..')

        return
