from ligand_param.stages.abstractstage import AbstractStage
from ligand_param.io.leapIO import LeapWriter
from ligand_param.interfaces import Leap

class StageLeap(AbstractStage):
    
        def __init__(self, name, base_cls=None, leaprc = []) -> None:
            self.name = name
            self.base_cls = base_cls
            self.leaprc = leaprc

            return
        
    
        def _append_stage(self, stage: "AbstractStage") -> "AbstractStage":
            return stage
        
    
        def execute(self, dry_run=False):
            print(f"Executing {self.name} with netcharge={self.base_cls.net_charge}")
            leapgen = LeapWriter("param")
            for rc in self.leaprc:
                leapgen.add_leaprc(rc)
            leapgen.add_line(f"loadamberparams {self.base_cls.base_name}.frcmod")
            leapgen.add_line(f"mol = loadmol2 {self.base_cls.base_name}.resp.mol2")
            leapgen.add_line(f"saveOff mol {self.base_cls.base_name}.off")
            leapgen.add_line("quit")
            leapgen.write()

            leap = Leap()
            leap.call(f="tleap.param.in", dry_run = dry_run)
            return