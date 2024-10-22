import MDAnalysis as mda

from ligand_param.stages.abstractstage import AbstractStage
from ligand_param.io.leapIO import LeapWriter
from ligand_param.interfaces import Leap

class StageLeap(AbstractStage):
    """ This is class to run a basic leap calculations on the ligand. """
    def __init__(self, name, base_cls=None) -> None:
        """ Initialize the StageLeap class.
        
        Parameters
        ----------
        name : str
            The name of the stage
        base_cls : Ligand
            The base class of the ligand
        """

        self.name = name
        self.base_cls = base_cls
        self.add_required(f"{self.base_cls.base_name}.frcmod")
        self.add_required(f"{self.base_cls.base_name}.resp.mol2")

        return
    

    def _append_stage(self, stage: "AbstractStage") -> "AbstractStage":
        """ Appends the stage. """
        return stage
    

    def _execute(self, dry_run=False):
        """ Setup and execute the leap lib file generation """
        # Generate the leap input file
        leapgen = LeapWriter("param")
        # Add the leaprc files
        for rc in self.base_cls.leaprc:
            leapgen.add_leaprc(rc)

        u = mda.Universe(f"{self.base_cls.base_name}.resp.mol2")
        resname = u.residues.resnames[0]
        # Add the leap commands
        leapgen.add_line(f"loadamberparams {self.base_cls.base_name}.frcmod")
        leapgen.add_line(f"{resname} = loadmol2 {self.base_cls.base_name}.resp.mol2")
        leapgen.add_line(f"saveOff mol {self.base_cls.base_name}.off")
        leapgen.add_line("quit")
        # Write the leap input file
        leapgen.write()
        # Call the leap program
        leap = Leap()
        leap.call(f="tleap.param.in", dry_run = dry_run)
        return
    
    def _clean(self):
        """ Clean the files generated during the stage. """
        raise NotImplementedError("clean method not implemented")