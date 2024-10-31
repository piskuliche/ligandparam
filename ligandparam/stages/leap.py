import MDAnalysis as mda

from ligandparam.stages.abstractstage import AbstractStage
from ligandparam.io.leapIO import LeapWriter
from ligandparam.interfaces import Leap

class StageLeap(AbstractStage):
    """ This is class to run a basic leap calculations on the ligand. """
    def __init__(self, name, inputoptions=None) -> None:
        """ Initialize the StageLeap class.
        
        Parameters
        ----------
        name : str
            The name of the stage
     : Ligand
            The base class of the ligand
        """

        self.name = name
        self._parse_inputoptions(inputoptions)
        self.add_required(f"{self.base_name}.frcmod")
        self.add_required(f"{self.base_name}.resp.mol2")

        return
    

    def _append_stage(self, stage: "AbstractStage") -> "AbstractStage":
        """ Appends the stage. """
        return stage
    

    def _execute(self, dry_run=False):
        """ Setup and execute the leap lib file generation """
        # Generate the leap input file
        leapgen = LeapWriter("param")
        # Add the leaprc files
        for rc in self.leaprc:
            leapgen.add_leaprc(rc)

        u = mda.Universe(f"{self.base_name}.resp.mol2")
        resname = u.residues.resnames[0]
        # Add the leap commands
        leapgen.add_line(f"loadamberparams {self.base_name}.frcmod")
        leapgen.add_line(f"{resname} = loadmol2 {self.base_name}.resp.mol2")
        leapgen.add_line(f"saveOff {resname} {self.base_name}.off")
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