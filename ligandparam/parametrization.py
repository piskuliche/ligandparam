import pathlib


from ligandparam.driver import Driver
from ligandparam.io.coordinates import Coordinates
from ligandparam.stages import *

class Parametrization(Driver):
    def __init__(self, inputoptions=None, **kwargs):
        """ This is the base class for all parametrizations, that is a sub class of the :class:`ligand_param.driver.Driver` class.

        The rough approach to using this class is to generate a new Parametrization class, and then generate self.stages as a list 
        of stages that you want to run.

        Parameters
        ----------
        pdb_file : str, optional
            The path to a PDB file containing the ligand structure.
        netcharge : int, optional
            The net charge of the ligand.
        atom_type : str, optional
            The atom type to use for the ligand. Default is 'gaff2'.

        """
        super().__init__()
        self.inputoptions = inputoptions
        if "base_name" in kwargs:
            self.base_name = kwargs['base_name']
        elif "base_name" in inputoptions:
            self.base_name = inputoptions['base_name']
        else:
            raise ValueError("Please provide the base name.")
            
        return
    
    def add_leaprc(self, leaprc):
        self.inputoptions['leaprc'].append(leaprc)
        return
    

class Recipe(Parametrization):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        return