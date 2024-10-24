import pathlib


from ligand_param.driver import Driver
from ligand_param.io.coordinates import Coordinates
from ligand_param.stages import *




class Parametrization(Driver):
    def __init__(self, pdb_file, netcharge=0.0, atom_type="gaff2", 
                 theory_low='HF/6-31G*', theory_high='PBE1PBE/6-31G*', 
                 nproc=6, mem='8GB',
                leaprc = [], target_pdb=None,
                force_gaussian_rerun=False):
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
        # Inputs
        self.pdb_filename = pdb_file
        self.net_charge = netcharge
        self.atom_type = atom_type
        self.theory={"low":theory_low, 
                     "high":theory_high}
        self.force_gaussian_rerun = force_gaussian_rerun
        self.target_pdb = target_pdb


        # Settings for Processor Environment
        self.nproc = nproc
        self.mem = mem

        self.leaprc = []
        if not leaprc:
            self.leaprc = ['leaprc.gaff2']
        else:
            self.leaprc = leaprc


        # Set the base name
        self.base_name = pathlib.Path(self.pdb_filename).stem

        # Generate the Coordinates
        self._generate_header(nproc, mem)
        self.coord_object = self.initial_coordinates()

        # Print out information for the user.
        self.print_info()

        return
    
    def add_leaprc(self, leaprc):
        self.leaprc.append(leaprc)
        return
    
    def initial_coordinates(self):
        try:
            coord_object = Coordinates(self.pdb_filename, filetype='pdb')
            return coord_object
        except FileExistsError:
            raise FileExistsError(f"ERROR: File {self.pdb_filename} does not exist.")
        
    def _generate_header(self, nproc, mem):
        self.header = [f'%NPROC={nproc}', f'%MEM={mem}']
        return
    
    def print_info(self):
        print(f"New Parametrization for {self.pdb_filename}")
        print("*******************************")
        print("User supplied parameters:")
        print(f"Net charge: {self.net_charge}")
        print(f"Atom type: {self.atom_type}")
        print(f"Low level QM Theory: {self.theory['low']}")
        print(f"High level QM Theory: {self.theory['high']}")
        print(f"Number of processors: {self.nproc}")
        print(f"Memory: {self.mem}")
        print("Gaussian header is:")
        print(self.header)
        if self.force_gaussian_rerun:
            print("Forcing Gaussian calculations to rerun.")
        else:
            print("Defaulting to NOT rerunning Gaussian calculations.")

class Recipe(Parametrization):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        return