import os

from antechamber_interface import Antechamber
from gaussianIO import GaussianInput, GaussianWriter, GaussianReader
from coordinates import Coordinates
from driver import Driver
from stageinitialize import StageInitialize
from stagegaussian import StageGaussian
from stagegausrotation import StageGaussianRotation
from stagegaussiantomol2 import StageGaussiantoMol2
from stagenormalizecharges import StageNormalizeCharges


class Parametrization(Driver):
    """An example class that represents a simple entity."""


    def __init__(self, pdb_file, netcharge=None, atom_type=None, theory_low='HF/6-31G*', theory_high='PBE1PBE/6-31G*', nproc=6, mem='8GB', interactive=False):
        """Initialize the class with a PDB file and a net charge.

        Parameters
        ----------
        pdb_file : str, optional
            The path to a PDB file containing the ligand structure.
        netcharge : int, optional
            The net charge of the ligand.
        """
        super().__init__()
        # Inputs
        self.pdb_filename = pdb_file
        self.net_charge = netcharge
        self.atom_type = atom_type
        self.theory={"low":theory_low, 
                     "high":theory_high}
        
        self.nproc = nproc
        self.mem = mem
        self.interactive=interactive

        # Set None behavior
        if self.net_charge is None:
            self.net_charge = 0.0

        if self.atom_type is None:
            self.atom_type = 'gaff2'

        # Set the base name
        self.base_name = self.pdb_filename.strip('.pdb')

        # Things set later
        self.init_gaus_run = None
        self.coord_object = None

        self._generate_header(nproc, mem)
        self.coord_object = self.initial_coordinates()

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

    

class FreeLigand(Parametrization):
    pass

class ProteinLigand(Parametrization):
    def new():
        pass

class RNALigand(Parametrization):
    def new():
        pass
    


if __name__ == "__main__":

    test = FreeLigand('FBKRP.pdb', netcharge=0)
    """
    print(test.pdb_filename)
    print(test.net_charge)
    print(test.atom_type)

    test.initialize()
    test.run_gaussian(dry_run=True)
    """
    #test.add_stage(StageInitialize("Initialize", base_name="FBKRP"))
    #test.add_stage(StageGaussian("Minimize", base_cls=test))
    #test.add_stage(StageGaussianRotation("Rotation1", base_cls=test))
    #test.add_stage(StageGaussiantoMol2("GaussianToMol2", base_cls=test))
    test.add_stage(StageNormalizeCharges("NormalizeCharges", mol2file=test.base_name+"log.mol2", netcharge=0.0))
    test.execute(dry_run=False)

