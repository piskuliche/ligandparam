import os


class Parametrization:
    """An example class that represents a simple entity."""


    def __init__(self, pdb_file, netcharge=None, atom_type=None, theory_low='HF/6-31G*', theory_high='PBE1PBE/6-31G*'):
        """Initialize the class with a PDB file and a net charge.

        Parameters
        ----------
        pdb_file : str, optional
            The path to a PDB file containing the ligand structure.
        netcharge : int, optional
            The net charge of the ligand.
        """
        # Inputs
        self.pdb_filename = pdb_file
        self.net_charge = netcharge
        self.atom_type = atom_type
        self.theory={"low":theory_low, 
                     "high":theory_high}

        # Set None behavior
        if self.net_charge is None:
            self.net_charge = 0.0

        if self.atom_type is None:
            self.atom_type = 'gaff2'

        # Set the base name
        self.base_name = self.pdb_filename.strip('.pdb')

        return
    
    def parametrize(self):
        """Parametrize the ligand."""
        self.gaussian_init()
    
    def gaussian_init(self):
        """Initialize the Gaussian calculations for the ligand.
        
        This method will generate the input files required to run Gaussian calculations
        on the ligand. The files will be generated in the current working directory.
        
        Parameters
        ----------
        None

        Returns
        -------
        None

        """
        # Create a mol2 file
        self._call_antechamber(i = self.pdb_filename,
                               fi = 'pdb',
                               o = f'{self.base_name}.mol2',
                               fo = 'mol2',
                               c = 'bcc',
                               nc = self.net_charge,
                               pf = 'y',
                               at = self.atom_type,
                               run = True)

        # Create a Gaussian Input File
        self._create_gaussian_input()

        return
    
    def _call_antechamber(self, **kwargs):
        """
        Call the Antechamber program to generate force field parameters for the ligand.
        
        Parameters
        ----------
        **kwargs
            Arbitrary keyword arguments to pass to Antechamber. Arguments should be specified
            using the same names as in the Antechamber manual. e.g. i='input.pdb'

        Returns
        -------
        None
            
        """
        import subprocess
        run=False
        if "run" in kwargs:
            run = kwargs['run']
            del kwargs['run']
        
        command = ['antechamber']
        for key, value in kwargs.items():
            if value is not None:
                command.extend([f'-{key}', str(value)])
        if run:
            subprocess.run(command)
        else:
            print(command)
        return
    
    def _create_gaussian_input(self):
        
    
    
    
    


if __name__ == "__main__":

    test = Parametrization('FBKRP.pdb', netcharge=0)

    print(test.pdb_filename)
    print(test.net_charge)
    print(test.atom_type)

    test.parametrize()

