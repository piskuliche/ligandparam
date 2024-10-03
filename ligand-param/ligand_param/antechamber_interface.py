class Antechamber:
    def __init__(self):
        print("Antechamber object created")
        return
    def call(self, **kwargs):
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