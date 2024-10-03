class GaussianInput:
    """ Class to represent a Gaussian LINK1 block """
    def __init__(self, command="# HF/6-31G* SCF=QC", elements = None, initial_coordinates=None, charge=0, multiplicity=1, **kwargs):
        self.kwargs = kwargs

        if initial_coordinates is None:
            raise ValueError("Initial coordinates must be provided")
        
        if elements is None:
            raise ValueError("Elements must be provided")

        return
    
    def __str__(self):
        print(self.kwargs)
        return
    
    def print_block(self):
        """ Print the Gaussian input block """
        return




if __name__ == "__main__":
    test = GaussianInput(FCCALC=True)
    test.__str__()
