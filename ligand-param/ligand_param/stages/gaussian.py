import os

import MDAnalysis as mda

from ligand_param.stages.abstractstage import AbstractStage
from ligand_param.io.coordinates import Coordinates
from ligand_param.io.gaussianIO import GaussianWriter, GaussianInput
from ligand_param.interfaces import Gaussian, Antechamber



class StageGaussian(AbstractStage):
    """ This is class to run a basic Gaussian calculations on the ligand. 
    
    This does three gaussian steps, one at a low level of theory, one at a higher level of theory, 
    and one for the resp calculation. 
    
    """
    def __init__(self, name, base_cls=None) -> None:
        """ Initialize the StageGaussian class.
        
        Parameters
        ----------
        name : str
            The name of the stage
        base_cls : Ligand
            The base class of the ligand
            
        Returns
        -------
        None
        """
        self.name = name
        self.base_cls = base_cls
        return
    
    def _append_stage(self, stage: "AbstractStage") -> "AbstractStage":
        """Appends the stage.

        Args:
            stage (AbstractStage): _description_

        Returns:
            AbstractStage: _description_
        """
        return stage

    def _execute(self, dry_run=False):
        """ Execute the Gaussian calculations. 

        Parameters
        ----------
        dry_run : bool, optional
            If True, the stage will not be executed, but the function will print the commands that would
        
        Returns
        -------
        None

        """
        stageheader = self.base_cls.header
        stageheader.append(f"%chk={self.base_cls.base_name}.antechamber.chk")

        # Set up the Gaussian Block
        gau = GaussianWriter(f'gaussianCalcs/{self.base_cls.base_name}.com')
        gau.add_block(GaussianInput(command=f"#P {self.base_cls.theory['low']} OPT(CalcFC)",
                                    initial_coordinates = self.base_cls.coord_object.get_coordinates(),
                                    elements = self.base_cls.coord_object.get_elements(),
                                    header=stageheader))
        gau.add_block(GaussianInput(command=f"#P {self.base_cls.theory['high']} OPT(CalcFC) GEOM(ALLCheck) Guess(Read)", 
                                    header=stageheader))
        gau.add_block(GaussianInput(command=f"#P {self.base_cls.theory['low']} GEOM(AllCheck) Guess(Read) NoSymm Pop=mk IOp(6/33=2) GFInput GFPrint", 
                                    header=stageheader))
        
        # Check if the path exists, and make if needed.
        if not os.path.exists(f'gaussianCalcs'):
            os.mkdir('gaussianCalcs')

        gau_complete = False    
        # Check if the Gaussian calculation has already been run
        if os.path.exists(f'gaussianCalcs/{self.base_cls.base_name}.log'):
            gau_complete = True

        # Check if the Gaussian calculation should be rerun
        if self.base_cls.force_gaussian_rerun:
            gau_complete = False

        # Run the Gaussian calculations in the gaussianCalcs directory
        os.chdir('gaussianCalcs')
        if not gau_complete:
            gau.write(dry_run=dry_run)
            gau_run = Gaussian()
            gau_run.call(inp_pipe=self.base_cls.base_name+'.com', 
                         out_pipe=self.base_cls.base_name+'.log',
                         dry_run=dry_run)
        os.chdir('..')

        return
    
    def _clean(self):
        """ Clean the files generated during the stage. """
        raise NotImplementedError("clean method not implemented")
    
class StageGaussianRotation(AbstractStage):
    """ This is class to rotate the ligand and run Gaussian calculations of the resp charges
    for each rotated ligand. """

    def __init__(self, name, alpha = [0.0], beta = [0.0], gamma = [0.0], base_cls=None) -> None:
        self.name = name
        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma

        if base_cls.coord_object is None:
            raise ValueError(f"Error (Stage {self.name}): Coordinate object not set")

        if base_cls.base_name is None:
            raise ValueError(f"Error (Stage {self.name}): Base name not set")

        if base_cls.header is None:
            raise ValueError(f"Error (Stage {self.name}): Header not set")

        self.base_cls = base_cls
        
        return
    
    def _append_stage(self, stage: "AbstractStage") -> "AbstractStage":
        """ Append the stage to the current stage. 
        
        Parameters
        ----------
        stage : AbstractStage
            The stage to append to the current stage
            
        """
        return stage

    def _execute(self, dry_run=False):
        """ Execute the Gaussian calculations for the rotated ligands.
        
        Parameters
        ----------
        dry_run : bool, optional
            If True, the stage will not be executed, but the function will print the commands that would
        
        Returns
        -------
        """

        run_apply = print

        for a in self.alpha:
            for b in self.beta:
                for g in self.gamma:
                    #TODO: add elements and header, and make sure they are consistent between steps. Probably initialized with class
                    newgau = GaussianWriter('gaussianCalcs/'+self.base_cls.base_name+f'_rot_{a}_{b}.com')
                    
                    newgau.add_block(GaussianInput(command=f"#P {self.base_cls.theory['low']} OPT(CalcFC)",
                                        initial_coordinates = self.base_cls.coord_object.rotate(alpha=a, beta=b),
                                        elements = self.base_cls.coord_object.get_elements(),
                                        header=self.base_cls.header))
                    newgau.write(dry_run=dry_run)
                    run_apply(newgau.get_run_command())
        
        return
    
    def _clean(self):
        return
    
class StageGaussiantoMol2(AbstractStage):
    """ Convert Gaussian output to mol2 format. 
    
    This class converts the Gaussian output to mol2 format, and assigns the charges to the mol2 file. """
    def __init__(self, name, base_cls=None, dry_run = None) -> None:
        """ Initialize the StageGaussiantoMol2 class 
        
        Parameters
        ----------
        name : str
            The name of the stage
        base_cls : Ligand
            The base class of the ligand
        dry_run : bool, optional
            If True, the stage will not be executed, but the function will print the commands that would
        
        Returns
        -------
        None
        
        """
        self.name = name
        self.base_cls = base_cls
        self.dry_run = dry_run

    def _append_stage(self, stage: "AbstractStage") -> "AbstractStage":
        """ Append the stage to the current stage. """
        return stage

    def _execute(self, dry_run=False):
        """ Execute the Gaussian to mol2 conversion.

        Parameters
        ----------
        dry_run : bool, optional
            If True, the stage will not be executed, but the function will print the commands that would
        
        Returns
        -------
        None

        """

        if self.dry_run is not None:
            dry_run = self

        # Convert from gaussian to mol2
        ante = Antechamber()
        ante.call(i=self.base_cls.base_name+'.log', fi='gout',
                  o=self.base_cls.base_name+'.tmp1.mol2', fo='mol2',
                  pf='y', at=self.base_cls.atom_type,
                  run=(not dry_run))

        # Assign the charges 
        if not dry_run:
            u1 = mda.Universe(self.base_cls.base_name+'.antechamber.mol2')
            u2 = mda.Universe(self.base_cls.base_name+'.tmp1.mol2')
            assert len(u1.atoms) == len(u2.atoms), "Number of atoms in the two files do not match"

            u2.atoms.charges = u1.atoms.charges

            ag = u2.select_atoms("all")
            ag.write(self.base_cls.base_name+'.tmp2.mol2')

        # Use antechamber to clean up the mol2 format
        ante = Antechamber()
        ante.call(i=self.base_cls.base_name+'.tmp2.mol2', fi='mol2',
                  o=self.base_cls.base_name+'.log.mol2', fo='mol2',
                  pf='y', at=self.base_cls.atom_type,
                  run=(not dry_run))
        
        return
    
    def _clean(self):
        return
