import os
import MDAnalysis as mda

from pathlib import Path

from ligand_param.stages.abstractstage import AbstractStage
from ligand_param.io.coordinates import Coordinates, SimpleXYZ, Mol2Writer
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

        # Set up the Gaussian Block - it does not yet write anything,
        # so this part can be set up before the Gaussian calculations are run.
        gau = GaussianWriter(f'gaussianCalcs/{self.base_cls.base_name}.com')
        gau.add_block(GaussianInput(command=f"#P {self.base_cls.theory['low']} OPT(CalcFC)",
                                    initial_coordinates = self.base_cls.coord_object.get_coordinates(),
                                    elements = self.base_cls.coord_object.get_elements(),
                                    charge = self.base_cls.net_charge,
                                    header=stageheader))
        gau.add_block(GaussianInput(command=f"#P {self.base_cls.theory['high']} OPT(CalcFC) GEOM(ALLCheck) Guess(Read)", 
                                    charge=self.base_cls.net_charge,
                                    header=stageheader))
        gau.add_block(GaussianInput(command=f"#P {self.base_cls.theory['low']} GEOM(AllCheck) Guess(Read) NoSymm Pop=mk IOp(6/33=2) GFInput GFPrint", 
                                    charge=self.base_cls.net_charge,
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
        
        if not gau_complete:
            gau.write(dry_run=dry_run)

        # Run the Gaussian calculations in the gaussianCalcs directory
        os.chdir('gaussianCalcs')
        if not gau_complete:
            gau_run = Gaussian()
            gau_run.call(inp_pipe=self.base_cls.base_name+'.com', 
                         out_pipe=self.base_cls.base_name+'.log',
                         dry_run=dry_run)
        os.chdir('..')

        return
    
    def _clean(self):
        """ Clean the files generated during the stage. """
        raise NotImplementedError("clean method not implemented")file_path

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

        # Check if the path exists, and make if needed
        orig_dir = Path.cwd()
        calc_dir = Path('gaussianCalcs')
        if not calc_dir.exists():
            calc_dir.mkdir()

        run_apply = print

        store_coords = []
        for a in self.alpha:
            for b in self.beta:
                for g in self.gamma:
                    test_rotation = self.base_cls.coord_object.rotate(alpha=a, beta=b, gamma=g)
                    store_coords.append(test_rotation)
                    # Write a guassian input file
                    newgau = GaussianWriter(f'gaussianCalcs/{self.base_cls.base_name}_rot_{a}_{b}_{g}.com')
                    newgau.add_block(GaussianInput(command=f"#P {self.base_cls.theory['low']} SCF(Conver=6) NoSymm Test Pop=mk IOp(6/33=2) GFInput GFPrint",
                                        initial_coordinates = test_rotation,
                                        elements = self.base_cls.coord_object.get_elements(),
                                        header=self.base_cls.header))
                    newgau.write(dry_run=dry_run)

        # Write the coordinates to a "trajectory" file
        self.write_rotation(store_coords)

        # Run the Gaussian calculations in the gaussianCalcs directory
        os.chdir(calc_dir)
        try:
            rot_count = 0
            for a in self.alpha:
                for b in self.beta:
                    for g in self.gamma:
                        self._print_rotation(a, b, g)
                        if dry_run:
                            print("Dry run: Gaussian calculations not run")
                            print("Would run the following file:")
                            print(f'-->{self.base_cls.base_name}_rot_{a}_{b}_{g}.com')
                        else:
                            gau_run = Gaussian()
                            gau_run.call(inp_pipe=f'{self.base_cls.base_name}_rot_{a}_{b}_{g}.com', 
                                    out_pipe=f'{self.base_cls.base_name}_rot_{a}_{b}_{g}.log',
                                    dry_run=dry_run)
                        rot_count += 1
                        self._print_status(rot_count, self.alpha, self.beta, self.gamma)
        finally:
            os.chdir(orig_dir)

        return
    
    def _print_rotation(self, alpha, beta, gamma):
        """ Print the rotation to the user. """
        print(f"---> Rotation: alpha={alpha}, beta={beta}, gamma={gamma}")
        return

    def _print_status(self, count, alphas, betas, gammas):
        """ Print the status of the stage. 
        
        Parameters
        ----------
        count : int
            The current count of the rotations
        alphas : list
            The list of alpha angles
        betas : list
            The list of beta angles
        gammas : list
            The list of gamma angles
        """
        total_count = len(alphas) * len(betas) * len(gammas)
        percent = count / total_count * 100
        print(f"Current Rotation Progress: {percent}")
        return
    
    def write_rotation(self, coords):
        """ Write the rotation to a file. """
        print(f"--> Writing rotations to file: gaussianCalcs/{self.base_cls.base_name}_rotations.xyz")
        with open(f'gaussianCalcs/{self.base_cls.base_name}_rotations.xyz', 'w') as file_obj:
            for frame in coords:
                SimpleXYZ(file_obj, frame)
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
        import warnings
        warnings.filterwarnings("ignore")
        
        if self.dry_run is not None:
            dry_run = self

        logfile = Path(f'gaussianCalcs/{self.base_cls.base_name}.log')
        if not logfile.exists():
            print(f"Problem with {logfile}")
            raise FileNotFoundError(f"Error (Stage {self.name}): Gaussian log file not found")
        # Convert from gaussian to mol2
        ante = Antechamber()
        ante.call(i=logfile, fi='gout',
                  o=self.base_cls.base_name+'.tmp1.mol2', fo='mol2',
                  pf='y', at=self.base_cls.atom_type,
                  dry_run = dry_run)

        # Assign the charges 
        if not dry_run:
            u1 = mda.Universe(self.base_cls.base_name+'.antechamber.mol2')
            u2 = mda.Universe(self.base_cls.base_name+'.tmp1.mol2')
            assert len(u1.atoms) == len(u2.atoms), "Number of atoms in the two files do not match"

            u2.atoms.charges = u1.atoms.charges
            """
            ag = u2.select_atoms("all")
            ag.write(self.base_cls.base_name+'.tmp2.mol2')
            # This exists because for some reason antechamber misinterprets
            # the mol2 file's blank lines in the atoms section.
            self.remove_blank_lines(self.base_cls.base_name+'.tmp2.mol2')
            """
            Mol2Writer(u2, self.base_cls.base_name+'.tmp2.mol2', selection="all").write()


        # Use antechamber to clean up the mol2 format
        ante = Antechamber()
        ante.call(i=self.base_cls.base_name+'.tmp2.mol2', fi='mol2',
                  o=self.base_cls.base_name+'.log.mol2', fo='mol2',
                  pf='y', at=self.base_cls.atom_type,
                  dry_run = dry_run)
        
        return
    
    def _clean(self):
        return

    def remove_blank_lines(self, file_path):
        """ Remove blank lines from a file.
        
        Parameters
        ----------
        file_path : str
            The path to the file to remove blank lines from
        
        Returns
        -------
        None
        
        """
        if Path(file_path).exists():
            # Read the file and filter out blank lines
            with open(file_path, 'r') as file:
                lines = file.readlines()
                non_blank_lines = [line for line in lines if line.strip()]

            # Write the non-blank lines back to the file
            with open(file_path, 'w') as file:
                file.writelines(non_blank_lines)
