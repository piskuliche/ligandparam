import os
from typing import Union
import logging

import MDAnalysis as mda

from pathlib import Path
import shutil as sh

from ligandparam.stages.abstractstage import AbstractStage
from ligandparam.io.coordinates import Coordinates, SimpleXYZ, Mol2Writer
from ligandparam.io.gaussianIO import GaussianWriter, GaussianInput, GaussianReader
from ligandparam.interfaces import Gaussian, Antechamber
from ligandparam.log import get_logger

#
logger = logging.getLogger("ligandparam.gaussian")


class StageGaussian(AbstractStage):
    """
    This is class to run a basic Gaussian calculations on the ligand.
    This does three gaussian steps, one at a low level of theory, one at a higher level of theory,
    and one for the resp calculation.

    Parameters
    ----------
    inputoptions : dict
        The input options for the stage

    Returns
    -------
    None
    """

    def __init__(self, stage_name: str, in_filename: Union[Path, str], cwd: Union[Path, str], *args, **kwargs) -> None:
        super().__init__(stage_name, in_filename, cwd, *args, **kwargs)
        self.in_pdb = Path(in_filename)
        self.out_gaussian_log = Path(kwargs["out_gaussian_log"])

        self._validate_input_paths(**kwargs)
        self.theory = kwargs.get("theory", {"low": "HF/6-31G*", "high": "PBE1PBE/6-31G*"})
        self.net_charge = kwargs.get("net_charge", 0.0)
        self.force_gaussian_rerun = kwargs.get("force_gaussian_rerun", False)
        self.gaussian_cwd = Path(self.cwd, "gaussianCalcs")
        self.out_log = self.gaussian_cwd / f"{self.out_gaussian_log.stem}.log"
        self._add_outputs(self.out_log)

        self.header = [f"%NPROC={self.nproc}', f'%MEM={self.mem}MB"]

        # No required files for this stage to execute.
        return

    def _validate_input_paths(self, **kwargs):
        for opt in ("gaussian_root", "gauss_exedir", "gaussian_binary", "gaussian_scratch"):
            try:
                setattr(self, opt, kwargs.get(opt, ""))
            except KeyError:
                raise ValueError(f"ERROR: Please provide {opt} option as a keyword argument.")
        if self.gaussian_binary is None:
            self.gaussian_binary = "g16"

        # if not self.gaussian_root.is_dir():
        #     raise ValueError(f"ERROR: input gaussian root is not an existing directory: {self.gaussian_root}")
        # for dir in str(self.gauss_exedir).split(":"):
        #     if not Path(dir).is_dir():
        #         raise ValueError(f"ERROR: input gaussian executable directory is not an existing directory: {dir}")
        # if not self.gaussian_scratch.is_dir():
        #     raise ValueError(f"ERROR: input gaussian scratch is not an existing directory: {self.gaussian_scratch}")
        # if not self.gaussian_binary.is_file():
        #     raise ValueError(f"ERROR: input gaussian binary is not an existing file: {self.gaussian_binary}")

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
        stageheader = self.header

        stageheader.append(f"%chk={self.in_pdb.stem}.antechamber.chk")

        # Set up the Gaussian Block - it does not yet write anything,
        # so this part can be set up before the Gaussian calculations are run.
        gau = GaussianWriter(Path(self.gaussian_cwd, f"{self.out_gaussian_log.stem}.com"))
        gau.add_block(GaussianInput(command=f"#P {self.theory['low']} OPT(CalcFC)",
                                    initial_coordinates=self.coord_object.get_coordinates(),
                                    elements=self.coord_object.get_elements(),
                                    charge=self.net_charge,
                                    header=stageheader))
        gau.add_block(GaussianInput(command=f"#P {self.theory['high']} OPT(CalcFC) GEOM(ALLCheck) Guess(Read)",
                                    charge=self.net_charge,
                                    header=stageheader))
        gau.add_block(GaussianInput(
            command=f"#P {self.theory['low']} GEOM(AllCheck) Guess(Read) NoSymm Pop=mk IOp(6/33=2) GFInput GFPrint",
            charge=self.net_charge,
            header=stageheader))

        self.gaussian_cwd.mkdir(exist_ok=True)

        gau_complete = False
        # Check if the Gaussian calculation has already been run
        if os.path.exists(self.out_gaussian_log):
            reader = GaussianReader(self.out_gaussian_log)
            if reader.check_complete():
                self.logger.info("Gaussian calculation already complete")
                gau_complete = True

        # Check if the Gaussian calculation should be rerun
        if self.force_gaussian_rerun:
            gau_complete = False

        if not gau_complete:
            gau.write(dry_run=dry_run)

        # Run the Gaussian calculations in the gaussianCalcs directory
        if not gau_complete:
            gau_run = Gaussian(cwd=self.gaussian_cwd, logger=self.logger, gaussian_root=self.gaussian_root,
                               gauss_exedir=self.gauss_exedir,
                               gaussian_binary=self.gaussian_binary, gaussian_scratch=self.gaussian_scratch)
            gau_run.call(inp_pipe=self.in_pdb.stem + '.com',
                         out_pipe=self.out_log.name,
                         dry_run=dry_run)

            # Move the Gaussian log file to the output location
            sh.move(self.out_log, self.out_gaussian_log)

        return

    def _clean(self):
        """ Clean the files generated during the stage. """
        raise NotImplementedError("clean method not implemented")


class StageGaussianRotation(AbstractStage):

    def __init__(self, stage_name, in_filename, cwd, alpha=[0.0], beta=[0.0], gamma=[0.0], inputoptions=None, *args,
                 **kwargs) -> None:
        """ This is class to rotate the ligand and run Gaussian calculations of the resp charges
        for each rotated ligand. 
        
        Parameters
        ----------
        name : str
            The name of the stage
        alpha : list
            The list of alpha angles to rotate the ligand
        beta : list
            The list of beta angles to rotate the ligand
        gamma : list
            The list of gamma angles to rotate the ligand
        inputoptions : dict
            The input options for the stage
        """

        super().__init__(stage_name, in_filename, cwd, *args, **kwargs)
        self.in_mol2 = Path(in_filename)
        self.out_gaussian_log = Path(kwargs["out_gaussian_log"])
        self._parse_inputoptions(inputoptions)

        self.alpha = [float(a) for a in alpha]
        self.beta = [float(b) for b in beta]
        self.gamma = [float(g) for g in gamma]

        self.header = [f"%NPROC={self.nproc}', f'%MEM={self.mem}MB"]

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
                    test_rotation = self.coord_object.rotate(alpha=a, beta=b, gamma=g)
                    store_coords.append(test_rotation)
                    # Write a guassian input file
                    newgau = GaussianWriter(
                        f'gaussianCalcs/{self.out_gaussian_log.stem}_rot_{a:0.2f}_{b:0.2f}_{g:0.2f}.com')
                    newgau.add_block(GaussianInput(
                        command=f"#P {self.theory['low']} SCF(Conver=6) NoSymm Test Pop=mk IOp(6/33=2) GFInput GFPrint",
                        initial_coordinates=test_rotation,
                        elements=self.coord_object.get_elements(),
                        header=self.header))
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
                            self.logger.info("Dry run: Gaussian calculations not run")
                            self.logger.info("Would run the following file:")
                            self.logger.info(f'-->{self.out_gaussian_log.stem}_rot_{a:.2f}_{b:.2f}_{g:.2f}.com')
                        else:
                            gau_run = Gaussian(cwd=self.cwd, logger=self.logger)
                            gau_run.call(inp_pipe=f'{self.out_gaussian_log.stem}_rot_{a:.2f}_{b:.2f}_{g:.2f}.com',
                                         out_pipe=f'{self.out_gaussian_log.stem}_rot_{a:.2f}_{b:.2f}_{g:.2f}.log',
                                         dry_run=dry_run)
                        rot_count += 1
                        self._print_status(rot_count, self.alpha, self.beta, self.gamma)
        finally:
            os.chdir(orig_dir)

        return

    def _print_rotation(self, alpha, beta, gamma):
        """ Print the rotation to the user. """
        self.logger.info(f"---> Rotation: alpha={alpha}, beta={beta}, gamma={gamma}")
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
        self.logger.info(f"Current Rotation Progress: {percent:.2f}%%")
        return

    def write_rotation(self, coords):
        """ Write the rotation to a file. """
        self.logger.info(f"--> Writing rotations to file: gaussianCalcs/{self.out_gaussian_log.stem}_rotations.xyz")
        with open(f'gaussianCalcs/{self.out_gaussian_log.stem}_rotations.xyz', 'w') as file_obj:
            for frame in coords:
                SimpleXYZ(file_obj, frame)
        return

    def _clean(self):
        return


class StageGaussiantoMol2(AbstractStage):

    def __init__(self, stage_name: str, in_filename: Union[Path, str], cwd: Union[Path, str], *args, **kwargs) -> None:
        """ Convert Gaussian output to mol2 format. 
    
        This class converts the Gaussian output to mol2 format, and assigns the charges to the mol2 file.
        
        Parameters
        ----------
        name : str
            The name of the stage
        inputoptions : dict
            The input options for the stage
        dry_run : bool, optional
            If True, the stage will not be executed, but the function will print the commands that would
        
        Returns
        -------
        None
        
        """

        super().__init__(stage_name, in_filename, cwd, *args, **kwargs)
        self.in_log = Path(in_filename)
        self.in_mol2 = Path(kwargs["in_mol2"])
        self.out_mol2 = Path(kwargs["out_mol2"])
        self.temp1_mol2 = Path(self.cwd, f"{self.out_mol2.stem}.tmp1.mol2")
        self.temp2_mol2 = Path(self.cwd, f"{self.out_mol2.stem}.tmp2.mol2")

        self._parse_inputoptions(kwargs)
        self.dry_run = kwargs.get("dry_run", False)
        self.add_required(self.in_log)

        self.header = [f"%NPROC={getattr(kwargs, 'nproc', 12)}', f'%MEM={getattr(kwargs, 'mem', 48)}GB"]

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

        if not self.in_log.exists():
            self.logger.info(f"Problem with {self.in_log}")
            raise FileNotFoundError(f"Error (Stage {self.stage_name}): Gaussian log file not found")
        # Convert from gaussian to mol2
        ante = Antechamber(cwd=self.cwd, logger=self.logger, nproc=self.nproc)
        ante.call(i=self.in_log, fi='gout',
                  o=self.temp1_mol2, fo='mol2',
                  pf='y', at=self.atom_type,
                  gn=f"%nproc={self.nproc}", gm=f"%mem={self.mem}MB",
                  dry_run=dry_run)

        # Assign the charges 
        if not dry_run:
            u1 = mda.Universe(self.in_mol2)
            u2 = mda.Universe(self.temp1_mol2)
            assert len(u1.atoms) == len(u2.atoms), "Number of atoms in the two files do not match"

            u2.atoms.charges = u1.atoms.charges
            """
            ag = u2.select_atoms("all")
            ag.write(self.name+'.tmp2.mol2')
            # This exists because for some reason antechamber misinterprets
            # the mol2 file's blank lines in the atoms section.
            self.remove_blank_lines(self.name+'.tmp2.mol2')
            """
            Mol2Writer(u2, self.temp2_mol2, selection="all").write()

        # Use antechamber to clean up the mol2 format
        ante = Antechamber(cwd=self.cwd, logger=self.logger, nproc=self.nproc)
        ante.call(i=self.temp2_mol2, fi='mol2',
                  o=self.out_mol2, fo='mol2',
                  pf='y', at=self.atom_type,
                  gn=f"%nproc={self.nproc}", gm=f"%mem={self.mem}MB",
                  dry_run=dry_run)

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
