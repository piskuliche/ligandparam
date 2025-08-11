import warnings
from pathlib import Path
from typing import Optional,  Union, Any
from typing_extensions import override

import MDAnalysis as mda

from ligandparam.stages.abstractstage import AbstractStage
from ligandparam.interfaces import Leap
from ligandparam.io.leapIO import LeapWriter
from ligandparam.log import get_logger


class StageBuild(AbstractStage):
    """
    A stage for building molecular systems using Antechamber and Leap.

    Attributes
    ----------
    target_pdb : str
        The target PDB file.
    build_type : int
        The type of build to perform (0 for aqueous, 1 for gas, 2 for target).
    concentration : float
        The concentration of ions.
    rbuffer : float
        The buffer radius.
    leaprc : list of str
        The Leaprc configuration files.
    """

    @override
    def __init__(self, stage_name: str, name: Union[Path, str], cwd: Union[Path, str], *args, **kwargs) -> None:
        """ This class is used to initialize from pdb to mol2 file using Antechamber.
        
        Parameters
        ----------
        name : str
            The name of the stage
        build_type : str
            The type of build to perform [aq, gas, or target]
        target_pdb : str
            The target pdb file
        concentration : float
            The concentration of the ions
        rbuffer : float
            The buffer radius
        inputoptions : dict
            The input options
        """
        super().__init__(stage_name, name, cwd, *args, **kwargs)
        
        self.target_pdb = kwargs['target_pdb']
        self.build_type = self._validate_build_type(kwargs.get('build_type', 'aq'))
        self.concentration = kwargs.get('concentration', 0.14)
        self.rbuffer = kwargs.get('rbuffer', 9.0)
        self.leaprc = kwargs.get('leaprc', None)
        if not self.leaprc:
            self.leaprc = ["leaprc.water.OPC"]



        self.add_required(f"{self.name}.resp.mol2")
        self.add_required(f"{self.name}.frcmod")
        self.add_required(f"{self.name}.lib")

    def _validate_build_type(self, build_type: str) -> int:
        """
        Validate the build type.

        Parameters
        ----------
        build_type : str
            The build type to validate ('aq', 'gas', or 'target').

        Returns
        -------
        int
            The validated build type as an integer.

        Raises
        ------
        ValueError
            If the build type is not valid.
        """
        if build_type.lower() == 'aq':
            return 0
        elif build_type.lower() == 'gas':
            return 1
        elif build_type.lower() == 'target':
            self.add_required(f"{self.target_pdb}")
            return 2
        else:
            raise ValueError("ERROR: Please provide a valid build type: aq, gas, or target")

    def _append_stage(self, stage: "AbstractStage") -> "AbstractStage":
        """
        Append a stage to the current stage.

        Parameters
        ----------
        stage : AbstractStage
            The stage to append.

        Returns
        -------
        AbstractStage
            The appended stage.
        """
        return stage

    def execute(self, dry_run=False, nproc: Optional[int]=None, mem: Optional[int]=None) -> Any:
        """
        Execute the Gaussian calculations.

        Parameters
        ----------
        dry_run : bool, optional
            If True, the stage will not be executed, but the function will print the commands that would be run.
        nproc : int, optional
            The number of processors to use (default is None).
        mem : int, optional
            The amount of memory to use (in GB, default is None).

        Returns
        -------
        None
        """
        if self.build_type == 0:
            self._aq_build(dry_run=dry_run)
        elif self.build_type == 1:
            self._gas_build(dry_run=dry_run)
        elif self.build_type == 2:
            self._target_build(dry_run=dry_run)

    def _clean(self):
        """
        Clean the files generated during the stage.

        Raises
        ------
        NotImplementedError
            If the method is not implemented.
        """
        raise NotImplementedError("clean method not implemented")

    def _aq_build(self, dry_run=False):
        """
        Build the ligand in an aqueous environment.

        Parameters
        ----------
        dry_run : bool, optional
            If True, the stage will not be executed, but the function will print the commands that would be run.
        """
        aqleap = LeapWriter("aq")
        # Add the leaprc files
        for rc in self.leaprc:
            aqleap.add_leaprc(rc)

        solvent = None
        for lrc in self.leaprc:
            if "OPC" in lrc:
                solvent = "OPCBOX"
            elif "tip3p" in lrc:
                solvent = "TIP3PBOX"
            elif "tip4pew" in lrc:
                solvent = "TIP4PEWBOX"
        if solvent is None:
            solvent = "TIP3PBOX"

        # Add the leap commands
        aqleap.add_line(f"loadamberparams {self.name}.frcmod")
        aqleap.add_line(f"loadoff {self.name}.lib")
        aqleap.add_line(f"mol = loadmol2 {self.name}.resp.mol2")
        aqleap.add_line("\n")
        # Add counter ions
        aqleap.add_line(f"addions mol NA 0")
        aqleap.add_line(f"solvateOct mol {solvent} {self.buffer}")
        aqleap.add_line("\n")
        aqleap.add_line(f"saveamberparm mol {self.name}_aq_noions.parm7 {self.name}_aq_noions.rst7")
        aqleap.add_line("quit")
        # Write the leap input file
        aqleap.write()
        # Call the leap program to run initial check
        leap = Leap()
        leap.call(f="tleap.aq.in", dry_run=dry_run)
        num_NA, num_Cl = self.Get_Num_Ions(self.name + "_aq_noions.parm7")
        # Call the leap program to add ions
        aqleap.remove_line("quit")
        if self.concentration > 0.0:
            aqleap.add_line(f"addionsrand mol NA {num_NA} CL {num_Cl} 6.0")
            aqleap.add_line(f"saveamberparm mol {self.name}_aq.parm7 {self.name}_aq.rst7")
            aqleap.add_line("quit")
            aqleap.write()
            leap = Leap()
            leap.call(f="tleap.aq.in", dry_run=dry_run)

    def _gas_build(self, dry_run=False):
        """
        Build the ligand in a gas environment.

        Parameters
        ----------
        dry_run : bool, optional
            If True, the stage will not be executed, but the function will print the commands that would be run.
        """
        gasleap = LeapWriter("gas")
        # Add the leaprc files
        for rc in self.leaprc:
            gasleap.add_leaprc(rc)
        # Add the leap commands
        gasleap.add_line(f"loadamberparams {self.name}.frcmod")
        gasleap.add_line(f"loadoff {self.name}.lib")
        gasleap.add_line(f"mol = loadmol2 {self.name}.resp.mol2")
        gasleap.add_line("\n")
        gasleap.add_line(f"saveamberparm mol {self.name}_gas.parm7 {self.name}_gas.rst7")
        gasleap.add_line("quit")
        # Write the leap input file
        gasleap.write()
        # Call the leap program
        leap = Leap()
        leap.call(f="tleap.gas.in", dry_run=dry_run)

    def _target_build(self, dry_run=False):
        """
        Build the ligand in the target environment.

        Parameters
        ----------
        dry_run : bool, optional
            If True, the stage will not be executed, but the function will print the commands that would be run.
        """
        self.check_target()
        targetleap = LeapWriter("target")
        # Add the leaprc files
        if len(self.leaprc) == 0:
            targetleap.add_leaprc("leaprc.water.OPC")

        for rc in self.leaprc:
            targetleap.add_leaprc(rc)

        solvent = None
        for lrc in self.leaprc:
            if "OPC" in lrc:
                solvent = "OPCBOX"
            elif "tip3p" in lrc:
                solvent = "TIP3PBOX"
            elif "tip4pew" in lrc:
                solvent = "TIP4PEWBOX"
        if solvent is None:
            solvent = "TIP3PBOX"

        # Add the leap commands
        targetleap.add_line(f"loadamberparams {self.name}.frcmod")
        targetleap.add_line(f"loadoff {self.name}.lib")
        # targetleap.add_line(f"mol = loadmol2 {self.name}.resp.mol2")
        targetleap.add_line(f"mol = loadpdb {self.target_pdb}")
        targetleap.add_line("\n")
        targetleap.add_line(f"savepdb mol {self.name}_in_target.pdb")
        # Add counter ions
        targetleap.add_line(f"addions mol NA 0")
        targetleap.add_line(f"solvateoct mol {solvent} {self.buffer}")
        targetleap.add_line("\n")
        targetleap.add_line(f"saveamberparm mol {self.name}_target_noions.parm7 {self.name}_target_noions.rst7")
        targetleap.add_line("quit")
        # Write the leap input file
        targetleap.write()
        # Call the leap program to run initial check
        leap = Leap()
        leap.call(f="tleap.target.in", dry_run=dry_run)
        num_NA, num_Cl = self.Get_Num_Ions(self.name + "_target_noions.parm7")
        # Call the leap program to add ions
        targetleap.remove_line("quit")
        if self.concentration > 0.0:
            targetleap.add_line(f"addionsrand mol NA {num_NA} CL {num_Cl} 6.0")
            targetleap.add_line(f"saveamberparm mol {self.name}_target.parm7 {self.name}_target.rst7")
            targetleap.add_line("quit")
            targetleap.write()
            leap = Leap()
            leap.call(f="tleap.target.in", dry_run=dry_run)

    def Get_Num_Ions(self, parm7, wat_resname="WAT"):
        """
        Get the number of ions needed for the system.

        Parameters
        ----------
        parm7 : str
            The parameter file for the system.
        wat_resname : str, optional
            The residue name for water molecules (default is 'WAT').

        Returns
        -------
        tuple of int
            The number of sodium (NA) and chloride (CL) ions needed.

        Raises
        ------
        ValueError
            If the concentration of ions is negative.
        """
        water_concentration = 55.
        u = mda.Universe(parm7)
        total_charge = sum(u.atoms.charges)
        num_waters = len(u.select_atoms("resname WAT").residues)
        num_NA = len(u.select_atoms("resname NA")) + len(u.select_atoms("resname NA+"))
        num_CL = len(u.select_atoms("resname CL")) + len(u.select_atoms("resname CL-"))
        non_ion_charge = total_charge - num_NA + num_CL

        conc_na = ((num_waters + num_NA + num_CL) * self.concentration / water_concentration) - num_NA - (
            non_ion_charge if non_ion_charge < 0 else 0)
        conc_cl = ((num_waters + num_NA + num_CL) * self.concentration / water_concentration) - num_CL - (
            non_ion_charge if non_ion_charge > 0 else 0)

        parmconc = 0
        if num_waters > 0:
            parmconc = min(num_NA, num_CL) * water_concentration / (num_waters + num_NA + num_CL)
        self.logger.info(f"-> Current system is {total_charge}")
        self.logger.info(f"-> Current system has {non_ion_charge} non-ion charge")
        self.logger.info(f"-> Current system has {num_waters} water molecules")
        self.logger.info(f"-> Current system has {num_waters} water molecules")
        self.logger.info(f"-> Current system has {num_NA} NA ions")
        self.logger.info(f"-> Current system has {num_CL} CL ions")
        self.logger.info(f"-> Current concentration is {parmconc}")
        if conc_na > 0:
            num_NA = int(conc_na)
        else:
            raise ValueError("ERROR: Negative concentration of NA ions")
        if conc_cl > 0:
            num_CL = int(conc_cl)
        else:
            raise ValueError("ERROR: Negative concentration of CL ions")
        return num_NA, num_CL

    def check_target(self):
        """
        Check that the target PDB file is correct.

        Raises
        ------
        ValueError
            If the ligand residue name is not in the target PDB file.
        """
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            u = mda.Universe(self.target_pdb)
            u2 = mda.Universe(self.name + ".resp.mol2")
            lig_resname = u2.residues.resnames[0]
            if lig_resname not in u.residues.resnames:
                raise ValueError(f"ERROR: The ligand residue name {lig_resname} is not in the target pdb file.")
