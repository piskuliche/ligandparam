import MDAnalysis as mda

from ligandparam.stages.abstractstage import AbstractStage
from ligandparam.interfaces import Leap
from ligandparam.io.leapIO import LeapWriter



class StageBuild(AbstractStage):
    """ 

    Parameters
    ----------
    name : str
        Name of the stage.
 : object
        Object of the base class.

    """
    def __init__(self, name, build_type='aq', target_pdb=None, concentration=0.14, rbuffer=9.0, inputoptions=None) -> None:
        """ This class is used to initialize from pdb to mol2 file using Antechamber.
        
        Parameters
        ----------
        name : str
            The name of the stage
     : Ligand
            The base class of the ligand
        build_type : str
            The type of build to perform [aq, gas, or target]
        target_pdb : str
            The target pdb file
        concentration : float
            The concentration of the ions
        rbuffer : float
            The buffer radius
        """
        self.name = name
        self._parse_inputoptions(inputoptions)
        self.concentration = concentration
        self.target_pdb = target_pdb
        self.buffer = rbuffer
        
        self.add_required(f"{self.base_name}.resp.mol2")
        self.add_required(f"{self.base_name}.frcmod")
        self.add_required(f"{self.base_name}.off")

        if build_type.lower() == 'aq':
            self.build_type = 0
        elif build_type.lower() == 'gas':
            self.build_type = 1
        elif build_type.lower() == 'target':
            self.build_type = 2
            self.add_required(f"{self.target_pdb}")
        else:
            raise ValueError("ERROR: Please provide a valid build type. [aq, gas, or target]")
        return
    
    def _append_stage(self, stage: "AbstractStage") -> "AbstractStage":
        """ Appends the stage. """
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
        if self.build_type == 0:
            self._aq_build(dry_run=dry_run)
        elif self.build_type == 1:
            self._gas_build(dry_run=dry_run)
        elif self.build_type == 2:
            self._target_build(dry_run=dry_run)

        
    def _clean(self):
        """ Clean the files generated during the stage. """
        raise NotImplementedError("clean method not implemented")
    
    def _aq_build(self, dry_run=False):
        """ Build the ligand in aqueous environment. """
        aqleap = LeapWriter("aq")
        # Add the leaprc files
        if len(self.leaprc) == 0:
            aqleap.add_leaprc("leaprc.water.OPC")

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
        aqleap.add_line(f"loadamberparams {self.base_name}.frcmod")
        aqleap.add_line(f"loadoff {self.base_name}.off")
        aqleap.add_line(f"mol = loadmol2 {self.base_name}.resp.mol2")
        aqleap.add_line("\n")
        # Add counter ions
        aqleap.add_line(f"addions mol NA 0")
        aqleap.add_line(f"solvateOct mol {solvent} {self.buffer}")
        aqleap.add_line("\n")
        aqleap.add_line(f"saveamberparm mol {self.base_name}_aq_noions.parm7 {self.base_name}_aq_noions.rst7")
        aqleap.add_line("quit")
        # Write the leap input file
        aqleap.write()
        # Call the leap program to run initial check
        leap = Leap()
        leap.call(f="tleap.aq.in", dry_run = dry_run)
        num_NA, num_Cl = self.Get_Num_Ions(self.base_name+"_aq_noions.parm7")
        # Call the leap program to add ions
        aqleap.remove_line("quit")
        if self.concentration > 0.0:
            aqleap.add_line(f"addionsrand mol NA {num_NA} CL {num_Cl} 6.0")
            aqleap.add_line(f"saveamberparm mol {self.base_name}_aq.parm7 {self.base_name}_aq.rst7")
            aqleap.add_line("quit")
            aqleap.write()
            leap = Leap()
            leap.call(f="tleap.aq.in", dry_run = dry_run)

    
    def _gas_build(self, dry_run=False):
        """ Build the ligand in gas environment. """
        gasleap = LeapWriter("gas")
        # Add the leaprc files
        for rc in self.leaprc:
            gasleap.add_leaprc(rc)
        # Add the leap commands
        gasleap.add_line(f"loadamberparams {self.base_name}.frcmod")
        gasleap.add_line(f"loadoff {self.base_name}.off")
        gasleap.add_line(f"mol = loadmol2 {self.base_name}.resp.mol2")
        gasleap.add_line("\n")
        gasleap.add_line(f"saveamberparm mol {self.base_name}_gas.parm7 {self.base_name}_gas.rst7")
        gasleap.add_line("quit")
        # Write the leap input file
        gasleap.write()
        # Call the leap program
        leap = Leap()
        leap.call(f="tleap.gas.in", dry_run = dry_run)
    
    def _target_build(self, dry_run=False): 
        """ Build the ligand in the target environment. """
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
        targetleap.add_line(f"loadamberparams {self.base_name}.frcmod")
        targetleap.add_line(f"loadoff {self.base_name}.off")
        #targetleap.add_line(f"mol = loadmol2 {self.base_name}.resp.mol2")
        targetleap.add_line(f"mol = loadpdb {self.target_pdb}")
        targetleap.add_line("\n")
        targetleap.add_line(f"savepdb mol {self.base_name}_in_target.pdb")
        # Add counter ions
        targetleap.add_line(f"addions mol NA 0")
        targetleap.add_line(f"solvateoct mol {solvent} {self.buffer}")
        targetleap.add_line("\n")
        targetleap.add_line(f"saveamberparm mol {self.base_name}_target_noions.parm7 {self.base_name}_target_noions.rst7")
        targetleap.add_line("quit")
        # Write the leap input file
        targetleap.write()
        # Call the leap program to run initial check
        leap = Leap()
        leap.call(f="tleap.target.in", dry_run = dry_run)
        num_NA, num_Cl = self.Get_Num_Ions(self.base_name+"_target_noions.parm7")
        # Call the leap program to add ions
        targetleap.remove_line("quit")
        if self.concentration > 0.0:
            targetleap.add_line(f"addionsrand mol NA {num_NA} CL {num_Cl} 6.0")
            targetleap.add_line(f"saveamberparm mol {self.base_name}_target.parm7 {self.base_name}_target.rst7")
            targetleap.add_line("quit")
            targetleap.write()
            leap = Leap()
            leap.call(f="tleap.target.in", dry_run = dry_run)
    
    def Get_Num_Ions(self, parm7, wat_resname="WAT"):
        """ Get the number of ions needed for the system. """
        water_concentration = 55.
        u = mda.Universe(parm7)
        total_charge = sum(u.atoms.charges)
        num_waters = len(u.select_atoms("resname WAT").residues)
        num_NA = len(u.select_atoms("resname NA"))+len(u.select_atoms("resname NA+"))
        num_CL = len(u.select_atoms("resname CL"))+len(u.select_atoms("resname CL-"))
        non_ion_charge = total_charge - num_NA + num_CL

        conc_na = ((num_waters + num_NA + num_CL) * self.concentration / water_concentration) - num_NA - ( non_ion_charge if non_ion_charge < 0 else 0)
        conc_cl = ((num_waters + num_NA + num_CL) * self.concentration / water_concentration) - num_CL - ( non_ion_charge if non_ion_charge > 0 else 0)

        parmconc = 0
        if num_waters > 0:
            parmconc = min(num_NA, num_CL) * water_concentration / (num_waters + num_NA + num_CL)
        print(f"-> Current system is {total_charge}")
        print(f"-> Current system has {non_ion_charge} non-ion charge")
        print(f"-> Current system has {num_waters} water molecules")
        print(f"-> Current system has {num_waters} water molecules")
        print(f"-> Current system has {num_NA} NA ions")
        print(f"-> Current system has {num_CL} CL ions")
        print(f"-> Current concentration is {parmconc}")
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
        """ Check that the target pdb file is correct. """
        u = mda.Universe(self.target_pdb)
        u2 = mda.Universe(self.base_name+".resp.mol2")
        lig_resname = u2.residues.resnames[0]
        if lig_resname not in u.residues.resnames:
            raise ValueError(f"ERROR: The ligand residue name {lig_resname} is not in the target pdb file.")

