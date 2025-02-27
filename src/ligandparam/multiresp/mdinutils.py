""" This code is used to generate input files for AMBER simulations. 

    The main classes are:
    - AmberCmd: used to generate the command line for sander.MPI   
    - Mdin: used to generate the mdin file
    - Disang1D: used to generate the disang file for distance restraints
    - Disang2D: used to generate the disang file for distance restraints
    - Computer: used to generate the bash script for running the simulation
    - BASH: used to generate the bash script for running the simulation
    - PERCEVAL: used to generate the bash script for running the simulation

This code was initially developed by Prof. Timothy Giese at Rutgers University, and was reproduced more or less in
its entirety here (with modifications to the code for readibility and documentation). 

"""
from collections import defaultdict as ddict
import os

class AmberCmd(object):
    """ This class is used to generate the command line for sander.MPI """
    exe    = "${AMBERHOME}/bin/sander.MPI"
    parm7  = "prmtop"
    mdin   = "mdin"
    mdout  = "mdout"
    crd7   = "inpcrd"
    rst7   = "restrt"
    refc   = "refc"
    nc     = "mdcrd"
    mdinfo = "mdinfo"
    mdfrc  = None
    dipout = None

    def __init__(self,base = None, inpcrd = None):
        """ Initialize the AmberCmd object.
        
        Parameters
        ----------
        base : str, optional
            The base name of the simulation
        inpcrd : str, optional
            The name of the inpcrd file
        
        """
        self.SetDefaults()
        if base != None: self.SetBaseName( base )
        if inpcrd != None: self.CRD7 = inpcrd

    def SetDefaults(self):
        """ Set the default values for the AmberCmd object. """
        self.reference = False
        self.EXE    = AmberCmd.exe
        self.PARM7  = AmberCmd.parm7
        self.MDIN   = AmberCmd.mdin
        self.MDOUT  = AmberCmd.mdout
        self.CRD7   = AmberCmd.crd7
        self.RST7   = AmberCmd.rst7
        self.REFC   = AmberCmd.refc
        self.NC     = AmberCmd.nc
        self.MDINFO = AmberCmd.mdinfo
        self.INPNC  = None
        self.gsize  = None
        self.gmmsize = None
        self.MDFRC  = None
        self.DIPOUT  = None

    def SetBaseName(self,base):
        """ Set the base name of the simulation. """
        self.BASE   = base
        self.MDIN   = base + ".mdin"
        self.MDOUT  = base + ".mdout"
        self.RST7   = base + ".rst7"
        self.NC     = base + ".nc"
        self.MDINFO = base + ".mdinfo"
        self.DISANG = base + ".disang"
        self.DUMPAVE= base + ".dumpave"

    def SetGroupSize(self,gsize,gmmsize=None):
        """ Set the group size. 
        
        Parameters
        ----------
        gsize : int
            The group size
        gmmsize : int, optional
            The gmm size
        """
        self.gsize = gsize
        if gmmsize is not None:
            self.gmmsize = gmmsize
            
    def CmdString(self):
        """ Generate the command line string for sander.MPI """
        crd = self.CRD7
        if crd == self.RST7:
            crd = self.RST7.replace(".rst7",".crd7")
        if self.gsize is not None:
            cmd = "-gsize %4i"%(self.gsize)
            if self.gmmsize is not None:
                cmd += " -gmmsize %4i"%(self.gmmsize)
            else:
                cmd += " -gmmsize %4i"%(self.gsize)
        else:
            cmd = ""
        cmd += " -O -p " + self.PARM7 \
               + " -i " + self.MDIN + " -o " + self.MDOUT \
               + " -c " + crd + " -r " + self.RST7 \
               + " -x " + self.NC \
               + " -inf " + self.MDINFO
        if self.REFC != "refc":
            cmd += " -ref " + self.REFC
        if self.INPNC is not None:
            cmd += " -y " + self.INPNC
        if self.MDFRC is not None:
            cmd += " -frc " + self.MDFRC
        if self.DIPOUT is not None:
            cmd += " -dipout " + self.DIPOUT
             
        return cmd

class Disang1D(object):
    """ This class is used to generate the disang file for distance restraints. """
    def __init__(self,template,r1,k1):
        """ Initialize the Disang1D object.
        
        Parameters
        ----------
        template : str
            The template file
        r1 : float
            The distance restraint
        k1 : float
            The force constant
        """
        self.TEMPLATE = template
        self.R1 = r1
        self.K1 = k1

    def WriteDisang(self,fname):
        """ Write the disang file. 

        Parameters
        ----------
        fname : str
            The name of the file
        """
        import re
        cout = open(fname,"w")
        cin = open(self.TEMPLATE,"r")
        for line in cin:
            line = re.sub(r'R1',"%.2f"%(self.R1),line)
            line = re.sub(r'K1',"%.2f"%(self.K1),line)
            cout.write(line)
        cin.close()
        cout.close()
    
class Disang2D(object):
    """ This class is used to generate the disang file for distance restraints. """
    def __init__(self,template,r1,k1,r2,k2):
        """ Initialize the Disang2D object.
        
        Parameters
        ----------
        template : str
            The template file
        r1 : float
            The distance restraint for the first atom
        k1 : float
            The force constant for the first atom
        r2 : float
            The distance restraint for the second atom
        k2 : float
            The force constant for the second atom
        """
        self.TEMPLATE = template
        self.R1 = r1
        self.K1 = k1
        self.R2 = r2
        self.K2 = k2

    def WriteDisang(self,fname):
        """ Write the disang file.
        
        Parameters
        ----------
        fname : str
            The name of the file
        """
        import re
        cout = open(fname,"w")
        cin = open(self.TEMPLATE,"r")
        for line in cin:
            line = re.sub(r'R1',"%.2f"%(self.R1),line)
            line = re.sub(r'K1',"%.2f"%(self.K1),line)
            line = re.sub(r'R2',"%.2f"%(self.R2),line)
            line = re.sub(r'K2',"%.2f"%(self.K2),line)
            cout.write(line)
        cin.close()
        cout.close()
    

class Mdin(AmberCmd):
    """ This class is used to generate the mdin file. """
    def __init__(self,base = None, inpcrd = None):
        """ Initialize the Mdin object.
        
        Parameters
        ----------
        base : str, optional
            The base name of the simulation
        inpcrd : str, optional
            The name of the inpcrd file
        """
        AmberCmd.__init__(self,base,inpcrd)
        self.title="title"
        self.cntrl = ddict( str )
        self.qmmm = ddict( str )
        self.ewald = ddict( str )
        self.dlfind = ddict( str )
        self.shake = False # used for dlfind
        self.restraints = None
        self.DUMPFREQ = 25
        self.cntrl["ntxo"]=1
        self.cntrl["cut"]=12.0
        self.cntrl["nstlim"]=1000000
        self.cntrl["dt"]=0.001
        self.cntrl["iwrap"]=1
        self.cntrl["ntf"]=2
        self.cntrl["ntc"]=2
        self.cntrl["ioutfm"]=1
        self.cntrl["ig"]=-1
        self.Set_Restart(True)
        self.Set_PrintFreq(1000)
        self.temp_changes=[]


    def Set_DumpFreq(self,dumpfreq):
        """ Set the dump frequency. 
        
        Parameters
        ----------
        dumpfreq : int
            The dump frequency
            
        """
        self.DUMPFREQ = dumpfreq
        
    def Set_PrintFreq(self,freq):
        """ Set the print frequency.

        Parameters
        ----------
        freq : int
            The print frequency
        """
        self.cntrl["ntpr"]=freq
        self.cntrl["ntwr"]=freq
        self.cntrl["ntwx"]=freq

    def Set_Restart(self,restart=True):
        """ Set the restart.

        Parameters
        ----------
        restart : bool, optional
            If True, the simulation will be restarted. Default is True
        """

        if not restart:
            #self.cntrl["imin"]=0
            self.cntrl["ntx"]=1
            self.cntrl["irest"]=0
        else:
            #self.cntrl["imin"]=0
            self.cntrl["ntx"]=5
            self.cntrl["irest"]=1

    def Add_TempChange(self,time0,time1,temp0,temp1):
        """ Add a temperature change.
        
        Parameters
        ----------
        time0 : int
            The initial time
        time1 : int
            The final time
        temp0 : float
            The initial temperature
        temp1 : float
            The final temperature
        """
        self.temp_changes.append( (time0,time1,temp0,temp1) )
        
            
    def Set_NVT(self,temp0=298.0):
        """ Set the NVT ensemble.
        
        Parameters
        ----------
        temp0 : float, optional
            The initial temperature. Default is 298.0
        
        """
        self.cntrl["ntt"]=3
        self.cntrl["ntb"]=1
        self.cntrl["ntp"]=0
        self.cntrl["temp0"]=temp0
        self.cntrl["gamma_ln"]=5.0
        self.cntrl.pop("pres0",None)
        self.cntrl.pop("taup",None)
        self.cntrl.pop("nscm",None)

    def Set_NPT(self,temp0=298.0):
        """ Set the NPT ensemble.
        
        Parameters
        ----------
        temp0 : float, optional
            The initial temperature. Default is 298.0
        """
        self.cntrl["ntt"] = 3
        self.cntrl["ntb"] = 2
        self.cntrl["ntp"] = 1
        self.cntrl["temp0"] = temp0
        self.cntrl["gamma_ln"] = 5.0
        self.cntrl["barostat"] = 2
        self.cntrl["pres0"] = 1.013
        self.cntrl["taup"] = 2.0
        self.cntrl.pop("nscm",None)

    def Set_NVE(self):
        """ Set the NVE ensemble. """
        self.cntrl["ntt"]=0
        self.cntrl["ntb"]=1
        self.cntrl["ntp"]=0
        self.cntrl.pop("temp0",None)
        self.cntrl.pop("gamma_ln",None)
        self.cntrl.pop("pres0",None)
        self.cntrl.pop("taup",None)
        self.cntrl["nscm"]=0

        
    def Set_QMMM_AM1(self,qmmask=None,qmcharge=None,tight=False):
        """ Set the QMMM AM1 method.
        
        Parameters
        ----------
        qmmask : str, optional
            The quantum mask. Default is None
        qmcharge : int, optional
            The quantum charge. Default is None
        tight : bool, optional
            If True, the simulation will be tight. Default is False
        
        """
        self.qmmm["qm_theory"] = "'AM1D'"
        self.qmmm["qmmm_switch"]=1
        self.qmmm["qm_ewald"]=1
        self.qmmm["qmshake"]=0
        if tight:
            self.qmmm["diag_routine"] = 6
            self.qmmm["tight_p_conv"] = 1
            self.qmmm["scfconv"]=1.e-11

        if qmmask is not None:
            self.qmmm["qmmask"] = qmmask
            self.cntrl["ifqnt"] = 1
        if qmcharge is not None:
            self.qmmm["qmcharge"] = qmcharge
            
        self.qmmm.pop("hfdf_theory",None)
        self.qmmm.pop("hfdf_basis",None)
        self.qmmm.pop("hfdf_mempercore",None)
        self.qmmm.pop("hfdf_ewald",None)
        self.qmmm.pop("hfdf_mm_percent",None)
        self.qmmm.pop("hfdf_qmmm_wswitch",None)


    def Set_QMMM_PBE0(self,qmmask=None,qmcharge=None,basis=None,tight=False):
        """ Set the QMMM PBE0 method.
        
        Parameters
        ----------
        qmmask : str, optional
            The quantum mask. Default is None
        qmcharge : int, optional
            The quantum charge. Default is None
        basis : str, optional
            The basis set. Default is None
        tight : bool, optional
            If True, the simulation will be tight. Default is False
        """
        if qmmask is not None:
            self.qmmm["qmmask"] = qmmask
            self.cntrl["ifqnt"] = 1
        if qmcharge is not None:
            self.qmmm["qmcharge"] = qmcharge
            
        self.qmmm["qmshake"]=0
        self.qmmm["qm_theory"] = "'HFDF'"
        self.qmmm["hfdf_theory"] = "'PBE0'"
        if basis is not None:
            self.qmmm["hfdf_basis"] = basis
        else:
            self.qmmm["hfdf_basis"] = "'6-31G*'"
        if "scfconv" not in self.qmmm:
            self.qmmm["scfconv"] = 1.e-7
        if tight:
            self.qmmm["scfconv"] = 1.e-8
            
        self.qmmm["hfdf_mempercore"] = 2000
        self.qmmm["hfdf_ewald"] = "T"
        self.qmmm["qm_ewald"] = 1
        self.qmmm["hfdf_mm_percent"] = 0.0
        self.qmmm["hfdf_qmmm_wswitch"] = 0.0
        self.qmmm["qmmm_switch"] = 0

        self.qmmm.pop("tight_p_conv",None)
        self.qmmm.pop("diag_routine",None)

    def Set_DLFIND_Minimize(self):
        """ Set the DLFind minimization. """
        self.dlfind["crdrep"] = "'HDLC'"
        self.dlfind["optalg"] = "'LBFGS'"
        self.dlfind["trustrad"] = "'GRAD'"
        self.dlfind["hessini"] = "'IDENTITY'"
        self.dlfind["hessupd"] = "'BFGS'"
        self.dlfind["neb"] = "''"
        self.dlfind["dimer"] = "''"
        self.dlfind["crdfile"] = "''"
        if "wtmask" not in self.dlfind:
            self.dlfind["wtmask"] = "''"
        self.dlfind["hessout"] = "''"
        self.dlfind["ihess"] = -1
        self.dlfind["nhess"] = -1
        self.dlfind["tol"] = 5.e-6
        self.dlfind["tole"] = 1.e-4
        self.dlfind["maxcycle"] = 1200
        self.dlfind["maxene"] = 2400
        self.dlfind["maxstep"] = -1.
        self.dlfind["minstep"] = 1.e-8
        self.dlfind["scalestep"] = 0.5
        self.dlfind["lbfgsmem"] = -1
        self.dlfind["maxupdate"] = -1
        self.dlfind["hessdel"] = 0.005
        self.dlfind["dimer_delta"] = 0.01
        self.dlfind["dimer_maxrot"] = 25
        self.dlfind["dimer_tolrot"] = 2.
        self.dlfind["dimer_mmpercent"] = 0.0

        self.dlfind.pop("crdfile",None)

        
    def Set_DLFIND_TransitionState(self,displaced_rst7):
        """ Set the DLFind transition state.
        
        Parameters
        ----------
        displaced_rst7 : str
            The name of the displaced restart file
        
        """
        self.Set_DLFIND_Minimize()
        self.dlfind["trustrad"] = "'CONST'"
        self.dlfind["hessupd"] = "'BOFILL'"
        self.dlfind["dimer"] = "'LN_NO_EXTRA'"
        self.dlfind["crdfile"] = "'%s'"%(displaced_rst7) 
        
            
    def Set_DLFIND(self,active = None):
        """ Set the DLFind.

        Parameters
        ----------
        active : str, optional
            The active atoms. Default is None
        """
        #self.cntrl = ddict( str )
        #self.qmmm = ddict( str )
        #self.ewald = ddict( str )
        self.dlfind = ddict( str )
        self.cntrl["imin"]   = 1  # minimize
        self.cntrl["ntmin"]  = 5  # read &dlfind
        self.cntrl["irest"]  = 0  # no need to read forces from rst
        self.cntrl["ntx"]    = 1
        #self.cntrl["ntwx"]   = 50 # traj freq
        #self.cntrl["ioutfm"] = 1
        #self.cntrl["ntb"]    = 1  # periodic
        #self.cntrl["iwrap"]  = 1  # wrap
        #self.cntrl["cut"]    = 12.0 # nonbond
        #self.cntrl["ig"]     = -1   # random number seed
        #self.cntrl["ifqnt"]  = 1    # read &qmmm
        #self.cntrl["nmropt"] = 1    # read DISANG,DUMPAVE
        #self.cntrl["ntf"]    = 2    # shake
        #self.cntrl["ntc"]    = 2    # shake
#        if "qmmask" not in self.qmmm:
#            self.qmmm["qmmask"] = "''"
#        if "qmcharge" not in self.qmmm:
#            self.qmmm["qmcharge"]  = 0
#        self.qmmm["writepdb"]  = 0
#        self.qmmm["qmshake"]   = 0
#        self.qmmm["verbosity"] = 0
#        self.Set_QMMM_AM1()
        self.dlfind["active"] = "''"
        self.dlfind["prefix"] = ""
        self.Set_DLFIND_Minimize()
        if active is not None:
            self.dlfind["active"] = active

    def WriteMdin(self,directory=None):
        """ Write the mdin file.

        Parameters
        ----------
        directory : str, optional
            The directory where the mdin file will be written. Default is None
        """ 
        import os.path
        fname = self.MDIN
        if directory is not None:
            fname = os.path.join(directory,self.MDIN)
            
        fh = open(fname,"w")
        fh.write("%s\n"%(self.title))
        fh.write("&cntrl\n")
        iokeys = ddict(str)
        iokeys["irest"]  = "! 0 = start, 1 = restart"
        iokeys["ntx"]    = "! 1 = start, 5 = restart"
        iokeys["ntxo"]   = "! read/write rst as formatted file"
        iokeys["iwrap"]  = "! wrap crds to unit cell"
        iokeys["ioutfm"] = "! write mdcrd as netcdf"
        iokeys["ntpr"]   = "! mdout print freq"
        iokeys["ntwx"]   = "! mdcrd print freq"
        iokeys["ntwr"]   = "! rst print freq"

        runkeys = ddict(str)
        runkeys["imin"]   = "! 0=dynamics, 1=internal minimizer"
        runkeys["ntmin"]  = "! minimization algo (5=dlfind)" 
        runkeys["nstlim"] = "! number of time steps"
        runkeys["dt"]     = "! ps/step"
        runkeys["numexch"]= "! number of replica exchange attempts"
        runkeys["ntb"]    = "! 1=NVT periodic, 2=NPT periodic, 0=no box"

        tempkeys = ddict(str)
        tempkeys["ntt"]      = "! thermostat (3=Langevin)"
        tempkeys["tempi"]    = "! initial temp for generating velocities"
        tempkeys["temp0"]    = "! target temp"
        tempkeys["gamma_ln"] = "! Langevin collision freq"

        preskeys = ddict(str)
        preskeys["ntp"]      = "! 0=no scaling, 1=isotropic, 2=anisotropic"
        preskeys["barostat"] = "! barostat (1=Berendsen, 2=MC)"
        preskeys["pres0"]    = "! pressure (bar), 1.013 bar/atm"
        preskeys["taup"]     = "! pressure relaxation time"

        shakekeys = ddict(str)
        shakekeys["ntf"] = "! 1=cpt all bond E, 2=ignore HX bond E, 3=ignore all bond E"
        shakekeys["ntc"] = "! 1=no shake, 2=HX constrained, 3=all constrained"
        shakekeys["noshakemask"] = "! do not shake these"

        has_restraints=False
        if "nmropt" in self.cntrl:
            if self.cntrl["nmropt"] == 1:
                has_restraints=True
        if len(self.temp_changes) > 0:
            self.cntrl["nmropt"] = 1
            
            
        namedkeys = []
        for ktypes in [ iokeys, runkeys, tempkeys, preskeys, shakekeys ]:
            for key in ktypes:
                namedkeys.append( key )
        
        fh.write("! IO =======================================\n")
        for key in ["irest","ntx","ntxo","iwrap","ioutfm","ntpr","ntwx","ntwr"]:
            if key in self.cntrl:
                fh.write("%11s = %-7s %s\n"%(key,str(self.cntrl[key]),iokeys[key]))
        fh.write("! DYNAMICS =================================\n")
        for key in runkeys:
            if key in self.cntrl:
                fh.write("%11s = %-7s %s\n"%(key,str(self.cntrl[key]),runkeys[key]))
        fh.write("! TEMPERATURE ==============================\n")
        for key in tempkeys:
            if key in self.cntrl:
                fh.write("%11s = %-7s %s\n"%(key,str(self.cntrl[key]),tempkeys[key]))
        fh.write("! PRESSURE  ================================\n")
        for key in preskeys:
            if key in self.cntrl:
                fh.write("%11s = %-7s %s\n"%(key,str(self.cntrl[key]),preskeys[key]))
        fh.write("! SHAKE ====================================\n")
        for key in shakekeys:
            if key in self.cntrl:
                fh.write("%11s = %-7s %s\n"%(key,str(self.cntrl[key]),shakekeys[key]))
        is_default_not_shaked = False
        if "ntf" in self.cntrl:
            if self.cntrl["ntf"] == 1:
                is_default_not_shaked = True
        if not self.shake:
            if ("noshakemask" not in self.cntrl) and (not is_default_not_shaked):
                if "active" in self.dlfind:
                    fh.write("%11s = %s\n"%("noshakemask",self.dlfind["active"]))
        fh.write("! MISC =====================================\n")
        for key in sorted(self.cntrl):
            if key not in namedkeys:
                fh.write("%11s = %s\n"%(key,str(self.cntrl[key])))

        fh.write("/\n\n")

        if "nmropt" in sorted(self.cntrl):
            if self.cntrl["nmropt"] == 1:
                if has_restraints:
                    fh.write("&wt type='DUMPFREQ', istep1 = %i &end\n"%(self.DUMPFREQ))
                for c in self.temp_changes:
                    fh.write("&wt type='TEMP0', istep1=%i, istep2=%i, value1=%.2f, value2=%.2f &end\n"%(c[0],c[1],c[2],c[3]))
                fh.write("&wt type='END' &end\n")
                if has_restraints:
                    fh.write("LISTOUT=POUT\n")
                    fh.write("DISANG=%s\n"%(self.DISANG))
                    fh.write("DUMPAVE=%s\n\n"%(self.DUMPAVE))
                    fname = self.DISANG
                    if directory is not None:
                        fname = os.path.join(directory,self.DISANG)
                    if self.restraints is not None:
                        self.restraints.WriteDisang( fname )
                    elif not os.path.isfile( fname ):
                        raise Exception("%s not found!\n"%(fname))
                

        ewald_keys = [ key for key in self.ewald ]
        if len(ewald_keys) > 0:
            fh.write("&ewald\n")
            for key in sorted(self.ewald):
                fh.write("%14s = %s\n"%(key,str(self.ewald[key])))
            fh.write("/\n\n")
                
        if "ifqnt" in self.cntrl:
            if self.cntrl["ifqnt"] > 0:
                fh.write("&qmmm\n")
                for key in sorted(self.qmmm):
                    fh.write("%14s = %s\n"%(key,str(self.qmmm[key])))
                fh.write("/\n\n")
        if "ntmin" in self.cntrl:
            if self.cntrl["ntmin"] == 5:
                self.dlfind["prefix"] = "'%s'"%(self.BASE)
                if "hessout" in self.dlfind:
                    self.dlfind["hessout"] = "'%s.hes'"%(self.BASE)
                else:
                    self.dlfind["hessout"] = "''"
                fh.write("&dlfind\n")
                for key in sorted(self.dlfind):
                    fh.write("%14s = %s\n"%(key,str(self.dlfind[key])))
                fh.write("/\n\n")
                
        fh.close()

