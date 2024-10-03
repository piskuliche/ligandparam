#!/usr/bin/env python3

import sys
import glob
import copy
import parmed
import parmutils
from collections import defaultdict as ddict
import os.path
from io import StringIO as StringIO
import subprocess

from . import parmutils
from . import mdinutils
import os
import os.path



def PrintPerturbedCharges(native_fit,model_fit,extra=[],exclude=[]):

    native = native_fit.states[0]
    model  = model_fit.states[0]
    
    mm_charges     = ddict( lambda: ddict( float ) )
    native_charges = ddict( lambda: ddict( float ) )
    model_charges  = ddict( lambda: ddict( float ) )
    for a in native.frag.atomsel:
        res  = native.parm.atoms[a].residue.idx+1
        name = native.parm.atoms[a].name

        q,stddev,stderr = GetAvgStdDevAndErr(native.charge_data[a])
        native_charges[res][name] = q
        mm_charges[res][name] = native.mmcharges[a]

    for a in model.frag.atomsel:
        res  = model.parm.atoms[a].residue.idx+1
        name = model.parm.atoms[a].name

        q,stddev,stderr = GetAvgStdDevAndErr(model.charge_data[a])
        model_charges[res][name] = q


    new_charges = ddict( lambda: ddict( float ) )
   

    print("# %-4s %-5s %11s %11s %11s %11s %11s"%("Res","Atom","New Q","delta Q","New Resp","Nat. Resp","Nat. Q"))
    for res in sorted( model_charges ):
        rname = model.parm.residues[ res-1 ].name

        atms = []
        for a in model.frag.atomsel:
            r  = model.parm.atoms[a].residue.idx+1
            name = model.parm.atoms[a].name
            if r == res:
                atms.append( name )
                
        for atm in sorted(mm_charges[res]):
            if atm not in model_charges[res]:
                if atm in native_charges[res] or mm_charges[res][atm] == 0:
                    atms.append( atm )

        amber_sum=0
        native_sum=0
        model_sum=0
        new_sum=0
        delta_sum=0
                    
        for atm in atms:
            amber_charge   = "%11s"%("")
            native_charge  = "%11s"%("")
            model_charge   = "%11s"%("")
            dq_charge   = "%11s"%("")
            new_charge  = "%11s"%("")
            if atm in mm_charges[res]:
                amber_charge="%11.6f"%(mm_charges[res][atm])
                amber_sum += mm_charges[res][atm]
            if atm in native_charges[res]:
                native_charge="%11.6f"%(native_charges[res][atm])
                native_sum += native_charges[res][atm]
            dq=0
            charge=0
            if res in model_charges:
                if atm in model_charges[res]:
                    model_charge="%11.6f"%(model_charges[res][atm])
                    model_sum += model_charges[res][atm]
                    dq = (model_charges[res][atm]-native_charges[res][atm])
                    charge = mm_charges[res][atm] + dq
                    dq_charge = "%11.6f"%(dq)
                    new_charge = "%11.6f"%(charge)
                    new_charges[res][atm] = charge
                else:
                    dq = -native_charges[res][atm]
                    dq_charge = "%11.6f"%(dq)
                new_sum += charge
                delta_sum += dq
                
            print("# %-4s %-5s %s %s %s %s %s"%( res,atm, new_charge, dq_charge, model_charge, native_charge, amber_charge ))
        print("# %-4s %-5s %11.6f %11.6f %11.6f %11.6f %11.6f"%( res,"SUM", new_sum, delta_sum, model_sum, native_sum, amber_sum ))
        print("")


            
    for res in sorted( model_charges ):
        rname = model.parm.residues[ res-1 ].name
        print("")
        print("function set_charge_%s_%i"%(rname,res))
        print("{")
        print("   selection=\"$1\"")
        print("   cat <<EOF")


        atms = []
        for a in model.frag.atomsel:
            r  = model.parm.atoms[a].residue.idx+1
            name = model.parm.atoms[a].name
            if r == res:
                atms.append( name )

                
        seen=[]
        for atm in atms:
            if atm in exclude:
                continue
            else:
                print("set ${selection}.%-5s charge %10.6f"%(atm,new_charges[res][atm]))
                seen.append(atm)
        for atm in sorted(mm_charges[res]):
            if atm not in model_charges[res]:
                if atm in native_charges[res] or mm_charges[res][atm] == 0:
                    if atm in exclude:
                        continue
                    else:
                        print("set ${selection}.%-5s charge   0.0"%(atm))
                        seen.append(atm)
        for atm in extra:
            if atm not in seen:
                print("set ${selection}.%-5s charge   0.0"%(atm))
                    
        print("EOF")
        print("}")
        print("")







def ReadGauEsp(fname):
    import re
    import os
    import sys

    #sys.stderr.write( "ReadGauEsp %s\n"%(fname))


    if not os.path.exists(fname):
        raise Exception("ReadGauEsp file not found: %s"%(fname))
    
    fh = open(fname,"r")
    
    crds=[]
    pts=[]
    esp=[]

    atomic_center_line = re.compile(r" +Atomic Center[ 0-9]+ is at +([\-0-9]{1,3}\.[0-9]{6}) *([\-0-9]{1,3}\.[0-9]{6}) *([\-0-9]{1,3}\.[0-9]{6})")
    fit_center_line = re.compile(r" +ESP Fit Center[ 0-9]+ is at +([\-0-9]{1,3}\.[0-9]{6}) *([\-0-9]{1,3}\.[0-9]{6}) *([\-0-9]{1,3}\.[0-9]{6})")
    esp_line = re.compile(r"[ 0-9]+ Fit +([\-\.0-9]+)")

    for line in fh:
        m = atomic_center_line.match(line)
        if m is not None:
            crds.extend( [ float( m.group(1) ), float( m.group(2) ), float( m.group(3) ) ] )
            continue
        m = fit_center_line.match(line)
        if m is not None:
            pts.extend( [ float( m.group(1) ), float( m.group(2) ), float( m.group(3) ) ] )
            continue
        m = esp_line.match(line)
        if m is not None:
            esp.append( float( m.group(1) ) )
            continue
    if len(pts) != 3*len(esp):
        raise Exception("# pts (%i) != # esp values (%i) in %s"%(len(pts)//3,len(esp),fname))
   
    #sys.stderr.write("crds=%s\n"%(str(len(crds))))
    #sys.stderr.write("pts=%s\n"%(str(len(pts))))
    #sys.stderr.write("esp=%s\n"%(str(len(esp))))
    #exit(1)

    return crds,pts,esp


def RunGauEsp(atn,crds,charge,mult,basename,nproc):
    import subprocess
    fh = open(basename+".com","w")
    fh.write("""%%MEM=2000MB
%%NPROC=%i

#P HF/6-31G* SCF(Conver=6) NoSymm Test 
   Pop=mk IOp(6/33=2) GFInput GFPrint 

RESP

%i %i
"""%(nproc,charge,mult))
    for i in range(len(atn)):
        fh.write("%4s %12.8f %12.8f %12.8f\n"%(str(atn[i]),crds[0+i*3],crds[1+i*3],crds[2+i*3]))
    fh.write("\n\n\n")
    fh.close()
    print("# Running g09 < %s.com > %s.log"%(basename,basename))
    subprocess.call("g09 < %s.com > %s.log"%(basename,basename), shell=True)


def ReadGauOutput(fname):
    atn=[]
    crds=[]
    charge=0
    mult=1
    fh=open(fname,"r")
    arc=""
    for line in fh:
        if '1\\1\\' in line:
            arc=line.strip()
            for line in fh:
                arc += line.strip()
                if '\\\\@' in arc:
                    break
    secs = arc.split("\\\\")

    try:
        data   = [ sub.split(",") for sub in secs[3].split("\\") ]
        charge = int( data[0][0] )
        mult   = int( data[0][1] )
        for i in range( len(data)-1 ):
            atn.append( data[i+1][0] )
            if len(data[i+1]) == 5:
                crds.append( float(data[i+1][2]) )
                crds.append( float(data[i+1][3]) )
                crds.append( float(data[i+1][4]) )
            else:
                crds.append( float(data[i+1][1]) )
                crds.append( float(data[i+1][2]) )
                crds.append( float(data[i+1][3]) )
    except:
        print("Could not process gaussian file '%s'"%(fname))
        print("This is the archive:")
        for i,sec in enumerate(secs):
            subs = sec.split("\\")
            for j,sub in enumerate(subs):
                vals = sub.split(",")
                print("%2i %2i %s"%(i,j,str(vals)))
    
    return atn,crds,charge,mult



def ReadOrMakeGauEsp( mdout, nproc=4 ):
    import os
    crds,pts,esp = ReadGauEsp(mdout)
    if len(esp) == 0:
        ele,crds,charge,mult = ReadGauOutput(mdout)
        espmdout = mdout.replace(".log","").replace(".out","") + ".resp"
        if not os.path.exists(espmdout+".log"):
            RunGauEsp(ele,crds,charge,mult,espmdout,nproc)
        crds,pts,esp = ReadGauEsp(espmdout+".log")
        if len(esp) == 0:
            raise Exception("Failed to generate ESP using g09 for %s"%(mdout))
    return crds,pts,esp



def WriteFitSh(base):
    sh = open("%s.resp.sh"%(base),"w")
    sh.write("""#!/bin/bash
out=%s.resp
cat <<EOF > ${out}.qwt
1
1.

EOF

resp -O -i ${out}.inp -o ${out}.out -p ${out}.punch -t ${out}.qout -w ${out}.qwt -e ${out}.esp

rm -f ${out}.punch ${out}.qout ${out}.qwt ${out}.esp

"""%(base))
    sh.close()
    subprocess.call(["bash","%s.resp.sh"%(base)])


def GetAvgStdDevAndErr(x):
    import math
    n = len(x)
    if n == 0:
        return 0,0,0
    avg = sum(x)/n
    var = 0.
    for a in x:
        var += math.pow(a-avg,2)
    if n > 1:
        var = var / ( n-1. )
    stddev = math.sqrt( var )
    stderr = math.sqrt( var / n )
    return avg,stddev,stderr

def WriteArray8(fh,arr):
    for istart in range(0,len(arr),16):
        for ioff in range(16):
            i = istart + ioff
            if i >= len(arr):
                break
            if ioff == 0:
                pass
            fh.write("%5i"%(arr[i]))
        fh.write("\n")


def GetResidueNameFromAtomIdx(parm,iat,unique_residues):
#    rname = parm.atoms[iat].residue.name[0].upper()
#    if len(parm.atoms[iat].residue.name) > 1 and rname == "D":
#        print  parm.atoms[iat].residue.name
#        rname = parm.atoms[iat].residue.name[1].upper()

        
#    for r in [("C5","C"),("G3","G"),("AFK","A"),("afk","A"),("aqm","A"),("cqm","C"),("GFK","G"),("gfk","G")]:
#        rname = rname.replace( r[0], r[1] )
    #if unique_residues:
    #    name = "%s_%i"%(parm.atoms[iat].residue.name,parm.atoms[iat].residue.idx+1)
    #else:
    name = parm.atoms[iat].residue.name
    return name

class ComboFH(object):
    def __init__(self,fh1,fh2):
        self.fh1=fh1
        self.fh2=fh2
    def write(self,a):
        self.fh1.write(a)
        self.fh2.write(a)
    def close(self):
        self.fh1.close()
        self.fh2.close()

def GetAtomsBondedToIdx(parm,idx):
    cats=[]
    for bond in parm.bonds:
        if bond.atom1.idx == idx:
            cats.append(bond.atom2.idx)
        elif bond.atom2.idx == idx:
            cats.append(bond.atom1.idx)
    return cats

def GetEquivHydrogens(parm,sele):
    hpairs = []
    for hvy in sele:
        if parm.atoms[hvy].atomic_number > 1:
            hbnds=[]
            bnds = GetAtomsBondedToIdx(parm,hvy)
            for bat in bnds:
                if parm.atoms[bat].atomic_number == 1:
                    hbnds.append(bat)
            if len(hbnds) > 1:
                hbnds.sort()
                hpairs.append( hbnds )
    return hpairs

def GetEquivNonbridgeOxygens(parm,sele):
    hpairs = []
    for hvy in sele:
        if parm.atoms[hvy].atomic_number == 15:
            hbnds=[]
            bnds = GetAtomsBondedToIdx(parm,hvy)
            for bat in bnds:
                if parm.atoms[bat].type == "O2" or parm.atoms[bat].type == "o2":
                    hbnds.append(bat)
            if len(hbnds) > 1:
                hbnds.sort()
                hpairs.append( hbnds )
    return hpairs


def ReadRespCharges(fname):
    qs = []
    fh = open( fname, "r" )
    for line in fh:
        if "no.  At.no." in line:
            for line in fh:
                c = line.strip().split()
                if len(c) == 6:
                    qs.append( float(c[3]) )
                else:
                    break
            break
    return qs


def ReadNextRespCharges(fh):
    import re
    import sys
    #prog = re.compile(r"^ {1}([ \d]{4})([ \d]{4}) {5}([ \d\.\-]{10}) {5}([ \d\.\-]{10})([ \d]{7})([ \d\.\-]{15})([ \d\.\-]{12})")
    prog = re.compile(r"^ {1}([ \d]{4})([ \d]{4}) {5}([ \d\.\-]{10}) {5}([ \d\.\-]{10})([ \d]{7})([ \d\.\-]{15})")

    qs = []
    for line in fh:
        #sys.stderr.write("%s\n"%(line))
        result = prog.match(line)
        if result is not None:
            #print "match ",result.group(4)
            qs.append( float( result.group(4) ) )
            for line in fh:
                #print line
                result = prog.match(line)
                if result is not None:
                    #print "match ",result.group(4)
                    qs.append( float( result.group(4) ) )
                else:
                    break
            break
    return qs



class RespMol(object):
    def __init__(self,state,mdout,unique_residues):
        self.state = state
        self.mdout = mdout
        self.atoms = []
        for fitidx in range(self.state.nquant):
            parmidx = self.state.fit2parm_map[ fitidx ]
            rname = GetResidueNameFromAtomIdx(self.state.parm,parmidx,unique_residues)
            ridx  = self.state.parm.atoms[parmidx].residue.idx+1
            aname = self.state.parm.atoms[parmidx].name
            if unique_residues:
                self.atoms.append( "%s%i:%s"%( rname, ridx, aname ) )
            else:
                self.atoms.append( "%s:%s"%( rname, aname ) )
        for linkidx in range(self.state.nlink):
            fitidx = linkidx + self.state.nquant
            parmidx = self.state.fit2parm_map[ fitidx ]
            rname = GetResidueNameFromAtomIdx(self.state.parm,parmidx,unique_residues)
            ridx  = self.state.parm.atoms[parmidx].residue.idx+1
            aname = self.state.parm.atoms[parmidx].name
            if unique_residues:
                self.atoms.append( "%s%i:%s:%i"%( rname, ridx, aname, linkidx+1 ) )
            else:
                self.atoms.append( "%s:%s:%i"%( rname, aname, linkidx+1 ) )


    def find_fitidx(self,name):
        fitidx=-1
        if name in self.atoms:
            fitidx = self.atoms.index( name )
        return fitidx

    def get_selected_atoms( self, sele ):
        parmidxs = parmutils.GetSelectedAtomIndices(self.state.parm,sele)
        atoms = []
        for parmidx in parmidxs:
            for fitidx in self.state.parm2fit_map[parmidx]:
                atoms.append( self.atoms[ fitidx ] )
        return atoms


class IntermolEquiv(object):
    def __init__( self, equiv_mask, unique_residues ):
        self.equiv_mask = equiv_mask
        self.mols = []
        self.state_map = ddict( str )
        self.states = []
        self.unique_residues=unique_residues

    def append( self, state, mdout ):
        mol = RespMol( state, mdout, self.unique_residues )
        istate = len(self.states)
        if state.base in self.state_map:
            istate = self.state_map[ state.base ]
        else:
            self.state_map[ state.base ] = istate
            self.states.append( [] )
        mol.istate = istate
        self.states[ istate ].append( mol )
        self.mols.append( mol )

    def get_unique_atoms(self,mask="@*"):
        uatoms = []
        if mask is not None:
            for mol in self.mols:
                atoms = mol.get_selected_atoms( mask )
                for atom in atoms:
                    if atom not in uatoms:
                        uatoms.append(atom)
        return uatoms

    def get_unique_atoms_from_state(self,state,mask="@*"):
        uatoms = []
        if mask is not None:
            for mol in state:
                atoms = mol.get_selected_atoms( mask )
                for atom in atoms:
                    if atom not in uatoms:
                        uatoms.append(atom)
        return uatoms
    

    def print_equiv( self, fh ):
        uatoms = self.get_unique_atoms( mask=self.equiv_mask )
        for u in uatoms:
            idxs = []
            for imol,mol in enumerate(self.mols):
                idx = mol.find_fitidx(u)
                if idx >= 0:
                    idxs.append( imol+1 )
                    idxs.append( idx+1 )
            if len(idxs) > 0:
                fh.write("%5i\n"%( len(idxs)//2 ) )
                WriteArray8(fh,idxs)
        for state in self.states:
            satoms = self.get_unique_atoms_from_state( state, "@*" )
            satoms = [ x for x in satoms if (x not in uatoms) ]
            for s in satoms:
                idxs = []
                for mol in state:
                    imol = self.mols.index(mol)
                    idx = mol.find_fitidx(s)
                    if idx >= 0:
                        idxs.append( imol+1 )
                        idxs.append( idx+1 )
                if len(idxs) > 0:
                    fh.write("%5i\n"%( len(idxs)//2 ) )
                    WriteArray8(fh,idxs)
        fh.write("\n")


class EndState(object):
    def __init__(self,parmfile,rstfiles,comp,ires,qmmask=None,theory="PBE0",basis="6-31G*",maxgrad=1.e-6,etol=1.e-4,fitgasphase=False):
        self.backbone_restraint=False
        self.sugar_restraint=False
        self.nucleobase_restraint=False
        self.multifit=False
        self.equiv_hydrogens=True
        self.equiv_nonbridge=True
        self.equiv_atoms=[]
        
        self.theory=theory
        self.basis=basis
        self.maxgrad=maxgrad
        self.etol=etol
        self.fitgasphase=fitgasphase
        
        self.parmfile = parmfile
        self.rstfiles = sorted(rstfiles)
        self.comp = comp
        self.ires = ires

        has_g09 = False
        has_hfdf = False
        for rst in rstfiles:
            if ".log" in rst:
                has_g09 = True
            elif ".rst7" in rst:
                has_hfdf = True
        if has_g09 and has_hfdf:
            raise Exception("EndState called with gaussian .log files and amber .rst7 files -- must use one or the other")
        
                
        if has_g09:
            self.parm = parmutils.OpenParm( parmfile, xyz=None )
        else:
            self.parm = parmutils.OpenParm( parmfile, xyz=rstfiles[0] )
        self.mmcharges = [ a.charge for a in self.parm.atoms ]
        import os.path
        topdir,parmfile = os.path.split( self.parmfile )
        self.topdir = topdir
        self.base = "mod.%i.%s"%(self.ires,parmfile.replace(".parm7",""))
        self.modified_parmfile = os.path.join(self.topdir,self.base + ".parm7")
        self.gasphase_parmfile = os.path.join(self.topdir,"gas.%i.%s.parm7"%(self.ires,parmfile.replace(".parm7","")))


        if qmmask is None:
            self.qmmask = "(:%i&!@P,OP*,*5'*)|(:%i&@P,OP*,*5'*)"%(ires,ires+1)
        else:
            self.qmmask = qmmask

        self.fsys = parmutils.FragmentedSys( self.parm, self.comp )
        self.fsys.add_fragment( self.qmmask, coef0=1, coef1=1, method="HF" )
        self.fsys.sort()
        #self.fsys.check_overlaps()
        #self.fsys.redistribute_residue_charges()
        #if not os.path.isfile( self.modified_parmfile ):
        #    parmutils.WriteParm( self.fsys.parmobj, self.modified_parmfile )

        self.frag = self.fsys.frags[0]

        self.nquant = len( self.frag.atomsel )
        linkatom_pairs = self.frag.GetLinkPairs()
        self.nlink = len( linkatom_pairs )
        self.nquant_nlink = self.nquant + self.nlink

        self.fit2parm_map = []
        self.fit_atomic_numbers = []
        for i,a in enumerate( self.frag.atomsel ):
            self.fit2parm_map.append( a )
            self.fit_atomic_numbers.append( self.parm.atoms[a].atomic_number )
        for ip in range( self.nlink ):
            iqm = self.frag.atomsel.index( linkatom_pairs[ip][0] )
            imm = self.nquant_nlink - self.nlink + ip
            self.fit2parm_map.append( linkatom_pairs[ip][0] )
            self.fit_atomic_numbers.append( 1 )

        self.parm2fit_map = ddict( list )
        for imm in range( self.nquant_nlink ):
            iqm = self.fit2parm_map[imm]
            self.parm2fit_map[iqm].append(imm)

        self.grp_restraints = []

        
    def set_equiv_atoms(self,sele):
        idxs = parmutils.GetSelectedAtomIndices(self.parm,sele)
        self.equiv_atoms.append(idxs)
        
        
    def add_group_restraint(self,cmask,q=None):
        if q is None:
            grp = parmutils.GetSelectedAtomIndices(self.parm,cmask)
            grp_charges = [ self.parm.atoms[a].charge for a in grp ]
            q  = sum( grp_charges )
        self.grp_restraints.append( (cmask,q) )
        

    def apply_backbone_restraint(self,value=True):
        self.backbone_restraint=value


    def apply_sugar_restraint(self,value=True):
        self.sugar_restraint=value


    def apply_nucleobase_restraint(self,value=True):
        self.nucleobase_restraint=value


    def multimolecule_fit(self,value=True):
        self.multifit = value


    def write_mdin(self):

        has_g09=False
        for rstfile in self.rstfiles:
            if ".log" in rstfile:
                has_g09=True
        if has_g09:
            return
            
        self.fsys.check_overlaps()
        self.fsys.redistribute_residue_charges()
        if not os.path.isfile( self.modified_parmfile ):
            parmutils.SaveParm( self.fsys.parmobj, self.modified_parmfile )
        if self.fitgasphase:
            if not os.path.isfile( self.gasphase_parmfile ):
                gp = parmutils.CopyParm( self.fsys.parmobj )
                qmidxs = parmutils.GetSelectedAtomIndices( gp, self.fsys.get_noshake_selection() )
                for a in gp.atoms:
                    if a.idx not in qmidxs:
                        a.charge = 0
                parmutils.SaveParm( gp, self.gasphase_parmfile )


            
        for rstfile in self.rstfiles:
            base="%s.%s"%(self.base,rstfile.replace(".rst7",""))
            mdin = mdinutils.Mdin()
            mdin.SetBaseName( base )
            mdin.Set_NVT()
            mdin.Set_Restart(False)
            mdin.Set_PrintFreq(1)
            mdin.cntrl["nstlim"] = 0
            mdin.cntrl["xpol_c"] = 0
            mdin.cntrl["ifqnt"]=1
            mdin.cntrl["ig"] = -1
            mdin.cntrl["ntf"] = 1
            mdin.cntrl["ntc"] = 2
            mdin.cntrl["ntwx"] = 0
            mdin.cntrl["ntwr"] = 0
            mdin.cntrl["noshakemask"] = '"' + self.fsys.get_noshake_selection() + '"'
            mdin.Set_QMMM_PBE0( qmmask='"'+self.fsys.get_noshake_selection()+'"', qmcharge=self.frag.qmcharge )
            mdin.qmmm["hfdf_theory"] = "'%s'"%(self.theory)
            mdin.qmmm["hfdf_basis"] = "'%s'"%(self.basis)
            mdin.qmmm["verbosity"] = 2
            mdin.PARM7 = self.modified_parmfile


            opt = copy.deepcopy(mdin)

            if self.fitgasphase:
                mdin.PARM7 = self.gasphase_parmfile
            
            opt.SetBaseName( "opt."+base )
            opt.qmmm["verbosity"] = 0
            opt.Set_DLFIND_Minimize()
            opt.dlfind["tol"] = self.maxgrad
            opt.dlfind["tole"] = self.etol

            opt.dlfind["active"]='"' + self.fsys.get_noshake_selection() + '"'
            opt.cntrl["imin"]   = 1  # minimize
            opt.cntrl["ntmin"]  = 5  # read &dlfind
            opt.cntrl["ntx"]    = 1

            opt.CRD7 = rstfile
            opt.WriteMdin()

            mdin.CRD7 = opt.RST7
            mdin.WriteMdin()

            sfile = self.comp.open( "%s.slurm"%(base) )
            sfile.write("if [ ! -e %s ]; then\n\n"%(mdin.MDOUT))
            sfile.write("%s %s %s\n\n"%(self.comp.mpirun,opt.EXE,opt.CmdString()))
            sfile.write("%s %s %s\n\n"%(self.comp.mpirun,mdin.EXE,mdin.CmdString()))
            sfile.write("rm -f %s %s %s %s %s %s opt.%s*xyz\n\n"%(opt.MDOUT,opt.MDINFO,opt.NC,mdin.MDINFO,mdin.RST7,mdin.NC,base))
            sfile.write("fi\n")
            sfile.close()
    

    def clear_charge_data(self):
        self.charge_data = ddict( list )

                
    def read_next_resp(self,fh):
        import sys
        if True:
            fit_charges = ReadNextRespCharges( fh )
            #print fit_charges
            qm_charges = []
            mm_charges = []
            for i,a in enumerate( self.frag.atomsel ):
                qqm  = 0.
                for imm in self.parm2fit_map[a]:
                    #sys.stderr.write("%s %s\n"%( str(a),str(imm) ))
                    qqm += fit_charges[imm]
                qm_charges.append( qqm )
                mm_charges.append( self.parm.atoms[a].charge )
            dq = ( sum(mm_charges) - sum(qm_charges) ) / len( mm_charges )
            qm_charges = [ q + dq for q in qm_charges ]
            for i,a in enumerate(self.frag.atomsel):
                self.charge_data[a].append( qm_charges[i] )

                
    def read_respfile(self):
        import os.path
        self.charge_data = ddict( list )
        mdouts = self.get_mdouts()
        if self.multifit:
            out = open(os.path.join(self.topdir,"%s.resp.out"%(self.base)),"r")
            #print "%s.resp.out"%(self.base)
            for mdout in mdouts:
                self.read_next_resp(out)
        else:
            for mdout in mdouts:
                base=mdout.replace(".mdout","")
                base=base.replace(".log","") # g09
                out = open("%s.resp.out"%(base),"r")
                #print "%s.resp.out"%(base)
                self.read_next_resp(out)

    def preserve_residue_charges_by_shifting(self):
        for res in self.parm.residues:
            mmq = 0.
            q = 0.
            n = 0
            for a in res.atoms:
                mmq += self.mmcharges[a.idx]
                if a.idx in self.frag.atomsel:
                    avg,std,err = GetAvgStdDevAndErr(self.charge_data[a.idx])
                    q += avg
                    if a.epsilon > 0.001:
                        n += 1
                else:
                    q += self.mmcharges[a.idx]
            if n > 0:
                dq = (mmq-q)/n
                for a in res.atoms:
                    if a.idx in self.frag.atomsel:
                        if a.epsilon > 0.001:
                            for i in range(len(self.charge_data[a.idx])):
                                self.charge_data[a.idx][i] += dq


    def preserve_mm_charges_by_shifting(self,mm_mask):
        reset_sele = parmutils.GetSelectedAtomIndices(self.parm,mm_mask)
        for res in self.parm.residues:
            mmq = 0.
            q   = 0.
            n   = 0
            for a in res.atoms:
                mmq += self.mmcharges[a.idx]
                if a.idx in self.frag.atomsel:
                    if a.idx in reset_sele:
                        self.charge_data[a.idx] = [ self.mmcharges[a.idx] ]
                        avg = self.mmcharges[a.idx]
                    else:
                        avg,std,err = GetAvgStdDevAndErr(self.charge_data[a.idx])
                    q += avg
                    if a.epsilon > 0.001 and a.idx not in reset_sele:
                        n += 1
                else:
                    q += self.mmcharges[a.idx]
            if n > 0:
                dq = (mmq-q)/n
                for a in res.atoms:
                    if a.idx in self.frag.atomsel:
                        if a.epsilon > 0.001 and a.idx not in reset_sele:
                            for i in range(len(self.charge_data[a.idx])):
                                self.charge_data[a.idx][i] += dq


                                
                            
    def print_resp(self,prefix="",fh=sys.stdout):
        for a in self.frag.atomsel:
            q,stddev,stderr = GetAvgStdDevAndErr(self.charge_data[a])
            mm = self.mmcharges[a]
            #res = GetResidueNameFromAtomIdx(self.parm,a,self.unique_residues)
            res = "%3i"%(self.parm.atoms[a].residue.idx+1)
            atm = self.parm.atoms[a].name
            dct = "QQS[\"%s\"][%s][\"%s\"]"%(prefix,res,atm)
            fh.write("    %-27s = %10.6f # %4s %10.6f %4i %8.4f %8.4f\n"%(dct,q,self.parm.atoms[a].residue.name,mm,len(self.charge_data[a]),stderr,stddev))
        fh.write("\n")


    def get_mdouts(self):
        mdouts = []
        for rstfile in self.rstfiles:
            topdir,rst = os.path.split( rstfile.replace(".rst7","") )
            base="%s.%s"%(self.base,rst)
            mdout=os.path.join(topdir,"%s.mdout"%(base))
            if os.path.isfile( mdout ):
                mdouts.append( mdout )
            elif os.path.isfile( os.path.join(topdir,"%s.log"%(base)) ): # g09
                mdouts.append( os.path.join(topdir,"%s.log"%(base)) ) # g09
            elif ".log" in rstfile: # g09
                mdouts.append( rstfile )
        return mdouts

    
    def perform_fit(self,unique_residues=True):
        if self.multifit:
            equiv = IntermolEquiv(None,unique_residues)
            mdouts = self.get_mdouts()
            body = self.get_resp_body()
            nmols = len(self.get_mdouts())

            inp = open(os.path.join(self.topdir,"%s.resp.inp"%(self.base)),"w")
            esp = open(os.path.join(self.topdir,"%s.resp.esp"%(self.base)),"w")
            inp.write( self.get_resp_header() )
            for imol,mdout in enumerate(mdouts):
                inp.write(body)
                equiv.append( self, mdout )
                self.append_esp( esp, mdout )
            esp.close()
            inp.write("\n")
            equiv.print_equiv( inp )
            inp.write("\n\n")
            inp.close()

            print("# EndState.perform_fit writing multifit %s.resp.inp"%(os.path.join(self.topdir,self.base)))
            WriteFitSh( os.path.join(self.topdir,self.base) )
        else:
            mdouts = self.get_mdouts()
            body = self.get_resp_body()
            header = self.get_resp_header()
            for mdout in mdouts:
                base=mdout.replace(".mdout","")
                base=base.replace(".log","") # g09
                
                print("# EndState.perform_fit writing singlefit %s.resp.inp"%(base))
                inp = open("%s.resp.inp"%(base),"w")
                inp.write( header )
                inp.write( body )
                inp.close()
                esp = open("%s.resp.esp"%(base),"w")
                self.append_esp( esp, mdout )
                esp.close()
                
                WriteFitSh( base )
        self.read_respfile()

    def get_resp_header(self):
        if self.multifit:
           nmol = len( self.get_mdouts() )
        else:
           nmol = 1
        return "title\n &cntrl inopt=0 ioutopt=0 iqopt=1 ihfree=1 irstrnt=1 iunits=1 qwt=0.0005 nmol=%i &end\n"%(nmol)


    def get_resp_body(self):
        import copy
        fh = StringIO()
        fh.write("%10.5f\n"%(1.0))
        fh.write(" molecule\n")
        fh.write("%5i%5i\n"%(self.frag.qmcharge,self.nquant_nlink))

        equiv_mask=[0]*self.nquant_nlink
        equiv_grps = copy.deepcopy(self.equiv_atoms)
        if self.equiv_hydrogens:
            equiv_grps += GetEquivHydrogens(self.parm,self.frag.atomsel)
        if self.equiv_nonbridge:
            equiv_grps += GetEquivNonbridgeOxygens(self.parm,self.frag.atomsel)
        for grp in equiv_grps:
            equiv_atms=[]
            for parmidx in grp:
                for fitidx in self.parm2fit_map[ parmidx ]:
                    equiv_atms.append(fitidx)
            equiv_atms.sort()
            first_atm = equiv_atms[0]+1
            for idx in equiv_atms[1:]:
                equiv_mask[idx] = first_atm

        for z,eq in zip( self.fit_atomic_numbers , equiv_mask ):
            fh.write("%5i%5i\n"%(z,eq))

        pmask = ":%i&(@P,OP*,*5'*)"%(self.ires+1)
        smask = ":%i&(@*')&!(@*5'*)"%(self.ires)
        bmask = ":%i&!(@P,OP*,*')"%(self.ires)

        cons_grps = []
        if self.backbone_restraint:
           cons_grps.append( pmask )
        if self.sugar_restraint:
           cons_grps.append( smask )
        if self.nucleobase_restraint:
           cons_grps.append( bmask )

        for cmask in cons_grps:
            grp = parmutils.GetSelectedAtomIndices(self.parm,cmask)
            grp_charges = [ self.parm.atoms[a].charge for a in grp ]
            grp_charge  = sum( grp_charges )
            grp_fit_idxs = []
            for parmidx in grp:
                for fitidx in self.parm2fit_map[ parmidx ]:
                    grp_fit_idxs.append( fitidx )
 
            fh.write("%5i%10.5f\n"%(len(grp_fit_idxs),grp_charge))
            arr=[]
            for idx in grp_fit_idxs:
                arr.append(1)
                arr.append(idx+1)
            WriteArray8(fh,arr)

        for cmask,grp_charge in self.grp_restraints:
            grp = parmutils.GetSelectedAtomIndices(self.parm,cmask)
            grp_fit_idxs = []
            for parmidx in grp:
                for fitidx in self.parm2fit_map[ parmidx ]:
                    grp_fit_idxs.append( fitidx )

            ##
            fh.write("%5i%10.5f\n"%(len(grp_fit_idxs),grp_charge))
            arr=[]
            for idx in grp_fit_idxs:
                arr.append(1)
                arr.append(idx+1)
            WriteArray8(fh,arr)
            
        fh.write("\n")

        return fh.getvalue()

    def append_esp(self,fh,mdout):
        out = open(mdout,"r")
        if ".mdout" in mdout:
            found_something = False
            for line in out:
                if "RESPESP " in line:
                    found_something = True
                    fh.write(line.replace("RESPESP ",""))
                    
            if not found_something:
                raise Exception("Could not find RESPESP data in %s"%(mdout))
        elif ".log" in mdout: # g09
            crds,pts,vals = ReadOrMakeGauEsp( mdout, max(4,self.comp.num_nodes) )
            nat=len(crds)//3
            npt=len(vals)
            fh.write("%5i%5i\n"%(nat,npt))
            au = 1.88972613373
            for i in range(nat):
                fh.write("%17s%16.7f%16.7f%16.7f\n"%("",au*crds[0+i*3],au*crds[1+i*3],au*crds[2+i*3]))
            for i in range(npt):
                fh.write(" %16.7f%16.7f%16.7f%16.7f\n"%(vals[i],au*pts[0+i*3],au*pts[1+i*3],au*pts[2+i*3]))

            
    def append_esps(self,fh):
        for mdout in self.get_mdouts():
            self.append_esp(fh,mdout)




class ResidueResp(object):
    def __init__(self,comp,ires,theory="PBE0",basis="6-31G*",maxgrad=1.e-6,etol=0.0001,fitgasphase=False):
        self.comp = comp
        self.ires = ires
        self.states = []
        self.base = "multistate.%i"%(self.ires)
        self.theory=theory
        self.basis=basis
        self.maxgrad=maxgrad
        self.etol=etol
        self.fitgasphase = fitgasphase
        self.multifit=False

    def add_state(self,prefix,parmfile,rstfiles,qmmask=None):
        self.states.append( EndState( parmfile, rstfiles, self.comp, self.ires, qmmask=qmmask, theory=self.theory,basis=self.basis,maxgrad=self.maxgrad,etol=self.etol,fitgasphase=self.fitgasphase ) )
        self.states[-1].prefix = prefix

    def clear_charge_data(self):
        for state in self.states:
            state.clear_charge_data()

    def write_mdin(self):
        for state in self.states:
            state.write_mdin()

    def read_respfile(self):
        if self.multifit:
            self.clear_charge_data()
            out = open("%s.resp.out"%(self.base),"r")
            for state in self.states:
                mdouts = state.get_mdouts()
                for mdout in mdouts:
                    state.read_next_resp( out )
                

    def preserve_residue_charges_by_shifting(self):
        for state in self.states:
            state.preserve_residue_charges_by_shifting()
              
    def preserve_mm_charges_by_shifting(self,mm_mask):
        for state in self.states:
            state.preserve_mm_charges_by_shifting(mm_mask)
                
    def print_resp(self,fh=sys.stdout):
        for state in self.states:
            state.print_resp(state.prefix,fh=fh)

    def apply_backbone_restraint(self,value=True):
        for state in self.states:
            state.backbone_restraint=value

    def apply_sugar_restraint(self,value=True):
        for state in self.states:
            state.sugar_restraint=value

    def apply_nucleobase_restraint(self,value=True):
        for state in self.states:
            state.nucleobase_restraint=value

    def apply_equiv_hydrogens(self,value=True):
        for state in self.states:
            state.equiv_hydrogens=value

    def apply_equiv_nonbridge(self,value=True):
        for state in self.states:
            state.equiv_nonbridge=value
            
    def multimolecule_fit(self,value=True):
        self.multifit=value
        for state in self.states:
            state.multifit = value

    def perform_fit(self,equiv_mask="@P,OP*,*'",unique_residues=True):
        if self.multifit:
            equiv = IntermolEquiv(equiv_mask,unique_residues)
            nmols = 0
            for state in self.states:
                nmols += len(state.get_mdouts())
                
            print("# ResidueResp.perform_fit writing multifit %s.resp.inp"%(self.base))
            
            inp = open("%s.resp.inp"%(self.base),"w")
            esp = open("%s.resp.esp"%(self.base),"w")
            inp.write( "title\n &cntrl inopt=0 ioutopt=0 iqopt=1 ihfree=1 irstrnt=1 iunits=1 qwt=0.0005 nmol=%i &end\n"%(nmols) )

            for state in self.states:
                mdouts = state.get_mdouts()
                body = state.get_resp_body()
                for mdout in mdouts:
                    inp.write(body)
                    equiv.append( state, mdout )
                    state.append_esp( esp, mdout )
            esp.close()
            inp.write("\n")
            equiv.print_equiv( inp )
            inp.write("\n\n")
            inp.close()
            WriteFitSh( self.base )
            self.read_respfile()
        else:
            equiv = IntermolEquiv(equiv_mask,unique_residues)
            for state in self.states:
                state.perform_fit(unique_residues=unique_residues)
                state.read_respfile()
                mdouts = state.get_mdouts()
                for mdout in [ mdouts[0] ]:
                    equiv.append( state, mdout )

            # for each unique atom, create a mega-array of all charges
            # and then reset each dependent atom to the mega-array
            uatoms = equiv.get_unique_atoms( mask=equiv.equiv_mask )
            for uat in uatoms:
                qdat = []
                for mol in equiv.mols:
                    for fitidx,name in enumerate(mol.atoms):
                        #parmidx = mol.state.fit2parm_map[ fitidx ]
                        if name == uat:
                            #print name,mol.state.base,len(mol.state.charge_data[fitidx])
                            qdat.extend( mol.state.charge_data[fitidx] )
                for mol in equiv.mols:
                    for fitidx,name in enumerate(mol.atoms):
                        if name == uat:
                            #print "changing from ",len(mol.state.charge_data[fitidx])," to ",len(qdat)
                            #del mol.state.charge_data[fitidx][:]
                            #mol.state.charge_data[fitidx].extend( qdat )
                            self.states[ mol.istate ].charge_data[fitidx] = qdat
                    


if __name__ == "__main__":
        
    #####################################################################################
    #####################################################################################


    comp = parmutils.CALIBURN(numnodes=2)
    comp.amberhome="${HOME}/caliburn/hfdfamber"

    labels=["GHAH","GHAX","GXAH","GXAX"]



    sys.stdout.write("""def OnlyBaseA7andBaseG45change_nocon():
    QQS = ddict( lambda: ddict( lambda: ddict(float) ) )
    # U6:  common-frames:         common-states: :6|(:7&@P,OP*,*')   charge-con: none
    # A7:  common-frames: OP1=OP2 common-states: :8|(:7&@P,OP*,*')   charge-con: none
    # G45: common-frames: OP1=OP2 common-states: :46|(:45&@P,OP*,*') charge-con: none

""")

    for ires in [6,7,45]:
        fit = ResidueResp(comp,ires)
        for prefix in labels:
            rsts = glob.glob(prefix+".00*.rst7")
            fit.add_state( prefix, prefix+".parm7",rsts )
        if False:
            fit.write_mdin()
        elif True:
            fit.multimolecule_fit()
            #fit.apply_backbone_restraint()
            if ires == 6:
                fit.apply_equiv_nonbridge(False)
                fit.perform_fit(":6|(:7&@P,OP*,*')")
            elif ires == 7:
                fit.perform_fit(":8|(:7&@P,OP*,*')")
            elif ires == 45:
                fit.perform_fit(":46|(:45&@P,OP*,*')")
            else:
                fit.perform_fit("@P,OP*,*'")
            fit.read_respfile()
            fit.preserve_residue_charges_by_shifting()
            fit.print_resp()

    sys.stdout.write("""
    return QQS

""")




    sys.stdout.write("""def OnlyResA7andBaseG45change_nocon():
    QQS = ddict( lambda: ddict( lambda: ddict(float) ) )
    # U6:  common-frames:         common-states: :6   charge-con: none
    # A7:  common-frames: OP1=OP2 common-states: :8   charge-con: none
    # G45: common-frames: OP1=OP2 common-states: :46|(:45&@P,OP*,*') charge-con: none

""")

    for ires in [6,7,45]:
        fit = ResidueResp(comp,ires)
        for prefix in labels:
            rsts = glob.glob(prefix+".00*.rst7")
            fit.add_state( prefix, prefix+".parm7",rsts )
        if False:
            fit.write_mdin()
        elif True:
            fit.multimolecule_fit()
            #fit.apply_backbone_restraint()
            if ires == 6:
                fit.apply_equiv_nonbridge(False)
                fit.perform_fit(":6")
            elif ires == 7:
                fit.perform_fit(":8")
            elif ires == 45:
                fit.perform_fit(":46|(:45&@P,OP*,*')")
            else:
                fit.perform_fit("@P,OP*,*'")
            fit.read_respfile()
            fit.preserve_residue_charges_by_shifting()
            fit.print_resp()

    sys.stdout.write("""
    return QQS

""")




    sys.stdout.write("""def OnlyResA7andResG45change_nocon():
    QQS = ddict( lambda: ddict( lambda: ddict(float) ) )
    # U6:  common-frames:         common-states: :6   charge-con: none
    # A7:  common-frames: OP1=OP2 common-states: :8   charge-con: none
    # G45: common-frames: OP1=OP2 common-states: :46 charge-con: none

""")

    for ires in [6,7,45]:
        fit = ResidueResp(comp,ires)
        for prefix in labels:
            rsts = glob.glob(prefix+".00*.rst7")
            fit.add_state( prefix, prefix+".parm7",rsts )
        if False:
            fit.write_mdin()
        elif True:
            fit.multimolecule_fit()
            #fit.apply_backbone_restraint()
            if ires == 6:
                fit.apply_equiv_nonbridge(False)
                fit.perform_fit(":6")
            elif ires == 7:
                fit.perform_fit(":8")
            elif ires == 45:
                fit.perform_fit(":46")
            else:
                fit.perform_fit("@P,OP*,*'")
            fit.read_respfile()
            fit.preserve_residue_charges_by_shifting()
            fit.print_resp()

    sys.stdout.write("""
    return QQS

""")




    
