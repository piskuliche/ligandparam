import parmed
import math
from collections import defaultdict as ddict
import random

from ligandparam.multiresp import mdinutils
from ligandparam.multiresp import functions

class Computer(object):
    """ Base class for computer objects """
    def __init__(self,numnodes):
        """ Initialize the computer object 
        
        Parameters
        ----------
        numnodes : int
            The number of nodes to use
        """
        self.mpirun=""
        self.num_nodes = numnodes
        self.exclude = ""
        self.cores_per_node = 1
        self.amberhome="${AMBERHOME}"
        self.gpu=False
        self.array = []
        self.array_max_running = 1

        
    def get_array(self):
        alist = list(sorted(set(self.array)))
        rs=[]
        if len(alist) > 0:
            rs = [ (alist[0],alist[0]) ]
            for a in alist[1:]:
                if a == rs[-1][1]+1:
                    rs[-1] = ( rs[-1][0], a )
                else:
                    rs.append( (a,a) )
        sarr = []
        for r in rs:
            if r[0] != r[1]:
                sarr.append( "%i-%i"%(r[0],r[1]) )
            else:
                sarr.append( "%i"%(r[0]) )
        return ",".join(sarr)


    def write_array(self,fh):
        """ Write the array to the file handle 
        
        Parameters
        ----------
        fh : file handle
            The file handle to write to
        """
        if len(self.array) > 0:
            fh.write("#SBATCH --array=%s"%(self.get_array()))
            if self.array_max_running > 0:
                fh.write("%%%i"%(self.array_max_running))
            fh.write("\n")


    def unset_amberhome(self):
        """ Unset the amberhome variable """
        self.amberhome = None

    def set_exclude(self,x):
        """ Set the exclude variable """
        self.exclude = x

    def use_gpu(self,x=True):
        """ Use the gpu 
        
        Parameters
        ----------
        x : bool, optional
            If True, use the gpu
        """
        self.gpu=x

    def get_num_cores(self):
        """ Get the number of cores """
        return self.num_nodes * self.cores_per_node



class BASH(Computer):
    """ Generates a bash script """
    def __init__(self,numnodes):
        """ Initialize the bash object

        Parameters
        ----------
        numnodes : int
            The number of nodes to use
        """

        Computer.__init__(self,numnodes)
        self.mpirun="mpirun -n %i"%(self.get_num_cores())
                    
    def open(self,fname):
        """ Open the bash file

        Parameters
        ----------
        fname : str
            The name of the file to open
        """
        fh = open(fname,"w")
        fh.write("#!/bin/bash\n\n")
        return fh
    

def OpenParm( fname, xyz=None ):
    """ Open a file with parmed.
    
    Parameters
    ----------
    fname : str
        The name of the file to open
    xyz : str, optional
        The name of the xyz file to open

    Returns
    -------
    parmed object
        The parmed object
    """
    import parmed

    try:
        from parmed.constants import IFBOX
    except:
        from parmed.constants import PrmtopPointers
        IFBOX = PrmtopPointers.IFBOX

    if ".mol2" in fname:
        param = parmed.load_file( fname, structure=True )
    else:
        param = parmed.load_file( fname, xyz=xyz )
        if xyz is not None:
            if ".rst7" in xyz:
                param.load_rst7(xyz)
    if param.box is not None:
        if abs(param.box[3]-109.471219)<1.e-4 and \
           abs(param.box[4]-109.471219)<1.e-4 and \
           abs(param.box[5]-109.471219)<1.e-4:
            param.parm_data["POINTERS"][IFBOX]=2
            param.pointers["IFBOX"]=2
    return param

def CopyParm( parm ):
    """ Copy the parmed object 
    
    Parameters
    ----------
    parm : parmed object
        The parmed object to copy
    
    Returns
    -------
    parmed object
        The copied parmed object
    """
    import copy
    try:
        parm.remake_parm()
    except:
        pass
    p = copy.copy( parm )
    p.coordinates = copy.copy( parm.coordinates )
    p.box = copy.copy( parm.box )
    try:
        p.hasbox = copy.copy( parm.hasbox )
    except:
        p.hasbox = False
    return p

def MakeUniqueBondParams( p, xlist, scale=1.0 ):
    """ Make unique bond parameters
    
    Parameters
    ----------
    p : parmed object
        The parmed object
    xlist : list
        The list of bonds
    scale : float, optional
        The scale factor. Default is 1.0
    """
    from collections import defaultdict as ddict
    byidx = ddict( list )
    for x in xlist:
        byidx[ x.type.idx ].append( x )
    for idx in byidx:
        x = byidx[idx][0].type
        p.bond_types.append( parmed.BondType( x.k*scale, x.req, p.bond_types ) )
        for x in byidx[idx]:
            #help(p.bonds[0])
            #p.bonds[ x.idx ].type = p.bond_types[-1]
            x.type = p.bond_types[-1]

def MakeUniqueAngleParams( p, xlist, scale=1.0 ):
    """ Make unique angle parameters
    
    Parameters
    ----------
    p : parmed object
        The parmed object
    xlist : list
        The list of angles
    scale : float, optional
        The scale factor. Default is 1.0
    """
    from collections import defaultdict as ddict
    byidx = ddict( list )
    for x in xlist:
        byidx[ x.type.idx ].append( x )
    for idx in byidx:
        x = byidx[idx][0].type
        p.angle_types.append( parmed.AngleType( x.k*scale, x.theteq, p.angle_types ) )
        for x in byidx[idx]:
            #p.angles[ x.idx ].type = p.angle_types[-1]
            x.type = p.angle_types[-1]

def MakeUniqueDihedralParams( p, xlist, scale=1.0 ):
    """ Make unique dihedral parameters
    
    Parameters
    ----------
    p : parmed object
        The parmed object
    xlist : list
        The list of dihedrals
    scale : float, optional
        The scale factor. Default is 1.0"""
    from collections import defaultdict as ddict
    byidx = ddict( list )
    for x in xlist:
        byidx[ x.type.idx ].append( x )
    for idx in byidx:
        x = byidx[idx][0].type
        p.dihedral_types.append( parmed.DihedralType( x.phi_k*scale, x.per, x.phase, x.scee, x.scnb, p.dihedral_types ) )
        for x in byidx[idx]:
            #p.dihedrals[ x.idx ].type = p.dihedral_types[-1]
            x.type = p.dihedral_types[-1]

                
def GetSelectedAtomIndices(param,maskstr):
    """ Get the selected atom indices
    
    Parameters
    ----------
    param : parmed object
        The parmed object
    maskstr : str
        The mask string
    
    """
    #param = parmed.load_file(parmfile)
    #mask = parmed.amber.mask.AmberMask( param, maskstr )
    #aidxs = mask.Selected()
    #for aidx in aidxs:
    #    atom = param.atoms[aidx]
    #    res  = atom.residue
    sele = []
    if len(maskstr) > 0:
        newmaskstr = maskstr.replace("@0","!@*")
        sele = [ param.atoms[i].idx for i in parmed.amber.mask.AmberMask( param, newmaskstr ).Selected() ]
    return sele


def GetSelectedResidueIndices(param,maskstr):
    """ Get the selected residue indices
    
    Parameters
    ----------
    param : parmed object
        The parmed object
    maskstr : str
        The mask string
    
    """
    a = GetSelectedAtomIndices(param,maskstr)
    b = list(set([ param.atoms[c].residue.idx for c in a ]))
    b.sort()
    return b

def ListToSelection(atomlist):
    """ Convert a list to a selection
    
    Parameters
    ----------
    atomlist : list
        The list of atoms
    
    Returns
    -------
    str
        The selection
    """
    alist = list(sorted(set(atomlist)))
    rs=[]
    if len(alist) > 0:
        rs = [ (alist[0],alist[0]) ]
        for a in alist[1:]:
            if a == rs[-1][1]+1:
                rs[-1] = ( rs[-1][0], a )
            else:
                rs.append( (a,a) )
    sarr = []
    for r in rs:
        if r[0] != r[1]:
            sarr.append( "%i-%i"%(r[0]+1,r[1]+1) )
        else:
            sarr.append( "%i"%(r[0]+1) )
    sele = "@0"
    if len(sarr) > 0:
        sele = "@" + ",".join(sarr)
    return sele

class Fragment(object):
    """ A fragment """
    def __init__(self,parmobj,ambmask,coef0=None,coef1=None,method="AM1D"):
        """ Initialize the fragment 
        
        Parameters
        ----------
        parmobj : parmed object
            The parmed object
        ambmask : str
            The ambmask
        coef0 : float, optional
            The coefficient for lambda=0. Default is None
        coef1 : float, optional
            The coefficient for lambda=1. Default is None
        method : str, optional
            The method. Default is AM1D
        
        """
        self.ambmask = ambmask
        self.parmobj = parmobj
        self.parmfilename = None
        # Check if coef0 and coef1 are None
        if coef0 is None:
            if coef1 is None:
                raise Exception("Fragment %s must have a coefficient"%(ambmask))
            else:
                self.coef0 = coef1
                self.coef1 = coef1
        else:
            if coef1 is None:
                self.coef0=coef0
                self.coef1=coef0
            else:
                self.coef0=coef0
                self.coef1=coef1

        self.method = method
        self.atomsel = GetSelectedAtomIndices(parmobj,ambmask)
        self.nat = len(self.atomsel)
        self.mmcharge = 0.

        for iat in self.atomsel:
            self.mmcharge += parmobj.atoms[iat].charge

        self.qmcharge = 0
        resmmcharges = self.get_selected_mmcharge_from_each_touched_residue()

        for residx in resmmcharges:
            self.qmcharge += int(round(resmmcharges[residx]))
        if self.qmcharge != int(round(self.mmcharge)):
            print("WARNING: qm charge (%i) of fragment '%s' is suspect (sum of mm charges: %.4f)"%(self.qmcharge,self.ambmask,self.mmcharge))
        for ires in self.get_touched_residues():
            self.funkify_residue_name(ires)

    def funkify_residue_name(self,ires):
        """ Funkify the residue name for the given residue index
        
        Parameters
        ----------
        ires : int
            The residue index
            
        """
        origname = self.parmobj.residues[ires].name
        newname = self.get_funkified_residue_name(origname)
        if newname is not None:
            self.parmobj.residues[ires].name = newname

    def get_funkified_residue_name(self,origname):
        """ Get the funkified residue name
        
        Parameters
        ----------
        origname : str
            The original name
        
        Returns
        -------
        str
            The funkified residue name
        """
        if origname[0].islower():
            return origname
        else:
            resnames = [ res.name for res in self.parmobj.residues ]
            for charoffset in range(20):
                firstchar = chr( ord(origname[0].lower()) + charoffset ).lower()
                for i in range(100):
                    name = "%s%02i"%(firstchar,i)
                    if name in resnames:
                        continue
                    return name
       
    def get_coef(self,lam=0):
        """ Get the coefficient
        
        Parameters
        ----------
        lam : int, optional
            The lambda value. Default is 0
        
        Returns
        -------
        float
            The coefficient
        """
        return self.coef0 + lam*(self.coef1-self.coef0)
 
    def get_touched_residues(self):
        """ Get the touched residues """
        residues = []
        for a in self.atomsel:
            res = self.parmobj.atoms[a].residue.idx
            residues.append( res )
        return list(set(residues))

    def get_selected_mmcharge_from_each_touched_residue(self):
        """ Get the selected mm charge from each touched residue """
        residues = self.get_touched_residues()
        charges = ddict(float)
        for idx in residues:
            res = self.parmobj.residues[idx]
            resselecharge = 0.
            for atom in res.atoms:
                if atom.idx in self.atomsel:
                    resselecharge += atom.charge
            charges[idx] = resselecharge
        return charges

    
    def redistribute_residue_charges(self):
        """ Redistribute the residue charges """
        TOL=1.e-4
        residues = self.get_touched_residues()

        #print self.mmcharge, self.qmcharge,abs( self.mmcharge - self.qmcharge )
        if abs( self.mmcharge - self.qmcharge ) > 0.001:
            #print "changing charges"
            
            initresq = 0.
            for idx in residues:
                res = self.parmobj.residues[idx]
                for atom in res.atoms:
                    initresq += atom.charge

            chargeable=[]
            for idx in residues:
                res = self.parmobj.residues[idx]
                resmmcharge = 0.
                selmmcharge = 0.
                num_nondummy_atoms = 0
                num_sele_atoms = 0
                for atom in res.atoms:
                    resmmcharge += atom.charge
                    if atom.idx in self.atomsel:
                        selmmcharge += atom.charge
                        if abs(atom.charge) > TOL:
                            num_sele_atoms += 1
                    elif abs(atom.charge) > TOL:
                        num_nondummy_atoms += 1
                        chargeable.append( atom.idx )
                selqmcharge = int(round(selmmcharge))
                seldq=0
                remsq=0
                if num_sele_atoms > 0:
                    seldq = (selqmcharge - selmmcharge) / num_sele_atoms
                    if num_nondummy_atoms > 0:
                        remdq = (selmmcharge - selqmcharge) / num_nondummy_atoms
                if num_sele_atoms > 0 and num_nondummy_atoms > 0:
                    for atom in res.atoms:
                        if abs(atom.charge) > TOL:
                            if atom.idx in self.atomsel:
                                atom.charge += seldq
                            else:
                                atom.charge += remdq
                        else:
                            pass

            postresq = 0.
            for idx in residues:
                res = self.parmobj.residues[idx]
                for atom in res.atoms:
                    postresq += atom.charge
            if len(chargeable) > 0:
                dq = (initresq-postresq) / len(chargeable)
                for idx in chargeable:
                    atom = self.parmobj.atoms[idx]
                    #print "was",atom.charge,
                    atom.charge += dq
                    #print "now",atom.charge
            postresq = 0.
            for idx in residues:
                res = self.parmobj.residues[idx]
                for atom in res.atoms:
                    postresq += atom.charge
            if abs(postresq-initresq) > 0.0001:
                print("WAS NOT ABLE TO PRESERVE CHARGE FOR FRAGMENT:",self.ambmask)


            
        cats=self.GetConnectionAtoms()
        for cat in cats:
             if self.parmobj.atoms[cat].residue.idx not in residues:
                 self.funkify_residue_name(self.parmobj.atoms[cat].residue.idx)
                 # If the connection atom is not one of the touched residues, then
                 # we don't want to monkey with the charge because we don't want
                 # to apply the change to ALL guanine residues (if the connection atom
                 # belonged to a GUA).  We can, however, monkey with it if it is one
                 # of the touched residues... well, unless we touch more than 1
                 # residue with the same name... hmmm. The tleap.sh script should
                 # manually set mol.residx.name CHARGE for every atom in the system
                 # or the residue names should be reset in a more intelligent way
             ats = []
             for at in self.GetAtomsBondedToIdx(cat):
                 if at in cats:
                    continue
                 if at in self.atomsel:
                    continue
                 if self.parmobj.atoms[at].residue.idx != self.parmobj.atoms[cat].residue.idx:
                    continue
                 ats.append(at)
             #print "cat=",cat,"ats=",ats
             if len(ats) > 0:
                 dq = self.parmobj.atoms[cat].charge / len(ats)
                 for at in ats:
                     self.parmobj.atoms[at].charge += dq
                 self.parmobj.atoms[cat].charge = 0.


    def GetMMBoundaryTerms(self):
        """ Get the MM boundary terms """
        linkpairs = self.GetLinkPairs()
        bonds  = []
        angles = []
        dihedrals  = []
        for x in self.parmobj.bonds:
            a = (x.atom1.idx,x.atom2.idx)
            b = (x.atom2.idx,x.atom1.idx)
            if a in linkpairs or b in linkpairs:
                bonds.append(x)
                for y in self.parmobj.angles:
                    a=(y.atom1.idx,y.atom2.idx)
                    b=(y.atom2.idx,y.atom3.idx)
                    c=(y.atom2.idx,y.atom1.idx)
                    d=(y.atom3.idx,y.atom2.idx)
                    if a in linkpairs or b in linkpairs or c in linkpairs or d in linkpairs:
                        angles.append( y )
                        for z in self.parmobj.dihedrals:
                            a=(z.atom1.idx,z.atom2.idx)
                            b=(z.atom2.idx,z.atom3.idx)
                            c=(z.atom3.idx,z.atom4.idx)
                            d=(z.atom2.idx,z.atom1.idx)
                            e=(z.atom3.idx,z.atom2.idx)
                            f=(z.atom4.idx,z.atom3.idx)
                            if a in linkpairs or b in linkpairs \
                               or c in linkpairs or d in linkpairs \
                               or e in linkpairs or f in linkpairs:
                                dihedrals.append( z )
        return bonds,angles,dihedrals



    def GetConnectionAtoms(self):
        """ Get the connection atoms """
        cats=[]
        for bond in self.parmobj.bonds:
            if bond.atom1.idx in self.atomsel:
                if bond.atom2.idx not in self.atomsel:
                    cats.append(bond.atom2.idx)
            elif bond.atom2.idx in self.atomsel:
                if bond.atom1.idx not in self.atomsel:
                    cats.append(bond.atom1.idx)
        #print "connection atoms:",cats
        return cats

    def GetLinkPairs(self):
        """ Get the link pairs """
        cats=[]
        for bond in self.parmobj.bonds:
            if bond.atom1.idx in self.atomsel:
                if bond.atom2.idx not in self.atomsel:
                    cats.append( (bond.atom1.idx,bond.atom2.idx) )
            elif bond.atom2.idx in self.atomsel:
                if bond.atom1.idx not in self.atomsel:
                    cats.append( (bond.atom2.idx,bond.atom1.idx) )
        cats.sort(key=lambda x: x[0])
        #print "connection atoms:",cats
        return cats


    def GetAtomsBondedToIdx(self,idx):
        """ Get the atoms bonded to the index
        
        Parameters
        ----------
        idx : int
            The index
        """
        
        cats=[]
        for bond in self.parmobj.bonds:
            if bond.atom1.idx == idx:
                cats.append(bond.atom2.idx)
            elif bond.atom2.idx == idx:
                cats.append(bond.atom1.idx)
        return cats

class FragmentedSys(object):
    """ A fragmented system """

    def __init__(self,parmobj,compobj):
        """ Initialize the fragmented system 
        
        Parameters
        ----------
        parmobj : parmed object
            The parmed object
        compobj : BASH
            The computer object
        """
        #self.parmfile = parmfile
        #self.rstfile = rstfile
        #self.parmobj = parmutils.OpenParm( parmfile, rstfile )
        self.parmobj = CopyParm( parmobj )
        self.compobj = compobj
        self.frags = []
        # self.nve = False
        # self.restart = True
        # self.cut=10.
        # self.nstlim=2000
        # self.ntpr=100
        # self.ntwr=100
        # self.ntwx=100
        self.mdin_templates = ddict()
        self.mdin = mdinutils.Mdin()
        self.parmfilename = None
        self.mm_parmfilename = None
        self.disang = None
        self.ig = random.randint(1,650663)
        
    def set_mm_parm(self,fname):
        """ Set the MM parm file"""
        self.mm_parmfilename = fname
        
    def add_fragment(self,ambmask,coef0,coef1=None,method="AM1D"):
        """ Add a fragment
        
        Parameters
        ----------
        ambmask : str
            The ambmask
        coef0 : float
            The coefficient for lambda=0
        coef1 : float, optional
            The coefficient for lambda=1. Default is None
        method : str, optional
            The method. Default is AM1D
        
        """
        self.frags.append( Fragment(self.parmobj,ambmask,coef0=coef0,coef1=coef1,method=method) )

    def GetQMAtoms(self):
        """ Get the QM atoms 
        
        Returns
        -------
        list
            The QM atoms
        """
        qmatoms=[]
        for f in self.frags:
            qmatoms.extend( f.atomsel )
        qmatoms = list(set(qmatoms))
        return qmatoms

    def GetMMTermsForQMRegion(self):
        """ Get the MM terms for the QM region 
        
        Returns
        -------
        bonds : list
            The bonds
        angles : list
            The angles
        dihes : list
            The dihedrals
        """

        bonds=[]
        angles=[]
        dihes=[]
        qmatoms=[]
        for f in self.frags:
            qmatoms.extend( f.atomsel )
        qmatoms = list(set(qmatoms))
        for x in self.parmobj.bonds:
            if x.atom1.idx in qmatoms or x.atom2.idx in qmatoms:
                bonds.append( x )
        for y in self.parmobj.angles:
            if y.atom1.idx in qmatoms or y.atom2.idx in qmatoms or y.atom3.idx in qmatoms:
                angles.append( y )
        for y in self.parmobj.dihedrals:
            if y.atom1.idx in qmatoms or y.atom2.idx in qmatoms or y.atom3.idx in qmatoms or y.atom4.idx in qmatoms:
                dihes.append( y )
        if len(bonds) > 0:
            bonds  = list(set(bonds))
        if len(angles) > 0:
            angles = list(set(angles))
        if len(dihes) > 0:
            dihes  = list(set(dihes))
        return bonds,angles,dihes

    
    def GetMMTermsForQMRegionAsDict(self):
        p = ddict( lambda:ddict( list ) )
        b,a,d = self.GetMMTermsForQMRegion()
        for x in b:
            p["bonds"][x.type.idx].append(x)
        for x in a:
            p["angles"][x.type.idx].append(x)
        for x in d:
            p["dihedrals"][x.type.idx].append(x)
        return p
    
    def MakeNewMMBoundaryTerms(self):
        """ Make new MM boundary terms """
        # p = self.GetMMBoundaryTermsAsDict()
        # for itype in p["bonds"]:
        #     a1=[ x.atom1.idx for x in p["bonds"][itype] ]
        #     a2=[ x.atom2.idx for x in p["bonds"][itype] ]
        #     AddNewBondType(self.parmobj,a1,a2)
        # for itype in p["angles"]:
        #     a1=[ x.atom1.idx for x in p["angles"][itype] ]
        #     a2=[ x.atom2.idx for x in p["angles"][itype] ]
        #     a3=[ x.atom3.idx for x in p["angles"][itype] ]
        #     AddNewAngleType(self.parmobj,a1,a2,a3)
        b,a,d = self.GetMMBoundaryTerms()
        MakeUniqueBondParams( self.parmobj, b )
        MakeUniqueAngleParams( self.parmobj, a )
        MakeUniqueDihedralParams( self.parmobj, d )

        sel = []
        for frag in self.frags:
            sel += frag.atomsel
        sel = list(set(sel))
        
        from collections import defaultdict as ddict
        ljs = ddict( list )
        for s in sel:
            ljs[ self.parmobj.atoms[s].nb_idx ].append( s )
        for nbidx in ljs:
            a = self.parmobj.atoms[ ljs[nbidx][0] ]
            rad   = a.rmin
            eps   = a.epsilon
            rad14=0
            eps14=0
            try:
                rad14 = a.rmin14
                eps14 = a.epsilon14
            except:
                pass
            from parmed.tools.addljtype import AddLJType
            print(nbidx)
            sel = [ 0 ]*len(self.parmobj.atoms)
            for a in ljs[nbidx]:
                sel[a] = 1
            AddLJType( self.parmobj, sel, rad, eps, rad14, eps14 )

        for i in range(len(self.parmobj.atoms)):
            self.parmobj.atoms[i].nb_idx = self.parmobj.parm_data['ATOM_TYPE_INDEX'][i]
#        for i in range( len(self.parmobj.atoms) ):
#            a = self.parmobj.atoms[i]
#            print "%5i %5i"%( i+1, self.parmobj.parm_data['ATOM_TYPE_INDEX'][i] )


    
    def add_mm(self):
        """ Add the MM fragment """
        self.frags = [ f for f in self.frags if f.method != "MM" ]
        s0 = 0.
        s1 = 0.
        for f in self.frags:
            s0 += f.get_coef(0)
            s1 += f.get_coef(1)
        self.frags.append( Fragment(self.parmobj,":0",coef0=(1-s0),coef1=(1-s1),method="MM") )

    def sort(self):
        """ Sort the fragments """
        qm  = []
        sqm = []
        mm  = []
        for f in self.frags:
            if f.method == "MM":
                mm.append(f)
            elif f.method == "AM1D" or f.method == "DFTB":
                sqm.append(f)
            else:
                qm.append(f)
        qm  = sorted( qm,  key=lambda x: x.nat, reverse=True )
        sqm = sorted( sqm, key=lambda x: x.nat, reverse=True )
        mm  = sorted( mm,  key=lambda x: x.nat, reverse=True )
        self.frags = qm + sqm + mm
        self.redistribute_cores()

    def get_noshake_selection(self):
        """ Get the no shake selection """
        alist = []
        for f in self.frags:
            alist.extend( f.atomsel )
        return ListToSelection(alist)
        
    def check_overlaps(self):
        """ Check for overlaps
        
        Raises
        ------
        Exception
            If an overlap is found
        """
        TOL=1.e-8
        methods = set( [ f.method for f in self.frags ] )
        
        found_error = False
        for method in methods:
            if "MM" in method or "mm" in method:
                continue
            csum0 = ddict(int)
            csum1 = ddict(int)
            for f in self.frags:
                if f.method == method:
                    for a in f.atomsel:
                        csum0[a] += f.coef0
                        csum1[a] += f.coef1
            for a in sorted(csum0):
                ok = False
                if (abs(csum0[a]) < TOL or abs(csum0[a]-1.) < TOL) and (abs(csum1[a]) < TOL or abs(csum1[a]-1.) < TOL):
                    ok = True
                if not ok:
                    found_error = True
                    print("Sum of %8s fragment coefs invalid for atom %7i : (lam0: %13.4e, lam1:%13.4e)"%\
                        (method,a+1,csum0[a],csum1[a]))
        if found_error:
            raise Exception("Unaccounted fragment overlap exists")

        s0 = 0.
        s1 = 0.
        for f in self.frags:
            s0 += f.coef0
            s1 += f.coef1
        if abs(s0-1.) > 1.e-8:
            raise Exception("Sum of lambda=0 coefficients != 1 (%.8f)"%(s0))
        if abs(s1-1.) > 1.e-8:
            raise Exception("Sum of lambda=1 coefficients != 1 (%.8f)"%(s1))

    def get_coefs(self,lam=0):
        """ Get the coefficients

        Parameters
        ----------
        lam : int, optional
            The lambda value. Default is 0
        """ 

        c=[]
        for f in self.frags:
            c.append( f.get_coef(lam) )
        return c

    def get_energy(self,fragenes,lam=0):
        """ Get the energy
        
        Parameters
        ----------
        fragenes : list
            The list of fragment energies
        lam : int, optional
            The lambda value. Default is 0
        """
        cs = self.get_coefs(lam)
        e = 0.
        for c,efrag in zip(cs,fragenes):
            e += c * efrag
        return e

    def get_dvdl(self,fragenes):
        """ Get the dvdl
        
        Parameters
        ----------
        fragenes : list
            The list of fragment energies
        
        Returns
        -------
        float
            The dvdl
        """
        return self.get_energy(fragenes,1.) - self.get_energy(fragenes,0.)

    def get_mbar(self,fragenes,nlam=11):
        """ Get the mbar
        
        Parameters
        ----------
        fragenes : list
            The list of fragment energies
        nlam : int, optional
            The number of lambdas. Default is 11
        
        Returns
        -------
        list
            The mbar
        """
        ene = []
        for i in range(nlam):
            lam = i/(nlam-1.)
            ene.append( self.get_energy(fragenes,lam) )
        return ene
    
    def get_dvdl_coefs(self):
        """ Get the dvdl coefficients """
        c0=self.get_coefs(0)
        c1=self.get_coefs(1)
        dc=[]
        for a,b in zip(c0,c1):
            dc.append( b-a )

    def redistribute_residue_charges(self):
        """ Redistribute the residue charges """
        self.sort()
        for f in reversed(self.frags):
            f.redistribute_residue_charges()


    def write_parm(self,parmname="frag.parm7",overwrite=True):
        """ Write the parameter file
        
        Parameters
        ----------
        parmname : str, optional
            The parameter name. Default is "frag.parm7"
        overwrite : bool, optional
            If True, overwrite the file. Default is True
        """
        import subprocess
        self.parmfilename = parmname
        aidxs = []
        for f in self.frags:
            aidxs.extend( f.atomsel )
        if len(aidxs) < 1:
            raise Exception("No fragments")

        base=parmname.replace(".parm7","")
        print("Writing %s.notsele.frcmod and %s.sele.frcmod"%(base,base))
        functions.WriteMaskedFrcmod(self.parmobj,aidxs,"%s.notsele.frcmod"%(base),"%s.sele.frcmod"%(base))

        print("Writing %s.pdb"%(base))
        parmed.tools.writeCoordinates(self.parmobj, "%s.pdb" % (base)).execute()
    
        print("Writing %s.lib"%(base))
        parmed.tools.writeOFF(self.parmobj, "%s.lib" % (base)).execute()

        print("Writing %s.sh"%(base))
        functions.WriteLeapSh \
            ("%s.sh"%(base),
             self.parmobj,
             ["%s.lib"%(base)],
             ["%s.notsele.frcmod"%(base),"%s.sele.frcmod"%(base)],
             "%s.pdb"%(base),
              base,overwrite=overwrite)
        print("Running tleap script %s.sh"%(base))
        subprocess.call("bash %s.sh"%(base),shell=True)


    def write_mm_optimization(self,reftraj,refmdin=None):
        import copy
        parmname    = self.parmfilename
        base        = parmname.replace(".parm7","")
        compobj     = copy.deepcopy( self.compobj )
        compobj.num_nodes = min( compobj.num_nodes, 1 )
        mdin        = copy.deepcopy( self.mdin )
        mdin.SetBaseName("trial")
        mdin.SetGroupSize(None)
        mdin.PARM7 = "trial.parm7"
        if refmdin is None:
           refmdin = "${REFTRAJ%.nc}.mdin"
        mdin.CRD7 = refmdin.replace(".mdin",".rst7")
        print("Writing %s.mmopt.slurm"%(base))
        slurm = compobj.open( "%s.mmopt.slurm"%(base) )
        slurm.write("function make_trial_parm()\n")
        slurm.write("{\n")
        slurm.write("REFBASE=\"$1\"\n")
        functions.WriteLeapSh \
            ("trial",
             self.parmobj,
             ["%s.lib"%("${REFBASE}")],
             ["%s.notsele.frcmod"%("${REFBASE}"),"${BASE}.frcmod"],
             "%s.pdb"%("${REFBASE}"),
             "trial",
             overwrite=True,
             fh=slurm)
        slurm.write("rm trial.rst7 trial.out trial.cmds\n")
        slurm.write("}\n\n\n\n")
        slurm.write("##############################################\n")
        slurm.write("REFBASE=%s\n"%(base))
        slurm.write("REFTRAJ=%s\n"%(reftraj))
        slurm.write("REFMDIN=%s\n"%(refmdin))
        slurm.write("##############################################\n")
        slurm.write("""
start=${REFBASE}.sele.frcmod
reqfiles=(${start} ${REFTRAJ} ${REFBASE}.lib ${REFBASE}.notsele.frcmod ${REFBASE}.pdb)
for f in ${reqfiles[@]}; do
    if [ ! -e "${f}" ]; then
        echo "Missing required file: ${f}"
        exit 1
    fi
done


if [ ! -e trial.mdin ]; then
    if [ -e ${REFMDIN} ]; then
        cp ${REFMDIN} trial.mdin
        sed -i 's|ifqnt *= *[0-9]*|ifqnt = 0|' trial.mdin
        sed -i 's|nmropt *= *[0-9]*|nmropt = 0|' trial.mdin
        sed -i 's|ntwx *= *[0-9]*|ntwx = 50|' trial.mdin
        sed -i 's|xpol_c|\!xpol_c|' trial.mdin
    else
        echo "trial.mdin not found"
        exit 1
    fi
fi

for iter in $(seq 10); do
    next=$(printf "trial.frcmod.%02i" ${iter})
    trial=$(printf "trial.frcmod.%02i" $(( ${iter} - 1 )))
    if [ -e "${next}" ]; then
        echo "${next} already exists; skipping"
        continue
    fi
    if [ "${iter}" == "1" ]; then
        cp ${start} ${trial}
    fi    
    if [ -e trial.frcmod ]; then
        echo "trial.frcmod exists; exiting"
        exit 1
    fi
    if [ ! -e ${trial} ]; then
        echo "File not found: ${trial}"
        exit 1
    fi

    nstlim=$(( ${iter} * 10000 ))
    sed -i "s|nstlim.*|nstlim = ${nstlim}|" trial.mdin


    echo "Making trial parm7"
    ln -s ${trial} trial.frcmod
    make_trial_parm "${REFBASE}"
    rm trial.frcmod


    echo "Running pmemd"\n""")
        slurm.write("    %s %s %s"%( compobj.mpirun, "pmemd.MPI", mdin.CmdString() ))
        slurm.write("""
    rm trial.rst7 trial.mdinfo trial.mdout
    
    
    echo "Generating next frcmod: ${next}"
    python2.7 `which parmutils-UpdateParamFromTraj.py` -p trial.parm7 --reftraj ${REFTRAJ} --trialtraj trial.nc --trialfrcmod ${trial} --newfrcmod ${next} > ${next}.out
    
done

if [ -e trial.nc ]; then
  rm -f trial.nc
fi

""")
        

        
    def has_ab_initio(self):
        """ Check if there is ab initio """
        has_abinitio = False
        for f in self.frags:
            if f.method == "AM1D" or f.method == "DFTB" or f.method == "MM":
                continue
            else:
                has_abinitio = True
                break
        return has_abinitio

        
    def redistribute_cores(self):
        """ Redistribute the cores """
        has_abinitio = self.has_ab_initio()
        ncores = self.compobj.get_num_cores()
        #print "has_abinitio?",has_abinitio,ncores
        for f in self.frags:
            f.ncores=1
        if not has_abinitio and ncores > 0:
           nfrag=len(self.frags)
           nrem = ncores - nfrag
           if nrem < 0:
               raise Exception("Must use at least %i cores"%(nfrag))
           i=0
           while nrem > 0:
              self.frags[ i % nfrag ].ncores += 1
              i += 1
              nrem -= 1
        elif ncores > 0:
            nfrag=len(self.frags)
            nrem = ncores - nfrag
            if nrem < 0:
                raise Exception("Must use at least %i cores"%(nfrag))
            wts=[0]*nfrag

            for i,f in enumerate(self.frags):
                wts[i] = 0
                if "MM" in f.method:
                    wts[i]=0
                elif "AM1D" in f.method or "DFTB" in f.method:
                    if not has_abinitio:
                        #wts[i] = math.pow(f.nat,2.25) 
                        wts[i] = 1. + f.nat + math.pow(f.nat/2.,1.2)
                else:
                    x0=1
                    x1=15
                    y0=0.0001
                    y1=1.
                    m = (y1-y0)/(x1-x0)
                    b = y0-m*x0
                    if f.nat < x1:
                        pref = m*f.nat+b
                    else:
                        pref = 1.0
                    wts[i] = pref*math.pow(f.nat,2.4)
            twt = sum(wts)
            tnc = 0
            for i,f in enumerate(self.frags):
                wts[i] = int(round(nrem * wts[i] / twt))
                f.ncores += wts[i]
                tnc += f.ncores
            for f in reversed(self.frags):
                if tnc > ncores:
                    if f.ncores > 1:
                        f.ncores -= 1
                        tnc -= 1
                else:
                    break
            for f in self.frags:
                if tnc < ncores:
                    f.ncores += 1
                    tnc += 1
                else:
                    break
                
    def write_mdin( self, prefix="frag", init="init", lam=0, init_from_same_rst=False, directory=None, dipout=False, same_ntpr=False ):
        """ Write the mdin file
        
        Parameters
        ----------
        prefix : str, optional
            The prefix. Default is "frag"
        init : str, optional
            The initialization. Default is "init"
        lam : int, optional
            The lambda value. Default is 0
        init_from_same_rst : bool, optional
            If True, initialize from the same restart. Default is False
        directory : str, optional
            The directory. Default is None
        dipout : bool, optional
            If True, output the dipole. Default is False
        same_ntpr : bool, optional
            If True, use the same ntpr. Default is False
        """
        import copy
        import os.path
        gfilename = "%s.%.8f.groupfile"%(prefix,lam)
        if directory is not None:
            gfile = open(os.path.join( directory, gfilename ),"w")
        else:
            gfile = open(gfilename,"w")
        for i,f in enumerate(self.frags):
            
            base="%s.%06i_%.8f"%(prefix,i+1,lam)
            if i+1 == len(self.frags):
                base = "%s.mm_%.8f"%(prefix,lam)
                
            if f.method in self.mdin_templates:
                mdin = copy.deepcopy(self.mdin_templates[f.method])
            else:
                mdin = copy.deepcopy(self.mdin)
                if f.method == "AM1D":
                    mdin.cntrl["ifqnt"]=1
                    mdin.Set_QMMM_AM1(qmmask='"'+f.ambmask+'"',qmcharge=f.qmcharge)
                    #mdin.qmmm["scfconv"] = 1.e-9
                    mdin.qmmm["diag_routine"] = 6
                elif f.method == "MM":
                    mdin.cntrl["ifqnt"]=0
                else:
                    mdin.cntrl["ifqnt"]=1
                    mdin.Set_QMMM_PBE0(qmmask='"'+f.ambmask+'"',qmcharge=f.qmcharge)
                    mdin.qmmm["hfdf_theory"] = "'%s'"%(f.method)

            if i+1 < len(self.frags):
                mdin.cntrl["ntwr"]=0
                mdin.cntrl["ntwx"]=0
                if not same_ntpr:
                    mdin.cntrl["ntpr"]=0


            mdin.title = f.ambmask
            mdin.cntrl["xpol_c"] = "%14.10f"%( f.get_coef(lam) )
            mdin.cntrl["ntf"] = 1
            #mdin.cntrl["ntc"] = 2
            mdin.cntrl["noshakemask"] = '"' + self.get_noshake_selection() + '"'
            mdin.cntrl["ig"] = self.ig
            mdin.SetGroupSize( f.ncores )
            mdin.SetBaseName( base )
            if self.disang is not None:
                mdin.DISANG=self.disang

            if dipout:
                mdin.DIPOUT = "%s.dipout"%( base )
                
            #mdin.CRD7 = "%s.%06i_%.8f.rst7"%(init,i+1,lam)
            #if i+1 == len(self.frags):
            mdin.CRD7 = "%s.mm_%.8f.rst7"%(init,lam)
            if mdin.cntrl["irest"] == 0 or init_from_same_rst:
                mdin.CRD7 = "%s.rst7"%(init)
            if self.parmfilename is not None:
                mdin.PARM7 = self.parmfilename
            if self.mm_parmfilename is not None and f.method == "MM":
                mdin.PARM7 = self.mm_parmfilename
            mdin.WriteMdin( directory=directory )
            gfile.write("%s\n"%( mdin.CmdString() ))

        f = "%s.%.8f.slurm"%(prefix,lam)
        if directory is not None:
            f = os.path.join( directory, f )
        sfile = self.compobj.open( f )
        sfile.write("%s %s -ng %i -groupfile %s\n\n"%(self.compobj.mpirun,mdin.EXE,len(self.frags),gfilename))





    def read_mdout( self, prefix="frag", lam=0, nlam=11 ):
        """ Read the mdout file
        
        Parameters
        ----------
        prefix : str, optional
            The prefix. Default is "frag"
        lam : int, optional
            The lambda value. Default is 0
        nlam : int, optional
            The number of lambdas. Default is 11

        Returns
        -------
        float
            The time step
        list
            The dvdl
        list
            The mbar
        """
        lams = [ i/(nlam-1.) for i in range(nlam) ]
        base = "%s.mm_%.8f"%(prefix,lam)
        fh=open(base+".mdout","r")
        dvdl=[]
        mbar=[]
        data=[]
        times=[]
        reading=False
        for line in fh:
            if not reading:
                if "A V E R A G E S" in line:
                    reading = False
                    break
                elif "TIME(PS)" in line:
                    reading = True
                    times.append( float(line.strip().split()[5]) )
                    data = []
            else:
                if "FragEne" in line:
                    data.append( float(line.strip().split()[2]) )
                elif "--------" in line:
                    dvdl.append( self.get_dvdl( data ) )
                    mbar.append( self.get_mbar( data, nlam ) )
                    reading=False
                    time=None
        if len(times) > 1:
            dt = times[1]-times[0]
        else:
            print("WARNING: %s.mdout is empty"%(base))
            dt=-1
            dvdl=[]
            mbar=[]
        return dt,dvdl,mbar


    def mdouts_to_pymbar( self, name="frag", nlam=11, datadir="mbar/data" ):
        """ Convert the mdouts to pymbar
        
        Parameters
        ----------
        name : str, optional
            The name. Default is "frag"
        nlam : int, optional
            The number of lambdas. Default is 11
        datadir : str, optional
            The data directory. Default is "mbar/data"
        
        Raises
        ------
        Exception
            If no files match the glob 'prod*.%s%s'
        """
        import os
        import glob
        
        lams = [ i/(nlam-1.) for i in range(nlam) ]

        try: 
            os.makedirs(datadir)
        except OSError:
            if not os.path.isdir(datadir):
                raise

        dvdl_data = ddict( lambda: ddict( float ) )
        efep_data = ddict( lambda: ddict( lambda: ddict( float ) ) )

        for lam in lams:
            post = ".mm_%.8f.mdout"%(lam)
            mdouts = sorted( glob.glob("prod*.%s%s"%(name,post)) )
            if len(mdouts) == 0:
                raise Exception("No files match glob 'prod*.%s%s'\n"%(name,post))
            dt = None
            dvdl = []
            mbar = []
            for mdout in mdouts:
                prefix = mdout.replace(post,"")
                mydt,dvdl_tmp,mbar_tmp = self.read_mdout(prefix=prefix,lam=lam,nlam=nlam)
                if mydt > 0:
                    dt=mydt
                    dvdl.extend(dvdl_tmp)
                    mbar.extend(mbar_tmp)
            if len(dvdl) > 5:
                g = functions.statisticalInefficiency( dvdl )
                print("dvdl_%s_%.8f.dat g=%.1f tau=%.3f"%(name,lam,g,0.5*(g-1.)*dt))
        
            for i,x in enumerate(dvdl):
                t = i * dt
                dvdl_data[lam][t] = x

            for iplam,plam in enumerate(lams):
                for i,v in enumerate(mbar):
                    t = i * dt
                    efep_data[lam][plam][t] = v[iplam]

        # ########################
        # alltimes = [ set(sorted(dvdl_data[0.6])) ]
        # alltimes.append( set(sorted(dvdl_data[0.8])) )
        # times = set.intersection( *alltimes )
        # lam=0.7
        # dvdl_data[lam] = ddict(float)
        # efep_data[lam] = ddict( lambda: ddict(float) )
        # for t in times:
        #     dvdl_data[lam][t] = 0.5*( dvdl_data[0.6][t] + dvdl_data[0.8][t] )
        #     for plam in lams:
        #         efep_data[lam][plam][t] = 0.5*( efep_data[0.6][plam][t] + efep_data[0.8][plam][t] )
        # ########################

        # prune the time series so all data has consistent times in the trajectory
        pruned_dvdl = ddict( lambda: ddict(float) )
        pruned_efep = ddict( lambda: ddict( lambda: ddict(float) ) )
        remove_time_gaps = False
        for tlam in dvdl_data:
            # get the unique set of times for this trajectory
            alltimes = [ set(sorted(efep_data[tlam][plam])) for plam in efep_data[tlam] ]
            alltimes.append( set(sorted(dvdl_data[tlam])) )
            times = set.intersection( *alltimes )
            
            for it,t in enumerate(sorted(times)):
                if remove_time_gaps:
                    time = times[0] + it * dt
                    pruned_dvdl[tlam][time] = dvdl_data[tlam][t]
                else:
                    pruned_dvdl[tlam][t] = dvdl_data[tlam][t]
                
            for plam in efep_data[tlam]:
                for it,t in enumerate(sorted(times)):
                    if remove_time_gaps:
                        time = times[0] + it * dt
                        pruned_efep[tlam][plam][time] = efep_data[tlam][plam][t]
                    else:
                        pruned_efep[tlam][plam][t] = efep_data[tlam][plam][t]

        dvdl_data = pruned_dvdl
        efep_data = pruned_efep

        for tlam in lams:
            dvdlfilename = "%s/dvdl_%s_%.8f.dat"%(datadir,name,tlam)
            fh = open(dvdlfilename,"w")
            for t in sorted(dvdl_data[tlam]):
                fh.write("%10.3f %23.14e\n"%(t,dvdl_data[tlam][t]))
            fh.close()

        for tlam in lams:
            for plam in lams:
                efepfilename = "%s/efep_%s_%.8f_%.8f.dat"%(datadir,name,tlam,plam)
                fh = open(efepfilename,"w")
                for t in sorted(efep_data[tlam][plam]):
                    fh.write("%10.3f %23.14e\n"%(t,efep_data[tlam][plam][t]))
                fh.close()

    