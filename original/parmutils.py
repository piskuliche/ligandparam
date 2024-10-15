#!/usr/bin/python

import parmed
import os
import math
from collections import defaultdict as ddict
import random



def OpenParm( fname, xyz=None ):
    import parmed
    #from parmed.constants import IFBOX
    
    try:
        from parmed.constants import IFBOX
    except:
        from parmed.constants import PrmtopPointers
        IFBOX = PrmtopPointers.IFBOX
        
    if ".mol2" in fname:
        param = parmed.load_file( fname, structure=True )
        #help(param)
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


def SaveParm( param, fname ):
    #from parmed.constants import IFBOX
    
    try:
        from parmed.constants import IFBOX
    except:
        from parmed.constants import PrmtopPointers
        IFBOX = PrmtopPointers.IFBOX
        
    for a in param.atoms:
        param.parm_data["CHARGE"][ a.idx ] = a.charge
    if param.box is not None:
       if abs(param.box[3]-109.471219)<1.e-4 and \
          abs(param.box[4]-109.471219)<1.e-4 and \
          abs(param.box[5]-109.471219)<1.e-4:
           param.parm_data["POINTERS"][IFBOX]=2
           param.pointers["IFBOX"]=2
    try:
        param.save( fname, overwrite=True )
    except:
        param.save( fname )

        
def SaveCrds( param, fname ):
    parmed.tools.writeCoordinates( param, fname ).execute()


        
def SaveLeapParm(parmobj,fname,aidxs=None):
    import subprocess
    base=fname.replace(".parm7","")
    
    if aidxs is None:
        print("Writing %s.frcmod"%(base))
        WriteFrcmod(parmobj,"%s.frcmod"%(base),uniqueparams=True)
    else:
        print("Writing %s.notsele.frcmod and %s.sele.frcmod"%(base,base))
        WriteMaskedFrcmod(parmobj,aidxs,"%s.notsele.frcmod"%(base),"%s.sele.frcmod"%(base))
        
    print("Writing %s.pdb"%(base))
    parmed.tools.writeCoordinates(parmobj, "%s.pdb"%(base)).execute()
    
    print("Writing %s.lib"%(base))
    parmed.tools.writeOFF(parmobj, "%s.lib"%(base)).execute()

    if aidxs is None:
        print("Writing %s.sh"%(base))
        WriteLeapSh \
            ("%s.sh"%(base),
             parmobj,
             ["%s.lib"%(base)],
             ["%s.frcmod"%(base)],
             "%s.pdb"%(base),
             base,overwrite=True)
    else:
        print("Writing %s.sh"%(base))
        WriteLeapSh \
            ("%s.sh"%(base),
             parmobj,
             ["%s.lib"%(base)],
             ["%s.notsele.frcmod"%(base),"%s.sele.frcmod"%(base)],
             "%s.pdb"%(base),
             base,overwrite=True)
        
    print("Running tleap script %s.sh"%(base))
    subprocess.call("bash %s.sh"%(base),shell=True)
        



def CopyParm( parm ):
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


def Strip( parm, mask ):
    p = CopyParm( parm )
    p.strip( "%s"%(mask) )
    return p


def Extract( parm, mask ):
    return Strip( parm, "!(%s)"%(mask) )


def Join( p, q ):
    import numpy
    import copy
    pq = p + q
    pq.coordinates = numpy.concatenate( (p.coordinates, q.coordinates) )
    pq.box = copy.copy( p.box )
    pq.hasbox = copy.copy( p.hasbox )
    return pq


def MakeUniqueBondParams( p, xlist, scale=1.0 ):
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

            
            
    
    
def FindDummyAtoms(param):
    dats=[]
    for a in param.atoms:
        if abs(a.charge) < 0.00001 and abs(a.epsilon) < 0.00001:
            dats.append( a.idx )
    return dats


def SaveParmWithoutDummyAtoms( param, fname ):
    import subprocess
    import os
    dats = FindDummyAtoms(param)
    SaveParm( param, fname )
    if len(dats) > 0:
        sele = ListToSelection(dats)
        ptrajinp = fname + ".ptraj.in"
        fh = open( ptrajinp, "w" )
        have_rst = False
        try:
            rst=fname.replace(".parm7","") + ".rst7"
            parmed.tools.writeCoordinates(param, rst).execute()
            fh.write("""
parm %s
trajin %s
strip %s
trajout %s
"""%( fname,rst,sele,rst ))
            have_rst = True
        except:
            pass

        fh.write("parm %s name sparm\n"%(fname))
        fh.write("parmstrip %s parm sparm\n"%(sele))
        fh.write("parmwrite out %s parm sparm\n"%(fname))
        fh.write("run\n")
        fh.write("quit\n\n")
        fh.close()
        subprocess.call("cpptraj -i %s"%(ptrajinp),shell=True)
        os.remove(ptrajinp)
        if have_rst:
            p = OpenParm( fname, xyz=rst )
            SaveParm( p, fname )

        

    
def GetSelectedAtomIndices(param,maskstr):
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
    a = GetSelectedAtomIndices(param,maskstr)
    b = list(set([ param.atoms[c].residue.idx for c in a ]))
    b.sort()
    return b

    
def ListToSelection(atomlist):
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




def WriteLeapSh(leapsh,param,lib,frcmod,pdb,base,fh=None,overwrite=False):
    """Writes a shell-script for running tleap

@param leapsh: name of shell-script to write
@param param: parm7 object (from parmed)
@param lib: list of off files to read
@param frcmod: list of frcmod files to read
@param pdb: pdb file to read
@param base: the output basnemae of the parm7 and rst7 files
@return none
    """
    if fh is None:
        fh = open(leapsh,"w")
        fh.write("#!/bin/bash\n\n")
    fh.write("BASE=\"%s\"\n\n"%(base))
    if not overwrite:
        fh.write("if [ -e \"%s.parm7\" ]; then BASE=\"%s.new\"; fi\n\n"%(base,base))

    for f in lib:
        fh.write("if grep -Eq 'A33|A55|G33|G55|C33|C55|U33|U55' %s; then\n"%(f))
        fh.write("   sed -i -e 's/A33/A3/' -e 's/A55/A5/' -e 's/G33/G3/' -e 's/G55/G5/'")
        fh.write(" -e 's/C33/C3/' -e 's/C55/C5/' -e 's/U33/U3/' -e 's/U55/U5/' %s\n"%(f))
        fh.write("fi\n")
    fh.write("echo \"Writing ${BASE}.parm7 and ${BASE}.rst7\"\n\n")
    fh.write("cat <<EOF > %s.cmds\n"%(leapsh))
    for f in lib:
        fh.write("loadOff %s\n"%(f))
    for f in frcmod:
        fh.write("loadAmberParams %s\n"%(f))

    if isinstance(pdb, str):
        pdbs = [ pdb ]
    else:
        pdbs = pdb
    combine = []
    for i,x in enumerate(pdbs):
        name = "x%i"%(i+1)
        fh.write("%s = loadPdb %s\n"%(name,x))
        combine.append( name )
    fh.write("x = combine { %s }\n"%(" ".join(combine)))
    
    if param.box is not None:
        fh.write("setbox x centers")
    fh.write("""
saveAmberParm x ${BASE}.parm7 ${BASE}.rst7
quit
EOF

tleap -s -f %s.cmds | grep -v "+---" | grep -v "+Currently" >> %s.out
    """%(leapsh,leapsh))
    
    if param.box is not None:
       fh.write("""
# Set the box dimensions
ChBox -X %.12f -Y %.12f -Z %.12f -al %.12f -bt %.12f -gm %.12f -c ${BASE}.rst7 -o ${BASE}.rst7.tmp; mv ${BASE}.rst7.tmp ${BASE}.rst7
"""%(param.box[0],param.box[1],param.box[2],param.box[3],param.box[4],param.box[5]))

       if abs(param.box[-1]-109.471219000000) < 1.e-4:
          fh.write("""
# Reset the ifbox flag
sed -i 's|0       0       1|0       0       2|' ${BASE}.parm7

""")


def GetLinkPairs(p,qmmask):
    if isinstance(qmmask,list):
        atomsel = qmmask
    else:
        atomsel = GetSelectedAtomIndices(p,qmmask)
    cats=[]
    for bond in p.bonds:
        if bond.atom1.idx in atomsel:
            if bond.atom2.idx not in atomsel:
                cats.append( (bond.atom1.idx,bond.atom2.idx) )
        elif bond.atom2.idx in atomsel:
            if bond.atom1.idx not in atomsel:
                cats.append( (bond.atom2.idx,bond.atom1.idx) )
    cats.sort(key=lambda x: x[0])
    return cats




def tisplit( parm, timask1, timask2, base="tisplit" ):
    import copy
    import os
    import subprocess
    
    parmed.tools.writeOFF( parm, "%s.lib"%(base) ).execute()
    WriteFrcmod( parm, "%s.frcmod"%(base), uniqueparams=False )

    p = CopyParm( parm )
    p.strip( "!(%s)"%( timask1 ) )
    parmed.tools.writeCoordinates(p, "%s.1.pdb"%(base)).execute()

    p = CopyParm( parm )
    p.strip( "!(%s)"%( timask2 ) )
    parmed.tools.writeCoordinates(p, "%s.2.pdb"%(base)).execute()

    p = CopyParm( parm )
    p.strip( "((%s)|(%s))"%( timask1, timask2 ) )
    parmed.tools.writeCoordinates(p, "%s.3.pdb"%(base)).execute()

    fh = open( "%s.sh"%(base), "w" )
    
    WriteLeapSh(
        "%s.sh"%(base),
        parm,
        [ "%s.lib"%(base) ],
        [ "%s.frcmod"%(base) ],
        [ "%s.1.pdb"%(base),
          "%s.1.pdb"%(base),
          "%s.3.pdb"%(base) ],
        "%s.1"%(base),
        fh=fh,
        overwrite=True )

    WriteLeapSh(
        "%s.sh"%(base),
        parm,
        [ "%s.lib"%(base) ],
        [ "%s.frcmod"%(base) ],
        [ "%s.2.pdb"%(base),
          "%s.2.pdb"%(base),
          "%s.3.pdb"%(base) ],
        "%s.2"%(base),
        fh=fh,
        overwrite=True )

    fh.close()

    subprocess.call("bash %s.sh"%(base),shell=True)

    os.remove( "%s.lib"%(base) )
    os.remove( "%s.frcmod"%(base) )
    os.remove( "%s.1.pdb"%(base) )
    os.remove( "%s.2.pdb"%(base) )
    os.remove( "%s.3.pdb"%(base) )
    os.remove( "%s.sh"%(base) )
    os.remove( "%s.sh.cmds"%(base) )
    os.remove( "%s.sh.out"%(base) )

    



def AddNewBondType(parm,sel1,sel2):
    from parmed.topologyobjects import BondType,AngleType,DihedralType
    from parmed.topologyobjects import Bond,Angle,Dihedral
    import copy
    
#    sel1 = parmed.amber.mask.AmberMask(parm,mask1).Selection()
#    sel2 = parmed.amber.mask.AmberMask(parm,mask2).Selection()
#    foo = GetSelectedAtomIndices(parm,mask1)
    if len(sel1) != len(sel2):
        raise Exception('parmutils.AddNewBondType: Each mask must select the same '
                        'number of atoms!')

    # If no atoms, nothing to do
    if len(sel1) == 0: return

    new_bnd_typ = None
    #atnum1, atnum2 = -1, -1
    # Loop through all of the selected atoms
    for atnum1,atnum2 in zip(sel1,sel2):
        # Collect the atoms involved
        #atnum1 = sel1.index(1, atnum1+1)
        #atnum2 = sel2.index(1, atnum2+1)
        atm1 = parm.atoms[atnum1]
        atm2 = parm.atoms[atnum2]
        # See if the bond exists in the first place, and if so, replace its
        # bond type with our new bond type (new_bnd)
        if atm2 in atm1.bond_partners and atm1 in atm2.bond_partners:
            for bond in atm1.bonds:
                if atm2 in bond:
                    # print "found %4s@%-4s - %4s@%-4s : %7.4f %7.4f"%\
                    #     ( atm1.residue.name,atm1.name,
                    #       atm2.residue.name,atm2.name,
                    #       bond.type.req,bond.type.k)
                    new_bnd_typ = copy.copy(bond.type)
                    break

    
    exists = False
    # If the bond is new, add it to the type list
    if not exists:
        parm.bond_types.append(new_bnd_typ)
        new_bnd_typ.list = parm.bond_types

    #atnum1, atnum2 = -1, -1
    # Loop through all of the selected atoms
    for atnum1,atnum2 in zip(sel1,sel2):
        # Collect the atoms involved
        #atnum1 = sel1.index(1, atnum1+1)
        #atnum2 = sel2.index(1, atnum2+1)
        atm1 = parm.atoms[atnum1]
        atm2 = parm.atoms[atnum2]

        # See if the bond exists in the first place, and if so, replace its
        # bond type with our new bond type (new_bnd)
        if atm2 in atm1.bond_partners and atm1 in atm2.bond_partners:
            for bond in atm1.bonds:
                if atm2 in bond:
                    bond.type = new_bnd_typ
                    parm.bonds.changed = True
                    break
            # Otherwise, it doesn't exist, so we just create a new one
        else:
            parm.bonds.append(Bond(atm1, atm2, new_bnd_typ))
    # Make sure we update 1-4 exception handling if we created any rings
    parm.update_dihedral_exclusions()

    


def AddNewAngleType(parm,sel1,sel2,sel3):
    from parmed.topologyobjects import BondType,AngleType,DihedralType
    from parmed.topologyobjects import Bond,Angle,Dihedral
    import copy
    
    #sel1 = parmed.amber.mask.AmberMask(parm,mask1).Selection()
    #sel2 = parmed.amber.mask.AmberMask(parm,mask2).Selection()
    #sel3 = parmed.amber.mask.AmberMask(parm,mask3).Selection()

    if len(sel1) != len(sel2):
        raise Exception('parmutils.AddNewAngleType: Each mask must select the same '
                        'number of atoms!')

    if len(sel1) != len(sel3):
        raise Exception('parmutils.AddNewAngleType: Each mask must select the same '
                        'number of atoms!')

    # If no atoms, nothing to do
    if len(sel1) == 0: return


    new_ang_typ = None
    #atnum1, atnum2, atnum3 = -1, -1, -1
    # Loop through all of the selections
    for ii in range(len(sel1)):
        # Collect the atoms involved
        atnum1 = sel1[ii]
        atnum2 = sel2[ii]
        atnum3 = sel3[ii]
        atm1 = parm.atoms[atnum1]
        atm2 = parm.atoms[atnum2]
        atm3 = parm.atoms[atnum3]

        # See if the angle exists in the first place, and if so, replace its
        # angle type with our new angle type (new_ang)
        found = False
        if atm1 in atm3.angle_partners:
            for ang in atm1.angles:
                if atm2 in ang and atm3 in ang:
                    new_ang_typ = copy.copy(ang.type)
                    break
    
    exists = False
    # If the angle is new, add it to the type list
    if not exists:
        parm.angle_types.append(new_ang_typ)
        new_ang_typ.list = parm.angle_types

    # Loop through all of the selections
    for ii in range(len(sel1)):
        # Collect the atoms involved
        atnum1 = sel1[ii]
        atnum2 = sel2[ii]
        atnum3 = sel3[ii]
        atm1 = parm.atoms[atnum1]
        atm2 = parm.atoms[atnum2]
        atm3 = parm.atoms[atnum3]

        # See if the angle exists in the first place, and if so, replace its
        # angle type with our new angle type (new_ang)
        found = False
        if atm1 in atm3.angle_partners:
            for ang in atm1.angles:
                if atm2 in ang and atm3 in ang:
                    ang.type = new_ang_typ
                    parm.angles.changed = True
                    found = True
                    break
        # If not found, create a new angle
        if not found:
            parm.angles.append(Angle(atm1, atm2, atm3, new_ang_typ))
    # Make sure we update 1-4 exception handling if we created any rings
    parm.update_dihedral_exclusions()
    





# def AddNewDihedType(parm,mask1,mask2,mask3,mask4,improper=False):
#     from parmed.topologyobjects import BondType,AngleType,DihedralType
#     from parmed.topologyobjects import Bond,Angle,Dihedral
#     import copy
    
#     sel1 = parmed.amber.mask.AmberMask(parm,mask1).Selection()
#     sel2 = parmed.amber.mask.AmberMask(parm,mask2).Selection()
#     sel3 = parmed.amber.mask.AmberMask(parm,mask3).Selection()
#     sel4 = parmed.amber.mask.AmberMask(parm,mask4).Selection()

#     if sum(sel1) != sum(sel2):
#         raise Exception('parmutils.AddNewAngleType: Each mask must select the same '
#                         'number of atoms!')

#     if sum(sel1) != sum(sel3):
#         raise Exception('parmutils.AddNewAngleType: Each mask must select the same '
#                         'number of atoms!')

#     if sum(sel1) != sum(sel4):
#         raise Exception('parmutils.AddNewAngleType: Each mask must select the same '
#                         'number of atoms!')

#     # If no atoms, nothing to do
#     if sum(sel1) == 0: return

#     new_dih_typ = None


#     exists = False

#     if not exists:
#         parm.dihedral_types.append(new_dih_typ)
#         new_dih_typ.list = parm.dihedral_types

#     # Loop through all of the atoms
#     atnum1, atnum2, atnum3, atnum4 = -1, -1, -1, -1

#     for it in range(sum(sel1)):
#         # Collect the atoms involved
#         atnum1 = sel1.index(1, atnum1+1)
#         atnum2 = sel2.index(1, atnum2+1)
#         atnum3 = sel3.index(1, atnum3+1)
#         atnum4 = sel4.index(1, atnum4+1)
#         atm1 = parm.atoms[atnum1]
#         atm2 = parm.atoms[atnum2]
#         atm3 = parm.atoms[atnum3]
#         atm4 = parm.atoms[atnum4]
#         if (atm1 is atm2 or atm1 is atm3 or atm1 is atm4 or
#             atm2 is atm3 or atm2 is atm4 or atm3 is atm4):
#             raise Exception('addDihedral: Duplicate atoms found!')
#         # Determine if end-group interactions need to be ignored
#         ignore_end = (atm1 in atm4.bond_partners or
#                       atm1 in atm4.angle_partners or
#                       atm1 in atm4.dihedral_partners)
#         parm.dihedrals.append(
#             Dihedral(atm1, atm2, atm3, atm4, improper=improper,
#                      ignore_end=ignore_end, type=new_dih_typ)
#         )






    

def FindParamTypeObj(param,atom_type_tuple):
    pt=[]
    ats=atom_type_tuple
    if len(ats) == 2:
        for x in param.bonds:
            atoms = ( x.atom1.type, x.atom2.type )
            if atoms == ats or reversed(atoms) == ats:
                ptappend( x.type )
    elif len(ats) == 3:
        for x in param.angles:
            atoms = ( x.atom1.type, x.atom2.type, x.atom3.type )
            if atoms == ats or reversed(atoms) == ats:
                pt.append( x.type )
    elif len(ats) == 4:
        for x in param.dihedrals:
            atoms = ( x.atom1.type, x.atom2.type, x.atom3.type, x.atom4.type )
            if atoms == ats or reversed(atoms) == ats:
                pt.append( x.type )
    return pt

def WriteFrcmod(param,native_frcmod,uniqueparams=False,with_mass=True,with_nonb=True):
    self = parmed.amber.parameters.AmberParameterSet.from_structure(param)
    WriteFrcmodObj(self,native_frcmod,angfact=0.9999995714245039,uniqueparams=uniqueparams,with_mass=with_mass,with_nonb=with_nonb)

   
def WriteMaskedFrcmod(param,aidxs,native_frcmod,changes_frcmod,with_mass=True,with_nonb=True):
    """
Converts selected atoms to new atom-types and
writes two frcmod files, one containing the
those parameters involving non-substituted
atoms, and the other for those involving the
substituted atoms.

@param param: Amber parm file object
@param aidxs: integer list of atoms that will use new atom-type parameters
@param native_frcmod: string filename of the nonsubstituted parameters
@param changes_frcmod: string filename listing the new atom-type parameters
@return: none
    """
    
    from copy import copy,deepcopy
    from parmed.utils.six import add_metaclass, string_types, iteritems
    from parmed.topologyobjects import BondType,AngleType,DihedralType
    from collections import defaultdict as ddict
    import re
    
    selected_names = {}
    for aidx in aidxs:
        atom = param.atoms[aidx]
        atom.residue.name = atom.residue.name.lower()
        oldname = atom.atom_type.name
        newname = oldname.lower()
        param.atoms[aidx].atom_type.name = newname
        param.atoms[aidx].type = newname
        selected_names[newname] = 1


    self = parmed.amber.parameters.AmberParameterSet.from_structure(param)
    WriteFrcmodObj(self,native_frcmod,angfact=0.9999995714245039,uniqueparams=False,selected_names=selected_names,changed_frcmod=changes_frcmod,with_mass=with_mass,with_nonb=with_nonb)
    return


    # # =============== BONDS ===================================
    # newtypes = {}
    # newtypemask = {}
    # # Get the list of type-parameter indices that pass through the selection
    # for bond in param.bonds:
    #     if bond.atom1.idx in aidxs or bond.atom2.idx in aidxs:
    #         newtypes[bond.type.idx] = bond.type

    # # Create new type-parameters in the param with the same values
    # for idx in newtypes:
    #     bondtype = newtypes[idx]
    #     param.bond_types.append( BondType(bondtype.k,bondtype.req,param.bond_types) )
    #     newtypemask[idx] = param.bond_types[-1].idx

    # # Point to the new parameters, where appropriate
    # for bond in param.bonds:
    #     if bond.atom1.idx in aidxs or bond.atom2.idx in aidxs:
    #         bond.type = param.bond_types[ newtypemask[bond.type.idx] ]

    # # =============== ANGLES ===================================
    # newtypes = {}
    # newtypemask = {}
    # # Get the list of type-parameter indices that pass through the selection
    # for angle in param.angles:
    #     if angle.atom1.idx in aidxs or angle.atom2.idx in aidxs:
    #         newtypes[angle.type.idx] = angle.type

    # # Create new type-parameters in the param with the same values
    # for idx in newtypes:
    #     angletype = newtypes[idx]
    #     param.angle_types.append( AngleType(angletype.k,angletype.theteq,param.angle_types) )
    #     newtypemask[idx] = param.angle_types[-1].idx

    # # Point to the new parameters, where appropriate
    # for angle in param.angles:
    #     if angle.atom1.idx in aidxs or angle.atom2.idx in aidxs:
    #         angle.type = param.angle_types[ newtypemask[angle.type.idx] ]

    # # =============== TORSIONS ===================================

            
    # parmed.tools.writeFrcmod(param, "foo.frcmod").execute()

        
    # if False:

    #     angfact = 0.9999995714245039 

    #     class combofile(object):
    #         def __init__(self,fh1,fh2):
    #             self.fh1 = fh1
    #             self.fh2 = fh2
    #         def write(self,s):
    #             self.fh1.write(s)
    #             self.fh2.write(s)
        
    #     nfile = file(native_frcmod,"w")
    #     if len(aidxs) > 0:
    #         cfile = file(changes_frcmod,"w")
    #         outfile = combofile(cfile,nfile)
    #     else:
    #         cfile = nfile
    #         outfile = nfile

    #     self = parmed.amber.parameters.AmberParameterSet.from_structure(param)
    #     outfile.write("modified parameters")
    #     outfile.write('\n')
    #     # Write the atom mass
    #     outfile.write('MASS\n')
    #     for atom, typ in iteritems(self.atom_types):
    #         fh=nfile
    #         if atom in selected_names:
    #             fh = cfile
    #         fh.write('%s%11.8f\n' % (atom.ljust(6), typ.mass))
                
    #     outfile.write('\n')
    #     # Write the bonds
    #     outfile.write('BOND\n')
    #     cdone = set()
    #     ndone = set()
    #     deltas = ddict( lambda: ddict(float) )
    #     for (a1, a2), typ in iteritems(self.bond_types):
    #         typ.k = float("%.8f"%(typ.k))
    #         fh=nfile
    #         delta = 0
    #         if a1 in selected_names or a2 in selected_names:
    #             fh=cfile
    #             qq = (a1,a2)
    #             if qq in cdone: continue
    #             qq = (a2,a1)
    #             if qq in cdone: continue
    #             cdone.add(qq)
    #             deltas[typ.k][typ.req] += 1.e-13
    #             delta = deltas[typ.k][typ.req]
    #         else:
    #             fh=nfile
    #             if id(typ) in ndone: continue
    #             ndone.add(id(typ))
    #         fh.write('%s-%s   %19.14f  %11.8f\n' %
    #                  (a1.ljust(2), a2.ljust(2), typ.k+delta, typ.req))
    #     outfile.write('\n')
    #     # Write the angles
    #     outfile.write('ANGLE\n')
    #     cdone = set()
    #     ndone = set()
    #     deltas = ddict( lambda: ddict(float) )
    #     for (a1, a2, a3), typ in iteritems(self.angle_types):
    #         typ.k = float("%.8f"%(typ.k))
    #         delta = 0.
    #         if a1 in selected_names or a2 in selected_names or \
    #            a3 in selected_names:
    #             fh=cfile
    #             qq = (a1,a2,a3)
    #             if qq in cdone: continue
    #             qq = (a3,a2,a1)
    #             if qq in cdone: continue
    #             cdone.add(qq)
    #             deltas[typ.k][typ.theteq] += 1.e-13
    #             delta = deltas[typ.k][typ.theteq]
    #         else:
    #             fh=nfile
    #             if id(typ) in ndone: continue
    #             ndone.add(id(typ))
    #         fh.write('%s-%s-%s   %19.14f  %17.3f\n' %
    #                  (a1.ljust(2), a2.ljust(2), a3.ljust(2), typ.k+delta,
    #                   typ.theteq * angfact))
    #     outfile.write('\n')
    #     # Write the dihedrals
    #     outfile.write('DIHE\n')
    #     cdone = set()
    #     ndone = set()
    #     deltas = ddict( lambda: ddict( lambda: ddict( float ) ) )
    #     for (a1, a2, a3, a4), typ in iteritems(self.dihedral_types):
    #         isnew = False
    #         if a1 in selected_names or a2 in selected_names or \
    #            a3 in selected_names or a4 in selected_names:
    #             fh=cfile
    #             qq = (a1,a2,a3,a4)
    #             if qq in cdone: continue
    #             qq = (a4,a3,a2,a1)
    #             if qq in cdone: continue
    #             cdone.add(qq)
    #             isnew = True
    #         else:
    #             fh=nfile
    #             if id(typ) in ndone: continue
    #             ndone.add(id(typ))
    #         if isinstance(typ, DihedralType) or len(typ) == 1:
    #             if not isinstance(typ, DihedralType):
    #                 typ = typ[0]
    #                 typ.phi_k = float("%.8f"%(typ.phi_k))
    #                 delta = 0
    #                 if isnew:
    #                     deltas[typ.phi_k][typ.phase][typ.per] += 1.e-13
    #                     delta = deltas[typ.phi_k][typ.phase][typ.per]
    #             if abs(typ.phase-180) < 0.0001:
    #                 fh.write('%s-%s-%s-%s %4i %20.14f %13.3f %5.1f    '
    #                          'SCEE=%s SCNB=%s\n' % (a1.ljust(2), a2.ljust(2),
    #                                                 a3.ljust(2), a4.ljust(2), 1, typ.phi_k+delta, typ.phase * angfact,
    #                                                 typ.per, typ.scee, typ.scnb))
    #             else:
    #                 fh.write('%s-%s-%s-%s %4i %20.14f %13.8f %5.1f    '
    #                          'SCEE=%s SCNB=%s\n' % (a1.ljust(2), a2.ljust(2),
    #                                                 a3.ljust(2), a4.ljust(2), 1, typ.phi_k+delta, typ.phase * angfact,
    #                                                 typ.per, typ.scee, typ.scnb))
    #         else:
    #             typ = sorted( typ, key=lambda x: x.per, reverse=False )
    #             for dtyp in typ[:-1]:
    #                 dtyp.phi_k = float("%.8f"%(dtyp.phi_k))
    #                 delta = 0
    #                 if isnew:
    #                     deltas[dtyp.phi_k][dtyp.phase][dtyp.per] += 1.e-13
    #                     delta = deltas[dtyp.phi_k][dtyp.phase][dtyp.per]
    #                 if abs(dtyp.phase-180) < 0.0001:
    #                     #print "%20.16f"%(180.0/dtyp.phase)
    #                     fh.write('%s-%s-%s-%s %4i %20.14f %13.3f %5.1f    '
    #                              'SCEE=%s SCNB=%s\n'%(a1.ljust(2), a2.ljust(2),
    #                                                   a3.ljust(2), a4.ljust(2), 1, dtyp.phi_k+delta,
    #                                                   dtyp.phase * angfact, -dtyp.per, dtyp.scee, dtyp.scnb))
    #                 else:
    #                     fh.write('%s-%s-%s-%s %4i %20.14f %13.8f %5.1f    '
    #                              'SCEE=%s SCNB=%s\n'%(a1.ljust(2), a2.ljust(2),
    #                                                   a3.ljust(2), a4.ljust(2), 1, dtyp.phi_k+delta,
    #                                                   dtyp.phase * angfact, -dtyp.per, dtyp.scee, dtyp.scnb))
    #             dtyp = typ[-1]
    #             dtyp.phi_k = float("%.8f"%(dtyp.phi_k))
    #             delta = 0
    #             if isnew:
    #                 deltas[dtyp.phi_k][dtyp.phase][dtyp.per] += 1.e-13
    #                 delta = deltas[dtyp.phi_k][dtyp.phase][dtyp.per]
    #             if abs(dtyp.phase-180) < 0.0001:
    #                 fh.write('%s-%s-%s-%s %4i %20.14f %13.3f %5.1f    '
    #                          'SCEE=%s SCNB=%s\n' % (a1.ljust(2), a2.ljust(2),
    #                                                 a3.ljust(2), a4.ljust(2), 1, dtyp.phi_k+delta,
    #                                                 dtyp.phase * angfact, dtyp.per, dtyp.scee, dtyp.scnb))
    #             else:
    #                 fh.write('%s-%s-%s-%s %4i %20.14f %13.8f %5.1f    '
    #                          'SCEE=%s SCNB=%s\n' % (a1.ljust(2), a2.ljust(2),
    #                                                 a3.ljust(2), a4.ljust(2), 1, dtyp.phi_k+delta,
    #                                                 dtyp.phase * angfact, dtyp.per, dtyp.scee, dtyp.scnb))
                    
    #     outfile.write('\n')
    #     # Write the impropers
    #     deltas = ddict( lambda: ddict( lambda: ddict( float ) ) )
    #     outfile.write('IMPROPER\n')
    #     for (a1, a2, a3, a4), typ in iteritems(self.improper_periodic_types):
    #         # Make sure wild-cards come at the beginning
    #         if a2 == 'X':
    #             assert a4 == 'X', 'Malformed generic improper!'
    #             a1, a2, a3, a4 = a2, a4, a3, a1
    #         elif a4 == 'X':
    #             a1, a2, a3, a4 = a4, a1, a3, a2

    #         typ.phi_k = float("%.8f"%(typ.phi_k))
    #         delta = 0
    #         if a1 in selected_names or a2 in selected_names or \
    #            a3 in selected_names or a4 in selected_names:
    #             fh=cfile
    #             deltas[typ.phi_k][typ.phase][typ.per] += 1.e-13
    #             delta = deltas[typ.phi_k][typ.phase][typ.per]
    #         else:
    #             fh=nfile
    #         if abs(typ.phase-180) < 0.0001:
    #             fh.write('%s-%s-%s-%s %20.14f %13.3f %5.1f\n' %
    #                      (a1.ljust(2), a2.ljust(2), a3.ljust(2), a4.ljust(2),
    #                       typ.phi_k+delta, typ.phase * angfact, typ.per))
    #         else:
    #             fh.write('%s-%s-%s-%s %20.14f %13.8f %5.1f\n' %
    #                      (a1.ljust(2), a2.ljust(2), a3.ljust(2), a4.ljust(2),
    #                       typ.phi_k+delta, typ.phase * angfact, typ.per))

                
    #     outfile.write('\n')
    #     # Write the LJ terms

    #     deltas = ddict( lambda: ddict( float ) )

    #     outfile.write('NONB\n')
    #     for atom, typ in iteritems(self.atom_types):
    #         #typ.rmin = float("%.8f"%(typ.rmin))
    #         typ.epsilon = float("%.9f"%(typ.epsilon))
    #         delta = 0.
    #         if atom in selected_names:
    #             fh=cfile
    #             deltas[typ.rmin][typ.epsilon] += 1.e-13
    #             delta = deltas[typ.rmin][typ.epsilon]
    #         else:
    #             fh=nfile
    #         if delta == 0.:
    #             fh.write('%-3s  %12.8f %18.9f\n' %
    #                      (atom.ljust(2), typ.rmin, typ.epsilon))
    #         else:
    #             fh.write('%-3s  %12.8f %18.14f\n' %
    #                      (atom.ljust(2), typ.rmin, typ.epsilon+delta))
    #     outfile.write('\n')
    #     # Write the NBFIX terms
    #     if self.nbfix_types:
    #         outfile.write('LJEDIT\n')
    #         for (a1, a2), (eps, rmin) in iteritems(self.nbfix_types):
    #             if a1 in selected_names or a2 in selected_names:
    #                 fh=cfile
    #             else:
    #                 fh=nfile
    #             fh.write('%s %s %13.8f %13.8f %13.8f %13.8f\n' %
    #                      (a1.ljust(2), a2.ljust(2), eps, rmin/2,
    #                       eps, rmin/2))
    #     cfile.close()
    #     nfile.close()
        

def statisticalInefficiency(A_n):
    """Compute the (cross) statistical inefficiency of (two) timeseries.

    Parameters
    ----------
    A_n : np.ndarray, float
        A_n[n] is nth value of timeseries A.  Length is deduced from vector.
    B_n : np.ndarray, float, optional, default=None
        B_n[n] is nth value of timeseries B.  Length is deduced from vector.
        If supplied, the cross-correlation of timeseries A and B will be estimated instead of the
        autocorrelation of timeseries A.  

    Returns
    -------
    g : np.ndarray,
        g is the estimated statistical inefficiency (equal to 1 + 2 tau, where tau is the correlation time).
        We enforce g >= 1.0.

    Notes
    -----
    The same timeseries can be used for both A_n and B_n to get the autocorrelation statistical inefficiency.
    The fast method described in Ref [1] is used to compute g.

    References
    ----------
    [1] J. D. Chodera, W. C. Swope, J. W. Pitera, C. Seok, and K. A. Dill. Use of the weighted
        histogram analysis method for the analysis of simulated and parallel tempering simulations.
        JCTC 3(1):26-41, 2007.


    """
    import numpy as np
    # Create numpy copies of input arguments.
    A_n = np.array(A_n)
    B_n = np.array(A_n)

    # Get the length of the timeseries.
    N = A_n.size

    # Initialize statistical inefficiency estimate with uncorrelated value.
    g = 1.0

    # Compute mean of each timeseries.
    mu_A = A_n.mean()
    mu_B = mu_A

    # Make temporary copies of fluctuation from mean.
    dA_n = A_n.astype(np.float64) - mu_A
    dB_n = dA_n

    # Compute estimator of covariance of (A,B) using estimator that will ensure C(0) = 1.
    sigma2_AB = (dA_n * dB_n).mean()  # standard estimator to ensure C(0) = 1

    # Trap the case where this covariance is zero, and we cannot proceed.
    if(sigma2_AB == 0):
        raise ParameterError('Sample covariance sigma_AB^2 = 0 -- cannot compute statistical inefficiency')

    # Accumulate the integrated correlation time by computing the normalized correlation time at
    # increasing values of t.  Stop accumulating if the correlation function goes negative, since
    # this is unlikely to occur unless the correlation function has decayed to the point where it
    # is dominated by noise and indistinguishable from zero.
    t = 1
    while (t < N - 1):

        # compute normalized fluctuation correlation function at time t
        C = np.sum(dA_n[0:(N - t)] * dB_n[t:N] + dB_n[0:(N - t)] * dA_n[t:N]) / (2.0 * float(N - t) * sigma2_AB)
        # Terminate if the correlation function has crossed zero and we've computed the correlation
        # function at least out to 'mintime'.
        if (C <= 0.0) and (t > 3):
            break

        # Accumulate contribution to the statistical inefficiency.
        g += 2.0 * C * (1.0 - float(t) / float(N)) * float(1)

        # Increment t and the amount by which we increment t.
        t += 1

    # g must be at least unity
    if (g < 1.0):
        g = 1.0

    # Return the computed statistical inefficiency.
    return g

        
def UpdateFrcmod(parmfile,ref_traj,trial_traj,trial_frcmod,next_frcmod):
    from copy import copy,deepcopy
    from parmed.utils.six import add_metaclass, string_types, iteritems
    from parmed.topologyobjects import BondType,AngleType,DihedralType
    from collections import defaultdict as ddict
    import pytraj as pt
    import numpy as np
    
    fmod = parmed.load_file(trial_frcmod)
    param = parmed.load_file(parmfile)
    modbonds = ddict( list )
    for bond in param.bonds:
        at1 = bond.atom1.atom_type.name
        at2 = bond.atom2.atom_type.name
        bpair = (at1,at2)
        for (a1,a2),btyp in iteritems(fmod.bond_types):
            if a2 > a1:
                continue
            if (a1,a2) == bpair or (a2,a1) == bpair:
                modbonds[ (a1,a2) ].append( bond )
            
    modangles = ddict( list )
    for angle in param.angles:
        at1 = angle.atom1.atom_type.name
        at2 = angle.atom2.atom_type.name
        at3 = angle.atom3.atom_type.name
        bpair = (at1,at2,at3)
        for (a1,a2,a3),btyp in iteritems(fmod.angle_types):
            if a3 > a1:
                continue
            if (a1,a2,a3) == bpair or (a3,a2,a1) == bpair:
                modangles[ (a1,a2,a3) ].append( angle )

    traj = pt.load( ref_traj, parmfile )

    refbond_cor = []
    refbondstats = ddict( list )
    for key in modbonds:
        v = []
        for bond in modbonds[key]:
            i = bond.atom1.idx+1
            j = bond.atom2.idx+1
            vals = pt.distance( traj, "@%i @%i"%(i,j) )
            refbond_cor.append( statisticalInefficiency(vals) )
            v.extend( vals )
        a = np.average(v)
        s = np.std(v)
        refbondstats[key] = (a,s)

    refangle_cor = []
    refanglestats = ddict( list )
    for key in modangles:
        v = []
        for angle in modangles[key]:
            i = angle.atom1.idx+1
            j = angle.atom2.idx+1
            k = angle.atom3.idx+1
            vals = pt.angle( traj, "@%i @%i @%i"%(i,j,k) )
            refangle_cor.append( statisticalInefficiency(vals) )
            v.extend( vals )
        a = np.average(v)
        s = np.std(v)
        refanglestats[key] = (a,s)

    a = np.average(refbond_cor)
    m = np.max(refbond_cor)
    print("Ref bond  stat. ineff: %6.1f (max: %6.1f)"%(a,m))
    a = np.average(refangle_cor)
    m = np.max(refangle_cor)
    print("Ref angle stat. ineff: %6.1f (max: %6.1f)"%(a,m))


    traj = pt.load( trial_traj, parmfile )
    
    bondstats = ddict( list )
    for key in modbonds:
        v = []
        for bond in modbonds[key]:
            i = bond.atom1.idx+1
            j = bond.atom2.idx+1
            v.extend( pt.distance( traj, "@%i @%i"%(i,j) ) )
        a = np.average(v)
        s = np.std(v)
        bondstats[key] = (a,s)
        if False:
            fname = "bond.%s-%s.dat"%(key[0],key[1])
            fh=open(fname,"w")
            n=len(v)
            for i in range(n):
                a = np.average(v[:i+1])
                s = np.std(v[:i+1])
                fh.write("%i6 %13.4e %13.4e\n"%(i+1,a,s))
            fh.close()

    anglestats = ddict( list )
    for key in modangles:
        v = []
        for angle in modangles[key]:
            i = angle.atom1.idx+1
            j = angle.atom2.idx+1
            k = angle.atom3.idx+1
            v.extend( pt.angle( traj, "@%i @%i @%i"%(i,j,k) ) )
        a = np.average(v)
        s = np.std(v)
        anglestats[key] = (a,s)

    seen = set()
    for (a1,a2),par in iteritems(fmod.bond_types):
        if a2 > a1:
            continue
        key=(a1,a2)
#        if (a1,a2) in seen or (a2,a1) in seen:
#            continue
#        else:
#            seen.add( key )

        kold = par.k
        eqold = par.req

        refavg,refstd = refbondstats[key]
        obsavg,obsstd = bondstats[key]

        knew = kold
        if obsstd > 0.002 and refstd > 0.002:
            knew = kold * ( 0.2 + 0.8 * math.pow( obsstd / refstd , 2 ) )
            knew = min(knew,800.)
        eqnew=eqold
        if obsavg > 0.002 and refavg > 0.002:
            eqnew = eqold * ( refavg / obsavg )

        par.k = knew
        par.req = eqnew
        
        print("%2s %2s     k: %8.3f <- %8.3f  eq: %8.4f <- %8.4f  avg: %8.3f/%8.3f  std: %8.3f/%8.3f"%(
            a1,a2,
            knew,kold,eqnew,eqold,
            refavg,obsavg,refstd,obsstd ))



    #seen=set()
    for (a1,a2,a3),par in iteritems(fmod.angle_types):
        if a3 > a1:
            continue
        key=(a1,a2,a3)
#        if (a1,a2,a3) in seen or (a3,a2,a1) in seen:
#            continue
#        else:
#            seen.add(key)
            
        kold = par.k
        eqold = par.theteq

        refavg,refstd = refanglestats[key]
        obsavg,obsstd = anglestats[key]

        knew = kold
        if obsstd > 0.002 and refstd > 0.002:
            knew = kold * ( 0.2 + 0.8 * math.pow( obsstd / refstd , 2 ) )
            knew = min(knew,800.)
        eqnew=eqold
        if obsavg > 0.002 and refavg > 0.002:
            eqnew = eqold * ( refavg / obsavg )
        par.k = knew
        par.theteq = eqnew
        
        print("%2s %2s %2s  k: %8.3f <- %8.3f  eq: %8.4f <- %8.4f  avg: %8.3f/%8.3f  std: %8.3f/%8.3f"%(
            a1,a2,a3,
            knew,kold,eqnew,eqold,
            refavg,obsavg,refstd,obsstd ))

    #SymmetrizeFrcmodObj(fmod)
    WriteFrcmodObj(fmod,next_frcmod,uniqueparams=True)


# def SymmetrizeFrcmodObj(fmod):
#     from parmed.utils.six import add_metaclass, string_types, iteritems
#     for (a1,a2),par in iteritems(fmod.bond_types):
#         if a2 > a1:
#             fmod.bond_types[ (a1,a2) ] = fmod.bond_types[ (a2,a1) ]
#     for (a1,a2,a3),par in iteritems(fmod.angle_types):
#         if a3 > a1:
#             fmod.angle_types[ (a1,a2,a3) ] = fmod.angle_types[ (a3,a2,a1) ]
#     for (a1,a2,a3,a4),par in iteritems(fmod.dihedral_types):
#         if a4 > a1 or (a4==a1 and a3>a2):
#             fmod.dihedral_types[ (a1,a2,a3,a4) ] = fmod.dihedral_types[ (a4,a3,a2,a1) ]

    
        
def MakeAtomTypesCommon( param0, param1, aidxs ):
    from copy import copy,deepcopy
    from parmed.utils.six import add_metaclass, string_types, iteritems
    from parmed.topologyobjects import BondType,AngleType,DihedralType,Dihedral
    from collections import defaultdict as ddict

    newtypes = []
    inew = 0
    for iat in range(len(param0.atoms)):
        a = param0.atoms[iat]
        b = param1.atoms[iat]
        if iat in aidxs:
            if a.atom_type.name != b.atom_type.name:
                t = "X%1i"%(inew)
                newtypes.append(t)
                inew += 1
                a.atom_type.name = t
                b.atom_type.name = t
                a.type = t
                b.type = t
            else:
                newtypes.append(a.type)
        else:
            if a.atom_type.name != b.atom_type.name:
                raise Exception("Atom %i has different atom_type's in the 2 parm7 files"%(iat+1))
            
    self0 = parmed.amber.parameters.AmberParameterSet.from_structure(param0)
    dtypes0 = ddict( list )
    cdone = set()
    for (a1, a2, a3, a4), typ0 in iteritems(self0.dihedral_types):
        if a1 in newtypes or a2 in newtypes or a3 in newtypes or a4 in newtypes:
            qq = (a1,a2,a3,a4)
            if qq in cdone: continue
            qq = (a4,a3,a2,a1)
            if qq in cdone: continue
            cdone.add(qq)
            #print a1,a2,a3,a4
            dtypes0[qq] = typ0
            #for t in typ0:
            #    t.phi_k = 10.
                #help(t)
                #exit(0)

    self1 = parmed.amber.parameters.AmberParameterSet.from_structure(param1)
    dtypes1 =ddict( list )
    cdone = set()
    for (a1, a2, a3, a4), typ1 in iteritems(self1.dihedral_types):
        if a1 in newtypes or a2 in newtypes or a3 in newtypes or a4 in newtypes:
            qq = (a1,a2,a3,a4)
            if qq in cdone: continue
            qq = (a4,a3,a2,a1)
            if qq in cdone: continue
            if qq not in dtypes0:
                qq = (a1,a2,a3,a4)
            cdone.add(qq)
            #print a1,a2,a3,a4
            dtypes1[qq] = typ1
            #print qq
            #for t in typ1:
            #    t.phi_k = 10.
            #help(typ1)
            #exit(0)

    missing0=ddict( list )
    for qq in dtypes1:
        for t in dtypes1[qq]:
            missing=True
            for s in dtypes0[qq]:
                if s.per == t.per:
                    missing = False
            if missing:
                missing0[qq].append( DihedralType( 0., t.per, t.phase, t.scee, t.scnb ) )
                print("adding dihe to lam0 from lam1: %s (k %7.3f phase %7.3f per %6.2f"%(str(qq),t.phi_k,t.phase,t.per))

    missing1=ddict( list )
    for qq in dtypes0:
        for t in dtypes0[qq]:
            missing=True
            for s in dtypes1[qq]:
                if s.per == t.per:
                    missing = False
            if missing:
                missing1[qq].append( DihedralType( 0., t.per, t.phase, t.scee, t.scnb ) )
                print("adding dihe to lam1 from lam0: %s (k %7.3f phase %7.3f per %6.2f"%(str(qq),t.phi_k,t.phase,t.per))

    missing0map = ddict( list )
    for qq in missing0:
        for t in missing0[qq]:
            param0.dihedral_types.append( \
                DihedralType(0., t.per, t.phase, t.scee, t.scnb, param0.dihedral_types ) )
            missing0map[qq].append( param0.dihedral_types[-1].idx )

    missing1map = ddict( list )
    for qq in missing1:
        for t in missing1[qq]:
            param1.dihedral_types.append( \
                DihedralType(0., t.per, t.phase, t.scee, t.scnb, param1.dihedral_types ) )
            missing1map[qq].append( param1.dihedral_types[-1].idx )
        
                
           

        
    # # =============== TORSIONS ===================================
    for it in range(2):
        if it == 0:
            p0 = param0
            mmap = missing0map
            p1 = param1
        else:
            p0 = param1
            mmap = missing1map
            p1 = param0
        seen_quartets = set()
        for dihedral in p0.dihedrals:
            a = dihedral.atom1
            b = dihedral.atom2
            c = dihedral.atom3
            d = dihedral.atom4
            i = a.idx
            j = b.idx
            k = c.idx
            l = d.idx

            if i in aidxs or j in aidxs or k in aidxs or l in aidxs:
                quartet = (i,j,k,l)
                if quartet in seen_quartets:
                    continue
                else:
                    seen_quartets.add( quartet )
                qq = (a.type,b.type,c.type,d.type)
                if qq not in mmap:
                    qq = (d.type,c.type,b.type,a.type)
                    if qq not in mmap:
                        continue
                for t in mmap[qq]:
                    p0.dihedrals.append( \
                        Dihedral( a,b,c,d,dihedral.improper,dihedral.ignore_end, \
                                  p0.dihedral_types[ t ] ) )
                
    # # =============== TORSIONS ===================================



class Fragment(object):
    def __init__(self,parmobj,ambmask,coef0=None,coef1=None,method="AM1D"):
        self.ambmask = ambmask
        self.parmobj = parmobj
        self.parmfilename = None
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
        origname = self.parmobj.residues[ires].name
        newname = self.get_funkified_residue_name(origname)
        if newname is not None:
            self.parmobj.residues[ires].name = newname

    def get_funkified_residue_name(self,origname):
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
        return self.coef0 + lam*(self.coef1-self.coef0)
 
    def get_touched_residues(self):
        residues = []
        for a in self.atomsel:
            res = self.parmobj.atoms[a].residue.idx
            residues.append( res )
        return list(set(residues))

    def get_selected_mmcharge_from_each_touched_residue(self):
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
        cats=[]
        for bond in self.parmobj.bonds:
            if bond.atom1.idx == idx:
                cats.append(bond.atom2.idx)
            elif bond.atom2.idx == idx:
                cats.append(bond.atom1.idx)
        return cats
        
class FragmentedSys(object):
    def __init__(self,parmobj,compobj):
        from . import mdinutils
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
        self.mm_parmfilename = fname
        
    def add_fragment(self,ambmask,coef0,coef1=None,method="AM1D"):
        self.frags.append( Fragment(self.parmobj,ambmask,coef0=coef0,coef1=coef1,method=method) )

    def GetQMAtoms(self):
        qmatoms=[]
        for f in self.frags:
            qmatoms.extend( f.atomsel )
        qmatoms = list(set(qmatoms))
        return qmatoms

    def GetMMTermsForQMRegion(self):
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
    
        
    def GetMMBoundaryTerms(self):
        bonds=[]
        angles=[]
        dihes=[]
        for frag in self.frags:
            b,a,d = frag.GetMMBoundaryTerms()
            bonds.extend(b)
            angles.extend(a)
            dihes.extend(d)
        if len(bonds) > 0:
            bonds  = list(set(bonds))
        if len(angles) > 0:
            angles = list(set(angles))
        if len(dihes) > 0:
            dihes  = list(set(dihes))
        return bonds,angles,dihes
        
    def GetMMBoundaryTermsAsDict(self):
        p = ddict( lambda:ddict( list ) )
        b,a,d = self.GetMMBoundaryTerms()
        for x in b:
            p["bonds"][x.type.idx].append(x)
        for x in a:
            p["angles"][x.type.idx].append(x)
        for x in d:
            p["dihedrals"][x.type.idx].append(x)
        return p
class FragmentedSys(object):
    def __init__(self,parmobj,compobj):
        from . import mdinutils
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
        self.ig = random.randint(1,650
    
    def MakeNewMMBoundaryTerms(self):
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
        self.frags = [ f for f in self.frags if f.method != "MM" ]
        s0 = 0.
        s1 = 0.
        for f in self.frags:
            s0 += f.get_coef(0)
            s1 += f.get_coef(1)
        self.frags.append( Fragment(self.parmobj,":0",coef0=(1-s0),coef1=(1-s1),method="MM") )

    def sort(self):
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
        alist = []
        for f in self.frags:
            alist.extend( f.atomsel )
        return ListToSelection(alist)
        
    def check_overlaps(self):
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
        c=[]
        for f in self.frags:
            c.append( f.get_coef(lam) )
        return c

    def get_energy(self,fragenes,lam=0):
        cs = self.get_coefs(lam)
        e = 0.
        for c,efrag in zip(cs,fragenes):
            e += c * efrag
        return e

    def get_dvdl(self,fragenes):
        return self.get_energy(fragenes,1.) - self.get_energy(fragenes,0.)

    def get_mbar(self,fragenes,nlam=11):
        ene = []
        for i in range(nlam):
            lam = i/(nlam-1.)
            ene.append( self.get_energy(fragenes,lam) )
        return ene
    
    def get_dvdl_coefs(self):
        c0=self.get_coefs(0)
        c1=self.get_coefs(1)
        dc=[]
        for a,b in zip(c0,c1):
            dc.append( b-a )

    def redistribute_residue_charges(self):
        self.sort()
        for f in reversed(self.frags):
            f.redistribute_residue_charges()


    def write_parm(self,parmname="frag.parm7",overwrite=True):
        import subprocess
        self.parmfilename = parmname
        aidxs = []
        for f in self.frags:
            aidxs.extend( f.atomsel )
        if len(aidxs) < 1:
            raise Exception("No fragments")

        base=parmname.replace(".parm7","")
        print("Writing %s.notsele.frcmod and %s.sele.frcmod"%(base,base))
        WriteMaskedFrcmod(self.parmobj,aidxs,"%s.notsele.frcmod"%(base),"%s.sele.frcmod"%(base))

        print("Writing %s.pdb"%(base))
        parmed.tools.writeCoordinates(self.parmobj, "%s.pdb"%(base)).execute()
    
        print("Writing %s.lib"%(base))
        parmed.tools.writeOFF(self.parmobj, "%s.lib"%(base)).execute()

        print("Writing %s.sh"%(base))
        WriteLeapSh \
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
        WriteLeapSh \
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
        has_abinitio = False
        for f in self.frags:
            if f.method == "AM1D" or f.method == "DFTB" or f.method == "MM":
                continue
            else:
                has_abinitio = True
                break
        return has_abinitio

        
    def redistribute_cores(self):
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
                g = statisticalInefficiency( dvdl )
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

    
    
        
            
def WriteFrcmodObj(self,native_frcmod,angfact=1.0,uniqueparams=False,selected_names=None,changed_frcmod=None,with_mass=True,with_nonb=True):
    
    from copy import copy,deepcopy
    from parmed.utils.six import add_metaclass, string_types, iteritems
    from parmed.topologyobjects import BondType,AngleType,DihedralType
    from collections import defaultdict as ddict
    import re

    if uniqueparams and selected_names is None:
        selected_names = {}
        for atom, typ in iteritems(self.atom_types):
            selected_names[atom]=1
    elif selected_names is None:
        selected_names = {}
        

    if True:

#        angfact = 0.9999995714245039

        class combofile(object):
            def __init__(self,fh1,fh2):
                self.fh1 = fh1
                self.fh2 = fh2
            def write(self,s):
                self.fh1.write(s)
                self.fh2.write(s)
        
        nfile = open(native_frcmod,"w")
        if changed_frcmod is None:
            cfile = nfile
            outfile = nfile
        else:
            cfile = open(changed_frcmod,"w")
            outfile = combofile( nfile, cfile )

#        self = parmed.amber.parameters.AmberParameterSet.from_structure(param)
        outfile.write("modified parameters")
        outfile.write('\n')
        # Write the atom mass
        outfile.write('MASS\n')
        if with_mass:
            for atom, typ in iteritems(self.atom_types):
                fh=nfile
                if atom in selected_names:
                    fh = cfile
                fh.write('%s%11.8f\n' % (atom.ljust(6), typ.mass))
                
        outfile.write('\n')
        # Write the bonds
        outfile.write('BOND\n')
        cdone = set()
        ndone = set()
        deltas = ddict( lambda: ddict(float) )
        for (a1, a2), typ in iteritems(self.bond_types):
            typ.k = float("%.8f"%(typ.k))
            fh=nfile
            delta = 0
            if a1 in selected_names or a2 in selected_names:
                fh=cfile
                qq = (a1,a2)
                if qq in cdone: continue
                qq = (a2,a1)
                if qq in cdone: continue
                cdone.add(qq)
                deltas[typ.k][typ.req] += 1.e-13
                delta = deltas[typ.k][typ.req]
            else:
                fh=nfile
                if id(typ) in ndone: continue
                ndone.add(id(typ))
            fh.write('%s-%s   %19.14f  %11.8f\n' %
                     (a1.ljust(2), a2.ljust(2), typ.k+delta, typ.req))
        outfile.write('\n')
        # Write the angles
        outfile.write('ANGLE\n')
        cdone = set()
        ndone = set()
        deltas = ddict( lambda: ddict(float) )
        for (a1, a2, a3), typ in iteritems(self.angle_types):
            typ.k = float("%.8f"%(typ.k))
            delta = 0.
            if a1 in selected_names or a2 in selected_names or \
               a3 in selected_names:
                fh=cfile
                qq = (a1,a2,a3)
                if qq in cdone: continue
                qq = (a3,a2,a1)
                if qq in cdone: continue
                cdone.add(qq)
                deltas[typ.k][typ.theteq] += 1.e-13
                delta = deltas[typ.k][typ.theteq]
            else:
                fh=nfile
                if id(typ) in ndone: continue
                ndone.add(id(typ))
            fh.write('%s-%s-%s   %19.14f  %17.3f\n' %
                     (a1.ljust(2), a2.ljust(2), a3.ljust(2), typ.k+delta,
                      typ.theteq * angfact))
        outfile.write('\n')
        # Write the dihedrals
        outfile.write('DIHE\n')
        cdone = set()
        ndone = set()
        deltas = ddict( lambda: ddict( lambda: ddict( float ) ) )
        for (a1, a2, a3, a4), typ in iteritems(self.dihedral_types):
            isnew = False
            if a1 in selected_names or a2 in selected_names or \
               a3 in selected_names or a4 in selected_names:
                fh=cfile
                qq = (a1,a2,a3,a4)
                if qq in cdone: continue
                qq = (a4,a3,a2,a1)
                if qq in cdone: continue
                cdone.add(qq)
                isnew = True
            else:
                fh=nfile
                if id(typ) in ndone: continue
                ndone.add(id(typ))
            if isinstance(typ, DihedralType) or len(typ) == 1:
                if not isinstance(typ, DihedralType):
                    typ = typ[0]
                    typ.phi_k = float("%.8f"%(typ.phi_k))
                    delta = 0
                    if isnew:
                        deltas[typ.phi_k][typ.phase][typ.per] += 1.e-13
                        delta = deltas[typ.phi_k][typ.phase][typ.per]
                if abs(typ.phase-180) < 0.0001:
                    fh.write('%s-%s-%s-%s %4i %20.14f %13.3f %5.1f    '
                             'SCEE=%s SCNB=%s\n' % (a1.ljust(2), a2.ljust(2),
                                                    a3.ljust(2), a4.ljust(2), 1, typ.phi_k+delta, typ.phase * angfact,
                                                    typ.per, typ.scee, typ.scnb))
                else:
                    fh.write('%s-%s-%s-%s %4i %20.14f %13.8f %5.1f    '
                             'SCEE=%s SCNB=%s\n' % (a1.ljust(2), a2.ljust(2),
                                                    a3.ljust(2), a4.ljust(2), 1, typ.phi_k+delta, typ.phase * angfact,
                                                    typ.per, typ.scee, typ.scnb))
            else:
                typ = sorted( typ, key=lambda x: x.per, reverse=False )
                for dtyp in typ[:-1]:
                    dtyp.phi_k = float("%.8f"%(dtyp.phi_k))
                    delta = 0
                    if isnew:
                        deltas[dtyp.phi_k][dtyp.phase][dtyp.per] += 1.e-13
                        delta = deltas[dtyp.phi_k][dtyp.phase][dtyp.per]
                    if abs(dtyp.phase-180) < 0.0001:
                        #print "%20.16f"%(180.0/dtyp.phase)
                        fh.write('%s-%s-%s-%s %4i %20.14f %13.3f %5.1f    '
                                 'SCEE=%s SCNB=%s\n'%(a1.ljust(2), a2.ljust(2),
                                                      a3.ljust(2), a4.ljust(2), 1, dtyp.phi_k+delta,
                                                      dtyp.phase * angfact, -dtyp.per, dtyp.scee, dtyp.scnb))
                    else:
                        fh.write('%s-%s-%s-%s %4i %20.14f %13.8f %5.1f    '
                                 'SCEE=%s SCNB=%s\n'%(a1.ljust(2), a2.ljust(2),
                                                      a3.ljust(2), a4.ljust(2), 1, dtyp.phi_k+delta,
                                                      dtyp.phase * angfact, -dtyp.per, dtyp.scee, dtyp.scnb))
                dtyp = typ[-1]
                dtyp.phi_k = float("%.8f"%(dtyp.phi_k))
                delta = 0
                if isnew:
                    deltas[dtyp.phi_k][dtyp.phase][dtyp.per] += 1.e-13
                    delta = deltas[dtyp.phi_k][dtyp.phase][dtyp.per]
                if abs(dtyp.phase-180) < 0.0001:
                    fh.write('%s-%s-%s-%s %4i %20.14f %13.3f %5.1f    '
                             'SCEE=%s SCNB=%s\n' % (a1.ljust(2), a2.ljust(2),
                                                    a3.ljust(2), a4.ljust(2), 1, dtyp.phi_k+delta,
                                                    dtyp.phase * angfact, dtyp.per, dtyp.scee, dtyp.scnb))
                else:
                    fh.write('%s-%s-%s-%s %4i %20.14f %13.8f %5.1f    '
                             'SCEE=%s SCNB=%s\n' % (a1.ljust(2), a2.ljust(2),
                                                    a3.ljust(2), a4.ljust(2), 1, dtyp.phi_k+delta,
                                                    dtyp.phase * angfact, dtyp.per, dtyp.scee, dtyp.scnb))
                    
        outfile.write('\n')
        # Write the impropers
        deltas = ddict( lambda: ddict( lambda: ddict( float ) ) )
        outfile.write('IMPROPER\n')
        for (a1, a2, a3, a4), typ in iteritems(self.improper_periodic_types):
            # Make sure wild-cards come at the beginning
            if a2 == 'X':
                assert a4 == 'X', 'Malformed generic improper!'
                a1, a2, a3, a4 = a2, a4, a3, a1
            elif a4 == 'X':
                a1, a2, a3, a4 = a4, a1, a3, a2

            typ.phi_k = float("%.8f"%(typ.phi_k))
            delta = 0
            if a1 in selected_names or a2 in selected_names or \
               a3 in selected_names or a4 in selected_names:
                fh=cfile
                deltas[typ.phi_k][typ.phase][typ.per] += 1.e-13
                delta = deltas[typ.phi_k][typ.phase][typ.per]
            else:
                fh=nfile
            if abs(typ.phase-180) < 0.0001:
                fh.write('%s-%s-%s-%s %20.14f %13.3f %5.1f\n' %
                         (a1.ljust(2), a2.ljust(2), a3.ljust(2), a4.ljust(2),
                          typ.phi_k+delta, typ.phase * angfact, typ.per))
            else:
                fh.write('%s-%s-%s-%s %20.14f %13.8f %5.1f\n' %
                         (a1.ljust(2), a2.ljust(2), a3.ljust(2), a4.ljust(2),
                          typ.phi_k+delta, typ.phase * angfact, typ.per))

                
        outfile.write('\n')
        # Write the LJ terms

        deltas = ddict( lambda: ddict( float ) )

        outfile.write('NONB\n')
        if with_nonb:
            for atom, typ in iteritems(self.atom_types):
                #typ.rmin = float("%.8f"%(typ.rmin))
                typ.epsilon = float("%.9f"%(typ.epsilon))
                delta = 0.
                if atom in selected_names:
                    fh=cfile
                    deltas[typ.rmin][typ.epsilon] += 1.e-13
                    delta = deltas[typ.rmin][typ.epsilon]
                else:
                    fh=nfile
                if delta == 0.:
                    fh.write('%-3s  %12.8f %18.9f\n' %
                             (atom.ljust(2), typ.rmin, typ.epsilon))
                else:
                    fh.write('%-3s  %12.8f %18.14f\n' %
                             (atom.ljust(2), typ.rmin, typ.epsilon+delta))
            outfile.write('\n')
            # Write the NBFIX terms
            if self.nbfix_types:
                outfile.write('LJEDIT\n')
                for (a1, a2), (eps, rmin) in iteritems(self.nbfix_types):
                    if a1 in selected_names or a2 in selected_names:
                        fh=cfile
                    else:
                        fh=nfile
                    fh.write('%s %s %13.8f %13.8f %13.8f %13.8f\n' %
                             (a1.ljust(2), a2.ljust(2), eps, rmin/2,
                              eps, rmin/2))
        cfile.close()
        nfile.close()
        


# def AvgFrcmods(frcmods,output):
#     from copy import copy,deepcopy
#     from parmed.utils.six import add_metaclass, string_types, iteritems
#     from parmed.topologyobjects import BondType,AngleType,DihedralType
#     from collections import defaultdict as ddict
#     import pytraj as pt
#     import numpy as np
    
#     fmods = [ parmed.load_file(fmod) for fmod in frcmods ]

    
#     kk=[]
#     for fmod in fmods:
#         kk.extend( [ (a1,a2) for (a1,a2) in fmod.bond_types if a1 <= a2 ] )
#     kk=set(kk)
#     nmod = parmed.load_file(frcmods[0])
#     for (a1,a2) in kk:
#         fkey=(a1,a2)
#         rkey=(a2,a1)
#         n=0.
#         k=0.
#         eq=0.
#         for fmod in fmods:
#             if fkey in fmod:
#                 k  += fmod[fkey].k
#                 eq += fmod[fkey].req
#                 n  += 1.
#         typ = BondType( k/n, eq/n )
#         nmod.bond_types[fkey] = typ
#         nmod.bond_types[rkey] = typ
        
#     kk=[]
#     for fmod in fmods:
#         kk.extend( [ (a1,a2,a3) for (a1,a2,a3) in fmod.angle_types if a3 <= a2 ] )
#     kk=set(kk)
#     nmod = parmed.load_file(frcmods[0])
#     for (a1,a2) in kk:
#         fkey=(a1,a2,a3)
#         rkey=(a3,a2,a1)
#         n=0.
#         k=0.
#         eq=0.
#         for fmod in fmods:
#             if fkey in fmod:
#                 k  += fmod[fkey].k
#                 eq += fmod[fkey].theteq
#                 n  += 1.
#         typ = AngleType( k/n, eq/n )
#         nmod.angle_types[fkey] = typ
#         nmod.angle_types[rkey] = typ

