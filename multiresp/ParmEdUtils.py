"""ParmEd open/copy/save/mask helpers for multi-RESP fragment work."""

from __future__ import annotations

import parmed

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

def CopyParm(parm):
    """Shallow-copy a ParmEd AmberParm, including coordinates and box."""
    import copy

    try:
        parm.remake_parm()
    except Exception:
        pass
    p = copy.copy(parm)
    p.coordinates = copy.copy(parm.coordinates)
    p.box = copy.copy(parm.box)
    try:
        p.hasbox = copy.copy(parm.hasbox)
    except Exception:
        p.hasbox = False
    return p

def MakeUniqueParams(p, xlist, *, type_attr: str, type_factory, scale: float = 1.0):
    """Duplicate selected bonded parameter types so they are unique to ``xlist``.

    Parameters
    ----------
    p : parmed Structure
        Parameter container (has ``bond_types`` / ``angle_types`` / ...).
    xlist : list
        Bond / angle / dihedral objects whose types should be uniquified.
    type_attr : str
        Name of the type list attribute on ``p`` (e.g. ``"bond_types"``).
    type_factory : callable
        ``type_factory(old_type, scale, type_list) -> new_type``.
    scale : float, optional
        Scale factor applied when cloning types.
    """
    from collections import defaultdict as ddict

    byidx = ddict(list)
    for x in xlist:
        byidx[x.type.idx].append(x)
    type_list = getattr(p, type_attr)
    for idx in byidx:
        old = byidx[idx][0].type
        type_list.append(type_factory(old, scale, type_list))
        for x in byidx[idx]:
            x.type = type_list[-1]


def MakeUniqueBondParams(p, xlist, scale=1.0):
    """Make unique bond parameters for the given bonds."""
    MakeUniqueParams(
        p,
        xlist,
        type_attr="bond_types",
        type_factory=lambda t, s, lst: parmed.BondType(t.k * s, t.req, lst),
        scale=scale,
    )


def MakeUniqueAngleParams(p, xlist, scale=1.0):
    """Make unique angle parameters for the given angles."""
    MakeUniqueParams(
        p,
        xlist,
        type_attr="angle_types",
        type_factory=lambda t, s, lst: parmed.AngleType(t.k * s, t.theteq, lst),
        scale=scale,
    )


def MakeUniqueDihedralParams(p, xlist, scale=1.0):
    """Make unique dihedral parameters for the given dihedrals."""
    MakeUniqueParams(
        p,
        xlist,
        type_attr="dihedral_types",
        type_factory=lambda t, s, lst: parmed.DihedralType(
            t.phi_k * s, t.per, t.phase, t.scee, t.scnb, lst
        ),
        scale=scale,
    )

                
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


def SaveParm(parm, fname, overwrite=True):
    """Write an Amber topology (parm7) from a ParmEd structure.

    EndState and gas-phase RESP paths call this after charge edits.
    """
    parm.save(fname, overwrite=overwrite)
