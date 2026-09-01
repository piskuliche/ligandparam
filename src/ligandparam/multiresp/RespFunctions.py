def ReadNextRespCharges(fh):
    """ Read the next set of RESP charges from the file handle
    
    Parameters
    ----------
    fh : file handle
        The file handle to read the charges from
    
    Returns
    -------
    qs : list of float
        The list of charges
    """
    import re
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
def GetAvgStdDevAndErr(x):
    """Get the average, standard deviation, and standard error of the mean.

    Parameters
    ----------
    x : list of float
        The list of numbers.

    Returns
    -------
    avg : float
        Arithmetic mean.
    stddev : float
        Sample standard deviation.
    stderr : float
        Standard error of the mean.
    """
    import math
    # Could use numpy, but kept pure-Python for minimal deps in this helper.
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

def BuildBondAdjacency(parm):
    """Build a bond adjacency map from a ParmEd structure.

    Parameters
    ----------
    parm : parmed.Structure
        The ParmEd structure.

    Returns
    -------
    dict of int to list of int
        Mapping from atom index to bonded neighbor indices.
    """
    adjacency = {atom.idx: [] for atom in parm.atoms}
    for bond in parm.bonds:
        i = bond.atom1.idx
        j = bond.atom2.idx
        adjacency[i].append(j)
        adjacency[j].append(i)
    return adjacency


def GetBondAdjacency(parm):
    """Return a bond adjacency map, building and caching it on ``parm`` if needed.

    Parameters
    ----------
    parm : parmed.Structure
        The ParmEd structure.

    Returns
    -------
    dict of int to list of int
        Mapping from atom index to bonded neighbor indices.
    """
    adjacency = getattr(parm, "_bond_adjacency", None)
    if adjacency is None:
        adjacency = BuildBondAdjacency(parm)
        try:
            parm._bond_adjacency = adjacency
        except (AttributeError, TypeError):
            pass
    return adjacency


def GetAtomsBondedToIdx(parm, idx, adjacency=None):
    """Get the atoms bonded to a given atom index.

    Parameters
    ----------
    parm : parmed.Structure
        The ParmEd structure.
    idx : int
        The atom index.
    adjacency : dict of int to list of int, optional
        Precomputed bond adjacency. If omitted, one is built/cached on ``parm``.

    Returns
    -------
    list of int
        Indices of atoms bonded to ``idx``.
    """
    if adjacency is None:
        adjacency = GetBondAdjacency(parm)
    return adjacency.get(idx, [])


def GetEquivNeighbors(
    parm,
    sele,
    *,
    center_pred,
    neighbor_pred,
):
    """Collect sorted neighbor groups bonded to selected centers.

    For each atom in ``sele`` that passes ``center_pred``, gather bonded atoms
    that pass ``neighbor_pred``. Groups with more than one neighbor are returned.
    """
    adjacency = GetBondAdjacency(parm)
    groups = []
    for hvy in sele:
        if not center_pred(parm.atoms[hvy]):
            continue
        neighbors = [
            bat
            for bat in GetAtomsBondedToIdx(parm, hvy, adjacency=adjacency)
            if neighbor_pred(parm.atoms[bat])
        ]
        if len(neighbors) > 1:
            neighbors.sort()
            groups.append(neighbors)
    return groups


def GetEquivHydrogens(parm, sele):
    """Get the equivalent hydrogens bonded to a given atom.

    Parameters
    ----------
    parm : parmed structure
        The ParmEd structure.
    sele : list of int
        The atom indices.

    Returns
    -------
    hpairs : list of list of int
        The indices of the equivalent hydrogens.
    """
    return GetEquivNeighbors(
        parm,
        sele,
        center_pred=lambda atom: atom.atomic_number > 1,
        neighbor_pred=lambda atom: atom.atomic_number == 1,
    )


def GetEquivNonbridgeOxygens(parm, sele):
    """Get the equivalent non-bridge oxygens bonded to a given atom.

    Parameters
    ----------
    parm : parmed structure
        The ParmEd structure.
    sele : list of int
        The atom indices.

    Returns
    -------
    hpairs : list of list of int
        The indices of the equivalent non-bridge oxygens.
    """
    return GetEquivNeighbors(
        parm,
        sele,
        center_pred=lambda atom: atom.atomic_number == 15,
        neighbor_pred=lambda atom: atom.type in ("O2", "o2"),
    )

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

def WriteArray8(fh,arr):
    """ Write an array of numbers to a file handle in 8 columns
    
    Parameters
    ----------
    fh : file handle
        The file handle to write to
    arr : list of float
        The array of numbers
    """

    for istart in range(0,len(arr),16):
        for ioff in range(16):
            i = istart + ioff
            if i >= len(arr):
                break
            if ioff == 0:
                pass
            fh.write("%5i"%(arr[i]))
        fh.write("\n")




# Compiled once - ReadGauEsp is called for every multi-RESP orientation log.
_ATOMIC_CENTER_LINE = None
_FIT_CENTER_LINE = None
_ESP_LINE = None


def _gau_esp_regexes():
    global _ATOMIC_CENTER_LINE, _FIT_CENTER_LINE, _ESP_LINE
    if _ATOMIC_CENTER_LINE is None:
        import re

        _ATOMIC_CENTER_LINE = re.compile(
            r" +Atomic Center[ 0-9]+ is at +([\-0-9]{1,3}\.[0-9]{6}) *([\-0-9]{1,3}\.[0-9]{6}) *([\-0-9]{1,3}\.[0-9]{6})"
        )
        _FIT_CENTER_LINE = re.compile(
            r" +ESP Fit Center[ 0-9]+ is at +([\-0-9]{1,3}\.[0-9]{6}) *([\-0-9]{1,3}\.[0-9]{6}) *([\-0-9]{1,3}\.[0-9]{6})"
        )
        _ESP_LINE = re.compile(r"[ 0-9]+ Fit +([\-\.0-9]+)")
    return _ATOMIC_CENTER_LINE, _FIT_CENTER_LINE, _ESP_LINE


def _coords_after_is_at(line: str):
    """Parse ``is at  x y z`` without a regex when the line shape is standard."""
    idx = line.find("is at")
    if idx < 0:
        return None
    try:
        vals = [float(tok) for tok in line[idx + 5 :].split()]
    except ValueError:
        return None
    if len(vals) != 3:
        return None
    return vals


def ReadGauEsp(fname):
    """ Read the ESP from a Gaussian output file
    
    Parameters
    ----------
    fname : str
        The Gaussian output file
    
    Returns
    -------
    crds : list of float
        The coordinates of the atoms
    
    """
    import os

    if not os.path.exists(fname):
        raise Exception("ReadGauEsp file not found: %s"%(fname))

    atomic_center_line, fit_center_line, esp_line = _gau_esp_regexes()

    crds=[]
    pts=[]
    esp=[]

    with open(fname, "r") as fh:
        for line in fh:
            # Substring gates avoid three regex matches on every log line
            # (x ~28 orientation files in FreeLigand multi-RESP).
            if "Atomic Center" in line:
                vals = _coords_after_is_at(line)
                if vals is None:
                    m = atomic_center_line.match(line)
                    if m is not None:
                        vals = [float(m.group(1)), float(m.group(2)), float(m.group(3))]
                if vals is not None:
                    crds.extend(vals)
                continue
            if "ESP Fit Center" in line:
                vals = _coords_after_is_at(line)
                if vals is None:
                    m = fit_center_line.match(line)
                    if m is not None:
                        vals = [float(m.group(1)), float(m.group(2)), float(m.group(3))]
                if vals is not None:
                    pts.extend(vals)
                continue
            if " Fit " in line:
                parts = line.split()
                if len(parts) >= 3 and parts[1] == "Fit":
                    try:
                        esp.append(float(parts[-1]))
                        continue
                    except ValueError:
                        pass
                m = esp_line.match(line)
                if m is not None:
                    esp.append(float(m.group(1)))
                    continue
    if len(pts) != 3*len(esp):
        raise Exception("# pts (%i) != # esp values (%i) in %s"%(len(pts)//3,len(esp),fname))

    return crds,pts,esp


def RunGauEsp(atn,crds,charge,mult,basename,nproc):
    """ Run Gaussian to calculate the ESP
    
    Parameters
    ----------
    atn : list of str
        The atomic symbols
    crds : list of float
        The atomic coordinates
    charge : int
        The charge
    mult : int
        The multiplicity
    basename : str
        The base name for the output files
    nproc : int
        The number of processors to use
    
    """
    import subprocess
    import warnings
    warnings.warn("Careful, running gaussian task as a subprocess!!!")
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
    """ Read the atomic coordinates, charge, and multiplicity from a Gaussian output file
    
    Parameters
    ----------
    fname : str
        The Gaussian output file
    
    Returns
    -------
    atn : list of str
        The atomic symbols
    crds : list of float
        The atomic coordinates
    charge : int
        The charge
    mult : int
        The multiplicity
    """
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
    except Exception as exc:
        print("Could not process gaussian file '%s'"%(fname))
        print("This is the archive:")
        for i,sec in enumerate(secs):
            subs = sec.split("\\")
            for j,sub in enumerate(subs):
                vals = sub.split(",")
                print("%2i %2i %s"%(i,j,str(vals)))
        raise ValueError(
            f"Could not parse Gaussian archive in {fname!r}: {exc}"
        ) from exc
    
    return atn,crds,charge,mult

def ReadOrMakeGauEsp( mdout, nproc=4 ):
    """ Read the ESP from a Gaussian output file, or generate it if it doesn't exist
    
    Parameters
    ----------
    mdout : str
        The Gaussian output file
    nproc : int
        The number of processors to use
    
    Returns
    -------
    crds : list of list of float
        The coordinates of the atoms
    pts : list of list of float
        The coordinates of the grid points
    esp : list of float
        The ESP values
    """
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