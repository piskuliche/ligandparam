import parmed

def WriteFitSh(base):
    import subprocess
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
          

def WriteMaskedFrcmod(param,aidxs,native_frcmod,changes_frcmod,with_mass=True,with_nonb=True):
    """
Converts selected atoms to new atom-types and
writes two frcmod files, one containing the
those parameters involving non-substituted
atoms, and the other for those involving the
substituted atoms.

Parameters
----------
param : Amber parm file object
    The parameter file object
aidxs : list of int
    List of atoms that will use new atom-type parameters
native_frcmod : str
    The filename of the nonsubstituted parameters
changes_frcmod : str
    The filename listing the new atom-type parameters
with_mass : bool, optional
    If True, include mass parameters. Default is True
with_nonb : bool, optional
    If True, include nonbonded parameters. Default is True
    """

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


def WriteFrcmodObj(self,native_frcmod,angfact=1.0,uniqueparams=False,selected_names=None,changed_frcmod=None,with_mass=True,with_nonb=True):

    """ Writes an frcmod file from a parmed object
    
    Parameters
    ----------
    self : parmed object
        The parmed object
    native_frcmod : str
        The filename of the nonsubstituted parameters
    angfact : float, optional
        The angle factor. Default is 1.0
    uniqueparams : bool, optional
        If True, only write unique parameters. Default is False
    selected_names : dict, optional
        The selected atom names. Default is None
    changed_frcmod : str, optional
        The filename listing the new atom-type parameters. Default is None
    with_mass : bool, optional
        If True, include mass parameters. Default is True
    with_nonb : bool, optional
        If True, include nonbonded parameters. Default is True
    
    """
    
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
    
    Raises
    ------
    ValueError
        If the sample covariance is zero, the statistical inefficiency cannot be computed.

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
        raise ValueError('Sample covariance sigma_AB^2 = 0 -- cannot compute statistical inefficiency')

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