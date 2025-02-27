import os
import copy
import sys

from collections import defaultdict as ddict
from io import StringIO as StringIO


from ligandparam.multiresp import parmhelper, mdinutils, respfunctions, functions
from ligandparam.multiresp.intermolequiv import IntermolEquiv

class EndState(object):
    """ Class to handle the end state of a residue."""
    def __init__(self,parmfile,rstfiles,comp,ires,qmmask=None,theory="PBE0",basis="6-31G*",maxgrad=1.e-6,etol=1.e-4,fitgasphase=False):
        """ Initialize the EndState object

        Parameters
        ----------
        parmfile : str
            The name of the parameter file (either parm or mol2)
        rstfiles : list of str
            A list of files (could be rst or gaussian log files, specifically with .log or .rst7 extensions)
        comp : BASH
            The computer object
        ires : int
            The residue number
        qmmask : str, optional
            The quantum mask. Default is None
        theory : str
            The quantum theory used. Default is PBE0
        basis : str
            The basis set used. Default is 6-31G*
        maxgrad : float
            The maximum gradient. Default is 1.e-6
        etol : float
            The energy tolerance. Default is 0.0001
        fitgasphase : bool
            If True, fit the gas phase. Default is False
        
        """
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
            self.parm = parmhelper.OpenParm( parmfile, xyz=None )
        else:
            self.parm = parmhelper.OpenParm( parmfile, xyz=rstfiles[0] )
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

        self.fsys = parmhelper.FragmentedSys( self.parm, self.comp )
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
        """ Set the equivalent atoms for the residue. """
        idxs = parmhelper.GetSelectedAtomIndices(self.parm,sele)
        self.equiv_atoms.append(idxs)
        
        
    def add_group_restraint(self,cmask,q=None):
        """ Add a group restraint to the residue. 
        
        Parameters
        ----------
        cmask : str
            The mask for the group
        q : float, optional
            The charge. Default is None
        
        """
        if q is None:
            grp = parmhelper.GetSelectedAtomIndices(self.parm,cmask)
            grp_charges = [ self.parm.atoms[a].charge for a in grp ]
            q  = sum( grp_charges )
        self.grp_restraints.append( (cmask,q) )
        

    def apply_backbone_restraint(self,value=True):
        """ Apply a backbone restraint to the residue."""
        self.backbone_restraint=value


    def apply_sugar_restraint(self,value=True):
        """ Apply a sugar restraint to the residue."""
        self.sugar_restraint=value


    def apply_nucleobase_restraint(self,value=True):
        """  Apply a nucleobase restraint to the residue."""
        self.nucleobase_restraint=value


    def multimolecule_fit(self,value=True):
        """ Perform a multi-molecule fit"""
        self.multifit = value


    def write_mdin(self):
        """ Write the mdin files for the residue."""

        has_g09=False
        for rstfile in self.rstfiles:
            if ".log" in rstfile:
                has_g09=True
        if has_g09:
            return
            
        self.fsys.check_overlaps()
        self.fsys.redistribute_residue_charges()
        if not os.path.isfile( self.modified_parmfile ):
            parmhelper.SaveParm( self.fsys.parmobj, self.modified_parmfile )
        if self.fitgasphase:
            if not os.path.isfile( self.gasphase_parmfile ):
                gp = parmhelper.CopyParm( self.fsys.parmobj )
                qmidxs = parmhelper.GetSelectedAtomIndices( gp, self.fsys.get_noshake_selection() )
                for a in gp.atoms:
                    if a.idx not in qmidxs:
                        a.charge = 0
                parmhelper.SaveParm( gp, self.gasphase_parmfile )


            
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
        """ Clear the charge data for the residue."""
        self.charge_data = ddict( list )

                
    def read_next_resp(self,fh):
        """ Read the next set of RESP charges from the file handle
        
        Parameters
        ----------
        fh : file handle
            The file handle to read the charges from
        
        """
        import sys
        if True:
            fit_charges = respfunctions.ReadNextRespCharges( fh )
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
        """ Read the resp file """
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
        """ Preserve the residue charges by shifting """
        for res in self.parm.residues:
            mmq = 0.
            q = 0.
            n = 0
            for a in res.atoms:
                mmq += self.mmcharges[a.idx]
                if a.idx in self.frag.atomsel:
                    avg,std,err = respfunctions.GetAvgStdDevAndErr(self.charge_data[a.idx])
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
        """ Preserve the MM charges by shifting 
        
        Parameters
        ----------
        mm_mask : str
            The mask for the MM charges
        
        """
        reset_sele = parmhelper.GetSelectedAtomIndices(self.parm,mm_mask)
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
                        avg,std,err = respfunctions.GetAvgStdDevAndErr(self.charge_data[a.idx])
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
        """ Print the RESP charges
        
        Parameters
        ----------
        prefix : str, optional
            The prefix. Default is ""
        fh : file handle, optional
            The file handle. Default is sys.stdout
        """
        for a in self.frag.atomsel:
            q,stddev,stderr = respfunctions.GetAvgStdDevAndErr(self.charge_data[a])
            mm = self.mmcharges[a]
            #res = GetResidueNameFromAtomIdx(self.parm,a,self.unique_residues)
            res = "%3i"%(self.parm.atoms[a].residue.idx+1)
            atm = self.parm.atoms[a].name
            dct = "QQS[\"%s\"][%s][\"%s\"]"%(prefix,res,atm)
            fh.write("    %-27s = %10.6f # %4s %10.6f %4i %8.4f %8.4f\n"%(dct,q,self.parm.atoms[a].residue.name,mm,len(self.charge_data[a]),stderr,stddev))
        fh.write("\n")


    def get_mdouts(self):
        """ Get the mdout files 
        
        Returns
        -------
        mdouts : list of str
            The list of mdout files
        
        """
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
        """ Perform the fit

        Parameters
        ----------
        unique_residues : bool, optional
            If True, use unique residues. Default is True
        
        """
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
            functions.WriteFitSh( os.path.join(self.topdir,self.base) )
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
                
                functions.WriteFitSh( base )
        self.read_respfile()

    def get_resp_header(self):
        """ Get the RESP header """
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
            equiv_grps += respfunctions.GetEquivHydrogens(self.parm,self.frag.atomsel)
        if self.equiv_nonbridge:
            equiv_grps += respfunctions.GetEquivNonbridgeOxygens(self.parm,self.frag.atomsel)
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
            grp = parmhelper.GetSelectedAtomIndices(self.parm,cmask)
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
            respfunctions.WriteArray8(fh,arr)

        for cmask,grp_charge in self.grp_restraints:
            grp = parmhelper.GetSelectedAtomIndices(self.parm,cmask)
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
            respfunctions.WriteArray8(fh,arr)
            
        fh.write("\n")

        return fh.getvalue()

    def append_esp(self,fh,mdout):
        """ Append the ESP data to the file handle
        
        Parameters
        ----------
        fh : file handle
            The file handle to append the ESP data to
        mdout : str
            The mdout file
        """
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
            crds,pts,vals = respfunctions.ReadOrMakeGauEsp( mdout, max(4,self.comp.num_nodes) )
            nat=len(crds)//3
            npt=len(vals)
            fh.write("%5i%5i\n"%(nat,npt))
            au = 1.88972613373
            for i in range(nat):
                fh.write("%17s%16.7f%16.7f%16.7f\n"%("",au*crds[0+i*3],au*crds[1+i*3],au*crds[2+i*3]))
            for i in range(npt):
                fh.write(" %16.7f%16.7f%16.7f%16.7f\n"%(vals[i],au*pts[0+i*3],au*pts[1+i*3],au*pts[2+i*3]))

            
    def append_esps(self,fh):
        """ Append the ESP data to the file handle"""
        for mdout in self.get_mdouts():
            self.append_esp(fh,mdout)