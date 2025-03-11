""" This module is used to perform resp fitting for a molecule, based on multiple gaussian output files.

This was initially developed by Prof. Timothy Giese at Rutgers University, and was reproduced more or less in
its entirety here (with modifications to the code for readibility and documentation). This code was initially meant
for more options than are currently used by the ligand-param project, so some of the code may be redundant or not used.

"""


import sys

from ligandparam.multiresp.endstate import EndState
from ligandparam.multiresp.intermolequiv import IntermolEquiv
from ligandparam.multiresp.functions import WriteFitSh

class ResidueResp(object):
    """ This class is used to perform a residue-based resp fit. """
    def __init__(self,comp,ires,theory="PBE0",basis="6-31G*",maxgrad=1.e-6,etol=0.0001,fitgasphase=False):
        """ Initialize the ResidueResp object
        
        Parameters
        ----------
        comp : BASH
            The computer object
        ires : int
            The residue number
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
        """ Add a state to the residue resp object 

        TODO: Confirm meaning with TIM
        
        Parameters
        ----------
        prefix : str
            The prefix of the state (e.g. the residue name)
        parmfile : str
            The name of the file (either parm or mol2)
        rstfiles : list of str
            A list of files (could be rst or gaussian log files)
        qmmask : str, optional
            The quantum mask. Default is None
        """
        # Add state to the residue resp object
        self.states.append( EndState( parmfile, rstfiles, self.comp, self.ires, qmmask=qmmask,
                                      theory=self.theory,basis=self.basis,maxgrad=self.maxgrad,
                                      etol=self.etol,fitgasphase=self.fitgasphase ) )
        self.states[-1].prefix = prefix

    def clear_charge_data(self):
        """ Clear the charge data from the states"""
        for state in self.states:
            state.clear_charge_data()

    def write_mdin(self):
        """ Write the mdin file """
        for state in self.states:
            state.write_mdin()

    def read_respfile(self):
        """ Read the resp file """
        if self.multifit:
            self.clear_charge_data()
            out = open("%s.resp.out"%(self.base),"r")
            for state in self.states:
                mdouts = state.get_mdouts()
                for mdout in mdouts:
                    state.read_next_resp( out )
                

    def preserve_residue_charges_by_shifting(self):
        """ Preserve the residue charges by shifting for each state """
        for state in self.states:
            state.preserve_residue_charges_by_shifting()
              
    def preserve_mm_charges_by_shifting(self,mm_mask):
        """ Preserve the mm charges by shifting for each state """
        for state in self.states:
            state.preserve_mm_charges_by_shifting(mm_mask)
                
    def print_resp(self,fh=sys.stdout):
        """ Print the resp file for each state"""
        for state in self.states:
            state.print_resp(state.prefix,fh=fh)

    def apply_backbone_restraint(self,value=True):
        """ Apply the backbone restraint to each state
        
        Parameters
        ----------
        value : bool, optional
            If True, apply the backbone restraint. Default is True
            
        """
        for state in self.states:
            state.backbone_restraint=value

    def apply_sugar_restraint(self,value=True):
        """ Apply the sugar restraint to each state 

        Parameters
        ----------
        value : bool, optional
            If True, apply the sugar restraint. Default is True
        
        """
        for state in self.states:
            state.sugar_restraint=value

    def apply_nucleobase_restraint(self,value=True):
        """ Apply the nucleobase restraint to each state
        
        Parameters
        ----------
        value : bool, optional
            If True, apply the nucleobase restraint. Default is True
        
        """
        for state in self.states:
            state.nucleobase_restraint=value

    def apply_equiv_hydrogens(self,value=True):
        """ Apply the equivalent hydrogens to each state
        
        Parameters
        ----------
        value : bool, optional
            If True, apply the equivalent hydrogens. Default is True"""
        for state in self.states:
            state.equiv_hydrogens=value

    def apply_equiv_nonbridge(self,value=True):
        """ Apply the equivalent nonbridge to each state
        
        Parameters
        ----------
        value : bool, optional
            If True, apply the equivalent nonbridge. Default is True
        """
        for state in self.states:
            state.equiv_nonbridge=value
            
    def multimolecule_fit(self,value=True):
        """ Perform a multi-molecule fit 
        
        Parameters
        ----------
        value : bool, optional
            If True, perform a multi-molecule fit. Default is True
        """
        self.multifit=value
        for state in self.states:
            state.multifit = value

    def perform_fit(self,equiv_mask="@P,OP*,*'",unique_residues=True):
        """ Perform the fit for each state

        Parameters
        ----------
        equiv_mask : str, optional
            The equivalent mask. Default is "@P,OP*,*'"
        unique_residues : bool, optional
            If True, use unique residues. Default is True
        """
        if self.multifit:
            equiv = IntermolEquiv(equiv_mask,unique_residues)
            nmols = 0
            for state in self.states:
                nmols += len(state.get_mdouts())
                
            # print("# ResidueResp.perform_fit writing multifit %s.resp.inp"%(self.base))
            
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


