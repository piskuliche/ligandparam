import sys

from ligand_param.multiresp.endstate import EndState
from ligand_param.multiresp.intermolequiv import IntermolEquiv
from ligand_param.multiresp.functions import WriteFitSh

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