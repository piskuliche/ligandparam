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