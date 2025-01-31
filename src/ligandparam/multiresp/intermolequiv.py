from collections import defaultdict as ddict

from ligandparam.multiresp import respfunctions, parmhelper

class IntermolEquiv(object):
    """ Class to handle the equivalent atoms between molecules """
    def __init__( self, equiv_mask, unique_residues ):
        """ Initialize the IntermolEquiv class.
        
        Parameters
        ----------
        equiv_mask : str
            The equivalent mask
        unique_residues : bool
            If True, use unique residues
        
        """
        self.equiv_mask = equiv_mask
        self.mols = []
        self.state_map = ddict( str )
        self.states = []
        self.unique_residues=unique_residues

    def append( self, state, mdout ):
        """ Append a state to the IntermolEquiv class
        
        Parameters
        ----------
        state : str
            The state
        mdout : str
            The mdout file
        """
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
        """ Get the unique atoms
        
        Parameters
        ----------
        mask : str, optional
            The mask. Default is "@*"
        
        Returns
        -------
        uatoms : list of int
            The unique atoms
        """
        uatoms = []
        if mask is not None:
            for mol in self.mols:
                atoms = mol.get_selected_atoms( mask )
                for atom in atoms:
                    if atom not in uatoms:
                        uatoms.append(atom)
        return uatoms

    def get_unique_atoms_from_state(self,state,mask="@*"):
        """ Get the unique atoms from a state
        
        Parameters
        ----------
        state : str
            The state
        mask : str, optional
            The mask. Default is "@*"
        """
        uatoms = []
        if mask is not None:
            for mol in state:
                atoms = mol.get_selected_atoms( mask )
                for atom in atoms:
                    if atom not in uatoms:
                        uatoms.append(atom)
        return uatoms
    

    def print_equiv( self, fh ):
        """ Print the equivalent atoms
        
        Parameters
        ----------
        fh : file
            The file handle
        """

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
                respfunctions.WriteArray8(fh,idxs)
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
                    respfunctions.WriteArray8(fh,idxs)
        fh.write("\n")


class RespMol(object):
    """ Class to handle the RESP molecule """
    def __init__(self,state,mdout,unique_residues):
        """ Initialize the RespMol class
        
        Parameters
        ----------
        state : str
            The state
        mdout : str
            The mdout file
        unique_residues : bool
            If True, use unique residues
        """
        self.state = state
        self.mdout = mdout
        self.atoms = []
        # get the atoms
        for fitidx in range(self.state.nquant):
            parmidx = self.state.fit2parm_map[ fitidx ]
            rname = respfunctions.GetResidueNameFromAtomIdx(self.state.parm,parmidx,unique_residues)
            ridx  = self.state.parm.atoms[parmidx].residue.idx+1
            aname = self.state.parm.atoms[parmidx].name
            if unique_residues:
                self.atoms.append( "%s%i:%s"%( rname, ridx, aname ) )
            else:
                self.atoms.append( "%s:%s"%( rname, aname ) )
        # Get the link atoms
        for linkidx in range(self.state.nlink):
            fitidx = linkidx + self.state.nquant
            parmidx = self.state.fit2parm_map[ fitidx ]
            rname = respfunctions.GetResidueNameFromAtomIdx(self.state.parm,parmidx,unique_residues)
            ridx  = self.state.parm.atoms[parmidx].residue.idx+1
            aname = self.state.parm.atoms[parmidx].name
            if unique_residues:
                self.atoms.append( "%s%i:%s:%i"%( rname, ridx, aname, linkidx+1 ) )
            else:
                self.atoms.append( "%s:%s:%i"%( rname, aname, linkidx+1 ) )


    def find_fitidx(self,name):
        """ Find the fit index
        
        Parameters
        ----------
        name : str
            The name
        
        Returns
        -------
        fitidx : int
            The fit index
        """
        fitidx=-1
        if name in self.atoms:
            fitidx = self.atoms.index( name )
        return fitidx

    def get_selected_atoms( self, sele ):
        """ Get the selected atoms 
        
        Parameters
        ----------
        sele : str
            The selection
        """
        parmidxs = parmhelper.GetSelectedAtomIndices(self.state.parm,sele)
        atoms = []
        for parmidx in parmidxs:
            for fitidx in self.state.parm2fit_map[parmidx]:
                atoms.append( self.atoms[ fitidx ] )
        return atoms
