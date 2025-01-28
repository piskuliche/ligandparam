import glob
from typing import Union
from pathlib import Path

from ligandparam.stages.abstractstage import AbstractStage
from ligandparam.interfaces import Antechamber

from ligandparam.multiresp import parmhelper
from ligandparam.multiresp.residueresp import ResidueResp


class StageLazyResp(AbstractStage):
    """ This class runs a 'lazy' resp calculation based on only
        a single gaussian output file. """

    def __init__(self, stage_name: str, name: Union[Path, str], cwd: Union[Path, str], *args, **kwargs) -> None:
        super().__init__(stage_name, name, cwd, *args, **kwargs)
        self.net_charge = getattr(kwargs, 'net_charge', 0.0)
        self.atom_type = getattr(kwargs, 'atom_type', "gaff2")
        self.in_gaussian_log = Path(kwargs["in_gaussian_log"])
        self.add_required(self.in_gaussian_log)
        self.out_mol2 = Path(kwargs["out_mol2"])


    def _append_stage(self, stage: "AbstractStage") -> "AbstractStage":
        """ Appends the stage. """
        return stage

    def _execute(self, dry_run=False):
        """ Execute antechamber to convert the gaussian output to a mol2 file. 
        
        Parameters
        ----------
        dry_run : bool, optional
            If True, the stage will not be executed, but the function will print the commands that would
        """
        print(f"Executing {self.name} with netcharge={self.net_charge}")
        ante = Antechamber(cwd=self.cwd)
        ante.call(i=self.in_gaussian_log, fi='gout',
                  o=self.out_mol2, fo='mol2',
                  gv=0, c='resp',
                  nc=self.net_charge,
                  at=self.atom_type, dry_run=dry_run)
        return

    def _clean(self):
        """ Clean the files generated during the stage. """
        raise NotImplementedError("clean method not implemented")


class StageMultiRespFit(AbstractStage):
    """ This class runs a multi-state resp fitting calculation, based on 
        multiple gaussian output files. 
        
        TODO: Implement this class.
        TODO: Implement the clean method.
        TODO: Implement the execute method.
        TODO: Add a check that a multistate resp fit is possible. 
        """

    def __init__(self, stage_name: str, name: Union[Path, str], cwd: Union[Path, str], *args, **kwargs) -> None:
        super().__init__(stage_name, name, cwd, *args, **kwargs)
        self.net_charge = getattr(kwargs, 'net_charge', 0.0)
        self.gauss_logmol2_fname = Path(self.cwd, f"{self.name.stem}.log.mol2")
        self.add_required(self.gauss_logmol2_fname)

    def _append_stage(self, stage: "AbstractStage") -> "AbstractStage":
        """ Appends the stage. """
        return stage

    def _execute(self, dry_run=False):
        """Execute a multi-state respfitting calculation.

        if __name__ == "__main__":
        comp = parmutils.BASH( 12 )
        model = rf.ResidueResp( comp, 1 )


        model.add_state( "$base", "$base.log.mol2", glob.glob("gaussianCalcs/$base_*.log"), qmmask="@*" )


        model.multimolecule_fit(True)
        model.perform_fit("@*",unique_residues=False)
        #model.preserve_residue_charges_by_shifting()
        model.print_resp()


        Parameters
        ----------
        dry_run : bool, optional
            If True, the stage will not be executed, but the function will print the commands that would

        """
        comp = parmhelper.BASH(12)
        model = ResidueResp(comp, 1)
        gaussian_out_files = Path(self.cwd, "gaussianCalcs", f"{self.name.stem}_*.log")
        model.add_state(self.name.name, self.gauss_logmol2_fname, glob.glob(gaussian_out_files), qmmask="@*")
        model.multimolecule_fit(True)
        model.perform_fit("@*", unique_residues=False)
        with open(self.cwd / "respfit.out", "w") as f:
            model.print_resp(fh=f)

        return

    def _clean(self):
        """ Clean the files generated during the stage. """
        raise NotImplementedError("clean method not implemented")
