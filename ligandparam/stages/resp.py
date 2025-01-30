import glob
from typing import Union
from pathlib import Path

from ligandparam.stages.abstractstage import AbstractStage
from ligandparam.interfaces import Antechamber

from ligandparam.multiresp import parmhelper
from ligandparam.multiresp.residueresp import ResidueResp
from ligandparam.log import get_logger


class StageLazyResp(AbstractStage):
    """ This class runs a 'lazy' resp calculation based on only
        a single gaussian output file. """

    def __init__(self, stage_name: str, in_filename: Union[Path, str], cwd: Union[Path, str], *args, **kwargs) -> None:
        super().__init__(stage_name, in_filename, cwd, *args, **kwargs)
        self.in_gaussian_log = Path(in_filename)
        self.add_required(self.in_gaussian_log)
        self.out_mol2 = Path(kwargs["out_mol2"])

        self.net_charge = getattr(kwargs, 'net_charge', 0.0)
        self.atom_type = getattr(kwargs, 'atom_type', "gaff2")


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
        ante = Antechamber(cwd=self.cwd, logger=self.logger)
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

    def __init__(self, stage_name: str, in_filename: Union[Path, str], cwd: Union[Path, str], *args, **kwargs) -> None:
        super().__init__(stage_name, in_filename, cwd, *args, **kwargs)
        # in_gaussian_log could be any log file? This is where the assumption of each stage taking a filename breaks down.
        # This should only take a working dir.
        self.in_gaussian_log = Path(in_filename)
        # self.add_required(self.in_gaussian_log)
        self.out_mol2 = Path(kwargs["out_mol2"])
        self.add_required(self.out_mol2)
        self.out_respfit = Path(kwargs["out_respfit"])
        
        self.net_charge = getattr(kwargs, 'net_charge', 0.0)


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
        gaussian_out_files = Path(self.cwd, f"{self.in_gaussian_log.stem}_*.log")
        model.add_state(self.out_mol2.name, self.gauss_logmol2_fname, glob.glob(gaussian_out_files), qmmask="@*")
        model.multimolecule_fit(True)
        model.perform_fit("@*", unique_residues=False)
        with open(self.out_respfit, "w") as f:
            model.print_resp(fh=f)

        return

    def _clean(self):
        """ Clean the files generated during the stage. """
        raise NotImplementedError("clean method not implemented")
