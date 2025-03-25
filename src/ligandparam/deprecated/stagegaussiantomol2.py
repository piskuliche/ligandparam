import warnings

import MDAnalysis as mda

from ligandparam.abstractstage import AbstractStage
from ligandparam.interfaces import Antechamber


class StageGaussiantoMol2(AbstractStage):
    def __init__(self, name, base_cls=None, dry_run=None) -> None:
        self.name = name
        self.base_cls = base_cls
        self.dry_run = dry_run

    def _append_stage(self, stage: "AbstractStage") -> "AbstractStage":
        return stage

    def execute(self, dry_run=False, nproc=1, mem=1) -> Any:
        if self.dry_run is not None:
            dry_run = self

        # Convert from gaussian to mol2
        ante = Antechamber()
        ante.call(
            i=self.base_cls.name + ".log",
            fi="gout",
            o=self.base_cls.name + ".tmp1.mol2",
            fo="mol2",
            pf="y",
            at=self.base_cls.atom_type,
            an="no",
            nc=self.net_charge,
            run=(not dry_run),
        )

        # Assign the charges
        if not dry_run:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                u1 = mda.Universe(self.base_cls.name + ".antechamber.mol2")
                u2 = mda.Universe(self.base_cls.name + ".tmp1.mol2")
                assert len(u1.atoms) == len(u2.atoms), "Number of atoms in the two files do not match"

                u2.atoms.charges = u1.atoms.charges

                ag = u2.select_atoms("all")
                ag.write(self.base_cls.name + ".tmp2.mol2")

        # Use antechamber to clean up the mol2 format
        ante = Antechamber(self.cwd)
        ante.call(
            i=self.base_cls.name + ".tmp2.mol2",
            fi="mol2",
            o=self.base_cls.name + ".log.mol2",
            fo="mol2",
            pf="y",
            at=self.base_cls.atom_type,
            an="no",
            nc=self.net_charge,
            run=(not dry_run),
        )

        return
