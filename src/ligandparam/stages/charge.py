from typing import Optional, Union, Any
from typing_extensions import override
from pathlib import Path
import warnings

import numpy as np
import MDAnalysis as mda

from ligandparam.stages.abstractstage import AbstractStage
from ligandparam.interfaces import Antechamber
from ligandparam.io.coordinates import Mol2Writer


class StageUpdateCharge(AbstractStage):
    """
    Create a new mol2 file with updated charges.

    Parameters
    ----------
    stage_name : str
        The name of the stage.
    main_input : Union[Path, str]
        Path to the input mol2 file.
    cwd : Union[Path, str]
        Current working directory.
    charge_source : str
        Path to the file containing charges.
    charge_column : int, optional
        Column index in charge_source to use for charges (default: 3).
    out_mol2 : str
        Path to the output mol2 file.
    net_charge : float, optional
        Net charge for the molecule (default: 0.0).
    atom_type : str, optional
        Atom type (default: 'gaff2').

    Attributes
    ----------
    in_mol2 : Path
        Path to the input mol2 file.
    charge_source : str
        Path to the file containing charges.
    charge_column : int
        Column index in charge_source to use for charges.
    out_mol2 : Path
        Path to the output mol2 file.
    tmp_mol2 : Path
        Path to the temporary mol2 file.
    net_charge : float
        Net charge for the molecule.
    atom_type : str
        Atom type.
    """

    @override
    def __init__(self, stage_name: str, main_input: Union[Path, str], cwd: Union[Path, str], *args, **kwargs) -> None:
        super().__init__(stage_name, main_input, cwd, *args, **kwargs)
        self.in_mol2 = Path(main_input)
        self.charge_source = kwargs["charge_source"]
        self.charge_column = kwargs.get("charge_column", 3)
        self.out_mol2 = Path(kwargs["out_mol2"])
        self.tmp_mol2 = self.cwd / f"{self.out_mol2.stem}_tmp_update.mol2"  # tmpresp
        self.net_charge = kwargs.get("net_charge", 0.0)
        self.atom_type = kwargs.get("atom_type", "gaff2")

        self.add_required(Path(self.in_mol2))
        self.add_required(Path(self.charge_source))

        return

    def _append_stage(self, stage: "AbstractStage") -> "AbstractStage":
        return stage

    def execute(self, dry_run=False, nproc: Optional[int] = None, mem: Optional[int] = None) -> Any:
        super()._setup_execution(dry_run=dry_run, nproc=nproc, mem=mem)
        # Supress the inevitable mol2 file warnings.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            if Path(self.charge_source).exists():
                charges = np.genfromtxt(self.charge_source, usecols=self.charge_column, unpack=True)
            else:
                raise FileNotFoundError(f"File {self.charge_source} not found.")

            if not dry_run:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    u = mda.Universe(self.in_mol2, format="mol2")
                if len(charges) != len(u.atoms):
                    raise ValueError("Error: Number of charges does not match the number of atoms.")
                u.atoms.charges = charges
                # Write the Mol2 temporary file
                Mol2Writer(u, self.tmp_mol2, selection="all").write()

            ante = Antechamber(cwd=self.cwd, logger=self.logger, nproc=self.nproc)
            ante.call(
                i=self.tmp_mol2, fi="mol2", o=self.out_mol2, fo="mol2", pf="y", at=self.atom_type, an="no",
                nc=self.net_charge, dry_run=dry_run
            )

        return

    def _clean(self):
        raise NotImplementedError("clean method not implemented")


class StageNormalizeCharge(AbstractStage):
    """
    Normalize the charges in a mol2 file to the specified net charge.

    This class works by calculating the charge difference, and then normalizing the charges
    based on the overall precision that you select, by adjusting each atom charge by the precision
    until the charge difference is zero.

    Parameters
    ----------
    stage_name : str
        The name of the stage.
    main_input : Union[Path, str]
        Path to the input mol2 file.
    cwd : Union[Path, str]
        Current working directory.
    out_mol2 : str
        Path to the output mol2 file.
    atom_type : str, optional
        Atom type (default: 'gaff2').
    net_charge : float, optional
        Net charge for the molecule (default: 0.0).
    precision : float, optional
        Precision for charge normalization (default: 0.0001).

    Attributes
    ----------
    in_mol2 : Path
        Path to the input mol2 file.
    out_mol2 : Path
        Path to the output mol2 file.
    tmp_mol2 : Path
        Path to the temporary mol2 file.
    atom_type : str
        Atom type.
    net_charge : float
        Net charge for the molecule.
    precision : float
        Precision for charge normalization.
    decimals : int
        Number of decimals for rounding charges.
    """

    def __init__(self, stage_name: str, main_input: Union[Path, str], cwd: Union[Path, str], *args, **kwargs) -> None:
        super().__init__(stage_name, main_input, cwd, *args, **kwargs)
        self.in_mol2 = Path(main_input)

        self.out_mol2 = Path(kwargs["out_mol2"])
        self.tmp_mol2 = self.cwd / f"{self.in_mol2.stem}_tmp_norm.mol2"

        self.atom_type = kwargs.get("atom_type", "gaff2")
        self.net_charge = kwargs.get("net_charge", 0.0)
        self.precision = kwargs.get("precision", 0.0001)
        try:
            self.decimals = len(str(self.precision).split(".")[1])
        except IndexError:
            raise ValueError(f"ERROR: Invalid precision: {self.precision}. It should be a float between 0 and 0.1")

        self.add_required(self.in_mol2)

    def _append_stage(self, stage: "AbstractStage") -> "AbstractStage":
        return stage

    def execute(self, dry_run=False, nproc: Optional[int] = None, mem: Optional[int] = None) -> Any:
        """
        Execute the charge normalization stage.

        Parameters
        ----------
        dry_run : bool, optional
            If True, the stage will not be executed, but the function will print the commands that would be run.
        nproc : int, optional
            Number of processors to use.
        mem : int, optional
            Amount of memory to use (in GB).

        Raises
        ------
        ValueError
            If the charge normalization fails.

        Notes
        -----
        TODO: Check what happens when netcharge is nonzero.
        TODO: Check what happens when charge difference is larger than the number of atoms.
        """
        super()._setup_execution(dry_run=dry_run, nproc=nproc, mem=mem)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.logger.debug("Checking charges")
            self.logger.debug(f"Normalizing charges to {self.net_charge}")
            self.logger.debug(f"Precision {self.precision} with {self.decimals} decimals")

            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                u = mda.Universe(self.in_mol2, format="mol2")
            rounded_charges, total_charge, charge_difference = self.check_charge(u.atoms.charges)

            if not np.isclose(total_charge, self.net_charge, rtol=1e-10):
                self.logger.info("Normalizing charges")
                new_charges = self.normalize(rounded_charges, charge_difference)
                _, new_total, new_diff = self.check_charge(new_charges)
                if np.isclose(new_total, self.net_charge, rtol=1e-10):
                    u.atoms.charges = new_charges
                else:
                    raise ValueError(f"Error: Charge normalization failed, new charge: {new_total}.")
            else:
                self.logger.info("Charges are already normalized")
            if not dry_run:
                Mol2Writer(u, self.tmp_mol2, selection="all").write()

                ante = Antechamber(cwd=self.cwd, logger=self.logger, nproc=self.nproc)
                ante.call(
                    i=self.tmp_mol2, fi="mol2", o=self.out_mol2, fo="mol2", pf="y", at=self.atom_type, an="no",
                    nc=self.net_charge, dry_run=dry_run
                )

    def _clean(self):
        raise NotImplementedError("clean method not implemented")

    def normalize(self, charges, charge_difference):
        """
        Normalize the charges to the net charge.

        Parameters
        ----------
        charges : np.ndarray
            Array of atomic charges.
        charge_difference : float
            The charge difference to be corrected.

        Returns
        -------
        np.ndarray
            The normalized charges.
        """

        count = np.round(np.abs(charge_difference) / self.precision)
        adjust = np.round(charge_difference / count, self.decimals)
        natoms = len(charges)
        # Choosing charges closest to zero.
        sorted_indices = np.argsort(np.abs(charges))
        # Flip the order to choose the largest charges first.
        sorted_indices = sorted_indices[::-1]
        for i in range(int(count)):
            atom_idx = i % natoms
            charges[sorted_indices[atom_idx]] += adjust
        return charges

    def check_charge(self, charges):
        """
        Check the total charge and the charge difference.

        Parameters
        ----------
        charges : np.ndarray
            Array of atomic charges.

        Returns
        -------
        tuple
            charges : np.ndarray
                Rounded charges.
            total_charge : float
                The total charge.
            charge_difference : float
                The charge difference.
        """
        charges = np.round(charges, self.decimals)
        total_charge = np.sum(charges)
        charge_difference = self.net_charge - total_charge
        self.logger.debug(f"-> Total Charge is {total_charge}")
        self.logger.debug(f"-> Charge difference is {charge_difference}")
        return charges, total_charge, charge_difference
