"""CLI entry point for batch ligand parameterization."""

from __future__ import annotations

from pathlib import Path

from ligandparam.Log import set_file_logger
from ligandparam.stages import StagePdbNameFixer


def worker(
    recipe_name: str,
    mol: str,
    resname: str,
    cwd: Path,
    net_charge: float,
    atom_type: str = "gaff2",
    charge_model: str = "bcc",
    model: str = None,
    sqm: str = True,
    data_cwd: str = "param",
    nprocs: int = 1,
    mem: int = 1,
    reference_pdb: str = None,
    force_gaussian_rerun: bool = False,
) -> Path:
    """Execute a ligand parameterization recipe for one ligand.

        After completion, ALPS / ``lig-dihed-correct`` can apply optional
        torsion corrections to the generated mol2/lib/frcmod.

    Returns
    -------
    Path
        Output directory for this ligand.
    """
    binder_dir = cwd / data_cwd / resname
    binder_dir.mkdir(parents=True, exist_ok=True)
    binder_pdb = cwd / mol
    logger = set_file_logger(
        binder_dir / f"{resname}.log",
        logname=resname,
        filemode="w",
        include_version=True,
    )

    logger.info("Working on ligand: %s", resname)
    if not binder_pdb.is_file():
        raise FileNotFoundError(
            f"Input file {binder_pdb} does not exist. Please provide a valid PDB file."
        )
    if not binder_dir.is_dir():
        raise NotADirectoryError(
            f"Output directory {binder_dir} does not exist. Please provide a valid directory."
        )

    logger.info(
        f"Starting ligand parameterization for {resname} using recipe '{recipe_name}'"
    )
    logger.info(f"Input file: {binder_pdb}")
    logger.info(f"Output directory: {binder_dir}")
    logger.info(f"Net charge: {net_charge}")
    logger.info(f"Atom type: {atom_type}")
    logger.info(f"Charge model: {charge_model}")
    logger.info(f"force_gaussian_rerun (-O): {force_gaussian_rerun}")
    if model is not None:
        logger.info(f"Using DeepMD model: {model}")
    if sqm:
        logger.info("Using SQM calculations for geometry optimization.")
    else:
        logger.info("Not using SQM calculations for geometry optimization.")
    logger.info("Starting recipe execution...")

    if reference_pdb is not None:
        logger.info(f"Reference PDB file: {reference_pdb}")
        fix_pdb_stage = StagePdbNameFixer(
            f"build_{resname}",
            binder_pdb,
            binder_dir,
            out_pdb=f"{binder_pdb.parent}/fix_{binder_pdb.name}",
            reference_pdb=reference_pdb,
            align=True,
            logger=logger,
        )
        fix_pdb_stage.execute(dry_run=False)
        logger.info("PDB name fixing complete.")
        out_pdb = f"{binder_pdb.parent}/fix_{binder_pdb.name}"
    else:
        out_pdb = binder_pdb

    recipe = recipe_selector(
        recipe_name,
        in_filename=f"{out_pdb}",
        cwd=binder_dir,
        atom_type=atom_type,
        charge_model=charge_model,
        net_charge=net_charge,
        logger=logger,
        molname=resname,
        model=model,
        sqm=sqm,
        nproc=nprocs,
        mem=mem,
        force_gaussian_rerun=force_gaussian_rerun,
    )
    logger.info(f"Recipe selected: {recipe_name}")
    recipe.setup()
    recipe.execute()
    logger.info("Recipe execution complete.")
    return binder_dir


def recipe_selector(recipe_name: str, **kwargs):
    """Select and return the appropriate recipe class based on the recipe name."""
    from ligandparam.recipes.Registry import available_recipes, get_recipe

    try:
        return get_recipe(recipe_name, **kwargs)
    except ValueError as exc:
        raise ValueError(
            f"Unknown recipe name: {recipe_name}. Available recipes: "
            + ", ".join(available_recipes())
        ) from exc


def main():
    """Parse command line arguments and execute the ligand parameterization worker."""
    import argparse

    parser = argparse.ArgumentParser(description="Ligand parameterization CLI")
    parser.add_argument("-i", "--input", type=str, required=True, help="Input PDB file with ligand")
    parser.add_argument("-r", "--resname", type=str, required=True, help="Residue name for the ligand")
    parser.add_argument("-d", "--data_cwd", type=Path, required=True, help="Directory to store output files")
    parser.add_argument("-a", "--atom_type", type=str, default="gaff2", help="Atom type for the ligand (default: gaff2)")
    parser.add_argument(
        "-cm",
        "--charge_model",
        type=str,
        default="bcc",
        choices=["bcc", "abcg2"],
        help="Charge model for the ligand (default: bcc, options: bcc, abcg2)",
    )
    parser.add_argument("-c", "--net_charge", type=float, default=0.0, help="Net charge of the ligand")
    parser.add_argument("-m", "--model", type=str, default=None, help="DeepMD model file path (optional)")
    parser.add_argument("--sqm", action="store_true", help="Use SQM calculations")
    parser.add_argument("-rn", "--recipe_name", type=str, required=True, help="Recipe name for the ligand processing")
    parser.add_argument("-n", "--nproc", type=int, default=1, help="Number of processes to use (default: 1)")
    parser.add_argument(
        "-mem",
        "--mem",
        type=int,
        default=1,
        help=(
            "Node memory budget in GB, split across concurrent orientation "
            "ESP jobs (default: 1). Do not set this to the per-job Gaussian "
            "request; n_workers * %%MEM <= this value."
        ),
    )
    parser.add_argument("-ref", "--reference_pdb", type=str, default=None, help="Reference PDB file for name fixing (optional)")
    parser.add_argument(
        "-O",
        "--force-gaussian-rerun",
        action="store_true",
        help=(
            "Force re-run of Gaussian stages even when logs already show "
            "Normal termination. By default, complete jobs are skipped and "
            "only incomplete orientation ESP jobs are re-run."
        ),
    )

    args = parser.parse_args()

    worker(
        recipe_name=args.recipe_name,
        mol=args.input,
        cwd=Path.cwd(),
        resname=args.resname,
        data_cwd=args.data_cwd,
        net_charge=args.net_charge,
        atom_type=args.atom_type,
        charge_model=args.charge_model,
        model=args.model,
        sqm=args.sqm,
        nprocs=args.nproc,
        mem=args.mem,
        reference_pdb=args.reference_pdb,
        force_gaussian_rerun=args.force_gaussian_rerun,
    )


if __name__ == "__main__":
    main()
