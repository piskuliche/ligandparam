import logging
import sys, shutil
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from ligandparam import __version__



def set_file_logger(
    logfilename: Path, logname: str = None, filemode: str = "a"
    ) -> logging.Logger:
    if logname is None:
        logname = Path(logfilename).stem
    logger = logging.getLogger(logname)
    logger.setLevel(logging.INFO)
    formatter = logging.Formatter(
        "{asctime} - {levelname} - {version} {message}",
        style="{",
        datefmt="%Y-%m-%d %H:%M:%S",
        defaults={"version": __version__},
    )
    file_handler = logging.FileHandler(filename=logfilename, mode=filemode)
    file_handler.setLevel(logging.INFO)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    return logger

def worker(recipe_name: str, mol: str, resname: str, cwd: Path, net_charge: float, atom_type: str = "gaff2", charge_model: str = "bcc", model: str = None, sqm: str = True, data_cwd: str = "param") -> Path:
    binder_dir = cwd / data_cwd / resname
    binder_dir.mkdir(parents=True, exist_ok=True)
    binder_pdb = cwd / mol
    logger = set_file_logger(
        binder_dir / f"{resname}.log", filemode="w"
    )

    print("Working on ligand:", resname)
    if not binder_pdb.is_file():
        raise FileNotFoundError(f"Input file {binder_pdb} does not exist. Please provide a valid PDB file.")
    if not binder_dir.is_dir():
        raise NotADirectoryError(f"Output directory {binder_dir} does not exist. Please provide a valid directory.")
    
    logger.info(f"Starting ligand parameterization for {resname} using recipe '{recipe_name}'")
    logger.info(f"Input file: {binder_pdb}")
    logger.info(f"Output directory: {binder_dir}")
    logger.info(f"Net charge: {net_charge}")
    logger.info(f"Atom type: {atom_type}")
    logger.info(f"Charge model: {charge_model}")
    if model is not None:
        logger.info(f"Using DeepMD model: {model}")
    if sqm:
        logger.info("Using SQM calculations for geometry optimization.") 
    else:
        logger.info("Not using SQM calculations for geometry optimization.")
    logger.info("Starting recipe execution...")
    

    recipe = recipe_selector(
        recipe_name,
        in_filename = f"{binder_pdb}",
        cwd        = binder_dir,
        atom_type  = atom_type,
        charge_model = charge_model,
        net_charge = net_charge,
        logger     = logger,
        molname    = resname,
        model      = model,
        sqm        = sqm,
    )
    logger.info(f"Recipe selected: {recipe_name}")
    recipe.setup()
    logger.info("Recipe setup complete. Executing recipe...")
    recipe.execute()
    logger.info("Recipe execution complete.")

def recipe_selector(recipe_name: str, **kwargs):
    if recipe_name == "lazyligand":
        from ligandparam.recipes.lazyligand import LazyLigand
        return LazyLigand(**kwargs)
    elif recipe_name == "lazierligand":
        from ligandparam.recipes.lazierligand import LazierLigand
        return LazierLigand(**kwargs)
    elif recipe_name == "freeligand":
        from ligandparam.recipes.freeligand import FreeLigand
        return FreeLigand(**kwargs)
    elif recipe_name == "dplazyligand":
        from ligandparam.recipes.dplazyligand import DPLigand
        return DPLigand(**kwargs)
    elif recipe_name == "dpfreeligand":
        from ligandparam.recipes.dpfreeligand import DPFreeLigand
        return DPFreeLigand(**kwargs)
    else:
        raise ValueError(f"Unknown recipe name: {recipe_name}. Available recipes: lazyligand, lazierligand, freeligand, dplazyligand, dpfreeligand.")


def main():
    import argparse


    parser = argparse.ArgumentParser(description="Ligand parameterization CLI")
    parser.add_argument("-i", "--input", type=str, required=True, help="Input PDB file with ligand")
    parser.add_argument("-r", "--resname", type=str, required=True, help="Residue name for the ligand")
    parser.add_argument("-d", "--data_cwd", type=Path, required=True, help="Directory to store output files")
    parser.add_argument("-a", "--atom_type", type=str, default="gaff2", help="Atom type for the ligand (default: gaff2)")
    parser.add_argument("-cm", "--charge_model", type=str, default="bcc", choices=["bcc", "abcg2"], help="Charge model for the ligand (default: bcc, options: bcc, abcg2)")
    parser.add_argument("-c", "--net_charge", type=float, default=0.0, help="Net charge of the ligand")
    parser.add_argument("-m", "--model", type=str, default=None, help="DeepMD model file path (optional)")
    parser.add_argument("--sqm", action='store_true', help="Use SQM calculations")
    parser.add_argument("-rn", "--recipe_name", type=str, required=True, help="Recipe name for the ligand processing")

    args = parser.parse_args()

    cwd = Path.cwd()

    worker(
        recipe_name=args.recipe_name,
        mol=args.input,
        resname=args.resname,
        data_cwd=args.data_cwd,
        net_charge=args.net_charge,
        atom_type=args.atom_type,
        charge_model=args.charge_model,
        model=args.model,
        sqm=args.sqm,
        cwd=cwd
    )

if __name__ == "__main__":
    main()