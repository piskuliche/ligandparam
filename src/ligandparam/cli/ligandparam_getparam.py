import logging
import sys, shutil
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from ligandparam import __version__
from ligandparam.stages import PDB_Name_Fixer



def set_file_logger(
    logfilename: Path, logname: str = None, filemode: str = "a"
    ) -> logging.Logger:
    """ Set up a file logger for the ligand parameterization process.
    
    Parameters
    ----------
    logfilename : Path
        The path to the log file where the logs will be written.
    logname : str, optional
        The name of the logger. If None, it will be derived from the log filename.
    filemode : str, optional
        The mode in which the log file will be opened. Default is 'a' (append
        mode). Use 'w' for write mode to overwrite the log file.
    Returns
    -------
    logger : logging.Logger
        A configured logger instance that writes logs to the specified file.
    
    """
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

def worker(recipe_name: str, mol: str, resname: str, cwd: Path, net_charge: float, atom_type: str = "gaff2", charge_model: str = "bcc", model: str = None, sqm: str = True, data_cwd: str = "param", nprocs: int = 1, mem: int = 1, reference_pdb: str = None) -> Path:
    binder_dir = cwd / data_cwd / resname
    binder_dir.mkdir(parents=True, exist_ok=True)
    binder_pdb = cwd / mol
    logger = set_file_logger(
        binder_dir / f"{resname}.log", filemode="w"
    )
    """ Worker function to execute the ligand parameterization recipe.
    
    Parameters
    ----------
    recipe_name : str
        The name of the recipe to be used for ligand parameterization.
    mol : str
        The input PDB file containing the ligand.
    resname : str
        The residue name for the ligand.
    cwd : Path
        The current working directory where the output files will be stored.
    net_charge : float
        The net charge of the ligand.
    atom_type : str, optional
        The atom type for the ligand (default is "gaff2").
    charge_model : str, optional
        The charge model for the ligand (default is "bcc"). Options are "bcc" or
        "abcg2".
    model : str, optional
        The path to the DeepMD model file (optional).
    sqm : bool, optional
        Whether to use SQM calculations for geometry optimization (default is True).
    data_cwd : str, optional
        The directory to store output files (default is "param").
    nprocs : int, optional
        The number of processes to use for parallel execution (default is 1).
    mem : int, optional
        The amount of memory in GB to allocate for the process (default is 1GB).
    
    Returns
    -------
    None
    """

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

    
    
    if reference_pdb is not None:
        logger.info(f"Reference PDB file: {reference_pdb}")
        fix_pdb_stage = PDB_Name_Fixer(f"build_{resname}", binder_pdb, binder_dir, out_pdb=f"{binder_pdb.parent}/fix_{binder_pdb.name}", reference_pdb=reference_pdb, align=True, logger=logger)
        fix_pdb_stage.execute(dry_run=False)
        logger.info("PDB name fixing complete.")
        out_pdb = f"{binder_pdb.parent}/fix_{binder_pdb.name}"
    else:
        out_pdb = binder_pdb


    recipe = recipe_selector(
        recipe_name,
        in_filename = f"{out_pdb}",
        cwd        = binder_dir,
        atom_type  = atom_type,
        charge_model = charge_model,
        net_charge = net_charge,
        logger     = logger,
        molname    = resname,
        model      = model,
        sqm        = sqm,
        nproc      = nprocs,
        mem        = mem
    )
    logger.info(f"Recipe selected: {recipe_name}")
    recipe.setup()
    recipe.execute()
    logger.info("Recipe execution complete.")

def recipe_selector(recipe_name: str, **kwargs):
    """ Selects and returns the appropriate recipe class based on the recipe name.
    
    Parameters
    ----------
    recipe_name : str
        The name of the recipe to be used for ligand parameterization.
    **kwargs : dict
        Additional keyword arguments to be passed to the recipe class constructor.
    
    Returns
    -------
    AbstractStage
        An instance of the selected recipe class.
    
    Raises
    ------
    ValueError
        If the recipe name is not recognized, a ValueError is raised with a message
        listing the available recipes.
    
    """
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
    elif recipe_name == "sqmligand":
        from ligandparam.recipes import SQMLigand
        return SQMLigand(**kwargs)
    else:
        raise ValueError(f"Unknown recipe name: {recipe_name}. Available recipes: lazyligand, lazierligand, freeligand, dplazyligand, dpfreeligand.")


def main():
    """ Main function to parse command line arguments and execute the ligand parameterization worker.
    
    This function uses argparse to handle command line arguments and calls the worker
    function with the parsed arguments. It sets up the current working directory and
    initializes the logger for the ligand parameterization process.
    
    Raises
    ------
    SystemExit
        If the command line arguments are not provided correctly, argparse will raise
        a SystemExit exception, which will terminate the program with an error message.
    
    """
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
    parser.add_argument("-n", "--nproc", type=int, default=1, help="Number of processes to use (default: 1)")
    parser.add_argument("-mem", "--mem", type=int, default=1, help="Memory in GB to allocate for the process (default: 1GB)")
    parser.add_argument("-ref", "--reference_pdb", type=str, default=None, help="Reference PDB file for name fixing (optional)")

    args = parser.parse_args()

    cwd = Path.cwd()

    worker(
        recipe_name=args.recipe_name,
        mol=args.input,
        cwd=cwd,
        resname=args.resname,
        data_cwd=args.data_cwd,
        net_charge=args.net_charge,
        atom_type=args.atom_type,
        charge_model=args.charge_model,
        model=args.model,
        sqm=args.sqm,
        nprocs=args.nproc,
        mem=args.mem,
        reference_pdb=args.reference_pdb
    )

if __name__ == "__main__":
    main()