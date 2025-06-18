import sys
from pathlib import Path
import argparse
import logging

from ligandparam.stages import LigHFix


def get_opts() -> dict:
    parser = argparse.ArgumentParser(description="Fix hydrogenation and bonding in ligands that come from PDB")
    parser.add_argument("-i", "--ligand_id", help="Path to the input PDB file.", required=True)
    parser.add_argument("-p", "--in_pdb", help="Path to the input PDB file.", required=True)
    parser.add_argument("-o", "--out_pdb", help="Path to the output file.", default="output.txt", required=True)

    args = parser.parse_args()
    opts = {}
    opts["ligand_id"] = args.ligand_id
    opts["in_pdb"] = Path(args.in_pdb)
    opts["out_pdb"] = Path(args.out_pdb)

    if not opts["in_pdb"].is_file() or opts["in_pdb"].suffix != ".pdb":
        raise ValueError(f"Bad input PDB: {opts['in_pdb']}")

    if not opts["out_pdb"].parent.is_dir():
        raise ValueError(f"Bad dir for out_pdb: {opts['out_pdb'].parent}")

    if opts["out_pdb"].suffix != ".pdb":
        raise ValueError(f"Bad output PDB: {opts['out_pdb']}")

    return opts


def ligfix():
    opts = get_opts()

    # Send output to stdout, though it probably won't print anything unless there's an error
    logger = logging.getLogger("mylog")
    logger.setLevel(logging.INFO)
    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setLevel(logging.INFO)
    logger.addHandler(stream_handler)

    s = LigHFix("liga", main_input=opts["ligand_id"], cwd=opts["out_pdb"].parent, in_pdb=opts["in_pdb"],
                out_pdb=opts["out_pdb"], logger=logger)

    s.execute()


if __name__ == "__main__":
    ligfix()
