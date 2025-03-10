import sys
from pathlib import Path
import logging

from ligandparam.stages import SDFToPDBBatch

cwd = Path(sys.argv[0]).resolve().parent
sdf = cwd / "ligands.sdf"
bad_sdf = cwd / "bad_ligands.sdf"

# Send output to stdout, though SDFToPDBBatch won't print unless there's an error
logger = logging.getLogger("mylog")
logger.setLevel(logging.INFO)
stream_handler = logging.StreamHandler(sys.stdout)
stream_handler.setLevel(logging.INFO)
logger.addHandler(stream_handler)

# Quickest way to convert an SDF file to PDB files. Filenames and resnames will be the same as the SDF molecule names
sdf_to_pdb = SDFToPDBBatch("sdf2pdb", sdf, cwd, logger=logger)
sdf_to_pdb.execute()

# This will trigger an error, as there are multiple molecules with the same name
sdf_to_pdb = SDFToPDBBatch("sdf2pdb", bad_sdf, cwd, logger=logger)
sdf_to_pdb.execute()

"""
There're many alternatives when calling SDFToPDBBatch:

# Set the same resname for all the molecules, but the filenames will be extracted from the SDF molecule name
sdf_to_pdb = SDFToPDBBatch("sdf2pdb", sdf, cwd, resname="LIG")

# Set the same resname for all the molecules, and the filenames to ligand_0.pdb, ligand_1.pdb, etc.
sdf_to_pdb = SDFToPDBBatch("sdf2pdb", sdf, cwd, resname="LIG",
                           out_pdb_template=Path("/home/pb777/ligandparam/examples/05_FromSDF/ligand.pdb"))
          
# molecules will keep their SDF resnames, but we set their filenames                 
sdf_to_pdb = SDFToPDBBatch("sdf2pdb", sdf, cwd, out_pdb_template=Path("/home/pb777/ligandparam/examples/05_FromSDF/ligand.pdb"))

# we set resnames and filenames
resnames = ("17", "1h1q", "1h1r", "1h1s", "1oi9", "1oiu", "1oiy", "20", "21", "22", "26", "28", "29", "30", "31", "32")
out_pdbs = [cwd / f"{pdb}.pdb" for pdb in resnames]
sdf_to_pdb = SDFToPDBBatch("sdf2pdb", sdf, cwd, out_pdbs=out_pdbs, resnames=resnames)
"""
