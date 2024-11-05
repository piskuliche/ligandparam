from ligandparam.stages import StageMapTarget
from ligandparam.io.coordinates import Remove_PDB_CONECT

reference_structure = "1Y27.pdb"


baseoptions = {
    "base_name": "F3G",
    "nproc": 12,
    "mem": "60GB",
    "net_charge": 0,
    "atom_type": "gaff2",
    "target_pdb": reference_structure,
    "force_gaussian_rerun": False
}

Remove_PDB_CONECT(reference_structure)

newstage = StageMapTarget("Map", inputoptions=baseoptions)

newstage.execute()
