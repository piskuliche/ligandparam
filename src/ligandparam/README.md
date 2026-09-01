# ligandparam

Stage-based parameterization of nonstandard ligands and residues for Amber.
Recipes assemble a pipeline (Antechamber → Gaussian ESP / RESP → `parmchk2`
→ LEaP). Each stage is a single tool step you can add, remove, or replace.

Documentation: https://ligandparam.readthedocs.io/en/latest/

## Install

AmberTools (`antechamber`, `parmchk2`, `tleap`) must be on `PATH`. Gaussian
recipes also need `g16` (or the binary you pass as `gaussian_binary`).

From this directory:

```bash
python3 -m pip install -e .
```

From the ALPS / RutgersLBSR monorepo root:

```bash
python3 -m pip install -e src/ligandparam
```

Optional extras:

| Extra | What it enables |
|---|---|
| `ml` | DeepMD recipes (`DPLigand`, `DPFreeLigand`) |
| `sage` | `lig-to-sage` (OpenFF Sage → Amber) |
| `misc` | Open Babel |
| `docs` | Sphinx |

```bash
python3 -m pip install -e ".[sage]"
```

`numpy<2` matches ParmEd. TensorFlow / DeepMD on HPC are usually easier from
conda-forge than from pip.

## Recipes

```python
from ligandparam.recipes import LazyLigand, FreeLigand

recipe = LazyLigand(
    "ligand.pdb",
    "work",
    net_charge=0,
    nproc=8,
    mem=16,
    logger="stream",
)
recipe.setup()
recipe.execute()
```

| Recipe | CLI name | Typical use |
|---|---|---|
| `LazyLigand` | `lazyligand` | Gaussian min + single-orientation RESP |
| `FreeLigand` | `freeligand` | Multi-orientation RESP (`so3_n28`, 28 ESP jobs) |
| `LazierLigand` | `lazierligand` | BCC / Antechamber charges, no Gaussian |
| `SQMLigand` | `sqmligand` | SQM geometry, BCC charges |
| `DPLigand` | `dplazyligand` | DeepMD minimize, then Lazy-style RESP |
| `DPFreeLigand` | `dpfreeligand` | DeepMD minimize, then Free-style RESP |

Outputs (working directory): `{label}.mol2`, `{label}.lib`, `{label}.frcmod`.

`FreeLigand` / `DPFreeLigand` default to a fixed 28-point quaternion pack
(`orientation_protocol="so3_n28"`). Pass `"legacy_euler"` for the historical
alpha/beta grid. Same job count either way.

Gaussian jobs that already show `Normal termination` are skipped on resume.
Pass `force_gaussian_rerun=True` (CLI `-O`) to re-run them.

## Command-line tools

```bash
lig-getparam -i ligand.pdb -r LIG -d ./param -rn lazyligand -c 0 -n 8 --mem 16
lighfix -i LIG -p ligand.pdb -o ligand.fixed.pdb
smiles-to-pdb -s "c1ccccc1O" -o phenol.pdb -rn LIG
lig-to-sage input.mol2 out_tag
```

`lig-getparam --recipe_name` accepts the CLI names in the table above.
`--mem` is the **node** memory budget in GB (split across concurrent
orientation ESP jobs), not the per-job Gaussian `%MEM`.

## Layout

```text
ligandparam/
  recipes/     LazyLigand, FreeLigand, ...
  stages/      Initialize, Gaussian, RESP, Leap, ...
  io/          PDB / mol2 / Gaussian / SMILES / orientations
  multiresp/   multi-orientation RESP helpers
  runtime/     Gaussian orientation job board, CPU split
  cli/         lig-getparam, lighfix, smiles-to-pdb, lig-to-sage
  Log.py       stream / file logging (GitLab main formats)
```

## Logging

`ligandparam.Log` is self-contained (no ALPS, no `runtime.Console`).

| Setup | Where | Line shape |
|---|---|---|
| `logger="stream"` | stdout | message only |
| `logger="file"` | `{label}.log` | `YYYY-mm-dd HH:MM:SS - LEVEL - message` |
| `lig-getparam` | `{resname}.log` | same, with the package version after LEVEL |

Pass `logger="stream"` or `logger="file"` on a recipe, or a `logging.Logger`.
The default is a null handler. Configuring a stream or file logger prints the
ligandparam ASCII banner and a random quote from `pkgdata/quotes.txt` (once
per process). Recipes log another quote after a successful run. ALPS can
parse these files later.

This package does **not** run torsion fitting or ligand fragmentation.
Those stay in ffpopt / scission (wired by ALPS on a later push).
Recipes still accept `dihed_correct=...` so old callers do not break;
ligandparam only records the flags and warns.
