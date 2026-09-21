"""CLI to generate AMBER parameters from a MOL2 via OpenFF Sage."""

import os

from ligandparam.stages import StageSageCreate, StageSageToAmber

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Generate AMBER parameters from a mol2 file.")
    parser.add_argument("input_mol2", type=str, help="Path to the input mol2 file.")
    parser.add_argument("output_tag", type=str, help="Tag for the output files.")
    parser.add_argument(
        "--use-mol2-charges",
        action="store_true",
        help="Use atom charges from the input MOL2 instead of OpenFF charge assignment.",
    )
    args = parser.parse_args()
    out_parm = f"{args.output_tag}.parm7"
    out_rst7 = f"{args.output_tag}.rst7"
    out_lib = f"{args.output_tag}.lib"
    out_frcmod = f"{args.output_tag}.frcmod"

    stage_create = StageSageCreate(
        "sage_creation",
        args.input_mol2,
        os.getcwd(),
        out_parm=out_parm,
        use_mol2_charges=args.use_mol2_charges,
    )
    stage_create.execute()

    stage_export = StageSageToAmber(
        "sage_to_amber",
        out_parm,
        os.getcwd(),
        in_rst7=out_rst7,
        out_lib=out_lib,
        out_frcmod=out_frcmod,
        resname=stage_create.resname,
    )
    stage_export.execute()

if __name__ == "__main__":
    main()
