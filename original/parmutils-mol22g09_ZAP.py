#!/usr/bin/env python3

import sys
import argparse
import parmed
import parmutils
from io import StringIO as StringIO
import subprocess

if __name__ == "__main__":

    parser = argparse.ArgumentParser \
    ( formatter_class=argparse.RawDescriptionHelpFormatter,
      description="""Writes a Gaussian input file, given a mol2 file""" )
    
    parser.add_argument \
        ("-n","--nproc",
         help="number of cores (default: 6)",
         type=int,
         default=6,
         required=False )

    parser.add_argument \
        ("-t_low","--theory_low",
         help="model chemistry (default: 'None')",
         type=str,
         default=None,
         required=False )

    parser.add_argument \
        ("-t","--theory",
         help="model chemistry (default: 'None')",
         type=str,
         default="PBE1PBE/6-31G*",
         required=False )

    parser.add_argument \
        ("--sp",
         help="if present, then do a single-point, not a geometry optimization",
         action='store_true',
         required=False )

    parser.add_argument \
        ("--tight",
         help="if present, then set tight geometry optimization tolerance",
         action='store_true',
         required=False )

    parser.add_argument \
        ("--nohess",
         help="if present, then skip the CalcFC option within the geometry optimization",
         action='store_true',
         required=False )

    parser.add_argument \
        ("--resp",
         help="if present, then perform a RESP evaluation after the geometry optimization (or as part of the single-point calculation)",
         action='store_true',
         required=False )
    
    parser.add_argument \
        ("--run",
         help="if present, then run the job",
         action='store_true',
         required=False )

    parser.add_argument \
        ('mol2',
         metavar='filename',
         type=str,
         nargs=1,
         help='mol2 file')

    args = parser.parse_args()
    

    p = parmutils.OpenParm( args.mol2[0], xyz=None )
    q = int(round(sum( [ a.charge for a in p.atoms ] )))
    z = sum( [ a.atomic_number for a in p.atoms ] )
    nelec = z-q
    m = 1
    if nelec % 2 != 0:
       m = 2



    if args.sp:
        opts = "SP"
        if args.resp:
            opts += " SCF(Conver=6) NoSymm Test\n   Pop=mk IOp(6/33=2) GFInput GFPrint"
    else:
        opts = "OPT"
        o=[]
        if not args.nohess:
           o.append("CalcFC")
        if args.tight:
           o.append("Tight")
        if len(o) > 0:
           opts += "(%s)"%( ",".join(o) )

           
    base = args.mol2[0].replace(".mol2","")
    com = base+".com"
    log = base+".log"
    chk = base+".chk"
       
    fh = StringIO()

    fh.write("%%NPROC=%i\n"%(args.nproc))
    fh.write("%MEM=2GB\n")
    fh.write("%%CHK=%s\n"%(chk))
    if args.theory_low is not None:# Writes the low level theory in place of the high level if provided.
        fh.write("#P %s %s\n\n"%(args.theory_low, opts))
    else:
        fh.write("#P %s %s\n\n"%(args.theory,opts))
    fh.write("%s\n\n%i %i\n"%(args.mol2[0],q,m))
    for a in p.atoms:
        fh.write("%4s %12.6f %12.6f %12.6f\n"%(a.element,a.xx,a.xy,a.xz))
    fh.write("\n\n")

    # Write the high level, if low level is provided.
    if args.theory_low is not None:
        fh.write("--Link1--\n")
        fh.write("%%NPROC=%i\n"%(args.nproc))
        fh.write("%MEM=2GB\n")
        fh.write("%%CHK=%s\n"%(chk))
        fh.write("#P %s %s GEOM(AllCheck) Guess(Read)\n\n"%(args.theory,opts))
        fh.write("\n\n")


    if not args.sp and args.resp:
        fh.write("--Link1--\n")
        fh.write("%%NPROC=%i\n"%(args.nproc))
        fh.write("%MEM=2GB\n")
        fh.write("%%CHK=%s\n"%(chk))
        fh.write("#P HF/6-31G* Geom(AllCheck) Guess(Read) NoSymm Test\n   Pop=mk IOp(6/33=2) GFInput GFPrint\n")
        fh.write("\n\n")

            
    s = fh.getvalue()
    fh.close()

    if not args.run:
       print(s)
    else:

       fh = open(com,"w")
       fh.write(s)
       fh.close()
       print("Running: g09 < %s > %s"%(com,log))
       print("You can track the progress by opening a window and run:")
       print("tail -n 10000 -f %s | grep -E 'Maximum Force|RMS     Force|Maximum Displacement|RMS     Displacement|Step number|SCF Done'"%(log))
       subprocess.call("g09 < %s > %s"%(com,log),shell=True)


