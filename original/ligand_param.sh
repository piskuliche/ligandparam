#!/bin/bash

thisScript=`basename "$0"`
if [ "$#" -eq 0 ];then
   cat << EOF
------------------------
Needs pdb name and charge when non-zero (i.e. default is 0)
EXAMPLE:
${thisScript} ligand.pdb -1
OR
${thisScript} ligand.pdb
------------------------
Also needs:
- amber/antechamber
- parmutils
- Gaussian16 to be callable by g16
EOF
exit
fi

###########################
# Assigning input arguments
pdb=$1
netcharge=$2
if  [ -z $netcharge ]; then
  netcharge=0
fi
atomtype=$3
if  [ -z $atomtype ]; then
  atomtype=gaff2
fi
if [[ $atomtype == "gaff"* ]]; then
  extraFF="source leaprc.gaff2"
else
  extraFF=""
fi
###########################

base=${pdb%.pdb}
ac_mol2=${base}.antechamber.mol2

echo "######################################################"
echo "STEP1: creating mol2 from antechamber with bcc charges"
echo ""
echo $pdb
echo $ac_mol2
echo $netcharge
echo $atomtype
antechamber -i $pdb -fi pdb -o $ac_mol2 -fo mol2 -c bcc -nc $netcharge -pf y -at $atomtype 

parmutils-mol22g09_ZAP.py $ac_mol2 -t_low HF/6-31G* -t PBE1PBE/6-31G* --resp > ${base}.com

#submit local
cat << "EOF" > submit.slurm
#!/bin/bash
#SBATCH --partition=main
#SBATCH --requeue
#SBATCH --job-name=param
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=10GB
#SBATCH --time=06:00:00
#SBATCH --output=slurm.out

module load gaussian

mem=10
core=8

base=${1%.com}

mkdir gaussianCalcs 2>/dev/null
cd gaussianCalcs
cp -p ../${base}.com .

sed -i "s/MEM=.*GB/MEM=${mem}GB/g" ${base}.com
sed -i "s/NPROC=.*/NPROC=${core}/" ${base}.com
g16 < ${base}.com > ${base}.log
for alp in 0 30 60 90 120 150 180; do
  for bet in 0 30 60 90; do
      parmutils-rotate-g09-resp.py --alp ${alp} --bet ${bet} ${base}.log
  done
done
for alp in 0 30 60 90 120 150 180; do
  for bet in 0 30 60 90; do
       f=$(printf "%s_%.2f_%.2f_%.2f" ${base} ${alp} ${bet} 0)
       g16 <  ${f}.com > ${f}.log
   done
done

EOF

echo "######################################################"
echo "STEP2: running g16 to optimize structure and calc ESP "
echo ""

if ! [[ -f "gaussianCalcs/${base}.log" ]]; then
  echo "Gaussian calculation has been submitted, please wait for it to complete."
  exit
else
  echo "gaussianCalcs/${base}.log exists, skipping Gaussian calculations"
fi

# mol2 from Gaussian file (through antechamber)
# this file will NOT have the correct atomnames, but we want to compare atomtypes so we need it.

echo "######################################################"
echo "STEP3: getting mol2 from gaussian optimized structure"
echo ""
antechamber -fi gout -fo mol2 -i gaussianCalcs/${base}.log -o ${base}.log.mol2 -pf y -at $atomtype

#respFit.py
cat << EOF > respFit.py
#/home/solen/Programs/Amber20/bin/amber.python

import glob
import parmutils
import parmutils.respfit as rf



if __name__ == "__main__":
    comp = parmutils.BASH( 12 )
    model = rf.ResidueResp( comp, 1 )


    model.add_state( "$base", "$base.log.mol2", glob.glob("gaussianCalcs/$base_*.log"), qmmask="@*" )


    model.multimolecule_fit(True)
    model.perform_fit("@*",unique_residues=False)
    #model.preserve_residue_charges_by_shifting()
    model.print_resp()

    #rf.PrintPerturbedCharges(model)

EOF

echo "######################################################"
echo "STEP4: running RESP using parmutils & Tim's scripts"
echo ""

python3 respFit.py > respFit_out.dat


echo "######################################################"
echo "STEP5: updating charges in mol2 according to RESP"
echo ""

# Update charges in mol2
resp_mol2=$base.resp.mol2
cp $ac_mol2 $resp_mol2
sed -i "s/bcc/RESP/" $resp_mol2


chargeCol=72 #charges go after this column
totatom=$(sed '3p;d' $ac_mol2 | awk '{print $1}')
totcharge=0
charge_list=""
truncate -s0 tmp
for atom in $(seq 1 $totatom); do

strt_linnum=$(grep -n "TRIPOS.*ATOM" $ac_mol2 | grep -Eo '^[^:]+')
mol2line=$(( ${strt_linnum} + ${atom} ))
respline=$(( 1 + ${atom} ))
charge_new=$(sed "${respline}q;d" respFit_out.dat | awk '{printf("%7.4f\n",$4)}')
charge_new=$(printf "%9.6f" $charge_new)
echo ${charge_new} >> tmp
sed -i ''"${mol2line}"'s/^\(.\{'${chargeCol}'\}\).\{9\}/\1'"${charge_new}"'/' $resp_mol2
totcharge=$(echo "${totcharge} + ${charge_new}" | bc -l)
done

# If total charge isn't the integer it's supposed to be
if ! [[ ${totcharge} == ${netcharge} ]]; then
  diff=$(echo "${totcharge} - ${netcharge}" | bc -l)

  echo "RESP total charges give ${totcharge}"
  echo "The extra $diff charge will be adjusted"

  count=$(echo "${diff} / 0.0001" | bc -l)
  count=$(printf "%1.f" ${count#-})
  adjust=$(echo "${diff} / ${count}" | bc -l)

  sort tmp > tmp1
  vals=(`sort tmp | tail -n$count`)
  for i in $(seq 0 $(( $count - 1 )) ); do
    charge=${vals[$i]}
    strng=${charge/./\\.}
    strng=${strng/-/\\-}
    chrg_linnum=$(grep -n "$strng" $resp_mol2 | grep -Eo '^[^:]+')
    atom=$(sed "${chrg_linnum}q;d" $resp_mol2 | cut -b 4-7)
    new_charge=$(echo "${charge} - ${adjust}" | bc -l)
    new_charge=$(printf "%9.6f" $new_charge)
    sed -i ''"${chrg_linnum}"'s/^\(.\{'${chargeCol}'\}\).\{9\}/\1'"${new_charge}"'/' $resp_mol2
    echo "charge of atom $atom is adjusted from $charge to $new_charge"
  done

   #check total charges again
   totcharge=0
   for atom in $(seq 1 $totatom); do
   
   strt_linnum=$(grep -n "TRIPOS.*ATOM" $ac_mol2 | grep -Eo '^[^:]+')
   mol2line=$(( ${strt_linnum} + ${atom} ))
   charge=$(sed "${mol2line}q;d" $resp_mol2 | awk '{printf("%7.4f\n",$9)}')
   charge=$(printf "%9.6f" $charge)
   totcharge=$(echo "${totcharge} + ${charge}" | bc -l)
   done
fi

if ! [[ ${totcharge} == ${netcharge} ]]; then
  echo "ERROR something's wrong"
  echo "total charge is: $totcharge"

else
  echo "Total charge in $resp_mol2 is: $totcharge"
  echo "same as provided net charge: $netcharge"
fi

echo "######################################################"
echo "STEP6: checking atomtypes between bcc and gaussian fit"
echo ""

# check atomtype matching
cnt=""
for atom in $(seq 1 $totatom); do

strt_linnum=$(grep -n "TRIPOS.*ATOM" $resp_mol2 | grep -Eo '^[^:]+')
mol2line=$(( ${strt_linnum} + ${atom} ))
type_bcc=$(sed "${mol2line}q;d" $resp_mol2  | cut -b 51-52)
type_gau=$(sed "${mol2line}q;d" ${base}.log.mol2 | cut -b 51-52)

if ! [[ "$type_bcc" == "$type_gau" ]]; then
  sed -i ''"${mol2line}"'s/^\(.\{'50'\}\).\{2\}/\1'"${type_gau}"'/' $resp_mol2
  echo "Atom type of #$atom is changed from bcc $type_bcc to $type_gau based on gaussian optimization"
  cnt=1
fi

done

if [ -z $cnt ]; then
  echo "All atom types match"
fi

# Could add to swap the coordinate information from the gaussian log to the bcc based mol2#

# build off & frcmod 
echo "######################################################"
echo "STEP7: parmchk to create frcmod file"
echo ""

parmchk2 -i $resp_mol2 -f mol2 -o ${base}.frcmod

cat << EOF > leap.saveOffamber.in
source leaprc.protein.ff19SB
source leaprc.RNA.OL3
source leaprc.DNA.OL15
$extraFF
loadamberparams ${base}.frcmod

$base = loadmol2 "$resp_mol2"


saveOff $base ${base}.off
quit
EOF

echo "######################################################"
echo "STEP8: building off file in tleap"
echo ""

tleap -f leap.saveOffamber.in

echo "######################################################"
echo "STEPS ALL DONE"
echo "if everything went well you're files ${base}.off and ${base}.frcmod should be ready"
