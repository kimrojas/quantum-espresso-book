#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -q xe1.q
#$ -pe x8 8
#$ -j y
#$ -N calc

module load qe/7.0
module load ase
source activate espresso

export I_MPI_PIN=1
export OMP_NUM_THREADS=1
export I_MPI_FABRICS=shm:ofi

cat $PE_HOSTFILE | awk '{ print $1":"$2/ENVIRON["OMP_NUM_THREADS"] }' > hostfile.$JOB_ID

# ----------------- MAIN PROCESS -----------------------------------------------

echo -e "========= Job started  at `date` ==========" >> run_calc.log
   # Export the environment
export ESPRESSO_PSEUDO='/home/krojas/QE/pseudo'
export ESPRESSO_TMPDIR='./TMPDIR'
   # Run SCF calculation
mpirun pw.x -in espresso_scf.pwi > espresso_scf.pwo
   # Run NSCF(BANDS) calculation
mpirun pw.x -in espresso_bands.pwi > espresso_bands.pwo
echo -e "========= Job finished at `date` ==========" >> run_calc.log


echo -e "\nCALCULATING BAND PROPERTIES" >> run_calc.log
mpirun bands.x -nk 1 -in espresso_bandx.inp > espresso_bandx.out
echo -e "DONE" >> run_calc.log

echo -e "\nCREATING BAND STRUCTURE PLOT"
python bandplot.py
echo -e "PLOT READ --- look for dosplot.png" >> run_calc.log

# ----------------- MAIN PROCESS -----------------------------------------------

rm -f hostfile.$JOB_ID
