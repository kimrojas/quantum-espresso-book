#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -q xh1.q
#$ -pe x24 24
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
echo "========= Job started  at `date` ==========" >> run_calc.log
   # Export the environment
export ESPRESSO_PSEUDO='/home/krojas/QE/pseudo'
export ESPRESSO_TMPDIR='./TMPDIR'
   # Run calculation
mpirun pw.x -in espresso.pwi > espresso.pwo
echo "========= Job finished at `date` ==========" >> run_calc.log
   # Quick analysis: Get BFGS steps
echo -e "\nBFGS Steps:" >> run_calc.log
steps=$(grep 'bfgs steps' espresso.pwo)
echo -e "$steps" >> run_calc.log
   # Quick analysis: Get Final energy
echo -e "\nFinal energy:" >> run_calc.log
energy=$(grep Final espresso.pwo)
echo -e "$energy" >> run_calc.log
   # Quick analysis: Get Final structure
echo -e "\nFinal structure:" >> run_calc.log
structure=$(awk '/Begin final coordinates/{flag=1; next} /End final coordinates/{flag=0} flag' espresso.pwo)
echo -e "$structure" >> run_calc.log
# ----------------- MAIN PROCESS -----------------------------------------------

rm -f hostfile.$JOB_ID
