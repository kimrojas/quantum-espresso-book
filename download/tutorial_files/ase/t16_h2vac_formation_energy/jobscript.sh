#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -q xb1.q
#$ -pe x32 32
#$ -j y
#$ -N calc
#$ -t 1-10

module load qe/7.0
module load ase
source activate espresso

export I_MPI_PIN=1
export OMP_NUM_THREADS=1
export I_MPI_FABRICS=shm:ofi

cat $PE_HOSTFILE | awk '{ print $1":"$2/ENVIRON["OMP_NUM_THREADS"] }' > hostfile.$JOB_ID

echo "========= Job started  at `date` =========="

time python run_calc.py -i $SGE_TASK_ID

echo "========= Job finished at `date` =========="

rm -f hostfile.$JOB_ID
