#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -q xs2.q
#$ -pe x16 16
#$ -j y
#$ -N dftb_au5ni_multirun
#$ -t 1-3

# --- LOAD MODULES ---
module load cmake/3.18.3
module load intel/2020.2.254
module load intelmpi/2020.2.254
module load python/3.8

# --- ACTIVATE PYTHON/CONDA ENVIRONMENT ---
source activate tutorial317

# --- ACTIVATE DFTB+ / GOFEE ---
scripts_dir=~/tutorial_files/test-suite/scripts
source ${scripts_dir}/activate-dftb.sh
source ${scripts_dir}/activate-gofee.sh

# --- SET NODE ENVIRONMENT ---
export I_MPI_PIN=1
export OMP_NUM_THREADS=1
export I_MPI_FABRICS=shm:ofi
cat $PE_HOSTFILE | awk '{ print $1":"$2/ENVIRON["OMP_NUM_THREADS"] }' > hostfile.$JOB_ID


# --- START CALCULATION ---
echo "========= Job started  at `date` =========="

PYTHONFILE='input_au5ni_multirun.py'
WORKDIR="RUNS/run_${SGE_TASK_ID}"

mkdir -p $WORKDIR
cp $PYTHONFILE $WORKDIR
cd $WORKDIR

mpirun -np $NSLOTS python ${PYTHONFILE} 

echo "========= Job finished at `date` =========="

rm -f hostfile.$JOB_ID