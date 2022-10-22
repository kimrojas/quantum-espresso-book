# C6H6 system - parallelized

A parallelized structure search calculation for a C6H6 system using GOFEE. Files are available in the `test-suite` provided.

The degree of parallelization is based on `MPI` and can be controlled using the 

```bash
mpirun -np $NSLOTS python input_c6h6_parallel.py
``` 

command. The input file is the same as in the previous case, but we change the jobscript.





**PARALLELIZED JOBSCRIPT**
```bash
#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -q xs2.q
#$ -pe x16 16
#$ -j y
#$ -N dftb_c6h6_parallel

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

mpirun -np $NSLOTS python input_c6h6_parallel.py 

echo "========= Job finished at `date` =========="

rm -f hostfile.$JOB_ID

```