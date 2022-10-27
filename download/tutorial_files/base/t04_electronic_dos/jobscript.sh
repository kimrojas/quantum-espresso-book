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

echo -e "========= Job started  at `date` ==========" >> run_calc.log
   # Export the environment
export ESPRESSO_PSEUDO='/home/krojas/QE/pseudo'
export ESPRESSO_TMPDIR='./TMPDIR'
   # Run SCF calculation
mpirun pw.x -in espresso_scf.pwi > espresso_scf.pwo
   # Run NSCF calculation
mpirun pw.x -in espresso_nscf.pwi > espresso_nscf.pwo
echo -e "========= Job finished at `date` ==========" >> run_calc.log

echo -e "\nCALCULATING DOS: Wavefunction projections" >> run_calc.log
mpirun -n 1 projwfc.x -in projwfc.in > projwfc.out
echo -e "DONE" >> run_calc.log

echo -e "\nLOWDIN CHARGE SUMMARY: " >> run_calc.log
charges=$(awk '/Lowdin Charges/{flag=1; next} /PROJWFC/{flag=0} flag' projwfc.out)
echo -e "$charges" >> run_calc.log

echo -e "\nCREATING DOS PLOT" >> run_calc.log
python dosplot.py
echo -e "PLOT READ --- look for dosplot.png">> run_calc.log 

# ----------------- MAIN PROCESS -----------------------------------------------

rm -f hostfile.$JOB_ID
