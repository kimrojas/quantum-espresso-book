#!/usr/bin/env python3
import os
import time
from ase.io import read, write
from argparse import ArgumentParser
from ase.calculators.espresso import Espresso
from ase.build import molecule
from ase.optimize import BFGS
from ase.units import Bohr, Rydberg

parser = ArgumentParser()
parser.add_argument("-i", "--run_idx", type=int, default=None)
args = parser.parse_args()

# SYSTEM - - - -
objdict = {
    1: (molecule("H2"), "H2"),
    2: (read("struct_files/h2v_alpha.poscar"), "h2v_alpha"),
    3: (read("struct_files/h2v_beta.poscar"), "h2v_beta"),
    4: (read("struct_files/h2v_beta_alt.poscar"), "h2v_beta_alt"),
    5: (read("struct_files/h2v_delta.poscar"), "h2v_delta"),
    6: (read("struct_files/h2v_delta_alt.poscar"), "h2v_delta_alt"),
    7: (read("struct_files/h2v_epsilon.poscar"), "h2v_epsilon"),
    8: (read("struct_files/h2v_epsilon_alt.poscar"), "h2v_epsilon_alt"),
    9: (read("struct_files/h2v_gamma.poscar"), "h2v_gamma"),
    10: (read("struct_files/h2v_gamma_alt.poscar"), "h2v_gamma_alt"),
}
obj, label = objdict.get(args.run_idx)

if label == "H2":
    obj.set_cell([15, 15, 15])
    obj.center()

obj.rattle(stdev=0.001, seed=0)


# CALCULATOR - - - -
command = f"mpirun pw.x -in PREFIX.pwi > PREFIX.pwo"  ### (CHANGE)
os.environ["ASE_ESPRESSO_COMMAND"] = command
os.environ["ESPRESSO_PSEUDO"] = "/home/krojas/QE/pseudo"  ### (CHANGE)
os.environ["ESPRESSO_TMPDIR"] = "./TMPDIR"

# -- Declare pseudopotentials --
pseudopotentials = {"B": "b_pbe_v1.4.uspp.F.UPF", "H": "h_pbe_v1.4.uspp.F.UPF"}
# -- Declare Kpoint mesh --
kpts = [5, 5, 1]

# -- Declare input settings --
input_data = {
    "control": {
        "calculation": "scf",
        "prefix": "pwscf",
        "verbosity": "high",
        "tprnfor": True,
        "tstress": True,
    },
    "system": {
        "ecutwfc": 60,
        "ecutrho": 480,
        "occupations": "smearing",
        "smearing": "gaussian",
        "degauss": 0.01,
        "input_dft": "vdw-df2-b86r",
    },
    "electrons": {
        "electron_maxstep": 100,
        "conv_thr": 1.0e-09,
        "mixing_beta": 0.7,
        "mixing_mode": "plain",
        "diagonalization": "david",
    },
}

# -- Declare calcutor
workdir = f"run_files/{label}_dir"
full_label = f"{workdir}/espresso"
calculator = Espresso(
    pseudopotentials=pseudopotentials,
    input_data=input_data,
    kpts=kpts,
    label=full_label,
)
obj.calc = calculator

# - - - CONVERGENCE THRESHOLD / FMAX - - -
forc_conv_thr = 1e-4
fmax = forc_conv_thr * Rydberg / Bohr
print("Force threshold:")
print(f"forc_conv_thr =   {forc_conv_thr} Ry/Bohr")
print(f"fmax =  {fmax} eV/A")  # 0.0025711033545240325
print()


# -- NEW IMPLEMENTATION --
# RELAXATION VIA ASE (BFGS)
os.makedirs(workdir, exist_ok=True)
logfile = f"{workdir}/relax.log"
trajfile = f"{workdir}/relax.traj"

start_time = time.time()
opt = BFGS(obj, logfile=logfile, trajectory=trajfile)
opt.run(fmax=fmax)
elapsed_time = time.time() - start_time

with open(logfile, "a") as log:
    log.write("\n")
    log.write(f"Elapsed time (s): {elapsed_time:.4f}\n")
    log.write(f"Final energy (eV): {obj.get_potential_energy()}\n")
    log.write("\n")

# -- Cleanly save the final structures
savedir = "optimized_struct_files"
os.makedirs(savedir, exist_ok=True)
write(os.path.join(savedir, label + ".poscar"), obj)
