import sys, os, subprocess
from pprint import pprint
import numpy as np
from argparse import ArgumentParser

from ase.calculators.espresso import Espresso
from ase.build import molecule
from ase.io import read, write


parser = ArgumentParser()
parser.add_argument("-i", "--run_idx", type=int, default=None)
args = parser.parse_args()


# REFERENCES
# Reference positions
def midpoint(a, b):
    return np.mean([a, b], axis=0)


ref = read("super_cell_rect.poscar").get_positions()
siteref = {
    "Btop": ref[7],
    "Htop": ref[40],
    "Bridge": midpoint(ref[7], ref[5]),
    "Hollow": midpoint(ref[20], ref[5]),
}

# Reference rotations
rotref = {
    "A": dict(phi=0, theta=0, psi=0, center="COP"),
    "B": dict(phi=0, theta=90, psi=0, center="COP"),
    "C": dict(phi=0, theta=90, psi=45, center="COP"),
}

# Reference heights
heightref = {
    "Btop_A": 3.50,
    "Btop_B": 3.50,
    "Btop_C": 3.50,
    "Htop_A": 2.75,
    "Htop_B": 2.75,
    "Htop_C": 2.75,
    "Bridge_A": 3.50,
    "Bridge_B": 3.50,
    "Bridge_C": 3.50,
    "Hollow_A": 3.25,
    "Hollow_B": 3.25,
    "Hollow_C": 3.25,
}


# -- SETTINGS ------------------------------------------------------------------
all_settings = []
for s in siteref.keys():
    for r in rotref.keys():
        all_settings.append((s, r))

site, rot = all_settings[args.run_idx - 1]
prefix = "run_files"
settings = {"site": site, "rot": rot}

command = f"mpirun pw.x -in PREFIX.pwi > PREFIX.pwo"
os.environ["ASE_ESPRESSO_COMMAND"] = command
os.environ["ESPRESSO_PSEUDO"] = "/home/krojas/QE/pseudo"
os.environ["ESPRESSO_TMPDIR"] = "./TMPDIR"
# ------------------------------------------------------------------------------


# --- FUNCTIONS PART -----------------------------------------------------------

# CREATE ADSORBATE FUNCTION
def create_adsorbate(site, height, rotation):
    adsorbate = molecule("H2")
    adsorbate.center(about=(siteref[site] + [0, 0, height]))
    adsorbate.euler_rotate(**rotation)
    adsorbate.rattle(0.07, seed=0)
    return adsorbate


# CREATE ADSORPTION STRUCTURE FUNCTION
def create_adsorption_structure(opt):
    # Read substrate
    atoms = read("super_cell_rect.poscar")

    # Combine substrate and adsorbate
    adsorbate = create_adsorbate(
        site=opt["site"], height=opt["height"], rotation=rotref[opt["rot"]]
    )
    atoms.extend(adsorbate)

    # Save
    strsite = f"{opt['site']}"
    strheight = f"h_{opt['height']:05.2f}"
    strrot = f"rot{opt['rot']}"

    dirname = f"{prefix}_{strsite}_{strrot}"
    fname = f"{strheight}"
    return dirname, fname, atoms


def run_calculation(dirname, fname, atoms, h):
    # -- Declare pseudopotentials --
    pseudopotentials = {"B": "b_pbe_v1.4.uspp.F.UPF", "H": "h_pbe_v1.4.uspp.F.UPF"}
    # -- Declare Kpoint mesh --
    kpts = [5, 5, 1]

    # -- Declare input settings --
    input_data = {
        "control": {
            "calculation": "relax",
            "prefix": "pwscf",
            "verbosity": "high",
            "tprnfor": True,
            "tstress": True,
            "nstep": 0, #1000,
            "etot_conv_thr": 1.0e-5,
            "forc_conv_thr": 1.0e-4,
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
        "ions": {"ion_dynamics": "bfgs"},
    }
    # -- Declare label --
    full_label = os.path.join(dirname, "espresso_" + fname)
    calculator = Espresso(
        pseudopotentials=pseudopotentials,
        input_data=input_data,
        kpts=kpts,
        label=full_label,
    )
    # Connect atoms object and calculator object
    atoms.calc = calculator

    try:
        calculator.calculate(atoms)
    except Exception as e:
        print(e)

    # SAVE ENERGY
    subprocess.call(
        f"echo -e 'Relaxation done' >> {dirname}/run_calc.dat", shell=True
    )

    # Read output 
    opt_atoms = read(full_label+'.pwo', index=':')
    write(f'{dirname}/all_{fname}.traj', opt_atoms)
    write(f'{dirname}/opt_{fname}.traj', opt_atoms[-1])
    write(f'{dirname}/opt_{fname}.poscar', opt_atoms[-1])
    subprocess.call(
                f"echo -e 'FINAL ENERGY :: {opt_atoms[-1].get_potential_energy()} eV' >> {dirname}/run_calc.dat", shell=True
            )

# --- RUN PART -----------------------------------------------------------------

# RUN ROTATION CONFIG
height_key = "_".join([settings["site"], settings["rot"]])
settings["height"] = heightref[height_key]
dirname, fname, atoms = create_adsorption_structure(opt=settings)
run_calculation(dirname, fname, atoms, settings["height"])

