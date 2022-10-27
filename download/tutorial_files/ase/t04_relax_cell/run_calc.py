import os
import sys
import shutil
import numpy as np
from ase.io import read, write
from ase.calculators.espresso import Espresso
from datetime import datetime

# LOAD ENVIRONMENT =============================================================
command = f"mpirun pw.x -in PREFIX.pwi > PREFIX.pwo"
os.environ['ASE_ESPRESSO_COMMAND'] = command
os.environ['ESPRESSO_PSEUDO'] ='/home/krojas/QE/pseudo'
os.environ['ESPRESSO_TMPDIR'] ='./TMPDIR'

# CALCULATOR PARTS =============================================================
# -- Declare pseudopotentials --
pseudopotentials = {'B': 'b_pbe_v1.4.uspp.F.UPF',
                    'H': 'h_pbe_v1.4.uspp.F.UPF'}
# -- Declare Kpoint mesh --
kpts = [16,16,1]
# -- Declare label --
label = "run_files/espresso"
# -- Declare input settings --
input_data = {'control': {'calculation': 'vc-relax',
                          'prefix': 'pwscf',
                          'verbosity': 'high',
                          'tprnfor':True,
                          'tstress':True,
                          'etot_conv_thr': 1.0e-05,
                          'forc_conv_thr': 1.0e-04,
                          'nstep': 1000
                          },
              'system': {'ecutwfc': 60,
                         'ecutrho': 480,
                         'occupations': 'smearing',
                         'smearing': 'gaussian',
                         'degauss': 0.01,
                         'input_dft': 'vdw-df2-b86r'
                         },
              'electrons': {'electron_maxstep': 100,
                            'conv_thr': 1.e-09,
                            'mixing_beta': 0.7,
                            'mixing_mode': 'plain',
                            'diagonalization': 'david'
                            },
              'ions': {'ion_dynamics': 'bfgs'},
              'cell': {'cell_dofree': '2Dxy'}
              }

# CREATE ASE ATOMS OBJECT ======================================================
    # Create atoms object from structure file
atoms = read('unit_cell.poscar')
    # Create calculator object
calculator = Espresso(pseudopotentials=pseudopotentials,
                      input_data=input_data,
                      kpts=kpts,
                      label=label)
    # Connect atoms object and calculator object
atoms.calc = calculator


# RUN CALCULATION ==============================================================
with open('run_calc.log', 'w') as f:
    print(f"Start {datetime.now()}", file=f)

    # Run calculation
    calculator.calculate(atoms)

    print(f"End {datetime.now()}", file=f)

    # Save structure
    opt_atoms = read(label+'.pwo', index=':')
    write('run_files/all_structures.traj', opt_atoms)
    write('run_files/final_structure.traj', opt_atoms[-1])
    write('run_files/final_structure.poscar', opt_atoms[-1])

    # Save Energy
    print(f"INITIAL ENERGY {opt_atoms[0].get_potential_energy()} eV", file=f)
    print(f"FINAL ENERGY {opt_atoms[-1].get_potential_energy()} eV", file=f)
