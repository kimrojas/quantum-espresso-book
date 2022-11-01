import os
import sys
import shutil
import subprocess
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
input_data = {'control': {'calculation': 'scf',
                          'prefix': 'pwscf',
                          'verbosity': 'high',
                          'tprnfor':True,
                          'tstress':True,
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
                            'diagonalization': 'cg'
                            },
              }

# CREATE ASE ATOMS OBJECT ======================================================
    # Create atoms object from structure file
atoms = read('unit_cell_rect.poscar')
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



    # Loop over many KPOINTS
    for k in range(1, 18+1):
        calculator.set(kpts=[k,16,1])
        calculator.set_label(label=label+f"_{k:02}x")

        calculator.calculate(atoms)
        energy = atoms.get_potential_energy()
        # SAVE ENERGY
        subprocess.call(f"echo -e '{k:<4} {energy}' >> run_calc_x.dat", shell=True)


    print(f"End {datetime.now()}", file=f)
    
with open('run_calc.log', 'w') as f:

    print(f"Start {datetime.now()}", file=f)



    # Loop over many KPOINTS
    for k in range(1, 18+1):
        calculator.set(kpts=[16,k,1])
        calculator.set_label(label=label+f"_{k:02}y")

        calculator.calculate(atoms)
        energy = atoms.get_potential_energy()
        # SAVE ENERGY
        subprocess.call(f"echo -e '{k:<4} {energy}' >> run_calc_y.dat", shell=True)


    print(f"End {datetime.now()}", file=f)

