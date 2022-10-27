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
                            'diagonalization': 'david'
                            },
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
with open('run_calc.log', 'w', buffering=1) as f:

    # SCF CALCULATION
    print(f"Start SCF {datetime.now()}", file=f)
    calculator.set_label(label+"_scf")
    calculator.calculate(atoms)
    print(f"End SCF {datetime.now()}", file=f)

    # NSCF CALCULATION
    print(f"Start NSCF {datetime.now()}", file=f)
    # (update kpts and calculation type)
    input_data['control']['calculation'] = 'nscf'
    calculator.set(kpts=[21, 21, 1],
                   input_data=input_data)
    calculator.set_label(label+"_nscf")
    calculator.calculate(atoms)
    print(f"End NSCF {datetime.now()}", file=f)

    # DOS CALCULATION
    print(f"Start DOS {datetime.now()}", file=f)
    os.chdir('run_files')
    # (create dos input file manually, ASE doesn't have this)
    with open('projwfc.in', 'w') as finp:
        lines = [" &PROJWFC",
                 "    prefix = 'pwscf'",
                 "    Emin = -20.0",
                 "    Emax = 20.0",
                 "    DeltaE = 0.01",
                 "    ngauss = 0",
                 "    degauss = 0.01",
                 " /"]
        for line in lines:
            finp.write(line+"\n")
    
    
    # (run dos calculation)
    subprocess.call("mpirun -n 1 projwfc.x -in projwfc.in > projwfc.out", shell=True)
    os.chdir('../')
    print(f"End DOS {datetime.now()}", file=f)