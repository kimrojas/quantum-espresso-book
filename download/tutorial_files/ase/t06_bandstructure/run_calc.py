import os
import sys
import shutil
import subprocess
import numpy as np
from ase.io import read, write
from ase.calculators.espresso import Espresso
from datetime import datetime
from ase.dft.kpoints import bandpath


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

    # NSCF (bands) CALCULATION
    print(f"Start NSCF {datetime.now()}", file=f)
    # (update kpts and calculation type)
    input_data['control']['calculation'] = 'bands'
    input_data['control']['restart_mode'] = 'restart'

    spath =    {'G' :[0.0000000000,     0.0000000000,     0.0000000000],
                'Y' :[0.0000000000,     0.5000000000,     0.0000000000],
                "K" :[0.3375000000,     0.6708333333,     0.0000000000],
                'G' :[0.0000000000,     0.0000000000,     0.0000000000],
                "S" :[0.5000000000,     0.5000000000,     0.0000000000],
                'K' :[0.6625000000,     0.3291666667,     0.0000000000],
                'G' :[0.0000000000,     0.0000000000,     0.0000000000]}

    xpath =[[0.0000000000,     0.0000000000,     0.0000000000],
            [0.0000000000,     0.5000000000,     0.0000000000],
            [0.3375000000,     0.6708333333,     0.0000000000],
            [0.0000000000,     0.0000000000,     0.0000000000],
            [0.5000000000,     0.5000000000,     0.0000000000],
            [0.6625000000,     0.3291666667,     0.0000000000],
            [0.0000000000,     0.0000000000,     0.0000000000]]

    kpts = bandpath(path=xpath, special_points=spath,
                npoints=7*50, cell=atoms.get_cell())

    calculator.set(kpts=kpts,
                   input_data=input_data)
    calculator.set_label(label+"_bands")
    calculator.calculate(atoms)
    print(f"End NSCF {datetime.now()}", file=f)

    # BANDS CALCULATION
    print(f"Start BANDS {datetime.now()}", file=f)
    os.chdir('run_files')
    # (create bands input file manually, ASE doesn't have this)
    with open('espresso_bandx.inp', 'w') as finp:
        lines = [" &BANDS",
                 "    prefix = 'pwscf'",
                 " /"]
        for line in lines:
            finp.write(line+"\n")


    # (run dos calculation)
    subprocess.call("mpirun bands.x inp-n 1 bands.x -in espresso_bandx.inp > espresso_bandx.out", shell=True)
    os.chdir('../')
    print(f"End BANDS {datetime.now()}", file=f)