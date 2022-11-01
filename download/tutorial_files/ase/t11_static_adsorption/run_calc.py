import sys, os, subprocess
from pprint import pprint
import numpy as np
from argparse import ArgumentParser

from ase.calculators.espresso import Espresso
from ase.build import molecule
from ase.io import read, write


parser = ArgumentParser()
parser.add_argument('-i', '--run_idx', type=int, default=None)
args = parser.parse_args()


# REFERENCES
# Reference positions
def midpoint(a,b): return np.mean([a,b], axis=0)
ref = read('super_cell_rect.poscar').get_positions()
siteref = {'Btop': ref[7],
           'Htop': ref[40],
           'Bridge': midpoint(ref[7], ref[5]),
           'Hollow': midpoint(ref[20], ref[5])}

# Reference rotations
rotref = {'A': dict(phi=0, theta=0, psi=0, center='COP'),
          'B': dict(phi=0, theta=90, psi=0, center='COP'),
          'C': dict(phi=0, theta=90, psi=45, center='COP')}


# -- SETTINGS ------------------------------------------------------------------
prefix = 'run_files'
settings = {'site': list(siteref.keys())[args.run_idx-1]}
hrange = np.arange(1.5, 5.25, 0.25)



command = f"mpirun pw.x -in PREFIX.pwi > PREFIX.pwo"
os.environ['ASE_ESPRESSO_COMMAND'] = command
os.environ['ESPRESSO_PSEUDO'] ='/home/krojas/QE/pseudo'
os.environ['ESPRESSO_TMPDIR'] ='./TMPDIR'
# ------------------------------------------------------------------------------


# --- FUNCTIONS PART -----------------------------------------------------------

# CREATE ADSORBATE FUNCTION
def create_adsorbate(site, height, rotation):
    adsorbate = molecule('H2')
    adsorbate.center(about=(siteref[site]+[0,0,height]))
    adsorbate.euler_rotate(**rotation)
    adsorbate.rattle(0.07, seed=0)
    return adsorbate

# CREATE ADSORPTION STRUCTURE FUNCTION
def create_adsorption_structure(opt):
    # Read substrate
    atoms = read('super_cell_rect.poscar')

    # Combine substrate and adsorbate
    adsorbate = create_adsorbate(site=opt['site'],
                                 height=opt['height'],
                                 rotation=rotref[opt['rot']])
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
    pseudopotentials = {'B': 'b_pbe_v1.4.uspp.F.UPF',
                        'H': 'h_pbe_v1.4.uspp.F.UPF'}
    # -- Declare Kpoint mesh --
    kpts = [5,5,1]

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
    # -- Declare label --
    full_label = os.path.join(dirname, 'espresso_'+fname)
    calculator = Espresso(pseudopotentials=pseudopotentials,
                      input_data=input_data,
                      kpts=kpts,
                      label=full_label)
    # Connect atoms object and calculator object
    atoms.calc = calculator
    
    try:
        calculator.calculate(atoms)
        energy = atoms.get_potential_energy()
        
        # SAVE ENERGY
        subprocess.call(f"echo -e '{h:<10} {energy}' >> {dirname}/run_calc.dat", shell=True)
    except Exception as e:
        print(e)


def make_plot(dirname):
    cwd = os.getcwd()
    os.chdir(dirname)
    minima = subprocess.check_output("python ../plotdat.py", shell=True).decode('ascii').strip()
    os.chdir(cwd)
    return minima
    



# --- RUN PART -----------------------------------------------------------------

# RUN ROTATION CONFIG -- A --
settings['rot'] = 'A'
for h in hrange:
    settings['height'] = h
    dirname, fname, atoms = create_adsorption_structure(opt=settings)
    run_calculation(dirname, fname, atoms, h)
minima = make_plot(dirname)
subprocess.call(f"echo -e '{dirname} {minima}' >> run_calc.log")

# RUN ROTATION CONFIG -- B --
settings['rot'] = 'B'
for h in hrange:
    settings['height'] = h
    dirname, fname, atoms = create_adsorption_structure(opt=settings)
    run_calculation(dirname, fname, atoms, h)
minima = make_plot(dirname)
subprocess.call(f"echo -e '{dirname} {minima}' >> run_calc.log")


# RUN ROTATION CONFIG -- C --
settings['rot'] = 'C'
for h in hrange:
    settings['height'] = h
    dirname, fname, atoms = create_adsorption_structure(opt=settings)
    run_calculation(dirname, fname, atoms, h)
minima = make_plot(dirname)
subprocess.call(f"echo -e '{dirname} {minima}' >> run_calc.log")

