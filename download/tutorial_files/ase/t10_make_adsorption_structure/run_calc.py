from ase.build import molecule
from ase.io import read, write
from ase.visualize import view
import numpy as np
import sys
from pprint import pprint
import os

from argparse import ArgumentParser

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
hrange = np.arange(0.75, 4.25, 0.25)
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
    fname = f"{strheight}.poscar"
    return dirname, fname, atoms



# --- RUN PART -----------------------------------------------------------------

# RUN ROTATION CONFIG -- A --
settings['rot'] = 'A'
for h in hrange:
    settings['height'] = h
    dirname, fname, atoms = create_adsorption_structure(opt=settings)
    os.makedirs(dirname, exist_ok=True)
    write(os.path.join(dirname, fname), atoms)

# RUN ROTATION CONFIG -- B --
settings['rot'] = 'B'
for h in hrange:
    settings['height'] = h
    dirname, fname, atoms = create_adsorption_structure(opt=settings)
    os.makedirs(dirname, exist_ok=True)
    write(os.path.join(dirname, fname), atoms)

# RUN ROTATION CONFIG -- C --
settings['rot'] = 'C'
for h in hrange:
    settings['height'] = h
    dirname, fname, atoms = create_adsorption_structure(opt=settings)
    os.makedirs(dirname, exist_ok=True)
    write(os.path.join(dirname, fname), atoms)
