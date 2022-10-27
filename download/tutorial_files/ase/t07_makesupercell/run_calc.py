import os
import numpy as np
from ase.io import read, write
from ase.build import sort


hb_unit = read('unit_cell.poscar')

# Expand
hb_super = hb_unit.repeat((3,3,1))

#Center (xy plane)
hb_super.center(vacuum=None, axis=(0,1))

# Center (z-axis)
hb_super.center(vacuum=18/2, axis=(2))
                
hb_super = sort(hb_super)

# Save
os.makedirs('run_files', exist_ok=True)
write('run_files/supercell.poscar', hb_super)

