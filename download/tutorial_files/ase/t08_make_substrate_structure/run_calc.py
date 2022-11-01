##
# Creates the unit cell and a supercell (rectangular cell)
##
import numpy as np
from ase.io import read, write
from pprint import pprint
from ase import Atoms
from ase.build import sort
import os

# SUPER CELL SETTINGS --------------
rep_x = 2
rep_y = 3
vacuum = 18
# ----------------------------------




def get_dist(a,b):
    return (np.linalg.norm(a-b))

def get_xydist(a,b):
    return (np.linalg.norm(a[:-1] - b[:-1]))



# CREATE UNIT CELL
base = read('base.poscar')
base_pos = base.get_positions()
base_sym = base.get_chemical_symbols()

pprint(base_pos)
pprint(base_sym)

# Bond length measurements
BB = get_dist(base_pos[0], base_pos[1])
BHxy = get_xydist(base_pos[1], base_pos[2])
HH = get_dist(base_pos[2], base_pos[3])
BBlong = get_dist(base_pos[8], base_pos[5])

cell_y = get_dist(base_pos[0], base_pos[8])
cell_x = (BHxy*2) + BBlong
cell_z = vacuum

print('BB', BB)
print('BHxy', BHxy)
print('HH', HH)
print('cell_x', cell_x)
print('cell_y', cell_y)
print('cell_z', cell_z)

# Creation
cell = [cell_x, cell_y, cell_z]
pos = []
sym = []

# atom 0 (B)
pos.append([BHxy, 0., cell_z/2])
sym.append('B')
# atom 1 (B)
pos.append([cell_x-BHxy, 0., cell_z/2])
sym.append('B')
# atom 2 (B)
pos.append([BBlong/2, cell_y/2, cell_z/2])
sym.append('B')
# atom 3 (B)
pos.append([cell_x-(BBlong/2), cell_y/2, cell_z/2])
sym.append('B')
# atom 4 (H)
pos.append([0., 0., (cell_z/2)+(HH/2)])
sym.append('H')
# atom 5 (H)
pos.append([0., 0., (cell_z/2)-(HH/2)])
sym.append('H')
# atom 6 (H)
pos.append([cell_x/2., cell_y/2., (cell_z/2)+(HH/2)])
sym.append('H')
# atom 7 (H)
pos.append([cell_x/2., cell_y/2., (cell_z/2)-(HH/2)])
sym.append('H')


os.makedirs('run_files', exist_ok=True)
unit = Atoms(symbols=sym, positions=pos, cell=cell)
write('run_files/unit_cell_rect.poscar', unit)


# CREATE SUPER CELL
super = unit.repeat((rep_x,rep_y,1))
# super.center(vacuum=vacuum/2, axis=(2))
super = sort(super)
write('run_files/super_cell_rect.poscar', super)


# hb_unit = read('unit_cell.poscar')

# # Expand
# hb_super = hb_unit.repeat((3,3,1))

# #Center (xy plane)
# hb_super.center(vacuum=None, axis=(0,1))

# # Center (z-axis)
# hb_super.center(vacuum=18/2, axis=(2))
                
# hb_super = sort(hb_super)

# # Save
# os.makedirs('run_files', exist_ok=True)
# write('run_files/supercell.poscar', hb_super)
