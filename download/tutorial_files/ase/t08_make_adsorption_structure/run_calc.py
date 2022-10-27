from ase.build import molecule
from ase.io import read, write
from ase.visualize import view

# -- SETTINGS --
settings = {'site': 'A',
            'height': 2,
            'rotation_angle': 90,
            'rotation_axis': 'x'}

# --------------


position = {'A': [5.783956581259004, 4.814392921230736, 9.96739506],
            'B': [5.777748202759003, 3.0822587927307366, 9.96739506],
            'C': [7.290650588759004, 4.809093124730737, 10.934790119999999]}

# Read substrate
atoms = read('supercell.poscar')


# Create adsorbate
def add_adsorbate(site, height, rotation_angle, rotation_axis):
    adsorbate = molecule('H2')
    adsorbate.set_cell(atoms.get_cell())
    x,y,z = position[site]
    adsorbate.center(about=(x, y, z+height))
    adsorbate.rotate(a=rotation_angle, v=rotation_axis, center='COP')
    
    return adsorbate

# Combine substrate and adsorbate
atoms.extend(add_adsorbate(**settings))

# Save
write(f"hb_h2_site{settings['site']}_h{settings['height']}_rot{settings['rotation_angle']}{settings['rotation_axis']}.poscar", atoms)