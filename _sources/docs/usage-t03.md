# HB sheet system - parallelized

A parallelized calculation for an Hydrogen boride (HB) sheet using DFTB+ with the ASE interface. Files are available in the `test-suite` provided.

We control the parallelization using the `OMP_NUM_THREADS` variable.


```python
#!/usr/bin/env python

import os
from ase.calculators.dftb import Dftb
from ase.io import read

# --- SET THREADS ---
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('-t', '--threads', type=int, default=1)
args = parser.parse_args()

os.environ['DFTB_PREFIX'] = '/home/krojas/share/lib/slakos/matsci'
os.environ['ASE_DFTB_COMMAND'] = f'OMP_NUM_THREADS={args.threads} dftb+ > PREFIX.out'


# --- DEFINE SYSTEM ---
atoms = read('structure.traj')

# --- DEFINE CALCULATOR ---
calc = Dftb(label='hb',
            Hamiltonian_SCC='Yes',
            Hamiltonian_MaxAngularMomentum_='',
            Hamiltonian_MaxAngularMomentum_B='"p"',
            Hamiltonian_MaxAngularMomentum_H='"s"',
            Hamiltonian_Charge='0.000000',
            Hamiltonian_Filling='Fermi {',
            Hamiltonian_Filling_empty='Temperature [Kelvin] = 0.000000',
            kpts=[8,8,1])

# --- ATTACH CALCULATOR AND ATOMS
atoms.calc =  calc

# --- CALCULATE ENERGY ---
en = atoms.get_potential_energy()

# --- OUTPUT ---
print(f"The energy calculated is: {en} eV")
```



