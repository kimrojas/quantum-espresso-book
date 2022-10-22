#!/usr/bin/env python

import os
from ase.calculators.dftb import Dftb
from ase.io import read


os.environ['DFTB_PREFIX'] = '/home/krojas/share/lib/slakos/matsci'
os.environ['ASE_DFTB_COMMAND'] = 'dftb+ > PREFIX.out'


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
