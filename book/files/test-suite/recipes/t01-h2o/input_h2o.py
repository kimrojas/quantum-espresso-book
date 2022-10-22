#!/usr/bin/env python

import os
from ase.calculators.dftb import Dftb
from ase.build import molecule

os.environ['DFTB_PREFIX'] = '/home/krojas/share/lib/slakos/matsci'
os.environ['ASE_DFTB_COMMAND'] = 'dftb+ > PREFIX.out'


# --- DEFINE SYSTEM ---
atoms = molecule('H2O')

# --- DEFINE CALCULATOR ---
calc = Dftb(label='h2o',
            Hamiltonian_SCC='Yes',
            Hamiltonian_MaxAngularMomentum_='',
            Hamiltonian_MaxAngularMomentum_H='"s"',
            Hamiltonian_MaxAngularMomentum_O='"p"',
            Hamiltonian_Charge='0.000000',
            Hamiltonian_Filling='Fermi {',
            Hamiltonian_Filling_empty='Temperature [Kelvin] = 0.000000')

# --- ATTACH CALCULATOR AND ATOMS
atoms.calc = calc

# --- CALCULATE ENERGY ---
en = atoms.get_potential_energy()

# --- OUTPUT ---
print(f"The energy calculated is: {en} eV")