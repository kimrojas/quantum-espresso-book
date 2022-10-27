# ------ DESCRIPTION---------
# This script is used to create a DOS plot

import numpy as np
import matplotlib.pyplot as plt
import subprocess
import os, sys

# FILES TO READ
total_file = 'pwscf.pdos_tot'
partial_files = ['pwscf.pdos_atm#1(B)_wfc#1(s)', 'pwscf.pdos_atm#1(B)_wfc#2(p)',
                 'pwscf.pdos_atm#2(B)_wfc#1(s)', 'pwscf.pdos_atm#2(B)_wfc#2(p)',
                 'pwscf.pdos_atm#3(H)_wfc#1(s)', 'pwscf.pdos_atm#4(H)_wfc#1(s)']

# READ FERMI ENERGY LEVEL
line = subprocess.check_output("grep 'Fermi energy' espresso_scf.pwo", shell=True)
efermi = float(line.decode('ascii').split()[-2])

# QUICK FUNCTION TO READ EACH FILE
# (We center the energy with respect to the fermi level)
def get_data(fname):
    en, dos = np.genfromtxt(fname, skip_header=1, usecols=(0,1), unpack=True)
    return en-efermi, dos #

# PLOT -------------------------
fig, ax = plt.subplots(figsize=(12,5), dpi=300)

# Plot TOTAL
en, dos = get_data(total_file)
ax.plot(en, dos, label=total_file)

# Plot partial
for file in partial_files:
    _en, _dos = get_data(file)
    ax.plot(_en, _dos, label=file)
    
# Plot fermi level
ax.axvline(x=0, ls='--', color='k', alpha=0.5)

# Plot details and formats
ax.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
ax.set_ylabel('Density of States', fontsize=12)
ax.set_xlabel('Energy (eV)', fontsize=12)
fig.tight_layout()

# Save figure
fig.savefig('dosplot.png', dpi=300)
