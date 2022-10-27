import numpy as np
import matplotlib.pyplot as plt
from pprint import pprint
from io import StringIO
import subprocess

# ------- SETTINGS ---------------
bandfile = 'bands.out.gnu'
pwo = 'espresso_scf.pwo'
symmetry_file = 'espresso_bandx.out'
labels = [r"$\Gamma$", "Y", "K\'", r"$\Gamma$", "S", "K", r"$\Gamma$"]
# (labels are high symmetry points labels)
# --------------------------------


# READ DATA --------------------------------------------------------------------

# Get Fermi energy level
line = subprocess.check_output(f"grep 'Fermi energy' {pwo}", shell=True)
line = line.decode('ascii').split()
efermi = np.float64(line[-2])

# Read band file
with open(bandfile, 'r') as f:
    contents = f.read()
    # Split each band
    contents = contents.split(" \n")
    # Remove empty line
    contents = [i for i in contents if i]

# Parse file contents into numpy arrays then store to a list
# (Energy is centered to fermi energy)
bands = []
for i, content in enumerate(contents):
    k, e = np.genfromtxt(StringIO(content), unpack=True)
    bands.append((k, e-efermi))

# Read Symmetry file and find kpath mapping to x-axis
with open(symmetry_file, 'r') as f:
    lines = f.read().splitlines()
    symstring = [line for line in lines if ('high-symmetry point' in line)]
    sympoints = np.array([float(i.split()[-1]) for i in symstring])

# PLOT -------------------------------------------------------------------------

fig, ax = plt.subplots(figsize=(6,5))

# Plot bands
for x,y in bands:
    ax.plot(x, y, ls='-', color='k', alpha=0.8)
    
# Plot Fermi level
ax.axhline(y=0, ls='--', color='r', alpha=0.5)

# Plot high-symmetry points
for x in sympoints[1:-1]:
    ax.axvline(x=x, ls=':', color='k', alpha=0.3)

# Plot information and labels
ax.set_xticks(sympoints, labels, fontsize=16)
ax.set_xlim(0, sympoints[-1])
ax.set_ylabel(r"$E-E_F$ (eV)", fontsize=16)
fig.tight_layout()

# SAVE FIGURE ------------------------------------------------------------------
fig.savefig('bandstructureplot.png', dpi=300)