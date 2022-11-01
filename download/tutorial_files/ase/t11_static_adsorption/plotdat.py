import numpy as np
import matplotlib.pyplot as plt

x, y = np.genfromtxt('run_calc.dat', unpack=True)

fig, ax = plt.subplots(figsize=(6,5))
ax.plot(x, y, ls='--', marker='o')


# Minima
imin = np.where(y == np.min(y))
xmin, ymin = x[imin], y[imin]
ax.scatter(xmin, ymin, marker='X', s=140, alpha=0.75, 
           color='r', label=f"Minima X = {xmin[0]}")
print(f"Minima at X = {xmin[0]}")


ax.set_ylabel('Energy (eV)', fontsize=14)
ax.set_xlabel('height (A)', fontsize=14)
ax.legend()


fig.tight_layout()
fig.savefig('DATAPLOT.png', dpi=300)

