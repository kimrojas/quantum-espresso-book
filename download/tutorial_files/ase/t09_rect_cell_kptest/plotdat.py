import numpy as np
import matplotlib.pyplot as plt

def plot_data(fname, plotname, label):
    # Read data
    x, y = np.genfromtxt(fname, unpack=True)

    # Transform data 
    # (we will plot energy difference instead of energy value)
    x = x[:-1]
    y = np.abs(y[1:] - y[:-1])



    # # Plot
    fig, ax = plt.subplots()
    ax.plot(x, y, ls='--', marker='*')
    ax.axhline(0.001, ls=':')
    ax.set_ylim(0, 0.02)
    ax.set_xticks(x)
    ax.set_ylabel('Energy difference (eV)', fontsize=14)
    ax.set_xlabel(label, fontsize=14)
    fig.tight_layout()

    # SAVE FIGURE
    fig.savefig(plotname, dpi=300)


# -------------------------------
plot_data('run_calc_x.dat', 'plot_kpX.png', label='KX')
plot_data('run_calc_y.dat', 'plot_kpY.png', label='KY')
