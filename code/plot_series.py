"""
Script.

How to use:
        argv[1] is the data file
    
    e.g. python plot_series series_data.txt

Requires correct params file used to generate input data.
"""

import numpy as np
import matplotlib.pyplot as plt
from sys import argv

from params import *

# Input file given as argument variable
input_file = argv[1]

# Load input data
input_data = np.loadtxt(input_file)

# Unpack
subdom = input_data[:,0]
mu = input_data[:,1]
weights = input_data[:,2]
E_active = input_data[:,3]

# Find rows of -1's denoting new, independent simulations
restart_rows = np.where( subdom == -1 )[0]
if len(restart_rows) > 1:
    stop_rows = list(restart_rows[1:]) + [None]
else:
    stop_rows = [None]


# Plot parameters
plt.rc('text',usetex=True)
font = {'family' : 'serif',
        'size' : 14}
plt.rc('font', **font)
plt.rcParams['axes.linewidth'] = 2

# Create plots
fig, ax = plt.subplots()
ax.set_title("Mu series")
ax.set_xlabel("Sweeps")
ax.set_ylabel(r"$\mu$")

fig2, ax2 = plt.subplots()
ax2.set_title("Active lattice energy series")
ax2.set_xlabel("Sweeps")
ax2.set_ylabel("Active lattice energy")

# Iterate over independent simulations and add to plot
N_prev = 0
for r in range(len(restart_rows)):
    ilo = restart_rows[r] + 1
    ihi = stop_rows[r]

    mu_r = mu[ilo:ihi]
    E_active_r = E_active[ilo:ihi]
    N_samples = len(mu_r)
    sweeps = (np.arange(N_samples)+N_prev) * sweeps_series / float(Natoms)

    ax.plot(sweeps, mu_r)
    #ax.plot([sweeps[-1],sweeps[-1]], [np.min(mu_r),np.max(mu_r)], 'k')
    ax2.plot(sweeps, E_active_r)
    #ax2.plot([sweeps[-1],sweeps[-1]], [np.min(E_active_r),np.max(E_active_r)], 'k')

    N_prev += N_samples


plt.show()
