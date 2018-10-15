"""
Script.

Plot the results of find_domain.py
"""
import numpy as np
import matplotlib.pyplot as plt
import sys

# Import params from parent directory
# DO NOT import from this directory as the path to the interaction potential will be wrong
path_to_params = '..'
sys.path.insert(0, path_to_params) # looks in parent directory first
from params import *
sys.path.remove(path_to_params)

# Set cutoff - throw away (1-cutoff) * sum(data) data points at the far edge
cutoff_list = np.array([0.999, 0.9999])

# Lines for cutoffs
lines = ['m-', 'r-']

# --------------------------------------------- #
#  Remind the user what parameters we're using  #
# --------------------------------------------- #
print "Lattice alpha: %s,  lattice beta: %s" %(alpha_type, beta_type)
print "Number of atoms per lattice: %d" %Natoms
print "Density: %.4f" %density
print "Temperature: %.4f" %B
print "Max step size: %.4f" %dr_max
print ""


# Load input files
lhs = np.loadtxt("mu_range_s0.out")
rhs = np.loadtxt("mu_range_s1.out")

# Get counts
lhist = lhs[:,0]
rhist = rhs[:,0]

# Get mu values for bins
lbins = lhs[:,1]
rbins = rhs[:,1]

# Tell the user what's going on
print "Mu = E(%s) - E(%s)" %(alpha_type, beta_type)
print "Most negative bin has %d counts at mu = %f" %(lhist[-1], lbins[-1])
print "Most positive bin has %d counts at mu = %f" %(rhist[-1], rbins[-1])

# Plot parameters
plt.rc('text',usetex=True)
font = {'family' : 'serif',
        'size' : 14}
plt.rc('font', **font)
plt.rcParams['axes.linewidth'] = 2

# Create plot
fig, ax = plt.subplots()
ax.set_title("Histogram generated with Boltzmann weights")
ax.set_xlabel(r"$\mu$")
ax.set_ylabel("Counts")
ax.tick_params(direction='in', top=True, right=True)

# Add data to plot
ax.plot(lbins, lhist, 'b.')
ax.plot(rbins, rhist, 'g.')

for i in range(len(cutoff_list)):
    
    cutoff = cutoff_list[i]
    cut_pct = 100*cutoff
    line = lines[i%len(lines)]

    # Find the percentiles for the lhs and 99.9th for the rhs
    lcut = lbins[ np.where(np.cumsum(lhist) > cutoff*np.sum(lhist))[0][0] ]
    rcut = rbins[ np.where(np.cumsum(rhist) > cutoff*np.sum(rhist))[0][0] ]

    print "%.2f percent of the data lies within the range  (%f, %f)" \
            %(cut_pct, lcut, rcut)

    ax.plot([lcut,lcut], [0,0.2*np.max(lhist)], line, label="%.2fth percentile" %cut_pct)
    ax.plot([rcut,rcut], [0,0.2*np.max(lhist)], line)

ax.legend(loc=1)
plt.show()
