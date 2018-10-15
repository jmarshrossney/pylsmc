""" 
Script.

Essential argument variables:
        File name
        Type of plot: w (weights), h (hist), 2d (matrix), log (log plot)
Optional argument variables:
        Subdomain index, e.g 1,2,3

    e.g.
        python plot_binned.py h my_global_hist.txt
        python plot_binned.py w weights_s3.txt 3
        python plot_binned.py 2d Cmat_scomb.txt
        python plot_binned.py log cuts_s2.txt

Requires correct params file used to generate input data.
"""

import numpy as np
import matplotlib.pyplot as plt
from sys import argv
from os.path import basename

from params import *
import domain as dom

# Name of this file
this_file = basename(__file__)


# Extract file name
# This will be the first arg var given with a "." (as in .out, .txt etc)
input_file = [v for v in argv if '.' in v]
if len(input_file) == 0:
    error(this_file, "No file name found in list of argument variables.")

# Load data
input_data = np.loadtxt(input_file[1]) # index 0 is the plot script itself

# Extract subdomain index if provided as argument variable 
# This will be the first int given in the argv list
s = [v for v in argv if v.isdigit()]

# If no subdomain given as argument variable, plot global data
if len(s) == 0:
    mu_bins = dom.get_global_mu_bins() * B / float(Natoms)
    s = 'all'
else:
    s = s[0]
    mu_bins = dom.get_local_mu_bins(int(s)) * B / float(Natoms)

# Check data and mu are the same size
if len(input_data) != len(mu_bins):
    error(this_file, "Length of input data does not match number of mu bins in subdomain: "+s)


# Plot parameters
plt.rc('text',usetex=True)
font = {'family' : 'serif',
        'size' : 14}
plt.rc('font', **font)
plt.rcParams['axes.linewidth'] = 2

# Create plot
fig, ax = plt.subplots()
ax.set_title("Lattice Alpha, 0, Lattice Beta")
ax.set_xlabel(r"$\mu \times \beta/N$")
ax.tick_params(direction='in', top=True, right=True)


########################
## Plotting functions ##
########################
if 'h' in argv:
    """ Plot a histogram """
    ax.set_ylabel("$h$")
    
    zeros = np.where(input_data == 0)[0]
    if len(zeros) > 0:
        print "No counts for bins ", zeros

    # Add input data
    ax.plot(mu_bins[:dom.I_cross], input_data[:dom.I_cross], 'bo-')
    ax.plot(mu_bins[dom.I_cross+1:], input_data[dom.I_cross+1:], 'go-')
    if dom.I_cross < len(mu_bins):
        ax.plot(mu_bins[dom.I_cross], input_data[dom.I_cross], 'ro')

    # Plot subdomain boundaries
    b_arr = np.array(boundaries) * B / float(Natoms)
    plt.plot( [b_arr[0],b_arr[-1]], [0,0], 'k')
    for k in range(1, len(b_arr)-1):
        plt.plot( [b_arr[k],b_arr[k]], [0,np.max(input_data)], 'k')
    
# This is currently the same as above except axis labels
elif 'w' in argv:
    """ Plot weights """
    ax.set_ylabel(r"$\eta(\mu)$")
    
    # Add input data
    ax.plot(mu_bins[:dom.I_cross], input_data[:dom.I_cross], 'bo-')
    ax.plot(mu_bins[dom.I_cross+1:], input_data[dom.I_cross+1:], 'go-')
    if dom.I_cross < len(mu_bins):
        ax.plot(mu_bins[dom.I_cross], input_data[dom.I_cross], 'ro')
    
    # Plot subdomain boundaries
    b_arr = np.array(boundaries) * B / float(Natoms)
    plt.plot( [b_arr[0],b_arr[-1]], [0,0], 'k')
    for k in range(1, len(b_arr)-1):
        plt.plot( [b_arr[k],b_arr[k]], [0,np.max(input_data)], 'k')


elif '2d' in argv:
    """ Plot a 2D matrix """
    ax.set_ylabel(r"$\mu \times \beta/N$")
    
    ax.pcolor(input_data, mu_bins, mu_bins, cmap='hot')
    #ax.imshow(mat, cmap='hot', interpolation='nearest')
    
    ax.colorbar()

elif 'log' in argv:
    """ log plot """
    
    plt.semilogy(mu_bins, input_data, 'ko')

print "Mu = E(%s) - E(%s)" %(alpha_type, beta_type)

plt.show()
