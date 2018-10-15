"""
Script.

Combines data from multiple subdomains into a single, global dataset.

Also creates a number of additional datasets which can be used by calc_dF.py to compute the uncertainty in dF due to subdomain joining.

Run automatically by run_parallel.sh,
but can also be run manually using the following arguments:
    
    argv[1] - 'w', 'c', 'h', 'e' for weights, Cmat, histogram, eigvec data

    argv[2] - Optional. Filename STEM (in this case, everything before _sX.out). If none given uses automatic filenames based on params.py

Can also plot before/after join as a check - make sure to comment the plot requests out before running a long simulation!
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate, linalg
from pymbar import mbar
import copy

from sys import argv, exit
from os.path import basename

from params import *
import domain as dom
import initialise as ini

# Import module corresponding to the type of algorithm being used
if algorithm == 'wang_landau':
    import wang_landau as alg
    stage_in = 'output'
elif algorithm == 'multicanonical':
    import multicanonical as alg
    stage_in = 'unfld'
elif algorithm == 'transition':
    import transition_matrix as alg
    stage_in = 'pcomb'

# Name of this file
this_file = basename(__file__)


# Plotting parameters
plt.rc('text',usetex=True)
font = {'family' : 'serif',
        'size' : 16}
plt.rc('font', **font)
plt.rcParams['axes.linewidth'] = 3


# Compatibility with single-domain simulations
if Ns == 1 or TRAP == False:
    # Check that we're not trying to run this via the command line for non-default files
    if len(argv) == 2: # then no file stem was passed as a command line argument
        print "Doing nothing..."
        print ""
        exit(0)


# ------------------------------------------- #
#  Compile a list of input/output file names  #
# ------------------------------------------- #
# Check for sensible argv[1]
if argv[1] not in ('w','c','h','e'):
    error(this_file, "argv[1] should be 'w', 'c', 'h' or 'e'")

# Input file names for each subdomain
input_file_list = []
for s in range(Ns):

    if len(argv) > 2:   # file stem passed as argument
        input_file = argv[2] + "_s" + str(s) + ".out"
    else:               # default file names
        input_file = alg.file_names(stage_in, s)[argv[1]]
    
    input_file_list.append(input_file)

# Output file name
if len(argv) > 2: # use argument file stem as output file stem
    output_file = argv[2] + "_scomb.out"
else:
    output_file = alg.file_names('scomb')[argv[1]]


# --------------------------------------------- #
#  Optional: plot all subdomains in one figure  #
# --------------------------------------------- #
def plot_combined(data_list, mu_list):
    ''' Plot combined subdomain data, either joined up or not '''

    mkrs = ('b-','g-')
    mkrs2 = ('mo-','co-')
    
    ax1 = plt.axes()
    ax1.tick_params(direction='in',top=True,right=True)
    ax1.set_yticklabels([])
    ax1.set_xticklabels([])

    #fig, ax1 = plt.subplots()
    #ax1.set_title("<- fcc, bcc ->")
    #ax1.set_xlabel(r'$\mu \times \beta / N $')
    #ax2 = ax1.twinx()
    
    #if argv[1] == 'w':
        #ax1.set_ylabel(r'$\eta$')
        #ax2.set_ylabel(r'$\eta$')
    #elif argv[1] in ('h','e'):
        #ax1.set_ylabel(r"$P(\mu)$")
        #ax2.set_ylabel("P")

    for s in range(len(data_list)):
        ax1.plot(mu_list[s]*float(B)/Natoms,data_list[s],mkrs[s%2],linewidth=6,markersize=10)
        
        #if argv[1] in ('h','e'):
        #    ax2.semilogy(mu_list[s]*float(B)/Natoms,data_list[s],mkrs2[s%2])

    plt.tight_layout()
    plt.show()
    return


'''
# ----------------------------- #
#  Joining transition matrices  #
# ----------------------------- #
# Why is this commented out? Can't one can simply stitch Cmat's from different subdomains
# together (without overlaps!) and diagonalise to find the eigenvector for the entire domain?
# Possibly. I need to look into this more.
# However, one problem is that the method of updating the transition matrix at subdomain 
# boundaries doesn't seem particularly good, so an overlap is currently necessary in order
# for these edge bins to be discarded.
# Even so, we shouldn't need an overlap as big as the required by MBAR for stitching.
# Actually, the main reason it's commented out is because it's not working. The eigenvector
# of the joined matrix has discontinuities at subdomain boundaries, even when the bins are
# the same thickness. 

if argv[1] == 'c':

    Cmat_list = []
    for s in range(Ns):
        Cmat_s = np.loadtxt(input_file_list[s])
        lo = dom.subdom[s]['l_bin_olap']
        hi = -dom.subdom[s]['h_bin_olap']
        if s == Ns-1: hi = None
        
        Cmat_interior = Cmat_s[lo:hi,lo:hi]
        Cmat_list.append(Cmat_interior)
        print "Subdomain %d, Cmat shape: %s -> %s" \
                %(s,np.shape(Cmat_s),np.shape(Cmat_interior))

    # Combine into one banded matrix
    Cmat_comb = linalg.block_diag(*Cmat_list)
    print "Combined Cmat shape: ",np.shape(Cmat_comb)
    print "Saving combined Cmat to file: ", output_file
    np.savetxt(output_file,Cmat_comb)
    
    # If combined eigvec already exists, initialise with that
    output_files = alg.file_names('scomb')
    eigvec_file = output_files['e']
    if exists(eigvec_file):
        eigvec_init = np.loadtxt(eigvec_file)
    else:
        eigvec_init = np.random.random(len(Cmat_comb))

    # Pull out eigenvector from the combined matrix
    eigvec, TMweights, v = alg.get_Pmat_eigenvector(Cmat_comb, eigvec_init)

    # Possibly need to scale by bin width here, as in calc_dF.py

    # Save eigvec, no need to re-evaluate naming system
    print "Saving combined eigenvector to ", eigvec_file
    np.savetxt(eigvec_file, eigvec)

    # Save weights
    weights_file = output_files['w']
    print "Saving combined weights to ", weights_file
    np.savetxt(weights_file, TMweights)

    print ""
    exit(0)

'''

# ----------------------------------------- #
#  Iterate over subdomains and import data  #
# ----------------------------------------- #
# We're going to store arrays as elements in lists
data_list = []
data_list_err = []
for j in range(Nsamples_join):
    data_list_err.append([]) # Necessary to keep copies of the data independent
mu_list = []
                   
# Import input data from different subdomains
for s in range(Ns):
    
    input_data = np.loadtxt(input_file_list[s])

    # Take log if histogram/probability data
    if argv[1] in ('h','e'):
        input_data = -np.log(input_data)

    data_list.append(input_data)

    # Make independent copies of the data for future error estimation
    for j in range(Nsamples_join):
        input_data_copy = np.copy(input_data) # copy to independent array
        data_list_err[j].append(input_data_copy)
                          
    # Generate array for bin mu positions
    mu_s = dom.get_local_mu_bins(s)
    mu_list.append(mu_s)

# End loop over subdomains


# !todo! Should have a check: if data_list contains zeros, warn

# Plot original
#plot_combined(data_list, mu_list)


# -------------------------------------------------- #
#  Iterate over subdomains and stitch together data  #
# -------------------------------------------------- #
data_min = np.min(data_list[0])
data_min_err = [data_min,]*Nsamples_join

for s in range(Ns-1):

    # Min/max mu of overlap region, then converted to bin index from perspective of the 'other' subdomain
    olap_min = dom.subdom[s+1]['min']
    olap_max = dom.subdom[s]['max']
    olap_min_bin = dom.get_local_index(olap_min,s)
    olap_max_bin = dom.get_local_index(olap_max,s+1) 

    # Interpolate, within overlap region, datasets from adjacent subdomains
    fs = interpolate.interp1d( mu_list[s][olap_min_bin:], data_list[s][olap_min_bin:],
            kind='linear')
    fspl1 = interpolate.interp1d( mu_list[s+1][:olap_max_bin+1], data_list[s+1][:olap_max_bin+1],
            kind='linear')
    
    # Pick a discrete set of mu values in overlap region, more fine than the bin mu's
    interp_min = mu_list[s][olap_min_bin+3] # the 3 is to make sure all points are within the boundaries of fs, fspl1
    interp_max = mu_list[s+1][olap_max_bin-3]
    mu_interp = np.linspace(interp_min,interp_max,Ninterp)
    
    # Find the index of this discrete set corresponding to the actual boundary between subdomains (used by MBAR)
    interp_i_cross = Ninterp # in case all points are left of boundary!
    for i in range(len(mu_interp)-1):
        if mu_interp[i] < boundaries[s+1] and mu_interp[i+1] > boundaries[s+1]:
            interp_i_cross = i

    # Evaluate the interpolation functions at this set of points
    ls = fs(mu_interp)
    lspl1 = fspl1(mu_interp)

    # -------------------------------------------------- #
    #  Use MBAR to shift subdomains by calculating free  #
    #  energy difference between data in overlap region  #
    # -------------------------------------------------- #
    print "Joining subdomains %d and %d using MBAR" %(s,s+1)
    
    # See pymbar documentation for detailed explanation
    U_kn = np.array( [ls, lspl1] )
    N_k = [interp_i_cross, Ninterp - interp_i_cross]
    mbar_obj = mbar.MBAR(U_kn, N_k)
    matrix = mbar_obj.getFreeEnergyDifferences(compute_uncertainty=True)

    # Shift and error, output by pymbar
    shift = matrix[0][1,0]
    sigma = matrix[1][1,0]
    print "shift: ", shift, "+/-", sigma

    # Minimise free energy difference between subdomains by shifting the higher-mu subdomain
    data_list[s+1] += shift
        
    # Sample shift error from a mean 0 Gaussian with standard deviation 'sigma'
    for j in range(Nsamples_join):
        uncertainty = np.random.normal(loc=0.0, scale=sigma)
        data_list_err[j][s+1] += (shift + uncertainty) # Data copies are each shifted by a different amount
    
    # Update minimum values after shift
    data_min = min(data_min, np.min(data_list[s+1]))
    for j in range(Nsamples_join):
        data_min_err[j] = min(data_min_err[j], np.min(data_list_err[j][s+1]))


# End loop over subdomains


# ---------------------- #
#  Tidy up and undo log  #
# ---------------------- #
# Weights
if argv[1] == 'w':
    # Shift so that the minimum is at 0
    for s in range(Ns):
        data_list[s][:] -= data_min
        for j in range(Nsamples_join):
            data_list_err[j][s][:] -= data_min_err[j]

# Histogram or eigvec
if argv[1] in ('h','e'):
    # Sum data over all subdomains for later normalisation
    Psum = 0
    Psum_err = [0,]*Nsamples_join
    
    # Shift so minimum as at 0, then undo log for histogram/probability data
    for s in range(Ns):
        data_list[s][:] -= data_min
        data_list[s][:] = np.exp(-data_list[s][:])
        Psum += np.sum(data_list[s][:])
        
        for j in range(Nsamples_join):
            data_list_err[j][s][:] -= data_min_err[j]
            data_list_err[j][s][:] = np.exp(-data_list_err[j][s][:])
            Psum_err[j] += np.sum(data_list_err[j][s][:])

    # Normalise the entire thing
    for s in range(Ns):
        data_list[s][:] /= Psum
        for j in range(Nsamples_join):
            data_list_err[j][s][:] /= Psum_err[j]

# Plot joined up
#plot_combined(data_list, mu_list)


# ---------------------------- #
#  Combine data into one list  #
# ---------------------------- #
# First create one big list for the data
data_combined = []
data_combined_err = []
for j in range(Nsamples_join): data_combined_err.append([])

# 
for s in range(Ns):

    # Ignore overlap bins
    lo = dom.subdom[s]['l_bin_olap']
    hi = dom.subdom[s]['bins'] - dom.subdom[s]['h_bin_olap']
    data_interior = list(data_list[s][lo:hi])

    # Append to big list
    data_combined = data_combined + data_interior

    for j in range(Nsamples_join):
        data_interior_err = list(data_list_err[j][s][lo:hi])
        data_combined_err[j] = data_combined_err[j] + data_interior_err

# Convert to arrays
data_combined = np.array(data_combined)
data_combined_err = np.array(data_combined_err)


# ------------------ #
#  Save joined data  #
# ------------------ #
# Weights
if argv[1] == 'w':

    # Minimum to 0 once again
    data_combined -= np.min(data_combined)
    for j in range(Nsamples_join):
        data_combined_err[j,:] -= np.min(data_combined_err[j,:])
    
    print "Saving combined weights to file: ",output_file
    np.savetxt(output_file,data_combined)

# Histogram
if argv[1] == 'h':

    # Normalise once again          !todo! Is this necessary?
    data_combined = data_combined / np.sum(data_combined)
    for j in range(Nsamples_join):
        data_combined_err[j,:] = data_combined_err[j,:] / np.sum(data_combined_err[j,:])
    
    print "Saving combined histogram to file: ",output_file
    np.savetxt(output_file,data_combined)

# Eigvec
if argv[1] == 'e':

    # Normalise once again
    data_combined = data_combined / np.sum(data_combined)
    for j in range(Nsamples_join):
        data_combined_err[j,:] = data_combined_err[j,:] / np.sum(data_combined_err[j,:])
    
    print "Saving combined Probability estimate to file: ",output_file
    np.savetxt(output_file,data_combined)


# Save the other joined data to an error file
err_output_file = output_file.replace('.out','.err')
print "Saving %d joins with errors to: %s" %(Nsamples_join,err_output_file)
np.savetxt(err_output_file,data_combined_err)

print ""
