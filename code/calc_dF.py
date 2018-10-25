"""
Script.

Computes free energy difference (dF) from an unweighted probability distribution, by partitioning distribution at mu=0 and summing each side.
Also computes errors from subdomain joining by repeating the calculation for each probability distribution in the error file (if it exists), and taking the standard error of the dF values.

Called automatically in run_parallel.sh.
Possible to run manually as   python calc_dF.py (input_file) (error_file)

"""

import numpy as np
import math as m
from os.path import exists, basename
from sys import argv

from params import *
import domain as dom
import initialise as ini

# Name of this file
this_file = basename(__file__)

if algorithm == 'multicanonical':
    import multicanonical as alg
    if Ns > 1 and TRAP == True:
        scale_bins = True
    else:
        scale_bins = False
    key = 'h'
if algorithm == 'transition':
    import transition_matrix as alg
    scale_bins = True # Always need to scale bins by their width
    key = 'w'
if algorithm == 'wang_landau':
    import wang_landau as alg
    scale_bins = True
    key = 'w'


def do_calc(prob_dist):
    """ Calculates free energy difference by partitioning probability
        distribution at mu = 0 and applying 
            dF = -k T log(P_alpha/P_beta) """    
  
    if scale_bins == True:
        # Scale bins by their width
        bin_low = 0
        bin_high = 0
        for s in range(Ns):
            bin_high += bins[s]
            prob_dist[bin_low:bin_high] *= dom.subdom[s]['bin_width']
            bin_low = bin_high

    # Sum bins below mu=0
    alpha_sum = np.sum(prob_dist[:dom.I_cross])

    # Sum bins above mu=0
    beta_sum = np.sum(prob_dist[dom.I_cross+1:])

    # Bin centered on mu=0
    alpha_sum += 0.5*prob_dist[dom.I_cross]
    beta_sum += 0.5*prob_dist[dom.I_cross]
    
    # Compute free energy difference
    deltaF = -kT * np.log(float(alpha_sum)/beta_sum)

    # Account for adjusted zero-point energy
    deltaF += adjust

    return deltaF


# ----------------------------- #
#  Input and output file names  #
# ----------------------------- #
sweeps = None

if len(argv) == 1:
    # Called by run.sh automaticall - default file names
    input_file = alg.file_names('scomb')[key]
    output_file = alg.file_names('dF')

elif len(argv) > 1: 
    # File names given as arguments
    input_file = argv[1]

    # Output file can be given as argument
    # If none, just print to stdout
    if len(argv) > 2:
        output_file = argv[2]
    else:
        output_file = None

    # Sweeps for this iter can be given as argument
    # Only useful for Wang-Landau
    if len(argv) > 3:
        sweeps = float(argv[3])

# Error file has the same name as input file
err_file = input_file.replace('.out','.err')


# ---------------------------------------------------- #
#  Do calculation and errors due to subdomain joining  #
# ---------------------------------------------------- #
print "Input file: ", input_file
input_data = np.loadtxt(input_file)

if key == 'w': # convert to estimate for probability dist
    input_data -= np.min(input_data)
    input_data = input_data * B
    input_data = np.exp(input_data)

# Calculate dF for main input data
deltaF = do_calc(input_data)

if os.path.exists(err_file) == True:
    # Calculate join errors for dF

    err_data = np.loadtxt(err_file)
    print "Join errors: counted %d datasets" %len(err_data[:,0])
    
    if key == 'w': # convert to estimate for probability dist
        err_data -= np.amin(err_data, axis=1)
        err_data = err_data * B
        err_data = np.exp(err_data)
    
    # Calculate dF for each distribution/dataset in err_data
    deltaF_err = np.zeros(Nsamples_join)
    for j in range(Nsamples_join):
        deltaF_err[j] = do_calc(err_data[j,:])

    # Calculate standard error on the dF's
    sigma = np.std(deltaF_err)
    err_est = sigma / Nsamples_join

else:
    # No join errors
    err_est = 0
    print "Join errors: None"


# -------------------------------- #
#  Print result and write to file  #
# -------------------------------- #
print "Free energy difference (alpha-beta): %1.2e +/- %1.2e" %(deltaF, err_est)

# In scaled, dimensionless units
deltaF = deltaF * B/Natoms
err_est = err_est * B/Natoms
print "In units of kT per molecule: %1.2e +/- %1.2e" %(deltaF, err_est)


# Write to output file if called by run.sh, or if given as argument
if output_file != None:

    # How many sweeps did this iteration take?
    if sweeps == None:
        sweeps = alg.Nsteps / Natoms # True for MC, TM simulations, but not WL!

    print "Writing to file: ", output_file
    with open(output_file, 'a+') as f:
        f.write(str(sweeps)+'   '+str(deltaF)+'   '+str(err_est)+'\n')

print ""
print "============================================================================="
print ""

