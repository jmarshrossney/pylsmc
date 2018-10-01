"""
Script.

Computes free energy difference (dF) from an unweighted probability distribution, by partitioning distribution at mu=0 and summing each side.
Also computes errors from subdomain joining by repeating the calculation for each probability distribution in the error file (if it exists), and taking the standard error of the dF values.

Called automatically in run_parallel.sh.
Possible to run manually as   python calc_dF.py (input_file) (error_file)

"""

import numpy as np
import math as m
from sys import argv

from params import *
import domain as dom
import initialise as ini

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
    key = 'e'


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


# --------------------------------- #
#  Import probability distribution  #
# --------------------------------- #
# Automatic file names
if len(argv) == 1:
   
    input_file = alg.file_names('scomb')[key]

    if Ns > 1 and TRAP == True:
        err_file = input_file.replace('.out','.err')
    else:
        err_file = None

# File name given as argv[]
if len(argv) > 1: 
    input_file = argv[1]

    if len(argv) > 2:
        err_file = argv[2]
    else:
        err_file = None

print "Input file: ", input_file
input_data = np.loadtxt(input_file)


# ---------------------------------------------------- #
#  Do calculation and errors due to subdomain joining  #
# ---------------------------------------------------- #
# Calculate dF for main input data
deltaF = do_calc(input_data)

if err_file != None:
    # Calculate join errors for dF

    err_data = np.loadtxt(err_file)
    print "Join errors: counted %d datasets" %len(err_data[:,0])
    
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


# Only write to file if called during simulation, rather than with singular use
if len(argv) == 1:

    print "Writing to file: ", deltaF_file
    with open(deltaF_file, 'a+') as f_handle:
        f_handle.write(str(deltaF)+'    '+str(err_est)+'\n')

print ""
print "============================================================================="
print ""

