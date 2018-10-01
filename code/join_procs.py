"""
Script.

Combines data (histogram/transition matrix) from multiple simulations within the same subdomain, for all subdomains, saving to new files.

Called automatically from run_parallel.sh
"""

import numpy as np
from sys import argv
from os.path import basename

import energy
from params import *
import domain as dom
import initialise as ini

if algorithm == 'wang_landau': import wang_landau as alg
elif algorithm == 'multicanonical': import multicanonical as alg
elif algorithm == 'transition': import transition_matrix as alg

# Name of this file
this_file = basename('__file__')


if TRAP == True:
    if Np <= Ns:
        print "Doing nothing..."
        print ""
        exit(0)
    else:
        Ns_eff = Ns # combining multiple procs in multiple subdoms
else:
    Ns_eff = 1 # there may be more than 1 subdom, but each process is global

p_list = np.arange(Np)


# ----------------------------------- #
#  Iterate over different subdomains  #
# ----------------------------------- #
for s in range(Ns_eff):

    # Size of binned data for this subdomain
    size = ini.get_size(s)

    # Initialise arrays for combined data within this subdomain
    if argv[1] == 'h':
        data_comb = np.zeros(size)
    elif argv[1] == 'c':    
        data_comb = np.zeros( (size, size) )
    elif argv[1] == 's':
        data_comb == np.array([])


    # ------------------------------------------------------------ #
    #  Sum over processors operating in the same subdomain/domain  #
    # ------------------------------------------------------------ #
    for p in p_list[s::Ns_eff]:

        # Check correct mapping of processor to subdomain
        s_check = ini.map_proc_to_subdom(p)
        if s_check != s:
            error(this_file, "Mapping p -> s has gone wrong \nExpected p%d -> s%d, but got s%d" %(p,s_check,s))

        # Load files (output by main.py)
        input_file = alg.file_names('output', s, p)[argv[1]]
        this_data = np.loadtxt(input_file)
        
        # Add to sum
        if argv[1] in ('h', 'c'):
            data_comb = data_comb + this_data
        elif argv[1] == 's':
            new_run_indicator = [np.ones(5)*-1,] # Need to distinguish between independent runs
            data_comb = np.concatenate( (data_comb, new_run_indicator, this_data) )


    # End sum over processors
    # Continue for this subdomain


    # Save for each subdomain
    savename = alg.file_names('pcomb')[argv[1]]
    print "Saving combined data to ", savename
    np.savetxt(savename, data_comb)


    # --------------------------------------- #
    #  Additional steps if Transition Matrix  #
    # --------------------------------------- #
    if argv[1] == 'c':
        # Initialise eigvec computation with one of the other individual eigenvectors (just slightly faster than starting from random)
        eigvec_file = alg.file_names('output')['e']
        eigvec_init = np.loadtxt(eigvec_file)

        # Pull out a new eigenvector from the combined Pmat
        eigvec, TMweights, v = alg.get_Pmat_eigenvector(data_comb, eigvec_init)

        # Save eigvec
        savename = alg.file_names('pcomb')['e']
        print "Saving combined eigenvector to ", savename
        np.savetxt(savename, eigvec)

        # Save weights
        savename = alg.file_names('pcomb')['w']
        print "Saving combined weights to ", savename
        np.savetxt(savename, TMweights)

print ""
