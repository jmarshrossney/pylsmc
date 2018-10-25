"""
Script.

Combines data (histogram/transition matrix) from multiple simulations within the same subdomain, for all subdomains, saving to new files.

Called automatically from run_parallel.sh
"""

import numpy as np
from sys import argv
from os.path import basename

from params import *
import domain as dom
import initialise as ini

if algorithm == 'wang_landau': import wang_landau as alg
elif algorithm == 'multicanonical': import multicanonical as alg
elif algorithm == 'transition': import transition_matrix as alg

# Name of this file
this_file = basename(__file__)


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
        data_comb = []


    # ------------------------------------------------------------ #
    #  Sum over processors operating in the same subdomain/domain  #
    # ------------------------------------------------------------ #
    for p in p_list[s::Ns_eff]:

        if TRAP == True:
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
            ldc = len(data_comb)
            tmp = np.zeros( (ldc+len(this_data), 4) )
            if ldc > 0:
                tmp[:ldc,:] = data_comb
            tmp[ldc:,:] = this_data
            data_comb = tmp


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
        # Pull out a new eigenvector from the combined Pmat
        TMweights = alg.get_Pmat_eigenvector(data_comb)

        # Save weights
        savename = alg.file_names('pcomb')['w']
        print "Saving combined weights to ", savename
        np.savetxt(savename, TMweights)

print ""
