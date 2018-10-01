"""
Script.

'Unfolds' weights to recover physical probability distribution from a histogram generated using a fixed set of weights.

How to use:
        argv[1] specifies the format of the to be unfolded - h (histogram) or s (series).
    
    The files are the automatic ones specified in params/initialise.py. Currently there is
    no option to unfold a user-defined file with user-defined weights.

    e.g.
        python unfold_weights.py h
        python unfold_weights.py s

Called automatically from run_parallel.sh to unfold histogram data, after joining 
within subdomain but before joining subdomains.
"""

import math as m
import numpy as np
from sys import argv

from params import *
import domain as dom
import initialise as ini
import multicanonical as alg

if TRAP == True:
    Ns_eff = Ns # multiple procs on same subdomain already combined in join_procs
else:
    Ns_eff = 1 # there may be more than 1 subdom, but each process is global - simulation has access to entire domain

# ----------------------------------- #
#  Iterate over different subdomains  #
# ----------------------------------- #
for s in range(Ns_eff):

    # Files to load
    weights_in = alg.file_names('input', s)['w']
    hist_pcomb = alg.file_names('input', s)['h']

    # Import weights
    weights = np.loadtxt(weights_in)
    
    # Import histogram
    hist = np.loadtxt(hist_pcomb)
        
    # in hist: 1 -> exp(B*weight)
    factor = np.exp(B*weights)
    hist_unfld = hist*factor
   
    # Save unfolded histogram
    save_name = alg.file_names('unfld', s)['h']
    print "Saving unfolded histogram to ", save_name
    np.savetxt(save_name, hist_unfld)


print ""
