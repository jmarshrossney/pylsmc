"""
Script.

Takes the raw outputs from a batch-processed simulation (where each iteration of 'sweeps' sweeps
has been saved, but the post-processing not carried out), and selects a particular iteration
of the data. This 'snapshot' is saved, to be processed by the regular post-processing functions,
which are executed by   post_proc_batch.sh

"""

import numpy as np
from sys import argv, exit

import initialise as ini
from params import *

elif algorithm == 'multicanonical':
    import multicanonical as alg
    d = 'h' # key for histogram data
elif algorithm == 'transition':
    import transition_matrix as alg
    d = 'c' # key for transition matrix data

# Start counting from 0
this_iter = int(argv[1])

# Keep a tally of the number of files skipped due to not enough iterations
skipped = 0

for p in range(Np):
    
    # Get subdomain and file names
    s = ini.map_proc_to_subdom(p)
    file_names = alg.file_names('output', s, p)
    
    # Input file name
    input_file = "allIt" + file_names[d]
    
    # Import histogram data from all iterations
    data_all_iters = np.loadtxt(input_file)

    # Check that enough iterations have actually been completed?
    if this_iter >= len(data_all_iters):
        print "Not enough iterations in file: %s. Skipping..." %input_file
        
        # Exit if none of the simulations have completed this many iterations
        skipped += 1
        if skipped == Np:
            print "All files skipped. Exiting..."
            exit(1)
        else:
            continue

    # Pull out a 'snapshot' of just this iteration
    snapshot = data_all_iters[this_iter]

    # Save to regular output file name
    output_file = file_names[d]
    np.savetxt(output_file, snapshot)


