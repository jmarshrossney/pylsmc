"""
Main Script.

Calls functions in initialise.py to set up lattices and naming system.
Calls functions in lattice_switch.py to perform simulation.
"""

import numpy as np
import math as m
import os.path
from sys import argv, exit

from params import *
import domain as dom
import initialise as ini
import lattice_switch as ls

if track_dynamics == True:
    import dynamics as dyn

if algorithm == 'wang_landau':
    import wang_landau as alg
elif algorithm == 'multicanonical':
    import multicanonical as alg
    d = 'h' # key for histogram data
elif algorithm == 'transition':
    import transition_matrix as alg
    d = 'c' # key for transition matrix data


# ------------------------------------------------------------------------ #
#               The process number/subdomain index for this simulation     #
#                                                                          #
if len(argv) > 1: p = int(argv[1])                                         #
else: p = 0                                                                #
print "Starting process no. ", p                                           #
s = ini.map_proc_to_subdom(p)                                              #
#                                                                          #
#                                                                          #
# ------------------------------------------------------------------------ #


# ------------------ #
#  Load input files  #
# ------------------ #

disp, binned_data = alg.load_inputs(s, p)

if track_dynamics == True:
    dyn_data = dyn.load_inputs(s, p)
else:
    dyn_data = None

# ------------------------- #
#  Set up the two lattices  #
# ------------------------- #
# Shuffle atomic indices of one if both lattice types are the same
if alpha_type == beta_type: shuffle = True
else: shuffle = False

# Supercells - these are ASE Atoms objects
atoms_alpha = ini.build_supercell(disp, alpha_vec, alpha_a, alpha_type)
atoms_beta  = ini.build_supercell(disp, beta_vec, beta_a, beta_type, shuffle)


# ------------------------ #
#  Drift to subdomain 's'  #
# ------------------------ #
if disp.any() == False: # if starting from ideal lattice positions
    # Drift to correct subdomain and save configurations
    print "Drifting to subdomain ", s
    mu = ls.drift(atoms_alpha, atoms_beta, disp, s)
    print "Final lattice configurations give mu = ", mu

    # Separate independent runs by saving this line to the series output file
    if track_series == True:
        series_file = alg.file_names('output', s, p)['s']
        series = open(series_file, 'a')
        series.write("-1 -1 -1 -1 \n")
        series.close()



###################################################################################
##                                                                               ##
##                   Run a Lattice-switch Monte Carlo simulation                 ##
##                      by calling lattice_switch.run()                          ##
##                                                                               ##
###################################################################################


if algorithm == 'wang_landau':
    print "Please use the regular main.py file, which will work just fine for batch jobs"

########################################
## Multicanonical / Transition Matrix ##
########################################
if algorithm in ('multicanonical', 'transition'):

    print "Running for %d iterations of %d sweeps" %(interations, sweeps_dF)
    
    # Names of files to which data is to be saved after each iteration
    output_files = alg.file_names('output', s, p)
    
    # Add entry to output_files for the list containing all iters
    output_files['allIt'] = "allIt_" + output_files[d]
    
    # Initialise list to hold the results (hist or Pseq) of each iteration
    data_all_iters = []
    
    # Loop over iterations (without resetting atom positions in between)
    for ITER in range(iterations):

        print "Iteration no.", ITER

        # Run the Lattice-switch algorithm for a given number of steps
        steps = ls.run(
                    atoms_alpha, atoms_beta,
                    disp,
                    binned_data,
                    dyn_data,
                    F=1, p=p, s=s)
        
        # Append important data (hist or Cmat) from this iteration to list
        data_all_iters.append(binned_data[d])

        # Add entry to binned_data so that data from the iterations completed
        # so far is saved
        binned_data['allIt'] = np.array(data_all_iters)

        # -------------------- #
        #  Save (& overwrite)  # 
        # -------------------- #
        # Add displacement array to binned_data so it gets saved
        binned_data['d'] = disp
        
        # (For each iteration, in case simulation doesn't complete)
        # Keys indicate type of data ('w', 'h' etc.)
        for key in output_files.keys():
            if key != 's' and output_files[key] != None:
                print "Saving to file: ", output_files[key]
                np.savetxt(output_files[key], binned_data[key])

        if track_dynamics == True:
            dyn_output_files = dyn.file_names('output', s, p)
            for key in dyn_output_files.keys():
                if key not in ('led', 'red'):
                    print "Saving to file: ", dyn_output_files[key]
                    np.savetxt(dyn_output_files[key], dyn_data[key])
    
        print ""
    
    print "Finished all iterations."

print ""

