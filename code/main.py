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
    import lsmc_dynamics as dyn

if algorithm == 'wang_landau':
    import wang_landau as alg
elif algorithm == 'multicanonical':
    import multicanonical as alg
elif algorithm == 'transition':
    import transition_matrix as alg


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
atoms_alpha = ini.build_supercell(disp, alpha_vec, alpha_type)
atoms_beta  = ini.build_supercell(disp, beta_vec, beta_type, shuffle)


# ------------------------ #
#  Drift to subdomain 's'  #
# ------------------------ #
if disp.any() == False: # if starting from ideal lattice positions
    # Drift to correct subdomain and save configuration
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


###################
##  Wang-Landau  ##
###################
if algorithm == 'wang_landau':

    # Size of binned data for this subdomain
    size = ini.get_size(s)
    
    print ""
    print "Running with Boltzmann weights for %d sweeps..." %sweeps_relax
        
    # Allow the model to 'relax' away from its initial configuration
    # Note that disp vector needs to be passed onto the next stage
    steps = ls.run(
                atoms_alpha, atoms_beta,
                disp, 
                {'w': np.zeros(size), # weights
                'h': np.zeros(size) }, # hist
                dyn_data,
                F=1, p=p, s=s)

    # ----------------------------------------------------- #
    #  Iterate over incremental factors of decreasing size  #
    # ----------------------------------------------------- #
    F = F_init
    while F > F_min:
        print ""
        print "Wang factor = ", F   # Not technically the correct name
        
        print "Building weights until flatness achieved..."
        binned_data['h'][:] = 0 # reset hist for each F
        
        # Run the Lattice-switch algorithm with Wang-Landau until the histogram is flat
        steps = ls.run(
                    atoms_alpha, atoms_beta,
                    disp,
                    binned_data,
                    dyn_data,
                    F=F, p=p, s=s)

        # Save for this F
        alg.save_F(F, steps, binned_data, s)

        # Update F
        F = np.sqrt(F)

    print ""


########################################
## Multicanonical / Transition Matrix ##
########################################
if algorithm in ('multicanonical', 'transition'):

    print "Running for %d sweeps" %sweeps_dF
    
    # Run the Lattice-switch algorithm for a given number of steps
    steps = ls.run(
                atoms_alpha, atoms_beta,
                disp, 
                binned_data,
                dyn_data,
                F=1, p=p, s=s)


# ------------------------------------------------------ #
#  Save outputs (for Wang Landau and MC/TM simulations)  #
# ------------------------------------------------------ #
# Names of files to which data is to be saved
output_files = alg.file_names('output', s, p)

# Add displacement array to binned_data so it gets saved
binned_data['d'] = disp

# Keys indicate type of data ('w', 'h' etc.)
for key in output_files.keys():
    if key != 's' and output_files[key] != None:
        print "Saving to file: ", output_files[key]
        np.savetxt(output_files[key], binned_data[key])

# Save the data on the simulation dynamics
if track_dynamics == True:
    dyn_output_files = dyn.file_names('output', s, p)
    for key in dyn_output_files.keys():
        if key not in ('led', 'red'):
            print "Saving to file: ", dyn_output_files[key]
            np.savetxt(dyn_output_files[key], dyn_data[key])


print ""

