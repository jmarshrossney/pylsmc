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
p = int(argv[1])                                                           #
print "Starting process no. ", p                                           #
s = ini.map_proc_to_subdom(p)                                              #
#                                                                          #
#                                                                          #
# ------------------------------------------------------------------------ #


# ------------------------- #
#  Set up the two lattices  #
# ------------------------- #
# Supercells - these are ASE Atoms objects
atoms_alpha = ini.build_supercell(alpha_vec, alpha_a, alpha_type)
atoms_beta  = ini.build_supercell(beta_vec, beta_a, beta_type)

# If two of the same lattice, shuffle positions
if alpha_type == beta_type:
    pos = atoms_beta.positions
    np.random.shuffle(pos)
    atoms_beta.set_positions(pos)

# Save the ideal lattice positions
ideal_positions_alpha = atoms_alpha.get_positions()
ideal_positions_beta = atoms_beta.get_positions()


# ---------------------------------------- #
#  Find an appropriate range of mu values  #
# ---------------------------------------- #
if algorithm == 'domain':

    print "Running for %d sweeps..." %sweeps_dF

    output = ls.find_mu_limits(atoms_alpha, atoms_beta, s)

    savename = "mu_range_s" + str(s) + ".out"
    print "Saving to file: ", savename
    np.savetxt(savename, output)

    # Exit happily
    exit(0)



###################################################################################
##                                                                               ##
##                   Run a Lattice-switch Monte Carlo simulation                 ##
##                      by calling lattice_switch.run()                          ##
##                                                                               ##
###################################################################################

# ------------------------ #
#  Drift to subdomain 's'  #
# ------------------------ #
# Initial positions are ideal
disp = np.zeros((Natoms,3))

# Drift to correct subdomain and save configuration
disp_init = ls.drift(atoms_alpha, atoms_beta, disp, s)
init_positions_alpha = atoms_alpha.get_positions()
init_positions_beta = atoms_beta.get_positions()


# ------------------ #
#  Load input files  #
# ------------------ #
binned_data = alg.load_inputs(s, p)


###################
##  Wang-Landau  ##
###################
if algorithm == 'wang_landau':

    # Size of binned data for this subdomain
    size = ini.get_size(s)

    print "Powering up the Wang Machine..."
    F = F_init

    # ----------------------------------------------------- #
    #  Iterate over incremental factors of decreasing size  #
    # ----------------------------------------------------- #
    while F > F_min:
        print ""
        print "Wang factor = ", F   # Not technically the correct name
        
        # Reset atoms to initial positions in correct subdom (not necessary, I suspect)
        disp = disp_init
        atoms_alpha.set_positions(init_positions_alpha)
        atoms_beta.set_positions(init_positions_beta)

        print "Running with Boltzmann weights for %d sweeps..." %sweeps_relax
        
        # Allow the model to relax to an equilibrium
        # Note that disp vector needs to be passed onto the next stage
        disp = ls.run(
                atoms_alpha, atoms_beta, disp,
                {'w': np.zeros(size), # weights
                 'h': np.zeros(size) }, # hist
                f=1, p=p, s=s)[0]
        
        print "Building weights until flatness achieved..."
        binned_data['h'][:] = 0 # reset hist for each F
        
        # Run the Lattice-switch algorithm with Wang-Landau until the histogram is flat
        steps = ls.run(
                atoms_alpha, atoms_beta, disp,
                binned_data,
                f=F, p=p, s=s)[1]

        # Save for this F
        alg.save_F(F, steps, binned_data, s)

        # Update F
        F = np.sqrt(F)

    print ""

    print "Wang has concluded"
    print ""

    # ------------------------------ #
    #  Save outputs for Wang Landau  #
    # ------------------------------ #
    # Names of files to which data is to be saved
    output_files = alg.file_names('output', s, p)

    # Keys indicate type of data ('w', 'h' etc.)
    for key in output_files.keys():
        if key != 's' and output_files[key] != None:
            print "Saving to file: ", output_files[key]
            np.savetxt(output_files[key], binned_data[key])


########################################
## Multicanonical / Transition Matrix ##
########################################
if algorithm == 'multicanonical':

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
        disp = ls.run(
            atoms_alpha, atoms_beta, disp,
            binned_data,
            f=1, p=p, s=s)[0]
        
        # Append important data (hist or Cmat) from this iteration to list
        data_all_iters.append(binned_data[d])

        # Add entry to binned_data so that data from the iterations completed
        # so far is saved
        binned_data['allIt'] = np.array(data_all_iters)

        # -------------------- #
        #  Save (& overwrite)  # 
        # -------------------- #
        # (For each iteration, in case simulation doesn't complete)
        # Keys indicate type of data ('w', 'h' etc.)
        for key in output_files.keys():
            if key != 's' and output_files[key] != None:
                print "Saving to file: ", output_files[key]
                np.savetxt(output_files[key], binned_data[key])
    
        print ""
    
    print "Finished all iterations."

print ""

