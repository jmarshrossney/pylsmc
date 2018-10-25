"""
Module.

Contains functions specific to the multicanonical algorithm.
"""

import numpy as np
from os.path import basename
import sys

# Import params from parent directory
path_to_params = '..'
sys.path.insert(0, path_to_params) # looks in parent directory first
from params import *
sys.path.remove(path_to_params)

# Import lsmc modules from this directory
import domain as dom
import initialise as ini

# Name of this file
this_file = basename(__file__)


# Number of steps to run for before returning to main.py and recording dF
Nsteps = Natoms * sweeps_dF

#################################
## File naming and i/o control ##
#################################
def file_names(stage, s=0, p=0):
    """ Retuns a dictionary of file names corresponding to each stage of the program """

    # Get the correct names for p, s
    sname, pname = ini.naming_system(s, p)

    if stage == 'input':
        # Input file names include appropriate '_pY_sX'
        disp_in = "disp" + pname[0] + sname[0] + ".out"
        weights_in = ini.rename_inputs(params_weights_file, sname[0], '')
        hist_in = "hist_MC" + pname[0] + sname[0] + ".out" # same as output
        return {'d': disp_in, 'w': weights_in, 'h': hist_in}
    
    elif stage == 'output':
        # Out file names include appropriate '_pY_sX'
        disp_out = "disp" + pname[0] + sname[0] + ".out"
        hist_out = "hist_MC" + pname[0] + sname[0] + ".out"
        series_out = "series_MC" + pname[0] + sname[0] + ".out"
        return {'d': disp_out, 'h': hist_out, 's': series_out}

    elif stage == 'pcomb':
        # After combined over processors
        hist_pcomb = "hist_MC" + pname[1] + sname[0] + ".out"
        series_pcomb = "series_MC" + pname[1] + sname[0] + ".out"
        return {'h': hist_pcomb, 's': series_pcomb}

    elif stage == 'unfld':
        # After unfolding
        hist_unfld = "hist_MC" + "_unfld" + sname[0] + ".out"
        series_unfld = "series_MC" + "_unfld" + sname[0] + ".out"
        return {'h': hist_unfld, 's': series_unfld}
    
    elif stage == 'scomb':
        # After subdomain join
        hist_scomb = "hist_MC" + "_unfld" + sname[1] + ".out"
        return {'h': hist_scomb}

    elif stage == 'dF':
        # dF estimate
        dF_file = "deltaF_MC.out"
        return dF_file
    
    else:
        error(this_file, "Invalid argument: "+stage)


def load_inputs(s, p):
    """ Load the correct input files for a multicanonical simulation """
   
    input_files = file_names('input', s, p)

    disp = ini.file_input(input_files['d'], (Natoms, 3))

    size = ini.get_size(s)
    weights = ini.file_input(input_files['w'], size)
    hist = ini.file_input(input_files['h'], size)

    return disp, {'w': weights, 'h': hist}


#######################################
## Functions called within main loop ##
#######################################
def refresh_func(binned_data):
    return


def update_func(step, binned_data, delta_local_energies,
            new_index, old_index, old_index_copy,
            old_mu, mu_bins,
            old_weight, get_weight, lendata, kTlogF, s):
    """ Update appropriate binned data after a single-particle move has been accepted or rejected """
  
    # Update histogram
    binned_data['h'][old_index] += 1
    
    # Second returned value is FLAT, for compatibility with Wang Landau simulation
    return old_weight, 1


def save_func(binned_data, mu_bins, step, s, p):
    """ Save appropriate binned data periodically during a simulation """
    print "Sweep %d.  Saving and continuing..." %(step/Natoms)

    savename = "_in_prog_p"+str(p)+"_s"+str(s)
    np.savetxt("hist"+savename+".out", binned_data['h'])

    return

