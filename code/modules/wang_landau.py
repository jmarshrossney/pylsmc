"""
Module.

Contains functions specific to the Wang Landau algorithm.
"""

import numpy as np
from os.path import exists, basename
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


# Number of steps to take before only focusing on flatness
Nsteps = Natoms * sweeps_relax

def extrapolate(array):
    """ Adds linearly extrapolated end bins to a numpy array, assuming the
        array has spaces at either end for this extrapolated value to go """
    array[-1] = 2*array[0] - array[1] # this is really to the left of bin 0
    array[-2] = 2*array[-3] - array[-4]
    return

def check_flat(hist):
    """ Check for histogram 'flatness': flat if all bins > flat_tol x mean """
    mean = np.mean(hist)
    if len(np.where(hist < flat_tol*mean)[0]) == 0:
        return 3
    else:
        return 1


#################################
## File naming and i/o control ##
#################################
def file_names(stage, s=0, p=0):
    """ Retuns a dictionary of file names corresponding to each stage of the program """

    # Get the correct names for p, s
    sname, pname = ini.naming_system(s, p)
    
    if stage == 'input':
        # Input file names include appropriate '_pY_sX'
        disp_in = "disp" + sname[0] + ".out"
        weights_in = ini.rename_inputs(params_weights_file, sname[0], '')
        return {'d': disp_in, 'w': weights_in}

    elif stage == 'output':
        # Output file names include appropriate '_pY_sX'
        disp_out = "disp" + sname[0] + ".out"
        weights_out = "weights_WL" + sname[0] + ".out" # probably don't need pname
        hist_out = "hist_WL" + sname[0] + ".out"
        series_out = "series_WL" + sname[0] + ".out"
        return {'d': disp_out, 'w': weights_out, 'h': hist_out, 's': series_out}

    elif stage == 'scomb':
        # After subdomain join
        weights_scomb = "weights_WL" + sname[1] + ".out"
        hist_scomb = "hist_WL" + sname[1] + ".out"
        return {'w': weights_scomb, 'h': hist_scomb}

    else:
        error(this_file, "Invalid argument: "+stage)


def load_inputs(s, p):
    """ Load the correct input files for a Wang Landau simulation """

    input_files = file_names('input', s, p)
    
    disp = ini.file_input(input_files['d'], (Natoms, 3))
    
    size = ini.get_size(s)
    weights = ini.file_input(input_files['w'], size)
    hist = ini.file_input(None, size)

    return disp, {'w': weights, 'h': hist}
    

############################################
## Functions called from within main loop ##
############################################
def refresh_func(binned_data):
    return


def update_func(step, binned_data, delta_local_energies, 
            new_index, old_index, old_index_copy,
            old_mu, mu_bins,
            old_weight, get_weight, lendata, kTlogF, s): 
    """ Update appropriate binned data after a single-particle move has been accepted or rejected """

    # Update histogram and weights
    binned_data['w'][old_index] += kTlogF 
    binned_data['h'][old_index] += 1 

    # Need to update old weight after incrementation
    old_weight = get_weight(old_mu, old_index, mu_bins, binned_data['w'])
    if use_interpolated_weights == True:
        # Must extrapolate each time an edge bin is updated
        if old_index == 0 or old_index == lendata-1:
            extrapolate(binned_data['w'])
   
    # Check for flatness every sweep
    if step % Natoms == 0:
        FLAT = check_flat(binned_data['h'])
    else:
        FLAT = 1

    return old_weight, FLAT


def save_func(binned_data, mu_bins, step, s, p):
    """ Save appropriate binned data periodically during a simulation """

    print "Sweep %d.  Saving and continuing..." %(step/Natoms)

    # If interpolated weights used, remove those extrapolated bins / 'ghost points'
    if use_interpolated_weights == True:
        binned_data['w'] = binned_data['w'][:-2]

    savename = "_in_prog_s"+str(s)
    np.savetxt("weights"+savename+".out", binned_data['w'])
    np.savetxt("hist"+savename+".out", binned_data['h'])

    return


#################################
## Save after each F iteration ##
#################################
def save_F(F, counters, binned_data, s):

    # Save to F, sweeps file
    steps = counters[0]
    Nsweeps = sweeps_relax + steps / Natoms
    with open("sweeps_allF_s"+str(s)+".out", 'a+') as f_handle:
        f_handle.write("%.10f %d\n" %(F, Nsweeps))

    if save_all_F == True:
        Fstr = str( round(F-1, -int(m.floor(m.log10(abs(F-1)))) ) + 1)
        savename = "_F" + Fstr.replace('.','-')
    else: # overwrite 'in_prog' file
        savename = "_in_prog"

    # Weights
    np.savetxt("weights"+savename+"_s"+str(s)+".out", binned_data['w']) 
    
    # Histogram
    np.savetxt("hist"+savename+"_s"+str(s)+".out", binned_data['h'])
  
    return
