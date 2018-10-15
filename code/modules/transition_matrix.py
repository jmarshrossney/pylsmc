"""
Module.

Contains functions specific to the transition matrix, e.g. normalisation of columns, diagonalisation, TTT identity (detailed balance) check.
"""

import math as m
import numpy as np
from scipy.sparse.linalg import eigs
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


def extrapolate(array):
    """ Adds linearly extrapolated end bins to a numpy array, assuming the
        array has spaces at either end for this extrapolated value to go """
    array[-1] = 2*array[0] - array[1] # this is really to the left of bin 0
    array[-2] = 2*array[-3] - array[-4]
    return


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
        weights_in = "weights_TM" + pname[1] + sname[0] + ".out" # pcomb: combined within each subdom
        hist_in = "hist_TM" + pname[0] + sname[0] + ".out"
        Cmat_in = "Cmat_TM" + pname[0] + sname[0] + ".out"
        return {'d': disp_in, 'w': weights_in, 'h': hist_in, 'c': Cmat_in}
    
    elif stage == 'output':
        # Output file names include appropriate '_pY_sX'
        disp_out = "disp" + pname[0] + sname[0] + ".out"
        weights_out = "weights_TM" + pname[0] + sname[0] + ".out"
        hist_out = "hist_TM" + pname[0] + sname[0] + ".out"
        Cmat_out = "Cmat_TM" + pname[0] + sname[0] + ".out"
        eigvec_out = "eigvec_TM" + pname[0] + sname[0] + ".out"
        series_out = "series_TM" + pname[0] + sname[0] + ".out"
        return {'d': disp_out, 'w': weights_out, 'h': hist_out, \
                'c': Cmat_out, 'e': eigvec_out, 's': series_out}

    elif stage == 'pcomb':
        # After combined over processors
        weights_pcomb = "weights_TM" + pname[1] + sname[0] + ".out"
        hist_pcomb = "hist_TM" + pname[1] + sname[0] + ".out"
        Cmat_pcomb = "Cmat_TM" + pname[1] + sname[0] + ".out"
        eigvec_pcomb = "eigvec_TM" + pname[1] + sname[0] + ".out"
        series_pcomb = "series_TM" + pname[1] + sname[0] + ".out"
        return {'w': weights_pcomb, 'h': hist_pcomb, \
                'c': Cmat_pcomb, 'e': eigvec_pcomb, 's': series_pcomb}

    elif stage == 'scomb':
        # After subdomain join
        weights_scomb = "weights_TM" + pname[1] + sname[1] + ".out"
        hist_scomb = "hist_TM" + pname[1] + sname[1] + ".out"
        Cmat_scomb = "Cmat_TM" + pname[1] + sname[1] + ".out"
        eigvec_scomb = "eigvec_TM" + pname[1] + sname[1] + ".out"
        return {'w': weights_scomb, 'h': hist_scomb, 'c': Cmat_scomb, 'e': eigvec_scomb}

    else:
        error(this_file, "Invalid argument: "+stage)


def load_inputs(s, p):
    """ Load the correct input files for a Transition Matrix simulation """

    input_files = file_names('input', s, p)

    disp = ini.file_input(input_files['d'], (Natoms, 3))

    size = ini.get_size(s)
    weights = ini.file_input(input_files['w'], size)
    hist = ini.file_input(input_files['h'], size)
    Cmat = ini.file_input(input_files['c'], (size,size) )

    # Initialise eigvec as a random vector
    eigvec = np.random.random(size)

    return disp, {'w': weights, 'h': hist, 'c': Cmat, 'e': eigvec}


############################################################
### Functions called periodically from lattice_switch.py ###
############################################################
def refresh_func(binned_data):

    # Find eigenvector and associated weights from transition matrix
    binned_data['e'], TMweights, v = get_Pmat_eigenvector(binned_data['c'], binned_data['e'])
    
    if use_interpolated_weights == True:
        # Need to add bins to the edges of mu and weights for extrapolation
        TMweights = np.concatenate( (TMweights, [0,0]) )
        extrapolate(TMweights)
        TMweights -= np.min(TMweights)

    print "Refreshing weights with log of eigenvector"
    for i in range(len(binned_data['w'])):
        binned_data['w'][i] = TMweights[i] # easy to switch off for testing
            
    return


def update_func(step, binned_data, delta_local_energies, 
            new_index, old_index, old_index_copy, 
            old_mu, mu_bins,
            old_weight, get_weight, lendata, kTlogF, s): 
    """ Update appropriate binned data after a single-particle move has been accepted or rejected 
        Old_index has been updated according to move success already - must use old_index_copy """
   
    # Update Cmat with unbiased probabilities (indices = before/after move)
    P_unbiased = np.exp( min(0, -B*delta_local_energies) )
    binned_data['c'][new_index, old_index_copy] += P_unbiased# * sampling_adjust
    binned_data['c'][old_index_copy, old_index_copy] += (1-P_unbiased)# * sampling_adjust

    # Update histogram (index = result of move)
    binned_data['h'][old_index] += 1

    # Second returned value is FLAT, for compatibility with Wang Landau simulation
    return old_weight, 1


def save_func(binned_data, mu_bins, step, s, p):
    """ Save appropriate binned data periodically during a simulation """

    print "Sweep %d.  Saving and continuing..." %(step/Natoms)

    # If interpolated weights used, remove those extrapolated bins / 'ghost points'
    if use_interpolated_weights == True:
        binned_data['w'] = binned_data['w'][:-2]

    savename = "_in_prog_p"+str(p)+"_s"+str(s)
    np.savetxt("weights"+savename+".out", binned_data['w'])
    np.savetxt("hist"+savename+".out", binned_data['h'])
    np.savetxt("Cmat"+savename+".out", binned_data['c'])
    np.savetxt("eigvec"+savename+".out", binned_data['e'])


################################
### Matrix-related functions ###
################################
def sequential(Pmat, eigvec_init):
    """ Uses sequential evaluation (first off-diagonals only) to compute eigenvector
        of P matrix """
    eigvec = np.zeros(len(eigvec_init))
    eigvec[0] = 1.0
    for i in range(len(eigvec)-1):
        eigvec[i+1] = np.exp( np.log(eigvec[i]) + np.log(Pmat[i+1,i] / Pmat[i,i+1]) )
    eigvec = eigvec - np.min(eigvec) + 0.000001 # assume pdf goes to ~0 at edges
    eigvec = eigvec / np.sum(eigvec) # normalise
    return eigvec

def arpack(Pmat, eigvec_init):
    """ Uses python wrapper to Arpack to compute eigenvector of P matrix """
    eigval, eigvec = eigs(Pmat, k=1, which='LM', v0 = eigvec_init, tol=eigs_tol)
    eigvec_abs = np.abs(eigvec)
    return eigvec_abs / np.sum(eigvec_abs)

def Cmat_to_Pmat(Cmat):
    """ Normalise each column to make a stochastic matrix """
    col = np.sum(Cmat, axis=0)
    col[ col==0 ] = 1
    Pmat = Cmat/col

    return Pmat

def get_Pmat_eigenvector(Cmat, eigvec):
    """ Returns the transition matrix eigenvector and associated
        multicanonical weights estimate """

    lendata = len(eigvec)
    undersampled = np.where(Cmat.diagonal()<1000)[0]
        
    # Make a fresh copy - Cmat shouldn't be changed here
    Cmat_copy = np.copy(Cmat)

    # Pretend there aren't zeros on diagonal if there are
    if len(undersampled) != 0:
        Cmat_copy += np.diag(np.ones(lendata))
        
    # Symmetrise zeros in second off-diagonal
    for j in range(lendata-3):
        if Cmat_copy[j,j+2] == 0: Cmat_copy[j+2,j] = 0
        if Cmat_copy[j+2,j] == 0: Cmat_copy[j,j+2] = 0
             
    # Separately normalise transitions from each bin
    Pmat = Cmat_to_Pmat(Cmat_copy)
                
    # Check TTT identity
    v = np.zeros(lendata)
    for i in range(lendata-2):
        left = Pmat[i+1,i] * Pmat[i+2,i+1] * Pmat[i,i+2]
        right = Pmat[i+2,i] * Pmat[i+1,i+2] * Pmat[i,i+1]
        if left != 0 and right != 0:
            v[i+1] = 1 - left / right
        else:
            v[i+1] = float('NaN')
    v_abs = np.abs(v)
    j = np.nanargmax(v_abs)
    print "TTT identity: detailed balance violation: max is %f in col %d" %(v_abs[j],j)
        
    # Pull out an estimate for the probability
    print "Attempting to compute eigenvector using ", eigvec_method
    if eigvec_method == 'arpack':
        eigvec = arpack(Pmat, eigvec)
    elif eigvec_method == 'sequential':
        eigvec = sequential(Pmat, eigvec)
        
    # Scale by bin width and renormalise
    if TRAP == False:
        mu = dom.get_global_mu_bins()
        for i in range(len(eigvec)):
            s = dom.get_subdomain(mu[i])
            w = dom.subdom[s]['bin_width']
            eigvec[i] = eigvec[i] / w
        eigvec = eigvec / np.sum(eigvec)

    # Compute weights
    TMweights = np.log(eigvec)
    TMweights -= np.min(TMweights)
    TMweights = kT * TMweights
    
    # Smooth out undersampled bits
    # I don't know what to put here right now

    return eigvec, TMweights, v




