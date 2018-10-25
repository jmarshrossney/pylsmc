import numpy as np
import math as m
from os.path import basename
import sys

# Import params from parent directory
path_to_params = '..'
sys.path.insert(0, path_to_params) # looks in parent directory first
from params import *
sys.path.remove(path_to_params)

# Import lsmc modules from this directory
import initialise as ini
import domain as dom

# Name of this file
this_file = basename(__file__)


#############################
## User-defined parameters ##
#############################

# Update cuts histogram after how many steps?
cuts_step = 8

# Update 'mini' matrix with transition probabilities after how many steps?
minimat_step = 8

# INT Number of mini-bins per bin. Odd -> has a mini-bin centred on the bin centre
minibins_per_bin = 10

# Number of regular bins per window. Remember the window is centred on a bin.
bins_per_window = 4

# Set integer number of mini-bins per window based on above
minibins_per_window = int(round(minibins_per_bin * bins_per_window))



#################################
## File naming and i/o control ##
#################################
def file_names(stage, s=0, p=0):
    """ Returns a dictionary of file names corresponding to each stage of the program """

    # Get the correct names for p, s
    sname, pname = ini.naming_system(s, p)


    if stage in ('input', 'output'):
        leaf = pname[0] + sname[0] + ".out"
    elif stage == 'scomb':
        leaf = pname[0] + sname[1] + ".out"
    else:
        error(this_file, "Invalid argument: "+stage)

    cuts = "cuts" + leaf
    accepts = "accepts" + leaf
    dmu = "dmu" + leaf
    dmu_acc = "dmu_acc" + leaf
    dE = "dE" + leaf
    dE_acc = "dE_acc" + leaf
    minimat = "minimat" + leaf
    round_trips = "round_trips" + leaf
        
    return_dict =  {'cu': cuts, 'ac': accepts, \
                    'dm': dmu, 'dma': dmu_acc, \
                    'de': dE, 'dea': dE_acc, \
                    'rt': round_trips}
    
    if TRAP == True:
        return_dict['mm'] = minimat

    return return_dict


def load_inputs(s, p):
    """ Load the correct input files for the extra info """

    input_files = file_names('input', s, p)

    size = ini.get_size(s)
    
    # Load/initialise binned data in the usual way
    cuts = ini.file_input(input_files['cu'], size)
    accepts = ini.file_input(input_files['ac'], size)
    dmu = ini.file_input(input_files['dm'], size)
    dmu_acc = ini.file_input(input_files['dma'], size)
    dE = ini.file_input(input_files['de'], size)
    dE_acc = ini.file_input(input_files['dea'], size)
    
    round_trips = ini.file_input(input_files['rt'], 5)
    if int(round_trips[-1] == 0): round_trips[-1] = 1 # initialise round trip flag

    return_dict =  {'cu': cuts, 'ac': accepts, \
                    'dm': dmu, 'dma': dmu_acc, \
                    'de': dE, 'dea': dE_acc, \
                    'rt': round_trips}
    
    # Fine-grained matrix only works for TRAP == True (due to indexing issues)
    if TRAP == True:
        mm_size = (minibins_per_window, size)
        minimat = ini.file_input(input_files['mm'], mm_size)
        return_dict['mm'] = minimat

    return return_dict


##############################
## Set up fine-grained grid ##
##############################
""" Sets up 'windows' centered on each regular bin, containing 'mini-bins'.
    This allows more precise measurement of local transition probabilities,
    while keeping relatively thick bins for everything else. """

# ------------------------- #
#  Iterate over subdomains  #
# ------------------------- #
subdom = []
for s in range(Ns):

    # These are the centres of the windows
    mu_bins = dom.get_local_mu_bins(s)
    
    # Widths in terms of mu
    minibin_width = dom.subdom[s]['bin_width'] / float(minibins_per_bin)
    window_width = minibin_width * minibins_per_window

    # Dictionary of info for this subdomain, including a list of dictionaries
    # containing info about individual windows
    windows = {'mbin_width': minibin_width,
               'win_width': window_width,
               'win_info': []}

    # ---------------------------------- #
    #  Iterate over regular bin indices  #
    # ---------------------------------- #
    for iloc in range(dom.subdom[s]['bins']):

        win_info = {'min': mu_bins[iloc] - 0.5*window_width,
                    'max': mu_bins[iloc] + 0.5*window_width}

        # Mu coordinates of the centres of the mini-bins
        win_info['bin_mus'] = np.ones(minibins_per_window) * win_info['min'] + \
                              np.arange(minibins_per_window) * \
                              (win_info['max']-win_info['min']) / float(minibins_per_window) \
                              + 0.5 * minibin_width
        
        # Add info about this window to the dictionary for this subdomain
        windows['win_info'].append(win_info)


    subdom.append(windows)


def get_minibin_index(mu, s, root_index):
    """ Returns index for minibin within a given window (bin index) """
    return dom.eq_width(mu,
                subdom[s]['win_info'][root_index]['min'],
                subdom[s]['win_info'][root_index]['max'],
                minibins_per_window)


###########################################
## Functions to update data at each step ##
###########################################


def update_cuts(cuts, mu_bins, new_mu, old_mu):
    """ Updates cuts histogram by 1 for each 'cut' (bin edge)
        that's been crossed """
    mu_max = max(new_mu, old_mu)
    mu_min = min(new_mu, old_mu)
    cuts[(mu_bins<mu_max) & (mu_bins>mu_min)] += 1
    return


def update_func(dyn_data, mu_bins, mu, bin_index, s, \
        cuts, rtrips, accepted, dmu, dE):
    """ Update all of the extra binned data """

    # Update cuts if necessary
    cuts['counter'] += 1
    if cuts['counter'] == cuts_step:
        update_cuts(dyn_data['cu'], mu_bins, cuts['new_mu'], cuts['old_mu'])
        cuts['old_mu'] = cuts['new_mu']
        cuts['counter'] = 0

    # Check for round trip completion
    if bin_index == rtrips['minmax'][rtrips['flag']]: # assume bin width > step size
        rtrips['counts'] += 1 # no. half round trips
        rtrips['flag'] *= -1 # swap rt_flag to point to other min/max

    # Update whether move was accepted or not
    dyn_data['dm'][bin_index] += dmu # mu step size
    dyn_data['de'][bin_index] += dE # E step size
    
    # Only update if move was accepted
    if accepted == True:
        dyn_data['ac'][bin_index] += 1 # accepts
        dyn_data['dma'][bin_index] += dmu # mu step size
        dyn_data['dea'][bin_index] += dE # E step size
    
    return


def update_minimat(dyn_data, mu_bins, mu, bin_index, s, \
        minimat, accepted, expnt):
    """ Update mini transition matrix if TRAP == True """   

    # Product of probabilities over several steps
    if accepted == False:
        P_this_step = 1 - np.exp( min(0, expnt) )
    else: # covers edge cases: expnt = 0 -> P multipled by 1
        P_this_step = np.exp( min(0, expnt) )

    # Update unbiased probabilities for mini transition matrix
    minimat['Pprod'] = minimat['Pprod'] * P_this_step
    
    # Get minibin index and enforce boundaries
    new_minibin_index = get_minibin_index(mu, s, minimat['root_index'])    
    if new_minibin_index < 0:
        new_minibin_index = 0
    elif new_minibin_index > minibins_per_window-1:
        new_minibin_index = minibins_per_window-1

    # Update mini transition matrix if necessary
    minimat['counter'] += 1
    if minimat['counter'] == minimat_step:
        # Update with the probability of being in this minibin after
        # 'minimat_step' steps, starting from 'root_index'
        dyn_data['mm'][new_minibin_index,    minimat['root_index']] += minimat['Pprod']
        dyn_data['mm'][minimat['old_index'], minimat['root_index']] += (1-minimat['Pprod'])
    
        # Reset counter, product of probabilities and root index
        minimat['counter'] = 0
        minimat['Pprod'] = 1
        minimat['root_index'] = bin_index

    # Update old minibin index
    minimat['old_index'] = new_minibin_index
    
    return 

