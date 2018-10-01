"""
Module.

Reads domain parameters from params.py and sets up a domain.
"""

import numpy as np
import math as m

from params import *


##################################
## Container for subdomain info ##
##################################

subdom = []

for s in range(Ns):

    # Populate subdomain info based on params.py
    sd_info = { 'min': boundaries[s],
                'max': boundaries[s+1],
                'bins': bins[s],
                'rules': rules[s]
            }
    
    # Bin width seems useful while we're only using eq widths
    mu_width = boundaries[s+1] - boundaries[s]
    bin_width = mu_width / bins[s] # only works for eq width!
    sd_info['bin_width'] = bin_width
        
    # Take overlap into account
    if TRAP == True:
        
        # Calculate using abs olap
        bin_olap_useabs = int(m.floor(abs_olap / bin_width)) + 1

        # Calculate using frac olap
        bin_olap_usefrac = int(m.floor((frac_olap*mu_width) / bin_width)) + 1
       
        # Use larger one
        bin_olap = max(bin_olap_useabs, bin_olap_usefrac)
        mu_olap = bin_olap * bin_width

        # Lower boundary
        if s == 0:
            sd_info['l_bin_olap'] = 0
        else:
            sd_info['min'] -= mu_olap
            sd_info['l_bin_olap'] = bin_olap
            sd_info['bins'] += bin_olap

        # Upper boundary
        if s == Ns - 1:
            sd_info['h_bin_olap'] = 0
        else:
            sd_info['max'] += mu_olap
            sd_info['h_bin_olap'] = bin_olap
            sd_info['bins'] += bin_olap

    subdom.append(sd_info)


# Compatible with 1 subdomain
if len(subdom) == 1:
    subdom.append(None)


########################
## Rules for indexing ##
########################

def eq_width(mu, mu_min, mu_max, N_bins):
    index = m.floor((mu-mu_min)/(mu_max-mu_min) * N_bins)
    return int(index)

# ------------------- Container for these functions ------------------- 
# This could eventually allow different binning systems to be used, but
# currently isn't needed so I've stopped working on it. As a result,
# later developments have assumed equally spaced bins, which would need
# to be changed before new binning systems could be introduced
indexfunc = {'eq': eq_width}
# --------------------------------------------------------------------


#####################################
## Functions for binning mu values ##
#####################################

def get_subdomain(mu):
    """ Returns subdomain index for a given mu, ignoring overlaps """
    for s in range(Ns):
        if boundaries[s] <= mu and boundaries[s+1] > mu:
            this_s = s
            break
    return this_s

def get_local_index(mu, s):
    """ Returns index of bin within subdomain, including overlaps """
    return indexfunc[rules[s]](mu,subdom[s]['min'],subdom[s]['max'],subdom[s]['bins'])

def loc2glob_index(iloc, s):
    iglob = iloc + np.sum(bins[:s])
    return int(iglob)

def get_global_index(mu, s):
    """ Returns global index of bin, ignoring any overlap bins"""
    iloc = indexfunc[rules[s]](mu,boundaries[s],boundaries[s+1],bins[s])
    iglob = loc2glob_index(iloc, s)
    return iglob


######################################################
## Arrays of bin mu positions (useful for plotting) ##
######################################################

def get_global_mu_bins():
    """ Return array for mu position of bins over entire domain """

    ind = np.arange(0, np.sum(bins), 1)
    mu_list = []

    for s in range(Ns):
        mu_min = boundaries[s]
        mu_max = boundaries[s+1]
        Nbins = bins[s]

        # Assumed equal widths
        mu_s = np.ones(Nbins)*mu_min + np.arange(0,Nbins,1)*(mu_max-mu_min)/float(Nbins)
        mu_s += 0.5*subdom[s]['bin_width'] # centered on bin
        mu_s = list(mu_s)
        mu_list = mu_list + mu_s

    mu_all_s = np.array(mu_list)

    return mu_all_s

def get_local_mu_bins(s):
    """ Return array for mu position of bins in one subdomain,
    INCLUDING OVERLAPS """

    mu_min = subdom[s]['min']
    mu_max = subdom[s]['max']
    Nbins = subdom[s]['bins']

    # Assumed equal widths
    mu_s = np.ones(Nbins)*mu_min + np.arange(0,Nbins,1)*(mu_max-mu_min)/float(Nbins)
    mu_s += 0.5*subdom[s]['bin_width'] # centered on bin

    return mu_s


############################################
## Crossover from negative to positive mu ##
############################################

def find_crossover():
    """ Find the subdomain and global index of mu=0 bin """
   
    s_cross = get_subdomain(1e-9) # 0 does nasty things - use a v small positive number
    i_cross = get_local_index(1e-9, s_cross)
    I_cross = get_global_index(1e-9, s_cross)

    return s_cross, i_cross, I_cross

s_cross, i_cross, I_cross = find_crossover()

