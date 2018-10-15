"""
Module.

Imported by pretty much everything.
For almost all applications of the code, this is the only file which the user should directly edit.
"""

import math as m
import numpy as np
import os.path
import sys

# Simulation ID: must be an integer
simID = 0


#################
## Input files ##
#################

# Existing weights file for restart or for final fixed weights computation
params_weights_file = "weights_WL_scomb.out"

# Existing collection matrix
params_Cmat_file = None

# Existing histogram
params_hist_file = None

# Existing series data
params_series_file = None


########################
## Type of simulation ##
########################

# Options: 'wang_landau', 'multicanonical' or 'transition'
algorithm = 'multicanonical'

# Do we want to save a data as a series (for histogram reweighting)?
track_series = True

# Do we want to record the dynamics of the simulation?
track_dynamics = True

# Do we want to use interpolated or discrete weights?
use_interpolated_weights = True


########################
## General parameters ##
########################

# (Reciprocal of) temperature (in units of k_B, i.e. beta)
B = 125.
kT = 1.0/B

# Maximum displacement for a single move (units of SUPERCELL vector)!!!!
dr_max = 0.0177


########################
## Lattice parameters ##
########################

# Number of atoms per lattice
Natoms = 72

# Density of atoms
density = 0.24

# Lattice Alpha
alpha_type = 'fcc'
alpha_vec = np.array([2, 3, 3])
alpha_a = m.pow(4.0/density, 1.0/3.0)

# Lattice Beta
beta_type = 'bcc'
beta_vec = np.array([3, 4, 3])
beta_a = m.pow(2.0/density, 1.0/3.0)
"""
# Uncomment for two fcc lattices
beta_type = alpha_type
beta_vec = alpha_vec
beta_a = alpha_a
"""
# Free energy difference between alpha and beta (get from lattice_energies_diff in initialise.py)
adjust = 0.0397714088

# Ideal lattice energy of alpha (and beta after adjustment)
E_ideal = 16.9104578991


############
## Domain ##
############

# Number of processes
Np = 3

# Number of subdomains
Ns = 3

# Subdomain boundaries: MUST have length of Ns+1
boundaries = (-3.5, -1., 1., 3.5)

# Number of bins for each subdomain
bins = (31,)*Ns

# Type of binning system for each subdomain (only equal width 'eq' currently available)
rules = ('eq',)*Ns

# Trap in one subdomain (True), or allow each random walk to access entire domain (False)
# Must be set to True for Wang Landau simulations
TRAP = False


# Mu overlap between subdomains - the max of these two for each subdomain will be applied
abs_olap = 0.1 # absolute mu overlap
frac_olap = 0.1 # fraction of subdomain width

# Number of points to interpolate within overlap region
Ninterp = 20

# Number of joins to average over for dF calculation
Nsamples_join = 10


#############################
## Step-related parameters ##
#############################

# Number of sweeps between dF evaluations
sweeps_dF = 10000

# How often do we want to save during the simulation (sweeps)
sweeps_save = 10000000

# How often do we want to sample E,mu,eta for our 'time'-series (sweeps)
sweeps_series = 1

# How often do we want to recalculate the whole lattice energy (sweeps)
sweeps_recalc = 1000


############################
## Wang-Landau parameters ##
############################

# Flatness tolerance
flat_tol = 0.8

# Initial F
F_init = 1.05

# Minimum F (lower -> better converged weights)
F_min = 1.000001

# Allow for system to relax to a steady state before starting weights+histogram build-up
sweeps_relax = 1000

# Save for each F value?
save_all_F = True


#################
## Convergence ##
#################

# Number of data points in running average for delta F stdev
window = 10

# Iterations of 'sweeps' for multicanonical/TM
# Set to -1 to run until stdev_converged
iterations = 10

# Standard devation of the last 'window' delta F values required for convergence
stdev_converged = 1e-6


##################################
## Transition matrix parameters ##
##################################

# What method to use to extract eigenvector from P matrix ('sequential', 'arpack')
eigvec_method = 'arpack'

# Tolerance for arpack method
eigs_tol = 1e-13

# How often to refresh TM weights 
# set to greater than 'sweeps_dF' to use input weights the entire time
sweeps_refresh = 100000


##########################
## File names and paths ##
##########################

# Name for file to contain free energy difference series
deltaF_file = "DELTA_F.out"

# Directory of Params.py and Python scripts
pwd = os.path.dirname(os.path.realpath(__file__))

# Set path to interaction potential
path_to_pot = os.path.join(pwd, 'interaction_pots', 'GCM') 
sys.path.insert(1, path_to_pot) # added to path behind pwd

# Set path to Python modules
path_to_mods = os.path.join(pwd, 'modules')
sys.path.insert(2, path_to_mods) # added to path behind pwd and path_to_pot


############################
## Generic error function ##
############################

# Since this module is imported by everything it seems like a decent place to
# put a generic error function
def error(origin, message):
    print "Error! %s" %origin
    print message
    print "Exiting..."
    sys.exit(1)
