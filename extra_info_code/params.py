"""
Module.

Imported by pretty much everything.
For almost all applications of the code, this is the only file which the user should directly edit.
"""

import math as m
import numpy as np
from sys import exit

# Simulation ID: must be an integer
simID = 0


#################
## Input files ##
#################

# Existing weights file for restart or for final fixed weights computation
params_weights_file = None

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
track_series = False

# Do we want to use interpolated or discrete weights?
use_interpolated_weights = False


########################
## General parameters ##
########################

# Number of atoms per lattice
Natoms = 72

# (Reciprocal of) temperature (in units of k_B, i.e. beta)
B = 333.33
kT = 1.0/B

# Density of atoms
density = 0.24

# Maximum displacement for a single move (units of SUPERCELL vector)
dr_max = 0.0105


########################
## Lattice parameters ##
########################

# Lattice Alpha
alpha_type = 'fcc'
alpha_vec = np.array([2, 3, 3])
alpha_a = m.pow(4.0/density, 1.0/3.0)
"""
# Lattice Beta
beta_type = 'bcc'
beta_vec = np.array([3, 4, 3]) 
beta_a = m.pow(2.0/density, 1.0/3.0)
"""
# Uncomment for two fcc lattices
beta_type = alpha_type
beta_vec = alpha_vec
beta_a = alpha_a

# Free energy difference between alpha and beta (get from lattice_energies_diff in initialise.py)
adjust = 0.0

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
boundaries = (-1., -0.25, 0.25, 1.)

# Number of bins for each subdomain
bins = (31,21,31)

# Type of binning system for each subdomain
rules = ('eq',)*Ns

# Trap in one subdomain (True), or allow each random walk to access entire domain (False)
# Must be set to True for Wang Landau simulations
TRAP = True


########################
## Joining subdomains ##
########################

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
sweeps_save = 1000000

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

# Minimum F (lower -> flatter histogram), must be less than F_init
F_min = 1.00001

# Allow for system to relax to a steady state before starting weights+histogram build-up
sweeps_relax = 1000

# Save for each F value?
save_all_F = False


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


################
## File names ##
################

deltaF_file = "DELTA_F.out"


############################
## Generic error function ##
############################

# Since this module is imported by everything it seems like a decent place to
# put a generic error function
def error(origin, message):
    print "Error! %s" %origin
    print message
    print "Exiting..."
    exit(1)
