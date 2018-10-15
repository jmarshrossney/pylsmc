"""
Script.

Contains checks which should be run before starting a simulation.
Run automatically in run_parallel.sh
"""

import numpy as np
import math as m
import os.path
from sys import argv, exit

from params import *
import energy

import domain as dom
import initialise as ini

# Import module corresponding to the type of algorithm being used
if algorithm == 'wang_landau':
    import wang_landau as alg
elif algorithm == 'multicanonical':
    import multicanonical as alg
elif algorithm == 'transition':
    import transition_matrix as alg

# Name of this file
this_file = os.path.basename(__file__)

############################################################################
##           Functions to check that the input parameters are OK          ##
############################################################################

def check_compatible(atoms_alpha, atoms_beta):
    """ Check for equal atom numbers and densities """
    
    Nalpha = atoms_alpha.get_number_of_atoms()
    Nbeta  = atoms_beta.get_number_of_atoms()

    if Nalpha != Nbeta:
        error(this_file, "Different numbers of atoms in each lattice.")

    if Nalpha != Natoms:
        error(this_file, "Number of atoms does not match 'Natoms' in params.py.")

    dens_alpha = Nalpha / atoms_alpha.get_volume()
    dens_beta  = Nbeta / atoms_beta.get_volume()

    if abs(dens_alpha - dens_beta) > 1e-10:
        error(this_file, "Different density for each lattice.")

    print "Lattice compatibility check OK." 
    return


def calibrate(atoms_alpha, atoms_beta):
    """ Find the difference between the ideal lattice energies, to be
    set as a constant parameter in params.py """
    
    alpha_energy = energy.compute_lattice_energy(
                            atoms_alpha.positions,atoms_alpha.get_cell(),1)
    beta_energy = energy.compute_lattice_energy(
                            atoms_beta.positions,atoms_beta.get_cell(),0)

    energy_ideal = alpha_energy
    energy_diff = alpha_energy - beta_energy

    if np.allclose(energy_ideal, E_ideal) == False: # rtol =1e-5, atol=1e-8
        error(this_file, "Please set E_ideal = %.10f in params.py." %energy_ideal)
    
    if np.allclose(energy_diff, adjust) == False: # rtol =1e-5, atol=1e-8
        error(this_file, "Please set adjust = %.10f in params.py." %energy_diff)

    print "Ideal lattice energies and adjustment OK."
    return


def check_domain():
    """ Check sensible subdomain strategy """
        
    if len(boundaries) != Ns + 1 or len(bins) != Ns or len(rules) != Ns:
        error(this_file, "Mismatch between number of subdomains and number of boundaries/bins/rules. \nMake sure there are Ns + 1 boundaries, Ns bins and Ns rules in params.py!")

    for i in range(len(boundaries)-1):
        if boundaries[i+1] <= boundaries[i]:
            error(this_file, "Mu boundaries do not increase monotonically from left to right.")

    # These are center-of-bin values
    mu_s_cross = dom.get_local_mu_bins(dom.s_cross)
    
    # Want a bin centered on mu=0
    mu_hi = mu_s_cross[dom.i_cross]
    mu_lo = mu_s_cross[dom.i_cross - 1]
    
    # If we've failed to get a bin centered on mu=0
    if np.allclose(mu_hi,0.0) == False: # rtol=1e-5, atol=1e-8
        
        print "Error:", this_file
        print "mu = 0 falls into a bin centered on mu = %1.8f, rather than 0." %mu_hi

        bound_hi = boundaries[dom.s_cross + 1]
        bound_lo = boundaries[dom.s_cross]

        rnd_bin_width = np.around(dom.subdom[dom.s_cross]['bin_width'],decimals=10)
        hbw = 0.5*rnd_bin_width
        bins_lo = dom.i_cross - dom.subdom[dom.s_cross]['l_bin_olap']
        bins_hi = bins[dom.s_cross] - bins_lo
        prop_bound_hi = bins_hi * rnd_bin_width
        prop_bound_lo = bins_lo * rnd_bin_width * -1
        
        print "Try %.10f --> %.10f,    %.10f --> %.10f" \
                %(bound_lo,prop_bound_lo+hbw,bound_hi,prop_bound_hi+hbw)
        print "Or %.10f --> %.10f,    %.10f --> %.10f" \
                %(bound_lo,prop_bound_lo-hbw,bound_hi,prop_bound_hi-hbw)
        
        print "Or change the number of bins and try again."
        exit(1)

    print "Subdomain strategy OK."
    return

def final_check():
    """ Just some final checks """

    if algorithm not in ('wang_landau','multicanonical','transition'):
        error(this_file, "algorithm should be either 'wang_landau', 'multicanonical', or 'transition'.")
    
    if type(track_series) != bool:
        error(this_file, "Track_series should be set to either True or False.")

    if type(use_interpolated_weights) != bool:
        error(this_file, "use_interpolated_weights should be set to either True or False.")

    if (abs_olap>=0) == False or (0<=frac_olap<1) == False:
        error(this_file, "Please set abs_olap >= 0 and 0 <= frac_olap < 1.")

    if type(TRAP) != bool:
        error(this_file, "TRAP should be set to either True or False.")

    if eigvec_method not in ('sequential', 'arpack'):
        error(this_file, "eigvec_method should be either 'sequential' or 'arpack'.")
    
    # TRAP should be True for a weights-buildup
    if algorithm == 'wang_landau':
        if Ns != 1:
            if TRAP == False:
                error(this_file, "Please set TRAP = True for a Wang Landau simulation.")

    print "Final params check OK."
    return

def gen_input_files():
    """ If starting a new fixed-weights simulation, empty files will need
    to be created to be imported during the first loop iteration. """
   
    if algorithm not in ('multicanonical', 'transition'):
        return

    print ""

    for p in range(Np):
        # Map to subdomain index
        s = ini.map_proc_to_subdom(p)

        # Input files as decided by program
        input_files = alg.file_names('input', s, p)

        sname, pname = ini.naming_system(s, p)

        # User input files from params
        params_input_files = {'d': None,
                           'w': params_weights_file,
                           'h': params_hist_file,
                           'c': params_Cmat_file,
                           's': params_series_file,}
        # Check if we need to add '_pY_sX'
        params_weights_file_alt = ini.rename_inputs(params_weights_file, sname[0], '')
        params_hist_file_alt = ini.rename_inputs(params_hist_file, sname[0], pname[0])
        params_Cmat_file_alt = ini.rename_inputs(params_Cmat_file, sname[0], pname[0])
        params_series_file_alt = ini.rename_inputs(params_series_file, sname[0], pname[0])
        params_input_files_alt = {'d': None,
                               'w': params_weights_file_alt, 
                               'h': params_hist_file_alt,
                               'c': params_Cmat_file_alt,
                               's': params_series_file_alt}

        # Correct sizes
        size = ini.get_size(s)
        start_files = {'d': np.zeros((Natoms, 3)), 
                'w': np.zeros(size), 'h': np.zeros(size), 
                'c': np.zeros((size,size)), 's': [np.ones(5)*-1,]}
      
        for key in params_input_files.keys():
            # If an input file has been given in params.py
            if params_input_files[key] != None:
                # If the default input file doesn't yet exist,
                # We need to save a copy of the params file to this name
                if os.path.exists(input_files[key]) == False:
                    # First check for the params file modified by '_pY_sX'
                    if os.path.exists(params_input_files_alt[key]) == True:
                        tmp = np.loadtxt(params_input_files_alt[key])
                        print "     Copying %s -> %s" %(params_input_files_alt[key], input_files[key])
                        np.savetxt(input_files[key], tmp)
                    # Then check for the exact file name given in params.py 
                    elif os.path.exists(params_input_files[key]) == True:
                        tmp = np.loadtxt(params_input_files[key])
                        print "     Copying %s -> %s" %(params_input_files[key], input_files[key])
                        np.savetxt(input_files[key], tmp)
                    # Error if params file hasn't been found
                    else:
                        error(this_file, "     Neither files %s or %s not found." %(params_input_files_alt[key],params_input_files[key]))

                # If the default input file does exist, leave it alone
                else: 
                    print "     File %s exists and will be read as input" %input_files[key]
        
    
    print ""

    print "File check completed. "
    return


def check_shuffle_file():
    """ Check for a file for shuffled indices. Create a new one if not. """

    # Check if shuffle file exists
    input_file = "shuffled_indices.out"
    if os.path.exists(input_file) != True:
        
        # Any existing disp files should be deleted if the shuffle file has been lost
        disp_file = alg.file_names('input', 0, 0)['d']
        if os.path.exists(disp_file) == True:
            error(this_file, "Unfortunately, since the shuffle file has been lost, the initial displacement vectors must be reset to zero. Please move or delete files %s etc." %disp_file)
    
        # Create a new array and file for shuffled indices
        shuffled_indices = np.arange(Natoms)
        np.random.shuffle(shuffled_indices)
        print "Saving new shuffled indices to ", input_file
        np.savetxt(input_file, shuffled_indices.astype(int), fmt='%i')


#################################################################
## Call checking script to check specific groups of parameters ##
#################################################################

if len(argv) > 1:
    if argv[1] == 'lattice':
        # Create temporary supercells to perform checks
        atoms_alpha = ini.build_supercell(np.zeros((Natoms,3)), alpha_vec, alpha_a, alpha_type)
        atoms_beta = ini.build_supercell(np.zeros((Natoms,3)), beta_vec, beta_a, beta_type)
        # Check that the supercells are compatible
        check_compatible(atoms_alpha, atoms_beta)
        
        # Compare measured ideal lattice energies with those in params.py
        calibrate(atoms_alpha, atoms_beta)

    elif argv[1] == 'domain':
        # Check that there is a bin boundary at mu = 0
        check_domain()

    elif argv[1] == 'final':
        # Check sensible params.py values
        final_check()
    
    exit(0)

############################################################################
##                              Execute Checks                            ##
############################################################################


# Create temporary supercells to perform checks
atoms_alpha = ini.build_supercell(np.zeros((Natoms,3)), alpha_vec, alpha_a, alpha_type)
atoms_beta = ini.build_supercell(np.zeros((Natoms,3)), beta_vec, beta_a, beta_type)

# Check that the supercells are compatible
check_compatible(atoms_alpha, atoms_beta)

# Compare measured ideal lattice energies with those in params.py
calibrate(atoms_alpha, atoms_beta)

# Check that there is a bin boundary at mu = 0
check_domain()

# Check sensible params.py values
final_check()

# Generate input files if starting new fixed-weights run
gen_input_files()

# Check for a file containing the order by which atoms have been shuffled
if alpha_type == beta_type:
    check_shuffle_file()


print "Initial checks completed."
print ""
print "============================================================================="
print ""

###########################################################################
##                        Print simulation info                          ##
###########################################################################

print "Path to interaction potential: ", path_to_pot
print ""
print "Lattice Alpha: ", alpha_type, " with supercell vectors:"
print atoms_alpha.get_cell()
print "Lattice Beta: ", beta_type, " with supercell vectors:"
print atoms_beta.get_cell()
print ""
print "E_ideal (Alpha ideal energy): ", E_ideal
print "Adjust (Ideal energy difference Alpha-Beta): ", adjust
print "Number of atoms: ", Natoms
print "Density: ", density
print "Thermodynamic Beta: ", B
print "Max step size: ", dr_max
print "Using interpolated weights is set to: ", use_interpolated_weights
print ""
print "Global domain: %f to %f with %d bins" \
        %(boundaries[0],boundaries[-1],np.sum(bins))
print "Number of subdomains: ", Ns
print "Trapping simulations in one subdomain is set to: ", TRAP
print "Subdomain      -olap  |   boundaries   |  +olap        bins    bin width"
for s in range(Ns):
    print "%02d             %06.2f | %06.2f  %06.2f | %06.2f         %d       %05.3f" \
            %(s,dom.subdom[s]['min'], boundaries[s],boundaries[s+1],dom.subdom[s]['max'],
                    dom.subdom[s]['bins'],dom.subdom[s]['bin_width'])
print ""
print "Simulation type: ",algorithm
if algorithm == 'wang_landau':
    print "Initial F: ", F_init
    print "Final F: ", F_min
    print "Flatness tolerance: ", flat_tol
if algorithm == 'transition':
    print "P matrix eigenvector computed using: ", eigvec_method
    print "Weights refreshed every ", min(sweeps_refresh,sweeps_dF), " sweeps"
if algorithm in ('multicanonical', 'transition'):
    print "dF computed after every iteration of ", sweeps_dF, " sweeps"
    if iterations < 0:
        print "Running until the standard deviation on dF is less than ", stdev_converged
        print "Standard deviation calculated with a running average of ", window, "dF's"
    else:
        print "Running for ", iterations, " iterations"
if track_series == True:
    print ""
    print "Tracking data series: saving every %d steps" %int(sweeps_series/Natoms)
print "Measuring simulation dynamics is set to: ", track_dynamics
print ""
print "Splitting into ", Np, "processes"
print ""
print "============================================================================="
print ""
