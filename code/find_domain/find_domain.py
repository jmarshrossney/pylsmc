"""
Script.

Runs a Boltzmann-weighted Metropolis-Hastings (no switches) simulation on one side of mu=0,
in order to determine appropriate boundaries for future LSMC simulations using the same parameters.

All you really need to do to use this is copy in a params.py file with your chosen parameters.

Make sure sweeps_dF in params.py is set to a reasonably large number so that the simulation has
time to explore regions of configuration space with large lattice energy differences (|mu|).

This is also useful for tuning the max step size so that 50% of moves are accepted.
Really this tuning should happen automatically in this script, but I haven't gotten round to it.

Run as      parallel python find_domain.py ::: 0 1
"""

import numpy as np
import sys
from os.path import basename
from ase.lattice.cubic import FaceCenteredCubic, BodyCenteredCubic

# Import params from parent directory
# DO NOT import from this directory as the path to the interaction potential will be wrong
path_to_params = '..'
sys.path.insert(0, path_to_params) # looks in parent directory first
from params import *
sys.path.remove(path_to_params)

# Import energy using path given in params.py
import energy

# Name of this file
this_file = basename(__file__)

# Subdomain index (should be 0 or 1 for each side of mu=0)
s = int(sys.argv[1])

# Choose the size of the bins for this simulation
bin_width = 0.05*kT*Natoms 


#####################
## Local functions ##
#####################
def build_supercell(disp, vec, a, Ltype, shuffle=False):
    """ Build a supercell using a given set of lattice vectors, lattice constant and type"""

    if Ltype == 'fcc':
        build = FaceCenteredCubic
    elif Ltype == 'bcc':
        build = BodyCenteredCubic
    else:
        error(this_file, "Only fcc and bcc supported right now. Please set alpha_type and beta_type to 'fcc' or 'bcc' in params.py")

    atoms = build(size=tuple(vec),
                    symbol='Ar', # inconsequential
                    pbc=(1,1,1),
                    latticeconstant=a)
    
    # Translate so that center of mass is at (0,0,0)
    com = atoms.get_center_of_mass()
    atoms.translate(-com)

    # Haven't yet worked out how to translate the cell..
    cell = np.array([[vec[0]*a,0.0,0.0],
                     [0.0,vec[1]*a,0.0],
                     [0.0,0.0,vec[2]*a]])
    atoms.set_cell(cell)

    # Shuffle atomic indices of one lattice if both lattices are same type
    if shuffle != False:
        pos = atoms.positions
        shuffled_indices = np.arange(len(pos))
        np.random.shuffle(shuffled_indices)
        atoms.set_positions(pos[shuffled_indices])

    # Set the atomic positions according to disp
    ideal_positions = atoms.get_positions()
    new_positions = ideal_positions + np.dot(disp, atoms.get_cell())
    atoms.set_positions(new_positions)

    return atoms
        
def eq_width(mu, mu_min, mu_max, N_bins):
    index = m.floor((mu-mu_min)/(mu_max-mu_min) * N_bins)
    return int(index)

def local_energy(atoms_alpha, atoms_beta, ipart):
    """ Compute energies of ipart interacting with its neighbours """
    alpha_local_energy = energy.compute_local_energy(ipart,
                                atoms_alpha.positions,atoms_alpha.get_cell(),1)
    beta_local_energy = energy.compute_local_energy(ipart,
                                atoms_beta.positions,atoms_beta.get_cell(),0)
    return np.array([0, alpha_local_energy, beta_local_energy]) # alpha: 1, beta: -1

def lattice_energies(atoms_alpha, atoms_beta):
    """ Compute energies of lattices, return difference """
    alpha_energy = energy.compute_lattice_energy(
                            atoms_alpha.positions,atoms_alpha.get_cell(),1)
    beta_energy = energy.compute_lattice_energy(
                            atoms_beta.positions,atoms_beta.get_cell(),0)
    beta_energy += adjust # make the ideal lattice energies the same
    return alpha_energy, beta_energy

def accept_move(expnt):
    """ Metropolis acceptance criteria """
    expnt = min(0, expnt)
    P = np.exp(expnt)
    if np.random.random() < P: return True
    else: return False


####################
## Code execution ##
####################

# ------------------------- #
#  Set up the two lattices  #
# ------------------------- #
# Shuffle positions of one if both lattice types are the same
if alpha_type == beta_type: shuffle = True
else: shuffle = False

# Supercells - these are ASE Atoms objects
disp = np.zeros((Natoms,3))
atoms_alpha = build_supercell(disp, alpha_vec, alpha_a, alpha_type)
atoms_beta  = build_supercell(disp, beta_vec, beta_a, beta_type, shuffle)


# ------------------- #
#  Initialise things  #
# ------------------- #

# Convert parameters from sweeps to steps
Nsteps = Natoms * sweeps_dF
recalc_step = sweeps_recalc * Natoms
    
hist = [0,]*10
abs_mu_max = len(hist) * bin_width
    
# Track individual energies
Ea_Eb = np.zeros((Nsteps,2))
    
get_index = eq_width

# Get initial lattice energies and difference between them
old_E_alpha, old_E_beta = lattice_energies(atoms_alpha, atoms_beta)
mu = old_E_alpha - old_E_beta
old_index = get_index(mu, 0, abs_mu_max, len(hist))

if s == 0:
    print "Subdomain 0 (-ve mu)"
    ACT = 1 # start with lattice alpha(1) active
else:
    print "Subdomain 1 (+ve mu)"
    ACT = -1 # start with lattice beta(-1) active

step = 1
moves = 0

stuck = 0
prevstep = 1


# ----------- #
#  Main loop  #
# ----------- #
print "Running for %d sweeps..." %sweeps_dF
while step < Nsteps:
        
    # Pick a random particle
    ipart = np.random.randint(0,Natoms)
            
    # Calculate old local energies before trial move
    old_local_energies = local_energy(atoms_alpha, atoms_beta, ipart)
            
    # Pick a random displacement and update lattice positions
    dr = (np.random.random(3)*2.0 - 1.0)*dr_max
    if step < 10000: dr *= 0.01
    dr_alpha = np.dot(dr, atoms_alpha.get_cell())
    dr_beta = np.dot(dr, atoms_beta.get_cell())
    atoms_alpha.positions[ipart,:] += dr_alpha
    atoms_beta.positions[ipart,:] += dr_beta

    # Calculate new local energies after trial move
    new_local_energies = local_energy(atoms_alpha, atoms_beta, ipart)

    # Calculate energy difference
    delta_local_energies = new_local_energies - old_local_energies

    # New energy difference between lattices
    if step % recalc_step == 1: # full calculation
        new_E_alpha, new_E_beta = lattice_energies(atoms_alpha, atoms_beta)
    else: # update with local energy difference
        new_E_alpha = old_E_alpha + delta_local_energies[1]
        new_E_beta = old_E_beta + delta_local_energies[-1]
    mu = new_E_alpha - new_E_beta

    # absolute value of mu
    if s == 0:
        mu *= -1

    # Sanity check - don't want to be below ideal lattice energy!
    if new_E_alpha < E_ideal or new_E_beta < E_ideal:
        error(this_file, "Oh no! Energy is less than that of an Ideal lattice.")

    # Check if mu within 'global' boundaries
    if mu < 0:
        atoms_alpha.positions[ipart,:] -= dr_alpha
        atoms_beta.positions[ipart,:] -= dr_beta
        if step == prevstep:
            stuck += 1
            if stuck == 1000:
                error(this_file, "Simulation has failed to escape ideal lattice region.")
        else:
            prevstep == step
        continue
    while mu > abs_mu_max:
        abs_mu_max += bin_width
        hist.append(0)

    # Calculate 'mu' and from it get index for arrays
    new_index = get_index(mu, 0, abs_mu_max, len(hist))
       
    # Accept/reject trial move with weighted energy difference
    expnt = -B*(delta_local_energies[ACT])
    result = accept_move(expnt)
        
    # Update things accordingly
    if result == True:
        disp[ipart,:] += dr
        old_E_alpha = new_E_alpha
        old_E_beta = new_E_beta
        old_index = new_index
        moves += 1
    else:
        atoms_alpha.positions[ipart,:] -= dr_alpha
        atoms_beta.positions[ipart,:] -= dr_beta

    # Update weights and histogram
    hist[old_index] += 1

    
    step += 1


# End main loop


move_pct = int(float(moves)/step * 100)

print "- Steps: ",step
print "- Moves: ",moves," (%d percent)" %move_pct

hist = np.array(hist)
mu_bins = np.arange(len(hist)) * bin_width

if s == 0:
    mu_bins = mu_bins * -1

output = np.zeros( (len(hist), 2) )
output[:,0] = hist
output[:,1] = mu_bins


# ------------- #
#  Save output  #
# ------------- #
savename = "mu_range_s" + str(s) + ".out"
print "Saving to file: ", savename
np.savetxt(savename, output)

print ""


