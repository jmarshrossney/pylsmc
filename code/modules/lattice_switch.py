"""
Module.

Contains 'run' which executes lattice-switch Monte Carlo algorithm, calling low-level energy calculations at each timestep.

Also contains 'drift' - a function which moves atoms around with Boltzmann weights until the system drifts into an appropriate subdomain.

"""

import numpy as np
import math as m
from os.path import basename

from params import *
import energy
import domain as dom
import lsmc_dynamics as dyn

# Import module corresponding to the type of algorithm being used
if algorithm == 'wang_landau': import wang_landau as alg
elif algorithm == 'multicanonical': import multicanonical as alg
elif algorithm == 'transition': import transition_matrix as alg

# Name of this file
this_file = basename(__file__)


#####################
## Local functions ##
#####################

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

def attempt_switch(mu, ACT, switches):
    """ Attempt switch move """
    expnt = ACT*B*mu
    result = accept_move(expnt)
    if result == True:
        ACT *= -1
        switches += 1
    return ACT, switches
    
def coinflip():
    """ Binary choice """
    return np.random.choice([0,1])

# -------------------- #
# Weighting strategies #
# -------------------- #
def interpolated_weight(mu, index, mu_bins, weights):
    """ Return interpolated weight """
    if mu > mu_bins[index]:
        ilow = index
        ihigh = index + 1
    else:
        ilow = index - 1
        ihigh = index
    frac = (mu-mu_bins[ilow]) / (mu_bins[ihigh]-mu_bins[ilow])
    grad = weights[ihigh] - weights[ilow]
    return weights[ilow] + frac*grad

def discrete_weight(mu, index, mu_bins, weights):
    """ Simply returns weight in current bin """
    return weights[index]

def extrapolate(array):
    """ Adds linearly extrapolated end bins to a numpy array, assuming the
        array has spaces at either end for this extrapolated value to go """
    array[-1] = 2*array[0] - array[1] # this is really to the left of bin 0
    array[-2] = 2*array[-3] - array[-4]
    return



########################################################################################
##                                                                                    ##
##                Main routine: runs the Lattice-switch Monte Carlo with              ##
##                   Wang Landau, multicanonical or transition matrix.                ##
##                                Called from main.py                                 ##
##                                                                                    ##
########################################################################################


def run(atoms_alpha, atoms_beta, disp, binned_data, dyn_data, F, p, s):
    
    # ------------------------------------------------- #
    #  Set simulation parameters and useful quantities  #
    # ------------------------------------------------- #
    lendata = len(binned_data['w']) # length of weights
    kTlogF = kT*np.log(F)

    # Convert frequency parameters from sweeps to steps
    Nsteps = alg.Nsteps
    save_step = sweeps_save * Natoms
    series_step = int(sweeps_series * Natoms)
    recalc_step = sweeps_recalc * Natoms
    refresh_step = sweeps_refresh * Natoms

    # Set boundaries for this simulation
    if TRAP == True: # restrict to subdomain 's'; local bin indices
        global_mu_min = dom.subdom[s]['min']
        global_mu_max = dom.subdom[s]['max']
        get_index = dom.get_local_index
        mu_bins = dom.get_local_mu_bins(s)
    else: # can access entire domain; global bin indices
        global_mu_min = boundaries[0]
        global_mu_max = boundaries[-1]
        get_index = dom.get_global_index
        mu_bins = dom.get_global_mu_bins()

    # Set weights function
    if use_interpolated_weights == True:
        get_weight = interpolated_weight
        
        # Need to add bins to the edges of mu and weights for extrapolation
        mu_bins = np.concatenate( (mu_bins, [0,0]) )
        extrapolate(mu_bins)
        binned_data['w'] = np.concatenate( (binned_data['w'], [0,0]) )
        extrapolate(binned_data['w'])
        
        # Need to do the same for cuts since it uses mu_bins
        if track_dynamics == True:
            dyn_data['cu'] = np.concatenate( (dyn_data['cu'], [0,0]) )

    else:
        get_weight = discrete_weight


    # ------------------------ #
    #  Initialise usual stuff  #
    # ------------------------ #
    # Get initial lattice energies and difference between them
    old_E_alpha, old_E_beta = lattice_energies(atoms_alpha, atoms_beta)
    old_mu = old_E_alpha - old_E_beta
    new_mu = old_mu 
    
    # For a global simulation, the initial subdomain (corresponding to the
    # input disp vector) needs to be computed
    if TRAP != True:
        s = dom.get_subdomain(old_mu)
   
    # Initial index and weight
    old_index = get_index(old_mu, s)
    old_weight = get_weight(old_mu, old_index, mu_bins, binned_data['w'])
    
    # Active lattice
    if s < dom.s_cross: # start with lattice alpha(1) active
        ACT = 1
    elif s > dom.s_cross: # start with lattice beta(-1) active
        ACT = -1
    else: # randomise
        ACT = np.random.choice([-1,1])
    
    # Flag for histogram flatness
    FLAT = 1 # start with 1 < f if wang-landau, otherwise f = 1 = FLAT
    
    # Initialise counters
    step = 1
    moves = 0
    retakes = 0
    switches = 0
    

    # ---------------------------- #
    #  Initialise 'dynamics' info  #
    # ---------------------------- #
    if track_dynamics == True:
        
        cuts = {'counter': 0,
                'old_mu': old_mu,
                'new_mu': old_mu}
        
        if TRAP == True:
            rt_min = 0
            rt_max = lendata-1
        else:
            rt_min = np.argmax( binned_data['w'][ mu_bins<0 ] )
            rt_max = np.argmax( binned_data['w'][ mu_bins>0 ] ) + len(mu_bins[ mu_bins<0 ])
        
        rtrips = {'minmax': (0, rt_min, rt_max),
                  'counts': 0,
                  'flag': int(dyn_data['rt'][4])}
                    
        # Mini transition matrix only works for TRAP == True (due to indexing)
        if TRAP == True:
            minimat = {'counter': 0,
                       'root_index': old_index,
                       'old_index': dyn.get_minibin_index(old_mu, s, old_index),
                       'Pprod': 1}

        # Open files for edge dmu's
        dyn_output_files = dyn.file_names('output', s, p)
        l_edge_dmu = open(dyn_output_files['led'], 'a')
        r_edge_dmu = open(dyn_output_files['red'], 'a')


    # --------------------------- #
    #  Open file for data series  #
    # --------------------------- #
    if track_series == True:
        series_file = alg.file_names('output', s, p)['s']
        series = open(series_file, 'a')

   
    # ------------------------------------------------------------ #
    #  Main loop over individual particle moves + switch attempts  #
    # ------------------------------------------------------------ #
    while FLAT < F or step <= Nsteps: # Flat condition only relevant for Wang Landau
        
        # Refresh weights with eigenvector of transition matrix (TM only)
        if step % refresh_step == 0:
            alg.refresh_func(binned_data)
            
        # Write values from this step to the series file. Note high precision needed for mu
        if track_series == True:
            if step % series_step == 0:
                E_active = (0,old_E_alpha,old_E_beta)[ACT]
                series.write("%d %.10f %f %f \n" \
                        %(s, old_mu, old_weight, E_active) )


        # --------------------------------------------------- #
        #  Move a single particle and compute energy changes  #
        # --------------------------------------------------- #
        # Pick a random particle
        ipart = np.random.randint(0,Natoms)
            
        # Calculate old local energies for the two lattices before trial move
        old_local_energies = local_energy(atoms_alpha, atoms_beta, ipart)
            
        # Pick a random displacement and update lattice positions
        dr = (np.random.random(3)*2.0 - 1.0)*dr_max
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
            
            # Sanity check - don't want to be below ideal lattice energy!
            if new_E_alpha < E_ideal or new_E_beta < E_ideal:
                error(this_file, "Oh no! Energy is less than that of an Ideal lattice.")
        
        else: # update with local energy difference
            new_E_alpha = old_E_alpha + delta_local_energies[1]
            new_E_beta = old_E_beta + delta_local_energies[-1]
        
        new_mu = new_E_alpha - new_E_beta

        
        # ---------------- #
        #  Switch attempt  #
        # ---------------- #
        # Attempt a switch here 50% of the time - before any move evaluation
        coin = coinflip()
        ACT, switches = ( (ACT, switches),
                attempt_switch(new_mu, ACT, switches) )[coin]
        

        # --------------------------------------- #
        #  Check if boundaries have been crossed  #
        # --------------------------------------- #
        # 'Global' boundaries are those that the simulation may not cross
        if new_mu < global_mu_min or new_mu > global_mu_max:

            # Update the bin from which the move was attempted
            old_weight, FLAT = alg.update_func(
                        step, binned_data, delta_local_energies[ACT],
                        old_index, old_index, old_index, # only update edge bin
                        old_mu, mu_bins,
                        old_weight, get_weight, lendata, kTlogF, s)
            
            if track_dynamics == True:
                # Record size of attempted move out of (sub)domain
                dmu = abs(new_mu - old_mu)
                if new_mu < global_mu_min:
                    l_edge_dmu.write("%f\n" %dmu)
                elif new_mu > global_mu_max:
                    r_edge_dmu.write("%f\n" %dmu)
                
                # Update data on simulation dynamics / mobility
                dyn.update_func(
                        dyn_data, mu_bins, old_mu, old_index, s,
                        cuts, rtrips,
                        False, dmu, abs(delta_local_energies[ACT]))
              
                if TRAP == True:
                    dyn.update_minimat(dyn_data, mu_bins, old_mu, s,
                            minimat, True, 0)

            # Undo lattice positions update
            atoms_alpha.positions[ipart,:] -= dr_alpha
            atoms_beta.positions[ipart,:] -= dr_beta
            retakes += 1
            step += 1
            
            # Go back to the start of the while loop
            continue

        # Check if subdomain boundary crossed
        old_s = s
        if new_mu < dom.subdom[s]['min']:
            s -= 1
        elif new_mu > dom.subdom[s]['max']:
            s += 1

        
        # ---------------------- #
        #  Evaluate move result  #
        # ---------------------- #
        # Get index for arrays using mu
        new_index = get_index(new_mu, s)

        # Difference in weights augments delta_local_energy
        new_weight = get_weight(new_mu, new_index, mu_bins, binned_data['w'])
        augment = new_weight - old_weight

        # Accept/reject trial move with weighted energy difference
        expnt = -B*(delta_local_energies[ACT] + augment)
        accepted = accept_move(expnt)
       
        """ #Difference between bin-averaged weight and interpolated
        intra_bin_weight = old_weight - weights[old_index]
        #sampling_adjust = np.exp( B*intra_bin_weight ) # doesn't work: centre of bins?"""
        
        old_index_copy = old_index # to update Cmat and extra info after move evaluation
        dmu = abs(new_mu-old_mu) # to update extra info

        # Update things according to whether move accepted
        if accepted == True:
            disp[ipart,:] += dr
            old_E_alpha = new_E_alpha
            old_E_beta = new_E_beta
            old_index = new_index
            old_weight = new_weight
            old_mu = new_mu
            moves += 1

            cuts_new_mu = new_mu

        else:
            atoms_alpha.positions[ipart,:] -= dr_alpha
            atoms_beta.positions[ipart,:] -= dr_beta
            s = old_s # Account for possible subdomain change

            cuts_new_mu = old_mu


        # ------------------------ #
        #  Update bin information  #
        # ------------------------ #
        old_weight, FLAT = alg.update_func(
                    step, binned_data, delta_local_energies[ACT],
                    new_index, old_index, old_index_copy,
                    old_mu, mu_bins,
                    old_weight, get_weight, lendata, kTlogF, s)


        # ------------------------------------ #
        #  Update data on simulation dynamics  #
        # ------------------------------------ #
        if track_dynamics == True:
            
            # Update data on simulation dynamics / mobility
            cuts['new_mu'] = cuts_new_mu
            dyn.update_func(
                    dyn_data, mu_bins, old_mu, old_index, s,
                    cuts, rtrips,
                    accepted, dmu, abs(delta_local_energies[ACT]))
            
            if TRAP == True:
                dyn.update_minimat(dyn_data, mu_bins, old_mu, s,
                        minimat, accepted, expnt)

        # ---------------- #
        #  Switch attempt  #
        # ---------------- #
        # 50/50 chance of attempting a switch after the move
        ACT, switches = ( attempt_switch(old_mu, ACT, switches),
                (ACT, switches) )[coin]

        
        # ------------------- #
        #  Save periodically  #
        # ------------------- #
        if step % save_step == 0:
            alg.save_func(binned_data, mu_bins, step, s, p)

            # Take this opportunity to reset centre of mass
            atoms_alpha.positions = atoms_alpha.positions - atoms_alpha.get_center_of_mass()
            atoms_beta.positions = atoms_beta.positions - atoms_beta.get_center_of_mass()

        step += 1


    # End main loop
 

    print "- Steps: ", step
    print "- Moves: ", moves
    print "- Retakes: ", retakes
    print "- Switches: ", switches

    # --------------------------------------------- #
    #  Finalise things before returning to main.py  #
    # --------------------------------------------- #
    # Close series file
    if track_series == True:
        series.close()

    if track_dynamics == True:
        # Close the edge dmu files
        l_edge_dmu.close()
        r_edge_dmu.close()

        # Round trip output file
        dyn_data['rt'][0] += step
        dyn_data['rt'][1] += switches
        dyn_data['rt'][2] += rtrips['counts'] # half round trips
        rt_rate = float(0.5 * dyn_data['rt'][2] * Natoms) / dyn_data['rt'][0] # rtrips per sweep
        dyn_data['rt'][3] = 1.0 / rt_rate # sweeps per rtrip
        dyn_data['rt'][4] = rtrips['flag'] # to be picked up by next iteration

    # Compute eigenvector and refresh weights
    alg.refresh_func(binned_data)

    # If interpolated weights used, remove those extrapolated bins / 'ghost points'
    if len(binned_data['w']) != len(binned_data['h']):
        binned_data['w'] = binned_data['w'][:-2]
        if track_dynamics == True:
            dyn_data['cu'] = dyn_data['cu'][:-2]
    
    # Only care about difference in weights, set minimum to zero
    binned_data['w'] -= np.min(binned_data['w'])


    return step



#########################################
##  Drift into a particular subdomain  ##
#########################################

def drift(atoms_alpha, atoms_beta, disp, target_s):

    if target_s < dom.s_cross:
        ACT = 1 # start with lattice alpha (fcc,negative mu) active
    else:
        ACT = -1 # start with lattice beta (bcc, positive mu) active

    recalc_step = sweeps_recalc * Natoms

    # Get initial lattice energies and difference between them
    old_E_alpha, old_E_beta = lattice_energies(atoms_alpha, atoms_beta)
    mu = old_E_alpha - old_E_beta
   
    step = 0
    while (dom.subdom[target_s]['min'] < mu < dom.subdom[target_s]['max']) == False:

        
        # Pick a random particle
        ipart = np.random.randint(0,Natoms)
            
        # Calculate old local energies before trial move
        old_local_energies = local_energy(atoms_alpha, atoms_beta, ipart)
            
        # Pick a random displacement and update lattice positions
        dr = (np.random.random(3)*2.0 - 1.0)*dr_max
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

        # Only accept moves in the right direction
        # Drift to more negative mu: old_mu-new_mu > 0, ACT=1
        # Drift to more positive mu: old_mu-new_mu < 0, ACT=-1
        # Therefore moves with (old_mu-new_mu)*ACT < 0 should be rejected
        if ((old_E_alpha - old_E_beta) - (new_E_alpha - new_E_beta))*ACT < 0:
            atoms_alpha.positions[ipart,:] -= dr_alpha
            atoms_beta.positions[ipart,:] -= dr_beta
            step += 1
            continue

        # Accept/reject trial move
        expnt = -B*(delta_local_energies[ACT])
        result = accept_move(expnt)

        # Update things accordingly
        if result == True:
            disp[ipart,:] += dr
            old_E_alpha = new_E_alpha
            old_E_beta = new_E_beta
            mu = old_E_alpha - old_E_beta
        else:
            atoms_alpha.positions[ipart,:] -= dr_alpha
            atoms_beta.positions[ipart,:] -= dr_beta

        step += 1

    # Once in the right subdomain

    return mu


