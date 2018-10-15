"""
Script.

Performs post-processing steps for a series of data:
    Computes decorrelation time;
    Unfolds weights to give a histogram;
    Reweights to different temperatures * (* not yet implemented).

Data is assumed to be in the following format (for each row):
    subdomain  mu  weights  E_active

Some rows will look like -1  -1  -1  -1. These indicate a change to a new, independently initiated simulation. 
Independent simulations are dealt with separately since things like round-trip time and autocorrelation will be miscalculated otherwise.

If multiple simulations have been run in the same subdomain, this script should be run after joining data within subdomains (but before joining subdomains!)
i.e. first run:  python join_procs.py s

No argument variables: Run as   python analyse_series.py

"""

import numpy as np
import matplotlib.pyplot as plt

from params import *
import domain as dom

if algorithm == 'wang_landau': import wang_landau as alg
elif algorithm == 'multicanonical': import multicanonical as alg
elif algorithm == 'transition': import transition_matrix as alg


if TRAP == True:
    Ns_eff = Ns # multiple procs on same subdomain already combined in join_procs
else:
    Ns_eff = 1 # there may be more than 1 subdom, but each process is global - simulation has access to entire domain


# ----------------------------------- #
#  Iterate over different subdomains  #
# ----------------------------------- #
for s in range(Ns_eff):
    
    print ""
    if TRAP == True: print "Series data from subdomain: ", s
    else: print "Series data from global domain."

    # Load the series data
    input_file = alg.file_names('pcomb', s)['s']
    input_data = np.loadtxt(input_file)

    # Unpack
    subdom = np.array(input_data[:,0], dtype=np.int)
    mu = input_data[:,1]
    weights = input_data[:,2]
    E_active = input_data[:,3]

    # Number of bins, mu positions of bins
    if TRAP == True:
        N_bins = dom.subdom[s]['bins']
        mu_bins = dom.get_local_mu_bins(s)
    else:
        N_bins = np.sum(bins)
        mu_bins = dom.get_global_mu_bins()
    
    # Array for unfolded histogram
    hist_unfld = np.zeros(N_bins)
    print "Nbins = ",N_bins

    # Round-trip stuff
    mu_lo = mu_bins[N_bins/4]
    mu_hi = mu_bins[3*(N_bins/4)]
    tracker = -1 # neutral (neither 0 or 1) initial value
    rt_rate_list = []

    # Decorrelation length
    decorr_list = []
    
    # Find all the rows of -1's which denote that the simulation has been restarted
    restart_rows = np.where( subdom == -1 )[0]

    # List of rows representing the last entry for this independent simulation
    # Note that indexing with [:None] indexes up to the final element
    if len(restart_rows) > 1:
        stop_rows = list(restart_rows[1:]) + [None]
    else:
        stop_rows = [None]
   

    # ------------------------------- # 
    #  Iterate over independent runs  #
    # ------------------------------- #
    for r in range(len(restart_rows)):
        print "Independent run no.", r

        ilo = restart_rows[r] + 1
        ihi = stop_rows[r]
        mu_r = mu[ilo:ihi]
        N_samples = len(mu_r)

        # Decorrelation stuff
        mean = np.mean(mu_r)
        auto = np.correlate(mu_r-mean, mu_r-mean, mode='full')
        mid = len(auto)/2
        auto = auto[mid:] # autocorrelation is symmetric
        x0 = np.sign(auto[0])
        x1 = np.sign(auto[1])
        x2 = np.sign(auto[2])
        decorrelated = False

        """
        # Plotting
        plt.figure(1)
        sw = np.arange(len(auto))*sweeps_series
        plt.plot(sw, auto)
        plt.figure(2)
        sw = np.arange(len(mu_r))*sweeps_series
        plt.plot(sw, mu_r)
        plt.show()
        """
        
        # Initialise round trip counter
        rt_count = 0


        # ---------------------------------- #
        #  Iterate over samples in this run  #
        # ---------------------------------- #
        for j in range(len(mu_r)):

            # Check for decorrelation
            if decorrelated == False:
                x3 = np.sign(auto[j+3])
                if x1 != x2 and x0 == x1 and x2 == x3:
                    decorr_list.append(j+2)
                    decorrelated = True
                x0 = x1
                x1 = x2
                x2 = x3
            
            # Check for round trip completion
            if mu_r[j] < mu_lo:
                if tracker == 1: rt_count += 1
                tracker = 0
            elif mu_r[j] > mu_hi: tracker = 1
            
            # Unfold and bin this sample
            I = dom.get_local_index(mu_r[j], subdom[ilo+j])
            hist_unfld[I] = np.exp(B*weights[ilo+j])

            # Reweight to different temperatures
            ###### To be added ######


        # End loop over samples in this run


        # Round trip rate per sample
        rt_rate = float(rt_count) / N_samples
        rt_rate_list.append(rt_rate)


    # End loop over independent runs

    
    # Save unfolded histogram
    savename = alg.file_names('unfld', s)['s']
    print "Saving unfolded hist to ", savename
    np.savetxt(savename, hist_unfld)

    # Convert to sweeps and take average of decorrelation times
    decorr_arr = np.array(decorr_list) * sweeps_series
    decorr_mean = np.mean(decorr_arr)
    decorr_std = np.std(decorr_arr)
    print "Decorrelation time (sweeps) -- Average: %f , Std dev: %f " %(decorr_mean, decorr_std)

    # Convert to sweeps and take average round trip time
    rt_rate_arr = 1.0/(np.array(rt_rate_list)/sweeps_series)
    rt_time_mean = np.mean(rt_rate_arr)
    rt_time_std = np.std(rt_rate_arr)
    print "Round trip time (sweeps) -- Average: %f , Std dev: %f " %(rt_time_mean, rt_time_std)
       



