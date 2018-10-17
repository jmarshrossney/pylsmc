"""
Script.

Takes a time-series of free energy differences and computes the mean, standard deviation and standard error on the mean.
Then, does various things depending on the value of argv[1], which can be one of the following:
    
    eval : evaluate (True/False) whether simulation has converged to desired tolerance (called automatically by run_parallel.sh)

    report : Compute mean, std dev and std err within a moving window, for all timesteps > the number of timesteps in the window.
             Conveniently print these values for use by some scripts.

    plot : Same as above but plot results instead of printing.

argv[2] is an optional argument which gives the size of the 'window' - the number of data points used to calculate the estimators. 
This must be smaller than the total number of data points.
The default is in params.py

Can also use argv[3] to provide a different data file from the one given in params.py

Run as      python calc_dF.py (window) (input_file)

"""

import numpy as np
import matplotlib.pyplot as plt
from sys import argv, exit
from os.path import basename

from params import *

# Name of this file
this_file = basename(__file__)

# Can specify a non-default window size via argument variable
if len(argv) > 2:
    window = int(argv[2])

    # Can specify non-default data file
    if len(argv) > 3:
        deltaF_file = argv[3]


######################################################
## Functions to compute mean and standard deviation ##
######################################################
'''
########### Not convinced that this is working! ##########
def weighted_estimators(deltaF_series, error_series):
    """ Takes into account unequal errors on individual dF's due to join. 
        Formulae from NIST datasheets """

    N = len(deltaF_series)

    # Compute weighted mean
    num = np.sum( deltaF_series / error_series**2 )
    denom = np.sum( 1.0 / error_series**2 )
    weighted_mean = num/denom
    err_weighted_mean = np.sqrt( 1.0 / np.sum( 1.0 / error_series**2 ) )

    # Compute weighted standard deviation
    num = np.sum( (deltaF_series - weighted_mean)**2 / error_series**2 )
    denom = denom * (N-1) / N
    weighted_var = num/denom
    weighted_stdev = np.sqrt(weighted_var)

    return weighted_mean, weighted_stdev, err_weighted_mean
'''


def unweighted_estimators(deltaF_series):
    """ Faster, assumes errors on individual dF's are equal """
    
    N = len(deltaF_series)

    unweighted_mean = np.mean(deltaF_series)
    unweighted_stdev = np.std(deltaF_series)
    err_unweighted_mean = unweighted_stdev / np.sqrt(N)
    
    return unweighted_mean,  unweighted_stdev, err_unweighted_mean


###################################################
## Evaluate whether the simulation has converged ##
###################################################
if argv[1] == 'eval':
    # Called by run_parallel.sh automatically
    # intended for use during a simulation rather than after.

    # Load the series data
    input_data = np.loadtxt(deltaF_file)

    # How many iterations of 'sweeps' have we done so far
    Niters = len(input_data)

    # If a target number of iterations has been set
    # Handled by the run script
    if iterations > 0:
        print 'True'

    # If the simulation is to continue until stdev_converged
    # (iterations < 0 in params)
    else:

        # Don't bother if there's not enough points
        if len(Niters) < window:
            print 'False' # print statement picked up by run_parallel.sh
            exit(0) # exit happily

        # Unpack
        deltaF_series = input_data[:,0]
        error_series = input_data[:,1]
    
        # Only use the last 'window' data points
        deltaF_series = deltaF_series[-window:]
        error_series = error_series[-window:]
        
        # Weighted estimators used if there are join errors
        # Currently switched off because
        # (a) doubts with weighted_estimators() and 
        # (b) join errors don't seem to be significant.
        if False: #TRAP == True and Ns > 1:
            mean, stdev, stderr = weighted_estimators( \
                        deltaF_series, error_series)

        # Unweighted estimators if no join errors -> all dF errors equal
        else:
            mean, stdev, stderr = unweighted_estimators( \
                        deltaF_series)

        # Print statement picked up by run_parallel.sh
        if stdev < stdev_converged:
            print 'True'
        else:
            print 'False'


###########################################
## Compute estimators in a moving window ##
###########################################
if argv[1] in ('plot', 'report'):
    # Either plot or conveniently print resuts.

    # Load the series data
    input_data = np.loadtxt(deltaF_file)
    
    # Unpack
    deltaF_series = input_data[:,0]
    error_series = input_data[:,1]
    
    # Cumulative sweeps corresponding to each dF value
    sweeps_series = np.arange(1,len(deltaF_series)+1)*sweeps_dF

    # Size of arrays to take computed values of estimators
    N_computed = len(deltaF_series) - window

    if N_computed < 1:
        error(this_file, "The number of data points is less than or equal to 'window'. Either generate more data,\n or specify a smaller window size as the second argument variable when running this script.")

    # Initialise arrays
    mean_series = np.ones(N_computed)
    stdev_series = np.ones(N_computed)
    stderr_series = np.ones(N_computed)
    
    # -------------------------------------------- #
    #  Iterate over positions of moving dF window  #
    # -------------------------------------------- #
    for i in range(N_computed):
        
        # Weighted estimators used if there are join errors
        # Currently switched off because
        # (a) doubts with weighted_estimators() and 
        # (b) join errors don't seem to be significant.
        if False: # TRAP == True and Ns > 1:
            mean, stdev, stderr = weighted_estimators( \
                deltaF_series[i:i+window], error_series[i:i+window])
        
        # Unweighted estimators if no join errors -> all dF errors equal
        else:
            mean, stdev, stderr = unweighted_estimators( \
                    deltaF_series[i:i+window])

        # Add results from this window to arrays        
        mean_series[i] = mean
        stdev_series[i] = stdev
        stderr_series[i] = stderr


    # --------------------------- #
    #  Plot a pretty time series  #
    # --------------------------- #
    if argv[1] == 'plot':
        
        plt.rc('text',usetex=True)
        font = {'family' : 'serif',
                'size' : 16}
        plt.rc('font', **font)
        plt.rcParams['axes.linewidth'] = 2

        # Print final values
        print "Mean: %1.2e" %mean
        print "Error on the mean: %1.2e" %stderr
        print "Std dev: %1.2e" %stdev
        
        # Also want to plot mean +/- 1 standard deviation
        dF_up = mean_series + stdev_series
        dF_down = mean_series - stdev_series
       
        fig, ax = plt.subplots()
        ax.set_title("Convergence of free energy differences")
        ax.set_xlabel("Sweeps")
        ax.set_ylabel("Free energy difference (kT/atom)")
        ax.tick_params(direction='in',top=True,right=True)
        #ax.set_yticklabels([])
        #ax.set_xticklabels([])

        ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

        ax.plot(sweeps_series,np.zeros(len(deltaF_series)),'r--',label=r"$\Delta F = 0$")
        #ax.set_xlim(0,np.max(sweeps_series)*1.4)
        #ax.set_ylim(np.min(deltaF_series)*1.2, np.max(deltaF_series)*3.2)
        ax.plot(sweeps_series[:-1:1],deltaF_series[:-1:1],'ko',markersize=4,label=r"raw $\Delta F$")
        ax.plot(sweeps_series[-1],deltaF_series[-1],'ro',markersize=4)

        #ax.plot(sweeps_series[window:],mean_series,'g',label=r"$\bar{\Delta F}$")
        #ax.plot(sweeps_series[window:],dF_up,'m',label=r"$\Delta F \pm 1\sigma$")
        #ax.plot(sweeps_series[window:],dF_down,'m')

        ax.legend()
        plt.tight_layout()
        plt.show()


    # -------------------------------------- # 
    #  Print/save things for use by scripts  #
    # -------------------------------------- #
    elif argv[1] == 'report':
    
        # Print final values
        # Useful for comparison between simulations run for same Nsweeps
        print simID, Ns, sweeps_series[-1], mean, stdev, stderr

        # Save a time series of standard deviations
        save_data = np.zeros( (len(stdev_series)+1, 2) )
        save_data[0,0] = simID
        save_data[0,1] = Ns
        save_data[1:,0] = sweeps_series[window:]
        save_data[1:,1] = stdev_series
        np.savetxt("stdev_series.txt", save_data)
        


