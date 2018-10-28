"""
Script.

Takes various quantities output by a simulation and condenses them into a few measures of simulation efficiency,
which can then be compared with simulations with different parameters.

Can also be used to create plots for easy visualision of simulation efficiency.

argv[1] is plot or report

if argv[1] is report, argv[2] should be an array index corresponding to the F value you want to report on, e.g. 0 for the first, -1 for the last...

File names are automatic (has to be so to allow calling from script)
These are "sweeps_allF_sX.out"

Run as e.g.      python convergence_wl.py plot
                 python convergence_wl.py report -1

"""

import numpy as np
import math as m
import matplotlib.pyplot as plt
import sys
from sys import argv

# Import params from parent directory
path_to_params = '..'
sys.path.insert(0, path_to_params) # looks in parent directory first
from params import *
sys.path.remove(path_to_params)


# File names
input_file_stem = "sweeps_allF_s"

# Load first / only file
input_data = np.loadtxt(input_file_stem+"0.out")
sweeps_all_s = []
cum_sweeps_all_s = []

# Load other files if there are > 1 subdomains
for s in range(Ns):
    input_file = input_file_stem+str(s)+".out"
    input_data = np.loadtxt(input_file)
    sweeps_all_s.append(input_data[:,1])

    # Cumulative sweeps for this subdomain
    cum_sweeps_all_s.append(np.cumsum(input_data[:,1]))

# F factor is the first column, same for all subdomains
F = input_data[:,0]

# Should usually be equal to Ns in params
Ns_input = len(sweeps_all_s)

# Need the transpose later: (row, col) = (sweeps F, subdomain)
cum_sweeps_all_s_T = np.array(cum_sweeps_all_s).T


# ------------------------------------------------ #
#  Compute sweeps per F, averaged over subdomains  #
# ------------------------------------------------ #
# Initialise arrays for mean + standard err
sweeps_smean = np.zeros(len(F))
err_sweeps_smean = np.zeros(len(F))
cum_sweeps_smean = np.zeros(len(F))
err_cum_sweeps_smean = np.zeros(len(F))

# Iterate over different F factors
for f in range(len(F)):

    # Create Ns-d array of (cumulative) sweeps for this F
    sweeps_this_f = np.zeros(Ns_input)
    cum_sweeps_this_f = np.zeros(Ns_input)
    
    for s in range(Ns_input):
        sweeps_this_f[s] = sweeps_all_s[s][f]
        cum_sweeps_this_f[s] = cum_sweeps_all_s[s][f]

    # Mean and std deviation over subdomains
    sweeps_smean[f] = np.mean(sweeps_this_f)
    err_sweeps_smean[f] = np.std(sweeps_this_f) / np.sqrt(Ns_input)
    cum_sweeps_smean[f] = np.mean(cum_sweeps_this_f)
    err_cum_sweeps_smean[f] = np.std(cum_sweeps_this_f) / np.sqrt(Ns_input)


####################################################################
## Simply save array of mean sweeps per F, for use by wl_to_dF.sh ##
####################################################################
if argv[1] == 'save':
    output_file = "sweeps_allF_smean.out"
    print "Saving sweeps to ", output_file
    np.savetxt(output_file, sweeps_smean)
    

##############################################
## Print for comparison between simulations ##
##############################################
if argv[1] == 'report':
    
    # Print the total number of sweeps required to converge to the final F for each subdom
    X = int(argv[2])
    printlist = [simID, "%.10f"%F[X], Ns] + list(cum_sweeps_all_s_T[X,:])
    for item in printlist:
        print item,


###################################################
## Plots to show work balance in this simulation ##
###################################################
elif argv[1] == 'plot':
    subdom = np.arange(Ns_input)

    # Plot sweeps per F
    fig, ax = plt.subplots()
    ax.set_xlabel("$F-1$")
    ax.set_ylabel("Sweeps")
    ax.set_title("Sweeps per F, averaged over subdomains")
    
    # Plot cumulative sweeps per F
    fig2, ax2 = plt.subplots()
    ax2.set_xlabel("$F-1$")
    ax2.set_ylabel("Cumulative sweeps")
    ax2.set_title("Cumulative sweeps per F, averaged over subdomains")

    # Plot each subdomain separately
    for s in range(Ns_input):
        ax.semilogx(F-1, sweeps_all_s[s], 'o', label="Subdomain %d" %s)
        ax2.semilogx(F-1, cum_sweeps_all_s[s], '--', label="Subdomain %d" %s)
    
    # Plot mean and standard deviation of subdomains
    ax.errorbar(F-1, sweeps_mean, yerr=sweeps_stdev, fmt='ks', label="Mean +/- stdev")
    fig.gca().invert_xaxis()
    ax2.errorbar(F-1, cum_sweeps_mean, yerr=cum_sweeps_stdev, fmt='ks', label="Mean +/- stdev")
    fig2.gca().invert_xaxis()
        
    ax.legend(loc=2)
    ax2.legend(loc=2)


    # Plot cumulative sweeps for each F, each subdomain
    plt.figure(3)
    plt.xlabel("Subdomain index")
    plt.ylabel("Cumulative sweeps")
    plt.title("Cumulative sweeps per F")
    
    # Plot for each F factor
    for f in range(len(F)):
        Fstr = str( round(F[f]-1, -int(m.floor(m.log10(abs(F[f]-1)))) ) + 1)
        plt.plot(subdom, cum_sweeps_all_s_T[f,:], label="F="+Fstr)
    plt.legend(loc=3)

    plt.show()


