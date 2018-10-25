import numpy as np
from sys import argv
from os.path import basename

from params import *
import dynamics as dyn

if algorithm == 'wang_landau': import wang_landau as alg
elif algorithm == 'multicanonical': import multicanonical as alg
elif algorithm == 'transition': import transition_matrix as alg

# Name of this file
this_file = basename(__file__)

if TRAP == True:
    Ns_eff = Ns
else:
    Ns_eff = 1
p_list = np.arange(Np)


rt_sweeps_pmean_alls = []
rt_sweeps_pstderr_alls = []

# Iterate over different subdomains
for s in range(Ns_eff):

    rt_sweeps_allp = []

    # Iterate over different processes operating in the same subdomain/domain
    for p in p_list[s::Ns_eff]

        # Get file names to load
        input_files = dyn.file_names('output', s, p)

        # Load input files
        this_rt_data = np.loadtxt(input_files['rt'])
        
        # Round trip time in sweeps
        rt_sweeps_allp.append(this_rt_data[3])

    # End loop over processes

    # Mean and standard error for this subdomain
    rt_sweeps_pmean = np.mean(rt_sweeps_allp)
    rt_sweeps_pstderr = np.std(rt_sweeps_allp) / np.sqrt(len(rt_sweeps_allp))
    
    # Append to lists containing data from all subdoms
    rt_sweeps_pmean_alls.append(rt_sweeps_pmean)
    rt_sweeps_pstderr_alls.append(rt_sweeps_pstderr)

# End loop over effective subdomains

# Print round trip time info for get_data.sh script
if argv[1] == 'rt':
    printlist = [simID, Ns] + rt_sweeps_pmean_alls + rt_sweeps_pstderr_alls
    for item in printlist:
        print item,


