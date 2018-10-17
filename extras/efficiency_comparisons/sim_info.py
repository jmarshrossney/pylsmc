"""
Module.

Each simulation has an integer 'simID' parameter set in it's params.py file.
Simulations with the same simID will be averaged over in these comparison scripts,
whereas simulations with different simID's will be compared.

Use this module to add appropriate descriptors for the simID's, so that the plots make
sense rather than just comparing simulations labelled by integers.
"""

# Label to go on categorical axis
axis_label = "Number of bins"

# dictionary of strings describing input files (plot legend)
file_labels = {'rtd.txt': 'disc. weights',
               'wld.txt': 'disc. weights',
               'dFd.txt': 'disc. weights',
               'stdevd.txt': 'disc. weights',
               'rti.txt': 'interp. weights',
               'wli.txt': 'interp. weights',
               'dFi.txt': 'interp. weights',
               'stdevi.txt': 'interp. weights'}

# List of strings describing type of simulations, indexed by simID (categorical axis)
sim_labels = ['301', '251', '201', '151', '101']

