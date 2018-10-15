"""
Module.

Each simulation has an integer 'simID' parameter set in it's params.py file.
Simulations with the same simID will be averaged over in these comparison scripts,
whereas simulations with different simID's will be compared.

Use this module to add appropriate descriptors for the simID's, so that the plots make
sense rather than just comparing simulations labelled by integers.
"""

# Label to go on categorical axis
axis_label = "Number of subdomains"

# List of strings describing type of simulations, indexed by simID
sim_labels = ['1 subdomain',
             '3 subdomains',
             ]

