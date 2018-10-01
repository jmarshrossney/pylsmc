"""
Module.

Contains some stuff nicked from domain.py and initialise.py which is needed by find_domain.py
"""


import numpy as np
import math as m
from sys import exit
from os.path import basename

from ase.lattice.cubic import FaceCenteredCubic, BodyCenteredCubic

from params import *

# Name of this file
this_file = basename(__file__)


############################
## Lattice initialisation ##
############################
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
        

##############
## Indexing ##
##############
def eq_width(mu, mu_min, mu_max, N_bins):
    index = m.floor((mu-mu_min)/(mu_max-mu_min) * N_bins)
    return int(index)

