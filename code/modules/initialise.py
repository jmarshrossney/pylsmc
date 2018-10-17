"""
Module.

Initialises and does sanity checks on two proposed lattices with which to perform lattice switching.

Links processors to subdomains and sets up appropriate naming system.

"""

import numpy as np
import math as m
from os.path import exists, basename
import sys
from ase.lattice.cubic import FaceCenteredCubic, BodyCenteredCubic
from ase.lattice.hexagonal import HexagonalClosedPacked

# Import params from parent directory
path_to_params = '..'
sys.path.insert(0, path_to_params) # looks in parent directory first
from params import *
sys.path.remove(path_to_params)

# Import lsmc modules from this directory
import domain as dom

# Name of this file
this_file = basename(__file__)


############################
## Lattice initialisation ##
############################
def load_shuffled_indices():
    """ Attempts to load file containing order by which atoms have been shuffled.
        Returns an error if it doesn't exist. """
    # If shuffle file exists, load and return array
    input_file = "shuffled_indices.out"
    if exists(input_file) == True:
        shuffled_indices = np.loadtxt(input_file, dtype=np.int)
    
    # If not, the user needs to run check_params.py
    else:
        error(this_file, "File %s not found. Please run check_params.py to generate this file." %input_file)
    
    return shuffled_indices


def build_supercell(disp, vec, Ltype, shuffle=False):
    """ Build a supercell using a given set of lattice vectors, lattice constant and type"""

    if Ltype == 'fcc':
        build = FaceCenteredCubic
        Lconst = m.pow(4.0/density, 1.0/3.0)

    elif Ltype == 'bcc':
        build = BodyCenteredCubic
        Lconst = m.pow(2.0/density, 1.0/3.0)

    elif Ltype == 'hcp':
        build = HexagonalClosedPacked
        a = m.pow(m.sqrt(2.0)/density, 1.0/3.0)
        Lconst = (a, m.sqrt(8.0/3.0)*a)

    else:
        error(this_file, "Only fcc, bcc and hcp supported right now. Please set alpha_type and beta_type to 'fcc', 'bcc' or 'hcp' in params.py")

    atoms = build(size=tuple(vec),
                    symbol='Ar', # inconsequential
                    pbc=(1,1,1),
                    latticeconstant=Lconst)
    
    # Translate so that center of mass is at (0,0,0)
    com = atoms.get_center_of_mass()
    atoms.translate(-com)

    # Shuffle atomic indices of one lattice if both lattices are same type
    if shuffle != False:
        shuffled_indices = load_shuffled_indices()
        pos = atoms.positions
        atoms.set_positions(pos[shuffled_indices])

    # Set the atomic positions according to disp
    ideal_positions = atoms.get_positions()
    new_positions = ideal_positions + np.dot(disp, atoms.get_cell())
    atoms.set_positions(new_positions)

    return atoms
        


###########################################
## Processors, subdomains, naming system ##
## General to all algorithm types        ##
###########################################

def map_proc_to_subdom(p):
    """ For a processor index p, returns the appropriate subdom index s """
    if TRAP == True:
        if Np > Ns: # more processors than subdomains
            s = p % Ns
        else: # 1 - 1 mapping processors to subdomains
            s = p
    else:
        s = p % Ns # global, but initialise processes evenly across subdomains
    return s


def get_size(s):
    """ Returns the size of binned data in a given subdomain """
    if Ns == 1:
        size = bins[0]
    else:
        if TRAP == False:
            size = np.sum(bins)
        else:
            size = dom.subdom[s]['bins']
    return size

def naming_system(s, p):
    """ Sets up the appropriate file naming system for a given processor and subdomain index """
    if Ns == 1:
        if Np == 1:
            pname = ('', '')
            sname = ('', '')
        else:
            pname = ('_p'+str(p), '_pcomb')
            sname = ('', '')
    else:
        if TRAP == False:
            pname = ('_p'+str(p), '_pcomb')
            sname = ('', '')
        else:
            if Np > Ns:
                pname = ('_p'+str(p), '_pcomb')
            else:
                pname = ('', '')
            sname = ('_s'+str(s), '_scomb')
    return sname, pname


def rename_inputs(fname, sname, pname):
    """ Just adds the '_pY' '_sX' onto file input names """
    if fname == None:
        return fname
    else:
        stringsplit = fname.split('.')
        stringsplit[0] = stringsplit[0] + pname + sname
        return stringsplit[0] + '.' + stringsplit[1]


def file_input(input_file, shape):
    """ Reads in data from a file or, if 'None' argument, initialises an empty array """

    input_zeros = np.zeros(shape)
    correct_shape = np.shape(input_zeros)

    if input_file == None:
        return input_zeros

    else:
        if exists(input_file) == False:
            print "Input file %s not found. Initialising with zeros." %input_file
            return input_zeros
        
        else:    
            print "Loading input file: ", input_file

            input_data = np.loadtxt(input_file)
            if np.shape(input_data) != correct_shape:
                errmsg = "Expected input file "+str(input_file)+" to be shape "+\
                        str(correct_shape)+",\nbut had shape "+str(np.shape(input_data))
                error(this_file, errmsg)

            return input_data




