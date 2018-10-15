/* -*- mode: C */

%define DOCSTRING
"
Implementation of the 6-12 Lennard-Jones potential in (moderately) optimised Fortran2003 code.
This is intended to be used in Monte Carlo projects in which the higher level sampling algorithm
is implemented in Python. 

The Lennard-Jones potential is implemented in the usual reduced units, is truncated (but not 
shifted) at 3.5 length units. Long range corrections assuming uniform radial density beyond
the cut off length are included.

"
%enddef


/* energy.i */
%module(docstring=DOCSTRING) energy
%{
#define SWIG_FILE_WITH_INIT

/* Function prototypes. First arg in typemat triplet used for python variable name */
double compute_lattice_energy(int positions, int d, double *pos, int cellmatrix, int dh1, double *hmatrix, int ilat);
double compute_local_energy(int i, int positions, int d, double *pos, int cellmatrix, int dh1, double *hmatrix, int ilat);      
void compute_neighbour_list(int positions, int d, double *pos, int cellmatrix, int dh1, double *hmatrix, int ilat);

%}

/* Documentation strings */
%feature("docstring","

Computes the total energy of the system of n interacting particles.
The number of particles n is inferred from the size of the positions
array. This must be kept constant between calls.

Parameters:

positions  : Particle positions as a numpy array of shape (n, 3)

cellmatrix : 3x3 numpy array specifying the 3 vectors which define  
             the simulation cell. Implementation for non-orthogonal 
             cells is not present so this must be diagonal.

") compute_model_energy;

%feature("docstring","

Computes the energy of the ith particle due to interactions with its
entries in the Verlet neighbour list, which must be up to date.
The number of particles n is inferred from the size of the positions
array. This must be kept constant between calls.

Parameters:

i          : Index of particle i

positions  : Particle positions as a numpy array of shape (n, 3)

cellmatrix : 3x3 numpy array specifying the 3 vectors which define  
             the simulation cell. Implementation for non-orthogonal 
             cells is not present so this must be diagonal.

") compute_local_real_energy;



%feature("docstring", "

Updates the Verlet neighbour list storing information on neighbours of each particle. 
Particles within a radius of 3.5 + skin are added to the list. The list is initially
populated on first call to compute_model_energy or compute_local_real_energy.

After initialisation, the Python programmer is responsible for ensuring that the neighbour
list is updated with sufficient frequency such that all two-body interactions within the 
cut-off are maintained in the list at all times.

Parameters:

positions  : Particle positions as a numpy array of shape (n, 3)

cellmatrix : 3x3 numpy array specifying the 3 vectors which define  
             the simulation cell. Implementation for non-orthogonal 
             cells is not present so this must be diagonal.

") compute_neighbour_list; 

/* Numpy array typemap */
%include "numpy.i"

%init %{
  import_array();
%}

/* Map array array onto arguments used in the C interface */

/* Atom positions */
%apply (int DIM1, int DIM2, double* IN_ARRAY2) {(int positions, int d, double* pos)};

/* Matrix of cell vectors */
%apply (int DIM1, int DIM2, double* IN_ARRAY2) {(int cellmatrix, int dh1, double* hmatrix)};

double compute_lattice_energy(int positions, int d, double *pos, int cellmatrix, int dh1, double *hmatrix, int ilat);
double compute_local_energy(int i, int positions, int d, double *pos, int cellmatrix, int dh1, double *hmatrix, int ilat);      
void compute_neighbour_list(int positions, int d, double *pos, int cellmatrix, int dh1, double *hmatrix, int ilat);
