/* Header file for C-compatible variables/functions in energy.F90 */

/* Initialises a new lattice (creates neighbour lists and image vectors) */
void create_lattice(int n, int d, double *pos, int dh2, int dh1, double *hmatrix,int ilat);

/* Calculates the total energy of the system of particles */
double compute_lattice_energy(int n, int d, double *pos, int dh2, int dh1, double *hmatrix,int ilat);

/* Calculates the energy of particle i interacting with its neighbours */
double compute_local_energy(int i, int n, int d, double *pos, int dh2, int dh1, double *hmatrix,int ilat);

/* Recalculates the Verlet neighbour list. Should be called if atoms have moved sufficiently far */
/* that the originally computed list is invalidated.                                             */
void compute_neighbour_list(int n, int d, double *pos, int dh2, int dh1, double *hmatrix,int ilat);

/* May need to expose additional functions here if cell changes (need to recompute image vectors) */
