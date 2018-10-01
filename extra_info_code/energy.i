/* -*- mode: C */

/* energy.i */
%module energy
%{
#define SWIG_FILE_WITH_INIT

/* This will be included in the wrapper code */
#include "energy.h"
%}

/* Numpy array typemap */
%include "numpy.i"

%init %{
  import_array();
%}

/* Map array array onto arguments used in the C interface */

/* Atom positions */
%apply (int DIM1, int DIM2, double* IN_ARRAY2) {(int n, int d, double* pos)};

/* Matrix of cell vectors */
%apply (int DIM1, int DIM2, double* IN_ARRAY2) {(int dh2, int dh1, double* hmatrix)};

/* This will be parsed to generate the wrapper */
%include "energy.h"
