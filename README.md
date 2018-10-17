# py-lsmc

*py-lsmc* stands for *"Python - Lattice Switch Monte Carlo"*.


The purpose of the Lattice-switch algorithm<sup>1</sup> is to enable efficient computation of relative thermodynamic properties of different crystalline forms (often referred to as polymorphs) of a system.

*py-lsmc* combines highly optimised low-level FORTRAN routines with a convenient Python front-end to arrive at accurate Free energy differences between polymorphs of model systems.
Parallelisation via domain decomposition allows systems containing thousands of atoms to be tackled using additional CPUs.

## Getting started

### Prerequisites

* [NumPy, SciPy and Matplotlib](https://www.scipy.org/install.html)
* [Atomic Simulation Environment (ASE)](https://wiki.fysik.dtu.dk/ase/install.html)
* [pymbar](https://github.com/choderalab/pymbar)
* [SWIG](http://www.swig.org/download.html)
* [GNU Parallel](https://www.gnu.org/software/parallel/)


### Installing

Clone the repository using
```
git clone https://github.com/marshrossney/py-lsmc.git
```
or just download and unzip.

### Documentation

For instructions on how to set up and run simulations, please refer to the [wiki](https://github.com/marshrossney/py-lsmc/wiki).

The files themselves also contain a brief explanation of their purpose and, if applicable, instructions for use.

## Overview of usage

Basic usage can be generally broken down into the following steps:

1. Setting up a simulation.
    * Choosing an interaction potential.
    * Specifying the two lattices.
    * Specifying the domain: limits, discretisation, subdomain boundaries, overlaps.


2. Building a set of weights.
    * Employing the Wang-Landau algorithm.
    * Parallelisation via domain decomposition.


3. Calculating the free energy difference.
    * Employing either the Multicanonical algorithm, or
    * The Transition Matrix algorithm.
    * Parallelisation via domain decomposition or multiple global simulations.

Additional functionality includes:

* Measurements of the simulation dynamics and efficiency on a per-bin basis: diffusivity, move acceptance rates...
* Measurements of dynamics and efficiency via series analysis: decorrelation time, round-trip rate...
* Efficiency comparisons between simulations employing different strategies.

## Credit

*py-lsmc* was created by myself (Joe Marsh Rossney) during my final-year undergraduate research project at Warwick University.
It uses energy calculation routines written by Dr David Quigley.
With massive thanks to Hayden Miller and David Quigley for many brilliant discussions which were instrumental in the making of *py-lsmc*.


## References
<sup>1</sup> A. D. Bruce, N. B. Wilding and G. J. Ackland, *Phys. Rev. Lett.* **79**, 3002 (1997)
