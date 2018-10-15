#!/bin/bash

# Should have the same path to working directory on both ORAC and here.
export PTH=$(pwd)

echo "Copying output files from orac..."

# insert your username here.
scp USERNAME@orac.csc.warwick.ac.uk:"$PTH/params.py $PTH/all*.out $PTH/Cmat*.out" .
