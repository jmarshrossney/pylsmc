#!/bin/bash

# Should have the same path to working directory both on ORAC and here.
export PTH=$(pwd)

echo "Sending new input files to orac..."

# Insert your username here.
scp $PTH/hist_fixedw_p* USERNAME@orac.csc.warwick.ac.uk:$PTH/


# Setting up a new working directory - copy everything from here to ORAC.
#echo "Sending everything to orac..."
#scp $PTH/* USERNAME@orac.csc.warwick.ac.uk:$PTH/
