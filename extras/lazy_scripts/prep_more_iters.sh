#!/bin/bash

# Send new input files (e.g. histogram) to ORAC
./send_orac.sh

# Move all the data from this set of iterations into a folder $1
mv Cmat* $1
mv all* $1
cp DELTA_F.out $1
mv hist* $1
