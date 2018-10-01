#!/bin/bash

trap "echo ' Aborting...'; exit;" SIGINT

# Declare arrays for folder hierarchy
declare -a LEVEL1=("71_bins" "81_bins" "91_bins" "101_bins")
declare -a LEVEL2=("interpolated_weights" "discrete_weights")
declare -a LEVEL3=("trial1" "trial2")

# Directory from where this script was executed
HOMEDIR=$(pwd)

# Go into all the directories
for i in "${LEVEL1[@]}"
do
    cd $i

    for j in "${LEVEL2[@]}"
    do
        cd $j
        
        for k in "${LEVEL3[@]}"
        do
            cd $k
        
            # --------------------------------------- #
            #  Stuff we want to do in each directory  #
            # --------------------------------------- #
        
            #make clean
            #make
        
            #python check_params.py || exit 1 

            #echo "Submitting COWscript in dir $(pwd)"
            #sbatch -p spacebat COWscript.sh
            
            
            # ---------------------------------------------- #
            #  End of stuff we want to do in each directory  #
            # ---------------------------------------------- #
            
            cd ..
        done
        cd ..
    done
    cd ..
done
