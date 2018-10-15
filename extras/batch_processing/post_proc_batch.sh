#!/bin/bash

# Takes the raw outputs, which have been saved periodically, and executes the post-processing
# scripts to give a dF time series.

trap "echo ' Aborting...'; exit;" SIGINT

# Type of simulation
TYPE=$(grep "algorithm = " params.py | awk '{print $3}')

MAXITER=$(grep "iterations = " params.py | awk '{print $3}')
MAXITER=$(expr $MAXITER - 1)

NITERS=$(seq 0 $MAXITER)

if [ "$TYPE" == "'multicanonical'"]
then

    for iter in $NITERS
    do
        echo "Iteration $iter"

        python snapshot.py $iter || exit 1
    
        python join_procs.py h

        python unfold_weights.py

        python join_subdoms.py h

        python calc_dF.py
    done
fi

if [ "$TYPE" == "'transition'"]
then

    for iter in $NITERS
    do
        echo "Iteration $iter"
        
        python snapshot.py $iter || exit 1
    
        python join_procs.py c

        python join_subdoms.py p

        python calc_dF.py
    done
fi


