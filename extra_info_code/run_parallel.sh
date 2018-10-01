#!/bin/bash
echo

trap "echo ' Aborting...'; exit;" SIGINT

# Number of processors/subdomains
NP=$(grep "Np = " params.py | awk '{print $3}')
PMAX=$(expr $NP - 1)
PLIST=$(seq 0 $PMAX)

# Type of simulation
TYPE=$(grep "algorithm = " params.py | awk '{print $3}')

# Number of iterations (for multicanonical/TM simulations)
NITERS=$(grep "iterations = " params.py | awk '{print $3}')
ITER=0

# Initial checks
echo "-----> check_params.py"; python check_params.py || exit 1


# Wang-Landau
if [ "$TYPE" == "'wang_landau'" ]
then
    echo "-----> main.py"; parallel  python main.py ::: $PLIST || exit 1

    echo "-----> join_subdoms.py w"; python join_subdoms.py w

    echo "Finished."
fi


# Multicanonical - run until converged
if [ "$TYPE" == "'multicanonical'" ] 
then
    CONVERGED="False"
    while [ "$CONVERGED" == "False" ]
    do
        ITER=$((ITER+1))
        echo "Iteration $ITER"; echo
        
        echo "-----> main.py"; parallel python main.py ::: $PLIST || exit 1
        
        echo "-----> join_procs.py h"; python join_procs.py h

        echo "-----> unfold_weights.py"; python unfold_weights.py h

        echo "-----> join_subdoms.py h"; python join_subdoms.py h
              
        echo "-----> calc_dF.py"; python calc_dF.py || exit 1
        
        if (( $ITER >= $NITERS ))
        then
            CONVERGED=$(python convergence_dF.py eval)
        fi

    done

    echo "Finished - $ITER iterations completed."
fi


# Transition matrix - run until converged
if [ "$TYPE" == "'transition'" ] 
then
    CONVERGED="False"
    while [ "$CONVERGED" == "False" ]
    do
        ITER=$((ITER+1))
        echo "Iteration $ITER"; echo
        
        echo "-----> main.py"; parallel python main.py ::: $PLIST || exit 1

        echo "-----> join_procs.py h"; python join_procs.py h
        echo "-----> join_procs.py c"; python join_procs.py c
        
        echo "-----> join_subdoms.py e"; python join_subdoms.py e
        echo "-----> join_subdoms.py w"; python join_subdoms.py w

        echo "-----> calc_dF.py"; python calc_dF.py || exit 1

        if (( $ITER >= $NITERS ))
        then
            CONVERGED=$(python convergence_dF.py eval)
        fi

    done

    echo "Finished - $ITER iterations completed."

fi

