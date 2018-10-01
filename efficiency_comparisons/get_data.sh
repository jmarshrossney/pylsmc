#!/bin/bash

trap "echo; echo 'Something went wrong. Moving on to next directory ...'; cd ..; continue;" ERR

HOMEDIR=$(pwd)

# Set which directories to enter and fetch data from
declare -a LEVEL1=("" "")
declare -a LEVEL2=("" "")
declare -a LEVEL3=("" "")

# Keep track of the number of lines added to files
Nl=0

for i in "${LEVEL1[@]}"
do
    echo "Entering directory $i ... "
    cd $i

    for j in "${LEVEL2[@]}"
    do
        echo "Entering directory $j ..."
        cd $j

        for k in "${LEVEL3[@]}"
        do
            echo -n "Entering directory $k ..."
            cd $k

            # ----------------------------------------- #
            #  Begin code to execute in each directory  #
            # ----------------------------------------- #
            
            if [ "$1" == wl ]
            then
                echo -n "running convergence_wl.py ..."
                python convergence_wl.py report >> $HOMEDIR/wl_data.txt
                Nl=$((Nl+1))
                echo "done."

            elif [ "$1" == dF ]
            then
                echo -n "running convergence_dF.py ..."
                python convergence_dF.py report >> $HOMEDIR/dF_data.txt
                Nl=$((Nl+1))
                
                echo -n "adding to stdev_series.txt ..."
                cat stdev_series.txt >> $HOMEDIR/stdev_series.txt
                echo "-1 -1" >> $HOMEDIR/stdev_series.txt
                echo "done."
            
            else
                echo "doing nothing ..."
            fi
            
            # -------------------------------------- #
            # End code to execute in each directory  #
            # -------------------------------------- #
            
            cd ..
        done
        cd ..
    done
    cd ..
done

echo "Finished"
echo
echo "$Nl lines added to ${1}_data.txt"

