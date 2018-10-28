#!/bin/bash

trap "echo; echo 'Something went wrong. Moving on to next directory ...'; cd ..; continue;" ERR

HOMEDIR=$(pwd)

# Set which directories to enter and fetch data from
declare -a LEVEL1=("" "")
declare -a LEVEL2=("" "")
declare -a LEVEL3=("" "")

HOMEDIR=$(pwd)

# Keep track of the number of lines added to files
Nl=0

# Warn the use if file already exists (data will be appended to previous data)
if [ -f $2 ]
then
    echo "WARNING: file $2 already exists! Data will be appended to file."
    echo
fi

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
                python convergence_wl/convergence_wl.py report >> $HOMEDIR/$2
                Nl=$((Nl+1))
                echo "done."

            elif [ "$1" == dF ]
            then
                echo -n "running convergence_dF.py ..."
                python convergence_dF.py report >> $HOMEDIR/$2
                Nl=$((Nl+1))
                
                echo -n "adding to series_$2 ..."
                cat deltaF_mean_series.txt >> $HOMEDIR/series_$2
                echo "-1 -1" >> $HOMEDIR/series_$2
                echo "done."
            
            elif [ "$1" == rt ]
            then
                echo -n "grabbing round trip time ..."
                cp $HOMEDIR/round_trips.py .
                python round_trips.py >> $HOMEDIR/$2
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
echo "$Nl lines added to $2"

