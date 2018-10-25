#!/bin/bash

Ns=$(grep "Ns = " params.py | awk '{print $3}')
sMAX=$(($Ns - 1))
sLIST=$(seq 0 $sMAX)

declare -a Fminus1=("05" "02" "01" "006" "003" "002" "0008" "0004" "0002" "0001" "00005" "00002" "00001" "000006" "000003" "000001" "0000007" "0000004" "0000002")

# Compute dF for each Wang-Landau iteration
l=1
for F in ${Fminus1[@]}
do
    stem="weights_F1-$F"

    #python join_subdoms.py w $stem

    #wfile=${stem}_scomb.out
    wfile="${stem}_s0.out"

    sweeps=$(sed -n "${l}p" < sweeps_allF_smean.out)

    python calc_dF.py $wfile DELTA_F.out $sweeps
    
    #rm $wfile

    l=$((l+1))

done
    

