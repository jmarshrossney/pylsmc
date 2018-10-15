#!/bin/bash

declare -a DIRS=("2fcc" "B167" "B263" "B303" "B333")

for i in {0..11}
do
    if [ "$1" == "$i" ]
    then
        cd 2fcc
        p=$1
        #echo "$p in dir $(pwd)"
        python main.py $p
        cd ..
    fi
done

for i in {12..23}
do
    if [ "$1" == "$i" ]
    then
        cd B167
        p=$(($1 - 12))
        #echo "$p in dir $(pwd)"
        python main.py $p
        cd ..
    fi
done

for i in {24..35}
do
    if [ "$1" == "$i" ]
    then
        cd B263
        p=$(($1 - 24))
        #echo "$p in dir $(pwd)"
        python main.py $p
        cd ..
    fi
done

for i in {36..47}
do
    if [ "$1" == "$i" ]
    then
        cd B303
        p=$(( $1 - 36))
        #echo "$p in dir $(pwd)"
        python main.py $p
        cd ..
    fi
done

for i in {48..59}
do
    if [ "$1" == "$i" ]
    then
        cd B333
        p=$(($1 - 48))
        #echo "$p in dir $(pwd)"
        python main.py $p
        cd ..
    fi
done
