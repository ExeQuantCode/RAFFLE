#!/bin/bash 

gaussify () {

cap=$1
bins=$2
rm "$3"_"$4"_bond_gauss
for i in $(seq 1 $bins);
do
    echo $(echo "$i*$cap/$bins" | bc -l) "0" >> "$3"_"$4"_bond_gauss
done 
sigma=0.1

} 

gaussify $1 $2 $3 $4
