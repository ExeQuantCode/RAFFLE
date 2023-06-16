#!/bin/bash 
var=$(<prevstructures.txt)
previous=$(echo $var|sed -n 's| .*||p')
last_run=$(echo $var|sed -n 's|.* ||p')
total=$(echo "$previous + $last_run"| bc -l)








rm similarity_combinations.txt
echo $total>>similarity_combinations.txt


for i in $(seq 1 $total)
do 
    reader=$(sed -n '7p' < ../DC_MgO/pos/POSCAR_$(printf "%03d" $i)/POSCAR)
    if [[ -z $reader ]]; then 
	continue 2 
    fi 
     el1=$(echo $reader | cut -d' ' -f1)
     tester=$(echo $reader | sed -n "s/ //p")
     echo $i
     element_names=$(sed -n '6p' < ../DC_MgO/pos/POSCAR_$(printf "%03d" $i)/POSCAR)
     
    status=0
    statusp=0
    stochio_i=0
    el1=0
    attot=0
    attotp=0
    elle=1
    if [[ -z $tester ]]; then
        attotp=$el1
    else
   
	while [[  $status -ne "-1" ]]
	
	do
            if [[ -z $el1 ]]; then
		elle=$(echo "$status-1" | bc -l)
		status=-1
            else
		
		statusp=$(echo "$status + 1" | bc -l)
		status=$statusp
		#echo "found an element"
		attot=$(echo "$el1 + $attotp" | bc -l)
		attotp=$attot
		el1=$(echo $reader | cut -d' ' -f$status)
		
	    fi
	done
    fi
    echo $i>>similarity_combinations.txt
    echo $elle>>similarity_combinations.txt
    atoms_total_expression=0


    for p in $(seq 1 $elle)
    do
	
	x=$(echo $element_names | sed -n 's|[a-zA-Z]* ||1p')
	step_1_stochio=$(echo $reader | sed -n 's|[0-9]* ||1p')
	if [[ -z $x ]]; then
	    y=$element_names
            stochio_i=$reader
	    echo $stochio_i>>similarity_combinations.txt 
	    echo $y>>similarity_combinations.txt
	    
	else
            y=$(echo $element_names | sed -n "s| ${x}||p")
            element_names=$x
            stochio_i=$(echo $reader | sed -n "s| ${step_1_stochio}||p")
	    reader=$step_1_stochio
	    
	    echo $stochio_i >>similarity_combinations.txt
            echo $y>>similarity_combinations.txt


	fi
	atoms_total_expression=$atoms_total_expression"+ "$stochio_i

    done
    atoms_total=$(echo $atoms_total_expression | bc -l)
    echo $atoms_total>>similarity_combinations.txt
done
#Key 
#First entry=total structures
#General format
#N+1=structure number
#N+2=ath stoichio
#N+3=ath el name 

#N+2a+1=total atoms 
