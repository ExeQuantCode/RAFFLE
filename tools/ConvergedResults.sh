#!/bin/bash


rm results0.txt
rm results1.txt
rm results2.txt


################################################################
for i in $(seq -w 001 001)
do 
    echo "_______________________________________________________________________________________________________________________"
    a=$(grep "stopping structural" pos/POSCAR_$i/RELAX/out.out)
    attotp="0"
    status=0

    DIR=pos/POSCAR_$i/POSCAR
    if test -f "$DIR"; then 
	: 
    else
	continue
    fi

    reader=$(sed -n '7p' < pos/POSCAR_$i/POSCAR)
    el1=$(echo $reader | cut -d' ' -f1)
    tester=$(echo $reader | sed -n "s/ //p")
 
    element_names=$(sed -n '6p ' < pos/POSCAR_$i/POSCAR)
    
    if [[ -z $tester ]]; then
        attotp=$el1
    else
        while [[  $status -ne -1 ]]
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

    energy_total=0
    energy=0
    atoms_total_expression=0
    
    
    
    #isolated_host_element_stochiometry=$(echo $(sed -n '7p ' pos/POSCAR_$i/POSCAR) | sed -n "s| .*||p")
    #sqrt_host_dble=$(echo "sqrt (($isolated_host_element_stochiometry)/4)" | bc -l)
    #host_size=$(echo $sqrt_host_dble | sed -n "s|\..*||p")
    
    #if [[ -z $host_size ]]; then 
#	host_size=$sqrt_host_dble
#    fi
#    FILE=pos/POSCAR_$i/host_size_override
#    if test -f "$FILE"; then 
   
#	read host_size1 host_size2 <pos/POSCAR_$i/host_size_override
#        host_energy_string=$(grep "e  e" "host/POSCAR_$host_size1""x""$host_size2/OUTCAR")

#    else 
#	host_energy_string=$(grep "e  e" "host/POSCAR_$host_size""x""$host_size/OUTCAR")
	host_size1=$host_size
	host_size2=$host_size
#    fi
#    host_energy_string_trim=$(echo $host_energy_string | sed -n "s|.* = ||p")
#    host_energy=$(echo $host_energy_string_trim | sed -n "s| eV||p")
	#    echo $host_energy $host_size "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
#    #start_point=$(echo "$host_size+1" | bc -l)
    start_point=1
    for p in $( seq 1 $elle )
    do 
	
	x=$(echo $element_names | sed -n 's|[a-zA-Z]* ||1p')
	step_1_stochio=$(echo $reader | sed -n 's|[0-9]* ||1p')

	if [[ -z $x ]]; then 
	    y=$element_names
	    stochio_i=$reader
	else
	    y=$(echo $element_names | sed -n "s| ${x}||p")
	    element_names=$x
	    stochio_i=$(echo $reader | sed -n "s|${step_1_stochio}||p")
	    
	    reader=$step_1_stochio
	    echo $stochio_i $y
	    

	fi
	location="iso/POSCAR_"$y"/OUTCAR"
	iso_energy=$(echo $(echo $(grep "e  e" $location) | sed -n 's|.* = *||p') | sed -n 's| eV||p')
	
	energy_total=$(echo $stochio_i"*"$iso_energy | bc -l)
        energy=$energy$energy_total

        atoms_total_expression=$atoms_total_expression"+ "$stochio_i
       
    done	
    isolated_binding_energies=$(echo $energy | bc -l)
    atoms_total=$(echo $atoms_total_expression | bc -l)
    

  
    ascf=$(grep "e  e" pos/POSCAR_$i/OUTCAR | tail -1)
    bscf=$(echo $ascf | sed -e "s/free energy TOTEN = //g")
    ascf=$(echo $bscf | sed -e "s/eV//g")
    echo $atoms_total $host_energy $ascf
   
    #evaluationscf="$ascf- $isolated_binding_energies - $host_energy"
    evaluationscf="$ascf- $isolated_binding_energies"

    apscf=$(echo $evaluationscf | bc -l)
    #host_squared=$(echo $host_size1"*"$host_size2 | bc -l)
    concentration=$(echo $atoms_total | bc -l)
    bpscf=$(echo $apscf/$concentration | bc -l)
    
     
    if [[ -z $bpscf ]]; then 
	:
    else
	echo "scf" $i $bpscf  $host_size $atoms_total >> results0.txt
    fi
    echo "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"









    a=$(grep "e  e" pos/POSCAR_$i/RELAX/OUTCAR)


    if [[ -z $a ]]; then
	echo "Empty"
	continue
    else 
	echo $i "is the number of the structure"
   fi

    a=$(grep "e  e" pos/POSCAR_$i/RELAX/OUTCAR | tail -1)
    b=$(echo $a | sed -e "s/free energy TOTEN = //g")
    a=$(echo $b | sed -e "s/eV//g")
    #evaluation="$a- $isolated_binding_energies - $host_energy"
    evaluation="$a- $isolated_binding_energies"
    ap=$(echo "$evaluation" | bc -l)
    bp=$(echo $ap/$concentration | bc -l)
   
    echo $bp "is the value of the energy, and undivided is" $ap "and read" $a "and evaluation expr" $evaluation
    echo "relax1" $i $bp $atoms_total >> results1.txt 


    a=$(grep "e  e" pos/POSCAR_$i/RELAX/RELAX2/OUTCAR)

    if [[ -z $a ]]; then    
	continue
    else 
	echo $i
    fi
    
    
    d=$(grep "e  e" pos/POSCAR_$i/RELAX/RELAX2/OUTCAR | tail -1)
    e=$(echo $d | sed -e "s/free energy TOTEN = //g")
    d=$(echo $e | sed -e "s/eV//g")
    #evaluation="$d - $isolated_binding_energies - $host_energy"
    evaluation="$d -$isolated_binding_energies"
    e=$(echo $evaluation | bc -l)
    ep=$(echo $e/$concentration | bc -l)
    
    echo $e $ep "README"
    echo "relax2 "$i $ep $atoms_total>> results2.txt
    
    echo "_______________________________________________________________________________________________________________________"
done
sort -k3 -g results0.txt > results0.tmp
sort -k3 -g results1.txt > results1.tmp 
sort -k3 -g results2.txt > results2.tmp 
nl results1.tmp > results1.txt
rm results1.tmp
nl results2.tmp > results2.txt
rm results2.tmp
nl results0.tmp > results0.txt
rm results0.tmp

xmgrace -nxy results1.txt -nxy results2.txt
