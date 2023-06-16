#!/bin/bash 


best_energy () {

while read line; do 
    if [[ -z $line ]]; then 
	break 
    fi
    if [[ -z $(echo $line | sed -n "s|[0-9]*.[0-9]* ||p") ]]; then 
	continue 
    fi
    #break line down into location delimiter and energy value 
    structure_location=$(echo $line | sed -n 's|\..*||p' )
    structure_stage=$(echo $(echo $line | sed -n 's|[0-9]*\.||p' )| sed -n 's| \-*[0-9]*\.*[0-9]*||p')
    
    echo $structure_location $structure_stage
    if [[ $structure_stage == "1" ]]; then

        reader=$(sed -n '7p' < pos/POSCAR_$structure_location/POSCAR)
        element_names=$(sed -n '6p ' < pos/POSCAR_$structure_location/POSCAR)
        #cell_volume_a=$(sed -n '3p' < pos/POSCAR_$structure_location/POSCAR)
        #cell_volume_b=$(sed -n '4p' < pos/POSCAR_$structure_location/POSCAR)
    elif [[ $structure_stage == "2" ]]; then
        reader=$(sed -n '7p' < pos/POSCAR_$structure_location/RELAX/CONTCAR)
        element_names=$(sed -n '6p ' < pos/POSCAR_$structure_location/RELAX/CONTCAR)
        #cell_volume_a=$(sed -n '3p' <pos/POSCAR_$structure_location/RELAX/CONTCAR)
        #cell_volume_b=$(sed -n '4p' <pos/POSCAR_$structure_location/RELAX/CONTCAR)
    else
        reader=$(sed -n '7p' < pos/POSCAR_$structure_location/RELAX/RELAX2/CONTCAR)
        element_names=$(sed -n '6p ' < pos/POSCAR_$structure_location/RELAX/RELAX2/CONTCAR)
        #cell_volume_a=$(sed -n '3p' <pos/POSCAR_$structure_location/RELAX/RELAX2/CONTCAR)
        #cell_volume_b=$(sed -n '4p' <pos/POSCAR_$structure_location/RELAX/RELAX2/CONTCAR)
    fi
    el1=$(echo $reader | cut -d' ' -f1)
    tester=$(echo $reader | sed -n "s/ //p")

    echo $el1
    atomsp=0
    

    if [[ $(echo $reader | sed -n 's|[0-9]* ||p') == "" ]]; then 
	atoms=$reader
    else 
    while [[ $reader != "" ]]; 
    do 
	trim_string=$(echo $reader | sed -n 's|[0-9]* ||p')
	el_count=$(echo $reader | sed -n 's| .*||p')
	reader=$trim_string


	if [[ $reader == "" ]]; then
            break
        fi
	atoms=$(echo "$atomsp+$el_count" | bc -l ) 
	atomsp=$atoms
    done
    fi
    #interface_area="1"
    #a_1=$(echo $cell_volume_a | sed -n 's| .*||p')
    #a_2=$(echo $(echo $cell_volume_a | sed -n "s|\-*[0-9]*\.*[0-9]* ||p") | sed -n 's| .*||p')
    #a_3=$(echo $cell_volume_a | sed -e 's|.* ||g')
    #b_1=$(echo $cell_volume_b | sed -n 's| .*||p')
    #b_2=$(echo $(echo $cell_volume_b | sed -n "s|\-*[0-9]*\.*[0-9]* ||p") | sed -n 's| .*||p')
    #b_3=$(echo $cell_volume_b | sed -e 's|.* ||g')
    
    #sq_1=$(echo "(($a_2*$b_3)-($a_3*$b_2))^2" | bc -l)
    #sq_2=$(echo "(($a_1*$b_3)-($a_3*$b_1))^2" | bc -l)
    #sq_3=$(echo "(($a_1*$b_2)-($a_2*$b_1))^2" | bc -l)
    #interface_area=$(echo "sqrt($sq_1+$sq_2+$sq_3)" | bc -l)
    #if echo "$interface_area < 0" | bc -l | grep -q 1
    #then
    #    interface_area=$(echo "- $interface_area" | bc -l)
    #fi
    readin_location=$(echo $line | sed -n 's| .*||p')
    readin_energy=$(echo $line | sed -n 's|.* ||p')
    if [[ -z ${readin_energy} ]]; then 
	continue 
    fi
    #if [[ -z $interface_area ]]; then 
	#continue 
    #fi
    
    echo "---------------------------------------------------------------------------------------"

    echo $line
    if [[ -z ${best_energy+x} ]]; then
	echo "$readin_energy/$atoms"
	best_energy=$(echo "$readin_energy/$atoms" | bc -l)
	best_location=$readin_location
    elif [[ -z $readin_energy ]]; then 
	:
    elif ((  $(bc -l <<< "($readin_energy/$atoms) < $best_energy" ) )); then 
	best_energy=$(echo "($readin_energy)/$atoms" | bc -l)
	
	#echo $best_energy $readin_energy $interface_area
	best_location=$readin_location
    fi
    echo $(echo "$readin_energy/$atoms" | bc -l)
    echo $best_energy "THIS IS THE BEXT ENERGY CHECK FR BROKEN"

    echo "---------------------------------------------------------------------------------------"
	

done < tmp_energies.txt
} 



create_weightings () {

rm weightings.txt
touch weightings.txt
while read line; do 
    if [[ -z $line ]]; then 
	break
    fi
    echo $line 
    structure_location=$(echo $line | sed -n 's|\..*||p' )
    structure_stage=$(echo $(echo $line | sed -n 's|[0-9]*\.||p' )| sed -n 's| \-*[0-9]*\.*[0-9]*||p')
    echo $structure_location $line "is the stage that this calculation should be compared to"
    if [[ $structure_stage == "1" ]]; then

        reader=$(sed -n '7p' < pos/POSCAR_$structure_location/POSCAR)
        element_names=$(sed -n '6p ' < pos/POSCAR_$structure_location/POSCAR)
        #cell_volume_a=$(sed -n '3p' < pos/POSCAR_$structure_location/POSCAR)
        #cell_volume_b=$(sed -n '4p' < pos/POSCAR_$structure_location/POSCAR)
    elif [[ $structure_stage == "2" ]]; then
        reader=$(sed -n '7p' < pos/POSCAR_$structure_location/RELAX/CONTCAR)
        element_names=$(sed -n '6p ' < pos/POSCAR_$structure_location/RELAX/CONTCAR)
        #cell_volume_a=$(sed -n '3p' <pos/POSCAR_$structure_location/RELAX/CONTCAR)
        #cell_volume_b=$(sed -n '4p' <pos/POSCAR_$structure_location/RELAX/CONTCAR)
    else
        reader=$(sed -n '7p' < pos/POSCAR_$structure_location/RELAX/RELAX2/CONTCAR)
        element_names=$(sed -n '6p ' < pos/POSCAR_$structure_location/RELAX/RELAX2/CONTCAR)
        #cell_volume_a=$(sed -n '3p' <pos/POSCAR_$structure_location/RELAX/RELAX2/CONTCAR)
        #cell_volume_b=$(sed -n '4p' <pos/POSCAR_$structure_location/RELAX/RELAX2/CONTCAR)
    fi
    el1=$(echo $reader | cut -d' ' -f1)
    tester=$(echo $reader | sed -n "s/ //p")

    #interface_area="1"
    #a_1=$(echo $cell_volume_a | sed -n 's| .*||p')
    #a_2=$(echo $(echo $cell_volume_a | sed -n "s|\-*[0-9]*\.*[0-9]* ||p") | sed -n 's| .*||p')
    #a_3=$(echo $cell_volume_a | sed -e 's|.* ||g')
    #b_1=$(echo $cell_volume_b | sed -n 's| .*||p')
    #b_2=$(echo $(echo $cell_volume_b | sed -n "s|\-*[0-9]*\.*[0-9]* ||p") | sed -n 's| .*||p')
    #b_3=$(echo $cell_volume_b | sed -e 's|.* ||g')

    #sq_1=$(echo "(($a_2*$b_3)-($a_3*$b_2))^2" | bc -l)

    #sq_2=$(echo "(($a_1*$b_3)-($a_3*$b_1))^2" | bc -l)
    #sq_3=$(echo "(($a_1*$b_2)-($a_2*$b_1))^2" | bc -l)
    #interface_area=$(echo "sqrt($sq_1+$sq_2+$sq_3)" | bc -l)

    atomsp=0
    while [[ $reader != "" ]];
    do
        trim_string=$(echo $reader | sed -n 's|[0-9]* ||p')
        el_count=$(echo $reader | sed -n 's| .*||p')
        reader=$trim_string
        if [[ $reader == "" ]]; then
            break
        fi
        atoms=$(echo "$atomsp+$el_count" | bc -l )
        atomsp=$atoms
    done



    #if echo "$interface_area < 0" | bc -l | grep -q 1
    #then
    #    interface_area=$(echo "- $interface_area" | bc -l)
    #fi

    readin_location=$(echo $line | sed -n 's| .*||p')
    readin_energy=$(echo $line | sed -n 's|.* ||p')
    if [[  -z $readin_energy ]]; then 
	continue 
    fi
    #if [[ -z $interface_area ]]; then 
#	echo "here"
#	continue 
 #   fi
  #  echo $interface_area
    echo $best_energy

    # INTERFACE AREA IS BAD, WHAT IF RELAX?
    executable="e(-($readin_energy/$atoms)+$best_energy)"
    
    weighting_output=$(bc -l <<< "$executable") 

     echo $readin_location $weighting_output
     echo "()*)(*)(*)(*)(*)(*)(*)(*)(*)(*)(*)(*)(*)(*)(*)(*)(*)(*)(*"
    echo $readin_location $weighting_output>> weightings.txt
    
done < tmp_energies.txt

} 






create_distribution () {
#element_string_parent=$(echo $(grep "elements=" Infile.txt) | sed -n "s|elements=||p")
#stochio_list_parent=$(echo $(grep "stochiometry=" Infile.txt) | sed -n "s|stochiometry=||p")
#eltot=$(echo $(grep "speciestotal=" Infile.txt) | sed -n "s|speciestotal=||p")

rm *_evolved_angle*
rm *_evolved_bond*
rm *_evolved_4*


while read weightings_file; do
    echo "----------------------------------------------------------------------------------------------------------"
    if [[ -z $weightings_file ]]; then
        break
    fi
    structure_location=$(echo $weightings_file | sed -n 's|\..*||p' )
    structure_stage=$(echo $(echo $weightings_file | sed -n 's|[0-9]*\.||p' )| sed -n 's| \-*[0-9]*\.*[0-9]*||p') 
    echo $structure_location $structure_stage "is the stage that this calculation should be compared to"
    if [[ $structure_stage == "1" ]]; then 
	
	reader=$(sed -n '7p' < pos/POSCAR_$structure_location/POSCAR)
	element_names=$(sed -n '6p ' < pos/POSCAR_$structure_location/POSCAR)
	#cell_volume_a=$(sed -n '3p' < pos/POSCAR_$structure_location/POSCAR)
        #cell_volume_b=$(sed -n '4p' < pos/POSCAR_$structure_location/POSCAR)
    elif [[ $structure_stage == "2" ]]; then 
	reader=$(sed -n '7p' < pos/POSCAR_$structure_location/RELAX/CONTCAR)
	element_names=$(sed -n '6p ' < pos/POSCAR_$structure_location/RELAX/CONTCAR)
	#cell_volume_a=$(sed -n '3p' <pos/POSCAR_$structure_location/RELAX/CONTCAR)
        #cell_volume_b=$(sed -n '4p' <pos/POSCAR_$structure_location/RELAX/CONTCAR)
    else
	reader=$(sed -n '7p' < pos/POSCAR_$structure_location/RELAX/RELAX2/CONTCAR)
	element_names=$(sed -n '6p ' < pos/POSCAR_$structure_location/RELAX/RELAX2/CONTCAR)
	#cell_volume_a=$(sed -n '3p' <pos/POSCAR_$structure_location/RELAX/RELAX2/CONTCAR)
        #cell_volume_b=$(sed -n '4p' <pos/POSCAR_$structure_location/RELAX/RELAX2/CONTCAR)
    fi
    el1=$(echo $reader | cut -d' ' -f1)
    tester=$(echo $reader | sed -n "s/ //p")
    
    #interface_area="1"
    #echo $cell_volume_a 
    #echo $cell_volume_b 
    
    #a_1=$(echo $cell_volume_a | sed -n 's| .*||p')
    #a_2=$(echo $(echo $cell_volume_a | sed -n "s|\-*[0-9]*\.*[0-9]* ||p") | sed -n 's| .*||p')
    #a_3=$(echo $cell_volume_a | sed -e 's|.* ||g')
    #b_1=$(echo $cell_volume_b | sed -n 's| .*||p')
    #b_2=$(echo $(echo $cell_volume_b | sed -n "s|\-*[0-9]*\.*[0-9]* ||p") | sed -n 's| .*||p')
    #b_3=$(echo $cell_volume_b | sed -e 's|.* ||g')

    #sq_1=$(echo "(($a_2*$b_3)-($a_3*$b_2))^2" | bc -l)

    #sq_2=$(echo "(($a_1*$b_3)-($a_3*$b_1))^2" | bc -l)
    #sq_3=$(echo "(($a_1*$b_2)-($a_2*$b_1))^2" | bc -l)
    #interface_area=$(echo "sqrt($sq_1+$sq_2+$sq_3)" | bc -l)

    #if echo "$interface_area < 0" | bc -l | grep -q 1
    #then
#	interface_area=$(echo "- $interface_area" | bc -l)
 #   fi

  #  echo $interface_area
    
    
    
    status=0
    if [[ -z $tester ]]; then
        attotp=$el1
    else
        while [[ $status -ne -1 ]]
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
    
    eltot=$elle
    stochio_string=$reader
    element_string=$element_names
    echo "Beginning loop. Initial strings:" $element_string $stochio_string
    
    for i in $(seq 1 $eltot) 
    do 
	echo "Loop step" $i "initialised. Strings:" $element_string $stochio_string
	    
	stochio_current=$(echo $stochio_string | sed -n "s| .*||p")

	
	stochio_trimmed=$(echo $stochio_string | sed -n 's/[0-9]* //p')
	element_current=$(echo $element_string | sed -n "s| .*||p")
	element_trimmed=$(echo $element_string | sed -n 's/[a-zA-Z]* //p')
	echo "Current atom/stochio obtained.  Trimmed strings are now" $element_trimmed $stochio_trimmed 
	if [[ -z $stochio_current ]]; then 
	    stochio_current=$(echo $stochio_string | sed -e 's/ //g')
	fi
	if [[ -z $element_current ]]; then 
	    echo $element_current":"$element_string":"
	    element_current=$(echo $element_string | sed -e 's/ //g')
	fi
	
	element_string=$element_trimmed 
	stochio_string=$stochio_trimmed 
	
	echo "WELCOME. Stochiometries being used are" $element_current ":" $stochio_current
	echo "////////////////////////////////////////////////////////////////"
	echo $weightings_file
        echo "////////////////////////////////////////////////////////////////"

	location=$(echo $weightings_file | sed -n "s/ .*//p")
	weighting=$(echo $weightings_file | sed -n "s/.* //p") 
	echo $weighting
	#weighting_division=$(echo "$weighting / $interface_area" | bc -l)
	#weighting=$weighting_division
	stage=$(echo $location | sed -n "s/.*\.//p") 
	if [[ $stage == "1" ]]; then 
	    bondir="bon"
	    baddir="bad"
	    fourboddir="4body"
	else 
	    bondir="bon"$stage
	    baddir="bad"$stage
	    fourboddir="4body"$stage
	fi
	echo "////////////////////////////////////////////////////////////////"
	echo $location "///" $(echo $location | sed -n "s/\..*//p") "---" 
	
	echo "////////////////////////////////////////////////////////////////"

	identity=$(echo $location | sed -n "s/\..*//p")
	full_path_bon=$bondir"/BON_"$identity
	full_path_bad=$baddir"/BAD_"$identity
        full_path_fourbody=$fourboddir"/4BOD_"$identity
	echo $full_path_bon $full_path_bad $stochio_current $identity "££££££££££££££££££££££££££"
	for j in $(seq 1 $stochio_current)
	do 
	    if (( ( "$j<10" | bc -l) )); then
	        file_number=$(echo 00$j)
	    elif (( ( "$j<100" | bc -l) )); then
		file_number=$(echo 0$j)
	    else
		file_number=$(echo $j)
	    fi
	    atom_path_four=$full_path_fourbody"/"$element_current$"_"$file_number  
	    atom_path_bon=$full_path_bon"/"$element_current$"_"$file_number
	    atom_path_bad=$full_path_bad"/"$element_current$"_"$file_number
	    echo $atom_path_bad $full_path_bad $element_current $file_number
	   
	    #echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
	    while read atom_details; do 
		if [[ -z $atom_details ]]; then 
		    break 
		else 
		    #echo $atom_details
		    if [[ $(echo $atom_details | sed -n "s/[a-zA-Z]*/NewEl/p") == "NewEl" ]]; then 
			target_element=$atom_details 
		    else
			echo $atom_details $weighting >> $element_current"_"$target_element"_evolved_bondlength"
			if [[ $atom_details == "NaN" ]]; then
			    echo $atom_path_bon
			    exit 
			fi
			echo $atom_details >> $atom_path_bon$element_current"_"$target_element"_bond_distribution"
		    fi 
		fi
	    done < "$atom_path_bon"
	    while read atom_angles; do 
		if [[ -z $atom_angles ]]; then 
		    break 
		else
		    #echo $atom_angles 
		    #echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
		    echo $atom_angles $weighting >> $element_current"_evolved_angles" 
		fi
	    done < "$atom_path_bad"
	    while read atom_four; do 
		if [[ -z $atom_four ]]; then 
		    break 
		else
		    echo $atom_four
		    echo $atom_four $weighting >> $element_current"_evolved_4body"
		fi
	    done < $atom_path_four
	done
    done 
    
done < weightings.txt

} 
best_energy 
create_weightings
create_distribution
