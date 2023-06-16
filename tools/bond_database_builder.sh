#!/bin/bash 


rm bond_element_tempfile
pairing_string=$(echo $(ls *_evolved_bondlength) | sed -e 's|_evolved_bondlength||g') 
echo $pairing_string
loop_control=1
while [[ $loop_control -eq 1 ]]
do 
    first_element=$(echo $pairing_string | sed -n 's| .*||p')

    if [[ -z $first_element ]]; then 
	first_element=$pairing_string 
	loop_control=0
    fi
    trimmed_string=$(echo $pairing_string | sed -n 's|[a-zA-Z]*_[a-zA-Z]* ||p')
    
    pairing_string=$trimmed_string
    echo $first_element >> bond_element_tempfile    
done

rm angle_element_tempfile
single_string=$(echo $(ls *_evolved_angles) | sed -e 's|_evolved_angles||g')
echo $single_string
loop_control=1
while [[ $loop_control -eq 1 ]]
do
    first_element=$(echo $single_string | sed -n 's| .*||p')

    if [[ -z $first_element ]]; then
        first_element=$single_string
        loop_control=0
    fi
    trimmed_string=$(echo $single_string | sed -n 's|[a-zA-Z]* ||p')

    single_string=$trimmed_string
    echo $first_element >> angle_element_tempfile
done

