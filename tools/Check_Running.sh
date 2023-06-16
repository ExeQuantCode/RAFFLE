#!/bin/bash 

for i in {001..070}
do 
    updated=$(echo $(find -cmin -600) | sed -n "s/.*POSCAR_.*${i}.*/$i/p")
    if [[ ${updated} == "" ]]; then
	echo $i "has not been updated recently. Checking"
	location="/gpshome/jp552/DPhD/DCoding/DLearning/DCodeStore/pos/POSCAR_"$i"/RELAX"
	if [[ -z $(grep "stopping structural" $location/out.out) ]]; then 
	    testvariable=$(echo $location "i")
	    echo $testvariable
	    
	    if [[ -z $(grep "$testvariable" /scratch/Jobq.txt) ]]; then 
		
		less /scratch/Jobq.txt
		less $location/e.out
		less $location/out.out
		ls -l -t $location
		
		echo "rerun, Y/N"
		read continue 
		if [[ $continue == "Y" ]]; then
		    echo "manual restart" >> $location/e.out
		fi
		echo "Next? Y/N "
		read next 
		if [[ $next == "N" ]]; then
		    exit 
		fi
	    else
		echo "Was running on SC. Will skip"
	    fi
	else 
	    echo "Finished relaxing, move along soldier" 
	fi 
	
    fi
    updated=""
done
