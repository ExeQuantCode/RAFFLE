#!/bin/bash 
### This scrpit manages the QUEUEING system.

vasprun () {
    
   echo $var
    
    ##Jobrunning->determines if vasp should run in that location
    jobrunning
    if [[ $dovasp == "FALSE" ]]; then 
	echo "Job already running, going next"
	echo "---------------------------------------------------------"
    elif [[ $dovasp == "TRUE" ]]; then 
	echo "Job initialised"
	echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++" 
	echo $posloc
	tmploc=$(echo $posloc | sed -e 's/ q//g')
	######OPTIONAL INCAR REWRITING. COMMENT OUT IF NO BLANKET CHANGE REQUIRE#####
	#sed -e "s/ALGO = N/ALGO = V/g" <$tmploc"/INCAR" > $tmploc"/INCAR1"
	#mv $tmploc"/INCAR1" $tmploc"/INCAR"
	
	H=$(hostname)
	vartmp=$(echo $H | sed -e 's/.ex.ac.uk//g')
	if [[ -z $vartmp ]]; then
            H=$HOSTNAME
	else
            H=$vartmp
	fi
	
	if [[ $supercomputer == "NO" ]]; then 
	    #if [[ $H == $var ]]; then  
#		cd $working_directory
#		rm e.out
#		echo "I'm doing this"
#	        module load mpi/openmpi-x86_64; 
#		module load vasp-modulefiles/vasp/5.4.4; 
#		mpirun -np 4 vasp_std >out.out 2>e.out & disown 
#		cd $workdir
#		
#		
	    if [[ $var != "phy-cdt41" ]]; then 
		ssh $var "cd $working_directory; pwd;  module load mpi/openmpi-x86_64; \
module load vasp/5.4.4; ulimit -s unlimited; mpirun -np 4 vasp_std >out.out 2>e.out & disown"
 	    else
		ssh $var "cd $working_directory; pwd; module load mpi/openmpi-x86_64; \
module load vasp-modulefiles/vasp/5.4.4; ulimit -s unlimited; mpirun -np 4 vasp_std >out.out 2>e.out & disown"
		echo $var $working_directory >> sublist.txt
	    fi
	elif [[ $supercomputer == "YES" ]]; then 
	    if [[ $is_isca == "ISCA" ]]; then 
		ssh jp552@login.isca.ex.ac.uk "cd $working_directory; cat e.out >>errorhis.out;rm e.out;sbatch job_vasp_isca.in"
		echo $var $working_directory >> sublist.txt
	    elif [[ $is_archer == "ARCHER" ]]; then 
		ssh jp552@login.archer2.ac.uk "cd $working_directory; cat e.out>>errorhis.out;rm e.out;sbatch job_vasp_archer2.in"
		echo "HERE"
		echo $var $working_directory >> sublist.txt
	    fi
	fi
    fi
}

resallowed () {
    is_isca=""
    is_archer=""
    if [[ $allowed == "" ]]; then 
	tmpvar=$(sed -n 'p' </gpshome/jp552/locations.txt)
	#Setting up "allowed" string, which is cut up to find the machine parameters 
	allowed=$tmpvar
    fi
    var=$(echo $allowed | sed -n 's| .*||p')
    is_isca=$(echo $allowed | sed -e 's|isca.*|ISCA|g')
    is_archer=$(echo $allowed | sed -e 's|archer.*|ARCHER|g')
    allowed=$(echo $allowed | sed -n 's|[0-9a-zA-Z\-]* ||p')
    perc_minfree_cpu=$(echo $allowed | sed -n 's| .*||p')
    allowed=$(echo $allowed | sed -n 's|[0-9]* ||p')
    if [[ $is_isca == "ISCA" ]]; then 
	max_jobs=$perc_minfree_cpu
	supercomputer="YES"
    elif [[ $is_archer == "ARCHER" ]]; then 
	max_jobs=$perc_minfree_cpu
	supercomputer="YES"
    else
	supercomputer="NO"
	max_jobs=$(echo $allowed | sed -n 's| .*||p')
	if [[ $max_jobs == "" ]]; then 
	    max_jobs=$allowed 
	fi
	allowed=$(echo $allowed | sed -n 's|[0-9]* ||p' )
    fi
    perc_minfree_ram=$perc_minfree_cpu
    
    if [[ $supercomputer == "YES" ]]; then 
	
	
	echo $max_jobs "is the total number of jobs permitted on" $var ". This machine is a supercomputer?:" $supercomputer
    else
	echo $max_jobs "is the total number of jobs permitted on" $var ", which may not exceed" $perc_minfree_cpu "% free space"
	echo "Is this machine a supercomputer?" $supercomputer 
    fi
}

previous_run_machine () {
    forced_location=""
#first character gives info on current job status, second character on restrictions on the calculation 
    echo "_________________________________"
    echo $posloc
    echo "_________________________________"
    if [[   $(echo $posloc | sed -e 's|.* q||g') == "" ]]; then 
	#job was queuing. 
	job_status="QUEUING"
	forced_location=""
    fi
    if [[   $(echo $posloc | sed -e 's|.* qi||g') == "" ]]; then
        #job was queuing and is only to be placed on ISCA
        job_status="QUEUING"
	forced_location="ISCA"
    fi
    if [[ $(echo $posloc | sed -e 's|.* qa||g') == "" ]]; then 
	#job queueing for archer
	job_status="QUEUING"
	forced_location="ARCHER"
	echo "________________________"
	echo $forced_location
	echo "------------------------"
   fi
    if [[   $(echo $posloc | sed -e 's|.* i||g') == "" ]]; then
	#job was queuing.
	job_status="ISCA"
	forced_location=""
    fi
    if [[   $(echo $posloc | sed -e 's|.* ii||g') == "" ]]; then
        #job was one isca, and should only be run on isca
        job_status="ISCA"
	forced_location="ISCA"
    fi
    if [[   $(echo $posloc | sed -e 's|.* a||g') == "" ]]; then
        #job was queuing.
        job_status="ARCHER"
	#jobs on archer stay on archer if they don't complete; perhaps they should be manually checked!
        forced_location="ARCHER"
    fi
    if [[   $(echo $posloc | sed -e 's|.* aa||g') == "" ]]; then
        #job was queuing.
        job_status="ARCHER"
        forced_location="ARCHER"
    fi


    if [[   $(echo $posloc | sed -e 's|.* r||g') == "" ]]; then
	#job was 
	job_status="LOCAL"
	forced_location=""
    fi
    if [[   $(echo $posloc | sed -e 's|.* ri||g') == "" ]]; then
        #job was queuing.
        job_status="LOCAL"
	forced_location="ISCA"
    fi
    if [[   $(echo $posloc | sed -e 's|.* ra||g') == "" ]]; then
        #job was queuing.
        job_status="LOCAL"
        forced_location="ACRHER"
    fi

    echo "The location of the last job was" $job_status
}

decision_matrix () {

dovasp="FALSE"
if [[ $errors == "YES" ]]; then
    dovasp="TRUE"
else
    dovasp="UNDEFINED"
fi

if [[ $finished == "YES" ]]; then 
    dovasp="FINISHED"
else 
    if [[ $dovasp == "UNDEFINED" ]]; then 
	: 
    else
	dovasp="TRUE"
    fi
    if [[ $limit_reached == "YES" ]]; then
	dovasp="TRUE"
    fi
    
    if [[ $job_status == "QUEUING" ]]; then
	dovasp="TRUE"
	finished="NO"
    fi
    
    
fi 


}



jobrunning () {
    echo "THIS IS THE JOBRUNNING FUNCTION" 
    dovasp="UNDEFINED"    
    finished="NO"
    if [[ $is_isca == "ISCA" ]]; then
	    determine_appropriate_directory ISCA
    elif [[ $is_archer = "ARCHER" ]]; then 
	    determine_appropriate_directory ARCHER
    else 
	    determine_appropriate_directory LOCAL
    fi
    determine_appropriate_directory NO
    echo $working_directory 
    previous_run_machine
    test_dir_errors
    test_dir_SCF

    if [[ $SCF == "NO" ]]; then
	test_dir_converged
	test_dir_limitreached
	test_dir_running
	decision_matrix
    else 
	limit_reached="NO"
	test_scf_done
	decision_matrix 
    fi
    if [[ $dovasp == "FINISHED" ]]; then 
	tmp=$supercomputer
	supercomputer="NO"
	determine_appropriate_directory NO
	copy_files
	supercomputer=$tmp
	ammend_queue DELETE
    fi
    if [[ $dovasp == "UNDEFINED" ]]; then 
	if [[ $limit_reached == "NO" ]]; then 
	    if [[ $errors == "YES" ]]; then 
		echo "This job may have crashed, but it also may not have crashed"
		poscar_number=$(echo $(echo $posloc | sed -n 's|.*/POSCAR_||p') | sed -n 's|/RELAX.*||p')
		echo $poscar_number 
		#check_if_running $poscar_number       
		ammend_queue CHANGE
	    else 
		ammend_queue MAINTAIN 
	    fi
	fi
    fi
    if [[ $dovasp == "TRUE" ]]; then
        ammend_queue CHANGE
	if [[ $hold_copy != "YES" ]]; then 
	copy_files
	fi
    fi
    if [[ $dovasp == "FALSE" ]]; then 
	ammend_queue MAINTAIN
    fi
    if [[ -z $posloc ]]; then 
	echo $posloc "IS THE POSCAR LOCATION"
	dovasp="FALSE"
    fi
    if [[ $dovasp == "UNDEFINED" ]]; then 
	dovasp="FALSE"
    fi
    echo $finished "is the job finished"
    echo $errors "were there any errors?"
    echo $limit_reached "was the relaxation cap reached"
    echo $dovasp "is the dovasp variable"
    echo $supercomputer "is the machine the desired job is to be run on a supercomputer"
    echo $job_status "is the last known status of the job"
}

ammend_queue () {
    mode=$1
    suffix=""
    if [[ $mode == "DELETE" ]]; then 
	sed -e "s|${posloc}||g" </scratch/Jobq.txt >>/scratch/Jobqp.txt
	mv /scratch/Jobqp.txt /scratch/Jobq.txt
    elif [[ $mode == "CHANGE" ]]; then 
	hold_copy="NO"
	if [[ $job_status == "ISCA" ]]; then  
	    if [[ $supercomputer == "YES" ]]; then  
		if [[ $forced_location == "ARCHER" ]]; then 
		    echo "WHAT DO YOU MEAN YOU WANT TO MOVE A JOB FROM ISCA TO ARCHER? ARE YOU MAD?"
		    exit
		elif [[ $forced_location == "ISCA" ]]; then  
		    prefix=$posloc
		    suffix=""
		elif [[ $is_archer == "ARCHER" ]]; then 
		    #THIS ROUTINE WOULD NEED TO COPY FROM ISCA TO ARCHER, BUT SEEING AS WE DONT NEED THIS, I WONT DO IT
		    echo "This is an absolute wrong.THIS ROUTINE WOULD NEED TO COPY FROM ISCA TO ARCHER, BUT SEEING AS WE DONT NEED THIS, I WONT DO IT" 
		    exit
		elif [[ $is_isca == "ISCA" ]]; then 
		    prefix=$posloc
		    suffix=""
		fi
	    elif [[ $supercomputer == "NO" ]]; then 
		if [[ $forced_location == "ISCA" ]]; then 
		    prefix=$(echo $posloc | sed -n 's| .*||p')
                    suffix=" ii"
		    hold_copy="YES"
		elif [[ $forced_location == "ARCHER" ]]; then 
		    echo "Jobs cant move from isca to archer. Fix your shit"
		    exit
		else
		    prefix=$(echo $posloc | sed -n 's| .*||p')
		    suffix=" r"
		fi
	    fi
	elif [[ $job_status == "ARCHER" ]]; then 
	    if [[ $supercomputer == "YES" ]]; then
                if [[ $forced_location == "ARCHER" ]]; then 
		    prefix=$(echo $posloc | sed -n 's| .*||p')
                    suffix=" aa"
                    if [[ $is_archer != "ARCHER" ]]; then 
			hold_copy="YES"
		    fi
		elif [[ $forced_location == "ISCA" ]]; then 
		    echo "this job thinks it wants to move from archer to isca. It shouldn't" 
		    exit
		    ####No job should ever have been on ARCHER and then moved to ISCA. Exiting, no cycle possible with force
		    #prefix=$(echo $posloc | sed -n 's| .*||p')
                    #suffix=" aa"
		    #if [[ $is_isca != "ISCA" ]]; then
                    #    hold_copy="YES"
                    #fi

                    #hold_copy="YES"
		elif [[ $is_archer == "ARCHER" ]]; then 
		    prefix=$(echo $posloc | sed -n 's| .*||p')
                    suffix=" aa"
		elif [[ $is_isca == "ISCA" ]]; then 
		    ### No job should ever have been on archer, and then be running on ISCA. Cycling machine 
		    prefix=$(echo $posloc | sed -n 's| .*||p')
                    suffix=" aa"
                    hold_copy="YES"
                fi
	    elif [[ $supercomputer == "NO" ]]; then
                if [[ $forced_location == "ARCHER" ]]; then
                    prefix=$(echo $posloc | sed -n 's| .*||p')
                    suffix=" aa"
                    hold_copy="YES"
                fi
		#There is no allowment for jobs to come from archer to local machines automatically
	    else 
		#It is probably that all jobs on archer should always be run only on archer. This forces all jobs to stay on archer 
                forced_location == "ARCHER"
		prefix=$(echo $posloc | sed -n 's| .*||p')
                suffix=" aa"
                hold_copy="YES"
            fi
	elif [[ $job_status == "LOCAL" ]]; then 
	    if [[ $supercomputer == "NO" ]]; then 
		prefix=$posloc
		suffix=""
	    elif [[ $supercomputer == "YES" ]]; then 
		if [[ $forced_location == "ISCA" ]]; then 
		    prefix=$(echo $posloc | sed -n 's| .*||p')
		    suffix=" ii"
		    hold_copy="YES"
		elif [[ $forced_location == "ARCHER" ]]; then 
		    prefix=$(echo $posloc | sed -n 's| .*||p')
                    suffix=" aa"
                    hold_copy="YES"
		else
		    prefix=$(echo $posloc | sed -n 's| .*||p')
                    suffix=" i"
		fi
	    fi
	elif [[ $job_status == "QUEUING" ]]; then 
	    if [[ $supercomputer == "NO" ]]; then
		if [[ $forced_location == "ISCA" ]]; then
                    prefix=$(echo $posloc | sed -n 's| .*||p')
                    suffix=" ii"
		    hold_copy="YES"
		elif [[ $forced_location == "ARCHER" ]]; then 
		    prefix=$(echo $posloc | sed -n 's| .*||p')
                    suffix=" aa"   
		    hold_copy="YES"
		else
		    prefix=$(echo $posloc | sed -n 's| .*||p')
		    suffix=" r"
		fi
	    elif [[ $supercomputer == "YES" ]]; then 
		if [[ $forced_location == "ISCA" ]]; then
                    prefix=$(echo $posloc | sed -n 's| .*||p')
                    suffix=" ii"
		    if [[$is_isca != "ISCA" ]]; then 
			hold_copy="YES"
		    fi
		elif [[ $forced_location == "ARCHER" ]]; then 
		    prefix=$(echo $posloc | sed -n 's| .*||p')
		    suffix=" aa"
		    if [[ $is_archer != "ARCHER" ]]; then 
			hold_copy="YES"
		    fi
		elif [[ $is_archer == "ARCHER" ]]; then 
		    echo "Cannot submit no authorised job to archer. Check your input file"
		elif [[ $is_isca == "ISCA" ]]; then 
		    prefix=$(echo $posloc | sed -n 's| .*||p')
                    suffix=" i"
		fi
	    fi
	fi
    
    
	queuing_variable=$prefix$suffix
	echo $queueing_variable
	
	sed -e "s|${posloc}|${queuing_variable}|g" </scratch/Jobq.txt >>/scratch/Jobqp.txt
        mv /scratch/Jobqp.txt /scratch/Jobq.txt
    elif [[ $mode == "MAINTAIN" ]]; then 
	sed -e "s|${posloc}|${posloc}|g" </scratch/Jobq.txt >>/scratch/Jobqp.txt
	mv /scratch/Jobqp.txt /scratch/Jobq.txt
    fi
}

test_scf_done () {

    if [[ $job_status == "QUEUING" ]]; then 
	dovasp="TRUE" 
    elif [[ $job_status == "LOCAL" ]]; then
	determine_appropriate_directory LOCAL
	scf_energy=$(grep "e  e" $output/OUTCAR)
	if [[ -z $scf_energy ]]; then 
	    dovasp="TRUE"
	    finished="NO"
	else
	    finished="YES"
	    dovasp="FINISHED"
	fi
    elif [[ $job_status == "ISCA" ]]; then
	determine_appropriate_directory ISCA
	echo $output "!!!!!!!!!!!!!!!!!!!!!!"
	scf_energy=$(ssh jp552@login.isca.ex.ac.uk "grep 'e  e' $output/OUTCAR")
	echo $scf_energy
	if [[ -z $scf_energy ]]; then 
	    dovasp="TRUE"
	else 
	    finished="YES"
	fi
     elif [[ $job_status == "ARCHER" ]]; then
        determine_appropriate_directory ARCHER
        echo $output "!!!!!!!!!!!!!!!!!!!!!!"
        scf_energy=$(ssh jp552@login.archer2.ac.uk "grep 'e  e' $output/OUTCAR")
        echo $scf_energy
        if [[ -z $scf_energy ]]; then
            dovasp="TRUE"
        else
            finished="YES"
        fi
    else
	echo "Fake jobs"
	exit
    fi
}



build_supercomputer_directory () {
    
    
    if [[ $is_isca == "ISCA" ]]; then 
	determine_appropriate_directory ISCA
	step1=$(echo $output | sed -n 's|/RELAX.*||p')

	if [[ -z $step1 ]]; then 
	    ssh jp552@login.isca.ex.ac.uk "mkdir $output" 
	else 
	    ssh jp552@login.isca.ex.ac.uk "mkdir $step1"
	    step2=$(echo $output | sed -n 's|/RELAX2.*||p')
	    if [[ -z $step2 ]]; then
		ssh jp552@login.isca.ex.ac.uk "mkdir $output"
	    else
		ssh jp552@login.isca.ex.ac.uk "mkdir $step2"
		ssh jp552@login.isca.ex.ac.uk "mkdir $output"
	    fi
	fi
    elif [[ $is_archer == "ARCHER" ]]; then 
	
	determine_appropriate_directory ARCHER
	step1=$(echo $output | sed -n 's|/RELAX.*||p')
	if [[ -z $step1 ]]; then
            ssh jp552@login.archer2.ac.uk "mkdir $output"
        else
            ssh jp552@login.archer2.ac.uk "mkdir $step1"
            step2=$(echo $output | sed -n 's|/RELAX2.*||p')
            if [[ -z $step2 ]]; then
                ssh jp552@login.archer2.ac.uk "mkdir $output"
            else
                ssh jp552@login.archer2.ac.uk "mkdir $step2"
                ssh jp552@login.archer2.ac.uk "mkdir $output"
            fi
        fi
    fi
}

copy_files () { 
    FILE=$working_directory"/CONTCAR"
    echo $FILE $job_status $dovasp $supercomputer
    if [[ $job_status == "ISCA" ]]; then 
	if [[ $is_isca == "ISCA" ]]; then 
	    :
	    #We don't need to copy anything here, but we do need to move the CONTCAR to the POSCAR
	    if [[ $dovasp == "FINISHED" ]]; then 
		determine_appropriate_directory ISCA
		scp -r jp552@login.archer2.ac.uk:$output/* $working_directory/.
	    else
		echo $FILE
		ssh jp552@login.isca.ex.ac.uk 'if [ -s "$FILE" ] ; then cd $working_directory; cat CONTCAR >> CONTCAR_HISTORY; cat POSCAR >> POSCAR_HISTORY; 
 mv CONTCAR POSCAR; fi' 	 
		
	    fi
	elif [[ $is_isca != "ISCA" ]]; then 
	    # We need to move from ISCA to the local machines
	    # Create hypothetical location for filepath source
	    # NO FACILITY TO MOVE FROM ISCA TO ARCHER
	    determine_appropriate_directory ISCA
	    
	    scp jp552@login.isca.ex.ac.uk:$output/* $working_directory/.
	    #delete_supercomputer_data
	    if [[ $SCF == "NO" ]]; then 
		if [ -s "$FILE" ]; then 
		    if [[ $dovasp == "FINISHED" ]]; then 
			 determine_appropriate_directory ISCA 
			 scp -r jp552@login.isca.ac.uk:$output/* $working_directory/.
		    else 
			cat $working_directory/CONTCAR >>  $working_directory/CONTCAR_HISTORY
			cat $working_directory/POSCAR >>  $working_directory/POSCAR_HISTORY
			mv $working_directory/CONTCAR $working_directory/POSCAR 
		    fi
		else
		    : 
		fi
	    fi
	fi
    elif [[ $job_status == "ARCHER" ]]; then
        if [[ $is_archer == "ARCHER" ]]; then
            :
            #We don't need to copy anything here, but we do need to move the CONTCAR to the POSCAR
            if [[ $dovasp == "FINISHED" ]]; then
                echo $working_directory 
		echo $job_status
		determine_appropriate_directory ARCHER 
		scp -r jp552@login.archer2.ac.uk:$output/* $working_directory/.
            else

                ssh jp552@login.archer2.ac.uk 'if [ -s "$FILE" ] ; then cd $working_directory; cat CONTCAR >> CONTCAR_HISTORY; cat POSCAR >> POSCAR_HISTORY;\

 mv CONTCAR POSCAR; fi'

            fi
        elif [[ $is_archer != "ARCHER" ]]; then
            # We need to move from ARCHER to the local machines
            # Create hypothetical location for filepath source
            determine_appropriate_directory ARCHER

           scp jp552@login.archer2.ac.uk:$output/* $working_directory/.
            #delete_supercomputer_data
            if [[ $SCF == "NO" ]]; then
                if [ -s "$FILE" ]; then
                    if [[ $dovasp == "FINISHED" ]]; then
                        :
                    else
                        cat $working_directory/CONTCAR >>  $working_directory/CONTCAR_HISTORY
                        cat $working_directory/POSCAR >>  $working_directory/POSCAR_HISTORY
                        mv $working_directory/CONTCAR $working_directory/POSCAR
                    fi
                else
                    :
                fi
            fi
        fi


    elif [[ $job_status == "LOCAL" ]]; then 
	if [[ $supercomputer == "NO" ]]; then 
	    : 
	    #No copy required#
	    if [[ $dovasp == "FINISHED" ]]; then 
		: 
	    else
		if [ -s "$FILE" ]; then
		    cat $working_directory/CONTCAR >>  $working_directory/CONTCAR_HISTORY
                    cat $working_directory/POSCAR >>  $working_directory/POSCAR_HISTORY
		    mv $working_directory/CONTCAR $working_directory/POSCAR
		fi
	    fi
	elif [[ $is_isca == "ISCA" ]]; then
	    if [[ $dovasp == "FINISHED" ]]; then 
		: 
	    else
		if [ -s "$FILE" ]; then
		    cat $working_directory/CONTCAR >>  $working_directory/CONTCAR_HISTORY
                    cat $working_directory/POSCAR >>  $working_directory/POSCAR_HISTORY
		    mv $working_directory/CONTCAR $working_directory/POSCAR
		    
		fi
	    fi
	    #Move from local machine to ISCA 
	    build_supercomputer_directory
	    determine_appropriate_directory LOCAL
	    scp $output/* jp552@login.isca.ex.ac.uk:$working_directory/.
	elif [[ $is_archer == "ARCHER" ]]; then
            if [[ $dovasp == "FINISHED" ]]; then
                :
            else
                if [ -s "$FILE" ]; then
                    cat $working_directory/CONTCAR >>  $working_directory/CONTCAR_HISTORY
                    cat $working_directory/POSCAR >>  $working_directory/POSCAR_HISTORY
                    mv $working_directory/CONTCAR $working_directory/POSCAR

                fi
            fi
            #Move from local machine to ISCA
            build_supercomputer_directory
            determine_appropriate_directory LOCAL
            scp $output/* jp552@login.archer2.ac.uk:$working_directory/.

	fi
    fi
    if [[ $job_status == "QUEUING" ]]; then 
	if [[ $supercomputer == "YES" ]]; then
	    build_supercomputer_directory
	    determine_appropriate_directory LOCAL
	    if [[ $is_archer == "ARCHER" ]]; then 
		scp $output/* jp552@login.archer2.ac.uk:$working_directory/.
	    elif [[ $is_isca == "ISCA" ]]; then
		scp $output/* jp552@login.isca.ex.ac.uk:$working_directory/.
	    fi
	elif [[ $supercomputer == "NO" ]]; then 
	    : 
	fi
    fi
    
}

test_dir_limitreached () {
    limit_reached="NO"
    if [[ $job_status == "QUEUING" ]]; then 
	limit_reached="NO"
    elif [[ $job_status == "ISCA" ]]; then 
	determine_appropriate_directory ISCA

	
	NSW=$(ssh jp552@login.isca.ex.ac.uk "sed -n 's|NSW = ||p' <$output/INCAR")
	echo $NSW 
	
	limit_reached=$(ssh jp552@login.isca.ex.ac.uk "sed -n \"s|${NSW} F|FOUNDLIMIT|p\" < $output/out.out")
	echo $limit_reached
	
	if [[ -z $limit_reached ]]; then 
	    # the limit of RELAX steps not reached. Could still be running
	    limit_reached="NO"
	else
	    # steps reached. we can rerun safely 
	    limit_reached="YES"
	    #HOWEVER, THE JOB COULD STILL BE IN THE QUEUE
	    #BUILD FIX FOR THIS
	fi
    elif [[ $job_status == "ARCHER" ]]; then 
	determine_appropriate_directory ARCHER

        NSW=$(ssh jp552@login.archer2.ex.ac.uk "sed -n 's|NSW = ||p' <$output/INCAR")
	limit_reached=$(ssh jp552@login.archer2.ex.ac.uk "sed -n \"s|${NSW} F|FOUNDLIMIT|p\" < $output/out.out")

        if [[ -z $limit_reached ]]; then
            # the limit of RELAX steps not reached. Could still be running
            limit_reached="NO"
        else
            # steps reached. we can rerun safely
            limit_reached="YES"
            #HOWEVER, THE JOB COULD STILL BE IN THE QUEUE
            #BUILD FIX FOR THIS
        fi
    elif [[ $job_status == "LOCAL" ]]; then 
	determine_appropriate_directory LOCAL 
	#NSW=$(grep 'NSW = ' $output/INCAR | sed -n 's|.* ||p')
	echo $(sed -n 's|NSW = ||p' <$output/INCAR) >temp_output.txt
	NSW=$(<temp_output.txt)
	echo $NSW
	rm temp_output.txt
	limit_reached=$(sed -n "s|${NSW} F|FOUNDLIMIT|p" < $output/out.out)
	echo $limit_reached
	#limit_reached=$(grep '$NSW F' $output/out.out)
	if [[ -z $limit_reached ]]; then 
            limit_reached="NO"
	else
            limit_reached="YES"
	fi
    fi
}

test_dir_running () {
    if [[ $job_status == "QUEUING" ]]; then 
	running="NO"
    else
	:
	#prototype function which will later be able to search for real time estimate of if a job is actually being run or not 
    fi
}

test_dir_SCF () {
    if [[ $(echo $posloc | sed -e 's|.*RELAX.*||g') == "" ]]; then 
	SCF="NO"
    else
	SCF="YES" 
    fi
}

test_dir_errors () {
    if [[ $job_status == "QUEUING" ]]; then
	errors="NO"
    elif [[ $job_status == "ISCA" ]]; then
	determine_appropriate_directory ISCA
	echo $output "IGVGFCYFCJVKYCYCUVUVC"
	errors=$(ssh jp552@login.isca.ex.ac.uk "sed -e 1p<$output/e.out")
	echo $errors
	
	if [[ -z $errors ]]; then
            errors="NO"
	else
            errors="YES"
	fi
    elif [[ $job_status == "ARCHER" ]]; then 
	determine_appropriate_directory ARCHER
	errors=$(ssh jp552@login.archer2.ac.uk "sed -e 1p<$output/e.out")
	if [[ -z $errors ]]; then
            errors="NO"
        else
            errors="YES"
        fi
    elif [[ $job_status == "LOCAL" ]]; then
	determine_appropriate_directory LOCAL
	errors=$(sed -e 1p<$output/e.out)
	if [[ -z $errors ]]; then
            errors="NO"
	else
            errors="YES"
	fi
    else
	echo "Your job shouldn't exist. Fix it. Exiting"
	exit
    fi
}

test_dir_converged () {

    if [[ $job_status == "QUEUING" ]]; then 
	finished="NO"
    elif [[ $job_status == "ISCA" ]]; then 
	determine_appropriate_directory ISCA
	finished=$(ssh jp552@login.isca.ex.ac.uk "grep 'stopping structural' $output/out.out")
	if [[ -z $finished ]]; then 
	    finished="NO"
	else
	    finished="YES"
	fi
    elif [[ $job_status == "ARCHER" ]]; then
        determine_appropriate_directory ARCHER
        finished=$(ssh jp552@login.archer2.ac.uk "grep 'stopping structural' $output/out.out")
        if [[ -z $finished ]]; then
            finished="NO"
        else
            finished="YES"
        fi

    elif [[ $job_status == "LOCAL" ]]; then 
	determine_appropriate_directory LOCAL
	finished=$(grep 'stopping structural' $output/out.out)
	if [[ -z $finished ]]; then
            finished="NO"
	else
            finished="YES"
	fi
    else 
	echo "Your job shouldn't exist. Fix it. Exiting"
	exit
    fi 

}

determine_appropriate_directory () {

    force=$1
    if [[ $force == "NO" ]]; then 
	if [[ $supercomputer == "NO" ]]; then 
	    if [[ $(echo $posloc | sed -e 's|.*POSCAR_[0-9]*.*||g') == "" ]]; then
		working_directory=$(echo $posloc | sed -n 's| [a-z]*||p')
	    fi 
	elif [[ $supercomputer == "YES" ]]; then
	    if [[ $is_isca == "ISCA" ]]; then
		sc_homedir="/gpfs/ts0/projects/Research_Project-172804/DJP552"
		sc_suffix=$(echo $(echo $posloc | sed -n 's|.*pos|pos|p') | sed -n 's| [a-z]*||p')
		working_directory=$sc_homedir/$sc_suffix
		
	    elif [[ $is_archer == "ARCHER" ]]; then
		sc_homedir="/mnt/lustre/a2fs-work2/work/e05/e05/jp552"
		sc_suffix=$(echo $(echo $posloc | sed -n 's|.*pos|pos|p') | sed -n 's| [a-z]*||p')
		working_directory=$sc_homedir/$sc_suffix
	    fi
	fi
    elif [[ $force == "ARCHER" ]]; then 
	
	    sc_homedir="/mnt/lustre/a2fs-work2/work/e05/e05/jp552"
            sc_suffix=$(echo $(echo $posloc | sed -n 's|.*pos|pos|p') | sed -n 's| [a-z]*||p')
            output=$sc_homedir/$sc_suffix
	    echo $sc_homedir
	    echo $sc_suffix 
	    echo $output
    elif [[ $force == "ISCA" ]]; then 
 
	    sc_homedir="/gpfs/ts0/projects/Research_Project-172804/DJP552"
	    sc_suffix=$(echo $(echo $posloc | sed -n 's|.*pos|pos|p') | sed -n 's| [a-z]*||p')
	    output=$sc_homedir/$sc_suffix
	    
    
    elif [[ $force == "LOCAL" ]]; then 
	if [[ $(echo $posloc | sed -e 's|.*POSCAR_[0-9]*.*||g') == "" ]]; then
	    output=$(echo $posloc | sed -n 's| [a-z]*||p')
	fi
    fi
} 


jobcount () {
    FILE="/scratch/Jobq.txt"
    if test -f "$FILE"; then 
	echo "Reading Job Q"
    else 
	touch /scratch/Jobq.txt
    fi
    t="0"
    while read jobs; do 
	if [[ -z $jobs ]]; then 
	    ### YOU WERE JUST ABOUT TO BUILD SOMETHING IN TO IGNORE BLANK SPACES, MAKE COUNTER EQUAL TO THE TOTAL JOBS BEFORE QUITTING
	    counter=$(echo "$counter +1" | bc -l)
	    if (( $(bc -l <<< "$counter" > 100) )); then 	    
		break
	    fi 
	fi
	t=$(echo "$t +1" | bc -l) 
    done </scratch/Jobq.txt 


}

SCFJobrun () {
    
    workdir=$(readlink -f .)
    jobcount
    p="1"
    
    
    H=$(hostname)
    vartmp=$(echo $H | sed -e 's/.ex.ac.uk//g')
    if [[ -z $vartmp ]]; then 
	var=$(hostname)
    else 
	var=$vartmp
    fi
    echo $var "test"
    resallowed
    for i in $(seq 1 $t)
    do
	logic="1"   
	while [[ $logic > 0 ]]; do
	    


	    if [[ $supercomputer == "NO" ]]; then 
		ssh $var "source ~/bin/Free.sh; n_free_resources"
		perc_free_cpu=$( sed -n '1p' < /gpshome/jp552/resources.txt)
		perc_free_ram=$( sed -n '2p' < /gpshome/jp552/resources.txt)
		user_jobs_running=$(ssh $var "top -b -n1 | grep 'vasp_std'"| grep 'jp552')
		total_jobs_running=$(ssh $var "top -b -n1 | grep 'vasp_std'")
		i_userjobs=0
		i_totaljobs=0
		while :
                do
                    if [[ $(echo $user_jobs_running | sed -n "s/vasp_std/vasp_std/p") == "" ]]; then
                        break
                    fi

                    tmpstring1=$(echo $user_jobs_running | sed -n "s/vasp_std//p")
                    user_jobs_running=$tmpstring1
                    i_userjobs_prime=$(echo "$i_userjobs+1" | bc -l)
                    i_userjobs=$i_userjobs_prime
                done
		while :
                do
                    if [[ $(echo $total_jobs_running | sed -n "s/vasp_std/vasp_std/p") == "" ]]; then
                        break
                    fi

                    tmpstring1=$(echo $total_jobs_running | sed -n "s/vasp_std//p")
                    total_jobs_running=$tmpstring1
                    i_totaljobs_prime=$(echo "$i_totaljobs+1" | bc -l)
                    i_totaljobs=$i_totaljobs_prime
                done
		i_freejobs=$(echo "$i_totaljobs-$i_userjobs" | bc -l)



	    elif [[ $supercomputer == "YES" ]]; then 
		previous_run_machine
		if [[ $is_isca == "ISCA" ]]; then
		    SCvar1=$(grep -oP " CD " <<< $(ssh jp552@login.isca.ex.ac.uk "squeue --me"))
		    SCvar2=$(grep -oP " PD " <<< $(ssh jp552@login.isca.ex.ac.uk "squeue --me"))
		    SCvar3=$(grep -oP " R " <<<$(ssh jp552@login.isca.ex.ac.uk "squeue --me"))
		elif [[ $is_archer == "ARCHER" ]]; then 
		    if [[ $forced_location != "ARCHER" ]]; then 
			echo $forced_location "(&&^)*&^*^%(&^%*)^"
			echo "skipping resource allocation"
		    else
			SCvar1=$(grep -oP " CD " <<< $(ssh jp552@login.archer2.ac.uk "squeue --me"))
			SCvar2=$(grep -oP " PD " <<< $(ssh jp552@login.archer2.ac.uk "squeue --me"))
			SCvar3=$(grep -oP " R " <<<$(ssh jp552@login.archer2.ac.uk "squeue --me"))
		    fi
		fi
		echo $SCvar1 "Completed job string"
		echo $SCvar2 "Pending string" 
		echo $SCvar3 "Running string"
		tmpstring="$SCvar1 $SCvar2 $SCvar3"
		cj=0
		rj=0
		qj=0
		while :
		do
		    if [[ $(echo $tmpstring | sed -n "s/CD/CD/p") == "" ]]; then
                        break
                    fi

		    tmpstring1=$(echo $tmpstring | sed -n "s/CD//p")
		    tmpstring=$tmpstring1
		    cjp=$(echo "$cj+1" | bc -l)
		    cj=$cjp 
		done 
		while :
                do
		    if [[ $(echo $tmpstring | sed -n "s/R/R/p") == "" ]]; then
                        break
                    fi
		    #echo $tmpstring
                    tmpstring1=$(echo $tmpstring | sed -n "s/R//p")
                    #echo $tmpstring1
		    tmpstring=$tmpstring1
                    rjp=$(echo "$rj+1" | bc -l)
                    rj=$rjp
                done
		while :
                do
		    if [[ $(echo $tmpstring | sed -n "s/PD/PD/p") == "" ]]; then
                        break
                    fi

		    tmpstring1=$(echo $tmpstring | sed -n "s/PD//p")
                    tmpstring=$tmpstring1
                    qjp=$(echo "$qj+1" | bc -l)
                    qj=$qjp
		done
		perc_free_cpu=$(echo "$qj+$rj" | bc -l)
	    fi 
	    echo $qj $rj $perc_free_cpu $supercomputer "!!!!!!!!!!!!!!!!!!!!"
	    if [[ $supercomputer == "NO" ]]; then 
		if (( $(bc -l <<<"$perc_free_cpu > $perc_minfree_cpu" ) )); then
		    if (( $(bc -l <<<"$perc_free_ram > $perc_minfree_ram" ) )); then
			
			if (( $(bc -l <<<"$max_jobs > $i_totaljobs" ) )); then 
			    
			    echo $posloc $i_totaljobs $max_jobs $perc_free_cpu $perc_minfree_cpu >> analysis.txt
			    tempo=$(echo $i"p")
			    echo $tempo
			    posloc=$(sed -n $tempo </scratch/Jobq.txt)
			    if [[ -z $posloc ]]; then
				continue 2
			    fi
			    echo $posloc
			    if [[ $(echo $posloc | sed -n "s|.* [a-z]||p") == "i" ]]; then 
				echo "THIS JOB IS FORCED TO SUPERCOMPUTER, CYCLING MACHINE"
			    elif [[ $(echo $posloc | sed -n "s|.* [a-z]||p") == "a" ]]; then
                                echo "THIS JOB IS FORCED TO SUPERCOMPUTER, CYCLING MACHINE"
                            else

			    vasprun
			    logic="0"
			    echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
			    
			    continue
			    fi
			fi
		    fi		
		fi
	    else 
		echo $perc_free_cpu $perc_minfree_cpu
		if (( $(bc -l<<<"$perc_free_cpu < $perc_minfree_cpu" ) )); then 
		    echo $perc_free_cpu $perc_minfree_cpu 
		    echo "YOU MANAGED TO CODE IT IN"
		    tempo=$(echo $i"p") 
		    echo $tempo 
		    posloc=$(sed -n $tempo </scratch/Jobq.txt)
		    if [[ -z $posloc ]]; then
                        continue 2
                    fi
		    if [[ $(echo $posloc | sed -n "s|.* [a-z]||p") == "a" ]]; then
                        echo "THIS JOB SHOULD BE FORCED TO ARCHER"
                            if [[ $is_archer == "ARCHER" ]]; then 
				vasprun
				logic="0"
				echo "SUPERCOMPUTINGPOWERHASGIVENYOUTHEOPPORTUNITYTOCALCULATEATAFASTERRATE"
			    else 
				:
			    fi
		    elif [[ $is_archer == "ARCHER" ]]; then 
			:
		    elif [[ $is_isca == "ISCA" ]]; then 
			if [[ $(echo $posloc | sed -n "s|.* [a-z]||p") == "i" ]]; then
			     vasprun
			     logic="0"
			     echo "SUPERCOMPUTINGPOWERHASGIVENYOUTHEOPPORTUNITYTOCALCULATEATAFASTERRATE"
			else
			    vasprun
                            logic="0"
                            echo "SUPERCOMPUTINGPOWERHASGIVENYOUTHEOPPORTUNITYTOCALCULATEATAFASTERRATE"
			fi
			
		    fi
		else 
		    echo "you fixed it"
		fi
	    fi
	    echo "cycling"
	    resallowed	
	    
	done 
    done
} 
cp /scratch/Jobq.txt REALQUEUE.txt
#while : 
#do 
    SCFJobrun 
#done


