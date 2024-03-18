module gen
  use constants, only: real12, pi
  use geom, only: get_volume, get_sphere_overlap
  use inputs, only: &
       vdW, volvar, minbond, maxbond,&
       sigma_bondlength, bins, vps_ratio, filename_host,&
       enable_self_bonding
  use help
  use atomtype
  use add_atom
  implicit none 
contains 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine generation(leng, atomlist, alistrep, spacelist, formula, structno, &
       options, eltot, elnames, stochio, elrad, c_cut, c_min) 
    integer, intent(in) :: c_cut, c_min
    type(unitcell), dimension(:), allocatable :: formula
    real(real12) :: bondmin,posneg, r,r2, meanvol, q, normvol, cellmultiplier, calc,sigma1, tmpval
    integer :: l,leng, i,b, j, k, x, y, z, m,p, structures, structno, prev_structures, modeselect, prevpos
    integer :: errorcounter, ecount, eltot, options, loopcounter, addedelements, addtest, bondminindex,tmpint
    integer, dimension(:), allocatable :: shapeA
    integer :: eltype, scan, bonding_number_correction, num_VOID
    !! box = initial untransformed cubic unit cell 
    !! a 0 0
    !! 0 b 0
    !! 0 0 c
    real(real12), allocatable, dimension(:,:) :: peakseparation, temppeaks
    real(real12), allocatable, dimension(:) :: peaksindividual

    real(real12), dimension(3,3) :: box
    real(real12), dimension(3) :: angle, spacelist, tmpvector, bestlocation, backuplocation
    real(real12), dimension(:,:,:,:), allocatable :: bondlist
    real(real12) ::  bondcutoff,connectivity,tmpangle, volmin, volmax, bondpro1, bondpro2, distribution, tmpvalue
    real(real12) :: anglecutoffupper, anglecutofflower, sigma2, bestlocationindex, agausssamp, peak1, tmpdistribution
    real(real12) :: angle_distribution, bond_distribution, normalisation_a, prob_void, prob_scan, prob_pseudo
    type (atom), dimension(:,:), allocatable :: atomlist, alistrep, alistrepp, copy_list
    character(1024) :: name, tmp, command,location
    character(3), dimension(:), allocatable :: elnames, sing_el, elnames_copy
    integer, dimension(:), allocatable :: elno, stochio, stochio_copy
    integer, dimension(:,:), allocatable :: nearneighbourmatrix
    logical :: dir_e, new_position_needed, placed
    real(real12), dimension(:,:,:), allocatable :: elrad
    real(real12), dimension(:,:), allocatable :: bondavg, bondminimum
    real(real12), dimension(:,:,:,:,:), allocatable :: results_matrix
    real(real12), dimension(:,:), allocatable :: tempmatrix
    character(1) :: equality_String

    !! The info file doesn't contain much of use yet. Could build in if relevant 
    open(11, file="Info")
    !! Atomlist contains the information about the positions of all atoms in ALL structures. This may cause issues 
    !! with memory when large numbers of structures are used, may consider breaking down into seperate iterations 
    !! (e.g paralyse) 
    loopcounter=0
    addtest=0
    modeselect=options
    !write(*,*) "elrad=",elrad
    !used to test a potential atomic location for interactions with itself
    allocate(atomlist(structno,leng))
    allocate(alistrep(structno,leng*27))
    prevpos=structurecounter("pos")
    addedelements=0
    if(modeselect.eq.2) then 
       !write(*,*) "Initialise host"
       deallocate(alistrep)
       deallocate(atomlist)
       !write(*,*) "Deallocation completed"
       call initialisehost(leng)
       allocate(atomlist(structno,leng))
       allocate(alistrep(structno,leng*27)) 

       call addhost(eltot,structures,formula,location,atomlist,stochio,&
            &elnames,addedelements,name,structno)


    end if
    !! alistrep contains positions of all atoms repeated in adjacent unit cells. Same point as above. 

    !! calls the function structurecounter, which provides information about the number of currently existing 
    !! structures in the directory
    prev_structures=structurecounter("pos")
    !! assigns the length of elno to eltot. NOT SURE WHY, SHOULD BE LENG?. UNLESS ELNO CONTAINS ALL MATERIAL SPECS
    allocate(elno(eltot))
    allocate(elrad(4,eltot,eltot))
    !if(modeselect.ne.2) then 
    call chemread(elnames,eltot,elrad)
    !end if

    !call chemread(elnames,eltot,elrad)

    !! bondlist is a list of ALL the bonds between all the atoms and each of it's neighbours in the first tier of recursive repeated unit cells
    allocate(bondlist(leng,leng*27,eltot,eltot))
    allocate(bondavg(eltot,eltot))
    !! modeselect=1 is a special option allowing a new poscar to be added in at user specification


    !! Could implement structno>1 in the future for large imports
    if(modeselect.eq.1) structno=1
    structures=1
    loopcounter=0

    !! elnames is a 1D array containing the symbol for each of the atoms (length=eltot). [generator;eltot~>main;elno]----------------------------------------------------------!                                                        !
    !-----------------------------------------------------------------------------------------!

    inquire(file="iso", exist=dir_e) 
    if(dir_e) then 
    else 
       write(command,*) "mkdir iso" 
       call execute_command_line(command) 
    end if

!!! This section prepares isolation calculations
    allocate(copy_list(2,2))
    do i=1, eltot 

       write(name,'(A11,A,A7)')"iso/POSCAR_",trim(adjustl(elnames(i))),"/POSCAR"
       !!Calculates the new structure number, and writes it to tmp
       write(tmp,'(A11,A3)')"iso/POSCAR_",elnames(i)
       !!Checks if a directory to contain that file exsts already (It should never exist, coul add warning) 
       inquire(file=tmp, exist=dir_e)
       if(dir_e) then 
       else
          !!Writes a command to create said directory
          write(command,'(A17,A3)')"mkdir iso/POSCAR_",elnames(i)
          call execute_command_line(command)

          open(10+i, file=name) 
          write(10+i,'(A)') "test"
          write(10+i, *) "1.00000000"
          write(10+i, *) "   ","20.0000000000000000","        ","0.0000000000000000","        ","0.0000000000000000"
          write(10+i, *) "   ","0.0000000000000000","        ","20.0000000000000000","        ","0.0000000000000000"
          write(10+i, *) "   ","0.0000000000000000","        ","0.0000000000000000","        ","20.0000000000000000"
          write(10+i,*) "     ", elnames(i) 
          write(10+i,*) "          ", "1"
          write(10+i,*) "Direct"
          write(10+i,*) "        ","0.5000000000000000","        ","0.5000000000000000","        ","0.5000000000000000"
          close(10+i)

          allocate(sing_el(1)) 
          sing_el(1)=elnames(i) 
          call potwrite(tmp,sing_el,1) 
          deallocate(sing_el)
          call Jobwrite(tmp,1,1,1)
          call Incarwrite(tmp,500, 20) !!The 500 here is nstep electronic
          !write(*,*) trim(adjustl(tmp))
          call chdir(trim(adjustl(tmp)))
          call execute_command_line("qsub.sh")
          write(name,*) "cp ../../job_vasp_isca.in ."
          call execute_command_line(name)
          call chdir("../../")

       end if
    end do











    !! Builds pos and don subfolders if they do not already exist 

    call touchpos()

    inquire(file="don", exist=dir_e) 
    if(dir_e) then 
    else 
       write(command,*) "mkdir don"
       call execute_command_line(command) 
    end if






!!! Create the directories for all of the POSCARS to be placed into !!!
    !! BIGLOOP is the parent loop for all procesess, generating one structure for each full completed iteration
    b=0
    BIGLOOP: do while(structures.le.structno)

       !! DO WE WANT STRUCTNO TO BE UPDATABLE LATER ON?
       !call invar(2,tmpdig,tmpels) 
       !eltot=tmpdig(1)
       !deallocate(tmpdig)
       !! DO WE WANT ELNAME TO BE UPDATABLE LATER ON?
       !call invar(4,tmpdig,tmpels)
       !do i=1, eltot 
       !   elnames(i)=tmpels(i)
       !end do
       !! DO WE WANT STOCHIO UPDATABLE LATER ON?
       !deallocate(stochio)
       !allocate(stochio(eltot))
!!! How many of each would you like
       !call invar(5,tmpdig,tmpels)
       !do i=1, eltot
       !   stochio(i)=tmpdig(i)
       !end do
       leng=0
       !! Total number of atoms 
       do i=1, eltot
          leng=leng+stochio(i)
       end do


       if(modeselect.eq.2) then 
          write(*,*) "Initialise host"
          deallocate(alistrep)
          deallocate(atomlist)
          write(*,*) "Deallocation completed"
          call initialisehost(leng)
          allocate(atomlist(structno,leng))
          allocate(alistrep(structno,leng*27)) 

          call addhost(eltot,structures,formula,location,atomlist,stochio,&
               &elnames,addedelements,name,structno)


       end if




       call execute_command_line("rm buildmap_testfile.txt") 
       b=structures
       do while (b.gt.7) 
          b=b-7
       end do
       !write(*,*) b
       call touchposdir(structures,prevpos)
       !------------------------------------------------------------------------------------------------!   
       !    !!This section will count the loop number efficiently if no other prints are used in loop  !
       loopcounter=loopcounter+1                                                                  !
       !     write(6, '(A1, 40X, A1)', advance='no') achar(13), achar(13)                               !
       !     write(6,'(I0.4)', advance='no') loopcounter                                                !
       !------------------------------------------------------------------------------------------------!


!!!-------------------------------------------------------------------------------------!!!
!!!Decides if pseudorandom volume will be larger or smaller than atomic volume summation!!!
!!!-------------------------------------------------------------------------------------!!!

       posneg=1
       call random_number(r)
       if(r.lt.0.5) then 
          posneg=-posneg
       else 
       end if

!!!-------------------------------------------------------------------------------------!!!
!!!Decides the total desired pseudorandom cell volume                                   !!! 
!!!-------------------------------------------------------------------------------------!!!


       !! Meanvol takes the atomic radius and calculates a guestimate for the total rough cell volume. NEED A BETTER METHOD
       !! FOR ACCOMPLISHING THIS
       !meanvol=4/3*2.2**3*pi*leng
       meanvol=0
       volmin=0
       volmax=0 
       k=0
       do i=1, eltot
          do j=1, eltot
             !      if(elrad(3,i,i).gt.elrad(3,i,j)) write(*,*) "This element is bonded to more of it's partners than they are to it"
             !      if(elrad(3,i,i).lt.elrad(3,i,j)) write(*,*) "This element is bonded to less of it's partners than they are to it"

          end do
       end do

       !write(*,*) elrad(2,i,i)
       k=0
       bonding_number_correction=0

       if(.not.enable_self_bonding) then 
          do k=1, eltot 
             bonding_number_correction=bonding_number_correction+stochio(k)**2
          end do
       end if

       normalisation_a=0

       do i=1, eltot 
          do j=1, eltot
             if(i.eq.j) then 
                normalisation_a=normalisation_a+dble(dble(stochio(i))/leng)**2
                write(*,*) dble(stochio(i)/leng)**2, elnames(i), elnames(j)

             else 
                if(j.lt.i) cycle
                normalisation_a=normalisation_a+dble(dble(stochio(i)+stochio(j))/leng)**2
                write(*,*) dble(dble(stochio(i)+stochio(j))/leng)**2, elnames(i), elnames(j)

             end if



          end do
       end do
       normalisation_a=1/normalisation_a




       k=0
       do i=1, eltot
          volmin=volmin+stochio(i)*(4.0/3.0)*pi*(elrad(2,i,i)**3)
          do j=1, eltot
             if(j.lt.i) cycle
             connectivity=real(vdW/100.0, real12)
             !write(*,*) connectivity 

             if(i.eq.j) then 
                if(bonding_number_correction.eq.0) then 
                   volmin=volmin-connectivity*0.5*normalisation_a*min((stochio(i)*elrad(3,i,j)),(stochio(j)*elrad(3,j,i)))&
                        &*get_sphere_overlap(elrad(2,i,i),elrad(2,j,j),elrad(1,i,j))

                else if(eltot.eq.1) then 
                   volmin=volmin-connectivity*0.5*normalisation_a*min((stochio(i)*elrad(3,i,j)),(stochio(j)*elrad(3,j,i)))&
                        &*get_sphere_overlap(elrad(2,i,i),elrad(2,j,j),elrad(1,i,j))

                end if
             else 
                !volmin=volmin-connectivity*(stochio(i)*elrad(3,i,j)+stochio(j)*elrad(4,i,j))*0.5*&
                !     &((dble(stochio(j)*stochio(i))/(leng**2))**(0.5)*&
                !     &get_sphere_overlap(elrad(2,i,i),elrad(2,j,j),elrad(1,i,j)))
                !write(*,*) leng**2-bonding_number_correction, stochio(i)*stochio(j)

                !volmin=volmin-(2*stochio(i)*stochio(j)/(leng**2-bonding_number_correction))*&
                !     &connectivity*(min(stochio(i)*elrad(3,i,j)&
                !     &,stochio(j)*elrad(3,j,i))*&
                !     &get_sphere_overlap(elrad(2,i,i),elrad(2,j,j),elrad(1,i,j)))

                volmin=volmin-connectivity*normalisation_a*min((stochio(i)*elrad(3,i,j)),(stochio(j)*elrad(3,j,i)))&
                     &*get_sphere_overlap(elrad(2,i,i),elrad(2,j,j),elrad(1,i,j))





             end if
             meanvol=volmin
          end do
          !  meanvol=meanvol+stochio(i)*((connectivity*elrad(1,i,i))+((1.0-connectivity)*elrad(2,i,i)))**3*(4/3)*pi 
       end do
       !!Adds or subtracts a small quantity from the calculated volume 
       call random_number(r)
       write(*,*) meanvol
       meanvol=meanvol+((volvar/100.0)*r*posneg*meanvol)
       write(*,*) "The allocated volume is", meanvol

!!!-------------------------------------------------------------------------------------!!!
       !!Define the random unit cell lengths                                                  !!!
!!!-------------------------------------------------------------------------------------!!!

       !! Initialises box, which is a cubic unit serving as the basis for the random unit cel
       box=0
       !! q keeps a running total of the "volume" in the loosest sense of the word  
       q=1

       call random_number(r)
       box(1,1)=0.75+r*2.25
       q=q*r
       call random_number(r)
       box(2,2)=0.75+r*2.25
       q=q*r
       call random_number(r)
       box(3,3)=0.75+r*2.25
       q=q*r*1000

!!!--------------------------------------------------------------!!!
!!!Sets the random angles between the unit vectors between 60-120!!!
!!!--------------------------------------------------------------!!!
       angle(:)=0
!!! BRAVAIS LATTICES 
       if (b.eq.1) then !!! Triclinic 
          do i=1, 3
             call random_number(r)
             r=r*60.0+60.0                                                            !   
             r=(r*pi)/(180.0)                                                         !     
             angle(i)=r 
          end do
       else if (b.eq.2) then !!! Cubic
          q=1
          call random_number(r) 
          box(1,1)=0.75+r*2.25
          q=q*(r**3)*1000.0
          box(2,2)=box(1,1) 
          box(3,3)=box(1,1)
          angle(:)=pi/2.0
       else if (b.eq.3) then !!! Monoclinic 
          angle(3)=pi/2
          angle(1)=pi/2
          call random_number(r) 
          r=r*60.0+60.0
          r=(r*pi)/(180.0)
          angle(2)=r
       else if (b.eq.4) then !!! Orthorhombic 
          angle(:)=pi/2.0
       else if (b.eq.5) then !!! Tetragonal B
          q=1 
          call random_number(r) 
          box(1,1)=0.75+r*2.25
          box(2,2)=box(1,1) 
          q=q*(r**2)
          call random_number(r) 
          box(3,3)=0.75+r*2.25 
          q=q*r*1000
          angle(:)=pi/2.0

       else if (b.eq.6) then !!! Rhombohedral very broken/ Trigonal :-(
          q=1 
          call random_number(r) 

          box(1,1)=0.75+r*2.25
          ! box(1,1)=5.0
          box(2,2)=box(1,1) 
          box(3,3)=box(1,1) 
          q=q*(r**3)*1000 
          !write(*,*) box(1,1)
          call random_number(r)
          r=(r*60.0)+60.0 
          r=(r*pi)/(180.0)
          angle(:)=r
          ! angle(:)= 1.75*pi/3.0
          !write(*,*) angle(:)
       else if (b.eq.7) then !!! hexagonal 
          angle(1)=pi/2.0
          angle(2)=pi/2.0
          angle(3)=2.0*pi/3.0
          call random_number(r) 
          box(1,1)=0.75+r*2.25
          box(2,2)=box(1,1) 
          q=r**2
          call random_number(r) 
          box(3,3)=0.75+r*2.25
          q=q*r*1000
       end if


!!!---------------------------------!!!
       !!Sets the new unit vectors         !!!
!!!---------------------------------!!!

       !! Creates a TYPE called formula that contains all the info for a unit cell that is random

       calc=sin(angle(1))**2-cos(angle(2))**2-cos(angle(3))**2
       calc=calc+cos(angle(1))*cos(angle(2))*cos(angle(3))*2
       calc=sqrt(calc)
       calc=calc/(sin(angle(1)))
       write(name,'(A11,I0.3,A7)')"pos/POSCAR_",structures+prevpos,"/POSCAR"
       open(structures+10000, file=name)
       write(structures+10000,*) "Test"
       write(structures+10000,*) 1.0

       if(modeselect.ne.2) then      
          formula(structures)%cell=0
          formula(structures)%cell(1,1)=box(1,1)*calc
          formula(structures)%cell(1,2)=box(1,1)*(cos(angle(3))-cos(angle(1))*cos(angle(2)))/(sin(angle(1)))
          formula(structures)%cell(1,3)=box(1,1)*cos(angle(2)) 
          formula(structures)%cell(2,2)=box(2,2)*sin(angle(1)) 
          formula(structures)%cell(2,3)=box(2,2)*cos(angle(1))      
          formula(structures)%cell(3,3)=box(3,3)

          !write(*,*) formula(structures)%cell

          !! Adjusts the the lengths of the unit cell vectors by the 
          normvol=get_volume(formula(structures)%cell)
          normvol=abs(normvol)/meanvol

          normvol=normvol**(1.0/3.0)
          if(modeselect.ne.2) then;
             do j=1, 3
                do i=1, 3

                   formula(structures)%cell(i,j)=formula(structures)%cell(i,j)/normvol      
                end do
                !        write(structures+10000,*) formula(structures)%cell(:,j)
             end do
          end if
       end if
!!! CHECKS IF THE UNIT CELL WILL FORCE ATOMS TO BE TOO CLOSE TOGETHER
       do i=1, eltot 
          do k=1, eltot
             do j=1, 3
                !write(*,*) (formula(structures)%cell(j,1)**2+formula(structures)%cell(j,2)**2+formula(structures)%cell(j,3)**2)**0.5, &
                !&elrad(1,k,i), k, i, elnames(k), elnames(i)
                if((formula(structures)%cell(j,1)**2+formula(structures)%cell(j,2)**2+formula(structures)%cell(j,3)**2)**0.5&
                     &.lt.0.9*(elrad(1,k,i))) then
                   write(*,*) "This unit cell would definitely cause recursive atoms to be closer together than 0.9 Elrad"
                   close(structures+10000)
                   cycle bigloop
                end if
             end do
          end do
       end do





!!!----------------------------------------------!!!
!!!Set everything to what it should be and places the atoms
!!!----------------------------------------------!!!

!!! Assigns all the atoms to the correct species label. 
       k=0
       z=0
       m=0
       l=0

       !! Old linear version
       !     do j=1, eltot
       !
       !        k=k+stochio(j)
       !        m=m+27*stochio(j)
       !        do i=1, leng
       !           if(modeselect.eq.1) exit         
       !           if((i.le.k).and.(i.gt.z)) then  
       !              atomlist(structures,i)%name=elnames(j)        
       !              atomlist(structures,i)%element_index=j
       !           end if
       !        end do
       !        do i=1, leng*27
       !           if(modeselect.eq.1) exit
       !           if((i.le.m).and.(i.gt.l)) then
       !              alistrep(structures,i)%name=elnames(j)
       !              alistrep(structures,i)%element_index=j
       !           end if
       !        end do!

       !        z=k
       !        l=m
       !     end do
       !     L=0
       !! New random version
       do j=1, eltot

          k=k+stochio(j)
          m=m+27*stochio(j)
          do i=1, leng
             if(modeselect.eq.1) exit
             if((i.le.k).and.(i.gt.z)) then
                atomlist(structures,i)%name=elnames(j)
                atomlist(structures,i)%element_index=j
             end if
          end do
          do i=1, leng*27
             if(modeselect.eq.1) exit
             if((i.le.m).and.(i.gt.l)) then
                alistrep(structures,i)%name=elnames(j)
                alistrep(structures,i)%element_index=j
             end if
          end do

          z=k
          l=m
       end do
       L=0
       !! Randomly swap 1000 pairs of elements. Tailor this number to suit individual needs, should never need to be much higher though  
       !! Do not swap host atoms, hence furst step 
       z=0
       do i=1, addedelements
          z=z+stochio(i)
          !write(*,*) addedelements, stochio(i) 

       end do
       !write(*,*) z

       do i=1, leng/2
          call random_number(r)
          call random_number(r2)
          r=r*(leng-z)+z
          r2=r2*(leng-z)+z
          !write(*,*) r, r2
          copy_list(1,1)=atomlist(structures,ceiling(r))
          copy_list(2,2)=atomlist(structures,ceiling(r2))

          atomlist(structures,ceiling(r))=copy_list(2,2) 
          atomlist(structures,ceiling(r2))=copy_list(1,1) 
       end do





       !L=addedelements
       !do i=1, eltot 
       !   do j=1, addedelements
       !      if(elnames(j).eq.(elnames(addedelements+i))) then 
       !         L=L-1
       !         write(*,*) L
       !      end if
       !   end do
       !end do
       !
       !addedelements=L


       allocate(stochio_copy(1+addedelements)) 
       allocate(elnames_copy(1+addedelements))
       stochio_copy=0
       elnames_copy=""

       do i=1, addedelements
          stochio_copy(i)=stochio(i)
          elnames_copy(i)=elnames(i)
       end do

       deallocate(stochio)
       deallocate(elnames)

       allocate(stochio(1+addedelements))
       allocate(elnames(1+addedelements))

       stochio=stochio_copy
       elnames=elnames_copy

       deallocate(stochio_copy)
       deallocate(elnames_copy)

       L=1

       do i=1+z, leng
          !write(*,*) i, z, "£"
          if(i.eq.1+z) then
             !write(*,*) "This is the first new atom, thus elnames", addedelements+1, "should be ", atomlist(structures,i)%name
             elnames(addedelements+1)=atomlist(structures,i)%name 
             !write(*,*) "As this is the first atom, stochio (", addedelements +1, ") = 1"
             stochio(addedelements+1)=1
          else
             if(atomlist(structures,i)%name.eq.atomlist(structures,i-1)%name) then 
                stochio(size(stochio))=stochio(size(stochio))+1
                !write(*,*) atomlist(structures,i)%name, "is the name of atom", i, " which is equal to", &
                !&atomlist(structures,i-1)%name, "which is name of atom", i-1
                !write(*,*) "Thus, the count of", size(stochio), " is incremented by", 1
             else 
                !write(*,*) "!"

                L=size(stochio)+1

                allocate(stochio_copy(L))
                allocate(elnames_copy(L))
                do j=1, L-1
                   stochio_copy(j)=stochio(j)
                   elnames_copy(j)=elnames(j)               
                end do
                stochio_copy(size(stochio_copy))=1
                elnames_copy(size(elnames_copy))=atomlist(structures,i)%name
                deallocate(stochio) 
                deallocate(elnames)
                allocate(stochio(size(stochio_copy)))
                allocate(elnames(size(elnames_copy)))
                elnames=elnames_copy
                stochio=stochio_copy
                deallocate(stochio_copy)
                deallocate(elnames_copy)

             end if
          end if
          !write(*,*) atomlist(structures,i)%name
       end do

       do i=l, leng
          do j=1, leng
             if(i.eq.j) cycle
             if(atomlist(structures,i)%name.eq.atomlist(structures,j)%name) then 
                if(atomlist(structures,i)%element_index.gt.atomlist(structures,j)%element_index) then 
                   atomlist(structures,i)%element_index=atomlist(structures,j)%element_index
                else 
                   atomlist(structures,j)%element_index=atomlist(structures,i)%element_index
                end if

             end if
          end do
       end do

       do i=l, leng*27
          do j=1, leng*27
             if(alistrep(structures,i)%name.eq.alistrep(structures,j)%name) then
                if(alistrep(structures,i)%element_index.gt.alistrep(structures,j)%element_index) then
                   alistrep(structures,i)%element_index=alistrep(structures,j)%element_index
                else
                   alistrep(structures,j)%element_index=alistrep(structures,i)%element_index
                end if

             end if
          end do
       end do

       !do i=1, leng 
       !   write(*,*)  atomlist(structures,i)%element_index,  atomlist(structures,i)%name
       !end do

       !do i=1, leng*27
       !   write(*,*) alistrep(structures,i)%element_index
       !end do


       !write(*,*) elnames
       !write(*,*) stochio
       eltot=maxval(atomlist(structures,:)%element_index)

       L=0




       if(modeselect.eq.2) then;
          call chemread(elnames,eltot,elrad)
          do i=1, addedelements
             L=L+stochio(i)
          end do
          !write(*,*) L
          do i=1,L 
             !write(*,*) atomlist(structures,i)%position(:) 
          end do

          do i=1, L
             call atomrepeater(structures,atomlist(structures,i)%position,&
                  &atomlist(structures,i)%name,&
                  &alistrep,formula,i,leng)
          end do
          do i=1, L*27
             !write(*,*) alistrep(structures,i)%position, L*27 
          end do

          !end if
       end if
!!!!!!!!!!!!! This section places the atoms via distribution !!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !i controls which atom is being placed. it will be reduced to 0 if the first atom is not to be randomly placed. Bad for small cells
       i=1
       !random_seed determines if for a host calculation the first atom should be randomly placed. 1 results in random seeding

       if(modeselect.ne.2) then 
          !write(*,*) formula(structures)%cell
          do j=1, 3
             call random_number(r) 
             atomlist(structures,1+L)%position(j)=r
             alistrep(structures,1+L)%position(j)=r
          end do
          errorcounter=0
          atomlist(structures,1+L)%position(:)=matmul(formula(structures)%cell,atomlist(structures,1+L)%position(:))
          call atomrepeater(structures,atomlist(structures,i)%position,&
               &atomlist(structures,i)%name,alistrep,formula,i,leng)
       else if(modeselect.eq.2) then         
          i=i-1
       end if


       b=minbond
       !write(*,*) dble(b/100.0)
       open(99,file="errorfile") 

       sigma1=sigma_bondlength
       ! This next line is intended to auto-tune the resolution of the guassian sampling. Needs testing
       !    sigma2=minval(peakseparation)/(2.0*sqrt(2.0*LOG(4.0)))
       sigma2=0.5
       i=i+L
       bondcutoff=2
       anglecutofflower=0
       anglecutoffupper=180
       scan=0
       bestlocationindex=0
       agausssamp=sigma2

       ! First pass is bin size, which should be tied to gaussian size of evolved functions

       !call execute_command_line("rm buildmap_testfile.txt") 
       num_VOID=1
       aloop: do while (i.le.leng-1)
          bondmin=10.0
          tmpvector=0
          write(*,*) i, "$$$$$"
          do j=1, 3
             call random_number(r)
             tmpvector(j)=r
          end do
          y=0
          j=0
          scan=1
          allocate(results_matrix(bins(1)+1,bins(2)+1,bins(3)+1,4,eltot))

          !if(i.ne.1) then 
          call random_number(r)
          !else 
          !   r=1.0
          !end if

          prob_void=real(vps_ratio(1)/(vps_ratio(1)+vps_ratio(2)+vps_ratio(3)),real12)
          prob_pseudo=prob_void+&
               real(vps_ratio(2)/(vps_ratio(1)+vps_ratio(2)+vps_ratio(3)),real12)
          prob_scan=prob_pseudo+&
               real(vps_ratio(3)/(vps_ratio(1)+vps_ratio(2)+vps_ratio(3)),real12)


          !        if(r.le.ratio_voidscan/100) then 
          !           write(*,*) "ADD ATOM VOID"

          !           call add_atom_void (bins,formula,i, sigma1,&
          !                &structures,sigma2,elnames,eltot,bondcutoff,atomlist,alistrep,tmpvector,elrad,leng)
          !           num_VOID=num_VOID+1
          !        else if(r.gt.ratio_voidscan/100) then 
          !           write(*,*) num_VOID, "VOID THING"
          !           call add_atom_scan_2 (bins, formula, atomlist, alistrep, i, structures, elrad,&
          !                &leng,results_matrix,eltot,elnames,placed,num_VOID)
          !           num_VOID=1
          !           if(placed.eqv..FALSE.) then
          !              write(*,*) "ADD ATOM VOID"
          !              call add_atom_void (bins,formula,i, sigma1,&
          !                   &structures,sigma2,elnames,eltot,bondcutoff,atomlist,alistrep,tmpvector,elrad,leng)
          !            end if
          !         end if

          write(*,*) prob_void, prob_pseudo, prob_scan 

          if(r.le.prob_void) then 
             write(*,*) "ADD ATOM VOID"!

             call add_atom_void (bins,formula,i, sigma1,&
                  &structures,sigma2,elnames,eltot,bondcutoff,atomlist,alistrep,&
                  tmpvector,elrad,leng,c_cut)
             num_VOID=num_VOID+1
          else if(r.le.prob_pseudo) then 
             write(*,*) num_VOID, "Add Atom Pseudo"
             call add_atom_pseudo (bins, formula, atomlist, alistrep, i, structures, elrad,&
                  &leng,results_matrix,eltot,elnames,placed,num_VOID)
             num_VOID=1
             placed=.TRUE.
             if(placed.eqv..FALSE.) then
                write(*,*) "ADD ATOM VOID"
                call add_atom_void (bins,formula,i, sigma1,&
                     &structures,sigma2,elnames,eltot,bondcutoff,atomlist,alistrep,&
                     tmpvector,elrad,leng,c_cut)
             end if

          else if(r.le.prob_scan) then 
             write(*,*) num_VOID, "Add Atom Scan"
             call add_atom_scan_2 (bins, formula, atomlist, alistrep, i, structures, elrad,&
                  &leng,results_matrix,eltot,elnames,placed,num_VOID,c_cut,c_min)
             num_VOID=1
             placed=.TRUE.
             if(placed.eqv..FALSE.) then
                write(*,*) "ADD ATOM VOID"
                call add_atom_void (bins,formula,i, sigma1,&
                     &structures,sigma2,elnames,eltot,bondcutoff,atomlist,alistrep,tmpvector,elrad,leng,c_cut)
             end if
          end if









          deallocate(results_matrix)
          distribution=1
          !!Ditch next line for scan
          tmpvector=matmul(formula(structures)%cell,tmpvector)
          !! Implement input integers to mark probabilities or breakpoints 
          !if (i-L.le.-12) then 
          !   write(*,*) "PLACING ATOM RANDOMLY"
          !   call add_atom_random(formula,i,sigma1,structures,sigma2,elnames,eltot,&
          !        bondcutoff, atomlist, alistrep,tmpvector,elrad)
          !   atomlist(structures,i+1)%position=tmpvector
          !   write(*,*) "RANDOM ATOM SEEDED"
          !else if (i-L.le.-1000) then  
          !   call add_atom_scan(5,formula,i,sigma1,structures,sigma2,elnames,eltot,&
          !        bondcutoff, atomlist, alistrep,tmpvector,elrad,bestlocation)
          !   atomlist(structures,i+1)%position=bestlocation
          !else
          !   call add_atom_void10,formula,i,sigma1,structures,sigma2,elnames,eltot,&
          !        bondcutoff, atomlist, alistrep,tmpvector,elrad,bestlocation)
          !   atomlist(structures,i+1)%position=bestlocation
          !
          !end if


          !call atomrepeater(structures,atomlist(structures,i+1)%position,alistrep,formula,i+1,leng)

          if(i.eq.leng-1) exit 
          i=i+1
          bestlocationindex=0
          bestlocation=0

       end do aloop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !Build in some failsafes - changing sigma iteratively or change the probability width
       !can also change the minimum allowed bondlength now!! Wouldn't use this too often

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !Generate the atomic bonding information files used in upcoming learning algorithm
       write(*,*) eltot, leng
       prev_structures=structurecounter("bon")
       !do i=1, leng
       !   call generatebondfiles(1,structures,prev_structures,atomlist,alistrep,eltot,stochio,i)
       !   call generateanglefiles(1,structures,prev_structures,atomlist,alistrep,eltot,stochio,i,bondcutoff)
       !   call generate4files(1,structures,prev_structures,atomlist,alistrep,eltot,stochio,i,bondcutoff)

       !end do
       !do i=1, leng
       !write(*,*) atomlist(structures,i)%position, i
       !end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       j=0
       k=0

       !call invar(8,tmpdig,tmpels)
       !do i=1, eltot 
       !   do k=1, eltot 
       !      do x=1, leng*27
       !         if(modeselect.eq.1) exit
       !         do y=1, leng*27        
       !            if(x.ge.y) cycle
       !            if(alistrep(structures,x)%name.ne.elnames(k)) cycle 
       !            if(alistrep(structures,y)%name.ne.elnames(i)) cycle
       !TURN ME BACK ON
       !if(bondlength(alistrep(structures,x)%position,alistrep(structures,y)%position)&
       !    &.lt.(elrad(1,k,i)*dble(tmpdig(1)/100.0))) then
       !  j=j+1        
       !  close(structures+10000)
       !  write(*,*) "Terminating. Bond lengths too short" 
       !  deallocate(tmpdig)
       !cycle BIGLOOP
       !end if
       !         end do
       !      end do
       !   end do
       !end do
       !deallocate(tmpdig)
       !m=j
       !write(*,*) "1!!!!!!!!!!!!"

       !! This is defined by the MAXBOND input parameter
       !        call invar(9,tmpdig,tmpels)
       !        do i=1, eltot 
       !           do k=1, eltot
       !              if(i.gt.k) cycle
       !              do x=1, leng
       !                 
       !         if(modeselect.eq.1) exit
       !         do y=1, leng*27 
       !            if(atomlist(structures,x)%name.ne.elnames(k)) cycle
       !            if(alistrep(structures,y)%name.ne.elnames(i)) cycle
       !            bondlist(x,y,i,k)=bondlength(atomlist(structures,x)%position,alistrep(structures,y)%position)
       !            !bondlist(y,x,k,i)=bondlist(x,y,i,k)
       !            if(bondlist(x,y,i,k).gt.(dble(tmpdig(1)/100.0)*elrad(1,i,k))) then;
       !               bondlist(x,y,i,k)=0
       !               !bondlist(y,x,k,i)=0
       !               cycle
       !            end if
       !         end do
       !      end do
       !   end do
       !end do
       !deallocate(tmpdig)
       !write(*,*) "2!!!!!!!!!!!"

       !! The following, for each species pairing, generates an average bondlength. At the end of
       !! the first eltot loop, the value m can be read out to be the total bonding for that pairing 
       !bondavg=0
       !inquire(file="bonddata.txt", exist=dir_e)
       !     if(dir_e) then
       !        open(81,status="old",file="bonddata.txt", access="append")
       !     else
       !        open(81,status="new",file="bonddata.txt", access="append")
       !     end if
       !write(*,*) "3!!!!!!!!!!!!"


       !do i=1, eltot 
       !   do k=1, eltot
       !      if(i.gt.k) cycle
       !      m=0
       !      do x=1, leng
       !         if(atomlist(structures,x)%name.ne.elnames(k)) cycle
       !         if(modeselect.eq.1) exit
       !         do y=1, leng*27
       !            if(bondlist(x,y,i,k).gt.0.001) then
       !               if(alistrep(structures,y)%name.ne.elnames(i)) cycle
       !               bondavg(i,k)=bondavg(i,k)+bondlist(x,y,i,k) 
       !               m=m+1
       !            end if
       !         end do
       !      end do
       !      
       !      bondavg(i,k)=(dble(bondavg(i,k))/m)
       !      !bondavg(k,i)=bondavg(i,k)
       !      write(81,*) bondavg(i,k)
       !   end do
       !end do!

       !close(81)

       !write(*,*) "4!!!!!!!!!!!!"


       !allocate(tempmatrix(leng,leng*27))
       !
       !inquire(file="bondsfile.txt", exist=dir_e)
       !     if(dir_e) then
       !        open(81,status="old",file="bondsfile.txt", access="append")
       !     else
       !        open(81,status="new",file="bondsfile.txt", access="append")
       !     end if

       !do x=1, eltot 
       !   do y=1, eltot
       !      if(x.gt.y) cycle
       !      do i=1, leng 
       !         m=0
       !         do j=1, leng*27
       !            if(bondlist(i,j,x,y).gt.0.001) then;
       !               m=m+1
       !               write(81,*) bondlist(i,j,x,y),i,j,x,y
       !               !write(*,*) "!!"
       !               !tempmatrix(1,m)=bondlist(i,j,x,y)
       !            end if
       !         end do
       !      !write(*,*) minval(tempmatrix,MASK=tempmatrix.gt.0.01)
       !      !deallocate(tempmatrix)
       !      end do

       !   end do
       !end do
       !close(81)
       !        write(*,*) "5!!!!!!!!!!!!"

       !! The folowing calculates the difference between all bonds and their corresponding average 
       !! If that value is greater than a tolerance, the bonds are deemed too varied.
       !! This section is primed for deletion, as maybe unimportant
       !call invar(10,tmpdig,tmpels)
       !do i=1, eltot
       !   do k=1, eltot 
       !      do x=1, leng
       !         if(modeselect.eq.1) exit
       !         do y=1, leng*27
       !            q=abs(bondlist(x,y,i,k)-bondavg(i,k))
       !            if(q.gt.dble(tmpdig(1)/100.0)*bondavg(i,k)) then 
       !               !close(structures+10000) 
       !               !write(*,*) "Terminating. Bond lengths too varied"
       !               !deallocate(tmpdig)
       !               !cycle BIGLOOP
       !            end if
       !         end do
       !      end do
       !   end do
       !end do
       q=0
       !        write(*,*) "6!!!!!!!!!!!!"

       !deallocate(tmpdig)


       !do z=1, eltot 
       !   do l=1, eltot
       !      do i=1, leng*27 
       !         if(modeselect.eq.1) exit
       !         do k=1, leng*27
       !            do j=1, leng 
       !               if(i.eq.k) cycle
       !               if(alistrep(structures,i)%name.ne.elnames(z)) cycle
       !               if(alistrep(structures,j)%name.ne.elnames(l)) cycle
       !               if(bondlength(alistrep(structures,i)%position,atomlist(structures,j)%position).gt.bondavg(z,l)*1.2) cycle
       !               if(bondlength(alistrep(structures,k)%position,atomlist(structures,j)%position).gt.bondavg(z,l)*1.2) cycle
       !               if(bondlength(alistrep(structures,i)%position,atomlist(structures,j)%position).lt.0.001) cycle
       !               if(bondlength(alistrep(structures,k)%position,atomlist(structures,j)%position).lt.0.001) cycle

       !!NEEDS CORRECTLY IMPLEMENTING WITH MULTIPLE SPECIES

       !write(*,*) bondangle(alistrep(structures,i)%position,&
       !     &atomlist(structures,j)%position,alistrep(structures,k)%position)
       !if(bondangle(alistrep(structures,i)%position,&
       !     &atomlist(structures,j)%position,alistrep(structures,k)%position).lt.10) then 
       !   close(structures+10000) 
       !   write(*,*) "Terminating, tiny angles in system"
       !   cycle BIGLOOP
       !end if

       !if(bondangle(alistrep(structures,i)%position,&
       !     &atomlist(structures,j)%position,alistrep(structures,k)%position).gt.170) then 
       !   close(structures+10000) 
       !   write(*,*) "Terminating, big angles in system"
       !   cycle BIGLOOP
       !end if
       !            end do
       !         end do
       !      end do
       !   end do
       !end do
       !       write(*,*) "7!!!!!!!!!!!!"

!!!THIS SECTION ASSUMES AN AVERAGE COORDINATION NUMBER. THIS IS LIKELY VERY UNPHYSICAL!!! ALSO ASSUMES 2.5A BONDLENGTH
       !k=0
       !do i=1, leng   
       ! if(modeselect.eq.1) exit
       !   do j=1, leng*27
       !      if(bondlength(atomlist(structures,i)%position,alistrep(structures,j)%position).gt.(bondavg+0.2)) cycle 
       !     k=k+1
       !  end do
       !end do
       !k=nint(dble(k/leng))
       !m=0

       !! HYPER SPECIFIC COORDINATION SECTION 
       !if(k.ne.4) then
       !   close(structures+10000) 
       !   write(*,*) "Terminating. Coordination number incorrect"
       !   cycle bigloop
       !end if
!!!COORDINATION NUMBER NEEDS TO BE MORE SOPHISTICATED.
       !do i=1, leng   
       ! if(modeselect.eq.1) exit
       !   do j=1, leng*27
       !      if(bondlength(atomlist(structures,i)%position,alistrep(structures,j)%position).gt.(bondavg+0.2)) cycle 
       !      m=m+1
       !   end do
       !write(*,*) m, k
       !  if(m.gt.k+1) then 
       !     close(structures+10000) 
       !     write(*,*) "Terminating. Coordination number incorrect"
       !     cycle bigloop
       !  end if
       !  if(m.lt.k-1) then 
       !     close(structures+10000) 
       !     cycle bigloop
       !  end if
       !  m=0
       !end do
       !write(*,*) q


!!!! This changes the output of the unit cell size to the repeated lattice
       !   leng=leng*27 
       !   do i=1, 3
       !      do j=1, 3
       !         formula(structures)%cell(j,i)=formula(structures)%cell(j,i)*3
       !      end do
       !   end do

       do j=1, 3
          write(structures+10000,*) formula(structures)%cell(j,1), &
               formula(structures)%cell(j,2), &
               formula(structures)%cell(j,3)
          !write(*,*) formula(structures)%cell(j,1), &
          !     formula(structures)%cell(j,2), &
          !     formula(structures)%cell(j,3)

       end do


       !! Reorganises the POSCAR so all like atoms are next to each other, cutting down on space in the input line
       !! Identifies the string of elnames with doubles reduced
       L=0
       k=1
       allocate(elnames_copy(size(elnames)))
       elnames_copy=""
       elnames_copy(1)=elnames(1)
       equality_string=""
       envelope : do i=1, size(elnames)
          jloop: do j=1, size(elnames)

             if(j.lt.i) cycle
             if(elnames(j).eq.elnames(i)) then 
                L=l+1
             else 
                do k=1, size(elnames)
                   if (elnames(j) == elnames_copy(k) ) cycle jloop
                end do
                do k=1, size(elnames)
                   if( elnames_copy(k) == "" ) then  
                      elnames_copy(k)=elnames(j) 
                      !p!rint*, elnames(j), "loop2"
                      exit
                   end if
                end do


             end if

          end do jloop
       end do envelope
       deallocate(elnames) 
       allocate(elnames(size(elnames_copy)))

       elnames=elnames_copy

       !! Swap atoms to the correct positions

       do i=1, size(elnames)
          do j=1, leng
             if(atomlist(structures,j)%name.ne.elnames(i)) then
                do k=1, leng 
                   if(atomlist(structures,k)%name.eq.elnames(i)) then 
                      copy_list(1,1)=atomlist(structures,j)
                      copy_list(2,2)=atomlist(structures,k)

                      atomlist(structures,k)=copy_list(1,1) 
                      atomlist(structures,j)=copy_list(2,2)
                   end if
                end do
             end if
          end do
       end do





       k=size(elnames)
       do i=1, size(elnames)
          if(elnames(i) == "") then 
             k=k-1
          end if
       end do


       deallocate(elnames)
       allocate(elnames(k))
       do i=1, k
          elnames(i)=elnames_copy(i)
       end do
       deallocate(elnames)
       allocate(elnames(k))
       elnames=elnames_copy
       eltot=k

       deallocate(elnames_copy)


       !write(*,*) elnames, eltot, leng
       !write(*,*) atomlist



       call poswrite(formula(structures)%cell,atomlist,leng, structures, structno, prevpos)

       !call poswrite(formula(structures)%cell,atomlist,leng, structures, structno, prev_structures)

       write(tmp,'(A11,I0.3)')"pos/POSCAR_",prevpos+structures
       call Incarwrite(adjustl(tmp),500, 20*leng)
       !write(*,*) structures, prevpos, "!!!!!!!!!!!!!!"
       call Jobwrite(tmp,3,3,3)

       write(11,*) "For structure number", prevpos+structures
       write(11,*) "The average bond value is", bondavg
       write(11,*) "The lower bound for allowed bonds is", bondavg-0.2
       write(11,*) "The upper bound for allowed bonds is", bondavg+0.2
       call potwrite(tmp, elnames, eltot)
       close(structures+10000)
       structures=structures+1





       addtest=0
       call sleep(10)
    end do BIGLOOP
  end subroutine generation


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine chemread(elnames, eltot, elrad)
    character(3), dimension(:), allocatable :: elnames
    character(3) :: read1, read2
    real(real12) :: r_vdw, r_cov, c1, c2
    integer :: increment, i, eltot, j, Reason
    real(real12), dimension(:,:,:), allocatable :: elrad


    deallocate(elrad)
    allocate(elrad(4,eltot,eltot))
    elrad=0
    open(77, file="chem.in", status="old") 
    do 
       read(77, *, IOSTAT=Reason) read1, read2, r_cov, r_vdw, c1, c2
       if (Reason.gt.0) then;
          stop
          write(*,*) "Something wrong with chem.in file" 
       else if(Reason.lt.0) then; 
          write(*,*) "Done"
          exit 
       else 
          write(*,*) trim(adjustl(read1)), " ",trim(adjustl(read2))," ", r_cov, " ", r_vdw
          do i=1, eltot 
             do j=1, eltot
                if(elnames(i).eq.trim(adjustl(read1))) then;
                   if(elnames(j).eq.trim(adjustl(read2))) then; 
                      elrad(1,i,j)=r_cov
                      elrad(2,i,j)=r_vdw
                      elrad(3,i,j)=c1
                      elrad(1,j,i)=r_cov
                      elrad(2,j,i)=r_vdw
                      elrad(3,j,i)=c1
                      if(elnames(i).ne.elnames(j)) then; 
                         elrad(3,i,j)=c1
                         elrad(4,i,j)=c2
                         elrad(4,j,i)=c1
                         elrad(3,j,i)=c2

                         write(*,*) elnames(i),c1,",",elnames(j),c2
                      end if
                      continue
                   end if
                end if
                if(elnames(j).eq.trim(adjustl(read1))) then;
                   if(elnames(i).eq.trim(adjustl(read2))) then;
                      elrad(1,j,i)=r_cov
                      elrad(2,j,i)=r_vdw
                      elrad(3,j,i)=c1
                      if(elnames(i).ne.elnames(j)) then;
                         elrad(3,j,i)=c1
                         elrad(4,j,i)=c2
                         elrad(4,j,i)=c1
                         elrad(3,j,i)=c2

                      end if
                      continue

                   end if
                end if
             end do

          end do
       end if
    end do

    close(77)
  end subroutine chemread

  subroutine addhost(eltot,structures,formula,location,atomlist,stochio,elnames,&
       &addedelements,name,structno)
    type(unitcell), dimension(:), allocatable :: formula
    character(1024) :: tmp, name, location
    character(3), dimension(:), allocatable :: elnames,tmpelnames
    real(real12) :: cellmultiplier,meanvol,tmpdble
    integer :: d,l,k,j,i,structno, ecount,addedelements, eltot, leng,structures, prev_structures
    type (atom), dimension(:,:), allocatable :: tmplist,atomlist,tmplist2
    integer, dimension(:), allocatable :: stochio, tmpstochio, tmpstochiotot
    real(real12), dimension(3) :: tmparray


    !! Wipes the randomly generated formula
    close(61)

    open(61, file=trim(adjustl(filename_host)))
    read(61, *) tmp
    !write(*,*) tmp
    read(61, *) cellmultiplier 
    !write(*,*) cellmultiplier
    do i=1, 3
       !do while(.true.)
       read(61,*) tmparray
       write(*,*) tmparray
       !end do
       do d=1, structno       
          formula(d)%cell(i,:)=tmparray(:)
       end do
    end do

    !write(*,*) formula(1)%cell
    !! Probably incorrect to do this step here
    !meanvol=get_volume(formula(structures)%cell)
!!!!!!!!!!!!!
    ecount=0
    allocate(tmpelnames(100))
    allocate(tmplist(1,1000))

    addedelements=0

    read(61,'(A)') tmp

    !write(*,*) tmp
    do i=1,len(tmp)-1
       if(i.eq.1) then 
          if((scan(tmp(i:i+1)," ").eq.0).or.&
               &((scan(tmp(i:i+1)," ").eq.2))) then
             addedelements=addedelements+1
             eltot=eltot+1
             !write(*,*) tmp(i:i+1)
             if(scan(tmp(i:i+1)," ").eq.0) then
                tmpelnames(addedelements)=tmp(i:i+1)
             end if
             if((scan(tmp(i:i+1)," ").eq.2)) then
                tmpelnames(addedelements)=tmp(i:i)
             end if
          else
          end if
       else 
          if((scan(tmp(i:i+1)," ").eq.0).or.&
               &((scan(tmp(i:i+1)," ").eq.2).and.(scan(tmp(i-1:i)," ")&
               &.eq.1))) then
             addedelements=addedelements+1
             eltot=eltot+1
             !write(*,*) tmp(i:i+1)
             if(scan(tmp(i:i+1)," ").eq.0) then
                tmpelnames(addedelements)=tmp(i:i+1)
             end if
             if((scan(tmp(i:i+1)," ").eq.2).and.(scan(tmp(i-1:i)," ").eq.1)) then
                tmpelnames(addedelements)=tmp(i:i)
             end if
          else
          end if
       end if
    end do
    j=0
    do i=addedelements+1,eltot
       j=j+1
       tmpelnames(i)=elnames(j)
    end do
    deallocate(elnames)
    allocate(elnames(eltot)) 
    do i=1, eltot
       elnames(i)=tmpelnames(i)
    end do
    deallocate(tmpelnames)

    allocate(tmpstochio(addedelements))
    read(61,*) tmpstochio
    allocate(tmpstochiotot(eltot))

    do i=1,addedelements
       tmpstochiotot(i)=tmpstochio(i)
    end do


    do i=1, eltot-addedelements
       tmpstochiotot(i+addedelements)=stochio(i) 
    end do

    deallocate(stochio) 
    deallocate(tmpstochio) 
    allocate(stochio(eltot))
    do i=1, eltot
       stochio(i)=tmpstochiotot(i)
    end do
    k=0

    do i=1, eltot
       do j=1, stochio(i)
          k=k+1
          tmplist(1,k)%name=elnames(i)
       end do
    end do
    write(*,*) elnames
    l=0
    read(61,*) tmp   





    do i=1, addedelements
       do j=1, stochio(i)
          l=l+1
          !read(61,'(6X,F18.16)', advance='no') tmplist(structures,L)%position(1) 
          !read(61,'(6X,F18.16)', advance='no') tmplist(structures,L)%position(2)
          !read(61,'(6X,F18.16)') tmplist(structures,L)%position(3)


          read(61,*) tmplist(1,L)%position
          !                      write(*,*) tmplist(1,L)%position

          tmplist(1,L)%position=matmul(tmplist(1,L)%position,formula(1)%cell)
          write(*,*) tmplist(1,L)%position 
       end do
    end do
    !write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"


    do d=1, structno    
       l=0
       !write(*,*) d
       do i=1, eltot
          do j=1, stochio(i) 
             l=l+1
             !write(*,*) l
             if(i.le.addedelements) atomlist(d,L)%position=tmplist(1,L)%position
             atomlist(d,L)%name=tmplist(1,L)%name
             !write(*,*) atomlist(d,L)%position
          end do
       end do
    end do


    write(*,*) "Host accepted"
    deallocate(tmplist)

    close(61)


    !open(50,file=name) 
    !write(50,*) "Test"
    !write(50,*) 1.0 
    !do i=1, 3
    !   write(50,*) formula(structures)%cell(1,i),&
    !        &formula(structures)%cell(2,i),formula(structures)%cell(3,i)
    !end do
    !close(50)


  end subroutine addhost

  subroutine addxyzfile()
    type(unitcell), dimension(:), allocatable :: formula
    character(1024) :: name, buffer
    character(1024), dimension(:), allocatable :: elnames,elnames_tmp, elnames_list
    real(real12) :: cellmultiplier
    real(real12), dimension(:), allocatable :: bondcutoff
    integer :: tmp,q,l,k,j,i,structno, ecount, eltot, leng,structures, prev_structures
    type (atom), dimension(:,:), allocatable :: tmplist,atomlist, alistrep
    integer, dimension(:), allocatable :: stochio
    !! Wipes the randomly generated formula
    structures=1
    structno=1
    prev_structures=structurecounter("pos")
    call touchpos()
    call touchposdir(structures,prev_structures)
    allocate(formula(structno))

    write(6,*) "Please enter the filename you wish to add to the database"
    read(*, *) name
    !name="asio2_06.xyz"
    open(50, file=name)
    read(50, *) buffer
    read(50, *) cellmultiplier
    formula(1)%cell=0
    read(50, *) formula(1)%cell(1,1),formula(1)%cell(2,2),formula(1)%cell(3,3)
    write(name,'(A11,I0.3,A7)')"pos/POSCAR_",structures+prev_structures,"/POSCAR"
    open(structures+10000, file=name,status="new")
    write(structures+10000,*) "Test"
    write(structures+10000,*) 1.0
    do i=1, 3
       write(structures+10000,*)  formula(1)%cell(1,i),formula(1)%cell(2,i),formula(1)%cell(3,i)
    end do


!!!!WARNING ONLY FOR ONE SYSTEM USE WITH CAUTION

    tmp=3


    allocate(elnames_list(tmp))
    do i=1, tmp 
       read(50,*) elnames_list(i) 
    end do
    do i=1, tmp
       if(i.eq.1) then 
          eltot=1
          allocate(elnames(1))
          elnames(1)=elnames_list(1)
       else 
          loop2 :do j=1, eltot
             if(elnames_list(i).eq.elnames(j)) exit loop2 
             if(j.eq.eltot) then 
                allocate(elnames_tmp(eltot))
                elnames_tmp=elnames
                deallocate(elnames)
                eltot=eltot+1
                allocate(elnames(eltot))
                do k=1, eltot-1
                   elnames(k)=elnames_tmp(k)
                end do
                elnames(eltot)=elnames_list(i)
             end if
          end do loop2
       end if
    end do
    allocate(stochio(eltot))
    stochio=0
    do i=1, eltot 
       do j=1, tmp 
          if(elnames(i).eq.elnames_list(j)) stochio(i)=stochio(i)+1
       end do
    end do
    leng=0
    do i=1, eltot 
       leng=leng+stochio(i) 
    end do

    allocate(atomlist(1,leng))
    allocate(bondcutoff(eltot))
    write(*,*) "MANUALLY CHANGE"
    stop
    k=0
    do i=1, eltot
       rewind(50)
       read(50,*) 
       read(50,*) tmp
       read(50,*)
       do j=1, tmp
          !write(*,*) trim(adjustl(elnames(i))), trim(adjustl(elnames_list(j)))
          if(elnames(i).eq.elnames_list(j)) then
             k=k+1
             read(50,*) buffer, atomlist(1,k)%position
             atomlist(1,k)%name=elnames(i)
             !read(50, *) atomlist(1,k)%position
             !write(*,*) atomlist(1,k)%position
             !atomlist(1,k)%position(
          else 
             read(50,*)
          end if
       end do
    end do
    close(50)


    call poswrite(formula(1)%cell,atomlist,k,1,1,prev_structures)
    allocate(alistrep(1,k*27))
    do i=1, k 
       call atomrepeater(1,atomlist(1,i)%position,atomlist(1,i)%name,alistrep,formula,i,k)
    end do

    prev_structures=structurecounter("bon")


    !do i=1, k
    !   call generatebondfiles(1,1,prev_structures,atomlist,alistrep,eltot,stochio,i,elnames)
    !end do
    prev_structures=structurecounter("bad")
    bondcutoff=2.0
    !do i=1, k
    !   call generateanglefiles(1,1,prev_structures,atomlist,alistrep,eltot,stochio,i,bondcutoff)
    !end do


  end subroutine addxyzfile



end module gen
