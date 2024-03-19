module gen
  use constants, only: real12, pi
  use misc, only: touch
  use rw_geom, only: bas_type, geom_read
  use edit_geom, only: bas_merge
  use isolated, only: generate_isolated_calculations
  use vasp_file_handler, only: &
       touchposdir, &
       poswrite, incarwrite, jobwrite, generate_potcar, &
       addposcar
  use geom, only: get_volume, get_sphere_overlap, get_random_unit_cell
  use inputs, only: &
       vdW, volvar, minbond, maxbond,&
       sigma_bondlength, bins, vps_ratio, filename_host,&
       enable_self_bonding
  use help
  use atomtype
  use add_atom
  use read_chem, only: get_element_radius
  implicit none


contains 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine generation(alistrep, num_structures, &
       option, elnames, stochio, c_cut, c_min)
    implicit none
    integer, intent(inout) :: num_structures !! SHOULD NOT EVEN BE AN ARGUEMNT
    !! MAKE AN INPUT ARGUMENT THAT IS MAX_NUM_STRUCTURES
    integer, intent(in) :: option
    integer, intent(in) :: c_cut, c_min
    integer, dimension(:), allocatable, intent(inout) :: stochio

    type(bas_type) :: host_basis
    real(real12), dimension(3,3) :: host_lattice

    !! lattice = initial untransformed cubic unit cell 
    !! a 0 0
    !! 0 b 0
    !! 0 0 c
    real(real12), dimension(3,3) :: lattice
    type(bas_type), dimension(:), allocatable :: basis_list


    integer :: unit, info_unit, structure_unit
    integer :: bravais_type
    integer :: num_species, num_atoms

    real(real12) :: bondmin,posneg, r,r2, meanvol, q, normvol, calc,sigma1

    integer :: l, i, j, k, x, y, z, m,p
    integer :: structures, prev_structures, option_, prevpos, addedelements
    integer :: bonding_number_correction, num_VOID

    real(real12), dimension(3) :: angle, tmpvector, bestlocation
    real(real12) ::  bondcutoff, connectivity, volmin, distribution, sigma2
    real(real12) :: angle_distribution, bond_distribution, normalisation_a, prob_void, prob_scan, prob_pseudo
    type(atom), dimension(2) :: copy_list

    !! atomlist has dimensions(structure, atom)
    !! there is no distinguishing factor for species
    real(real12), dimension(:,:,:), allocatable :: elrad
    type (atom), dimension(:,:), allocatable :: atomlist, alistrep
    type(unitcell), dimension(:), allocatable :: formula
    character(1024) :: name, tmp, command,location
    character(3), dimension(:), allocatable :: elnames, elnames_copy
    integer, dimension(:), allocatable :: stochio_copy
    logical :: placed
    real(real12), dimension(:,:,:,:,:), allocatable :: results_matrix


    !! The info file doesn't contain much of use yet. Could build in if relevant 
    open(newunit=info_unit, file="Info")
    !! Atomlist contains the information about the positions of all atoms in ALL structures. This may cause issues 
    !! with memory when large numbers of structures are used, may consider breaking down into seperate iterations 
    !! (e.g paralyse)
    option_=option

    !!! THINK OF SOME WAY TO HANDLE THE HOST SEPARATELY
    !!! THAT CAN SIGNIFICANTLY REDUCE DATA USAGE
    num_species = size(elnames)
    num_atoms = sum(stochio)
    allocate(basis_list(num_structures))
    do i = 1, num_structures
       allocate(basis_list(i)%spec(num_species))
       basis_list(i)%spec(:)%name = elnames
       basis_list(i)%spec(:)%num = stochio
       basis_list(i)%natom = num_atoms
    end do
    select case(option_)
    case(2)
       call initialisehost(num_atoms, filename_host)
       allocate(atomlist(num_structures,num_atoms))
       allocate(alistrep(num_structures,num_atoms*27))
       open(newunit = unit, file = trim(adjustl(filename_host)))
       call geom_read(unit,host_lattice, host_basis)
       close(unit)
       basis_list = bas_merge(host_basis,basis_list(1))

       call addhost(num_species,structures,formula,location,atomlist,stochio,&
            &elnames,addedelements,name,num_structures)
    case default
       allocate(atomlist(num_structures,num_atoms))
       allocate(alistrep(num_structures,num_atoms*27))
       prevpos=structurecounter("pos")
       addedelements=0
    end select

    !! alistrep contains positions of all atoms repeated in adjacent unit cells. Same point as above. 

    !! calls the function structurecounter, which provides information about the number of currently existing 
    !! structures in the directory
    prev_structures=structurecounter("pos")
    allocate(elrad(4,num_species,num_species))
    elrad = get_element_radius(elnames)

    !! option_=1 is a special option allowing a new poscar to be added in at user specification


    !! Could implement num_structures>1 in the future for large imports
    if(option_.eq.1) num_structures=1
    structures = 1

    !! elnames is a 1D array containing the symbol for each of the atoms (length=num_species).

    !!--------------------------------------------------------------------------
    !! set up isolated element calculations
    !!--------------------------------------------------------------------------
    call generate_isolated_calculations(elnames)


    !!--------------------------------------------------------------------------
    !! create the 'pos' and 'don' directories
    !!--------------------------------------------------------------------------
    call touch("pos")
    call touch("don")


!!! Create the directories for all of the POSCARS to be placed into !!!
    !! BIGLOOP is the parent loop for all procesess, generating one structure for each full completed iteration
    bravais_type = 0
    BIGLOOP: do while(structures.le.num_structures)
!!!#############################################################################
!!!#############################################################################
!!! THIS SHOULD BE ITS OWN PROCEDURE THEN CALLED IN THE MAIN LOOP
!!!#############################################################################

       !! Total number of atoms 
       num_atoms = sum(stochio)

       !! WHY DEALLOCATE EVERY TIME?
       !! HANDLES THE CASE WHERE THE HOST IS TO BE USED
       if(option_.eq.2) then 
          write(*,*) "Initialise host"
          deallocate(alistrep)
          deallocate(atomlist)
          write(*,*) "Deallocation completed"
          call initialisehost(num_atoms, filename_host)
          allocate(atomlist(num_structures,num_atoms))
          allocate(alistrep(num_structures,num_atoms*27)) 

          call addhost(num_species,structures,formula,location,atomlist,stochio,&
               &elnames,addedelements,name,num_structures)
       end if


       call execute_command_line("rm buildmap_testfile.txt") 
       call touchposdir(structures,prevpos)


!!!-------------------------------------------------------------------------------------!!!
!!! Decides the total desired pseudorandom cell volume                                   !!! 
!!!-------------------------------------------------------------------------------------!!!
       !! Meanvol takes the atomic radius and calculates a guestimate for ...
       !! ... the total rough cell volume.
       !! NEED A BETTER METHOD
       !! FOR ACCOMPLISHING THIS
       !meanvol= 4/3 * pi * 2.2**3 * num_atoms
       meanvol = 0._real12

       bonding_number_correction = 0

       if(.not.enable_self_bonding) then 
          do k=1, num_species 
             bonding_number_correction=bonding_number_correction+stochio(k)**2
          end do
       end if
       
       !! calculate the normalisation factor
       normalisation_a = sum(stochio**2)
       do i=1, num_species 
          do j = i + 1, num_species, 1
             normalisation_a = normalisation_a + ( stochio(i) + stochio(j) )**2
          end do
       end do
       normalisation_a = ( num_atoms ** 2._real12 ) / normalisation_a


       !! calculate the minimum volume
       volmin = 0._real12
       do i = 1, num_species
          volmin = volmin + stochio(i) * (4._real12/3._real12) * pi * &
               ( elrad(2,i,i) ** 3._real12 )
          do j = i, num_species, 1
             connectivity = vdW / 100._real12
             if( ( i .eq. j ) .and. &
                  (bonding_number_correction.eq.0 .or. num_species .eq. 1)) then 
                volmin = volmin - connectivity * 0.5 * normalisation_a * &
                     min((stochio(i)*elrad(3,i,j)),(stochio(j)*elrad(3,j,i))) * &
                     get_sphere_overlap(elrad(2,i,i),elrad(2,j,j),elrad(1,i,j))
             else
                volmin = volmin - connectivity * normalisation_a * &
                     min((stochio(i)*elrad(3,i,j)),(stochio(j)*elrad(3,j,i))) * &
                     get_sphere_overlap(elrad(2,i,i),elrad(2,j,j),elrad(1,i,j))
             end if
             meanvol = volmin
          end do
          !  meanvol=meanvol+stochio(i)*((connectivity*elrad(1,i,i))+((1.0-connectivity)*elrad(2,i,i)))**3*(4/3)*pi 
       end do

       !!-----------------------------------------------------------------------
       !! sets pseudorandom volume to larger or smaller than atomic volume summation
       !!-----------------------------------------------------------------------
       posneg = 1
       call random_number(r)
       if(r.lt.0.5_real12) posneg = -posneg

       !!-----------------------------------------------------------------------
       !!Adds or subtracts a small quantity from the calculated volume
       !!-----------------------------------------------------------------------
       call random_number(r)
       write(*,*) meanvol
       meanvol = meanvol + ((volvar/100.0)*r*posneg*meanvol)
       write(*,*) "The allocated volume is", meanvol


       !!-----------------------------------------------------------------------
       !! generate random unit cell lengths
       !!-----------------------------------------------------------------------
       !!! ARE YOU SERIOUSLY SAYING THAT THE BRAVAIS LATTICE JUST CYCLES ...
       !!! ... FROM 1 TO 7 CONSTANTLY?
       !!! WHY IS IT NOT RANDOMLY GENERATED?
       bravais_type = structures
       do while (bravais_type.gt.7) 
          bravais_type = bravais_type - 7
       end do
       lattice = get_random_unit_cell(bravais_type, angle, q)


       !! Creates a TYPE called formula that contains all the info for a unit cell that is random

       calc = sin(angle(1))**2 - cos(angle(2))**2 - cos(angle(3))**2
       calc = calc + cos(angle(1)) * cos(angle(2)) * cos(angle(3)) * 2
       calc = sqrt(calc)
       calc = calc / sin(angle(1))

       if(option_.ne.2) then      
          formula(structures)%cell = 0._real12
          formula(structures)%cell(1,1)=lattice(1,1)*calc
          formula(structures)%cell(1,2)=lattice(1,1)*(cos(angle(3))-cos(angle(1))*cos(angle(2)))/(sin(angle(1)))
          formula(structures)%cell(1,3)=lattice(1,1)*cos(angle(2))
          formula(structures)%cell(2,2)=lattice(2,2)*sin(angle(1))
          formula(structures)%cell(2,3)=lattice(2,2)*cos(angle(1))
          formula(structures)%cell(3,3)=lattice(3,3)

          !! Adjusts the the lengths of the unit cell vectors by the 
          normvol=get_volume(formula(structures)%cell)
          normvol=abs(normvol)/meanvol

          normvol=normvol**(1._real12/3._real12)
          formula(structures)%cell = formula(structures)%cell / normvol
       end if


!!! CHECKS IF THE UNIT CELL WILL FORCE ATOMS TO BE TOO CLOSE TOGETHER
       do i=1, num_species 
          do k=1, num_species
             do j=1, 3
                !write(*,*) (formula(structures)%cell(j,1)**2+formula(structures)%cell(j,2)**2+formula(structures)%cell(j,3)**2)**0.5, &
                !&elrad(1,k,i), k, i, elnames(k), elnames(i)
                if((formula(structures)%cell(j,1)**2+formula(structures)%cell(j,2)**2+formula(structures)%cell(j,3)**2)**0.5&
                     &.lt.0.9*(elrad(1,k,i))) then
                   write(*,*) "This unit cell would definitely cause recursive atoms to be closer together than 0.9 Elrad"
                   cycle bigloop
                end if
             end do
          end do
       end do


!!!----------------------------------------------!!!
!!!Set everything to what it should be and places the atoms
!!!----------------------------------------------!!!

!!! Assigns all the atoms to the correct species label. 
       k = 0; z = 0; m = 0; l = 0
       do j = 1, num_species
          k=k+stochio(j)
          m=m+27*stochio(j)
          do i=1, num_atoms
             if(option_.eq.1) exit
             if((i.le.k).and.(i.gt.z)) then
                atomlist(structures,i)%name=elnames(j)
                atomlist(structures,i)%element_index=j
             end if
          end do
          do i=1, num_atoms*27
             if(option_.eq.1) exit
             if((i.le.m).and.(i.gt.l)) then
                alistrep(structures,i)%name=elnames(j)
                alistrep(structures,i)%element_index=j
             end if
          end do

          z=k
          l=m
       end do
       l = 0
       !! Randomly swap 1000 pairs of elements.
       !! Tailor this number to suit individual needs, should never need to be much higher though  
       !! Do not swap host atoms, hence furst step 
       z = sum(stochio(:addedelements))

       do i=1, num_atoms/2
          call random_number(r)
          call random_number(r2)
          r=r*(num_atoms-z)+z
          r2=r2*(num_atoms-z)+z
          !write(*,*) r, r2
          copy_list(1)=atomlist(structures,ceiling(r))
          copy_list(2)=atomlist(structures,ceiling(r2))

          atomlist(structures,ceiling(r))=copy_list(2) 
          atomlist(structures,ceiling(r2))=copy_list(1) 
       end do

       allocate(stochio_copy(1+addedelements)) 
       allocate(elnames_copy(1+addedelements))
       stochio_copy=0
       elnames_copy=""

       do i=1, addedelements
          stochio_copy(i)=stochio(i)
          elnames_copy(i)=elnames(i)
       end do

       deallocate(stochio, elnames)

       allocate(stochio(1+addedelements))
       allocate(elnames(1+addedelements))

       stochio=stochio_copy
       elnames=elnames_copy

       deallocate(stochio_copy, elnames_copy)

       do i=1+z, num_atoms
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

                L=size(stochio)+1

                allocate(stochio_copy(L))
                allocate(elnames_copy(L))
                do j=1, L-1
                   stochio_copy(j)=stochio(j)
                   elnames_copy(j)=elnames(j)               
                end do
                stochio_copy(size(stochio_copy))=1
                elnames_copy(size(elnames_copy))=atomlist(structures,i)%name
                deallocate(stochio, elnames)
                allocate(stochio(size(stochio_copy)))
                allocate(elnames(size(elnames_copy)))
                elnames = elnames_copy
                stochio = stochio_copy
                deallocate(stochio_copy, elnames_copy)

             end if
          end if
       end do

       do i=l, num_atoms
          do j=1, num_atoms
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

       do i=l, num_atoms*27
          do j=1, num_atoms*27
             if(alistrep(structures,i)%name.eq.alistrep(structures,j)%name) then
                if(alistrep(structures,i)%element_index.gt.alistrep(structures,j)%element_index) then
                   alistrep(structures,i)%element_index=alistrep(structures,j)%element_index
                else
                   alistrep(structures,j)%element_index=alistrep(structures,i)%element_index
                end if

             end if
          end do
       end do

       !! repeat atoms in adjacent unit cells
       num_species=maxval(atomlist(structures,:)%element_index)
       if(option_.eq.2) then;
          elrad = get_element_radius(elnames)
          do i=1, sum(stochio(:addedelements))
             call atomrepeater(structures,atomlist(structures,i)%position,&
                  &atomlist(structures,i)%name,&
                  &alistrep,formula,i,num_atoms)
          end do
       end if
!!!!!!!!!!!!! This section places the atoms via distribution !!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !i controls which atom is being placed. it will be reduced to 0 if the first atom is not to be randomly placed. Bad for small cells
       i=1
       !random_seed determines if for a host calculation the first atom should be randomly placed. 1 results in random seeding

       if(option_.ne.2) then 
          !write(*,*) formula(structures)%cell
          do j=1, 3
             call random_number(r) 
             atomlist(structures,1+L)%position(j)=r
             alistrep(structures,1+L)%position(j)=r
          end do
          atomlist(structures,1+L)%position(:)=matmul(formula(structures)%cell,atomlist(structures,1+L)%position(:))
          call atomrepeater(structures,atomlist(structures,i)%position,&
               &atomlist(structures,i)%name,alistrep,formula,i,num_atoms)
       else if(option_.eq.2) then         
          i=i-1
       end if

       sigma1=sigma_bondlength
       ! This next line is intended to auto-tune the resolution of the guassian sampling. Needs testing
       !    sigma2=minval(peakseparation)/(2.0*sqrt(2.0*LOG(4.0)))
       sigma2=0.5
       i=i+L
       bondcutoff=2

       ! First pass is bin size, which should be tied to gaussian size of evolved functions

       !call execute_command_line("rm buildmap_testfile.txt") 
       num_VOID=1
       aloop: do while (i.le.num_atoms-1)
          bondmin=10.0
          tmpvector=0
          write(*,*) i, "$$$$$"
          do j=1, 3
             call random_number(r)
             tmpvector(j)=r
          end do
          y=0
          j=0
          allocate(results_matrix(bins(1)+1,bins(2)+1,bins(3)+1,4,num_species))

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
          !                &structures,sigma2,elnames,num_species,bondcutoff,atomlist,alistrep,tmpvector,elrad,num_atoms)
          !           num_VOID=num_VOID+1
          !        else if(r.gt.ratio_voidscan/100) then 
          !           write(*,*) num_VOID, "VOID THING"
          !           call add_atom_scan_2 (bins, formula, atomlist, alistrep, i, structures, elrad,&
          !                &num_atoms,results_matrix,num_species,elnames,placed,num_VOID)
          !           num_VOID=1
          !           if(placed.eqv..FALSE.) then
          !              write(*,*) "ADD ATOM VOID"
          !              call add_atom_void (bins,formula,i, sigma1,&
          !                   &structures,sigma2,elnames,num_species,bondcutoff,atomlist,alistrep,tmpvector,elrad,num_atoms)
          !            end if
          !         end if

          write(*,*) prob_void, prob_pseudo, prob_scan 

          if(r.le.prob_void) then 
             write(*,*) "ADD ATOM VOID"!

             call add_atom_void (bins,formula,i, sigma1,&
                  &structures,sigma2,elnames,num_species,bondcutoff,atomlist,alistrep,&
                  tmpvector,elrad,num_atoms,c_cut)
             num_VOID=num_VOID+1
          else if(r.le.prob_pseudo) then 
             write(*,*) num_VOID, "Add Atom Pseudo"
             call add_atom_pseudo (bins, formula, atomlist, alistrep, i, structures, elrad,&
                  &num_atoms,results_matrix,num_species,elnames,placed,num_VOID)
             num_VOID=1
             placed=.TRUE.
             if(placed.eqv..FALSE.) then
                write(*,*) "ADD ATOM VOID"
                call add_atom_void (bins,formula,i, sigma1,&
                     &structures,sigma2,elnames,num_species,bondcutoff,atomlist,alistrep,&
                     tmpvector,elrad,num_atoms,c_cut)
             end if

          else if(r.le.prob_scan) then 
             write(*,*) num_VOID, "Add Atom Scan"
             call add_atom_scan_2 (bins, formula, atomlist, alistrep, i, structures, elrad,&
                  &num_atoms,results_matrix,num_species,elnames,placed,num_VOID,c_cut,c_min)
             num_VOID=1
             placed=.TRUE.
             if(placed.eqv..FALSE.) then
                write(*,*) "ADD ATOM VOID"
                call add_atom_void (bins,formula,i, sigma1,&
                     &structures,sigma2,elnames,num_species,bondcutoff,atomlist,alistrep,tmpvector,elrad,num_atoms,c_cut)
             end if
          end if



          deallocate(results_matrix)
          distribution=1
          !!Ditch next line for scan
          tmpvector=matmul(formula(structures)%cell,tmpvector)
          !! Implement input integers to mark probabilities or breakpoints 
          !if (i-L.le.-12) then 
          !   write(*,*) "PLACING ATOM RANDOMLY"
          !   call add_atom_random(formula,i,sigma1,structures,sigma2,elnames,num_species,&
          !        bondcutoff, atomlist, alistrep,tmpvector,elrad)
          !   atomlist(structures,i+1)%position=tmpvector
          !   write(*,*) "RANDOM ATOM SEEDED"
          !else if (i-L.le.-1000) then  
          !   call add_atom_scan(5,formula,i,sigma1,structures,sigma2,elnames,num_species,&
          !        bondcutoff, atomlist, alistrep,tmpvector,elrad,bestlocation)
          !   atomlist(structures,i+1)%position=bestlocation
          !else
          !   call add_atom_void10,formula,i,sigma1,structures,sigma2,elnames,num_species,&
          !        bondcutoff, atomlist, alistrep,tmpvector,elrad,bestlocation)
          !   atomlist(structures,i+1)%position=bestlocation
          !
          !end if


          !call atomrepeater(structures,atomlist(structures,i+1)%position,alistrep,formula,i+1,num_atoms)

          if(i.eq.num_atoms-1) exit 
          i=i+1
          bestlocation=0

       end do aloop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !Build in some failsafes - changing sigma iteratively or change the probability width
       !can also change the minimum allowed bondlength now!! Wouldn't use this too often

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !Generate the atomic bonding information files used in upcoming learning algorithm
       write(*,*) num_species, num_atoms
       prev_structures = structurecounter("bon")
       !do i=1, num_atoms
       !   call generatebondfiles(1,structures,prev_structures,atomlist,alistrep,num_species,stochio,i)
       !   call generateanglefiles(1,structures,prev_structures,atomlist,alistrep,num_species,stochio,i,bondcutoff)
       !   call generate4files(1,structures,prev_structures,atomlist,alistrep,num_species,stochio,i,bondcutoff)

       !end do
       !do i=1, num_atoms
       !write(*,*) atomlist(structures,i)%position, i
       !end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       !call invar(8,tmpdig,tmpels)
       !do i=1, num_species 
       !   do k=1, num_species 
       !      do x=1, num_atoms*27
       !         if(option_.eq.1) exit
       !         do y=1, num_atoms*27        
       !            if(x.ge.y) cycle
       !            if(alistrep(structures,x)%name.ne.elnames(k)) cycle 
       !            if(alistrep(structures,y)%name.ne.elnames(i)) cycle
       !TURN ME BACK ON
       !if(bondlength(alistrep(structures,x)%position,alistrep(structures,y)%position)&
       !    &.lt.(elrad(1,k,i)*dble(tmpdig(1)/100.0))) then
       !  j=j+1 
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

       
       write(name,'(A11,I0.3,A7)')"pos/POSCAR_",structures+prevpos,"/POSCAR"
       open(newunit = structure_unit, file=name)
       write(structure_unit,*) "Test"
       write(structure_unit,*) 1.0
       do j=1, 3
          write(structure_unit,*) &
               formula(structures)%cell(j,1), &
               formula(structures)%cell(j,2), &
               formula(structures)%cell(j,3)
       end do


       !! Reorganises the POSCAR so all like atoms are next to each other, ...
       !! ... cutting down on space in the input line
       !! Identifies the string of elnames with doubles reduced
       L=0
       allocate(elnames_copy(size(elnames)))
       elnames_copy=""
       elnames_copy(1)=elnames(1)
       k = 1
       elnames_copy = [elnames(1)]
       envelope : do i=2, size(elnames)
          jloop: do j = 1, i - 1, 1
             if(trim(elnames(j)).eq.trim(elnames(i))) cycle envelope
          end do jloop
          k = k + 1
          elnames_copy = [ elnames(i), elnames_copy ]
       end do envelope
       deallocate(elnames)
       elnames = elnames_copy


       !! Swap atoms to the correct positions
       do i = 1, size(elnames)
          do j = 1, num_atoms
             if(atomlist(structures,j)%name.ne.elnames(i)) then
                do k=1, num_atoms 
                   if(atomlist(structures,k)%name.eq.elnames(i)) then 
                      copy_list(1)=atomlist(structures,j)
                      copy_list(2)=atomlist(structures,k)

                      atomlist(structures,k)=copy_list(1) 
                      atomlist(structures,j)=copy_list(2)
                   end if
                end do
             end if
          end do
       end do


       !!-----------------------------------------------------------------------
       !! write files for the structure
       !!-----------------------------------------------------------------------
       call poswrite(structure_unit,formula(structures)%cell,atomlist,num_atoms, structures, num_structures, prevpos)
       write(tmp,'(A11,I0.3)')"pos/POSCAR_",prevpos+structures
       call Incarwrite(adjustl(tmp),500, 20*num_atoms)
       !write(*,*) structures, prevpos, "!!!!!!!!!!!!!!"
       call Jobwrite(tmp,3,3,3)

       call generate_potcar(tmp, elnames)
       close(structure_unit)

       structures=structures+1

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


  subroutine addhost(num_species,structures,formula,location,atomlist,stochio,elnames,&
       addedelements,name,num_structures)
    type(unitcell), dimension(:), allocatable :: formula
    character(1024) :: tmp, name, location
    character(3), dimension(:), allocatable :: elnames,tmpelnames
    real(real12) :: cellmultiplier,meanvol,tmpdble
    integer :: d,l,k,j,i,num_structures, ecount,addedelements, num_species, leng,structures, prev_structures
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
       read(61,*) tmparray
       do d=1, num_structures       
          formula(d)%cell(i,:) = tmparray(:)
       end do
    end do

    !write(*,*) formula(1)%cell
    !! Probably incorrect to do this step here
    !meanvol=get_volume(formula(structures)%cell)
!!!!!!!!!!!!!
    ecount=0
    allocate(tmpelnames(100))
    allocate(tmplist(1,1000))


    !!! READS SPECIES LINE
    read(61,'(A)') tmp
    addedelements=0
    do i=1,len(tmp)-1
       if(i.eq.1) then 
          if((scan(tmp(i:i+1)," ").eq.0).or.&
               &((scan(tmp(i:i+1)," ").eq.2))) then
             addedelements=addedelements+1
             num_species=num_species+1
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
             num_species=num_species+1
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
    !!! APPENDS THE NEW SPECIES TO THE END OF THE OLD SPECIES
    do i=addedelements+1,num_species
       j=j+1
       tmpelnames(i)=elnames(j)
    end do
    deallocate(elnames)
    allocate(elnames(num_species)) 
    do i=1, num_species
       elnames(i)=tmpelnames(i)
    end do
    deallocate(tmpelnames)

    !!! READS STOICHIOMETRY LINE
    allocate(tmpstochio(addedelements))
    read(61,*) tmpstochio
    allocate(tmpstochiotot(num_species))

    do i=1,addedelements
       tmpstochiotot(i)=tmpstochio(i)
    end do
    do i=1, num_species-addedelements
       tmpstochiotot(i+addedelements)=stochio(i) 
    end do
    deallocate(stochio) 
    deallocate(tmpstochio) 
    allocate(stochio(num_species))
    do i=1, num_species
       stochio(i)=tmpstochiotot(i)
    end do


    !!! not sure what tmplist is
    !!! clearly, this is appending elnames to each atom in the list
    k=0
    do i=1, num_species
       do j=1, stochio(i)
          k=k+1
          tmplist(1,k)%name=elnames(i)
       end do
    end do
    write(*,*) elnames
    l=0

    !!! READ CARTESIAN/DIRECT LINE
    read(61,*) tmp   

    !!! READ ATOM POSITIONS
    do i=1, addedelements
       do j=1, stochio(i)
          l=l+1
          read(61,*) tmplist(1,L)%position
          tmplist(1,L)%position=matmul(tmplist(1,L)%position,formula(1)%cell)
       end do
    end do

    !!! APPENDS THE NEW ATOMS TO THE OLD ATOMS
    do d=1, num_structures    
       l=0
       do i=1, num_species
          do j=1, stochio(i) 
             l=l+1
             if(i.le.addedelements) atomlist(d,L)%position=tmplist(1,L)%position
             atomlist(d,L)%name=tmplist(1,L)%name
          end do
       end do
    end do


    write(*,*) "Host accepted"
    deallocate(tmplist)

    close(61)

  end subroutine addhost


  subroutine addxyzfile()
    type(unitcell), dimension(:), allocatable :: formula
    character(1024) :: name, buffer
    character(1024), dimension(:), allocatable :: elnames,elnames_tmp, elnames_list
    real(real12) :: cellmultiplier
    real(real12), dimension(:), allocatable :: bondcutoff
    integer :: tmp,q,l,k,j,i,num_structures, ecount, num_species, num_atoms,structures, prev_structures
    type (atom), dimension(:,:), allocatable :: tmplist,atomlist, alistrep
    integer, dimension(:), allocatable :: stochio

    integer :: unit
    !! Wipes the randomly generated formula
    structures=1
    num_structures=1
    prev_structures=structurecounter("pos")
    call touch("pos")
    call touchposdir(structures,prev_structures)
    allocate(formula(num_structures))

    write(6,*) "Please enter the filename you wish to add to the database"
    read(*, *) name
    !name="asio2_06.xyz"
    open(50, file=name)
    read(50, *) buffer
    read(50, *) cellmultiplier
    formula(1)%cell=0
    read(50, *) formula(1)%cell(1,1),formula(1)%cell(2,2),formula(1)%cell(3,3)
    write(name,'(A11,I0.3,A7)')"pos/POSCAR_",structures+prev_structures,"/POSCAR"
    
    open(newunit = unit, file=name,status="new")
    write(unit,*) "Test"
    write(unit,*) 1.0
    do i=1, 3
       write(unit,*)  formula(1)%cell(:,i)
    end do


!!!!WARNING ONLY FOR ONE SYSTEM USE WITH CAUTION

    tmp=3


    allocate(elnames_list(tmp))
    do i=1, tmp 
       read(50,*) elnames_list(i) 
    end do
    do i=1, tmp
       if(i.eq.1) then 
          num_species=1
          allocate(elnames(1))
          elnames(1)=elnames_list(1)
       else 
          loop2 :do j=1, num_species
             if(elnames_list(i).eq.elnames(j)) exit loop2 
             if(j.eq.num_species) then 
                allocate(elnames_tmp(num_species))
                elnames_tmp=elnames
                deallocate(elnames)
                num_species=num_species+1
                allocate(elnames(num_species))
                do k=1, num_species-1
                   elnames(k)=elnames_tmp(k)
                end do
                elnames(num_species)=elnames_list(i)
             end if
          end do loop2
       end if
    end do
    allocate(stochio(num_species))
    stochio=0
    do i=1, num_species 
       do j=1, tmp 
          if(elnames(i).eq.elnames_list(j)) stochio(i)=stochio(i)+1
       end do
    end do
    num_atoms=0
    do i=1, num_species 
       num_atoms=num_atoms+stochio(i) 
    end do

    allocate(atomlist(1,num_atoms))
    allocate(bondcutoff(num_species))
    write(*,*) "MANUALLY CHANGE"
    stop
    k=0
    do i=1, num_species
       rewind(50)
       read(50,*) 
       read(50,*) tmp
       read(50,*)
       do j=1, tmp
          if(elnames(i).eq.elnames_list(j)) then
             k=k+1
             read(50,*) buffer, atomlist(1,k)%position
             atomlist(1,k)%name=elnames(i)
          else 
             read(50,*)
          end if
       end do
    end do
    close(50)


    call poswrite(unit, formula(1)%cell,atomlist,k,1,1,prev_structures)
    allocate(alistrep(1,k*27))
    do i=1, k 
       call atomrepeater(1,atomlist(1,i)%position,atomlist(1,i)%name,alistrep,formula,i,k)
    end do

    prev_structures=structurecounter("bon")


    !do i=1, k
    !   call generatebondfiles(1,1,prev_structures,atomlist,alistrep,num_species,stochio,i,elnames)
    !end do
    prev_structures=structurecounter("bad")
    bondcutoff=2.0
    !do i=1, k
    !   call generateanglefiles(1,1,prev_structures,atomlist,alistrep,num_species,stochio,i,bondcutoff)
    !end do

    close(unit)

  end subroutine addxyzfile



end module gen
