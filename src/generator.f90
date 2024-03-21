module gen
  use constants, only: real12, pi
  use misc, only: touch, shuffle
  use misc_linalg, only: get_spheres_overlap
  use rw_geom, only: bas_type, geom_read, geom_write, clone_bas
  use edit_geom, only: bas_merge
  use isolated, only: generate_isolated_calculations
  use vasp_file_handler, only: &
       touchposdir, &
       incarwrite, kpoints_write, generate_potcar, &
       get_num_atoms_from_poscar
  use inputs, only: &
       vdW, volvar, minbond, maxbond,&
       sigma_bondlength, bins, vps_ratio, filename_host
  use help
  use atomtype
  use add_atom
  use read_chem, only: get_element_radius
  implicit none


  private

  public :: generation

contains 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine generation(alistrep, num_structures, task, &
       element_list, stoichiometry_list)
    implicit none
    integer, intent(inout) :: num_structures !! SHOULD NOT EVEN BE AN ARGUEMNT
    !! MAKE AN INPUT ARGUMENT THAT IS MAX_NUM_STRUCTURES
    integer, intent(in) :: task
    !! element_list is a 1D array containing the symbol for each of the atoms (length=num_species).
    integer, dimension(:), allocatable, intent(inout) :: stoichiometry_list
    character(3), dimension(:), allocatable, intent(inout) :: element_list

    type(bas_type) :: host_basis
    real(real12), dimension(3,3) :: host_lattice

    real(real12), dimension(3,3) :: lattice
    type(bas_type), allocatable :: basis, basis_store

    integer, dimension(:,:), allocatable :: placement_list, placement_list_shuffled

    integer :: unit, info_unit, structure_unit
    integer :: task_
    integer :: num_species, num_atoms
    integer :: num_insert_species
    integer :: num_host_species, num_host_atoms

    real(real12) :: rtmp1, rtmp2


    real(real12) :: posneg, meanvol, q, normvol, calc,sigma1

    integer :: l, i, j, k, x, y, z, m,p
    integer :: structures, prev_structures, prevpos
    integer :: num_VOID

    real(real12), dimension(3) :: angle, tmpvector, bestlocation
    real(real12) ::  bondcutoff, connectivity, volmin, distribution, sigma2
    real(real12) :: angle_distribution, bond_distribution, normalisation_a, prob_void, prob_scan, prob_pseudo
    type(atom), dimension(2) :: copy_list

    real(real12), dimension(:,:,:), allocatable :: radius_arr
    !! Atomlist contains the information about the positions of all atoms in ALL structures. This may cause issues 
    !! with memory when large numbers of structures are used, may consider breaking down into seperate iterations 
    !! (e.g paralyse)
    !! atomlist has dimensions(structure, atom)
    !! there is no distinguishing factor for species
    !! alistrep contains positions of all atoms repeated in adjacent unit cells. Same point as above. 
    type (atom), dimension(:,:), allocatable :: atomlist, alistrep
    character(1024) :: name, tmp, command,location
    character(3), dimension(:), allocatable :: element_list_copy
    integer, dimension(:), allocatable :: stoichiometry_list_copy
    logical :: placed


    !! The info file doesn't contain much of use yet. Could build in if relevant 
    open(newunit=info_unit, file="Info")
    task_=task

    !!! THINK OF SOME WAY TO HANDLE THE HOST SEPARATELY
    !!! THAT CAN SIGNIFICANTLY REDUCE DATA USAGE
    num_insert_species = size(element_list)
    num_atoms = sum(stoichiometry_list)
    allocate(basis_store%spec(num_insert_species))
    basis_store%spec(:)%name = element_list
    basis_store%spec(:)%num = stoichiometry_list
    basis_store%natom = num_atoms
    allocate(placement_list(num_atoms,2))
    k = 0
    do i = 1, basis_store%nspec
       do j = 1, basis_store%spec(i)%num
          k = k + 1
          placement_list(k,1) = i
          placement_list(k,2) = j
       end do
    end do
    select case(task_)
    case(2)
       num_host_atoms = get_num_atoms_from_poscar(filename_host)
       allocate(atomlist(num_structures,num_atoms))
       allocate(alistrep(num_structures,num_atoms*27))
       open(newunit = unit, file = trim(adjustl(filename_host)))
       call geom_read(unit,host_lattice, host_basis)
       close(unit)
       basis_store = bas_merge(host_basis,basis_store)
       num_host_species = host_basis%nspec
    case default
       allocate(atomlist(num_structures,num_atoms))
       allocate(alistrep(num_structures,num_atoms*27))
       prevpos=structurecounter("pos")
       num_host_species=0
    end select
    num_species        = basis_store%nspec
    element_list       = basis_store%spec(:)%name
    stoichiometry_list = basis_store%spec(:)%num


    !! calls the function structurecounter, which provides information about the number of currently existing 
    !! structures in the directory
    prev_structures = structurecounter("pos")
    radius_arr = get_element_radius(element_list)

    !! task_=1 is a special task allowing a new poscar to be added in at user specification


    !! Could implement num_structures>1 in the future for large imports
    if(task_.eq.1) num_structures=1
    structures = 1


    !!--------------------------------------------------------------------------
    !! set up isolated element calculations
    !!--------------------------------------------------------------------------
    call generate_isolated_calculations(element_list)


    !!--------------------------------------------------------------------------
    !! create the 'pos' and 'don' directories
    !!--------------------------------------------------------------------------
    call touch("pos")
    call touch("don")


!!! Create the directories for all of the POSCARS to be placed into !!!
    !! BIGLOOP is the parent loop for all procesess, generating one structure for each full completed iteration
    BIGLOOP: do while(structures.le.num_structures)
!!!#############################################################################
!!!#############################################################################
!!! THIS SHOULD BE ITS OWN PROCEDURE THEN CALLED IN THE MAIN LOOP
!!!#############################################################################

       call clone_bas(basis_store, basis)

       !! Total number of atoms 
       num_atoms = sum(stoichiometry_list)

       !! WHY DEALLOCATE EVERY TIME?
       !! HANDLES THE CASE WHERE THE HOST IS TO BE USED
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

       
       !! calculate the normalisation factor
       normalisation_a = sum(stoichiometry_list**2)
       do i=1, num_species 
          do j = i + 1, num_species, 1
             normalisation_a = normalisation_a + ( stoichiometry_list(i) + stoichiometry_list(j) )**2
          end do
       end do
       normalisation_a = ( num_atoms ** 2._real12 ) / normalisation_a


       !! calculate the minimum volume
       volmin = 0._real12
       connectivity = vdW / 100._real12
       do i = 1, num_species
          volmin = volmin + stoichiometry_list(i) * (4._real12/3._real12) * pi * &
               ( radius_arr(2,i,i) ** 3._real12 )
          do j = i, num_species, 1
             total_volume = ( 4._real12 / 3._real12 ) * pi * &
                            ( radius_arr(2,i,i) ** 3._real12 + &
                              radius_arr(2,j,j) ** 3._real12 ) - &
                            get_spheres_overlap(&
                                 radius_arr(2,i,i),&
                                 radius_arr(2,j,j),&
                                 radius_arr(1,i,j))
             rtmp1 = connectivity * normalisation_a * &
                  min( &
                       (stoichiometry_list(i)*radius_arr(3,i,j)), &
                       (stoichiometry_list(j)*radius_arr(3,j,i))
                  ) * total_volume
             if( i .eq. j )then ! I think this is to reduce significance of same-species bonding
                rtmp1 = 0.5_real12 * rtmp1
             else ! Ned introduced this to account for loop now ignoring half the triangle, TEST!!!
                rtmp1 = 2._real12 * rtmp1
             end if
             volmin = volmin - rtmp1
          end do
          !  meanvol=meanvol+stoichiometry_list(i)*((connectivity*radius_arr(1,i,i))+((1.0-connectivity)*radius_arr(2,i,i)))**3*(4/3)*pi 
       end do
       meanvol = volmin

       !!-----------------------------------------------------------------------
       !! sets pseudorandom volume to larger or smaller than atomic volume summation
       !!-----------------------------------------------------------------------
       posneg = 1
       call random_number(rtmp1)
       if(rtmp1.lt.0.5_real12) posneg = -posneg

       !!-----------------------------------------------------------------------
       !!Adds or subtracts a small quantity from the calculated volume
       !!-----------------------------------------------------------------------
       call random_number(rtmp1)
       write(*,*) meanvol
       meanvol = meanvol + ((volvar/100.0)*rtmp1*posneg*meanvol)
       write(*,*) "The allocated volume is", meanvol




!!!!!!!!!!!!! This section places the atoms via distribution !!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !i controls which atom is being placed. it will be reduced to 0 if the first atom is not to be randomly placed. Bad for small cells
       !random_seed determines if for a host calculation the first atom should be randomly placed. 1 results in random seeding

       !!-----------------------------------------------------------------------
       !! place the first atom in the cell
       !! if a host exists, this step is skipped
       !!-----------------------------------------------------------------------
       placement_list_shuffled = placement_list
       call shuffle(placement_list_shuffled,1) !!! NEED TO SORT OUT RANDOM SEED
       select case(task_)
       case(2)   
          i = 0
       case default
          do j = 1, 3
             call random_number(rtmp2)
             basis % &
                  spec(placement_list_shuffled(1,1)) % &
                  atom(placement_list_shuffled(1,2),j) = rtmp2
          end do
          basis % &
                  spec(placement_list_shuffled(1,1)) % &
                  atom(placement_list_shuffled(1,2),:) = &
                  matmul( host_lattice, &
                  basis % &
                  spec(placement_list_shuffled(1,1)) % &
                  atom(placement_list_shuffled(1,2),:) )
          i = 1
       end select

       sigma1 = sigma_bondlength
       ! This next line is intended to auto-tune the resolution of the guassian sampling. Needs testing
       !    sigma2=minval(peakseparation)/(2.0*sqrt(2.0*LOG(4.0)))
       sigma2 = 0.5
       i=i+L
       bondcutoff=2

       ! First pass is bin size, which should be tied to gaussian size of evolved functions

       num_VOID=1
       aloop: do while (i.lt.num_atoms)
          tmpvector = 0._real12
          do j = 1, 3
             call random_number(rtmp1)
             tmpvector(j) = rtmp1
          end do

          call random_number(rtmp1)

          prob_void=real(vps_ratio(1)/(vps_ratio(1)+vps_ratio(2)+vps_ratio(3)),real12)
          prob_pseudo=prob_void+&
               real(vps_ratio(2)/(vps_ratio(1)+vps_ratio(2)+vps_ratio(3)),real12)
          prob_scan=prob_pseudo+&
               real(vps_ratio(3)/(vps_ratio(1)+vps_ratio(2)+vps_ratio(3)),real12)


          write(*,*) prob_void, prob_pseudo, prob_scan 

          i = i + 1
          if(rtmp1.le.prob_void) then 
             write(*,*) "ADD ATOM VOID"

             !!! PROVIDE SOMETHING THAT TELLS IT WHAT ATOMS TO IGNORE IN THE CHECK
             call add_atom_void(bins, &
                  host_lattice, basis, &
                  placement_list_shuffled(i:,:))
             num_VOID=num_VOID+1
          else if(rtmp1.le.prob_pseudo) then 
             write(*,*) num_VOID, "Add Atom Pseudo"
             call add_atom_pseudo (bins, &
                  host_lattice, basis, &
                  placement_list_shuffled(i:,:), radius_arr, placed)
             num_VOID=1
          else if(rtmp1.le.prob_scan) then 
             write(*,*) num_VOID, "Add Atom Scan"
             call add_atom_pseudo (bins, &
                  host_lattice, basis, &
                  placement_list_shuffled(i:,:), radius_arr, placed)
             num_VOID=1
          end if


          distribution=1
          !!Ditch next line for scan
          tmpvector=matmul(host_lattice,tmpvector)

          bestlocation=0

       end do aloop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !Build in some failsafes - changing sigma iteratively or change the probability width
       !can also change the minimum allowed bondlength now!! Wouldn't use this too often

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !Generate the atomic bonding information files used in upcoming learning algorithm
       prev_structures = structurecounter("bon")
       !do i=1, num_atoms
       !   call generatebondfiles(1,structures,prev_structures,atomlist,alistrep,num_species,stoichiometry_list,i)
       !   call generateanglefiles(1,structures,prev_structures,atomlist,alistrep,num_species,stoichiometry_list,i,bondcutoff)
       !   call generate4files(1,structures,prev_structures,atomlist,alistrep,num_species,stoichiometry_list,i,bondcutoff)

       !end do
       !do i=1, num_atoms
       !write(*,*) atomlist(structures,i)%position, i
       !end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       !call invar(8,tmpdig,tmpels)
       !do i=1, num_species 
       !   do k=1, num_species 
       !      do x=1, num_atoms*27
       !         if(task_.eq.1) exit
       !         do y=1, num_atoms*27        
       !            if(x.ge.y) cycle
       !            if(alistrep(structures,x)%name.ne.element_list(k)) cycle 
       !            if(alistrep(structures,y)%name.ne.element_list(i)) cycle
       !TURN ME BACK ON
       !if(bondlength(alistrep(structures,x)%position,alistrep(structures,y)%position)&
       !    &.lt.(radius_arr(1,k,i)*dble(tmpdig(1)/100.0))) then
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

       
       !!-----------------------------------------------------------------------
       !! write generated POSCAR
       !!-----------------------------------------------------------------------
       write(name,'(A11,I0.3,A7)')"pos/POSCAR_",structures+prevpos,"/POSCAR"
       open(newunit = structure_unit, file=name)
       call geom_write(structure_unit, host_lattice, basis)
       close(structure_unit)

       !!-----------------------------------------------------------------------
       !! write additional VASP files
       !!-----------------------------------------------------------------------
       write(tmp,'(A11,I0.3)')"pos/POSCAR_",prevpos+structures
       call Incarwrite(adjustl(tmp),500, 20*num_atoms)
       call kpoints_write(tmp,3,3,3)
       call generate_potcar(tmp, element_list)

       structures = structures + 1

       call sleep(10)
    end do BIGLOOP
  end subroutine generation

end module gen
