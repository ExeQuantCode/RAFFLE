module gen
  use constants, only: real12, pi
  use misc, only: touch, shuffle
  use misc_linalg, only: get_spheres_overlap
  use rw_geom, only: bas_type, geom_read, geom_write, clone_bas
  use edit_geom, only: bas_merge
  ! use isolated, only: generate_isolated_calculations
  use vasp_file_handler, only: incarwrite, kpoints_write, generate_potcar
  use inputs, only: vdW, volvar, bins, filename_host
  use add_atom, only: add_atom_void, add_atom_pseudo, add_atom_scan
  use read_chem, only: get_element_radius
  use evolver, only: gvector_container_type
  implicit none


  private

  public :: generation

  !!! MOVE HOST STRUCTURE TO A BASIS THAT IS AN OPTIONAL ARGUMENT FOR THE ...
  !!! ... GENERATION PROCEDURE
  !!! that way, remove it from the inputs use
  !!! move bins, vdw, volvar


contains 

!!!#############################################################################
!!! 
!!!#############################################################################
  subroutine generation(gvector_container, num_structures, task, &
       element_list, stoichiometry_list, method_probab, output_dir)
    implicit none
    integer, intent(inout) :: num_structures !! SHOULD NOT EVEN BE AN ARGUEMNT
    !! MAKE AN INPUT ARGUMENT THAT IS MAX_NUM_STRUCTURES
    integer, intent(in) :: task
    type(gvector_container_type), intent(in) :: gvector_container
    character(len=1024), intent(in), optional :: output_dir
    integer, dimension(:), allocatable, intent(inout) :: stoichiometry_list
    character(3), dimension(:), allocatable, intent(inout) :: element_list
    real(real12), dimension(3), intent(in), optional :: method_probab

    type(bas_type) :: basis_host
    real(real12), dimension(3,3) :: lattice_host

    real(real12), dimension(3,3) :: lattice
    type(bas_type) :: basis, basis_store

    integer, dimension(:,:), allocatable :: placement_list, placement_list_shuffled

    integer :: i, j, k
    integer :: istructure
    integer :: unit, info_unit, structure_unit
    integer :: task_
    integer :: num_species, num_atoms, num_insert_atoms
    integer :: num_insert_species

    real(real12) :: rtmp1
    real(real12) :: posneg, meanvol, connectivity, volmin, total_volume
    real(real12) :: normalisation_a
    logical :: placed
    character(1024) :: buffer, output_dir_ = "iteration1"

    real(real12), dimension(3) :: method_probab_ = [0.33_real12, 0.66_real12, 1.0_real12]
    real(real12), dimension(:,:,:), allocatable :: radius_arr


    write(*,*) "HERE0"
    task_=task
    if(present(method_probab)) method_probab_ = method_probab
    write(*,*) "HERE1"


    !!! THINK OF SOME WAY TO HANDLE THE HOST SEPARATELY
    !!! THAT CAN SIGNIFICANTLY REDUCE DATA USAGE
    num_insert_species = size(element_list)
    num_insert_atoms = sum(stoichiometry_list)
    allocate(basis_store%spec(num_insert_species))
    basis_store%spec(:)%name = element_list
    basis_store%spec(:)%num = stoichiometry_list
    basis_store%natom = num_insert_atoms
    basis_store%nspec = num_insert_species
    basis_store%sysname = "inserts"
    write(*,*) "HERE2"
    allocate(placement_list(num_insert_atoms,2))
    k = 0
    write(*,*) "HERE3"
    do i = 1, basis_store%nspec
       allocate(basis_store%spec(i)%atom(basis_store%spec(i)%num,3), source = 0._real12)
       do j = 1, basis_store%spec(i)%num
          k = k + 1
          placement_list(k,1) = i
          placement_list(k,2) = j
       end do
    end do
    select case(task_)
    case(2)
       open(newunit = unit, file = trim(adjustl(filename_host)))
       call geom_read(unit,lattice_host, basis_host)
       close(unit)
       basis_store = bas_merge(basis_host,basis_store)
    end select
    num_species        = basis_store%nspec
    element_list       = basis_store%spec(:)%name
    stoichiometry_list = basis_store%spec(:)%num


    !! calls the function structurecounter, which provides information about the number of currently existing 
    !! structures in the directory
    radius_arr = get_element_radius(element_list)

    ! !!--------------------------------------------------------------------------
    ! !! set up isolated element calculations
    ! !!--------------------------------------------------------------------------
    ! call generate_isolated_calculations(element_list)


    !!--------------------------------------------------------------------------
    !! create the output directory
    !!--------------------------------------------------------------------------
    if(present(output_dir)) output_dir_ = output_dir
    call touch(output_dir_)


    !!--------------------------------------------------------------------------
    !! calculate the expected cell volume
    !!--------------------------------------------------------------------------
    !! Meanvol takes the atomic radius and calculates a guestimate for ...
    !! ... the total rough cell volume.

    !! calculate the normalisation factor
    normalisation_a = sum(stoichiometry_list**2)
    do i = 1, num_species 
       do j = i + 1, num_species, 1
          normalisation_a = normalisation_a + ( stoichiometry_list(i) + stoichiometry_list(j) )**2
       end do
    end do
    normalisation_a = ( basis_store%natom ** 2._real12 ) / normalisation_a

    !!! YOU CAN GET EXACT VOLUME TAKEN UP BY HOST STRUCTURE
    !!! get volume of all atoms associated with chem.in radius
    !!! then subtract all the overlaps
    !!! this is calculated by checking for nearest neighbours
    !!! apply a packing fraction (0.74 for FCC, 0.68 for BCC, 0.52 for SC)
    !!! then work out the estimated volume needed for the inserts

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
                    (stoichiometry_list(j)*radius_arr(3,j,i)) &
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

    posneg = 1
    call random_number(rtmp1)
    if(rtmp1.lt.0.5_real12) posneg = -posneg

    call random_number(rtmp1)
    write(*,*) meanvol
    meanvol = meanvol + ((volvar/100.0)*rtmp1*posneg*meanvol)
    write(*,*) "The allocated volume is", meanvol


    !!--------------------------------------------------------------------------
    !! generate the structures
    !!--------------------------------------------------------------------------
    BIGLOOP: do istructure = 1, num_structures

       basis = generate_structure(basis_store, basis_host, lattice, &
            placement_list, radius_arr, method_probab_)

       !!-----------------------------------------------------------------------
       !! write generated POSCAR
       !!-----------------------------------------------------------------------
       write(buffer,'(A,"/struc",I0.3)') trim(output_dir_),istructure
       call touch(buffer)
       open(newunit = structure_unit, file=trim(buffer)//"/POSCAR")
       call geom_write(structure_unit, lattice_host, basis)
       close(structure_unit)
    
       !!-----------------------------------------------------------------------
       !! write additional VASP files
       !!-----------------------------------------------------------------------
       call Incarwrite(adjustl(buffer),500, 20*num_atoms)
       call kpoints_write(buffer,3,3,3)
       call generate_potcar(buffer, element_list)

    end do BIGLOOP
    write(*,*) "Finished generating structures"

  end subroutine generation
!!!#############################################################################


!!!#############################################################################
!!! 
!!!#############################################################################
  function generate_structure(basis_initial, basis_host, lattice, &
       placement_list, radius_arr, method_probab) result(basis)
    implicit none
    type(bas_type), intent(in) :: basis_initial, basis_host
    real(real12), dimension(3,3), intent(in) :: lattice
    integer, dimension(:,:), intent(in) :: placement_list
    real(real12), dimension(3) :: method_probab
    type(bas_type) :: basis

    integer :: i, j, iplaced
    integer :: num_insert_atoms
    real(real12) :: rtmp1
    logical :: placed
    integer, dimension(size(placement_list,1),size(placement_list,2)) :: &
         placement_list_shuffled
    real(real12), dimension(:,:,:) :: radius_arr


    call clone_bas(basis_initial, basis)
    num_insert_atoms = basis%natom - basis_host%natom

    placement_list_shuffled = placement_list
    call shuffle(placement_list_shuffled,1) !!! NEED TO SORT OUT RANDOM SEED

    iplaced = 0
    placement_loop: do while (iplaced.lt.num_insert_atoms)

       call random_number(rtmp1)
       if(rtmp1.le.method_probab(1)) then 
          write(*,*) "ADD ATOM VOID"
          call add_atom_void( bins, &
                lattice, basis, &
                placement_list_shuffled(iplaced+1:,:), placed)
       else if(rtmp1.le.method_probab(2)) then 
          write(*,*) "Add Atom Pseudo"
          call add_atom_pseudo( bins, &
                lattice, basis, &
                placement_list_shuffled(iplaced+1:,:), radius_arr, placed)
       else if(rtmp1.le.method_probab(3)) then 
          write(*,*) "Add Atom Scan"
          call add_atom_scan( bins, &
                lattice, basis, &
                placement_list_shuffled(iplaced+1:,:), radius_arr, placed)
       end if
       write(*,*) "placed", placed
       write(*,*) "iplaced", iplaced
       write(*,*) "method", rtmp1, method_probab
       if(.not. placed) cycle placement_loop
       iplaced = iplaced + 1

    end do placement_loop

  end function generate_structure
!!!#############################################################################

end module gen