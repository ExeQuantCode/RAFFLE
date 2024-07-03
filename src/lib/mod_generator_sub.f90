submodule(generator) generator_submodule
  use constants, only: verbose
  use misc_raffle, only: shuffle
  use rw_geom, only: geom_read, geom_write, clone_bas
  use edit_geom, only: bas_merge
  use add_atom, only: add_atom_void, add_atom_pseudo, add_atom_scan, &
       get_viable_gridpoints, update_viable_gridpoints

#ifdef ENABLE_ATHENA
  use read_structures, only: get_graph_from_basis
  use machine_learning, only: network_predict_graph
  use athena, only: graph_type
#endif

  implicit none



contains


  module function init_raffle_generator( &
       lattice_host, basis_host, width, sigma, cutoff_min, cutoff_max ) &
       result(generator)
    !! Initialise an instance of the raffle generator.
    !! Set up run-independent parameters.
    implicit none
    ! Arguments
    real(real12), dimension(3,3), intent(in) :: lattice_host
    !! Lattice vectors of the host structure.
    type(bas_type), intent(in) :: basis_host
    !! Basis of the host structure.
    real(real12), dimension(3), intent(in), optional :: width
    !! Width of the gaussians used in the 2-, 3-, and 4-body 
    !! distribution functions.
    real(real12), dimension(3), intent(in), optional :: sigma
    !! Width of the gaussians used in the 2-, 3-, and 4-body
    !! distribution functions.
    real(real12), dimension(3), intent(in), optional :: cutoff_min
    !! Minimum cutoff for the 2-, 3-, and 4-body distribution functions.
    real(real12), dimension(3), intent(in), optional :: cutoff_max
    !! Maximum cutoff for the 2-, 3-, and 4-body distribution functions.

    type(raffle_generator_type) :: generator


    generator%lattice_host = lattice_host
    generator%basis_host = basis_host

    if( present(width) ) &
         call generator%distributions%set_width(width)
    if( present(sigma) ) &
         call generator%distributions%set_sigma(sigma)
    if( present(cutoff_min) ) &
         call generator%distributions%set_cutoff_min(cutoff_min)
    if( present(cutoff_max) ) &
         call generator%distributions%set_cutoff_max(cutoff_max)


  end function init_raffle_generator



  module subroutine generate(this, num_structures, &
       stoichiometry, method_probab)
    !! Generate random structures.
    implicit none
    ! Arguments
    class(raffle_generator_type), intent(inout) :: this
    !! Instance of the raffle generator.
    integer, intent(in) :: num_structures
    !! Number of structures to generate.
    type(stoichiometry_type), dimension(:), intent(in) :: stoichiometry
    !! Stoichiometry of the structures to generate.
    real(real12), dimension(:), intent(in), optional :: method_probab
    !! Probability of each placement method.

    type(bas_type) :: basis, basis_store

    integer, dimension(:,:), allocatable :: placement_list, placement_list_shuffled

    integer :: i, j, k
    integer :: istructure
    integer :: unit, info_unit, structure_unit
    integer :: num_insert_atoms, num_insert_species

    logical :: placed, success
    character(1024) :: buffer

    real(real12), dimension(3) :: method_probab_ = [0.33_real12, 0.66_real12, 1.0_real12]

#ifdef ENABLE_ATHENA
    type(graph_type), dimension(1) :: graph
#endif

    if(present(method_probab)) method_probab_ = method_probab


    !!! THINK OF SOME WAY TO HANDLE THE HOST SEPARATELY
    !!! THAT CAN SIGNIFICANTLY REDUCE DATA USAGE
    num_insert_species = size(stoichiometry)
    num_insert_atoms = sum(stoichiometry(:)%num)
    allocate(basis_store%spec(num_insert_species))
    basis_store%spec(:)%name = stoichiometry(:)%element
    basis_store%spec(:)%num = stoichiometry(:)%num
    basis_store%natom = num_insert_atoms
    basis_store%nspec = num_insert_species
    basis_store%sysname = "inserts"

    do i = 1, basis_store%nspec
       allocate(basis_store%spec(i)%atom(basis_store%spec(i)%num,3), source = 0._real12)
    end do
    basis_store = bas_merge(this%basis_host,basis_store)

    allocate(placement_list(num_insert_atoms,2))
    k = 0
    spec_loop1: do i = 1, basis_store%nspec
       success = .false.
       do j = 1, size(stoichiometry)
          if(trim(basis_store%spec(i)%name).eq.trim(stoichiometry(j)%element)) &
               success = .true.
       end do
       if(.not.success) cycle
       if(i.gt.this%basis_host%nspec)then
          do j = 1, basis_store%spec(i)%num
             k = k + 1
             placement_list(k,1) = i
             placement_list(k,2) = j
          end do
       else
          do j = 1, basis_store%spec(i)%num
             if(j.le.this%basis_host%spec(i)%num) cycle
             k = k + 1
             placement_list(k,1) = i
             placement_list(k,2) = j
          end do
       end if
    end do spec_loop1


    !!--------------------------------------------------------------------------
    !! generate the structures
    !!--------------------------------------------------------------------------
    structure_loop: do istructure = 1, num_structures

       basis = this%generate_structure( basis_store, &
            placement_list, method_probab_ )
       
#ifdef ENABLE_ATHENA
       !!-----------------------------------------------------------------------
       !! predict energy using ML
       !!-----------------------------------------------------------------------
       graph(1) = get_graph_from_basis(this%lattice_host, basis)
       write(*,*) "Predicted energy", network_predict_graph(graph(1:1))
#endif

    end do structure_loop
    write(*,*) "Finished generating structures"

  end subroutine generate



  module function generate_structure( &
       this, &
       basis_initial, &
       placement_list, method_probab ) result(basis)
    !! Generate a single random structure.
    implicit none
    ! Arguments
    class(raffle_generator_type), intent(in) :: this
    !! Instance of the raffle generator.
    type(bas_type), intent(in) :: basis_initial
    !! Initial basis to build upon.
    integer, dimension(:,:), intent(in) :: placement_list
    !! List of possible placements.
    real(real12), dimension(3) :: method_probab
    !! Probability of each placement method.
    type(bas_type) :: basis
    !! Generated basis.

    integer :: i, j, iplaced, void_ticker
    integer :: num_insert_atoms
    real(real12) :: rtmp1
    logical :: placed
    integer, dimension(size(placement_list,1),size(placement_list,2)) :: &
         placement_list_shuffled
    real(real12), dimension(3) :: method_probab_
    real(real12), dimension(:,:), allocatable :: viable_gridpoints



    call clone_bas(basis_initial, basis)
    num_insert_atoms = basis%natom - this%basis_host%natom

    placement_list_shuffled = placement_list
    call shuffle(placement_list_shuffled,1) !!! NEED TO SORT OUT RANDOM SEED

    viable_gridpoints = get_viable_gridpoints( this%bins, &
         this%lattice_host, basis, &
         [ this%distributions%bond_info(:)%radius_covalent ], &
         placement_list_shuffled )

    method_probab_ = method_probab

    iplaced = 0
    void_ticker = 0
    placement_loop: do while (iplaced.lt.num_insert_atoms)

    !!! CHANGE THESE PLACEMENT SUBROUTINES TO FUNCTIONS THAT OUTPUT THE COORDINATE
    !!! THEN, THIS LOOP ACTUALLY PLACES IT AT THE END
       call random_number(rtmp1)
       if(rtmp1.le.method_probab_(1)) then 
          if(verbose.gt.0) write(*,*) "Add Atom Void"
          call add_atom_void( this%bins, &
                this%lattice_host, basis, &
                placement_list_shuffled(iplaced+1:,:), placed)
       else if(rtmp1.le.method_probab_(2)) then 
          if(verbose.gt.0) write(*,*) "Add Atom Pseudo"
          call add_atom_pseudo( this%bins, &
                this%distributions, &
                this%lattice_host, basis, &
                placement_list_shuffled(iplaced+1:,:), &
                [ this%distributions%bond_info(:)%radius_covalent ], &
                placed )
          if(.not. placed) void_ticker = void_ticker + 1
       else if(rtmp1.le.method_probab_(3)) then 
          if(verbose.gt.0) write(*,*) "Add Atom Scan"
          call add_atom_scan( viable_gridpoints, &
                this%distributions, &
                this%lattice_host, basis, &
                placement_list_shuffled(iplaced+1:,:), &
                [ this%distributions%bond_info(:)%radius_covalent ], &
                placed)
       end if
       if(.not. placed) then
          if(void_ticker.gt.10) &
               call add_atom_void( this%bins, this%lattice_host, basis, &
                                  placement_list_shuffled(iplaced+1:,:), placed)
          void_ticker = 0
          if(.not.placed) cycle placement_loop
       end if
       if(verbose.gt.0)then
          write(*,'(A)',ADVANCE='NO') achar(13)
          write(*,*) "placed", placed
       end if
       iplaced = iplaced + 1
       if(allocated(viable_gridpoints)) &
            call update_viable_gridpoints( viable_gridpoints, &
                             this%lattice_host, basis, &
                             [ placement_list_shuffled(iplaced,:) ], &
                             this%distributions%bond_info( &
                                  ( basis%nspec - &
                                     placement_list_shuffled(iplaced,1)/2 ) * &
                                  ( placement_list_shuffled(iplaced,1) - 1 ) + &
                                  placement_list_shuffled(iplaced,1) &
                             )%radius_covalent )
       if(.not.allocated(viable_gridpoints).and. &
            abs( method_probab_(3) - method_probab_(2) ) .gt. 1.E-3) then
          write(*,*) "WARNING: No more viable gridpoints"
          write(*,*) "Suppressing SCAN method"
          method_probab_ = method_probab_ / method_probab_(2)
          method_probab_(3) = method_probab_(2)
       end if

    end do placement_loop

  end function generate_structure

end submodule generator_submodule