module generator
  !! Module for generating random structures from host structures.
  !!
  !! This module contains the raffle generator type, which is used to generate
  !! random structures from a host structure. The raffle generator uses
  !! distribution functions to determine the placement of atoms in the
  !! provided host structure.
  use constants, only: real12
  use misc_raffle, only: strip_null
  use rw_geom, only: bas_type
  use evolver, only: gvector_container_type

  use constants, only: verbose_global => verbose
  use misc_raffle, only: shuffle
  use rw_geom, only: clone_bas
  use edit_geom, only: bas_merge
  use add_atom, only: add_atom_void, add_atom_walk, add_atom_min, &
       get_viable_gridpoints, update_viable_gridpoints

#ifdef ENABLE_ATHENA
  use read_structures, only: get_graph_from_basis
  use machine_learning, only: network_predict_graph
  use athena, only: graph_type
#endif

  implicit none


  private
  public :: raffle_generator_type, stoichiometry_type


  type :: stoichiometry_type
    !! Type for storing the stoichiometry of atoms to be placed in the host
    !! structure.
    character(len=3) :: element
    !! Element symbol.
    integer :: num
    !! Number of atoms.
  end type stoichiometry_type


  type :: raffle_generator_type
    !! Type for instance of raffle generator.
    !!
    !! This type contains the parameters and methods for generating random
    !! structures from a host structure, using the RAFFLE method.
    integer :: num_structures = 0
    !! Number of structures generated. Initialised to zero.
    type(bas_type) :: host
    !! Host structure.
    integer, dimension(3) :: bins
    !! Number of bins to divide the host structure into along each axis.
    type(gvector_container_type) :: distributions
    !! Distribution function container for the 2-, 3-, and 4-body interactions.
    real(real12), dimension(3) :: method_probab
    !! Probability of each placement method.
    type(bas_type), dimension(:), allocatable :: structures
    !! Generated structures.
   contains
    procedure, pass(this) :: set_host
    !! Procedure to set the host structure.
    procedure, pass(this) :: generate
    !! Procedure to generate random structures.
    procedure, pass(this), private :: generate_structure
    !! Procedure to generate a single random structure.
    procedure, pass(this) :: get_structures
    !! Procedure to return the generated structures.
    procedure, pass(this) :: evaluate
    !! Procedure to evaluate the viability of a structure.
  end type raffle_generator_type

  interface raffle_generator_type
    !! Constructor for the raffle generator type.
    module function init_raffle_generator( &
         host, &
         width, sigma, cutoff_min, cutoff_max) result(generator)
      type(bas_type), intent(in), optional :: host
      real(real12), dimension(3), intent(in), optional :: width
      real(real12), dimension(3), intent(in), optional :: sigma
      real(real12), dimension(3), intent(in), optional :: cutoff_min
      real(real12), dimension(3), intent(in), optional :: cutoff_max
      type(raffle_generator_type) :: generator
    end function init_raffle_generator
  end interface raffle_generator_type


contains

!###############################################################################
  module function init_raffle_generator( &
       host, width, sigma, cutoff_min, cutoff_max ) &
       result(generator)
    !! Initialise an instance of the raffle generator.
    !! Set up run-independent parameters.
    implicit none

    ! Arguments
    type(bas_type), intent(in), optional :: host
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

    ! Local variables
    type(raffle_generator_type) :: generator
    !! Instance of the raffle generator.

    ! Handle optional arguments
    ! Set up the host structure
    if(present(host)) call generator%set_host(host)

    ! Set up the distribution function parameters
    if( present(width) ) &
         call generator%distributions%set_width(width)
    if( present(sigma) ) &
         call generator%distributions%set_sigma(sigma)
    if( present(cutoff_min) ) &
         call generator%distributions%set_cutoff_min(cutoff_min)
    if( present(cutoff_max) ) &
         call generator%distributions%set_cutoff_max(cutoff_max)

  end function init_raffle_generator
!###############################################################################


!###############################################################################
  subroutine set_host(this, host)
    !! Set the host structure.
    implicit none
    ! Arguments
    class(raffle_generator_type), intent(inout) :: this
    !! Instance of the raffle generator.
    type(bas_type), intent(in) :: host
    !! Basis of the host structure.

    ! Local variables
    integer :: i


    this%host = host
    do i = 1, this%host%nspec
       this%host%spec(i)%name = strip_null(this%host%spec(i)%name)
    end do
  end subroutine set_host
!###############################################################################


!###############################################################################
  subroutine generate(this, num_structures, &
       stoichiometry, method_probab, seed, verbose)
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
    integer, intent(in), optional :: seed
    !! Seed for the random number generator.
    integer, intent(in), optional :: verbose
    !! Verbosity level.

    ! Local variables
    integer :: i, j, k, istructure, num_structures_old, num_structures_new
    !! Loop counters.
    integer :: num_seed
    !! Number of seeds for the random number generator.
    integer :: num_insert_atoms, num_insert_species
    !! Number of atoms and species to insert (from stoichiometry).
    real(real12) :: total_probab
    !! Total probability of the placement methods.
    logical :: success
    !! Boolean comparison of element symbols.
    integer :: verbose_ = 0
    !! Verbosity level.
    type(bas_type) :: basis_template
    !! Basis of the structure to generate (i.e. allocated species and atoms).
    real(real12), dimension(3) :: &
         method_probab_ = [1.0_real12, 1.0_real12, 1.0_real12]
    !! Default probability of each placement method.

    integer, dimension(:), allocatable :: seed_arr
    !! Array of seeds for the random number generator.
    type(bas_type), dimension(:), allocatable :: tmp_structures
    !! Temporary array of structures (for memory reallocation).

    integer, dimension(:,:), allocatable :: placement_list
    !! List of possible atoms to place in the structure.


#ifdef ENABLE_ATHENA
    type(graph_type), dimension(1) :: graph
    !! Graph for machine learning.
#endif


    if(present(verbose)) verbose_ = verbose

    !---------------------------------------------------------------------------
    ! set the placement method probabilities
    !---------------------------------------------------------------------------
    if(verbose_.gt.0) write(*,*) "Setting method probabilities"
    if(present(method_probab)) method_probab_ = method_probab
    total_probab = real(sum(method_probab_), real12)
    method_probab_ = method_probab_ / total_probab
    method_probab_(2) = method_probab_(2) + method_probab_(1)
    method_probab_(3) = method_probab_(3) + method_probab_(2)
    if(verbose_.gt.0) write(*,*) "Method probabilities (void, walk, min): ", &
         method_probab_


    !---------------------------------------------------------------------------
    ! set the random seed
    !---------------------------------------------------------------------------
    if(present(seed))then
       call random_seed(size=num_seed)
       allocate(seed_arr(num_seed))
       seed_arr = seed 
       call random_seed(put=seed_arr)
    end if


    !---------------------------------------------------------------------------
    ! allocate memory for structures
    !---------------------------------------------------------------------------
    if(verbose_.gt.0) write(*,*) "Allocating memory for structures"
    if(.not.allocated(this%structures))then
       allocate(this%structures(num_structures))
    else
       allocate(tmp_structures(this%num_structures + num_structures))
       tmp_structures(:this%num_structures) = this%structures(:this%num_structures)
       call move_alloc(tmp_structures, this%structures)
    end if


    !---------------------------------------------------------------------------
    ! set up the template basis for generated structures
    !---------------------------------------------------------------------------
    if(verbose_.gt.0) write(*,*) "Setting up basis store"
    num_insert_species = size(stoichiometry)
    num_insert_atoms = sum(stoichiometry(:)%num)
    allocate(basis_template%spec(num_insert_species))
    do i = 1, size(stoichiometry)
       basis_template%spec(i)%name = strip_null(stoichiometry(i)%element)
    end do
    basis_template%spec(:)%num = stoichiometry(:)%num
    basis_template%natom = num_insert_atoms
    basis_template%nspec = num_insert_species
    basis_template%sysname = "inserts"

    do i = 1, basis_template%nspec
       allocate( &
            basis_template%spec(i)%atom(basis_template%spec(i)%num,3), &
            source = 0._real12 &
       )
    end do
    if(.not.allocated(this%host%spec)) stop "Host structure not set"
    basis_template = bas_merge(this%host,basis_template)
    basis_template%lat = this%host%lat


    !---------------------------------------------------------------------------
    ! generate the placement list
    ! placement list is the list of number of atoms of each species that can be
    ! placed in the structure
    ! ... the second dimension is the index of the species and atom in the
    ! ... basis_template
    !---------------------------------------------------------------------------
    if(verbose_.gt.0) write(*,*) "Generating placement list"
    allocate(placement_list(num_insert_atoms,2))
    k = 0
    spec_loop1: do i = 1, basis_template%nspec
       success = .false.
       do j = 1, size(stoichiometry)
          if( &
               trim(basis_template%spec(i)%name) .eq. &
               trim(strip_null(stoichiometry(j)%element))) &
               success = .true.
       end do
       if(.not.success) cycle
       if(i.gt.this%host%nspec)then
          do j = 1, basis_template%spec(i)%num
             k = k + 1
             placement_list(k,1) = i
             placement_list(k,2) = j
          end do
       else
          do j = 1, basis_template%spec(i)%num
             if(j.le.this%host%spec(i)%num) cycle
             k = k + 1
             placement_list(k,1) = i
             placement_list(k,2) = j
          end do
       end if
    end do spec_loop1


    !---------------------------------------------------------------------------
    ! generate the structures
    !---------------------------------------------------------------------------
    if(verbose_.gt.0) write(*,*) "Entering structure generation loop"
    num_structures_old = this%num_structures
    num_structures_new = this%num_structures + num_structures
    structure_loop: do istructure = num_structures_old + 1, num_structures_new
    
       if(verbose_.gt.0) write(*,*) "Generating structure", istructure
       this%structures(istructure) = this%generate_structure( basis_template, &
            placement_list, method_probab_, verbose_ )
       this%num_structures = istructure
       
#ifdef ENABLE_ATHENA
       !------------------------------------------------------------------------
       ! predict energy using ML
       !------------------------------------------------------------------------
       graph(1) = get_graph_from_basis(this%structures(istructure))
       write(*,*) "Predicted energy", network_predict_graph(graph(1:1))
#endif

    end do structure_loop
    write(*,*) "Finished generating structures"

  end subroutine generate
!###############################################################################


!###############################################################################
  module function generate_structure( &
       this, &
       basis_initial, &
       placement_list, method_probab, verbose ) result(basis)
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
    integer, intent(in) :: verbose
    !! Verbosity level.

    ! Local variables
    integer :: i, j, iplaced, void_ticker
    !! Loop counters.
    integer :: num_insert_atoms
    !! Number of atoms to insert.
    real(real12) :: rtmp1
    !! Random number.
    logical :: viable
    !! Boolean for viable placement.
    integer, dimension(size(placement_list,1),size(placement_list,2)) :: &
         placement_list_shuffled
    !! Shuffled placement list.
    real(real12), dimension(3) :: point
    !! Coordinate of the atom to place.
    real(real12), dimension(3) :: method_probab_
    !! Temporary probability of each placement method.
    !! This is used to update the probability of the global minimum method if
    !! no viable gridpoints are found.
    real(real12), dimension(:,:), allocatable :: viable_gridpoints
    !! Viable gridpoints for placing atoms.


    !---------------------------------------------------------------------------
    ! initialise the basis
    !---------------------------------------------------------------------------
    call clone_bas(basis_initial, basis)
    num_insert_atoms = basis%natom - this%host%natom


    !---------------------------------------------------------------------------
    ! shuffle the placement list
    !---------------------------------------------------------------------------
    placement_list_shuffled = placement_list
    call shuffle(placement_list_shuffled,1) !!! NEED TO SORT OUT RANDOM SEED


    !---------------------------------------------------------------------------
    ! check for viable gridpoints
    !---------------------------------------------------------------------------
    method_probab_ = method_probab
    if(abs( method_probab_(3) - method_probab_(2) ) .gt. 1.E-3)then
       viable_gridpoints = get_viable_gridpoints( this%bins, &
            basis, &
            [ this%distributions%bond_info(:)%radius_covalent ], &
            placement_list_shuffled )
    end if


    !---------------------------------------------------------------------------
    ! place the atoms
    !---------------------------------------------------------------------------
    iplaced = 0
    void_ticker = 0
    viable = .false.
    placement_loop: do while (iplaced.lt.num_insert_atoms)
       if(viable)then
          if(allocated(viable_gridpoints)) &
               call update_viable_gridpoints( viable_gridpoints, &
                     basis, &
                     [ placement_list_shuffled(iplaced,:) ], &
                     this%distributions%bond_info( &
                        nint( ( basis%nspec - &
                          placement_list_shuffled(iplaced,1)/2._real12 ) * &
                        ( placement_list_shuffled(iplaced,1) - 1 ) + &
                        placement_list_shuffled(iplaced,1) ) &                       
                     )%radius_covalent )
          if(.not.allocated(viable_gridpoints))then
             if(abs(method_probab_(2)).lt.1.E-3)then
                write(0,*) "ERROR: No viable gridpoints"
                write(0,*) "No placement methods available"
                stop 0
             else if(abs( method_probab_(3) - method_probab_(2) ) .gt. 1.E-3) then
                write(*,*) "WARNING: No more viable gridpoints"
                write(*,*) "Suppressing global minimum method"
                method_probab_ = method_probab_ / method_probab_(2)
                method_probab_(3) = method_probab_(2)
             end if
          end if
       end if
       call random_number(rtmp1)
       if(rtmp1.le.method_probab_(1)) then 
          if(verbose.gt.0) write(*,*) "Add Atom Void"
          point = add_atom_void( this%bins, &
                basis, &
                placement_list_shuffled(iplaced+1:,:), viable )
       else if(rtmp1.le.method_probab_(2)) then 
          if(verbose.gt.0) write(*,*) "Add Atom Walk"
          point = add_atom_walk( &
                this%distributions, &
                basis, &
                placement_list_shuffled(iplaced+1:,:), &
                [ this%distributions%bond_info(:)%radius_covalent ], &
                viable )
          if(.not. viable) void_ticker = void_ticker + 1
       else if(rtmp1.le.method_probab_(3)) then 
          if(verbose.gt.0) write(*,*) "Add Atom Min"
          point = add_atom_min( viable_gridpoints, &
                this%distributions, &
                basis, &
                placement_list_shuffled(iplaced+1:,:), &
                [ this%distributions%bond_info(:)%radius_covalent ], &
                viable )
          if(.not. viable .and. abs(method_probab_(2)).lt.1.E-3)then
             write(0,*) "ERROR: No viable gridpoints"
             write(0,*) "  Min method is the only method, but cannot place another atom"
             write(0,*) "  Species to place now: ", basis%spec(placement_list_shuffled(iplaced+1,1))%name
             stop 0
          end if
       end if
       if(.not. viable) then
          if(void_ticker.gt.10) &
               point = add_atom_void( this%bins, basis, &
                                  placement_list_shuffled(iplaced+1:,:), viable)
          void_ticker = 0
          if(.not.viable) cycle placement_loop
       end if
       iplaced = iplaced + 1
       basis%spec(placement_list_shuffled(iplaced,1))%atom( &
            placement_list_shuffled(iplaced,2),:3) = point(:3)
       if(verbose.gt.0)then
          write(*,'(A)',ADVANCE='NO') achar(13)
          write(*,*) "placed", viable
       end if

    end do placement_loop
    if(allocated(viable_gridpoints)) deallocate(viable_gridpoints)

  end function generate_structure
!###############################################################################


!###############################################################################
  function get_structures(this) result(structures)
    !! Get the generated structures.
    implicit none
    ! Arguments
    class(raffle_generator_type), intent(in) :: this
    !! Instance of the raffle generator.
    type(bas_type), dimension(:), allocatable :: structures
    !! Generated structures.

    structures = this%structures
  end function get_structures
!###############################################################################


!###############################################################################
  function evaluate(this, basis) result(viability)
    !! Evaluate the viability of the generated structures.
    implicit none
    ! Arguments
    class(raffle_generator_type), intent(in) :: this
    !! Instance of the raffle generator.
    type(bas_type), intent(in) :: basis
    !! Basis of the structure to evaluate.
    real(real12) :: viability
    !! Viability of the generated structures.

    viability = 0.0_real12
    stop "Not yet set up"
  end function evaluate
!###############################################################################


!###############################################################################
  subroutine allocate_structures(this, num_structures)
    !! Allocate memory for the generated structures.
    implicit none
    ! Arguments
    class(raffle_generator_type), intent(inout) :: this
    !! Instance of the raffle generator.
    integer, intent(in) :: num_structures
    !! Number of structures to allocate memory for.

    if(allocated(this%structures)) deallocate(this%structures)
    allocate(this%structures(num_structures))
    this%num_structures = num_structures
  end subroutine allocate_structures
!###############################################################################

end module generator