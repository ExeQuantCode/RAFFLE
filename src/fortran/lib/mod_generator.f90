module raffle__generator
  !! Module for generating random structures from host structures.
  !!
  !! This module contains the raffle generator type, which is used to generate
  !! random structures from a host structure. The raffle generator uses
  !! distribution functions to determine the placement of atoms in the
  !! provided host structure.
  use raffle__io_utils, only: stop_program
  use raffle__constants, only: real32
  use raffle__misc_linalg, only: modu
  use raffle__misc, only: strip_null, set, shuffle, sort1D
  use raffle__geom_rw, only: basis_type
  use raffle__geom_extd, only: extended_basis_type
  use raffle__distribs_container, only: distribs_container_type
  use raffle__geom_utils, only: basis_merge
  use raffle__place_methods, only: &
       place_method_void, place_method_rand, &
       place_method_growth, place_method_walk, &
       place_method_min
  use raffle__viability, only: &
       get_gridpoints_and_viability, update_gridpoints_and_viability

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
     type(basis_type) :: host
     !! Host structure.
     integer, dimension(3) :: grid = [0, 0, 0]
     !! Grid to divide the host structure into along each axis.
     real(real32), dimension(3) :: &
          grid_offset = [0.5_real32, 0.5_real32, 0.5_real32]
     !! Offset of the gridpoints.
     real(real32) :: grid_spacing = 0.1_real32
     !! Spacing of the gridpoints.
     type(distribs_container_type) :: distributions
     !! Distribution function container for the 2-, 3-, and 4-body interactions.
     integer :: max_attempts = 10000
     !! Limit for the number of attempts to place an atom.
     real(real32) :: &
          walk_step_size_coarse = 1._real32, &
          walk_step_size_fine = 0.1_real32
     !! Step size for the walk and grow methods.
     real(real32), dimension(5) :: method_probab
     !! Probability of each placement method.
     type(basis_type), dimension(:), allocatable :: structures
     !! Generated structures.
   contains
     procedure, pass(this) :: set_host
     !! Procedure to set the host structure.
     procedure, pass(this) :: set_grid
     !! Procedure to set the grid for the raffle generator.
     procedure, pass(this) :: reset_grid
     !! Procedure to reset the grid for the raffle generator.
     procedure, pass(this) :: generate
     !! Procedure to generate random structures.
     procedure, pass(this), private :: generate_structure
     !! Procedure to generate a single random structure.
     procedure, pass(this) :: get_structures
     !! Procedure to return the generated structures.
     procedure, pass(this) :: set_structures
     !! Procedure to set the array of generated structures.
     procedure, pass(this) :: remove_structure
     !! Procedure to remove a structure from the array of generated structures.
     procedure, pass(this) :: evaluate
     !! Procedure to evaluate the viability of a structure.
  end type raffle_generator_type

  interface raffle_generator_type
     !! Constructor for the raffle generator type.
     module function init_raffle_generator( &
          host, &
          width, sigma, cutoff_min, cutoff_max) result(generator)
       type(basis_type), intent(in), optional :: host
       real(real32), dimension(3), intent(in), optional :: width
       real(real32), dimension(3), intent(in), optional :: sigma
       real(real32), dimension(3), intent(in), optional :: cutoff_min
       real(real32), dimension(3), intent(in), optional :: cutoff_max
       type(raffle_generator_type) :: generator
     end function init_raffle_generator
  end interface raffle_generator_type


contains

!###############################################################################
  module function init_raffle_generator( &
       host, width, sigma, cutoff_min, cutoff_max &
  ) result(generator)
    !! Initialise an instance of the raffle generator.
    !!
    !! Set up run-independent parameters.
    implicit none

    ! Arguments
    type(basis_type), intent(in), optional :: host
    !! Basis of the host structure.
    real(real32), dimension(3), intent(in), optional :: width
    !! Width of the gaussians used in the 2-, 3-, and 4-body 
    !! distribution functions.
    real(real32), dimension(3), intent(in), optional :: sigma
    !! Width of the gaussians used in the 2-, 3-, and 4-body
    !! distribution functions.
    real(real32), dimension(3), intent(in), optional :: cutoff_min
    !! Minimum cutoff for the 2-, 3-, and 4-body distribution functions.
    real(real32), dimension(3), intent(in), optional :: cutoff_max
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
    !!
    !! This procedure sets the host structure for the raffle generator.
    implicit none

    ! Arguments
    class(raffle_generator_type), intent(inout) :: this
    !! Instance of the raffle generator.
    class(basis_type), intent(in) :: host
    !! Basis of the host structure.

    ! Local variables
    integer :: i
    !! Loop index.


    call this%host%copy(host)
    call this%distributions%host_system%set(this%host)

    call this%set_grid()
  end subroutine set_host
!###############################################################################


!###############################################################################
  subroutine set_grid(this, grid, grid_spacing, grid_offset)
    !! Set the grid for the raffle generator.
    !!
    !! This procedure sets the grid for the raffle generator. The grid is used
    !! to divide the host structure into bins along each axis on which
    !! atom placement viability will be evaluated
    implicit none

    ! Arguments
    class(raffle_generator_type), intent(inout) :: this
    !! Instance of the raffle generator.
    integer, dimension(3), intent(in), optional :: grid
    !! Number of bins to divide the host structure into along each axis.
    real(real32), intent(in), optional :: grid_spacing
    !! Spacing of the bins.
    real(real32), dimension(3), intent(in), optional :: grid_offset
    !! Offset of the gridpoints.

    ! Local variables
    integer :: i
    !! Loop index.


    if(present(grid).and.present(grid_spacing)) then
       call this%reset_grid()
       call stop_program("Cannot set grid and grid spacing simultaneously")
       return
    elseif(present(grid_spacing)) then
       this%grid_spacing = grid_spacing
       this%grid = 0
    elseif(present(grid)) then
       this%grid = grid
    end if

    if(present(grid_offset)) this%grid_offset = grid_offset

    if(all(this%grid.eq.0))then
       if(allocated(this%host%spec))then
          do i = 1, 3
             this%grid(i) = nint( modu(this%host%lat(i,:)) / this%grid_spacing )
          end do
       end if
    end if

  end subroutine set_grid
!###############################################################################


!###############################################################################
  subroutine reset_grid(this)
    !! Reset the grid for the raffle generator.
    implicit none

    ! Arguments
    class(raffle_generator_type), intent(inout) :: this
    !! Instance of the raffle generator.

    this%grid = 0
  end subroutine reset_grid
!###############################################################################


!###############################################################################
  subroutine generate(this, num_structures, &
       stoichiometry, method_probab, seed, verbose)
    !! Generate random structures.
    !!
    !! This procedure generates random structures from the contained host
    !! structure and the stoichiometry argument. The number of structures to
    !! generate is specified by the num_structures argument.
    !! The ratio of placement methods to be sampled is defined by method_probab.
    implicit none

    ! Arguments
    class(raffle_generator_type), intent(inout) :: this
    !! Instance of the raffle generator.
    integer, intent(in) :: num_structures
    !! Number of structures to generate.
    type(stoichiometry_type), dimension(:), intent(in) :: stoichiometry
    !! Stoichiometry of the structures to generate.
    real(real32), dimension(5), intent(in), optional :: method_probab
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
    real(real32) :: total_probab
    !! Total probability of the placement methods.
    logical :: success
    !! Boolean comparison of element symbols.
    integer :: verbose_ = 0
    !! Verbosity level.
    type(basis_type) :: basis_template
    !! Basis of the structure to generate (i.e. allocated species and atoms).
    real(real32), dimension(5) :: &
         method_probab_ = &
         [1.0_real32, 0.1_real32, 0.5_real32, 0.5_real32, 1.0_real32]
    !! Default probability of each placement method.

    integer, dimension(:), allocatable :: seed_arr
    !! Array of seeds for the random number generator.
    type(basis_type), dimension(:), allocatable :: tmp_structures
    !! Temporary array of structures (for memory reallocation).

    integer, dimension(:,:), allocatable :: placement_list
    !! List of possible atoms to place in the structure.


    if(present(verbose)) verbose_ = verbose

    !---------------------------------------------------------------------------
    ! set the placement method probabilities
    !---------------------------------------------------------------------------
    if(verbose_.gt.0) write(*,*) "Setting method probabilities"
    if(present(method_probab)) method_probab_ = method_probab
    total_probab = real(sum(method_probab_), real32)
    method_probab_ = method_probab_ / total_probab
    do i = 2, 5, 1
       method_probab_(i) = method_probab_(i) + method_probab_(i-1)
    end do
    if(verbose_.gt.0) write(*,*) &
         "Method probabilities (void, rand, walk, grow, min): ", &
         method_probab_


    !---------------------------------------------------------------------------
    ! set the random seed
    !---------------------------------------------------------------------------
    if(present(seed))then
       call random_seed(size=num_seed)
       allocate(seed_arr(num_seed))
       seed_arr = seed 
       call random_seed(put=seed_arr)
    else
       call random_seed(size=num_seed)
       allocate(seed_arr(num_seed))
       call random_seed(get=seed_arr)
    end if


    !---------------------------------------------------------------------------
    ! allocate memory for structures
    !---------------------------------------------------------------------------
    if(verbose_.gt.0) write(*,*) "Allocating memory for structures"
    if(.not.allocated(this%structures))then
       allocate(this%structures(num_structures))
    else
       allocate(tmp_structures(this%num_structures + num_structures))
       tmp_structures(:this%num_structures) = &
            this%structures(:this%num_structures)
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
            source = 0._real32 &
       )
    end do
    if(.not.allocated(this%host%spec))then
       call stop_program("Host structure not set")
       return
    end if
    basis_template = basis_merge(this%host,basis_template)
    basis_template%lat = this%host%lat


    !---------------------------------------------------------------------------
    ! ensure host element map is set
    !---------------------------------------------------------------------------
    call this%distributions%host_system%set_element_map( &
         this%distributions%element_info &
    )


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
               trim(strip_null(stoichiometry(j)%element)) &
          ) success = .true.
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
       call this%structures(istructure)%copy( basis = &
            this%generate_structure( &
                 basis_template, &
                 placement_list, &
                 method_probab_, &
                 verbose_ &
            ) &
       )
       this%num_structures = istructure

    end do structure_loop
    if(verbose_.gt.0) write(*,*) "Finished generating structures"

  end subroutine generate
!###############################################################################


!###############################################################################
  function generate_structure( &
       this, &
       basis_initial, &
       placement_list, method_probab, verbose &
  ) result(basis)
    !! Generate a single random structure.
    !!
    !! This function generates a single random structure from a host structure
    !! by placing atoms according to the ratio of placement methods.
    !! The input host structure will already have all host and insert species
    !! and atoms allocated. The placement list specifies the atoms in the
    !! host structure to be replaced by insert atoms.
    implicit none

    ! Arguments
    class(raffle_generator_type), intent(in) :: this
    !! Instance of the raffle generator.
    type(basis_type), intent(in) :: basis_initial
    !! Initial basis to build upon.
    integer, dimension(:,:), intent(in) :: placement_list
    !! List of possible placements.
    real(real32), dimension(5) :: method_probab
    !! Probability of each placement method.
    type(extended_basis_type) :: basis
    !! Generated basis.
    integer, intent(in) :: verbose
    !! Verbosity level.

    ! Local variables
    integer :: iplaced, void_ticker
    !! Loop counters.
    integer :: num_insert_atoms
    !! Number of atoms to insert.
    real(real32) :: rtmp1
    !! Random number.
    logical :: viable
    !! Boolean for viable placement.
    integer, dimension(size(placement_list,1),size(placement_list,2)) :: &
         placement_list_shuffled
    !! Shuffled placement list.
    real(real32), dimension(3) :: point
    !! Coordinate of the atom to place.
    real(real32), dimension(5) :: method_probab_
    !! Temporary probability of each placement method.
    !! This is used to update the probability of the global minimum method if
    !! no viable gridpoints are found.
    integer, dimension(:), allocatable :: species_index_list
    !! List of species indices to add.
    real(real32), dimension(:,:), allocatable :: gridpoint_viability
    !! Viable gridpoints for placing atoms.
    character(len=256) :: stop_msg
    !! Error message.


    !---------------------------------------------------------------------------
    ! initialise the basis
    !---------------------------------------------------------------------------
    call basis%copy(basis_initial)
    call basis%create_images( &
         max_bondlength = this%distributions%cutoff_max(1), &
         atom_ignore_list = placement_list &
    )
    num_insert_atoms = basis%natom - this%host%natom


    !---------------------------------------------------------------------------
    ! shuffle the placement list
    !---------------------------------------------------------------------------
    placement_list_shuffled = placement_list
    call shuffle(placement_list_shuffled,1)


    !---------------------------------------------------------------------------
    ! generate species index list to add
    !---------------------------------------------------------------------------
    species_index_list = placement_list_shuffled(:,1)
    call set(species_index_list)


    !---------------------------------------------------------------------------
    ! check for viable gridpoints
    !---------------------------------------------------------------------------
    method_probab_ = method_probab
    if(abs( method_probab_(5) - method_probab_(4) ) .gt. 1.E-3)then
       gridpoint_viability = get_gridpoints_and_viability( &
            this%distributions, &
            this%grid, &
            basis, &
            species_index_list, &
            [ this%distributions%bond_info(:)%radius_covalent ], &
            placement_list_shuffled, &
            this%grid_offset &
       )
    end if


    !---------------------------------------------------------------------------
    ! place the atoms according to the method probabilities
    !---------------------------------------------------------------------------
    iplaced = 0
    void_ticker = 0
    viable = .false.
    placement_loop: do while (iplaced.lt.num_insert_atoms)
       !------------------------------------------------------------------------
       ! check if there are any viable gridpoints remaining
       !------------------------------------------------------------------------
       if(viable)then
          if(allocated(gridpoint_viability)) &
               call update_gridpoints_and_viability( &
                    gridpoint_viability, &
                    this%distributions, &
                    basis, &
                    species_index_list, &
                    [ placement_list_shuffled(iplaced,:) ], &
                    [ this%distributions%bond_info(:)%radius_covalent ], &
                    placement_list_shuffled(iplaced+1:,:) &
               )
          if(.not.allocated(gridpoint_viability))then
             if(abs(method_probab_(4)).lt.1.E-6)then
                call stop_program( &
                     "No viable gridpoints" // &
                     achar(13) // achar(10) // &
                     "No placement methods available" &
                )
                return
             else if( &
                  abs( method_probab_(5) - method_probab_(4) ) .gt. 1.E-6 &
             ) then
                write(*,*) "WARNING: No more viable gridpoints"
                write(*,*) "Suppressing global minimum method"
                method_probab_ = method_probab_ / method_probab_(4)
                method_probab_(5) = method_probab_(4)
             end if
          end if
       end if
       viable = .false.
       !------------------------------------------------------------------------
       ! choose a placement method
       ! call a random number and query the method probabilities
       !------------------------------------------------------------------------
       call random_number(rtmp1)
       if(rtmp1.le.method_probab_(1)) then
          if(verbose.gt.0) write(*,*) "Add Atom Void"
          point = place_method_void( this%grid, &
               this%grid_offset, &
               basis, &
               placement_list_shuffled(iplaced+1:,:), viable &
          )
       else if(rtmp1.le.method_probab_(2)) then
          if(verbose.gt.0) write(*,*) "Add Atom Random"
          point = place_method_rand( &
               this%distributions, &
               basis, &
               placement_list_shuffled(iplaced+1:,:), &
               [ this%distributions%bond_info(:)%radius_covalent ], &
               this%max_attempts, &
               viable &
          )
          if(.not. viable) cycle placement_loop
       else if(rtmp1.le.method_probab_(3)) then
          if(verbose.gt.0) write(*,*) "Add Atom Walk"
          point = place_method_walk( &
               this%distributions, &
               basis, &
               placement_list_shuffled(iplaced+1:,:), &
               [ this%distributions%bond_info(:)%radius_covalent ], &
               this%max_attempts, &
               this%walk_step_size_coarse, this%walk_step_size_fine, &
               viable &
          )
          if(.not. viable) void_ticker = void_ticker + 1
       else if(rtmp1.le.method_probab_(4)) then
          if(iplaced.eq.0)then
             if(verbose.gt.0) write(*,*) "Add Atom Random (growth seed)"
             point = place_method_rand( &
                  this%distributions, &
                  basis, &
                  placement_list_shuffled(iplaced+1:,:), &
                  [ this%distributions%bond_info(:)%radius_covalent ], &
                  this%max_attempts, &
                  viable &
             )
          else
             if(verbose.gt.0) write(*,*) "Add Atom Growth"
             point = place_method_growth( &
                  this%distributions, &
                  basis%spec(placement_list_shuffled(iplaced,1))%atom( &
                       placement_list_shuffled(iplaced,2),:3 &
                  ), &
                  placement_list_shuffled(iplaced,1), &
                  basis, &
                  placement_list_shuffled(iplaced+1:,:), &
                  [ this%distributions%bond_info(:)%radius_covalent ], &
                  this%max_attempts, &
                  this%walk_step_size_coarse, this%walk_step_size_fine, &
                  viable &
             )
          end if
          if(.not. viable) void_ticker = void_ticker + 1
       else if(rtmp1.le.method_probab_(5)) then
          if(verbose.gt.0) write(*,*) "Add Atom Minimum"
          point = place_method_min( gridpoint_viability, &
               placement_list_shuffled(iplaced+1,1), &
               species_index_list, &
               viable &
          )
          if(.not. viable .and. abs(method_probab_(4)).lt.1.E-6)then
             write(stop_msg,*) &
                  "No viable gridpoints" // &
                  achar(13) // achar(10) // &
                  "Min method is the only method, but cannot place another &
                  &atom" // &
                  achar(13) // achar(10) // &
                  "Species to place now: ", &
                  basis%spec(placement_list_shuffled(iplaced+1,1))%name
             call stop_program(stop_msg)
             return
          elseif(.not. viable)then
             deallocate(gridpoint_viability)
             write(*,*) "WARNING: No more viable gridpoints"
             write(*,*) "Suppressing global minimum method"
             method_probab_ = method_probab_ / method_probab_(4)
             method_probab_(5) = method_probab_(4)
          end if
       end if
       !------------------------------------------------------------------------
       ! check if the placement method returned a viable point
       ! if not, cycle the loop
       !------------------------------------------------------------------------
       if(.not. viable) then
          if(void_ticker.gt.10) &
               point = place_method_void( &
                    this%grid, this%grid_offset, basis, &
                    placement_list_shuffled(iplaced+1:,:), viable &
               )
          void_ticker = 0
          if(.not.viable) cycle placement_loop
       end if
       !------------------------------------------------------------------------
       ! place the atom and update the image atoms in the basis
       !------------------------------------------------------------------------
       iplaced = iplaced + 1
       basis%spec(placement_list_shuffled(iplaced,1))%atom( &
            placement_list_shuffled(iplaced,2),:3) = point(:3)
       call basis%update_images( &
            max_bondlength = this%distributions%cutoff_max(1), &
            is = placement_list_shuffled(iplaced,1), &
            ia = placement_list_shuffled(iplaced,2) &
       )
       if(verbose.gt.0)then
          write(*,'(A)',ADVANCE='NO') achar(13)
          write(*,*) "placed", viable
       end if

    end do placement_loop
    if(allocated(gridpoint_viability)) deallocate(gridpoint_viability)

  end function generate_structure
!###############################################################################


!###############################################################################
  function get_structures(this) result(structures)
    !! Get the generated structures.
    implicit none
    ! Arguments
    class(raffle_generator_type), intent(in) :: this
    !! Instance of the raffle generator.
    type(basis_type), dimension(:), allocatable :: structures
    !! Generated structures.

    structures = this%structures
  end function get_structures
!###############################################################################


!###############################################################################
  subroutine set_structures(this, structures)
    !! Set the generated structures.
    !!
    !! This procedure overwrites the array of generated structures with the
    !! input array.
    !! This can be useful for removing structures that are not viable from the
    !! array.
    implicit none
    ! Arguments
    class(raffle_generator_type), intent(inout) :: this
    !! Instance of the raffle generator.
    type(basis_type), dimension(..), allocatable, intent(in) :: structures
    !! Array of structures to set.

    select rank(structures)
    rank(0)
       this%structures = [ structures ]
    rank(1)
       this%structures = structures
    rank default
       call stop_program("Invalid rank for structures")
    end select
    this%num_structures = size(this%structures)
  end subroutine set_structures
!###############################################################################


!###############################################################################
  subroutine remove_structure(this, index)
    !! Remove structures from the generated structures.
    !!
    !! This procedure removes structures from the array of generated structures
    !! at the specified indices.
    implicit none
    ! Arguments
    class(raffle_generator_type), intent(inout) :: this
    !! Instance of the raffle generator.
    integer, dimension(..), intent(in) :: index
    !! Indices of the structures to remove.

    ! Local variables
    integer :: i
    !! Loop index.
    integer, dimension(:), allocatable :: index_
    !! Indices of the structures to keep.

    select rank(index)
    rank(0)
       index_ = [ index ]
    rank(1)
       index_ = index
    rank default
       call stop_program("Invalid rank for index")
    end select

    if(any(index_.lt.1) .or. any(index_.gt.this%num_structures))then
       call stop_program("Invalid index")
       return
    end if

    call sort1D(index_, reverse=.true.)

    do i = 1, size(index_)
       this%structures = [ &
            this%structures(:index_(i)-1:1), &
            this%structures(index_(i)+1:this%num_structures:1) &
       ]
       this%num_structures = this%num_structures - 1
    end do

   end subroutine remove_structure
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


!###############################################################################
  function evaluate(this, basis) result(viability)
    !! Evaluate the viability of the generated structures.
    use raffle__evaluator, only: evaluate_point
    implicit none
    ! Arguments
    class(raffle_generator_type), intent(in) :: this
    !! Instance of the raffle generator.
    type(basis_type), intent(in) :: basis
    !! Basis of the structure to evaluate.
    real(real32) :: viability
    !! Viability of the generated structures.

    ! Local variables
    integer :: is, ia, species
    !! Loop indices.
    integer, dimension(1,2) :: atom_ignore
    !! Atom to ignore.
    type(extended_basis_type) :: basis_extd
    !! Extended basis for the structure to evaluate.


    call basis_extd%copy(basis)
    call basis_extd%create_images( &
         max_bondlength = this%distributions%cutoff_max(1) &
    )
    viability = 0.0_real32
    do is = 1, basis%nspec
       species = this%distributions%get_element_index( basis%spec(is)%name )
       if(species.eq.0)then
          call stop_program("Species not found in distribution functions")
          return
       end if
       do ia = 1, basis%spec(is)%num
          atom_ignore(1,:) = [is,ia]
          viability = viability + &
               evaluate_point( this%distributions, &
                    basis%spec(is)%atom(ia,1:3), &
                    species, basis_extd, &
                    atom_ignore, &
                    [ this%distributions%bond_info(:)%radius_covalent ] &
               )
       end do
    end do

    viability = viability / real(basis%natom, real32)
  end function evaluate
!###############################################################################

end module raffle__generator