module raffle__generator
  !! Module for generating random structures from host structures.
  !!
  !! This module contains the raffle generator type, which is used to generate
  !! random structures from a host structure. The raffle generator uses
  !! distribution functions to determine the placement of atoms in the
  !! provided host structure.
  use raffle__io_utils, only: stop_program, print_warning, suppress_warnings
  use raffle__constants, only: real32
  use raffle__tools_infile, only: assign_val, assign_vec
  use raffle__misc_linalg, only: modu
  use raffle__misc, only: strip_null, set, shuffle, sort1D, to_upper
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
     real(real32), dimension(2,3) :: bounds = reshape( &
          (/ &
               0.0_real32, 1.0_real32, &
               0.0_real32, 1.0_real32, &
               0.0_real32, 1.0_real32 &
          /), [2,3] &
     )
     !! Bounds for atom placement.
     type(distribs_container_type) :: distributions
     !! Distribution function container for the 2-, 3-, and 4-body interactions.
     integer :: max_attempts = 10000
     !! Limit for the number of attempts to place an atom.
     real(real32) :: &
          walk_step_size_coarse = 1._real32, &
          walk_step_size_fine = 0.1_real32
     !! Step size for the walk and grow methods.
     real(real32), dimension(5) :: method_ratio
     !! Ratio of each placement method.
     type(basis_type), dimension(:), allocatable :: structures
     !! Generated structures.
   contains
     procedure, pass(this) :: set_host
     !! Procedure to set the host structure.
     procedure, pass(this) :: set_grid
     !! Procedure to set the grid for the raffle generator.
     procedure, pass(this) :: reset_grid
     !! Procedure to reset the grid for the raffle generator.
     procedure, pass(this) :: set_bounds
     !! Procedure to set the bounds for the raffle generator.
     procedure, pass(this) :: reset_bounds
     !! Procedure to reset the bounds for the raffle generator.
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

     procedure, pass(this) :: print_settings => print_generator_settings
     !! Procedure to print the raffle generator settings.
     procedure, pass(this) :: read_settings => read_generator_settings
     !! Procedure to read the raffle generator settings.
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
             this%grid(i) = nint( &
                  ( this%bounds(2,i) - this%bounds(1,i) ) * &
                  modu(this%host%lat(i,:)) / this%grid_spacing &
             )
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
  subroutine set_bounds(this, bounds)
    !! Set the bounds for the raffle generator.
    !!
    !! This procedure sets the bounds for the raffle generator. The bounds are
    !! used to determine the placement of atoms in the host structure.
    implicit none

    ! Arguments
    class(raffle_generator_type), intent(inout) :: this
    !! Instance of the raffle generator.
    real(real32), dimension(2,3), intent(in) :: bounds
    !! Bounds for atom placement.

    this%bounds = bounds
    call this%set_grid()

  end subroutine set_bounds
!###############################################################################


!###############################################################################
  subroutine reset_bounds(this)
    !! Reset the grid for the raffle generator.
    implicit none

    ! Arguments
    class(raffle_generator_type), intent(inout) :: this
    !! Instance of the raffle generator.

    this%bounds(1,:) = 0.0_real32
    this%bounds(2,:) = 1.0_real32
  end subroutine reset_bounds
!###############################################################################



!###############################################################################
  subroutine generate(this, num_structures, &
       stoichiometry, method_ratio, seed, settings_out_file, &
       verbose, exit_code &
  )
    !! Generate random structures.
    !!
    !! This procedure generates random structures from the contained host
    !! structure and the stoichiometry argument. The number of structures to
    !! generate is specified by the num_structures argument.
    !! The ratio of placement methods to be sampled is defined by method_ratio.
    implicit none

    ! Arguments
    class(raffle_generator_type), intent(inout) :: this
    !! Instance of the raffle generator.
    integer, intent(in) :: num_structures
    !! Number of structures to generate.
    type(stoichiometry_type), dimension(:), intent(in) :: stoichiometry
    !! Stoichiometry of the structures to generate.
    real(real32), dimension(5), intent(in), optional :: method_ratio
    !! Ratio of each placement method.
    integer, intent(in), optional :: seed
    !! Seed for the random number generator.
    character(*), intent(in), optional :: settings_out_file
    !! File to print the settings to.
    integer, intent(in), optional :: verbose
    !! Verbosity level.
    integer, intent(out), optional :: exit_code
    !! Exit code.

    ! Local variables
    integer :: i, j, k, istructure, num_structures_old, num_structures_new
    !! Loop counters.
    integer :: exit_code_
    !! Exit code.
    integer :: num_seed
    !! Number of seeds for the random number generator.
    integer :: num_insert_atoms, num_insert_species
    !! Number of atoms and species to insert (from stoichiometry).
    real(real32) :: ratio_norm
    !! Normalisation factor for the method ratios.
    logical :: success
    !! Boolean comparison of element symbols.
    integer :: verbose_
    !! Verbosity level.
    logical :: suppress_warnings_store
    !! Boolean to store the suppress_warnings value.
    type(basis_type) :: basis_template
    !! Basis of the structure to generate (i.e. allocated species and atoms).
    real(real32), dimension(5) :: &
         method_rand_limit = &
              [1.0_real32, 0.1_real32, 0.5_real32, 0.5_real32, 1.0_real32]
    !! Default ratio of each placement method.

    integer, dimension(:), allocatable :: seed_arr
    !! Array of seeds for the random number generator.
    type(basis_type), dimension(:), allocatable :: tmp_structures
    !! Temporary array of structures (for memory reallocation).

    integer, dimension(:,:), allocatable :: placement_list
    !! List of possible atoms to place in the structure.


    !---------------------------------------------------------------------------
    ! Set the verbosity level
    !---------------------------------------------------------------------------
    exit_code_ = 0
    verbose_ = 0
    if(present(verbose)) verbose_ = verbose
    if(verbose_ .eq. 0)then
       suppress_warnings_store = suppress_warnings
       suppress_warnings = .true.
    end if


    !---------------------------------------------------------------------------
    ! Set the placement method selection limit numbers
    !---------------------------------------------------------------------------
    if(verbose_.gt.0) write(*,*) "Setting method ratios"
    if(present(method_ratio))then
       method_rand_limit = method_ratio
    else
       method_rand_limit = this%method_ratio
    end if
    this%method_ratio = method_rand_limit
    ratio_norm = real(sum(method_rand_limit), real32)
    method_rand_limit = method_rand_limit / ratio_norm
    do i = 2, 5, 1
       method_rand_limit(i) = method_rand_limit(i) + method_rand_limit(i-1)
    end do
    if(verbose_.gt.0) write(*,*) &
         "Method random limits (void, rand, walk, grow, min): ", &
         method_rand_limit


    !---------------------------------------------------------------------------
    ! Print the settings to a file
    !---------------------------------------------------------------------------
    if(present(settings_out_file))then
       if(trim(settings_out_file).ne."")then
          call this%print_settings(settings_out_file)
       end if
    end if

    !---------------------------------------------------------------------------
    ! Set the random seed
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
    call this%distributions%set_element_map( &
         [ basis_template%spec(:)%name ] &
    )
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
    allocate(placement_list(2, num_insert_atoms))
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
             placement_list(1,k) = i
             placement_list(2,k) = j
          end do
       else
          do j = 1, basis_template%spec(i)%num
             if(j.le.this%host%spec(i)%num) cycle
             k = k + 1
             placement_list(1,k) = i
             placement_list(2,k) = j
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
                 method_rand_limit, &
                 verbose_, &
                 exit_code_ &
            ) &
       )
       this%num_structures = istructure

    end do structure_loop
    if(verbose_ .gt. 0 .and. exit_code_ .eq. 0) &
         write(*,*) "Finished generating structures"

    if(verbose_ .eq. 0)then
       suppress_warnings = suppress_warnings_store
    end if

    if(present(exit_code))then
       exit_code = exit_code_
    elseif(exit_code_ .ne. 0)then
       call stop_program("Error generating structures", exit_code_)
    end if

  end subroutine generate
!###############################################################################


!###############################################################################
  function generate_structure( &
       this, &
       basis_initial, &
       placement_list, method_rand_limit, verbose, &
       exit_code &
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
    real(real32), dimension(5) :: method_rand_limit
    !! Upper limit of the random number to call each placement method.
    type(extended_basis_type) :: basis
    !! Generated basis.
    integer, intent(in) :: verbose
    !! Verbosity level.
    integer, intent(inout) :: exit_code
    !! Exit code.

    ! Local variables
    integer :: iplaced, void_ticker
    !! Loop counters.
    integer :: num_insert_atoms
    !! Number of atoms to insert.
    real(real32) :: rtmp1
    !! Random number.
    logical :: skip_min
    !! Boolean for skipping the minimum method.
    logical :: viable
    !! Boolean for viable placement.
    logical :: placement_aborted
    !! Boolean for aborted placement.
    integer, dimension(size(placement_list,1),size(placement_list,2)) :: &
         placement_list_shuffled
    !! Shuffled placement list.
    real(real32), dimension(3) :: point
    !! Coordinate of the atom to place.
    real(real32), dimension(5) :: method_rand_limit_, method_rand_limit_store
    !! Temporary random limit of each placement method.
    !! This is used to update the contribution of the global minimum method if
    !! no viable gridpoints are found.
    integer, dimension(:), allocatable :: species_index_list
    !! List of species indices to add.
    real(real32), dimension(:,:), allocatable :: gridpoint_viability
    !! Viable gridpoints for placing atoms.
    character(len=256) :: stop_msg, warn_msg
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
    call shuffle(placement_list_shuffled,2)


    !---------------------------------------------------------------------------
    ! generate species index list to add
    !---------------------------------------------------------------------------
    species_index_list = placement_list_shuffled(1,:)
    call set(species_index_list)


    !---------------------------------------------------------------------------
    ! check for viable gridpoints
    !---------------------------------------------------------------------------
    method_rand_limit_ = method_rand_limit
    gridpoint_viability = get_gridpoints_and_viability( &
         this%distributions, &
         this%grid, &
         this%bounds, &
         basis, &
         species_index_list, &
         [ this%distributions%bond_info(:)%radius_covalent ], &
         placement_list_shuffled, &
         this%grid_offset &
    )


    !---------------------------------------------------------------------------
    ! place the atoms according to the method ratios
    !---------------------------------------------------------------------------
    iplaced = 0
    void_ticker = 0
    viable = .false.
    skip_min = .false.
    placement_aborted = .false.
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
                    [ placement_list_shuffled(:,iplaced) ], &
                    [ this%distributions%bond_info(:)%radius_covalent ], &
                    placement_list_shuffled(:,iplaced+1:) &
               )
       end if
       if(.not.allocated(gridpoint_viability))then
          write(warn_msg, '("No more viable gridpoints")')
          warn_msg = trim(warn_msg) // &
               achar(13) // achar(10) // &
               "Stopping atom placement for this structure"
          call print_warning(warn_msg)
          placement_aborted = .true.
          exit placement_loop
       end if
       viable = .false.
       !------------------------------------------------------------------------
       ! Choose a placement method
       ! call a random number and query the method ratios
       !------------------------------------------------------------------------
       call random_number(rtmp1)
       if(rtmp1.le.method_rand_limit_(1)) then
          if(verbose.gt.0) write(*,*) "Add Atom Void"
          point = place_method_void( &
               gridpoint_viability, &
               basis, &
               placement_list_shuffled(:,iplaced+1:), viable &
          )
       elseif(rtmp1.le.method_rand_limit_(2)) then
          if(verbose.gt.0) write(*,*) "Add Atom Random"
          point = place_method_rand( &
               this%distributions, &
               this%bounds, &
               basis, &
               placement_list_shuffled(:,iplaced+1:), &
               [ this%distributions%bond_info(:)%radius_covalent ], &
               this%max_attempts, &
               viable &
          )
       elseif(rtmp1.le.method_rand_limit_(3)) then
          if(verbose.gt.0) write(*,*) "Add Atom Walk"
          point = place_method_walk( &
               this%distributions, &
               this%bounds, &
               basis, &
               placement_list_shuffled(:,iplaced+1:), &
               [ this%distributions%bond_info(:)%radius_covalent ], &
               this%max_attempts, &
               this%walk_step_size_coarse, this%walk_step_size_fine, &
               viable &
          )
       elseif(rtmp1.le.method_rand_limit_(4)) then
          if(iplaced.eq.0)then
             if(verbose.gt.0) write(*,*) "Add Atom Random (growth seed)"
             point = place_method_rand( &
                  this%distributions, &
                  this%bounds, &
                  basis, &
                  placement_list_shuffled(:,iplaced+1:), &
                  [ this%distributions%bond_info(:)%radius_covalent ], &
                  this%max_attempts, &
                  viable &
             )
          else
             if(verbose.gt.0) write(*,*) "Add Atom Growth"
             point = place_method_growth( &
                  this%distributions, &
                  basis%spec(placement_list_shuffled(1,iplaced))%atom( &
                       placement_list_shuffled(2,iplaced),:3 &
                  ), &
                  placement_list_shuffled(1,iplaced), &
                  this%bounds, &
                  basis, &
                  placement_list_shuffled(:,iplaced+1:), &
                  [ this%distributions%bond_info(:)%radius_covalent ], &
                  this%max_attempts, &
                  this%walk_step_size_coarse, this%walk_step_size_fine, &
                  viable &
             )
          end if
       elseif(rtmp1.le.method_rand_limit_(5)) then
          if(verbose.gt.0) write(*,*) "Add Atom Minimum"
          point = place_method_min( gridpoint_viability, &
               placement_list_shuffled(1,iplaced+1), &
               species_index_list, &
               viable &
          )
          if(.not. viable .and. abs(method_rand_limit_(4)).lt.1.E-6)then
             write(warn_msg, &
                  '("Minimum method failed, no other methods available")' &
             )
             warn_msg = trim(warn_msg) // &
                  achar(13) // achar(10) // &
                  "Stopping atom placement for this structure"
             call print_warning(warn_msg)
             placement_aborted = .true.
             exit placement_loop
          elseif(.not.viable)then
             skip_min = .true.
             method_rand_limit_store = method_rand_limit_
             method_rand_limit_ = method_rand_limit_ / method_rand_limit_(4)
             method_rand_limit_(5) = method_rand_limit_(4)
          end if
       end if
       !------------------------------------------------------------------------
       ! check if the placement method returned a viable point
       ! if not, cycle the loop
       !------------------------------------------------------------------------
       if(.not. viable)then
          void_ticker = void_ticker + 1
          if(void_ticker.gt.10.and..not.allocated(gridpoint_viability))then
             write(warn_msg, '("No more viable gridpoints")')
             warn_msg = trim(warn_msg) // &
                  achar(13) // achar(10) // &
                  "Stopping atom placement for this structure"
             call print_warning(warn_msg)
             placement_aborted = .true.
             exit placement_loop
          elseif(void_ticker.gt.10)then
             point = place_method_void( &
                  gridpoint_viability, basis, &
                  placement_list_shuffled(:,iplaced+1:), viable &
             )
             void_ticker = 0
          end if
          if(.not.viable) cycle placement_loop
       end if
       !------------------------------------------------------------------------
       ! place the atom and update the image atoms in the basis
       !------------------------------------------------------------------------
       if(skip_min)then
          method_rand_limit_ = method_rand_limit_store
          skip_min = .false.
       end if
       iplaced = iplaced + 1
       basis%spec(placement_list_shuffled(1,iplaced))%atom( &
            placement_list_shuffled(2,iplaced),:3) = point(:3)
       call basis%update_images( &
            max_bondlength = this%distributions%cutoff_max(1), &
            is = placement_list_shuffled(1,iplaced), &
            ia = placement_list_shuffled(2,iplaced) &
       )
       if(verbose.gt.0)then
          write(*,'(A)',ADVANCE='NO') achar(13)
          write(*,'(2X,"placed atom ",I0," [",I0,",",I0,"] at",3(1X,F6.3))') &
               iplaced, placement_list_shuffled(1:2,iplaced), point(:3)
          write(*,*)
       end if

    end do placement_loop

    if(placement_aborted)then
       call stop_program( &
            "Placement routine aborted, not all atoms placed", &
            block_stop = .true. &
       )
       exit_code = 1
       call basis%remove_atoms(placement_list_shuffled(:,iplaced+1:))
    end if

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
    class(raffle_generator_type), intent(inout) :: this
    !! Instance of the raffle generator.
    type(basis_type), intent(in) :: basis
    !! Basis of the structure to evaluate.
    real(real32) :: viability
    !! Viability of the generated structures.

    ! Local variables
    integer :: is, ia, species
    !! Loop indices.
    integer, dimension(2,1) :: atom_ignore
    !! Atom to ignore.
    type(extended_basis_type) :: basis_extd
    !! Extended basis for the structure to evaluate.


    call basis_extd%copy(basis)
    call basis_extd%create_images( &
         max_bondlength = this%distributions%cutoff_max(1) &
    )
    viability = 0.0_real32
    call this%distributions%set_element_map( &
         [ basis_extd%spec(:)%name ] &
    )
    do is = 1, basis%nspec
       species = this%distributions%get_element_index( &
            basis_extd%spec(is)%name &
       )
       if(species.eq.0)then
          call stop_program( &
               "Species "//&
               trim(basis_extd%spec(is)%name)//&
               " not found in distribution functions" &
          )
          return
       end if
       do ia = 1, basis%spec(is)%num
          atom_ignore(:,1) = [is,ia]
          viability = viability + &
               evaluate_point( this%distributions, &
                    [ basis%spec(is)%atom(ia,1:3) ], &
                    species, basis_extd, &
                    atom_ignore, &
                    [ this%distributions%bond_info(:)%radius_covalent ] &
               )
       end do
    end do

    viability = viability / real(basis%natom, real32)
  end function evaluate
!###############################################################################


!###############################################################################
  subroutine print_generator_settings(this, file)
    !! Print the raffle generator settings.
    implicit none

    ! Arguments
    class(raffle_generator_type), intent(in) :: this
    !! Instance of the raffle generator.
    character(*), intent(in) :: file
    !! Filename to write the settings to.

    ! Local variables
    integer :: i
    !! Loop index.
    integer :: unit
    !! Unit number for the file.

    ! Open the file
    open(newunit=unit, file=file)

    write(unit,'("# RAFFLE Generator Settings")')
    write(unit,'("# GENERATOR SETTINGS")')
    write(unit,'("HOST_LATTICE # not a setting, just for reference")')
    write(unit,'("  ",3(1X,F5.2))') this%host%lat(1,:)
    write(unit,'("  ",3(1X,F5.2))') this%host%lat(2,:)
    write(unit,'("  ",3(1X,F5.2))') this%host%lat(3,:)
    write(unit,'("END HOST_LATTICE")')

    write(unit,'("GRID =",3(1X,I0))') this%grid
    write(unit,'("GRID_OFFSET =",3(1X,F15.9))') this%grid_offset
    write(unit,'("GRID_SPACING = ",F15.9)') this%grid_spacing
    write(unit,'("BOUNDS_LW =",3(1X,F15.9))') this%bounds(1,:)
    write(unit,'("BOUNDS_UP =",3(1X,F15.9))') this%bounds(2,:)

    write(unit,'("MAX_ATTEMPTS =",I0)') this%max_attempts
    write(unit,'("WALK_STEP_SIZE_COARSE = ",F15.9)') this%walk_step_size_coarse
    write(unit,'("WALK_STEP_SIZE_FINE = ",F15.9)') this%walk_step_size_fine
    write(unit,'("METHOD_VOID = ",F15.9)') this%method_ratio(1)
    write(unit,'("METHOD_RANDOM = ",F15.9)') this%method_ratio(2)
    write(unit,'("METHOD_WALK = ",F15.9)') this%method_ratio(3)
    write(unit,'("METHOD_GROW = ",F15.9)') this%method_ratio(4)
    write(unit,'("METHOD_MIN = ",F15.9)') this%method_ratio(5)

    write(unit,'("# DISTRIBUTION SETTINGS")')
    write(unit,'("KBT = ",F5.2)') this%distributions%kbt
    write(unit,'("SIGMA =",3(1X,F15.9))') this%distributions%sigma
    write(unit,'("WIDTH =",3(1X,F15.9))') this%distributions%width
    write(unit,'("CUTOFF_MIN =",3(1X,F15.9))') this%distributions%cutoff_min
    write(unit,'("CUTOFF_MAX =",3(1X,F15.9))') this%distributions%cutoff_max
    write(unit,'("RADIUS_DISTANCE_TOLERANCE =",4(1X,F15.9))') &
         this%distributions%radius_distance_tol
    write(unit,'("ELEMENT_INFO # element : energy")')
    do i = 1, size(this%distributions%element_info)
       write(unit,'("  ",A," : ",F15.9)') &
            this%distributions%element_info(i)%name, &
            this%distributions%element_info(i)%energy
    end do
    write(unit,'("END ELEMENT_INFO")')
    write(unit,'("BOND_INFO # element1 element2 : radius")')
    do i = 1, size(this%distributions%bond_info)
       write(unit,'("  ",A," ",A," : ",F15.9)') &
            this%distributions%bond_info(i)%element(1), &
            this%distributions%bond_info(i)%element(2), &
            this%distributions%bond_info(i)%radius_covalent
    end do
    write(unit,'("END BOND_INFO")')

    close(unit)

  end subroutine print_generator_settings
!###############################################################################


!###############################################################################
  subroutine read_generator_settings(this, file)
    !! Read the raffle generator settings.
    implicit none

    ! Arguments
    class(raffle_generator_type), intent(inout) :: this
    !! Instance of the raffle generator.
    character(*), intent(in) :: file
    !! Filename to read the settings from.

    ! Local variables
    integer :: i
    !! Loop index.
    integer :: itmp1, status
    !! Temporary integer.
    integer :: unit
    !! Unit number for the file.
    logical :: exist
    !! Boolean for file existence.
    character(len=256) :: line, buffer, tag
    !! Line from the file.
    character(3), dimension(2) :: elements
    !! Element symbols.
    real(real32) :: rtmp1
    !! Temporary real number.

    ! Check if the file exists
    inquire(file=file, exist=exist)
    if(.not.exist)then
       call stop_program("File does not exist")
       return
    end if

    ! Open the file
    open(newunit=unit, file=file)

    if(allocated(this%distributions%element_info)) &
         deallocate(this%distributions%element_info)
    if(allocated(this%distributions%bond_info)) &
         deallocate(this%distributions%bond_info)
    itmp1 = 0
    do
       read(unit, '(A)', iostat = status) line
       ! encounter end of line
       if(status.ne.0) exit

       if(index(line,'#').gt.0) line = line(1:index(line,'#')-1)
       line = to_upper(trim(adjustl(line)))
       if(len(trim(line)).eq.0) cycle

       tag=trim(adjustl(line))
       if(scan(line,"=").ne.0) tag=trim(tag(:scan(tag,"=")-1))

       select case(trim(adjustl(tag)))
       case("HOST_LATTICE")
          do i = 1, 4
             read(unit,*)
          end do
       case("GRID")
          call assign_vec(line, this%grid, itmp1)
       case("GRID_OFFSET")
          call assign_vec(line, this%grid_offset, itmp1)
       case("GRID_SPACING")
          call assign_val(line, this%grid_spacing, itmp1)
       case("BOUNDS_LW")
          call assign_vec(line, this%bounds(1,:), itmp1)
       case("BOUNDS_UP")
          call assign_vec(line, this%bounds(2,:), itmp1)
       case("MAX_ATTEMPTS")
          call assign_val(line, this%max_attempts, itmp1)
       case("WALK_STEP_SIZE_COARSE")
          call assign_val(line, this%walk_step_size_coarse, itmp1)
       case("WALK_STEP_SIZE_FINE")
          call assign_val(line, this%walk_step_size_fine, itmp1)
       case("METHOD_VOID")
          call assign_val(line, this%method_ratio(1), itmp1)
       case("METHOD_RANDOM")
          call assign_val(line, this%method_ratio(2), itmp1)
       case("METHOD_WALK")
          call assign_val(line, this%method_ratio(3), itmp1)
       case("METHOD_GROW")
          call assign_val(line, this%method_ratio(4), itmp1)
       case("METHOD_MIN")
          call assign_val(line, this%method_ratio(5), itmp1)
       case("KBT")
          call assign_val(line, this%distributions%kbt, itmp1)
       case("SIGMA")
          call assign_vec(line, this%distributions%sigma, itmp1)
       case("WIDTH")
          call assign_vec(line, this%distributions%width, itmp1)
       case("CUTOFF_MIN")
          call assign_vec(line, this%distributions%cutoff_min, itmp1)
       case("CUTOFF_MAX")
          call assign_vec(line, this%distributions%cutoff_max, itmp1)
       case("RADIUS_DISTANCE_TOLERANCE")
          call assign_vec(line, this%distributions%radius_distance_tol, itmp1)
       case("ELEMENT_INFO")
          do
             read(unit,'(A)') line
             if(index(line,'#').gt.0) line = line(1:index(line,'#')-1)
             line = to_upper(trim(adjustl(line)))
             if(len(trim(line)).eq.0) exit
             if(index(line,'END').gt.0) exit
             read(line(:scan(line,":")-1),*) elements(1)
             read(line(scan(line,":")+1:),*) rtmp1
             call this%distributions%set_element_energy(elements(1), rtmp1)
          end do
       case("BOND_INFO")
          do
             read(unit,'(A)') line
             if(index(line,'#').gt.0) line = line(1:index(line,'#')-1)
             line = to_upper(trim(adjustl(line)))
             if(len(trim(line)).eq.0) exit
             if(index(line,'END').gt.0) exit
             read(line(:scan(line,":")-1),*) elements(1), elements(2)
             read(line(scan(line,":")+1:),*) rtmp1
             call this%distributions%set_bond_radius(elements, rtmp1)
          end do
       end select
    end do

    close(unit)

  end subroutine read_generator_settings
!###############################################################################

end module raffle__generator
