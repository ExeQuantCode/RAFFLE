program test_place_methods
  use raffle__io_utils
  use raffle__place_methods
  use raffle__distribs_container, only: distribs_container_type
  use raffle__constants, only: real32
  use raffle__geom_rw, only: basis_type
  use raffle__geom_extd, only: extended_basis_type
  use raffle__generator, only: raffle_generator_type
  implicit none

  integer :: num_seed, seed
  logical :: viable = .true.
  type(basis_type) :: basis
  type(extended_basis_type) :: basis_extd
  type(raffle_generator_type) :: generator
  real(real32), dimension(3) :: point
  real(real32), dimension(2, 3) :: bounds
  character(3), dimension(1) :: element_symbols
  real(real32), dimension(1) :: element_energies

  integer, dimension(:), allocatable :: seed_arr
  type(basis_type), allocatable :: database(:)
  integer, dimension(:,:), allocatable :: atom_ignore_list
  
  logical :: success = .true.


  test_error_handling = .true.


  !-----------------------------------------------------------------------------
  ! set up database
  !-----------------------------------------------------------------------------
  allocate(database(1))
  database(1)%nspec = 1
  database(1)%natom = 8
  allocate(database(1)%spec(database(1)%nspec))
  database(1)%spec(1)%num = 8
  database(1)%spec(1)%name = 'C'
  allocate(database(1)%spec(1)%atom(database(1)%spec(1)%num, 3))
  database(1)%spec(1)%atom(1, :3) = [0.0, 0.0, 0.0]
  database(1)%spec(1)%atom(2, :3) = [0.5, 0.5, 0.0]
  database(1)%spec(1)%atom(3, :3) = [0.5, 0.0, 0.5]
  database(1)%spec(1)%atom(4, :3) = [0.0, 0.5, 0.5]
  database(1)%spec(1)%atom(5, :3) = [0.25, 0.25, 0.25]
  database(1)%spec(1)%atom(6, :3) = [0.75, 0.75, 0.25]
  database(1)%spec(1)%atom(7, :3) = [0.75, 0.25, 0.75]
  database(1)%spec(1)%atom(8, :3) = [0.25, 0.75, 0.75]

  database(1)%lat(1,:) = [3.5607451090903233, 0.0, 0.0]
  database(1)%lat(2,:) = [0.0, 3.5607451090903233, 0.0]
  database(1)%lat(3,:) = [0.0, 0.0, 3.5607451090903233]
  database(1)%energy = -72.213492


  ! Initialise basis
  basis%nspec = 1
  allocate(basis%spec(basis%nspec))
  basis%spec(1)%name = 'C'
  basis%spec(1)%num = 2
  allocate(basis%spec(1)%atom(basis%spec(1)%num,3))
  basis%spec(1)%atom(1,:) = [0.0_real32, 0.0_real32, 0.0_real32]
  basis%spec(1)%atom(2,:) = [0.5_real32, 0.5_real32, 0.5_real32]
  basis%lat = 0.0_real32
  basis%lat(1,1) = 5.0_real32
  basis%lat(2,2) = 5.0_real32
  basis%lat(3,3) = 5.0_real32


  !-----------------------------------------------------------------------------
  ! set up distribution functions
  !-----------------------------------------------------------------------------
  bounds(1,:) = 0.25_real32
  bounds(2,:) = 0.75_real32
  call generator%set_host( basis )
  element_symbols(1) = 'C'
  element_energies(1) = -9.0266865
  call generator%distributions%set_element_energies( &
       element_symbols, &
       element_energies &
  )
  call generator%distributions%create( &
       basis_list = database, &
       deallocate_systems = .true. &
  )
  call generator%distributions%set_element_map( &
       [ "C  "] &
  )
  call generator%distributions%host_system%set_element_map( &
       generator%distributions%element_info &
  )
  allocate(atom_ignore_list(2,1))
  atom_ignore_list(:,1) = [1,2]
  seed = 0
  call random_seed(size=num_seed)
  allocate(seed_arr(num_seed))
  seed_arr = seed 
  call random_seed(put=seed_arr)
  call basis_extd%copy(basis)
  call basis_extd%create_images( &
       max_bondlength = 6._real32, &
       atom_ignore_list = atom_ignore_list &
  )


  !-----------------------------------------------------------------------------
  ! test place_method_void
  !-----------------------------------------------------------------------------
  call test_place_method_void(generator%distributions, basis, success)


  !-----------------------------------------------------------------------------
  ! test place_method_rand
  !-----------------------------------------------------------------------------
  viable = .true.
  point = place_method_rand( &
       generator%distributions, &
       bounds, &
       basis_extd, &
       atom_ignore_list, &
       radius_list = [ 0.5_real32 ], &
       max_attempts = 1000, &
       viable = viable &
  )
  if(.not. viable) then
     write(*,*) "test_place_method_rand failed"
     success = .false.
  end if
  if( any( point .gt. 0.75_real32 ) .or. any( point .lt. 0.25_real32 ) ) then
     write(*,*) "test_place_method_rand failed"
     success = .false.
  end if


  !-----------------------------------------------------------------------------
  ! test place_method_walk
  !-----------------------------------------------------------------------------
  viable = .true.
  point = place_method_walk( &
       generator%distributions, &
       bounds, &
       basis_extd, &
       atom_ignore_list, &
       radius_list = [ 0.5_real32 ], &
       max_attempts = 1000, &
       step_size_fine = 0.1_real32, &
       step_size_coarse = 0.5_real32, &
       viable = viable &
  )
  if(.not. viable) then
     write(*,*) "test_place_method_walk failed, viable = ", viable
     success = .false.
  end if
  if( any( point .gt. 0.75_real32 ) .or. any( point .lt. 0.25_real32 ) ) then
     write(*,*) "test_place_method_walk failed, point = ", point
     success = .false.
  end if


  !-----------------------------------------------------------------------------
  ! test place_method_growth
  !-----------------------------------------------------------------------------
  viable = .true.
  point = place_method_growth( &
       generator%distributions, &
       prior_point = [0.45_real32, 0.45_real32, 0.45_real32], &
       prior_species = 1, &
       bounds = bounds, &
       basis = basis_extd, &
       atom_ignore_list = atom_ignore_list, &
       radius_list = [ 0.5_real32 ], &
       max_attempts = 1000, &
       step_size_fine = 0.1_real32, &
       step_size_coarse = 0.5_real32, &
       viable = viable &
  )
  if(.not. viable) then
     write(*,*) "test_place_method_growth failed, viable = ", viable
     success = .false.
  end if
  if( any( point .gt. 0.75_real32 ) .or. any( point .lt. 0.25_real32 ) ) then
     write(*,*) "test_place_method_growth failed, point = ", point
     success = .false.
  end if




  !-----------------------------------------------------------------------------
  ! check for any failed tests
  !-----------------------------------------------------------------------------
  write(*,*) "----------------------------------------"
  if(success)then
     write(*,*) 'test_place_methods passed all tests'
  else
     write(0,*) 'test_place_methods failed one or more tests'
     stop 1
  end if

contains

  subroutine test_place_method_void(distributions, basis, success)
    use raffle__viability, only: get_gridpoints_and_viability
    implicit none
    logical, intent(inout) :: success
    type(basis_type), intent(in) :: basis
    type(distribs_container_type), intent(in) :: distributions

    integer :: i
    type(extended_basis_type) :: basis_copy
    logical :: viable
    integer, dimension(3) :: grid
    real(real32), dimension(3) :: grid_offset
    real(real32), dimension(3) :: point
    integer, dimension(:,:), allocatable :: atom_ignore_list
    real(real32), dimension(3) :: tolerance
    real(real32), dimension(2,3) :: bounds
    real(real32), dimension(:,:), allocatable :: gridpoints, viability_grid

    ! Initialise test data
    grid = [10, 10, 10]
    bounds(1,:) = 0.0_real32
    bounds(2,:) = 1.0_real32
    allocate(atom_ignore_list(2,1))  ! No atoms to ignore
    atom_ignore_list(:,1) = [1,2]
    grid_offset = [0.5_real32, 0.5_real32, 0.5_real32]

    ! Initialise basis
    call basis_copy%copy(basis)
    call basis_copy%create_images( &
         max_bondlength = 6._real32, &
         atom_ignore_list = atom_ignore_list &
    )

    !---------------------------------------------------------------------------
    ! set up gridpoints
    !---------------------------------------------------------------------------
    gridpoints = get_gridpoints_and_viability( &
         distributions, &
         grid, &
         bounds, &
         basis_copy, &
         [ 1 ], &
         [ distributions%bond_info(:)%radius_covalent ], &
         atom_ignore_list, &
         grid_offset = grid_offset &
    )

    ! Call the void subroutine
    point = place_method_void( &
         gridpoints, basis_copy, &
         atom_ignore_list, &
         viable &
    )

    ! Check if viable
    call assert(viable, "No viable gridpoints found.", success)

    do i = 1, 3
       tolerance(i) = 1._real32 / real(grid(i),real32) / 2._real32
    end do
    ! Check point is correct
    call assert( &
         all( abs( point - 0.5_real32) .lt. tolerance + 1.E-6_real32 ), &
         "Incorrect gridpoint found.", &
         success &
    )

  end subroutine test_place_method_void


!###############################################################################

  subroutine assert(condition, message, success)
    implicit none
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    logical, intent(inout) :: success
    if (.not. condition) then
       write(0,*) "Test failed: ", message
       success = .false.
    end if
  end subroutine assert

end program test_place_methods