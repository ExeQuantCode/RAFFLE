program test_place_methods
  use raffle__io_utils
  use raffle__viability
  use raffle__distribs_container, only: distribs_container_type
  use raffle__constants, only: real32
  use raffle__geom_rw, only: basis_type
  use raffle__geom_extd, only: extended_basis_type
  implicit none

  type(basis_type) :: basis
  logical :: success = .true.

  test_error_handling = .true.

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


  call test_get_gridpoints_and_viability(basis, success)
  call test_update_gridpoints_and_viability(basis, success)


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

  subroutine test_get_gridpoints_and_viability(basis, success)
    implicit none
    logical, intent(inout) :: success
    type(basis_type), intent(in) :: basis

    integer :: i
    type(extended_basis_type) :: basis_copy
    type(distribs_container_type) :: distribs_container
    integer, dimension(3) :: grid
    integer, dimension(:,:), allocatable :: atom_ignore_list
    real(real32), dimension(:), allocatable :: radius_list
    real(real32) :: lowtol
    real(real32), dimension(:,:), allocatable :: points
    real(real32), dimension(3) :: grid_offset
    real(real32), dimension(2,3) :: bounds

    ! Initialise test data
    grid = [10, 10, 10]
    bounds(1,:) = 0.0_real32
    bounds(2,:) = 1.0_real32
    allocate(atom_ignore_list(2,1))  ! No atoms to ignore
    atom_ignore_list(:,1) = [1,2]
    allocate(radius_list(1))
    radius_list = 1.0_real32
    lowtol = 0.5_real32
    grid_offset = [0.5_real32, 0.5_real32, 0.5_real32]

    ! Initialise basis
    call basis_copy%copy(basis)
    call basis_copy%create_images( &
         max_bondlength = distribs_container%cutoff_max(1), &
         atom_ignore_list = atom_ignore_list &
    )

    ! Initialise gvector container
    call distribs_container%set_element_energies( &
         [basis%spec(:)%name], &
         [ ( 0.0_real32, i = 1, basis%nspec ) ] &
    )
    call distribs_container%create([basis])
    call distribs_container%set_element_map( &
         [ basis%spec(:)%name ] &
    )
    call distribs_container%host_system%set(basis)
    call distribs_container%host_system%set_element_map( &
         distribs_container%element_info &
    )

    ! Call the function to test
    points = get_gridpoints_and_viability( &
         distribs_container, &
         grid, bounds, &
         basis_copy, &
         [ 1 ], &
         radius_list, &
         atom_ignore_list, &
         grid_offset &
    )

    ! Check points exist
    call assert(size(points, 2) .gt. 0, "No viable gridpoints found.", success)

    ! Check number of points
    call assert( &
         size(points, 2) .lt. 1000, &
         "Incorrect number of gridpoints found.", &
         success &
    )
    ! Check number of points
    call assert( &
         size(points, 2) .eq. 907, &
         "Incorrect number of gridpoints found.", &
         success &
    )

  end subroutine test_get_gridpoints_and_viability

  subroutine test_update_gridpoints_and_viability(basis, success)
    implicit none
    logical, intent(inout) :: success
    type(basis_type), intent(in) :: basis

    integer :: i
    type(extended_basis_type) :: basis_copy
    type(distribs_container_type) :: distribs_container
    integer, dimension(3) :: grid
    integer, dimension(:,:), allocatable :: atom_ignore_list
    real(real32), dimension(:), allocatable :: radius_list
    real(real32) :: lowtol
    real(real32), dimension(:,:), allocatable :: points
    real(real32), dimension(3) :: grid_offset
    real(real32), dimension(2,3) :: bounds

    ! Initialise test data
    grid = [10, 10, 10]
    bounds(1,:) = 0.0_real32
    bounds(2,:) = 1.0_real32
    allocate(atom_ignore_list(2,1))  ! No atoms to ignore
    atom_ignore_list(:,1) = [1,2]
    allocate(radius_list(1))
    radius_list = 1.0_real32 !!! NO!!! USING CARBON RADIUS
    lowtol = 0.5_real32
    grid_offset = [0.5_real32, 0.5_real32, 0.5_real32]

    ! Initialise basis
    call basis_copy%copy(basis)
    call basis_copy%create_images( &
         max_bondlength = distribs_container%cutoff_max(1), &
         atom_ignore_list = atom_ignore_list &
    )

    ! Initialise gvector container
    call distribs_container%set_element_energies( &
         [basis%spec(:)%name], &
         [ ( 0.0_real32, i = 1, basis%nspec ) ] &
    )
    call distribs_container%create([basis])
    call distribs_container%set_element_map( &
         [ basis%spec(:)%name ] &
    )
    call distribs_container%host_system%set(basis)
    call distribs_container%host_system%set_element_map( &
         distribs_container%element_info &
    )

    ! Call the function to test
    points = get_gridpoints_and_viability( &
         distribs_container, &
         grid, bounds, &
         basis_copy, &
         [ 1 ], &
         radius_list, &
         atom_ignore_list, &
         grid_offset &
    )

    ! Call the update subroutine
    call update_gridpoints_and_viability( &
         points, distribs_container, basis_copy, &
         [1], &
         [1,2], &
         radius_list, &
         atom_ignore_list &
    )

    ! Check points exist
    call assert(size(points, 2) .gt. 0, "No viable gridpoints found.", success)

    ! Check number of points
    call assert( &
         size(points, 2) .lt. 1000, &
         "Incorrect number of gridpoints found.", &
         success &
    )

    ! Check number of points
    call assert( &
         size(points, 2) .eq. 802, &
         "Incorrect number of gridpoints found.", &
         success &
    )

    ! Call the update subroutine
    distribs_container%radius_distance_tol(1) = 100._real32
    call update_gridpoints_and_viability( &
         points, distribs_container, basis_copy, &
         [1], &
         [1,2], &
         radius_list, &
         atom_ignore_list &
    )

    ! Check all points have been removed
    call assert(.not.allocated(points), "Some grid points remain.", success)

  end subroutine test_update_gridpoints_and_viability


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
