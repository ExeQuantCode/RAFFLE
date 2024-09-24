program test_atom_adder
  use error_handling
  use add_atom
  use evolver, only: gvector_container_type
  use constants, only: real12
  use rw_geom, only: basis_type
  use extended_geom, only: extended_basis_type
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
  basis%spec(1)%atom(1,:) = [0.0_real12, 0.0_real12, 0.0_real12]
  basis%spec(1)%atom(2,:) = [0.5_real12, 0.5_real12, 0.5_real12]
  basis%lat = 0.0_real12
  basis%lat(1,1) = 5.0_real12
  basis%lat(2,2) = 5.0_real12
  basis%lat(3,3) = 5.0_real12


  call test_get_gridpoints_and_viability(basis, success)
  call test_update_gridpoints_and_viability(basis, success)
  call test_add_atom_void(basis, success)


  !-----------------------------------------------------------------------------
  ! check for any failed tests
  !-----------------------------------------------------------------------------
  write(*,*) "----------------------------------------"
  if(success)then
     write(*,*) 'test_add_atom passed all tests'
  else
     write(0,*) 'test_add_atom failed one or more tests'
     stop 1
  end if

contains

  subroutine test_get_gridpoints_and_viability(basis, success)
    implicit none
    logical, intent(inout) :: success
    type(basis_type), intent(in) :: basis

    integer :: i
    type(extended_basis_type) :: basis_copy
    type(gvector_container_type) :: gvector_container
    integer, dimension(3) :: grid
    integer, dimension(:,:), allocatable :: atom_ignore_list
    real(real12), dimension(:), allocatable :: radius_list
    real(real12) :: lowtol
    real(real12), dimension(:,:), allocatable :: points
    real(real12), dimension(3) :: grid_offset

    ! Initialise test data
    grid = [10, 10, 10]
    allocate(atom_ignore_list(1, 2))  ! No atoms to ignore
    atom_ignore_list(1,:) = [1,2]
    allocate(radius_list(1))
    radius_list = 1.0_real12
    lowtol = 0.5_real12
    grid_offset = [0.5_real12, 0.5_real12, 0.5_real12]

    ! Initialise basis
    call basis_copy%copy(basis)
    call basis_copy%create_images( &
         max_bondlength = gvector_container%cutoff_max(1), &
         atom_ignore_list = atom_ignore_list &
    )

    ! Initialise gvector container
    call gvector_container%set_element_energies( &
         [basis%spec(:)%name], &
         [ ( 0.0_real12, i = 1, basis%nspec ) ] &
    )
    call gvector_container%create([basis])

    ! Call the function to test
    points = get_gridpoints_and_viability( &
         gvector_container, &
         grid, basis_copy, &
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
         size(points, 2) .eq. 864, &
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
    type(gvector_container_type) :: gvector_container
    integer, dimension(3) :: grid
    integer, dimension(:,:), allocatable :: atom_ignore_list
    real(real12), dimension(:), allocatable :: radius_list
    real(real12) :: lowtol
    real(real12), dimension(:,:), allocatable :: points
    real(real12), dimension(3) :: grid_offset

    ! Initialise test data
    grid = [10, 10, 10]
    allocate(atom_ignore_list(1, 2))  ! No atoms to ignore
    atom_ignore_list(1,:) = [1,2]
    allocate(radius_list(1))
    radius_list = 1.0_real12 !!! NO!!! USING CARBON RADIUS
    lowtol = 0.5_real12
    grid_offset = [0.5_real12, 0.5_real12, 0.5_real12]

    ! Initialise basis
    call basis_copy%copy(basis)
    call basis_copy%create_images( &
         max_bondlength = gvector_container%cutoff_max(1), &
         atom_ignore_list = atom_ignore_list &
    )

    ! Initialise gvector container
    call gvector_container%set_element_energies( &
         [basis%spec(:)%name], &
         [ ( 0.0_real12, i = 1, basis%nspec ) ] &
    )
    call gvector_container%create([basis])

    ! Call the function to test
    points = get_gridpoints_and_viability( &
         gvector_container, &
         grid, basis_copy, &
         [ 1 ], &
         radius_list, &
         atom_ignore_list, &
         grid_offset &
    )

    ! Call the update subroutine
    call update_gridpoints_and_viability( &
         points, gvector_container, basis_copy, &
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
         size(points, 2) .eq. 728, &
         "Incorrect number of gridpoints found.", &
         success &
    )

    ! Call the update subroutine
    gvector_container%radius_distance_tol(1) = 100._real12
    call update_gridpoints_and_viability( &
         points, gvector_container, basis_copy, &
         [1], &
         [1,2], &
         radius_list, &
         atom_ignore_list &
    )

    ! Check all points have been removed
    call assert(.not.allocated(points), "Some grid points remain.", success)

  end subroutine test_update_gridpoints_and_viability

  subroutine test_add_atom_void(basis, success)
    implicit none
    logical, intent(inout) :: success
    type(basis_type), intent(in) :: basis

    integer :: i
    type(extended_basis_type) :: basis_copy
    logical :: viable
    integer, dimension(3) :: grid
    real(real12), dimension(3) :: grid_offset
    real(real12), dimension(3) :: point
    integer, dimension(:,:), allocatable :: atom_ignore_list
    real(real12), dimension(3) :: tolerance

    ! Initialise test data
    grid = [10, 10, 10]
    allocate(atom_ignore_list(1, 2))  ! No atoms to ignore
    atom_ignore_list(1,:) = [1,2]
    grid_offset = [0.5_real12, 0.5_real12, 0.5_real12]

    ! Initialise basis
    call basis_copy%copy(basis)

    ! Call the void subroutine
    point = add_atom_void( &
         grid, grid_offset, basis_copy, &
         atom_ignore_list, &
         viable &
    )

    ! Check if viable
    call assert(viable, "No viable gridpoints found.", success)

    do i = 1, 3
       tolerance(i) = 1._real12 / real(grid(i),real12) / 2._real12
    end do
    ! Check point is correct
    call assert( &
         all( abs( point - 0.5_real12) .lt. tolerance + 1.E-6_real12 ), &
         "Incorrect gridpoint found.", &
         success &
    )

  end subroutine test_add_atom_void


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

end program test_atom_adder