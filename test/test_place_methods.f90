program test_place_methods
  use raffle__io_utils
  use raffle__place_methods
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


  call test_place_method_void(basis, success)


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

  subroutine test_place_method_void(basis, success)
    implicit none
    logical, intent(inout) :: success
    type(basis_type), intent(in) :: basis

    integer :: i
    type(extended_basis_type) :: basis_copy
    logical :: viable
    integer, dimension(3) :: grid
    real(real32), dimension(3) :: grid_offset
    real(real32), dimension(3) :: point
    integer, dimension(:,:), allocatable :: atom_ignore_list
    real(real32), dimension(3) :: tolerance
    real(real32), dimension(2,3) :: bounds

    ! Initialise test data
    grid = [10, 10, 10]
    bounds(1,:) = 0.0_real32
    bounds(2,:) = 1.0_real32
    allocate(atom_ignore_list(1, 2))  ! No atoms to ignore
    atom_ignore_list(1,:) = [1,2]
    grid_offset = [0.5_real32, 0.5_real32, 0.5_real32]

    ! Initialise basis
    call basis_copy%copy(basis)

    ! Call the void subroutine
    point = place_method_void( &
         grid, grid_offset, bounds, basis_copy, &
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