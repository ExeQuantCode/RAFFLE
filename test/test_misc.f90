program test_misc
  use error_handling
  use misc_raffle
  use constants, only: real12
  implicit none

  logical :: success = .true.

  test_error_handling = .true.

  call test_sort_str(success)
  call test_sort_str_order(success)
  call test_isort1D(success)
  call test_rsort1D(success)
  call test_iset(success)
  call test_rset(success)
  call test_cset(success)


  !-----------------------------------------------------------------------------
  ! check for any failed tests
  !-----------------------------------------------------------------------------
  write(*,*) "----------------------------------------"
  if(success)then
     write(*,*) 'test_misc passed all tests'
  else
     write(0,*) 'test_misc failed one or more tests'
     stop 1
  end if

contains

  subroutine test_sort_str(success)
    implicit none
    logical, intent(inout) :: success
    character(len=20), dimension(5) :: list
    character(len=20), dimension(5) :: expected_list
    list = [ &
         'banana    ', 'apple     ', 'cherry    ', 'date      ', 'elderberry' &
    ]
    expected_list = [ &
         'apple     ', 'banana    ', 'cherry    ', 'date      ', 'elderberry' &
    ]
    call sort_str(list)
    call assert( &
         all(list .eq. expected_list), &
         'test_sort_str failed', success &
    )
  end subroutine test_sort_str

  subroutine test_sort_str_order(success)
    implicit none
    logical, intent(inout) :: success
    character(len=20), dimension(5) :: list
    integer, dimension(5) :: expected_order = (/ 2, 1, 3, 4, 5 /)
    integer, dimension(:), allocatable :: order
    list = [ &
         'banana    ', 'apple     ', 'cherry    ', 'date      ', 'elderberry' &
    ]
    order = sort_str_order(list)
    call assert( &
         all(order .eq. expected_order), &
         'test_sort_str_order failed', success &
    )
  end subroutine test_sort_str_order

  subroutine test_isort1D(success)
    implicit none
    logical, intent(inout) :: success
    integer, dimension(5) :: arr = [5, 3, 4, 1, 2]
    integer, dimension(5) :: expected_arr = [1, 2, 3, 4, 5]
    call sort1D(arr)
    call assert( &
         all(arr .eq. expected_arr), &
         'test_isort1D failed', success &
    )
  end subroutine test_isort1D

  subroutine test_rsort1D(success)
    implicit none
    logical, intent(inout) :: success
    real(real12), dimension(5) :: arr = [5.0_real12, 3.0_real12, 4.0_real12, 1.0_real12, 2.0_real12]
    real(real12), dimension(5) :: expected_arr = [1.0_real12, 2.0_real12, 3.0_real12, 4.0_real12, 5.0_real12]
    call sort1D(arr)
    call assert( &
         all( abs(arr - expected_arr) .lt. 1.E-6), &
         'test_rsort1D failed', success &
    )
  end subroutine test_rsort1D

  subroutine test_iset(success)
    implicit none
    logical, intent(inout) :: success
    integer, dimension(:), allocatable :: arr
    integer, dimension(:), allocatable :: expected_arr
    allocate(arr(6))
    arr = [1, 2, 2, 3, 3, 3]
    allocate(expected_arr(3))
    expected_arr = [1, 2, 3]
    call set(arr)
    call assert( &
         all(arr .eq. expected_arr), &
         'test_iset failed', success &
    )
  end subroutine test_iset

  subroutine test_rset(success)
    implicit none
    logical, intent(inout) :: success
    real(real12), dimension(:), allocatable :: arr
    real(real12), dimension(:), allocatable :: expected_arr
    allocate(arr(6))
    arr = [1.0_real12, 2.0_real12, 2.0_real12, 3.0_real12, 3.0_real12, 3.0_real12]
    allocate(expected_arr(3))
    expected_arr = [1.0_real12, 2.0_real12, 3.0_real12]
    call set(arr)
    call assert( &
         all( abs(arr - expected_arr) .lt. 1.E-6), &
         'test_rset failed', success &
    )
  end subroutine test_rset

  subroutine test_cset(success)
    implicit none
    logical, intent(inout) :: success
    character(len=20), dimension(:), allocatable :: arr
    character(len=20), dimension(:), allocatable :: expected_arr
    allocate(arr(6))
    arr(:) = [ 'apple ', 'banana', 'banana', 'cherry', 'cherry', 'cherry' ]
    allocate(expected_arr(3))
    expected_arr(:) = [ 'apple ', 'banana', 'cherry' ]
    call set(arr)
    call assert( &
         all(arr .eq. expected_arr), &
         'test_cset failed', success &
    )
  end subroutine test_cset

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

end program test_misc
