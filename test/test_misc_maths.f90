program test_misc_maths
  use raffle__io_utils, only: test_error_handling
  use raffle__misc_maths
  use raffle__constants, only: real32
  implicit none

  logical :: success = .true.

  test_error_handling = .true.


  call test_lnsum(success)
  call test_triangular_number(success)
  call test_set_difference(success)


  !-----------------------------------------------------------------------------
  ! check for any failed tests
  !-----------------------------------------------------------------------------
  write(*,*) "----------------------------------------"
  if(success)then
     write(*,*) 'test_misc_maths passed all tests'
  else
     write(0,*) 'test_misc_maths failed one or more tests'
     stop 1
  end if

contains

  subroutine test_lnsum(success)
    implicit none
    logical, intent(inout) :: success
    integer :: n
    real(real32) :: result

    n = 5
    result = lnsum(n)
    call assert( &
         abs( &
              result - &
              ( &
                   log(1.0_real32) + &
                   log(2.0_real32) + &
                   log(3.0_real32) + &
                   log(4.0_real32) + &
                   log(5.0_real32) &
              ) &
         ) .lt. 1.E-6_real32, &
         'lnsum failed', &
         success &
    )
  end subroutine test_lnsum

  subroutine test_triangular_number(success)
    implicit none
    logical, intent(inout) :: success
    integer :: n, result

    n = 5
    result = triangular_number(n)
    call assert( &
         result .eq. 15, &
         'Triangular number failed', &
         success &
    )
  end subroutine test_triangular_number

  subroutine test_set_difference(success)
    implicit none
    logical, intent(inout) :: success
    real(real32), dimension(3) :: a, b, result, expected
    real(real32), dimension(4) :: c

    a = [1.0_real32, 2.0_real32, 3.0_real32]
    b = [1.0_real32, 1.0_real32, 1.0_real32]
    expected = [0.0_real32, 1.0_real32, 2.0_real32]
    result = set_difference(a, b)

    call assert( &
         all( abs(result - expected) .lt. 1.E-6_real32 ), &
         'Set difference failed', &
         success &
    )

    b = [0.0_real32, 1.0_real32, 4.0_real32]
    expected = [1.0_real32, 1.0_real32, 0.0_real32]
    result = set_difference(a, b, set_min_zero=.true.)

    call assert( &
         all( abs(result - expected) .lt. 1.E-6_real32 ), &
         'Set difference min zero failed', &
         success &
    )

    c = [1.0_real32, 2.0_real32, 3.0_real32, 4.0_real32]
    write(*,*) "Testing set_difference error handling"
    result = set_difference(a, c)
    write(*,*) "Handled error: set difference of arrays of different lengths"


  end subroutine test_set_difference

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

end program test_misc_maths