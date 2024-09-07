program test_misc_maths
  use error_handling, only: test_error_handling
  use misc_maths
  use constants, only: real12
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
     write(*,*) 'test_rw_geom passed all tests'
  else
     write(0,*) 'test_rw_geom failed one or more tests'
     stop 1
  end if

contains

  subroutine test_lnsum(success)
    implicit none
    logical, intent(inout) :: success
    integer :: n
    real(real12) :: result

    n = 5
    result = lnsum(n)
    call assert( &
         abs( &
              result - &
              ( &
                   log(1.0_real12) + &
                   log(2.0_real12) + &
                   log(3.0_real12) + &
                   log(4.0_real12) + &
                   log(5.0_real12) &
              ) &
         ) .lt. 1.E-6_real12, &
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
    real(real12), dimension(3) :: a, b, result, expected
    real(real12), dimension(4) :: c

    a = [1.0_real12, 2.0_real12, 3.0_real12]
    b = [1.0_real12, 1.0_real12, 1.0_real12]
    expected = [0.0_real12, 1.0_real12, 2.0_real12]
    result = set_difference(a, b)

    call assert( &
         all( abs(result - expected) .lt. 1.E-6_real12 ), &
         'Set difference failed', &
         success &
    )

    b = [0.0_real12, 1.0_real12, 4.0_real12]
    expected = [1.0_real12, 1.0_real12, 0.0_real12]
    result = set_difference(a, b, set_min_zero=.true.)

    call assert( &
         all( abs(result - expected) .lt. 1.E-6_real12 ), &
         'Set difference min zero failed', &
         success &
    )

    c = [1.0_real12, 2.0_real12, 3.0_real12, 4.0_real12]
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