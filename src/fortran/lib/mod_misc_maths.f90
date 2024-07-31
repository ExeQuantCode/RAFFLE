module misc_maths
  !! Module for miscellaneous mathematical functions.
  implicit none


  private

  public :: lnsum, triangular_number



contains

!###############################################################################
double precision function lnsum(n) 
    !! Return the sum of the logs of the integers from 1 to n.
    implicit none

    ! Arguments
    integer :: n
    !! The upper limit of the range.

    ! Local variables
    integer :: i
    !! Loop index.

    lnsum=0
    do i=1,n
       lnsum=lnsum+log(real(i))
    end do

    return
  end function lnsum
!###############################################################################


!###############################################################################
  integer function triangular_number(n)
    !! Return the nth triangular number.
    implicit none

    ! Arguments
    integer :: n
    !! The index of the triangular number to return.
    triangular_number = n * ( n + 1 ) / 2
  end function triangular_number
!###############################################################################

end module misc_maths
