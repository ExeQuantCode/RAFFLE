!!!#############################################################################
!!! Code written by Ned Thaddeus Taylor and Francis Huw Davies
!!! Code part of the ARTEMIS group (Hepplestone research group).
!!! Think Hepplestone, think HRG.
!!!#############################################################################
!!! module contains various miscellaneous maths functions and subroutines.
!!! module includes the following functionsand subroutines:
!!! lnsum            (sum a set of log(i), where i=1,n)
!!! triangular_number (returns the nth triangular number)
!!!#############################################################################
module misc_maths
  use constants, only: real12, pi
  implicit none


  private

  public :: lnsum, triangular_number


!!!updated 2020/02/03


contains

!!!#####################################################
!!! Sum of logs of range from 1 to n
!!!#####################################################
  double precision function lnsum(n) 
    implicit none
    integer :: i,n
    lnsum=0
    do i=1,n
       lnsum=lnsum+log(real(i))
    end do

    return
  end function lnsum
!!!#####################################################


!!!#####################################################
!!! returns the nth triangular number
!!!#####################################################
  integer function triangular_number(n)
    implicit none
    integer :: n
    triangular_number = n * ( n + 1 ) / 2
  end function triangular_number
!!!#####################################################

end module misc_maths
