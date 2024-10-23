module raffle__misc_maths
  !! Module for miscellaneous mathematical functions.
  use raffle__io_utils, only: stop_program
  use raffle__constants, only: real32
  implicit none


  private

  public :: lnsum, triangular_number, set_difference



contains

!###############################################################################
  function lnsum(n)
    !! Return the sum of the logs of the integers from 1 to n.
    implicit none

    ! Arguments
    integer :: n
    !! The upper limit of the range.
    real(real32) :: lnsum
    !! The sum of the logs of the integers from 1 to n.

    ! Local variables
    integer :: i
    !! Loop index.

    lnsum = 0._real32
    do i = 1, n
       lnsum = lnsum + log( real(i, real32) )
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


!###############################################################################
  function set_difference(a, b, set_min_zero)
    !! Return the set difference of two arrays.
    implicit none

    ! Arguments
    real(real32), dimension(:), intent(in) :: a
    !! The first array.
    real(real32), dimension(:), intent(in) :: b
    !! The second array.
    logical, optional :: set_min_zero
    !! Boolean to set the maximum value of the output array to zero.
    real(real32), dimension(size(a)) :: set_difference
    !! The set difference of the two arrays.

    ! Local variables
    integer :: i
    !! Loop indices.
    logical :: set_min_zero_
    !! Boolean to set all values below zero to zero.


    if(present(set_min_zero)) then
       set_min_zero_ = set_min_zero
    else
       set_min_zero_ = .false.
    end if

    if(size(a,1) .ne. size(b,1)) then
       call stop_program('Arrays must be the same size.')
       return
    end if

    if(set_min_zero_)then
       do i = 1, size(a,1)
          set_difference(i) = max(0.0_real32, a(i) - b(i))
       end do
    else
       set_difference = a - b
    end if

  end function set_difference
!###############################################################################

end module raffle__misc_maths
