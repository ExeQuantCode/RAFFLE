module raffle__constants
  !! Module with global constants
  !!
  !! This module contains global constants that may be used throughout the
  !! library.
  implicit none
  integer, parameter, public :: real32 = Selected_real_kind(6,37)!(15,307)
  real(real32), parameter, public :: pi = 4._real32 * atan(1._real32)
  real(real32), parameter, public :: c = 0.26246582250210965422_real32
  real(real32), parameter, public :: c_vasp = 0.262465831_real32
  real(real32), parameter, public :: INF = huge(0._real32)
  complex(real32), parameter, public :: imag=(0._real32, 1._real32)
end module raffle__constants
