module raffle__constants
  !! Module with global constants
  !!
  !! This module contains global constants that may be used throughout the
  !! library.
  implicit none
  integer, parameter, public :: real12 = Selected_real_kind(6,37)!(15,307)
  real(real12), parameter, public :: pi = 4._real12 * atan(1._real12)
  real(real12), parameter, public :: c = 0.26246582250210965422_real12
  real(real12), parameter, public :: c_vasp = 0.262465831_real12
  real(real12), parameter, public :: INF = huge(0._real12)
  complex(real12), parameter, public :: imag=(0._real12, 1._real12)
  integer, public :: verbose = 0
end module raffle__constants
