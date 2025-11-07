module raffle__io_utils
  !! Module for handling errors and io calls in the program.
  !!
  !! This module provides the expected procedure for stopping a program.
  !! If in testing mode, the stop can be suppressed.
  implicit none

  character(len=*), parameter :: raffle__version__ = "1.1.1"

  private

  public :: raffle__version__
  public :: print_version, print_build_info


contains

!###############################################################################
  subroutine print_version()
    !! Print the version number of the program.
    implicit none

    write(*,'("version: ",A)') raffle__version__
  end subroutine print_version
!###############################################################################


!###############################################################################
  subroutine print_build_info()
    !! Print the build information of the program.
    implicit none

    write(*,'("RAFFLE: pseudoRandom Approach For Finding Local Energy minima")')
    call print_version()
    write(*,'(" (build ",A,1X,A,")")') __DATE__, __TIME__

  end subroutine print_build_info
!###############################################################################

end module raffle__io_utils
