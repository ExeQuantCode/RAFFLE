module raffle__io_utils
  !! Module for handling errors and io calls in the program.
  !!
  !! This module provides the expected procedure for stopping a program.
  !! If in testing mode, the stop can be suppressed.
  implicit none
  logical :: test_error_handling = .false.

  logical :: suppress_warnings = .false.
  character(len=*), parameter :: raffle__version__ = "1.0.1"

  private

  public :: raffle__version__
  public :: test_error_handling, suppress_warnings
  public :: stop_program, print_warning
  public :: print_version, print_build_info


contains

!###############################################################################
  subroutine stop_program(message, exit_code, block_stop)
    !! Stop the program and print an error message.
    implicit none
    character(len=*), intent(in) :: message
    integer, intent(in), optional :: exit_code
    logical, intent(in), optional :: block_stop

    integer :: exit_code_
    logical :: block_stop_

    if(present(exit_code)) then
       exit_code_ = exit_code
    else
       exit_code_ = 1
    end if
    if(present(block_stop)) then
       block_stop_ = block_stop
    else
       block_stop_ = .false.
    end if

    write(0,*) 'ERROR: ', trim(message)
    if(.not.block_stop_)then
       if(.not.test_error_handling) then
          stop exit_code_
       end if
    end if
  end subroutine stop_program
!###############################################################################


!###############################################################################
  subroutine print_warning(message)
    !! Print a warning message
    implicit none
    character(len=*), intent(in) :: message

    if(.not.suppress_warnings) then
       write(0,*) 'WARNING: ', trim(message)
    end if
  end subroutine print_warning
!###############################################################################


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
