module raffle__error_handling
  !! Module for handling errors in the program.
  !!
  !! This module provides the expected procedure for stopping a program.
  !! If in testing mode, the stop can be suppressed.
  implicit none
  logical :: test_error_handling = .false.
  
  private

  public :: test_error_handling
  public :: stop_program


contains

!###############################################################################
  module subroutine stop_program(message, exit_code)
    !! Stop the program and print an error message.
    implicit none
    character(len=*), intent(in) :: message
    integer, intent(in), optional :: exit_code

    integer :: exit_code_

    if(present(exit_code)) then
       exit_code_ = exit_code
    else
       exit_code_ = 1
    end if

    write(0,*) 'ERROR: ', trim(message)
    if(.not.test_error_handling)then
       stop exit_code_
    end if
  end subroutine stop_program
!###############################################################################

end module raffle__error_handling