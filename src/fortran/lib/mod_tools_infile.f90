module raffle__tools_infile
  !! This module contains a collection of tools for reading input files.
  !!
  !! Code written by Ned Thaddeus Taylor and Francis Huw Davies
  !! Code part of the ARTEMIS group (Hepplestone research group).
  !! Think Hepplestone, think HRG.
  !! Original distribution: https://github.com/ExeQuantCode/ARTEMIS
  !! This module is distributed under the CC-BY-3.0 license.
  !! License: http://creativecommons.org/licenses/by/3.0/
  !! This module has been copied and modified with permission from the
  !! original authors.
  use raffle__constants, only: real32
  use raffle__misc, only: grep,icount

  implicit none


  private
  public :: getline, rm_comments
  public :: assign_val, assign_vec


  interface assign_val
     procedure assignI, assignR, assignS, assignL
  end interface assign_val

  interface assign_vec
     procedure assignIvec, assignRvec
  end interface assign_vec



contains

!###############################################################################
  function val(buffer)
    !! Return the section of buffer that occurs after an "=".
    implicit none

    ! Arguments
    character(*), intent(in) :: buffer
    !! The input buffer.

    ! Local variables
    character(100) :: val
    !! The output value.

    val = trim( adjustl( buffer((scan(buffer,"=",back=.false.)+1):) ) )

  end function val
!###############################################################################


!###############################################################################
  subroutine getline(unit, pattern, buffer)
    !! Get the line from a grep and assign it to buffer.
    implicit none

    ! Arguments
    integer, intent(in) :: unit
    !! The unit to read from.
    character(*), intent(in) :: pattern
    !! The pattern to grep for.
    character(*), intent(out) :: buffer
    !! The buffer to assign the line to.

    ! Local variables
    integer :: iostat
    !! input output status

    call grep(unit,pattern)
    backspace(unit)
    read(unit,'(A)',iostat=iostat) buffer

  end subroutine getline
!###############################################################################


!###############################################################################
  subroutine assignI(buffer, variable, found, keyword)
    !! Assign an integer to a variable.
    implicit none

    ! Arguments
    integer, intent(inout) :: found
    !! The number of variables found. External counter
    character(*), intent(inout) :: buffer
    !! The buffer to read from.
    integer, intent(out) :: variable
    !! The variable to assign to.
    character(*), optional, intent(in) :: keyword
    !! The keyword to search for.

    character(1024) :: buffer2

    if(present(keyword))then
       buffer = buffer(index(buffer,keyword):)
    end if
    if(scan("=",buffer).ne.0) buffer2 = val(buffer)
    if(trim(adjustl(buffer2)).ne.'') then
       found = found + 1
       read(buffer2,*) variable
    end if
  end subroutine assignI
!###############################################################################


!###############################################################################
  subroutine assignIvec(buffer, variable, found, keyword)
    !! Assign an arbitrary length vector of integers to a variable.
    implicit none

    ! Arguments
    integer, intent(inout) :: found
    !! The number of variables found. External counter
    character(*), intent(inout) :: buffer
    !! The buffer to read from.
    integer, dimension(:) :: variable
    !! The variable to assign to.
    character(*), optional, intent(in) :: keyword
    !! The keyword to search for.

    ! Local variables
    integer :: i
    !! Loop counter
    character(1024) :: buffer2
    !! Temporary buffer


    if(present(keyword))then
       buffer = buffer(index(buffer,keyword):)
    end if
    if(scan("=",buffer).ne.0) buffer2 = val(buffer)
    if(trim(adjustl(buffer2)).ne.'') then
       found = found + 1
       if(icount(buffer2).eq.1.and.&
            icount(buffer2).ne.size(variable))then
          read(buffer2,*) variable(1)
          variable = variable(1)
       else
          read(buffer2,*) ( variable(i), i = 1, size(variable) )
       end if
    end if
  end subroutine assignIvec
!###############################################################################


!###############################################################################
  subroutine assignR(buffer, variable, found, keyword)
    !! Assign a float value to a variable.
    implicit none

    ! Arguments
    integer, intent(inout) :: found
    !! The number of variables found. External counter
    character(*), intent(inout) :: buffer
    !! The buffer to read from.
    real(real32), intent(out) :: variable
    !! The variable to assign to.
    character(*), optional, intent(in) :: keyword
    !! The keyword to search for.

    ! Local variables
    character(1024) :: buffer2
    !! Temporary buffer

    if(present(keyword))then
       buffer = buffer(index(buffer,keyword):)
    end if
    if(scan("=",buffer).ne.0) buffer2 = val(buffer)
    if(trim(adjustl(buffer2)).ne.'') then
       found = found + 1
       read(buffer2,*) variable
    end if
  end subroutine assignR
!###############################################################################


!###############################################################################
  subroutine assignRvec(buffer, variable, found, keyword)
    !! Assign an arbitrary length float vector to a variable.
    implicit none

    ! Arguments
    integer, intent(inout) :: found
    !! The number of variables found. External counter
    character(*), intent(inout) :: buffer
    !! The buffer to read from.
    real(real32), dimension(:), intent(out) :: variable
    !! The variable to assign to.
    character(*), optional, intent(in) :: keyword
    !! The keyword to search for.

    ! Local variables
    integer :: i
    !! Loop counter
    character(1024) :: buffer2
    !! Temporary buffer

    if(present(keyword))then
       buffer = buffer(index(buffer,keyword):)
    end if
    if(scan("=",buffer).ne.0) buffer2=val(buffer)
    if(trim(adjustl(buffer2)).ne.'') then
       found = found + 1
       if(icount(buffer2).eq.1.and.&
            icount(buffer2).ne.size(variable))then
          read(buffer2,*) variable(1)
          variable = variable(1)
       else
          read(buffer2,*) (variable(i),i=1,size(variable))
       end if
    end if
  end subroutine assignRvec
!###############################################################################


!###############################################################################
  subroutine assignS(buffer, variable, found, keyword)
    !! Assign a string to a variable.
    implicit none

    ! Arguments
    integer, intent(inout) :: found
    !! The number of variables found. External counter
    character(*), intent(inout) :: buffer
    !! The buffer to read from.
    character(*), intent(out) :: variable
    !! The variable to assign to.
    character(*), optional, intent(in) :: keyword
    !! The keyword to search for.

    ! Local variables
    character(1024) :: buffer2
    !! Temporary buffer

    if(present(keyword))then
       buffer = buffer(index(buffer,keyword):)
    end if
    if(scan("=",buffer).ne.0) buffer2 = val(buffer)
    if(trim(adjustl(buffer2)).ne.'') then
       found = found + 1
       read(buffer2,'(A)') variable
    end if
  end subroutine assignS
!###############################################################################


!###############################################################################
  subroutine assignL(buffer, variable, found, keyword)
    !! Assign a logical to a variable.
    !!
    !! This subroutine will assign a logical value to a variable. The
    !! logical can take the form of a string or an integer. The following
    !! are all valid logical values:
    !! T, F, t, f, 1, 0
    implicit none

    ! Arguments
    integer, intent(inout) :: found
    !! The number of variables found. External counter
    character(*), intent(inout) :: buffer
    !! The buffer to read from.
    logical, intent(out) :: variable
    !! The variable to assign to.
    character(*), optional, intent(in) :: keyword
    !! The keyword to search for.

    ! Local variables
    character(1024) :: buffer2
    !! Buffer to read from

    if(present(keyword))then
       buffer=buffer(index(buffer,keyword):)
    end if
    if(scan("=",buffer).ne.0) buffer2 = val(buffer)
    if(trim(adjustl(buffer2)).ne.'') then
       found=found+1
       if(index(buffer2,"T").ne.0.or.&
            index(buffer2,"t").ne.0.or.&
            index(buffer2,"1").ne.0) then
          variable = .TRUE.
       end if
       if(index(buffer2,"F").ne.0.or.&
            index(buffer2,"f").ne.0.or.&
            index(buffer2,"0").ne.0) then
          variable = .FALSE.
       end if
    end if
  end subroutine assignL
!###############################################################################


!###############################################################################
  subroutine rm_comments(buffer, iline)
    !! Remove comments from a buffer.
    implicit none

    ! Arguments
    character(*), intent(inout) :: buffer
    !! Buffer to remove comments from.
    integer, optional, intent(in) :: iline
    !! Line number.

    ! Local variables
    integer :: lbracket,rbracket
    !! Bracketing variables
    integer :: iline_
    !! Line number

    iline_ = 0
    if(present(iline)) iline_ = iline

    if(scan(buffer,'!').ne.0) buffer = buffer(:(scan(buffer,'!')-1))
    if(scan(buffer,'#').ne.0) buffer = buffer(:(scan(buffer,'#')-1))
    do while(scan(buffer,'(').ne.0.or.scan(buffer,')').ne.0)
       lbracket = scan( buffer, '(', back=.true. )
       rbracket = scan( buffer(lbracket:), ')' )
       if( lbracket .eq. 0 .or. rbracket .eq. 0 )then
          write(6,'(A,I0)') &
               ' NOTE: a bracketing error was encountered on line ',iline_
          buffer = ""
          return
       end if
       rbracket = rbracket + lbracket - 1
       buffer = buffer(:(lbracket-1)) // buffer((rbracket+1):)
    end do

  end subroutine rm_comments
!###############################################################################

end module raffle__tools_infile
