module rw_vasprun
  use constants, only: real12
  implicit none

  private


contains
  recursive subroutine find_section(unit, section, found, endsection, depth)
    implicit none
    integer, intent(in) :: unit
    character(len=*), dimension(:), intent(in) :: section
    logical, intent(out) :: found
    character(len=*), intent(in), optional :: endsection
    integer, intent(in), optional :: depth

    character(len=100) :: line
    integer :: ierror, depth_
    character(len=:), allocatable :: &
         section_string, enclosing_section_end_string

    found = .false.

    if(present(depth)) then
       depth_ = depth
    else
       depth_ = 0
    end if
   
    !!! write depth_ number of spaces to the beginning of section_string
    !write(section_string, '(A)') repeat(' ', depth_), section(1)
    section_string = repeat(' ', depth_)//'<'//trim(section(1))//'>'
    if(present(endsection)) then
       enclosing_section_end_string = &
            repeat(' ', max(depth_-1,1))//'</'//trim(endsection)//'>'
    end if
    do
       read(unit, '(A)', iostat=ierror) line
       if(ierror /= 0) then
          write(0,*) 'Error reading vasprun.xml'
          stop
       end if
       if(index(line, trim(section_string)) .eq. 1) then
          found = .true.
          if(size(section) .gt. 1) then
             call find_section(unit, section(2:), found, section(1), depth_+1)
          end if
          exit
       elseif(present(endsection)) then
          if(index(line, trim(enclosing_section_end_string)) .eq. 1) then
             found = .false.
             exit
          end if
       end if
    end do

  end subroutine find_section

  function get_energy(unit, rewind_file) result(energy)
    implicit none
    integer, intent(in) :: unit
    real(real12) :: energy
    character(len=100) :: line
    logical, intent(in), optional :: rewind_file

    integer :: ierror
    logical :: found = .false.
    character(len=32), dimension(3) :: section_list

    if(present(rewind_file)) then
       if(rewind_file) rewind(unit)
    end if

    !!! read until inside <modeling>, but not encountered </modeling>
    !!! then read until <calculation> but not encountered </calculation>
    !!! then read until <energy> but not encountered </energy>
    !!! then read <i name="e_fr_energy"> VALUE </i> and extract VALUE as energy

    !!! make a recursive subroutine that checks for the open and close
    !!! pass in a list of strings to check for open and close
    !!! and then a value to return at the end
    !!! well, best to give it the part to, and then it exits the recursive
    section_list(1) = 'modeling'
    section_list(2) = 'calculation'
    section_list(3) = 'energy'
    call find_section(unit, section_list, found)
    !!! then, call a new procedure that gets the value from a value named <i name="NAME">
    !!! there should also be similar ones that handle vectors/arrays

    do
       read(unit, '(A)', iostat=ierror) line
       if(ierror .ne. 0) then
          write(0,*) 'Error reading energy from vasprun.xml'
          stop
       end if
       if(index(line, 'energy') /= 0) exit

    end do

  end function get_energy

  !!! for atomic structure:
  !!! read until inside <modeling>, but not encountered </modeling>
  !!! then read until <structure name="finalpos"> but not encountered </structure>
  !!! then read until <crystal> but not encountered </crystal>
  !!! THEN STOP AND READ THE CRYSTAL LATTICE
  !!! then read until <varray name="positions"> but not encountered </varray>
  !!! THEN STOP AND READ THE POSITIONS


  !!! BEFORE THAT:
  !!! read until inside <modeling>, but not encountered </modeling>
  !!! then read until <atominfo> but not encountered </atominfo>
  !!! then read <atoms> VALUE </atoms>
  !!! then read <types> VALUE </types>
  !!! then read until <array name="atoms" > but not encountered </array>
  !!! then read until <set> but not encountered </set>
  !!! then read <rc><c>NAME</c><c> NUMBER</c></rc>
  !!!  ... repeat previous line for VALUE (types) number of times
  !!! NAH that is each atom listed
  !!! instead, read until <array name="atomtypes" > but not encountered </array>
  !!! then read until <set> but not encountered </set>
  !!! then read <rc><c>NUMBER</c><c> NAME</c><c> MASS</c><c> VALENCY</c><c> PSEUDO_NAME</c></rc>

end module