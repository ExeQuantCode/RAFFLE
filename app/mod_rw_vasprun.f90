module rw_vasprun
  use constants, only: real12
  use rw_geom, only: basis_type
  implicit none

  private

  public :: get_energy_from_vasprun
  public :: get_structure_from_vasprun


contains
  recursive subroutine find_section( &
       unit, section, found, name, end_section, depth)
    implicit none
    integer, intent(in) :: unit
    character(len=*), dimension(:), intent(in) :: section
    logical, intent(out) :: found
    character(len=*), intent(in), optional :: end_section
    character(len=*), dimension(:), intent(in), optional :: name
    integer, intent(in), optional :: depth

    character(len=100) :: line
    integer :: ierror, depth_
    character(len=:), allocatable :: &
         section_string, enclosing_section_end_string, name_string
    character(len=:), dimension(:), allocatable :: name_

    found = .false.

    if(present(depth)) then
       depth_ = depth
    else
       depth_ = 0
    end if
    if(present(name)) then
       if(size(name) .ne. size(section)) then
          write(0,*) 'Error in find_section: name and section must be same size'
          stop
       end if
       name_ = name
    else
       allocate(character(len=1) :: name_(size(section)))
       name_ = ""
    end if
    if(trim(name_(1)).ne."")then
       name_string = ' name="'//trim(adjustl(name_(1)))//'" >'
    else
      name_string = ">"
    end if

    !!! write depth_ number of spaces to the beginning of section_string
    section_string = repeat(' ', depth_)//'<'//trim(section(1))//trim(name_string)
    if(present(end_section)) then
       enclosing_section_end_string = &
            repeat(' ', max(depth_-1,1))//'</'//trim(end_section)//'>'
    end if
    do
       read(unit, '(A)', iostat=ierror) line
       if(is_iostat_end(ierror)) exit
       if(ierror .ne. 0) then
          write(0,*) 'Error reading vasprun.xml'
          stop
       end if
       if(index(line, trim(section_string)) .eq. 1) then
          found = .true.
          if(size(section) .gt. 1) then
             call find_section(unit, [ section(2:) ], found, &
                  end_section = section(1), depth = depth_+1, name = [ name_(2:) ] )
          end if
          exit
       elseif(present(end_section)) then
          if(index(line, trim(enclosing_section_end_string)) .eq. 1) then
             found = .false.
             exit
          end if
       end if
    end do

  end subroutine find_section

  function get_energy_from_vasprun(unit, found, rewind_file) result(energy)
    implicit none
    integer, intent(in) :: unit
    real(real12) :: energy
    character(len=100) :: line, buffer
    logical, intent(out) :: found
    logical, intent(in), optional :: rewind_file

    integer :: ierror
    logical :: found_ = .false.
    real(real12), dimension(:), allocatable :: energy_list
    character(len=32), dimension(3) :: section_list


    found = .false.
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
    call find_section(unit, section_list(1:1), found_)
    !!! then, call a new procedure that gets the value from a value named <i name="NAME">
    !!! there should also be similar ones that handle vectors/arrays

    allocate(energy_list(0))
    do
       call find_section(unit, section_list(2:), found_, depth=1)
       if (.not. found_) exit
       read(unit, '(A)', iostat=ierror) line
       if(ierror .ne. 0) then
          write(0,*) 'Error reading energy from vasprun.xml'
          stop
       end if
       read(line, '(3X, A22, F16.8)') buffer, energy
       energy_list = [ energy_list, energy ]
       found = .true.
    end do

  end function get_energy_from_vasprun

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

  subroutine get_structure_from_vasprun(unit, basis, found)
    implicit none
    integer, intent(in) :: unit
    type(basis_type) :: basis

    integer :: ierror, i, is, ia
    logical, intent(out), optional :: found
    logical :: found_ = .false.
    character(len=100) :: line, buffer
    character(len=32), dimension(4) :: section_list, name_list

    integer :: number
    character(len=3) :: element
    real(real12) :: mass, valency
    character(len=40) :: pseudo
    integer, dimension(:), allocatable :: number_list
    character(len=3), dimension(:), allocatable :: element_list
    real(real12), dimension(:), allocatable :: mass_list, valency_list
    character(len=40), dimension(:), allocatable :: pseudo_list
    
    
    found = .false.
    section_list(1) = 'modeling'
    section_list(2) = 'incar'
    call find_section(unit, section_list(:2), found_)
    if(.not. found_) then
       write(0,*) 'Error finding incar in vasprun.xml'
       stop
    end if
    read(unit, '(A)', iostat=ierror) line
    if(ierror .ne. 0) then
       write(0,*) 'Error reading incar from vasprun.xml'
       stop
    end if
    read(line, '(2X,A31,A)') buffer, basis%sysname
    basis%sysname = basis%sysname(:index(basis%sysname, '<')-1)


    section_list(1) = 'modeling'
    section_list(2) = 'atominfo'
    section_list(3) = 'array'
    section_list(4) = 'set'

    name_list(1) = ""
    name_list(2) = ""
    name_list(3) = "atomtypes"
    name_list(4) = ""

    call find_section(unit, section_list(2:), found_, name=name_list(2:), depth=1)
    if(.not. found_) then
       write(0,*) 'Error finding set in vasprun.xml'
       stop
    end if

    i = 0
    do
       i = i + 1
       read(unit, '(A)', iostat=ierror) line
       if(ierror .ne. 0) then
          write(0,*) 'Error reading atomtypes from vasprun.xml'
          stop
       end if
       if(index(line, '   </set>') .eq. 1) exit
       read( line, '(4X,A7, I4, A7, A2, A7, F16.8, A7, F16.8, A7, A40)' ) &
            buffer, number, &
            buffer, element, &
            buffer, mass, &
            buffer, valency, &
            buffer, pseudo
       if(i .eq. 1) then
          number_list  = [ number ]
          element_list = [ element ]
          mass_list    = [ mass ]
          valency_list = [ valency ]
          pseudo_list  = [ pseudo ]
       else
          number_list  = [ number_list, number ]
          element_list = [ element_list, element ]
          mass_list    = [ mass_list, mass ]
          valency_list = [ valency_list, valency ]
          pseudo_list  = [ pseudo_list, pseudo ]
       end if
    end do
    basis%natom = sum(number_list)
    basis%nspec = size(element_list)
    allocate(basis%spec(basis%nspec))
    basis%spec(:)%name = element_list
    basis%spec(:)%num  = number_list
    basis%spec(:)%mass = mass_list
    do is = 1, basis%nspec
       allocate(basis%spec(is)%atom(basis%spec(is)%num,3))
    end do



    section_list(2) = 'structure'
    section_list(3) = 'crystal'
    section_list(4) = 'varray'
    name_list(2) = "finalpos"
    name_list(3) = ""
    name_list(4) = "basis"

    call find_section(unit, section_list(2:), found_, name=name_list(2:), depth=1)
    if(.not. found_) then
       write(0,*) 'Error finding structure in vasprun.xml'
       stop
    end if
    do i = 1, 3
       read(unit, '(A)', iostat=ierror) line
       if(ierror .ne. 0) then
          write(0,*) 'Error reading lattice from vasprun.xml'
          stop
       end if
       read( line, '(4X,A3,3(1X,F16.8))' ) buffer, basis%lat(i,:)
    end do


    section_list(3) = 'varray'
    name_list(3) = "positions"
    call find_section(unit, section_list(3:3), found_, name=name_list(3:3), depth=2)
    if(.not. found_) then
       write(0,*) 'Error finding structure in vasprun.xml'
       stop
    end if
    ia = 0
    is = 1
    do i = 1, basis%natom
       read(unit, '(A)', iostat=ierror) line
       if(ierror .ne. 0) then
          write(0,*) 'Error reading positions from vasprun.xml'
          stop
       end if
       ia = ia + 1
       if(ia .gt. basis%spec(is)%num) then
          ia = 1
          is = is + 1
       end if
       read( line, '(3X,A3,3(1X,F16.8))' ) buffer, &
            basis%spec(is)%atom(ia,:3)
    end do
    found = .true.

  end subroutine get_structure_from_vasprun

end module