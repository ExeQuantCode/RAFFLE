module rw_vasprun
  !! Module for reading VASP vasprun.xml files.
  !!
  !! This module provides the procedures for reading the energy and atomic
  !! structure from a VASP vasprun.xml file.
  use raffle__constants, only: real32
  use raffle__io_utils, only: stop_program
  use raffle__geom_rw, only: basis_type
  implicit none


  private

  public :: get_energy_from_vasprun
  public :: get_structure_from_vasprun


contains

!###############################################################################
  recursive subroutine find_section( &
       unit, section, found, name, end_section, depth &
  )
    !! Find a section in a vasprun.xml file.
    !!
    !! This subroutine reads a vasprun.xml file and searches for an embedded
    !! section. The parent sections are the earlier elements of the section
    !! array.
    implicit none

    ! Arguments
    integer, intent(in) :: unit
    !! The unit number of the file.
    character(len=*), dimension(:), intent(in) :: section
    !! The section to find.
    logical, intent(out) :: found
    !! Whether the section was found.
    character(len=*), intent(in), optional :: end_section
    !! The end section to find.
    character(len=*), dimension(:), intent(in), optional :: name
    !! The optional name of the section.
    integer, intent(in), optional :: depth
    !! The current depth of the recursion.

    ! Local variables
    character(len=100) :: line
    !! A buffer for reading lines.
    integer :: ierror, depth_
    !! I/O status, depth.
    character(len=:), allocatable :: &
         section_string, enclosing_section_end_string, name_string
    !! Strings for the section and end section.
    character(len=:), dimension(:), allocatable :: name_
    !! Array of names.


    found = .false.

    !---------------------------------------------------------------------------
    ! handle optional arguments
    !---------------------------------------------------------------------------
    if(present(depth)) then
       depth_ = depth
    else
       depth_ = 0
    end if
    if(present(name)) then
       if(size(name) .ne. size(section)) then
          call stop_program( &
               'name and section must be same size in find_section' &
          )
          return
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


    !---------------------------------------------------------------------------
    ! write the indentation to the section string, based on the depth
    !---------------------------------------------------------------------------
    section_string = &
         repeat(' ', depth_)//'<'//trim(section(1))//trim(name_string)
    if(present(end_section)) then
       enclosing_section_end_string = &
            repeat(' ', max(depth_-1,1))//'</'//trim(end_section)//'>'
    end if


    !---------------------------------------------------------------------------
    ! read the file and search for the section
    !---------------------------------------------------------------------------
    do
       read(unit, '(A)', iostat=ierror) line
       if(is_iostat_end(ierror)) exit
       if(ierror .ne. 0) then
          call stop_program( &
               'Issue encountered when reading vasprun.xml' &
          )
          return
       end if
       if(index(line, trim(section_string)) .eq. 1) then
          found = .true.
          if(size(section) .gt. 1) then
             call find_section(unit, [ section(2:) ], found, &
                  end_section = section(1), &
                  depth = depth_+1, &
                  name = [ name_(2:) ] &
             )
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
!###############################################################################


!###############################################################################
  function get_energy_from_vasprun(unit, found, rewind_file) result(energy)
    !! Get the energy from a vasprun.xml file.
    implicit none

    ! Arguments
    integer, intent(in) :: unit
    !! The unit number of the file.
    logical, intent(out) :: found
    !! Whether the energy was found.
    logical, intent(in), optional :: rewind_file
    !! Whether to rewind the file after reading the energy.
    real(real32) :: energy
    !! The energy of the structure.

    ! Local variables
    integer :: ierror
    !! I/O status.
    logical :: found_ = .false.
    !! Whether the section was found.
    character(len=100) :: line, buffer
    !! A buffer for reading lines.
    real(real32), dimension(:), allocatable :: energy_list
    !! List of energies.
    character(len=32), dimension(3) :: section_list
    !! List of sections.


    found = .false.

    !---------------------------------------------------------------------------
    ! handle optional arguments
    !---------------------------------------------------------------------------
    if(present(rewind_file)) then
       if(rewind_file) rewind(unit)
    end if


    !---------------------------------------------------------------------------
    ! set up and find embedded sections
    !---------------------------------------------------------------------------
    section_list(1) = 'modeling'
    section_list(2) = 'calculation'
    section_list(3) = 'energy'
    call find_section(unit, section_list(1:1), found_)


    !---------------------------------------------------------------------------
    ! read the energy from the file
    !---------------------------------------------------------------------------
    allocate(energy_list(0))
    do
       call find_section(unit, section_list(2:), found_, depth=1)
       if (.not. found_) exit
       read(unit, '(A)', iostat=ierror) line
       if(ierror .ne. 0) then
          call stop_program( &
               'Issue encountered when reading energy from vasprun.xml' &
          )
          return
       end if
       read(line, '(3X, A22, F16.8)') buffer, energy
       energy_list = [ energy_list, energy ]
       found = .true.
    end do

  end function get_energy_from_vasprun
!###############################################################################


!###############################################################################
  subroutine get_structure_from_vasprun(unit, basis, found)
    !! Get the atomic structure from a vasprun.xml file.
    implicit none

    ! Arguments
    integer, intent(in) :: unit
    !! The unit number of the file.
    logical, intent(out), optional :: found
    !! Whether the structure was found.
    type(basis_type) :: basis
    !! The basis of the structure.

    ! Local variables
    integer :: ierror, i, is, ia
    !! I/O status, loop indices.
    logical :: found_ = .false.
    !! Whether the section was found.
    character(len=100) :: line, buffer
    !! A buffer for reading lines.
    character(len=32), dimension(4) :: section_list, name_list
    !! List of sections and names.

    integer :: number
    !! The number of atoms of a given type.
    character(len=3) :: element
    !! The element symbol.
    real(real32) :: mass, valency
    !! The mass and valency of the element.
    character(len=40) :: pseudo
    !! The pseudopotential file name.
    integer, dimension(:), allocatable :: number_list
    !! List of numbers of atoms.
    character(len=3), dimension(:), allocatable :: element_list
    !! List of element symbols.
    real(real32), dimension(:), allocatable :: mass_list, valency_list
    !! List of masses and valencies.
    character(len=40), dimension(:), allocatable :: pseudo_list
    !! List of pseudopotential file names.
    
    
    found = .false.

    !---------------------------------------------------------------------------
    ! set up and find embedded sections for the incar
    !---------------------------------------------------------------------------
    section_list(1) = 'modeling'
    section_list(2) = 'incar'
    call find_section(unit, section_list(:2), found_)
    if(.not. found_) then
       call stop_program('incar section not found in vasprun.xml')
       return
    end if
    read(unit, '(A)', iostat=ierror) line
    if(ierror .ne. 0) then
       call stop_program( &
            'Error encountered when reading incar from vasprun.xml' &
       )
       return
    end if
    read(line, '(2X,A31,A)') buffer, basis%sysname
    basis%sysname = basis%sysname(:index(basis%sysname, '<')-1)


    !---------------------------------------------------------------------------
    ! set up and find embedded sections for the atomtypes
    !---------------------------------------------------------------------------
    section_list(1) = 'modeling'
    section_list(2) = 'atominfo'
    section_list(3) = 'array'
    section_list(4) = 'set'

    name_list(1) = ""
    name_list(2) = ""
    name_list(3) = "atomtypes"
    name_list(4) = ""
    call find_section( &
         unit, section_list(2:), found_, name=name_list(2:), depth=1 &
    )
    if(.not. found_) then
       call stop_program('Section "set" not found in vasprun.xml')
       return
    end if

    !---------------------------------------------------------------------------
    ! read atomtypes data from file
    !---------------------------------------------------------------------------
    i = 0
    do
       i = i + 1
       read(unit, '(A)', iostat=ierror) line
       if(ierror .ne. 0) then
          call stop_program( &
               'Error encountered when reading atomtypes from vasprun.xml' &
          )
          return
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

    !---------------------------------------------------------------------------
    ! allocate the basis from the atomtypes data
    !---------------------------------------------------------------------------
    basis%natom = sum(number_list)
    basis%nspec = size(element_list)
    allocate(basis%spec(basis%nspec))
    basis%spec(:)%name = element_list
    basis%spec(:)%num  = number_list
    basis%spec(:)%mass = mass_list
    do is = 1, basis%nspec
       allocate(basis%spec(is)%atom(basis%spec(is)%num,3))
    end do


    !---------------------------------------------------------------------------
    ! set up and find embedded sections for the lattice
    !---------------------------------------------------------------------------
    section_list(2) = 'structure'
    section_list(3) = 'crystal'
    section_list(4) = 'varray'
    name_list(2) = "finalpos"
    name_list(3) = ""
    name_list(4) = "basis"
    call find_section( &
         unit, section_list(2:), found_, name=name_list(2:), depth=1 &
    )
    if(.not. found_) then
       call stop_program( &
            'Error encountered when finding structure in vasprun.xml' &
       )
       return
    end if
    do i = 1, 3
       read(unit, '(A)', iostat=ierror) line
       if(ierror .ne. 0) then
          call stop_program( &
               'Error encountered when reading lattice from vasprun.xml' &
          )
          return
       end if
       read( line, '(4X,A3,3(1X,F16.8))' ) buffer, basis%lat(i,:)
    end do


    !---------------------------------------------------------------------------
    ! set up and find embedded sections for the atomic structure
    !---------------------------------------------------------------------------
    section_list(3) = 'varray'
    name_list(3) = "positions"
    call find_section( &
         unit, section_list(3:3), found_, name=name_list(3:3), depth=2 &
    )
    if(.not. found_) then
       call stop_program( &
            'Error encountered when finding structure in vasprun.xml' &
       )
       return
    end if
    ia = 0
    is = 1
    do i = 1, basis%natom
       read(unit, '(A)', iostat=ierror) line
       if(ierror .ne. 0) then
          call stop_program( &
               'Error encountered when reading positions from vasprun.xml' &
          )
          return
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
!###############################################################################

end module