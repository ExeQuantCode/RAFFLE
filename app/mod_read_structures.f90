module read_structures
  !! Module for reading structures from various file formats.
  !!
  !! This module takes a list of directories and reads in the structures from
  !! the contained files. The structures are then converted to a set of
  !! generalised vectors (gvectors, aka distribution functions).
  !! NOTE: This does not currently handle relaxations in POSCAR/CONTCAR format.
  !! NOTE: This does not currently restart from a printed distribs_container.
  use raffle__constants, only: real32
  use raffle__io_utils, only: stop_program
  use raffle__misc, only: grep
  use raffle__misc_linalg, only: modu
  use raffle__geom_rw, only: basis_type, geom_read, geom_write, igeom_input
  use rw_vasprun, only: get_energy_from_vasprun, get_structure_from_vasprun
  use raffle__distribs_container, only: distribs_container_type
  implicit none


  private

  public :: get_gdfs_from_data


contains

!###############################################################################
  function get_gdfs_from_data(input_dir, &
       file_format, distribs_container_template &
  ) result(distribs_container)
    !! Read structures from the input directories and generate the generalised
    !! distribution functions.
    implicit none

    ! Arguments
    character(1024), dimension(..), intent(in) :: input_dir
    !! List of directories containing the structures to be read.
    type(distribs_container_type), intent(in), optional :: &
         distribs_container_template
    !! Optional. A template distribs_container to be used.
    type(distribs_container_type), allocatable :: distribs_container
    !! The distribs_container containing the distribution functions.
    character(*), intent(in), optional :: file_format
    !! Optional. The format of the input files. Default is vasprun.xml.

    ! Local variables
    integer :: i, j
    !! Loop indices.
    integer :: ifile_format
    !! The format of the input files.
    integer :: num_structures
    !! The number of structures read.
    real(real32) :: energy
    !! The energy of the structure.
    character(50) :: buffer
    !! A buffer for reading strings.
    character(256) :: format
    !! A format string for writing output.
    logical :: success
    !! Boolean for success of file operations.

    integer :: unit
    !! File units.
    integer :: iostat
    !! I/O status.
    type(basis_type) :: basis
    !! The basis of the structure.
    type(basis_type), dimension(:), allocatable :: basis_list
    !! The list of bases.
    character(256) :: stop_msg
    !! The stop message.
    character(256), dimension(:), allocatable :: structure_list
    !! The list of structure files.


    !---------------------------------------------------------------------------
    ! handle optional arguments
    !---------------------------------------------------------------------------
    if(present(distribs_container_template)) then
       distribs_container = distribs_container_template
    else
       distribs_container = distribs_container_type()
    end if

    if(present(file_format)) then
       select case(trim(adjustl(file_format)))
       case('vasprun.xml','xml','vasprun')
          ifile_format = 0
       case('POSCAR','OUTCAR')
          ifile_format = 1
          igeom_input = 1
       case('xyz','extxyz')
          ifile_format = 2
          igeom_input = 6
       case default
          write(stop_msg,'("Unknown file format: ",A)') file_format
          call stop_program(stop_msg)
          return
       end select
    else
       ifile_format = 0
    end if


    !---------------------------------------------------------------------------
    ! read the list of structure files from the input directories
    !---------------------------------------------------------------------------
    ! inquire(file='distribs_container.dat', exist=success)
    ! if(success) then
    !    call distribs_container%read('distribs_container.dat')
    !    goto 100
    ! end if
    select rank(dir_ptr => input_dir)
    rank(0)
       structure_list = [ get_structure_list( dir_ptr, ifile_format ) ]
    rank(1)
       do i = 1, size(dir_ptr)
          if(i.eq.1)then
             structure_list = [ &
                  get_structure_list( dir_ptr(i), ifile_format ) &
             ]
          else
             structure_list = [ &
                  structure_list, &
                  get_structure_list( dir_ptr(i), ifile_format ) &
             ]
          end if
       end do
    end select


    !---------------------------------------------------------------------------
    ! read the structures from the list of files
    !---------------------------------------------------------------------------
    num_structures = 0
    allocate(basis_list(0))
    do i = 1, size(structure_list)
       write(*,*) "Reading structure: ", trim(adjustl(structure_list(i)))
       select case(ifile_format)
       case(0) ! vasprun.xml
          open( &
               newunit=unit, &
               file=trim(adjustl(structure_list(i)))//"/vasprun.xml" &
          )
          write(*,*) "Reading structures from xml"
          basis%energy = get_energy_from_vasprun(unit, success)
          if(.not.success) cycle
          rewind(unit)
          call get_structure_from_vasprun(unit, basis, success)
          if(.not.success) cycle
          close(unit)
       case(1) ! POSCAR/OUTCAR
          open(newunit=unit, file=trim(adjustl(structure_list(i)))//"/POSCAR")
          write(*,*) "Reading structures from POSCAR"
          call geom_read(unit, basis)
          close(unit)
          open(newunit=unit, file=trim(adjustl(structure_list(i)))//"/OUTCAR")
          call grep( &
               unit, 'free  energy   TOTEN  =', &
               lline=.false., success=success &
          )
          if(.not.success) cycle
          backspace(unit)
          read(unit,*) buffer, buffer, buffer, buffer, energy
          close(unit)
          basis%energy = energy
       case(2) ! extended xyz
          open(newunit=unit, file=trim(adjustl(structure_list(i))))
          write(*,*) "Reading structures from xyz"
          do
             read(unit,'(A)',iostat=iostat) buffer
             if(iostat .ne. 0) exit
             if(trim(buffer).eq."") cycle
             backspace(unit)
             call geom_read(unit, basis)

             num_structures = num_structures + 1
             write(format,'&
                  &("(""Found structure: "",I4,"" with energy: "",&
                  &F13.7,"" and formula:"",",I0,"(1X,A,I0))")' &
             ) basis%nspec
             write(*,format) &
                  num_structures, &
                  basis%energy, &
                  ( &
                       trim(basis%spec(j)%name), basis%spec(j)%num, &
                       j=1, basis%nspec &
                  )
             basis_list = [ basis_list, basis ]
          end do
          cycle
       end select

       write(*,*) &
            "Found structure: ", trim(adjustl(structure_list(i))), &
            " with energy: ", basis%energy
       basis_list = [ basis_list, basis ]
       num_structures = num_structures + 1
    end do


    !---------------------------------------------------------------------------
    ! generate the generalised distribution functions
    !---------------------------------------------------------------------------
    call distribs_container%create(basis_list)

    igeom_input = 1

  end function get_gdfs_from_data
!###############################################################################


!###############################################################################
  function get_structure_list(input_dir, ifile_format) result(structure_list)
    !! Get a list of structures from a directory.
    implicit none

    ! Arguments
    character(*), intent(in) :: input_dir
    !! The directory containing the structures.
    integer, intent(in) :: ifile_format
    !! The format of the input files.
    character(256), dimension(:), allocatable :: structure_list
    !! The list of structure files.

    ! Local variables
    integer :: unit
    !! File unit.
    integer :: iostat
    !! I/O status.
    character(256) :: name
    !! The name of the structure file.
    character(256) :: stop_msg
    !! The stop message.
    logical :: file_exists, addit_file_exists, structures_found
    !! Booleans for file existence.


    !---------------------------------------------------------------------------
    ! get the list of structure files using the command line
    !---------------------------------------------------------------------------
    structures_found = .false.
    call execute_command_line( &
         "ls "//trim(adjustl(input_dir))//"/. >structures.txt", wait = .TRUE. )
    open(newunit=unit, file="structures.txt", status="old")
    do
       read(unit,'(A)',iostat=iostat) name
       if(is_iostat_end(iostat)) exit
       name = trim(adjustl(input_dir))//"/"//trim(adjustl(name))
       select case(ifile_format)
       case(0)
          inquire(file=trim(name)//"/vasprun.xml", exist=file_exists)
       case(1)
          inquire(file=trim(name)//"/POSCAR", exist=file_exists)
          inquire(file=trim(name)//"/OUTCAR", exist=addit_file_exists)
          file_exists = file_exists.and.addit_file_exists
       case(2)
          inquire(file=trim(name), exist=file_exists)
          file_exists = file_exists.and.(index(name, ".xyz") .gt. 0)
          if(.not.file_exists)then
             inquire(file=trim(name)//"/data.xyz", exist=file_exists)
             if(file_exists) name = trim(name)//"/data.xyz"
          end if
       end select
       if(.not.file_exists) cycle
       if(structures_found) then
          structure_list = [ structure_list, name ]
       else
          structure_list = [ name ]
          structures_found = .true.
       end if
    end do
    close(unit)

    if(.not.structures_found) then
       write(stop_msg,'("No structures found in directory: """,A, &
            &""" with format: ",I0)' &
       ) trim(input_dir), ifile_format
       call stop_program(stop_msg)
       return
    end if

  end function get_structure_list
!###############################################################################

end module read_structures