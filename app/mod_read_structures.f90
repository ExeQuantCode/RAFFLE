module read_structures
  !! Module for reading structures from various file formats.
  !!
  !! This module takes a list of directories and reads in the structures from
  !! the contained files. The structures are then converted to a set of
  !! generalised vectors (gvectors, aka distribution functions).
  use constants, only: real12
  use misc_raffle, only: grep
  use misc_linalg, only: modu
  use rw_geom, only: basis_type, geom_read, geom_write, igeom_input
  use rw_vasprun, only: get_energy_from_vasprun, get_structure_from_vasprun
  use evolver, only: gvector_container_type
#ifdef ENABLE_ATHENA
  use machine_learning, only: network_setup, &
       network_train, network_train_graph, &
       network_predict, network_predict_graph
  use athena, only: graph_type
#endif
  implicit none


  private

  public :: get_evolved_gvectors_from_data


contains

!###############################################################################
  function get_evolved_gvectors_from_data(input_dir, &
       file_format, gvector_container_template) &
       result(gvector_container)
    !! Read structures from the input directories and evolve them to gvectors.
    implicit none

    ! Arguments
    character(*), dimension(..), intent(in) :: input_dir
    !! List of directories containing the structures to be read.
    type(gvector_container_type), intent(in), optional :: &
         gvector_container_template
    !! Optional. A template gvector_container to be used.
    type(gvector_container_type), allocatable :: gvector_container
    !! The gvector_container containing the evolved gvectors.
    character(*), intent(in), optional :: file_format
    !! Optional. The format of the input files. Default is vasprun.xml.

    ! Local variables
    character(256) :: name
    !! The name of the structure file.
    integer :: i, j
    !! Loop indices.
    integer :: ifile_format
    !! The format of the input files.
    integer :: num_structures
    !! The number of structures read.
    real(real12) :: energy
    !! The energy of the structure.
    character(50) :: buffer
    !! A buffer for reading strings.
    character(256) :: format
    !! A format string for writing output.
    logical :: success
    !! Boolean for success of file operations.

    integer :: xml_unit, unit
    !! File units.
    integer :: iostat
    !! I/O status.
    integer :: num_files
    !! The number of files in the input directories.
    type(basis_type) :: basis
    !! The basis of the structure.
    character(256), dimension(:), allocatable :: structure_list
    !! The list of structure files.
#ifdef ENABLE_ATHENA
    type(graph_type), dimension(:), allocatable :: graphs
    !! Graph representations of the structures.
    real(real12), dimension(:), allocatable :: &
         labels, labels_train, labels_validate
    !! The labels for the structures.
    real(real12), dimension(:,:), allocatable :: &
         dataset, data_train, data_validate
    !! The dataset for the structures.
#endif


    if(present(gvector_container_template)) then
       gvector_container = gvector_container_template
    else
       gvector_container = gvector_container_type()
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
          write(*,*) "Unknown file format: ", file_format
          stop
       end select
    else
       ifile_format = 0
    end if

    ! inquire(file='gvector_container.dat', exist=success)
    ! if(success) then
    !    call gvector_container%read('gvector_container.dat')
    !    goto 100
    ! end if
    ! ! For each new run of the code, it should populate a new directory and ...
    ! ! ... add to the existing ones.
    ! ! And it should check that output_dir never equals any of the input_dirs (or database_dirs)
    select rank(input_dir)
    rank(0)
       structure_list = [ get_structure_list( input_dir, ifile_format ) ]
    rank(1)
       do i = 1, size(input_dir)
         if(i.eq.1)then
            structure_list = [ get_structure_list( input_dir(i), ifile_format ) ]
         else
            structure_list = [ structure_list, &
                               get_structure_list( input_dir(i), ifile_format )]
         end if
       end do
    end select


    num_structures = 0
#ifdef ENABLE_ATHENA
    allocate(graphs(0))
    allocate(labels(0))
#endif
    do i = 1, size(structure_list)

       write(*,*) "Reading structure: ", trim(adjustl(structure_list(i)))
       select case(ifile_format)
       case(0) ! vasprun.xml
          open(newunit=unit, file=trim(adjustl(structure_list(i)))//"/vasprun.xml")
          write(*,*) "Reading structures from xml"
          basis%energy = get_energy_from_vasprun(unit, success)
          if(.not.success) cycle
          rewind(unit)
          call get_structure_from_vasprun(unit, basis, success)
          if(.not.success) cycle
          close(unit)
       case(1)
          open(newunit=unit, file=trim(adjustl(structure_list(i)))//"/POSCAR")
          write(*,*) "Reading structures from POSCAR"
          call geom_read(unit, basis)
          close(unit)
          open(newunit=unit, file=trim(adjustl(structure_list(i)))//"/OUTCAR")
          call grep(unit, 'free  energy   TOTEN  =', lline=.false., success=success)
          if(.not.success) cycle
          backspace(unit)
          read(unit,*) buffer, buffer, buffer, buffer, energy
          close(unit)
          basis%energy = energy
       case(2)
          open(newunit=unit, file=trim(adjustl(structure_list(i))))
          write(*,*) "Reading structures from xyz"
          do
             read(unit,'(A)',iostat=iostat) buffer
             if(iostat .ne. 0) exit
             if(trim(buffer).eq."") cycle
             backspace(unit)
             call geom_read(unit, basis)
#ifdef ENABLE_ATHENA
             graphs = [ graphs, get_graph_from_basis(basis) ]
             labels = [ labels, basis%energy ]
#endif

             num_structures = num_structures + 1
             write(format,'("(""Found structure: "",I4,"" with energy: "",&
                  &F13.7,"" and formula:"",",I0,"(1X,A,I0))")') basis%nspec
             write(*,format) &
                  num_structures, &
                  basis%energy, &
                  ( trim(basis%spec(j)%name), basis%spec(j)%num, &
                  j=1, basis%nspec )
             call gvector_container%add(basis)
          end do
          cycle
       end select

       write(*,*) &
            "Found structure: ", trim(adjustl(structure_list(i))), &
            " with energy: ", basis%energy
       call gvector_container%add(basis)
       num_structures = num_structures + 1
      
       ! ! STORE THE ENERGY IN AN ARRAY
       ! probably new structure format of crystal
       ! where crystal contains lattice, basis, and energy
       ! ! DO SOMETHING ABOUT NESTED RELAXATIONS
    end do



100 call gvector_container%evolve(deallocate_systems_after_evolve=.false.)


#ifdef ENABLE_ATHENA
    call network_setup(num_inputs = 2, num_outputs = 1)
    call network_train_graph(graphs(:), labels(:), num_epochs = 200)
#endif

    igeom_input = 1

  end function get_evolved_gvectors_from_data
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
    logical :: file_exists, addit_file_exists, structures_found
    !! Booleans for file existence.


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
       write(*,'("No structures found in directory: """,A, &
            &""" with format: ",I0)') &
            trim(input_dir), ifile_format
       stop 0
    end if

  end function get_structure_list
!###############################################################################

end module read_structures