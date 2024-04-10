module read_structures
  use constants, only: real12
  use misc, only: grep
  use rw_geom, only: bas_type, geom_read, geom_write, igeom_input
  use rw_vasprun, only: get_energy_from_vasprun, get_structure_from_vasprun
  use evolver, only: gvector_container_type, gvector_type
  use machine_learning, only: network_setup, network_train
  implicit none


  private

  public :: get_evolved_gvectors_from_data


contains

!!!#############################################################################
!!! read in the structures from the input directories and generate the gvectors
!!!#############################################################################
  function get_evolved_gvectors_from_data(input_dir, &
       element_file, bond_file, element_list, &
       file_format, gvector_container_template) &
       result(gvector_container)
    implicit none
    character(*), dimension(..), intent(in) :: input_dir
    character(*), intent(in), optional :: element_file, bond_file
    type(gvector_container_type), intent(in), optional :: gvector_container_template
    type(gvector_container_type), allocatable :: gvector_container
    character(len=3), dimension(:), intent(in), optional :: element_list
    character(*), intent(in), optional :: file_format

    character(256) :: name
    integer :: i, ifile_format, num_structures
    real(real12) :: energy
    character(50) :: buffer
    logical :: success, file_exists

    integer :: xml_unit, unit, ierror
    integer :: num_files
    type(bas_type) :: basis
    type(gvector_type) :: gvector
    real(real12), dimension(3,3) :: lattice
    character(256), dimension(:), allocatable :: structure_list

    real(real12), dimension(:,:), allocatable :: dataset, labels


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


    !!! For each new run of the code, it should populate a new directory and ...
    !!! ... add to the existing ones.
    !!! And it should check that output_dir never equals any of the input_dirs (or database_dirs)
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
    do i = 1, size(structure_list)

       select case(ifile_format)
       case(0) ! vasprun.xml
          open(newunit=unit, file=trim(adjustl(structure_list(i)))//"/vasprun.xml")
          basis%energy = get_energy_from_vasprun(unit, success)
          if(.not.success) cycle
          rewind(unit)
          call get_structure_from_vasprun(unit, lattice, basis, success)
          if(.not.success) cycle
          close(unit)
       case(1)
          open(newunit=unit, file=trim(adjustl(structure_list(i)))//"/POSCAR")
          call geom_read(unit, lattice, basis)
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
          do
             read(unit,'(A)',iostat=ierror) buffer
             if(ierror .ne. 0) exit
             if(trim(buffer).eq."") cycle
             backspace(unit)
             call geom_read(unit, lattice, basis)

             write(*,*) &
                  "Found structure: ", trim(adjustl(basis%sysname)), &
                  " in file: ", trim(adjustl(structure_list(i))), &
                  " with energy: ", basis%energy
             call gvector_container%add(basis, lattice)
          end do
          cycle
       end select

       write(*,*) &
            "Found structure: ", trim(adjustl(structure_list(i))), &
            " with energy: ", basis%energy
       call gvector_container%add(basis, lattice)
       num_structures = num_structures + 1
      
       !!! STORE THE ENERGY IN AN ARRAY
       ! probably new structure format of crystal
       ! where crystal contains lattice, basis, and energy
       !!! DO SOMETHING ABOUT NESTED RELAXATIONS
 
    end do

    if(present(element_file).and.present(element_list)) then
       call gvector_container%set_element_info(element_file, element_list)
    end if
    if(present(bond_file)) then
       call gvector_container%set_bond_info(bond_file)
    else
      call gvector_container%set_bond_info()
    end if
    call gvector_container%evolve(deallocate_systems_after_evolve=.false.)

    

    !!! do not deallocate structures
    !!! then load the athena library
    !!! set up the network
    !!! append 2, 3, and 4 body potentials
    !!! HOW DO WE HANDLE SPECIES?
    !!! A network for each species?
    !!! split dataset into train and test sets
    !!! train the network

    write(*,*) "LOOKY", gvector_container%nbins
    allocate(dataset(sum(gvector_container%nbins), num_structures))
    do i = 1, num_structures
       dataset(1:gvector_container%nbins(1),i) = &
            sum(gvector_container%system(i)%df_2body,dim=2)
       dataset(gvector_container%nbins(1)+1:&
            sum(gvector_container%nbins(1:2)),i) = &
            sum(gvector_container%system(i)%df_3body,dim=2)
       dataset(sum(gvector_container%nbins(1:2))+1:&
            sum(gvector_container%nbins(1:3)),i) = &
            sum(gvector_container%system(i)%df_4body,dim=2)
    end do
    allocate(labels(1,num_structures))
    labels(1,:) = [ gvector_container%system(:)%energy ]

    call network_setup(num_inputs = sum(gvector_container%nbins), &
         num_outputs = 1)
    call network_train(dataset, labels, num_epochs = 100)


  end function get_evolved_gvectors_from_data
!!!#############################################################################


!!!#############################################################################
!!! get the list of structures from the input directories
!!!#############################################################################
  function get_structure_list(input_dir, ifile_format) result(structure_list)
    implicit none
    character(*), intent(in) :: input_dir
    integer, intent(in) :: ifile_format
    character(256), dimension(:), allocatable :: structure_list

    integer :: unit, ierror
    character(256) :: name
    logical :: file_exists, addit_file_exists, structures_found


    structures_found = .false.
    call execute_command_line( &
         "ls "//trim(adjustl(input_dir))//"/. >structures.txt", wait = .TRUE. )
    open(newunit=unit, file="structures.txt", status="old")
    do
       read(unit,'(A)',iostat=ierror) name
       if(is_iostat_end(ierror)) exit
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
!!!#############################################################################

end module read_structures