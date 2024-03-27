module read_structures
  use constants, only: real12
  use misc, only: grep
  use rw_geom, only: bas_type, geom_read
  use evolver, only: gvector_container_type, gvector_type
  implicit none


  private

  public :: bond_evolution


contains

!!!#############################################################################
!!! read in the structures from the input directories and generate the gvectors
!!!#############################################################################
  function bond_evolution(input_dir, element_file, element_list) &
       result(gvector_container)
    implicit none
    character(*), dimension(..), intent(in) :: input_dir
    character(*), intent(in), optional :: element_file  
    type(gvector_container_type) :: gvector_container
    character(len=3), dimension(:), intent(in), optional :: element_list

    character(256) :: name
    integer :: i
    real(real12) :: energy
    character(50) :: buffer
    logical :: success, file_exists

    integer :: xml_unit, unit, ierror
    integer :: num_directories
    type(bas_type) :: basis
    type(gvector_type) :: gvector
    real(real12), dimension(3,3) :: lattice
    character(256), dimension(:), allocatable :: structure_list


    !!! SCRAP ALL OF THIS
    !!! Code should look in a directory called initial/seed/database/data (user defined)
    !!! ... it should scrape that for structures and, for each one it encounters, ...
    !!! ... get the energy and structure, then put it through the ...
    !!! ... gvector_type%calculate procedure.
    !!! Once done, it should then generate new structures in the iteration/ directory ...
    !!! ... (user defined again).
    !!! For each new run of the code, it should populate a new directory and ...
    !!! ... add to the existing ones.
    !!! Would the user give it a list of directories to search for structures?
    !!! And it should check that output_dir never equals any of the input_dirs (or database_dirs)
    select rank(input_dir)
    rank(0)
       num_directories = 1
       call execute_command_line("ls "//trim(adjustl(input_dir))//"/. >structures.txt",wait=.TRUE.)
       open(newunit=unit, file="structures.txt", status="old")
       read(unit,*,iostat=ierror) name
       name = trim(adjustl(input_dir))//trim(adjustl(name))
       structure_list = [ name ]
       do
          read(unit,*,iostat=ierror) name
          name = trim(adjustl(input_dir))//trim(adjustl(name))
          structure_list = [ structure_list, name ]
          if(is_iostat_end(ierror)) exit
       end do
       close(unit)
    rank(1)
       num_directories = size(input_dir)
       do i = 1, num_directories
         call execute_command_line("ls "//trim(adjustl(input_dir(i)))//">structures.txt",wait=.TRUE.)
         open(newunit=unit, file="structures.txt", status="new")
         read(unit,*,iostat=ierror) name
         name = trim(adjustl(input_dir(i)))//trim(adjustl(name))
         structure_list = [ name ]
         do
            read(unit,*,iostat=ierror) name
            name = trim(adjustl(input_dir(i)))//trim(adjustl(name))
            structure_list = [ structure_list, name ]
            if(is_iostat_end(ierror)) exit
         end do
         close(unit)
       end do
    end select


    do i = 1, size(structure_list)
       write(name,'(A,"/vasprun.xml")') trim(adjustl(structure_list(i)))

       inquire(file=trim(adjustl(structure_list(i)))//"/POSCAR", exist=file_exists)
       if(.not.file_exists) cycle
       inquire(file=trim(adjustl(structure_list(i)))//"/OUTCAR", exist=file_exists)
       if(.not.file_exists) cycle

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

       write(*,*) &
            "Found structure: ", trim(adjustl(structure_list(i))), &
            " with energy: ", energy
       call gvector_container%add(basis, lattice)


       !open(newunit=xml_unit, file=trim(adjustl(structure_list(i)))//"/vasprun.xml")
       !call grep(xml_unit,'   <i name="e_fr_energy">',lline=.true., success=success)
       !if(.not.success) cycle
       !backspace(xml_unit)
       !read(xml_unit,*) buffer, buffer, energy
       !close(xml_unit)
      
       !!! STORE THE ENERGY IN AN ARRAY
       ! either energy[], or store it alongside the basis, probably the latter
       ! probably new structure format of crystal
       ! where crystal contains lattice, basis, and energy
       !!! DO SOMETHING ABOUT NESTED RELAXATIONS 
 
    end do

    if(present(element_file).and.present(element_list)) then
       call gvector_container%set_species_list(element_file, element_list)
    end if
    call gvector_container%evolve()


    ! bondcut=5
    ! nbin=1000
    ! sigma=0.1 ! 2-body
    ! sigma=0.05 ! 3-body
    ! sigma=0.05 ! 4-body

  end function bond_evolution
!!!#############################################################################

end module read_structures