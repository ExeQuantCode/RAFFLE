module evolve
  use constants, only: real12
  use misc, only: grep
  use rw_geom, only: bas_type, geom_read
  use evolver, only: gvector_container_type, gvector_type
  use file_generator
  use contributions
  implicit none


  private

  public :: bond_evolution


contains

!!!#############################################################################
!!! 
!!!#############################################################################
  subroutine bond_evolution(input_dir)
    implicit none
    character(1024), dimension(..), intent(in) :: input_dir
    character(1024) :: name, read_element_pairing, read_element
    integer :: prev_structures, i, nbin, stat, exitst, exitst2, exitst3
    logical :: dir_e
    real(real12), dimension(:,:), allocatable :: gaussian
    real(real12), dimension(2) :: read_in, norma_vector
    real(real12) :: sigma, bondcut, dist_height, energy
    character(50) :: buffer1, buffer2
    logical :: success

    integer :: previous_structures_unit, xml_unit, unit, ierror
    integer :: num_directories
    type(bas_type) :: basis
    type(gvector_container_type) :: gvector_container
    type(gvector_type) :: gvector
    real(real12), dimension(3,3) :: lattice
    character(100), dimension(:), allocatable :: structure_list


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
       call execute_command_line("ls "//trim(adjustl(input_dir))//">structures.txt",wait=.TRUE.)
       open(newunit=unit, file="structures.txt", status="new")
       do
          read(unit,*,iostat=ierror) name
          structure_list = [ structure_list, input_dir//trim(name) ]
          if(is_iostat_end(ierror)) exit
       end do
       close(unit)
    rank(1)
       num_directories = size(input_dir)
       do i = 1, num_directories
         call execute_command_line("ls "//trim(adjustl(input_dir(i)))//">structures.txt",wait=.TRUE.)
         open(newunit=unit, file="structures.txt", status="new")
         do
            read(unit,*,iostat=ierror) name
            structure_list = [ structure_list, input_dir(i)//trim(name) ]
            if(is_iostat_end(ierror)) exit
         end do
         close(unit)
       end do
    end select


    do i = 1, size(structure_list)
       write(name,'(A,"/vasprun.xml")') trim(adjustl(structure_list(i)))

       open(newunit=unit, file=trim(adjustl(structure_list(i)))//"/POSCAR")
       call geom_read(unit, lattice, basis)
       close(unit)

       open(newunit=unit, file=trim(adjustl(structure_list(i)))//"/OUTCAR")
       call grep(unit, 'free  energy   TOTEN  =', lline=.true., success=success)
       backspace(unit)
       read(unit,*) buffer1, buffer1, buffer1, buffer1, energy
       close(unit)

       call gvector%calculate(lattice, basis)
       gvector%energy = energy

       !!! NO ADD SET UP YET, WORK ON THIS !!!
       !!! Better yet, make it also make calculate as well, just by providing a basis !!!
       !call gvector_container%add(gvector)


       !open(newunit=xml_unit, file=trim(adjustl(structure_list(i)))//"/vasprun.xml")
       !call grep(xml_unit,'   <i name="e_fr_energy">',lline=.true., success=success)
       !if(.not.success) cycle
       !backspace(xml_unit)
       !read(xml_unit,*) buffer1, buffer2, energy
       !close(xml_unit)
      
       !!! STORE THE ENERGY IN AN ARRAY
       ! either energy[], or store it alongside the basis, probably the latter
       ! probably new structure format of crystal
       ! where crystal contains lattice, basis, and energy
       !!! DO SOMETHING ABOUT NESTED RELAXATIONS 
 
    end do


    !!! BEFORE HERE, NEED TO CALCULATE FORMATION ENERGY FOR EACH STRUCTURE !!!
    ! need to read in an isolated energy file (or bulks, it can be user choice)

    call gvector_container%evolve()
    ! once the total has been calculated, deallocate the rest


    !!!! JUST WORK OUT THE GVECTORS HERE !!!!
    !! we have other codes for this
    !! read in the POSCAR (eventually from the xml file)
    !! calculate the neighbour tables


    ! bondcut=5
    ! nbin=1000
    ! sigma=0.1 ! 2-body
    ! sigma=0.05 ! 3-body
    ! sigma=0.05 ! 4-body

  end subroutine bond_evolution
!!!#############################################################################

end module evolve