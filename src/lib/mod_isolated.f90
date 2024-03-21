module isolated
  use constants, only: real12
  use rw_geom, only: bas_type, geom_write
  use vasp_file_handler, only: generate_potcar, kpoints_write, Incarwrite
  implicit none

  private

  public :: generate_isolated_calculations

contains

!!!#############################################################################
!!! make isolated atom directories and set up enclosing calculations
!!!#############################################################################
  subroutine generate_isolated_calculations(element_list)
    implicit none
    character(len=3), dimension(:), intent(in) :: element_list

    integer :: unit
    integer :: i, num_species
    type(bas_type) :: basis
    real(real12), dimension(3,3) :: lattice = 0._real12
    logical :: exists
    character(len=1024) :: tmp


    num_species = size(element_list)

    !!--------------------------------------------------------------------------
    !! inquires if the directory 'iso' exists, if not, create it
    !!--------------------------------------------------------------------------
    inquire(file="iso", exist=exists) 
    if(.not.exists)then
       call execute_command_line("mkdir iso") 
    end if

    lattice(1,1) = 20._real12
    lattice(2,2) = 20._real12
    lattice(3,3) = 20._real12
    basis%nspec = 1
    basis%natom = 1
    allocate(basis%spec(1))
    basis%spec(1)%num = 1
    allocate(basis%spec(1)%atom(1,3), source = 0.5_real12)
    !!--------------------------------------------------------------------------
    !!! prepare isolation calculations for each element
    !!--------------------------------------------------------------------------
    do i = 1, num_species 
       basis%sysname = trim(adjustl(element_list(i)))//" isolated"
       basis%spec(1)%name = trim(adjustl(element_list(i)))

       !! write the name of the directory to tmp
       write(tmp,'("iso/POSCAR_",A3)') trim(adjustl(element_list(i)))

       !! check if directory already exists, if not, create it
       !! COMPLAIN IF EXISTS
       inquire(file=tmp, exist=exists)
       if(exists) then
          write(*,'("ERROR: Directory ",A," already exists")') &
               trim(adjustl(tmp))
          stop 1
       else
          !! make directory
          call execute_command_line("mkdir " // trim(adjustl(tmp)))
    
          !! write POSCAR file
          open(newunit=unit, file=trim(adjustl(tmp))//"/POSCAR")
          call geom_write(unit, lattice, basis)
          close(unit)
    
          call generate_potcar(tmp,[element_list(i)]) 
          call kpoints_write(tmp,1,1,1)
          call Incarwrite(tmp,500, 20) !!The 500 here is nstep electronic
       end if
    end do

  end subroutine generate_isolated_calculations
!!!#############################################################################

end module isolated