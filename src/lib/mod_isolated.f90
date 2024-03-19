module isolated
  use constants, only: real12
  use vasp_file_handler, only: generate_potcar, Jobwrite, Incarwrite
  implicit none

  private

  public :: generate_isolated_calculations

contains

  subroutine generate_isolated_calculations(element_list)
    implicit none
    character(len=3), dimension(:), intent(in) :: element_list

    integer :: unit
    integer :: i, num_species
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


    !!--------------------------------------------------------------------------
    !!! prepare isolation calculations for each element
    !!--------------------------------------------------------------------------
    do i = 1, num_species 

    
       !! write the name of the directory to tmp
       write(tmp,'("iso/POSCAR_",A3)') element_list(i)

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
    
          open(newunit=unit, &
               file="iso/POSCAR_"//trim(adjustl(element_list(i)))//"/POSCAR") 
          write(unit,'(A)') "test"
          write(unit, *) "1.00000000"
          write(unit, '(" 20.0000000000000000  0.0000000000000000  0.0000000000000000")')
          write(unit, '("  0.0000000000000000 20.0000000000000000  0.0000000000000000")')
          write(unit, '("  0.0000000000000000  0.0000000000000000 20.0000000000000000")')
          write(unit,*) element_list(i) 
          write(unit,*) 1
          write(unit,'("Direct")')
          write(unit,'("0.5000000000000000 0.5000000000000000 0.5000000000000000")')
          close(unit)
    
          call generate_potcar(tmp,[element_list(i)]) 
          call Jobwrite(tmp,1,1,1)
          call Incarwrite(tmp,500, 20) !!The 500 here is nstep electronic
          !write(*,*) trim(adjustl(tmp))
          call chdir(trim(adjustl(tmp)))
          call execute_command_line("qsub.sh")
          call execute_command_line("cp ../../job_vasp_isca.in .")
          call chdir("../../")
    
       end if
    end do

  end subroutine generate_isolated_calculations

end module isolated