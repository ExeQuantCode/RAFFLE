module vasp_file_handler
  use constants, only: real12
  use misc_raffle, only: touch, icount
  implicit none

  private

  public :: Incarwrite, kpoints_write, generate_potcar


contains

!!!#############################################################################
!!! write INCAR file
!!!#############################################################################
  subroutine Incarwrite(filepath,nstep,bandno)
    implicit none
    character(1024) ::  name
    integer :: nstep, bandno
    character(*) :: filepath

    integer :: unit
    
    write(name,'(A,A6)') trim(filepath), "/INCAR"
    write(*,*) "The name of the new file will be", trim(adjustl(name))
    
    !write(*,*) trim(adjustl(name))
    open(newunit=unit,file=trim(adjustl(name)))
    
    write(unit, *)"SYSTEM = RSS"
    write(unit, *)"### sys ###"
    write(unit, *)"GGA = PE "
    write(unit, *)"ISYM = 0"
    write(unit, *)"ENCUT = 500"
    write(unit, *)"ENAUG = 500"
    write(unit, *)"PREC = Accurate"
    write(unit, *)"ISTART = 1"
    write(unit, *)"ICHARG = 2"
    write(unit, *)"LWAVE = .FALSE."
    write(unit, *)"LAECHG =.FALSE."
    write(unit, *)"LMAXMIX = 6"
    write(unit, *)"LASPH = .TRUE."
    write(unit, *)"LMIXTAU = .TRUE."
    write(unit, *)"LREAL = .FALSE."
    write(unit, *)"NWRITE = 3"
    write(unit, *)"ALGO = N"
    write(unit, *)"LHAR = .FALSE."
    write(unit, *)"LVTOT = .FALSE."
    write(unit, *)
    
    write(unit, *)"### vdw ###"
    write(unit, *)"#IVDW = 11"
    write(unit, *)"#VDW_RADIUS = 80"
    write(unit, *)"#VDW_CNRADIUS = 50.0"
    write(unit, *)"#VDW_S8 = 0.722"
    write(unit, *)"#VDW_SR = 1.217"
    write(unit, *) 
    
    write(unit, *)"### elc ###"
    write(unit, *)"#AMIX = 0.6"
    write(name,'(I0.3)') nstep
    write(unit, *)"NELM = ", adjustL(trim(name))
    write(unit, *)"NELMIN = 5"
    write(unit, *)"NELMDL = -5"
    write(unit, '(1X,A,1X,I0.3)') "NBANDS =", bandno
    write(unit, *)"EDIFF = 10d-8"
    write(unit, *)"ISMEAR = 0"
    
    !!! THIS SHOULD BE CHANGED BASED ON INTUITION
    
    write(unit, *)"SIGMA = 0.2"
    write(unit, *) 
    
    write(unit, *)"### Mag ###"
    write(unit, *)"#ISPIN = 2"
    write(unit, *)"#MAGMOM = 1 -1 1 1 -1 -1"
    write(unit, *) 
    
    write(unit, *)"### ncl ###"
    write(unit, *)"#LNONCOLLINEAR = .TRUE."
    write(unit, *)"#LSORBIT = .TRUE."
    write(unit, *)"#GGA_COMPAT = .FALSE."
    write(unit, *) 
    
    write(unit, *)"### MPI ###"
    write(unit, *)"#NCORE = 24"
    write(unit, *)"#KPAR = 2"
    write(unit, *) 
    
    write(unit, *)"### Rlx ###"
    write(unit, *)"#LMAXPAW =-1"
    write(unit, *)"ADDGRID = .TRUE."
    write(unit, *)"#POTIM = 0.1"
    write(unit, *)"#NFREE = 15"
    write(unit, *)"#NSW = 150"
    write(unit, *)"#ISIF = 2"
    write(unit, *)"#IBRION = 1"
    write(unit, *)"#EDIFFG = -0.001"
    write(unit, *) 
    
    write(unit, *)"### dos ###"
    write(unit, *)"#EMIN = -3.0"
    write(unit, *)"#EMAX =  6.0"
    write(unit, *)"#NEDOS = 5000"
    write(unit, *)"#LORBIT = 11"
    
    close(unit)
  end subroutine Incarwrite
!!!#############################################################################


!!!#############################################################################
!!! write POTCAR file
!!!#############################################################################
  subroutine generate_potcar(filepath, element_list) 
    implicit none
    character(20), intent(in) :: filepath
    character(3), dimension(:), intent(in) :: element_list

    integer :: i

    do i=1, size(element_list) 
       if(i.eq.1) then
          call execute_command_line( &
               "cat potcar"//trim(adjustl(element_list(i)))//" > "//&
               &trim(adjustl(filepath))//"/POTCAR")
       else
          call execute_command_line( &
               "cat "//trim(adjustl(filepath))//"/POTCAR"//&
               trim(adjustl(element_list(i)))//" potcar"//&
               trim(adjustl(element_list(i)))//" >> "//&
               trim(adjustl(filepath))//"/POTCAR1")
          call execute_command_line( &
               "mv "//trim(adjustl(filepath))//"/POTCAR1 "//&
               trim(adjustl(filepath))//"/POTCAR")
       end if
    end do

  end subroutine generate_potcar
!!!#############################################################################


!!!#############################################################################
!!! write job and KPOINTS files
!!!#############################################################################
  subroutine kpoints_write(filepath, a,b,c)
    implicit none
    integer :: a,b,c
    character(1024) ::  name
    character(20) :: filepath
    integer :: unit

    name=" "
    !write(*,*) filepath
    write(name,'(A,A8)') trim(filepath), "/KPOINTS"
    open(newunit=unit,file=trim(adjustl(name)), status='new')
    write(unit, '(A7)') "KPOINTS"
    write(unit, '(A1)') "0" 
    write(unit, '(A1)') "G" 
    write(unit, '(I0,X,I0,X,I0)') a, b, c
    write(unit, '(I0,X,I0,X,I0)') 0, 0, 0
    close(unit) 
  end subroutine kpoints_write
!!!#############################################################################

end module vasp_file_handler