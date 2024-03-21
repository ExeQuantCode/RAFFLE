module vasp_file_handler
  use constants, only: real12
  use misc, only: touch, icount
  use atomtype
  implicit none

  private

  public :: unitcell
  public :: structurecounter
  public :: atomrepeater
  public :: Incarwrite, kpoints_write, generate_potcar, poscar_read
  public :: touchposdir
  public :: get_num_atoms_from_poscar




  type unitcell 
     real(real12), dimension(3,3) :: cell 
  end type unitcell


contains

!!!#############################################################################
!!! repeat the atoms in the unit cell out to a -1,0,1 grid
!!!#############################################################################
  subroutine atomrepeater(structures,position,element,array,unit,atomnumber,length)
    use atomtype
    implicit none
    type (atom), dimension(:,:), allocatable :: array
    real(real12), dimension(3) :: position
    type(unitcell), dimension(:), allocatable :: unit

    integer :: j,atomnumber,length,x,y,z,structures,m
    character(3) :: element
    m=((atomnumber-1)*27)
    do x=-1,1
       do y=-1,1
          do z=-1,1
             m = m + 1
             do j=1, 3
                array(structures,m)%position(j)=position(j) + &
                     &(x*unit(structures)%cell(1,j)) + &
                     &(y*unit(structures)%cell(2,j)) + &
                     &(z*unit(structures)%cell(3,j))
             end do
             array(structures,m)%name=element
          end do
       end do
    end do
  end subroutine atomrepeater
!!!#############################################################################


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
    open(newunit=unit,file=trim(adjustl(name)))!"pos/POSCAR_001")!trim(name))
    
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
    open(newunit=unit,file=trim(adjustl(name)), status='new')!"pos/POSCAR_001")!trim(name))
    write(unit, '(A7)') "KPOINTS"
    write(unit, '(A1)') "0" 
    write(unit, '(A1)') "G" 
    write(unit, '(I0,X,I0,X,I0)') a, b, c
    write(unit, '(I0,X,I0,X,I0)') 0, 0, 0
    close(unit) 
  end subroutine kpoints_write
!!!#############################################################################


!!!#############################################################################
!!! touch POSCAR subdirectories within pos/ directory
!!!#############################################################################
  subroutine touchposdir(structure)
    implicit none
    integer, intent(in) :: structure 
    character(1024) :: buffer

    write(buffer,'("pos/POSCAR_",I0.3)') structure
    call touch(trim(adjustl(buffer)))
  end subroutine touchposdir
!!!#############################################################################


!!!#############################################################################
!!! read POSCAR file
!!!#############################################################################
  subroutine poscar_read(pathway,unit_cell,structure_elements,structure_stochiometry, structure_factor,&
    &coordinate_type,atomic_positions, header)
   use atomtype 
   implicit none 
   character(1024), intent(out) :: coordinate_type, header 
   character(1024) :: read_in, tmp
   character(1024), intent(in) :: pathway
   character(3), dimension(:), allocatable, intent(out) :: structure_elements
   character(3), dimension(:), allocatable :: temp_names
   integer, dimension(:), allocatable, intent(out) :: structure_stochiometry 
   integer :: loop, status, size,eltot, i, j,k,l,m
   integer, parameter :: line_buf_len= 1024*4
   character(LEN=line_buf_len) :: InS
   type(unitcell), dimension(1), intent(out) :: unit_cell 
   real(real12), dimension(:,:), intent(out), allocatable :: atomic_positions
   real(real12), dimension(:,:), allocatable :: temp_positions
   real(real12), intent(out) :: structure_factor 
   logical :: file_e  
   
   eltot=0
   allocate(structure_elements(1))
   inquire(file=trim(adjustl(pathway)), exist=file_e)
   write(*,*) file_e
   if(file_e) then 
    open(102,file=pathway)
   
    read(102,*) header
    read(102,*) structure_factor
    do loop=1, 3
       read(102,*) unit_cell(1)%cell(loop,1), unit_cell(1)%cell(loop,2), unit_cell(1)%cell(loop,3)
    end do
    read(102,'(A)') tmp
    write(*,*) tmp
    write(read_in,'(X,A)') trim(adjustl(tmp))
    write(*,*) read_in
    do i=1,len(read_in)-1
       if(i.eq.1) then
          if((scan(read_in(i:i+1)," ").eq.0).or.&
               &((scan(read_in(i:i+1)," ").eq.2))) then
             eltot=eltot+1
             if(scan(read_in(i:i+1)," ").eq.0) then
                allocate(temp_names(eltot))
   
                do j=1, eltot-1
                   temp_names(j)=structure_elements(j)
                end do
                temp_names(eltot)=read_in(i:i+1)
                deallocate(structure_elements)
                allocate(structure_elements(eltot))
                structure_elements=temp_names
                deallocate(temp_names)
   
             end if
             if((scan(read_in(i:i+1)," ").eq.2).and.(scan(read_in(i-1:i)," ").eq.1)) then
                allocate(temp_names(eltot))
   
                do j=1, eltot-1
                   temp_names(j)=structure_elements(j)
                end do
                temp_names(eltot)=read_in(i:i)
                deallocate(structure_elements)
                allocate(structure_elements(eltot))
                structure_elements=temp_names
                deallocate(temp_names)
             end if
          else
          end if
       else
          if((scan(read_in(i:i+1)," ").eq.0).or.&
               &((scan(read_in(i:i+1)," ").eq.2).and.(scan(read_in(i-1:i)," ")&
               &.eq.1))) then
             eltot=eltot+1
             !write(*,*) tmp(i:i+1)
             if(scan(read_in(i:i+1)," ").eq.0) then
                allocate(temp_names(eltot))
   
                do j=1, eltot-1
                   temp_names(j)=structure_elements(j)
                end do
                temp_names(eltot)=read_in(i:i+1)
                deallocate(structure_elements)
                allocate(structure_elements(eltot))
                structure_elements=temp_names
                deallocate(temp_names)
   
   
             end if
             if((scan(read_in(i:i+1)," ").eq.2).and.(scan(read_in(i-1:i)," ").eq.1)) then
                allocate(temp_names(eltot))
   
                do j=1, eltot-1
                   temp_names(j)=structure_elements(j)
                end do
                temp_names(eltot)=read_in(i:i)
                deallocate(structure_elements)
                allocate(structure_elements(eltot))
                structure_elements=temp_names
                deallocate(temp_names)
   
             end if
          else
          end if
       end if
    end do
    allocate(structure_stochiometry(eltot))
   
    read(102,*) structure_stochiometry 
    read(102,*) coordinate_type
    allocate(temp_positions(sum(structure_stochiometry),3))
    allocate(atomic_positions(sum(structure_stochiometry),3))
    if(coordinate_type.eq."Cartesian") then 
       do i=1, sum(structure_stochiometry)
          read(102,*) atomic_positions(i,1),atomic_positions(i,2),atomic_positions(i,3)
       end do
    else
       write(*,*) coordinate_type
       do i=1, sum(structure_stochiometry)
   
          read(102,*) temp_positions(i,1),temp_positions(i,2),temp_positions(i,3)
   
       end do
       do i=1, sum(structure_stochiometry)
          atomic_positions(i,:)=matmul(temp_positions(i,:),unit_cell(1)%cell) 
          !temp_positions(i,1)*unit_cell(1,1)+temp_positions(i,2)*unit_cell() 
          !write(*,*) atomic_positions(i,:)
          !write(*,*) temp_positions(i,:) 
          !write(*,*) unit_cell(i)%cell 
          !call sleep(5) 
       end do
    end if
   else
    write(*,*) "ERROR: NOT FOUND"
   end if
   
   end subroutine poscar_read
!!!#############################################################################


!!!#############################################################################
!!! count number of structures in a directory
!!!#############################################################################
  function structurecounter(dir_t)  result(n)
    logical :: file_exists
    character(1024) :: name, dir_name, dir_append
    character(3) :: dir_t
    integer :: n, gap_test
    n=1
    select case(dir_t)
    case("don")
       dir_name="/DON_"
       dir_append="/DON"
    case("pos")
       dir_name="/POSCAR_"
       dir_append="/POSCAR"
    case("bon")
       dir_name="/BON_"
       dir_append=""
    case("bad")
       dir_name="/BAD_"
       dir_append=""
    case default
       write(*,*) "ERROR: INVALID DIRECTORY TYPE"
       stop 1
    end select
    mainloop : do while(n.ne.0)
     write(name,'(A3,A,I0.3,A7)') dir_t,trim(adjustl(dir_name)), n,dir_append
     name=trim(adjustl(name))
    
     INQUIRE(FILE=name, EXIST=file_exists)
     if(file_exists) then 
        n=n+1
        cycle 
     else 
        do gap_test=1, 1000
           write(name,'(A3,A,I0.3,A7)') dir_t,trim(adjustl(dir_name)), n+gap_test,dir_append
           INQUIRE(FILE=name, EXIST=file_exists)
           if(file_exists) then
              n=n+gap_test
              cycle mainloop 
           end if
        end do
        n=n-1
    
    
        exit
     end if
    end do mainloop
  end function structurecounter
!!!#############################################################################
   

!!!#############################################################################
!!! get number of atoms from POSCAR file
!!!#############################################################################
  function get_num_atoms_from_poscar(filename_host) result(num_atoms)
    implicit none
    character(1024), intent(in) :: filename_host
    integer :: num_atoms

    integer :: unit, i
    character(1024) :: buffer
    integer, dimension(:), allocatable :: num_species_list


    open(newunit=unit, file=trim(adjustl(filename_host)))
    do i = 1, 6
       read(unit, * )
    end do
    read(unit,'(A)') buffer
    close(unit)
    allocate(num_species_list(icount(buffer)))
    read(buffer,*) num_species_list
    num_atoms = sum(num_species_list)
  end function get_num_atoms_from_poscar
!!!#############################################################################

end module vasp_file_handler