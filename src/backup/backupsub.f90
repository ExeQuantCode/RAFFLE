module help
use atomtype
implicit none


TYPE densitymatrix 
   double precision, dimension(3) :: position 
   double precision :: density 
   integer :: checked
end type densitymatrix

TYPE unitcell 
   double precision, dimension(3,3) :: cell 
end type unitcell

contains
subroutine touchpos() 
  
  character(1024) :: file, command 
  logical :: dir_e

  inquire(file="pos", exist=dir_e)
  if(dir_e) then
  else
     write(command,*) "mkdir pos"
     CALL execute_command_line(command)
  end if
end subroutine touchpos

subroutine touchposdir(structures,prev_structures)
  character(1024) :: name, tmp, command  
  logical :: dir_e
  integer :: structures, prev_structures 
  write(name,'(A11,I0.3,A7)')"pos/POSCAR_",structures+prev_structures,"/POSCAR"
   !!Calculates the new structure number, and writes it to tmp
  write(tmp,'(A11,I0.3)')"pos/POSCAR_",structures+prev_structures
  !!Checks if a directory to contain that file exsts already (It should never exist, could add warning)
  inquire(file=tmp, exist=dir_e)
  if(dir_e) then
  else
     !!Writes a command to create said directory
     write(command,'(A17,I0.3)')"mkdir pos/POSCAR_",structures+prev_structures
     CALL execute_command_line(trim(adjustl(command)))
  end if
end subroutine touchposdir

       

  
  recursive subroutine invar(a,b,c)   
    implicit none 
    integer :: a,i, tmpvar
    character(1024) :: buffer, command
    integer, dimension(:), allocatable :: b
    character(3), dimension(:), allocatable, intent(out) :: c 
    
    
    open(71, file="Infile.txt") 
    rewind(71)
    do i=1, a 
       read(71,'(A,X,A,X,A,X,A,X,A,X,A,X,A,X,A)') buffer
    end do
    close(71) 
   ! print*, buffer
    write(command,'(4A,X,A,X,A,X,A)') "(echo",' "',trim(adjustl(buffer)),'" | ','sed -e "s/^.*=//g")>tmp.txt'
    CALL execute_command_line(command)
    !print*, command
    open(72,file="tmp.txt")
    
    if(a.le.3) then
       allocate(b(1))
       read(72,*) b
       close(72)
    else if(a.eq.4) then 
        close(72)
       CALL invar(2,b,c)

       open(71, file="Infile.txt")
       rewind(71)
       do i=1, a
          read(71,'(A,X,A,X,A,X,A,X,A,X,A,X,A,X,A)') buffer
       end do
       close(71)
      ! print*, buffer
       write(command,'(4A,X,A,X,A,X,A)') "(echo",' "',trim(adjustl(buffer)),'" | ','sed -e "s/^.*=//g")>tmp.txt'
       CALL execute_command_line(command)
      ! print*, command
       open(72,file="tmp.txt")

       tmpvar=b(1)
       deallocate(b)
       allocate(c(tmpvar))
       read(72,*) c

       close(72)
   
    else if(a.eq.5) then
       close(72)
       CALL invar(2,b,c)
       
       open(71, file="Infile.txt")
       rewind(71)
       do i=1, a
          read(71,'(A,X,A,X,A,X,A,X,A,X,A,X,A,X,A)') buffer
       end do
       close(71)
      ! print*, buffer
       write(command,'(4A,X,A,X,A,X,A)') "(echo",' "',trim(adjustl(buffer)),'" | ','sed -e "s/^.*=//g")>tmp.txt'
       CALL execute_command_line(command)
      ! print*, command
       open(72,file="tmp.txt")
       
       tmpvar=b(1)
       deallocate(b)
       allocate(b(tmpvar))
       read(72,*) b
       print*, b
       close(72)
    else if(a.le.10) then
       allocate(b(1))
       read(72,*) b
       close(72)
    end if
    
    
    write(command,*) "rm tmp.txt"
    CALL execute_command_line(command)
  end subroutine invar
!!!-------------------------------------------------------------!!!
!!!This function provides the cross product of two vectors      !!!
!!!-------------------------------------------------------------!!!

function cross(a,b) result(axb)
implicit none
integer,parameter :: wp=selected_real_kind(15, 307) !double precision
double precision,dimension(3) :: axb
double precision,dimension(3),intent(in) :: a
double precision,dimension(3),intent(in) :: b 
axb(1) = a(2)*b(3) - a(3)*b(2)
axb(2) = a(3)*b(1) - a(1)*b(3)
axb(3) = a(1)*b(2) - a(2)*b(1)
end function cross

!!!---------------------------------------------------------------------------!!!
!!!This function writes a list of atomic positions to the POSCAR in cartesian !!!
!!!---------------------------------------------------------------------------!!!

subroutine  poswrite(box,atomlist,len, structures, structno, prev_structures)

integer :: i,j, len, reason, structures, structno, prev_structures
type(atom), dimension(:,:) :: atomlist
double precision, dimension(3,3) :: box

print*, box
do i=1, len 
   if(i.eq.len) then 
      write(structures+10000,'(5X,A2)', IOStat=reason, advance="no")&
           &atomlist(structures,i)%name
    else 
       !if(i.eq.(len-1)) cycle
      if (atomlist(structures,i)%name.eq.atomlist(structures,i+1)%name) cycle
      write(structures+10000,'(5X,A2)', IOStat=reason, advance="no") &
           &atomlist(structures,i)%name      
   end if
end do

write(structures+10000, '(10X)')
!write(structures+10000,*) "HERE"
j=1
do i=1, len 
   if(i.eq.len) then 
      write(structures+10000,'(10X,I3)', IOStat=reason, advance="no") j
   else
      if (atomlist(structures,i)%name.eq.atomlist(structures,i+1)%name) then 
         j=j+1
         cycle
      else 
         write(structures+10000,'(9X,I3)', IOStat=reason, advance="no") j
         j=1
  
      end if
   end if
end do


write(structures+10000,'(/,A)') "Cartesian"

do i=1, len 
   write(structures+10000,'(6X,F0.16,6X,F0.16,6X,F0.16)') atomlist(structures,i)%position(:)
end do

end subroutine poswrite

!!!-------------------------------------------------!!!
!!!This function returns normvol as the cell volume !!!
!!!-------------------------------------------------!!!

function sphereoverlap(r,rp,b, pi) result(volume) 
  double precision :: r,rp,volume, b, step1, step2, pi, step3
  step1=pi*(r+rp-b)**2
  step2=(b**2)+(2*b*rp)-(3*(rp**2))+(2*b*r)+(6*rp*r)-3*(r**2)
  step3=1.0/(12.0*b)
!  print*, "Sub.f90 ->" ,step1*step2*step3
  volume=step1*step2*step3

end function sphereoverlap
  
function  cellvol(x) result(normvol)
  double precision, dimension(3,3) :: x
  double precision, dimension(3) :: a,b,c, tmp
  integer :: i
  double precision :: normvol
  a=x(:,1)
  b=x(:,2)
  c=x(:,3)
  tmp=cross(b,c)
  normvol=0
  do i=1,3
     normvol=normvol+a(i)*tmp(i)
  end do
  !print*, "Cell volume is", abs(normvol)
end function cellvol


!!!Bondlength!!!
function bondlength(A,B) result(C)
  double precision, dimension(3), intent(in) :: A,B
  double precision :: C
  C= (((A(1)-B(1))**2)+((A(2)-B(2))**2)+((A(3)-B(3))**2))**0.5
end function bondlength

function structurecounter(dir_t)  result(n)
  logical :: file_exists
  character(1024) :: name, dir_name, dir_append
  character(3) :: dir_t
  integer :: n
  n=1
  do while(n.ne.0)
     if(dir_t.eq."don") dir_name="/DON_"     
     if(dir_t.eq."pos") dir_name="/POSCAR_"
     if(dir_t.eq."don") dir_append="/DON"
     if(dir_t.eq."pos") dir_append="/POSCAR"
     if(dir_t.eq."bon") dir_name="/BON_"
     if(dir_t.eq."bon") dir_append="/BON"
     write(name,'(A3,A,I0.3,A7)') dir_t,trim(adjustl(dir_name)), n,dir_append
     name=trim(adjustl(name))
     INQUIRE(FILE=name, EXIST=file_exists)
     if(file_exists) then 
        n=n+1
        cycle 
     else 
        n=n-1

        exit
     end if
  end do
end function structurecounter

function bondangle(A,B,C) result(theta)
double precision :: theta, x
double precision, dimension(3) :: A,B,C, bond1, bond2
integer :: i
bond1(:)=A(:)-B(:)
bond2(:)=C(:)-B(:)
x=0
do i=1,3
x=x+bond1(i)*bond2(i)
end do
x=x/(bondlength(A,B)*bondlength(B,C))
theta=(180.0/3.141592654)*acos(x)


!print*, A,B,C,bond1, bond2,"theta is", theta*180.0/3.141592654
end function bondangle

subroutine Incarwrite(filepath,nstep,bandno) 
character(1024) ::  name
integer :: nstep, bandno
character(1024) :: filepath


write(name,'(A,A6)') trim(filepath), "/INCAR"
!print*, trim(adjustl(name))
open(unit=11,file=trim(adjustl(name)), status='new')!"pos/POSCAR_001")!trim(name))

write(11, *)"SYSTEM = RSS"
write(11, *)"### sys ###"
write(11, *)"GGA = PE "
write(11, *)"ISYM = 0"
write(11, *)"ENCUT = 500"
write(11, *)"ENAUG = 500"
write(11, *)"PREC = Accurate"
write(11, *)"ISTART = 1"
write(11, *)"ICHARG = 2"
write(11, *)"LWAVE = .FALSE."
write(11, *)"LAECHG =.FALSE."
write(11, *)"LMAXMIX = 6"
write(11, *)"LASPH = .TRUE."
write(11, *)"LMIXTAU = .TRUE."
write(11, *)"LREAL = .FALSE."
write(11, *)"NWRITE = 3"
write(11, *)"ALGO = N"
write(11, *)"LHAR = .FALSE."
write(11, *)"LVTOT = .FALSE."
write(11, *)

write(11, *)"### vdw ###"
write(11, *)"#IVDW = 11"
write(11, *)"#VDW_RADIUS = 80"
write(11, *)"#VDW_CNRADIUS = 50.0"
write(11, *)"#VDW_S8 = 0.722"
write(11, *)"#VDW_SR = 1.217"
write(11, *) 

write(11, *)"### elc ###"
write(11, *)"#AMIX = 0.6"
write(name,'(I0.3)') nstep
write(11, *)"NELM = ", adjustL(trim(name))
write(11, *)"NELMIN = 5"
write(11, *)"NELMDL = -5"
write(11, '(1X,A,1X,I0.3)') "NBANDS =", bandno
write(11, *)"EDIFF = 10d-8"
write(11, *)"ISMEAR = 0"

!!! THIS SHOULD BE CHANGED BASED ON INTUITION

write(11, *)"SIGMA = 0.2"
write(11, *) 

write(11, *)"### Mag ###"
write(11, *)"ISPIN = 2"
write(11, *)"#MAGMOM = 1 -1 1 1 -1 -1"
write(11, *) 

write(11, *)"### ncl ###"
write(11, *)"#LNONCOLLINEAR = .TRUE."
write(11, *)"#LSORBIT = .TRUE."
write(11, *)"#GGA_COMPAT = .FALSE."
write(11, *) 

write(11, *)"### MPI ###"
write(11, *)"#NCORE = 24"
write(11, *)"#KPAR = 2"
write(11, *) 

write(11, *)"### Rlx ###"
write(11, *)"#LMAXPAW =-1"
write(11, *)"ADDGRID = .TRUE."
write(11, *)"#POTIM = 0.1"
write(11, *)"#NFREE = 15"
write(11, *)"#NSW = 150"
write(11, *)"#ISIF = 2"
write(11, *)"#IBRION = 1"
write(11, *)"#EDIFFG = -0.001"
write(11, *) 

write(11, *)"### dos ###"
write(11, *)"#EMIN = -3.0"
write(11, *)"#EMAX =  6.0"
write(11, *)"#NEDOS = 5000"
write(11, *)"#LORBIT = 11"

close(11)
end subroutine Incarwrite

subroutine Jobwrite(filepath, a,b,c) 

 integer :: a,b,c
character(1024) ::  name
character(20) :: filepath
name=" "
!print*, filepath
write(name,'(A,A8)') trim(filepath), "/KPOINTS"
open(unit=11,file=trim(adjustl(name)), status='new')!"pos/POSCAR_001")!trim(name))
write(11, '(A7)') "KPOINTS"
write(11, '(A1)') "0" 
write(11, '(A1)') "G" 
write(11, '(I1,1X,I1,1X,I1)') a, b, c
write(11, '(I1,1X,I1,1X,I1)') 0,0,0
close(11) 
end subroutine Jobwrite





subroutine potwrite(filepath, elnames, eltot) 

  character(3), dimension(:), allocatable :: elnames
  character(1024) ::  name
  character(20) :: filepath
  integer :: eltot, i
  name=" "
  !print*, filepath
  !print*, eltot
  do i=1, eltot 
     
     !print*, i
     if(i.eq.1) then
        !print*, "I reached here"
        write(name, '(A3,1X,A6,A2,1X,A1,A,A)') "cat", "potcar", elnames(i),">",trim(adjustl(filepath)),"/POTCAR"
        !print*, name
        call execute_command_line(name)
     else       
        write(name, '(A,1X,A,A,1X,A,A,1X,A,1X,A,A)') "cat ",trim(adjustl(filepath)),"/POTCAR", "potcar", elnames(i)," >>",&
             & trim(adjustl(filepath)),"/POTCAR1"   
        
        call execute_command_line(name)
        write(name,'(A,1X,A,A8,1X,A,A)') "mv ",trim(adjustl(filepath)),"/POTCAR1", trim(adjustl(filepath)),"/POTCAR"
        call execute_command_line(name)
     end if

  end do
end subroutine potwrite
end module help
