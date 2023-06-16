module gen
use help
use atomtype
implicit none 
contains 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine generation(leng, atomlist, alistrep, spacelist, formula, structno, options, eltot, elnames, stochio, elrad) 
  type(unitcell), dimension(:), allocatable :: formula
  double precision :: posneg, r, pi, meanvol, q, normvol, cellmultiplier, calc,sigma1
  integer :: l,b,leng, i, j, k, x, y, z, m, structures, structno, prev_structures, modeselect
  integer :: errorcounter, ecount, eltot, options, loopcounter

  integer :: eltype
  !! box = initial untransformed cubic unit cell 
  !! a 0 0
  !! 0 b 0
  !! 0 0 c
  double precision, dimension(3,3) :: box
  double precision, dimension(3) :: angle, spacelist, tmpvector
  double precision, dimension(:,:,:,:), allocatable :: bondlist
  double precision ::  connectivity, volmin, volmax, bondpro1, bondpro2, distribution, tmpvalue
  type (atom), dimension(:,:), allocatable :: atomlist, alistrep, alistrepp
  character(1024) :: name, tmp, command
  character(3), dimension(:), allocatable :: elnames, sing_el
  integer, dimension(:), allocatable :: elno, stochio 
  logical :: dir_e
  double precision, dimension(:,:,:), allocatable :: elrad
  double precision, dimension(:,:), allocatable :: bondavg, bondminimum

  character(3), dimension(:), allocatable :: tmpels
  integer, dimension(:), allocatable :: tmpdig
  double precision, dimension(:,:), allocatable :: tempmatrix
 

  !! The info file doesn't contain much of use yet. Could build in if relevant 
  open(11, file="Info")
  !! Atomlist contains the information about the positions of all atoms in ALL structures. This may cause issues 
  !! with memory when large numbers of structures are used, may consider breaking down into seperate iterations 
  !! (e.g paralyse) 
  allocate(atomlist(structno,leng))
  !! alistrep contains positions of all atoms repeated in adjacent unit cells. Same point as above. 
  allocate(alistrep(structno,leng*27))
  !! Calls the function structurecounter, which provides information about the number of currently existing 
  !! structures in the directory
  prev_structures=structurecounter("pos")
  !! assigns the length of elno to eltot. NOT SURE WHY, SHOULD BE LENG?. UNLESS ELNO CONTAINS ALL MATERIAL SPECS
  allocate(elno(eltot))



  !! bondlist is a list of ALL the bonds between all the atoms and each of it's neighbours in the first tier of recursive repeated unit cells
  allocate(bondlist(leng,leng*27,eltot,eltot))
  allocate(bondavg(eltot,eltot))
  !! modeselect=1 is a special option allowing a new poscar to be added in at user specification
  modeselect=options
  
  !! Could implement structno>1 in the future for large imports
  if(modeselect.eq.1) structno=1
  structures=1
  loopcounter=0

   !! elnames is a 1D array containing the symbol for each of the atoms (length=eltot). [generator;eltot~>main;elno]----------------------------------------------------------!
!Sets the value of pi. Will eventually come from a dedicated constants and parameters page! 
  pi=3.14159265358979323846                                                              !
!-----------------------------------------------------------------------------------------!

inquire(file="iso", exist=dir_e) 
if(dir_e) then 
else 
write(command,*) "mkdir iso" 
CALL execute_command_line(command) 
end if

!!! This section prepares isolation calculations
do i=1, eltot 

  write(name,'(A11,A,A7)')"iso/POSCAR_",trim(adjustl(elnames(i))),"/POSCAR"
   !!Calculates the new structure number, and writes it to tmp
   write(tmp,'(A11,A3)')"iso/POSCAR_",elnames(i)
   !!Checks if a directory to contain that file exsts already (It should never exist, coul add warning) 
   inquire(file=tmp, exist=dir_e)
   if(dir_e) then 
   else
      !!Writes a command to create said directory
      write(command,'(A17,A3)')"mkdir iso/POSCAR_",elnames(i)
      CALL execute_command_line(command)
      open(10+i, file=name) 
      write(10+i,'(A)') "test"
      write(10+i, *) "1.00000000"
      write(10+i, *) "   ","20.0000000000000000","        ","0.0000000000000000","        ","0.0000000000000000"
      write(10+i, *) "   ","0.0000000000000000","        ","20.0000000000000000","        ","0.0000000000000000"
      write(10+i, *) "   ","0.0000000000000000","        ","0.0000000000000000","        ","20.0000000000000000"
      write(10+i,*) "     ", elnames(i) 
      write(10+i,*) "          ", "1"
      write(10+i,*) "Direct"
      write(10+i,*) "        ","0.5000000000000000","        ","0.5000000000000000","        ","0.5000000000000000"
      close(10+i)

      allocate(sing_el(1)) 
      sing_el(1)=elnames(i) 
      call potwrite(tmp,sing_el,1) 
      deallocate(sing_el)
      call Jobwrite(tmp,1,1,1)
      call Incarwrite(tmp,500, 20) !!The 500 here is nstep electronic
      !print*, trim(adjustl(tmp))
      CALL chdir(trim(adjustl(tmp)))
      CALL execute_command_line("qsub.sh")
      CALL chdir("../../")
   end if
end do












!! Builds pos and don subfolders if they do not already exist 

call touchpos()

inquire(file="don", exist=dir_e) 
if(dir_e) then 
else 
   write(command,*) "mkdir don"
   CALL execute_command_line(command) 
end if






!!! Create the directories for all of the POSCARS to be placed into !!!
     !! BIGLOOP is the parent loop for all procesess, generating one structure for each full completed iteration
b=0
     BIGLOOP: do while(structures.le.structno)
        
        b=structures
        do while (b.gt.7) 
           b=b-7
        end do
        print*, b
        call touchposdir(structures,prev_structures)
!------------------------------------------------------------------------------------------------!   
!    !!This section will count the loop number efficiently if no other prints are used in loop   ! 
!     loopcounter=loopcounter+1                                                                  !
!     write(6, '(A1, 40X, A1)', advance='no') achar(13), achar(13)                               !
!     write(6,'(I0.4)', advance='no') loopcounter                                                !
!------------------------------------------------------------------------------------------------!


!!!-------------------------------------------------------------------------------------!!!
!!!Decides if pseudorandom volume will be larger or smaller than atomic volume summation!!!
!!!-------------------------------------------------------------------------------------!!!
     
     posneg=1
     call random_number(r)
     if(r.lt.0.5) then 
        posneg=-posneg
     else 
     end if

!!!-------------------------------------------------------------------------------------!!!
!!!Decides the total desired pseudorandom cell volume                                   !!! 
!!!-------------------------------------------------------------------------------------!!!


!! Meanvol takes the atomic radius and calculates a guestimate for the total rough cell volume. NEED A BETTER METHOD
!! FOR ACCOMPLISHING THIS
!meanvol=4/3*2.2**3*pi*leng
meanvol=0
volmin=0
volmax=0 
k=0
do i=1, eltot
   do j=1, eltot
!      if(elrad(3,i,i).gt.elrad(3,i,j)) print*, "This element is bonded to more of it's partners than they are to it"
!      if(elrad(3,i,i).lt.elrad(3,i,j)) print*, "This element is bonded to less of it's partners than they are to it"
      
   end do
end do

k=0
do i=1, eltot
   volmin=volmin+stochio(i)*(4.0/3.0)*pi*(elrad(2,i,i)**3)
   do j=1, eltot
      if(j.lt.i) cycle
      if(i.eq.j) then; 
         connectivity=0.5
      else 
         CALL invar(6,tmpdig,tmpels)
         connectivity=dble(tmpdig(1)/100.0)
         !print*, connectivity 
         deallocate(tmpdig)
      end if
 
      if(i.eq.j) then 
      volmin=volmin-connectivity*stochio(i)*elrad(3,i,j)*((dble(stochio(j))/leng)*&
           &sphereoverlap(elrad(2,i,i),elrad(2,j,j),elrad(1,i,j),pi))
      else 
         volmin=volmin-connectivity*(stochio(i)*elrad(3,i,j)+stochio(j)*elrad(4,i,j))*0.5*((dble(stochio(j))/leng)*&
           &sphereoverlap(elrad(2,i,i),elrad(2,j,j),elrad(1,i,j),pi))

      end if
 !"     print*, volmin, connectivity*stochio(i)*elrad(3,i,j)*((dble(stochio(j))/leng)*&
!           &sphereoverlap(elrad(2,i,i),elrad(2,j,j),elrad(1,i,j),pi))
      
         !print*,"bonding between elements",(min((elrad(3,i,j)*stochio(i)),(elrad(4,i,j)*stochio(j))))
         
            !volmax=volmax+stochio(i)*(4.0/3.0)*pi*(elrad(2,i,i)**3)-connectivity*((dble(stochio(j))/leng)*&
      !     &sphereoverlap(elrad(2,i,i),elrad(2,j,j),elrad(1,i,j),pi))
      !print*, volmax
      
      
      
      meanvol=volmin
   end do
   !  meanvol=meanvol+stochio(i)*((connectivity*elrad(1,i,i))+((1.0-connectivity)*elrad(2,i,i)))**3*(4/3)*pi 
end do 

!!Adds or subtracts a small quantity from the calculated volume 
call random_number(r) 
call invar(7,tmpdig,tmpels)
meanvol=meanvol+(dble(tmpdig(1)/100.0)*r*posneg*meanvol)
deallocate(tmpdig)

!!!-------------------------------------------------------------------------------------!!!
!!Define the random unit cell lengths                                                  !!!
!!!-------------------------------------------------------------------------------------!!!

!! Initialises box, which is a cubic unit serving as the basis for the random unit cel
box=0
!! q keeps a running total of the "volume" in the loosest sense of the word  
q=1

call random_number(r)
box(1,1)=0.75+r*2.25
q=q*r
call random_number(r)
box(2,2)=0.75+r*2.25
q=q*r
call random_number(r)
box(3,3)=0.75+r*2.25
q=q*r*1000

!!!--------------------------------------------------------------!!!
!!!Sets the random angles between the unit vectors between 60-120!!!
!!!--------------------------------------------------------------!!!
 angle(:)=0
!!! BRAVAIS LATTICES 
 if (b.eq.1) then !!! Triclinic 
    do i=1, 3
       call random_number(r)
       r=r*60.0+60.0                                                            !   
       r=(r*pi)/(180.0)                                                         !     
       angle(i)=r 
    end do
 else if (b.eq.2) then !!! Cubic
    q=1
    call random_number(r) 
    box(1,1)=0.75+r*2.25
    q=q*(r**3)*1000.0
    box(2,2)=box(1,1) 
    box(3,3)=box(1,1)
    angle(:)=pi/2.0
 else if (b.eq.3) then !!! Monoclinic 
    angle(3)=pi/2
    angle(1)=pi/2
    call random_number(r) 
    r=r*60.0+60.0
    r=(r*pi)/(180.0)
    angle(2)=r
 else if (b.eq.4) then !!! Orthorhombic 
    angle(:)=pi/2.0
 else if (b.eq.5) then !!! Tetragonal B
    q=1 
    call random_number(r) 
    box(1,1)=0.75+r*2.25
    box(2,2)=box(1,1) 
    q=q*(r**2)
    call random_number(r) 
    box(3,3)=0.75+r*2.25 
    q=q*r*1000
    angle(:)=pi/2.0

 else if (b.eq.6) then !!! Rhombohedral very broken/ Trigonal :-(
    q=1 
    call random_number(r) 

    box(1,1)=0.75+r*2.25
    ! box(1,1)=5.0
    box(2,2)=box(1,1) 
    box(3,3)=box(1,1) 
    q=q*(r**3)*1000 
    !print*, box(1,1)
    call random_number(r)
    r=(r*60.0)+60.0 
    r=(r*pi)/(180.0)
    angle(:)=r
    ! angle(:)= 1.75*pi/3.0
    !print*, angle(:)
 else if (b.eq.7) then !!! hexagonal 
    angle(1)=pi/2.0
    angle(2)=pi/2.0
    angle(3)=2.0*pi/3.0
    call random_number(r) 
    box(1,1)=0.75+r*2.25
    box(2,2)=box(1,1) 
    q=r**2
    call random_number(r) 
    box(3,3)=0.75+r*2.25
    q=q*r*1000
 end if
 

!!!---------------------------------!!!
!!Sets the new unit vectors         !!!
!!!---------------------------------!!!

!! Creates a TYPE called formula that contains all the info for a unit cell that is random

calc=sin(angle(1))**2-cos(angle(2))**2-cos(angle(3))**2
calc=calc+cos(angle(1))*cos(angle(2))*cos(angle(3))*2
calc=sqrt(calc)
calc=calc/(sin(angle(1)))

formula(structures)%cell=0
formula(structures)%cell(1,1)=box(1,1)*calc
formula(structures)%cell(2,1)=box(1,1)*(cos(angle(3))-cos(angle(1))*cos(angle(2)))/(sin(angle(1)))
formula(structures)%cell(3,1)=box(1,1)*cos(angle(2)) 
formula(structures)%cell(2,2)=box(2,2)*sin(angle(1)) 
formula(structures)%cell(3,2)=box(2,2)*cos(angle(1))      
formula(structures)%cell(3,3)=box(3,3)


write(name,'(A11,I0.3,A7)')"pos/POSCAR_",structures+prev_structures,"/POSCAR"
open(structures+10000, file=name)
write(structures+10000,*) "Test"
write(structures+10000,*) 1.0


!! Adjusts the the lengths of the unit cell vectors by the 
normvol=cellvol(formula(structures)%cell)
normvol=abs(normvol)/meanvol

normvol=normvol**(1.0/3.0)

do j=1, 3
   do i=1, 3
      formula%cell(i,j)=formula(structures)%cell(i,j)/normvol      
  
   end do
   write(structures+10000,*) formula(structures)%cell(:,j)
end do

!!!----------------------------------------------!!!
!!!Set everything to what it should be and places the atoms
!!!----------------------------------------------!!!

!!! Assigns all the atoms to the correct species label. 
k=0
z=0
m=0
l=0
do j=1, eltot

k=k+stochio(j)
m=m+27*stochio(j)
   do i=1, leng
      if(modeselect.eq.1) exit         
      if((i.le.k).and.(i.gt.z)) then  
         atomlist(structures,i)%name=elnames(j)          
      end if
   end do
    do i=1, leng*27
      if(modeselect.eq.1) exit
      if((i.le.m).and.(i.gt.l)) then
         alistrep(structures,i)%name=elnames(j)
      end if
   end do

z=k
l=m
end do


!!!!!!!!!!!!! This section places the atoms via distribution !!!!!!!!!!!!!!!!!!!!!!!!!!!!
do j=1, 3
   call random_number(r) 
   atomlist(structures,1)%position(j)=r
   alistrep(structures,1)%position(j)=r
end do
errorcounter=0
atomlist(structures,1)%position(:)=matmul(formula(structures)%cell,atomlist(structures,1)%position(:))
call atomrepeater(structures,atomlist(structures,1)%position,alistrep,formula,1,leng)
call invar(8,tmpdig,tmpels) 
b=tmpdig(1) 
deallocate(tmpdig) 
print*, dble(b/100.0)
open(99,file="errorfile") 
sigma1=0.1
i=1
aloop: do while (i.le.leng-1)
   if(errorcounter.gt.1000) then
      !sigma1=sigma1*1.01
      !print*, "Sigma being increased", sigma1

      errorcounter=0
   end if
   tmpvector=0
   do j=1, 3
      call random_number(r)
      tmpvector(j)=r
   end do
   distribution=0
   tmpvector=matmul(formula(structures)%cell,tmpvector)
   
   do y=1, i*27
      if(y.eq.i+1) cycle
      posneg=1
      do j=1, eltot
         do x=1, eltot
            if(atomlist(structures,i+1)%name.ne.elnames(j)) cycle 
            if(alistrep(structures,y)%name.ne.elnames(x)) cycle 
            if(bondlength(tmpvector,alistrep(structures,y)%position).lt.&
                 &elrad(1,j,x)*dble(b/100.0)) then 
               write(99,*) "Atoms too close together"
               errorcounter=errorcounter+1
               cycle aloop
            end if
            distribution=distribution*(y-1)
            tmpvalue=bondlength(tmpvector,alistrep(structures,y)%position)
            distribution=distribution+&
            &(1.0/4.0)*abs(erf((tmpvalue&
                 &-elrad(1,j,x)+0.5*sigma1)/(sigma1*sqrt(2.0))))-&
                 &(1.0/4.0)*abs(erf((tmpvalue&
                 &-elrad(1,j,x)-0.5*sigma1)/(sigma1*sqrt(2.0))))
            if(tmpvalue.lt.10) then 
                 distribution=distribution+0.5*sigma1*0.1
                  
               end if
            distribution=dble(distribution/(1.0*y))
         end do
      end do
   end do
   call random_number(r)
   if(r.gt.distribution) then 
      errorcounter=errorcounter+1
      cycle aloop
   end if
   print*,r, distribution, i
   atomlist(structures,i+1)%position=tmpvector
   call atomrepeater(structures,atomlist(structures,i+1)%position,alistrep,formula,i+1,leng)
   if(i.eq.leng-1) exit  
   i=i+1
   errorcounter=0
end do aloop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Build in some failsafes - changing sigma iteratively or change the probability width
!can also change the minimum allowed bondlength now!! Wouldn't use this too often

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Generate the atomic bonding information files used in upcoming learning algorithm

do i=1, leng
call generatebondfiles(structures,atomlist,alistrep,eltot,stochio,i)
end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
j=0
k=0


CALL invar(8,tmpdig,tmpels)
do i=1, eltot 
   do k=1, eltot 
      do x=1, leng*27
         if(modeselect.eq.1) exit
         do y=1, leng*27        
            if(x.ge.y) cycle
            if(alistrep(structures,x)%name.ne.elnames(k)) cycle 
            if(alistrep(structures,y)%name.ne.elnames(i)) cycle
            if(bondlength(alistrep(structures,x)%position,alistrep(structures,y)%position)&
                 &.lt.(elrad(1,k,i)*dble(tmpdig(1)/100.0))) then
               j=j+1        
               close(structures+10000)
               print*, "Terminating. Bond lengths too short" 
               deallocate(tmpdig)
               cycle BIGLOOP
            end if
         end do
      end do
   end do
end do
deallocate(tmpdig)
m=j
!! The following determines the total number of bonds that each atom has made
!! This is defined by the MAXBOND input parameter
CALL invar(9,tmpdig,tmpels)
do i=1, eltot 
   do k=1, eltot
      if(i.gt.k) cycle
      do x=1, leng
           
         if(modeselect.eq.1) exit
         do y=1, leng*27 
            if(atomlist(structures,x)%name.ne.elnames(k)) cycle
            if(alistrep(structures,y)%name.ne.elnames(i)) cycle
            bondlist(x,y,i,k)=bondlength(atomlist(structures,x)%position,alistrep(structures,y)%position)
            bondlist(y,x,k,i)=bondlist(x,y,i,k)
            if(bondlist(x,y,i,k).gt.(dble(tmpdig(1)/100.0)*elrad(1,i,k))) then;
               bondlist(x,y,i,k)=0
               bondlist(y,x,k,i)=0
               cycle
            end if
         end do
      end do
   end do
end do
deallocate(tmpdig)

!! The following, for each species pairing, generates an average bondlength. At the end of
!! the first eltot loop, the value m can be read out to be the total bonding for that pairing 
bondavg=0
inquire(file="bonddata.txt", exist=dir_e)
     if(dir_e) then
        open(81,status="old",file="bonddata.txt", access="append")
     else
        open(81,status="new",file="bonddata.txt", access="append")
     end if

 
do i=1, eltot 
   do k=1, eltot
      if(i.gt.k) cycle
      m=0
      do x=1, leng
         if(atomlist(structures,x)%name.ne.elnames(k)) cycle
         if(modeselect.eq.1) exit
         do y=1, leng*27
            if(bondlist(x,y,i,k).gt.0.001) then
               if(alistrep(structures,y)%name.ne.elnames(i)) cycle
               bondavg(i,k)=bondavg(i,k)+bondlist(x,y,i,k) 
               m=m+1
            end if
         end do
      end do
      bondavg(i,k)=(dble(bondavg(i,k))/m)
      bondavg(k,i)=bondavg(i,k)
   end do
end do
write(81,*) bondavg(1,2)
close(81)

!allocate(tempmatrix(leng,leng*27))

inquire(file="bondsfile.txt", exist=dir_e)
     if(dir_e) then
        open(81,status="old",file="bondsfile.txt", access="append")
     else
        open(81,status="new",file="bondsfile.txt", access="append")
     end if
do x=1, eltot 
   do y=1, eltot
      if(x.gt.y) cycle
      do i=1, leng 
         m=0
         do j=1, leng*27
            if(bondlist(i,j,x,y).gt.0.001) then;
               m=m+1
               write(81,*) bondlist(i,j,x,y),i,j,x,y
               !print*, "!!"
               !tempmatrix(1,m)=bondlist(i,j,x,y)
            end if
         end do
      !print*, minval(tempmatrix,MASK=tempmatrix.gt.0.01)
      !deallocate(tempmatrix)
      end do

   end do
end do
close(81)
!! The folowing calculates the difference between all bonds and their corresponding average 
!! If that value is greater than a tolerance, the bonds are deemed too varied.
!! This section is primed for deletion, as maybe unimportant
CALL invar(10,tmpdig,tmpels)
do i=1, eltot
   do k=1, eltot 
      do x=1, leng
         if(modeselect.eq.1) exit
         do y=1, leng*27
            q=abs(bondlist(x,y,i,k)-bondavg(i,k))
            if(q.gt.dble(tmpdig(1)/100.0)*bondavg(i,k)) then 
               close(structures+10000) 
               print*, "Terminating. Bond lengths too varied"
               deallocate(tmpdig)
               cycle BIGLOOP
            end if
         end do
      end do
   end do
end do
q=0
deallocate(tmpdig)


do z=1, eltot 
   do l=1, eltot
      do i=1, leng*27 
         if(modeselect.eq.1) exit
         do k=1, leng*27
            do j=1, leng 
               if(i.eq.k) cycle
               if(alistrep(structures,i)%name.ne.elnames(z)) cycle
               if(alistrep(structures,j)%name.ne.elnames(l)) cycle
               if(bondlength(alistrep(structures,i)%position,atomlist(structures,j)%position).gt.bondavg(z,l)*1.2) cycle
               if(bondlength(alistrep(structures,k)%position,atomlist(structures,j)%position).gt.bondavg(z,l)*1.2) cycle
               if(bondlength(alistrep(structures,i)%position,atomlist(structures,j)%position).lt.0.001) cycle
               if(bondlength(alistrep(structures,k)%position,atomlist(structures,j)%position).lt.0.001) cycle
               
               !!NEEDS CORRECTLY IMPLEMENTING WITH MULTIPLE SPECIES

               !print*, bondangle(alistrep(structures,i)%position,&
               !     &atomlist(structures,j)%position,alistrep(structures,k)%position)
               !if(bondangle(alistrep(structures,i)%position,&
               !     &atomlist(structures,j)%position,alistrep(structures,k)%position).lt.10) then 
               !   close(structures+10000) 
               !   print*, "Terminating, tiny angles in system"
               !   cycle BIGLOOP
               !end if
               
               !if(bondangle(alistrep(structures,i)%position,&
               !     &atomlist(structures,j)%position,alistrep(structures,k)%position).gt.170) then 
               !   close(structures+10000) 
               !   print*, "Terminating, big angles in system"
               !   cycle BIGLOOP
               !end if
            end do
         end do
      end do
   end do
end do

!!!THIS SECTION ASSUMES AN AVERAGE COORDINATION NUMBER. THIS IS LIKELY VERY UNPHYSICAL!!! ALSO ASSUMES 2.5A BONDLENGTH
!k=0
!do i=1, leng   
! if(modeselect.eq.1) exit
!   do j=1, leng*27
!      if(bondlength(atomlist(structures,i)%position,alistrep(structures,j)%position).gt.(bondavg+0.2)) cycle 
!     k=k+1
!  end do
!end do
!k=nint(dble(k/leng))
!m=0

!! HYPER SPECIFIC COORDINATION SECTION 
!if(k.ne.4) then
!   close(structures+10000) 
!   print*, "Terminating. Coordination number incorrect"
!   cycle bigloop
!end if
!!!COORDINATION NUMBER NEEDS TO BE MORE SOPHISTICATED.
!do i=1, leng   
! if(modeselect.eq.1) exit
!   do j=1, leng*27
!      if(bondlength(atomlist(structures,i)%position,alistrep(structures,j)%position).gt.(bondavg+0.2)) cycle 
!      m=m+1
!   end do
  !print*, m, k
!  if(m.gt.k+1) then 
!     close(structures+10000) 
!     print*, "Terminating. Coordination number incorrect"
!     cycle bigloop
!  end if
!  if(m.lt.k-1) then 
!     close(structures+10000) 
!     cycle bigloop
!  end if
!  m=0
!end do
!print*, q


if(q.lt.0.0001) then 
call poswrite(formula(structures)%cell,atomlist,leng, structures, structno, prev_structures)

write(tmp,'(A11,I0.3)')"pos/POSCAR_",structures+prev_structures
call Incarwrite(adjustl(tmp),500, 20*leng)
call Jobwrite(tmp,3,3,3)

write(11,*) "For structure number", structures+prev_structures
write(11,*) "The average bond value is", bondavg
write(11,*) "The lower bound for allowed bonds is", bondavg-0.2
write(11,*) "The upper bound for allowed bonds is", bondavg+0.2
call potwrite(tmp, elnames, eltot)
close(structures+10000)
structures=structures+1

end if
close(structures+10000)



end do BIGLOOP
end subroutine generation


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



subroutine densitylineplot2(spacelist, formula, atomlist, alistrep, nbin, len, sigma, nbin2, nbinf, structno, atomlistt)
  type(unitcell), dimension(:), allocatable :: formula
  double precision, dimension(3,3) :: box
  double precision, dimension(3) :: test,spacelist
  integer :: i,j,k,m,len, n, nbin, nbin2,x,y,z, nbinf, structno, structures, prev_structures,p
  type (atom), dimension(:,:), allocatable :: atomlist, alistrep, atomlistt
  double precision :: r, meanvol, alpha, beta, gamma,normvol,posneg, pi,q,a,b,c, sigma
  integer, dimension(:), allocatable :: seed
  double precision, dimension(3) :: angle
  type (densitymatrix), dimension(:), allocatable :: density
  character(1024) :: name, tmp, command
  logical :: dir_e
  allocate(atomlistt(structno,len))
  prev_structures=structurecounter("don")
  open(71, file='prevstructures.txt')
  write(71, '(I0.3)') prev_structures
  write(71, '(I0.3)') structno
  close(71)


  
  do structures=1, structno
     
     do i=1, len
        allocate(atomlist(structures,i)%radden(2,nbinf*nbin2))!**3))
        atomlist(structures,i)%radden=0
        allocate(atomlistt(structures,i)%radden(2,nbinf*nbin2))!**3))
        atomlistt(structures,:)=atomlist(structures,:)
     end do
     

     write(name, '(A8,I0.3,A4)') "don/DON_",structures+prev_structures,"/DON"
     write(tmp,'(A8,I0.3)') "don/DON_",structures+prev_structures
     inquire(file=tmp, exist=dir_e)
     if(dir_e) then 
     else 
        write(command,'(A14,I0.3)')"mkdir don/DON_",structures+prev_structures
        call execute_command_line(command)
     end if
     
     
     open(structures+1000, file=name)
     
     pi=3.14159265358979323846
     j=1
     do k=1, len
        do n=1, len*27
           if (n.eq.k) cycle
           do i=1, nbin2*nbinf
              !atomlist(structures,k)%radden(1,i)=(i-1)*(maxval(formula(1)%cell)/nbin2)!
              !!WE ARE HERE CHANGIN TO A 10Angstorm/Nbin circle radius!!
              atomlist(structures,k)%radden(1,i)=(i-1)*(10.0/nbin2)
             ! if(bondlength(atomlist(structures,k)%position,alistrep(structures,n)%position)&
              !     &.ge.((i+0.5)*maxval(formula(1)%cell)/nbin2)) cycle
              !print*, dble(i+0.5)*dble(5.0/nbin2), (dble(i-0.5)*dble(5.0/nbin2)),bondlength(atomlist(structures,k)%&
              !&position,alistrep(structures,n)%position)
              if(bondlength(atomlist(structures,k)%position,alistrep(structures,n)%position)&
                   &.ge.(dble(i+0.5)*dble(5.0/nbin2))) then 
              cycle
           end if
              !if(bondlength(atomlist(structures,k)%position,alistrep(structures,n)%position)&
              !     &.lt.((i-0.5)*maxval(formula(1)%cell)/nbin2)) cycle
              if(bondlength(atomlist(structures,k)%position,alistrep(structures,n)%position)&
                   &.lt.(dble(i-0.5)*dble(5.0/nbin2))) then 
              cycle
           end if
              atomlist(structures,k)%radden(2,i)=atomlist(structures,k)%radden(2,i)+1.0
              !print*, atomlist(structures,k)%radden(2,i)
              !print*,atomlist(k)%radden(2,j)
              !print*, "larger"   
           end do
        end do
     end do
     !print*, atomlist(structures,1)%radden(2,:)
     atomlistt=atomlist
     do i=1, nbin2*nbinf
        do j=1, nbin2*nbinf
           do k=1, len
              if(k.eq.n) cycle 
              if(i.eq.j) cycle
!              print*, atomlistt(structures,k)%radden(2,i)
              atomlistt(structures,k)%radden(2,i)=atomlistt(structures,k)%radden(2,i)+atomlist(structures,k)%radden(2,j)*&
                   &exp(-((atomlist(structures,k)%radden(1,i)-atomlist(structures,k)%radden(1,j))**2)/(2*sigma**2))
              !print*, atomlistt(structures, k)%radden(2,i)
              !!!Shifts down far field into lower sig
              !print*, dble((i*5.0)/nbin2)
              
              !if((dble((i*5.0)/nbin2).gt.2.0)) &
              !     &atomlistt(structures,k)%radden(2,i)=(atomlistt(structures,k)%radden(2,i))/(dble(i*5.0/nbin2)**3)
              !print*, atomlistt(structures,k)%radden(2,i)
           end do
        end do
     end do
     do k=1, len
        do i=1, nbin2*nbinf
           
           !if((dble((i*5.0)/nbin2).gt.3.0)) then 
           !   atomlistt(structures,k)%radden(2,i)=(atomlistt(structures,k)%radden(2,i))/(dble(i*5.0/nbin2)**3)
           !   write(structures+1000,*) atomlistt(structures,k)%radden(:,i)
           !else
              write(structures+1000, *) atomlistt(structures,k)%radden(:,i)
           !end if
        end do
        !write(structures+1000, *)
        write(structures+1000, *) 
     end do
     close(structures+1000)
  end do

end subroutine densitylineplot2

subroutine lineplotcomparison(atomlistt, nbin2, formula, structno, nbinf, len, sigma, pi)
  type(unitcell), dimension(:), allocatable :: formula
  integer :: i,j,k,m,len, n, nbin, nbin2,x,y,z, nbinf, structno, structures, io,l, jp, p
  type (atom), dimension(:,:), allocatable :: alistrep, atomlistt, atomlist
  double precision ::  pi,q,a,b,c, sigma,r, tmp
  logical :: file_exists
  character(1024) :: name 
  
  name='hello'
  
  !allocate(r(nbin2*nbinf))
  n=1
  open(96, file="test")
   
  do i=1, structurecounter("don")
     do k=1, len 
        jp=j
        j=1
        write(name,'(A8,I0.3,A4)') "don/DON_", i,"/DON"
        open(1000+i, file=name, status='old')
        read(1000+i, *) r
        do 
           read(1000+i,*, iostat=io) 
           if (io/=0) exit 
           j=j+1           
        end do
        if(i.ne.1) then
           if(jp.ne.j) print*, "Warning your lengths are inconsistent"
        end if
        close(1000+i)
     end do
     write(6, '(A1, 40X, A1)', advance='no') achar(13),  achar(13)

     write(6,'(A13,I0.3,A10)', advance='no')"Process 1 is ",nint(dble(100.0*i/structurecounter("don"))), "% complete"
  end do
  
  do i=1, structno 
     do l=1, structno
        atomloop:   do m=1, len
           atomloop2: do p=1, len
           r=0
           if(p.gt.m) cycle 
           if(i.eq.l) cycle
           if(i.gt.l) cycle 
           do k=1, nbin2*nbinf
              r=r+5.0*min(atomlistt(i,m)%radden(2,k),(atomlistt(l,p)%radden(2,k)))/nbin2
           end do
           r=1-(r/(sqrt(2.0*pi)*27*len*sigma))
           write(96, *) i,r
           if(r.gt.0.2) then 
              exit Atomloop
           else if (m.eq.len) then 
              print*, "I have identified two similar structures",&
                   &structurecounter("don")-structno+i,structurecounter("don")-structno+l, m
           end if

        end do atomloop2
     end do atomloop
     end do
  end do
  
  print*, "The analysis of the created structures reveals the above matches"
  
  if(structurecounter("don").le.structno) then 
     print*, "No structures exist in databse"
  else
     allocate(atomlist(1,len)) 
     do i=1, len
        allocate(atomlist(1,i)%radden(2,j))
     end do

     do i=1, structurecounter("don")-structno
        write(name,'(A8,I0.3,A4)') "don/DON_", i,"/DON"
        open(1000+i, file=name, status='old')

        do l=1, structno
           atomloopexistingloop:   do m=1, len
              do k=1, (j-2)/len
                 if(l.eq.1) read(1000+i, *) atomlist(1,m)%radden(:,k)
              end do
              r=0
              if(i.eq.l) cycle
              if(i.gt.l) cycle 
              do k=1,j
                 r=r+5.0*min(atomlistt(l,m)%radden(2,k),(atomlist(1,m)%radden(2,k)))/nbin2
              end do
              r=1-(r/(sqrt(2.0*pi)*27*len*sigma))
              !write(96, *) i,r
              if(r.gt.0.2) then 
                 exit Atomloopexistingloop
              else if (m.eq.len) then 
                 print*, "One of the structures generated matches a database structure", &
                      &structurecounter("don")-structno+l,i
              end if
           end do atomloopexistingloop          
           close(1000+i)
        end do
     end do
  end if
  
  print*, "Comparison of these structures with existing structures reveals the above matches"

  close(96)
end subroutine lineplotcomparison

subroutine symcalc(atomlist, nbin2, nbinf, len, structures, sigma, structno)
  integer :: i,j,k,m,len, n, nbin, nbin2,x,y,z, nbinf, structno, structures, io,l, jp
  type (atom), dimension(:,:), allocatable :: alistrep, atomlistt, atomlist 
  double precision ::  pi,q,a,b,c, sigma,r, tmp
  character(1024) :: name 
  pi=3.141592564
  open(99, file="test3") 
  a=0

  do i=1, structno
     do l=1, len
        do j=1, nbin2*nbinf
           a=a+(dble((5.0)/(nbin2))*atomlist(i,l)%radden(2,j))**2
           a=sqrt(a)
           !a=dble(a/(sqrt(2.0*pi)*len*sigma))
           
        end do
     end do
     write(99, *) i,a
  end do
end subroutine symcalc

subroutine chemread(elnames, eltot, elrad)
character(3), dimension(:), allocatable :: elnames
character(3) :: read1, read2
double precision :: r_vdw, r_cov, c1, c2
integer :: increment, i, eltot, j, Reason
double precision, dimension(:,:,:), allocatable :: elrad



open(77, file="chem.in", status="old") 
do 
   read(77, *, IOSTAT=Reason) read1, read2, r_cov, r_vdw, c1, c2
   if (Reason.gt.0) then;
      stop
     ! print*, "Something wrong with chem.in file" 
   else if(Reason.lt.0) then; 
      print*, "Done"
      exit 
   else 
      print*, trim(adjustl(read1)), " ",trim(adjustl(read2))," ", r_cov, " ", r_vdw
      do i=1, eltot 
         do j=1, eltot
            if(elnames(i).eq.trim(adjustl(read1))) then;
               if(elnames(j).eq.trim(adjustl(read2))) then; 
                  elrad(1,i,j)=r_cov
                  elrad(2,i,j)=r_vdw
                  elrad(3,i,j)=c1
                  elrad(1,j,i)=r_cov
                  elrad(2,j,i)=r_vdw
                  elrad(3,j,i)=c1
                  if(elnames(i).ne.elnames(j)) then; 
                     elrad(3,i,j)=c1
                     elrad(4,i,j)=c2
                     elrad(4,j,i)=c1
                     elrad(3,j,i)=c2

                     print*, elnames(i),c1,",",elnames(j),c2
                  end if
                  continue
               end if
            end if
            if(elnames(j).eq.trim(adjustl(read1))) then;
               if(elnames(i).eq.trim(adjustl(read2))) then;
                  elrad(1,j,i)=r_cov
                  elrad(2,j,i)=r_vdw
                  elrad(3,j,i)=c1
                  if(elnames(i).ne.elnames(j)) then;
                     elrad(3,j,i)=c1
                     elrad(4,j,i)=c2
                     elrad(4,j,i)=c1
                     elrad(3,j,i)=c2

                  end if
                  continue
                  
               end if
            end if
         end do
         
      end do
   end if
end do

print*, elrad(1,1,1), elrad(1,2,1) 
print*, elrad(1,1,2), elrad(1,2,2)


print*, elrad(1,1,1), elrad(1,2,1)


 
end subroutine chemread

subroutine atomrepeater(structures,position,array,unit,atomnumber,length)
use atomtype
implicit none 
type (atom), dimension(:,:), allocatable :: array
double precision, dimension(3) :: position 
integer :: j,atomnumber,length,x,y,z,structures,m
type(unitcell), dimension(:), allocatable :: unit
 

m=((atomnumber-1)*27)+1
do x=-1,1
   do y=-1,1
      do z=-1,1
         do j=1, 3
            !print*, array(structures,m)%position(j),position(j)
            array(structures,m)%position(j)=position(j)+&
                 &(x*unit(structures)%cell(j,1))+&
                 &(y*unit(structures)%cell(j,2))+&
                 &(z*unit(structures)%cell(j,3))
            !print*,   array(structures,m)%position(j)
         end do
         m=m+1
      end do
   end do
end do
end subroutine atomrepeater

subroutine generatebondfiles(structures,array,repeatedarray,eltot,stochio,atomnumber)
use atomtype
implicit none 
character(1024) :: command,name,tmp 
integer :: l,i,structures,prev_structures, eltot, atomnumber, tmpint,m,j 
type (atom), dimension(:,:), allocatable :: array, repeatedarray
integer, dimension(:), allocatable :: stochio
logical :: dir_e

inquire(file="bon",exist=dir_e)
if(dir_e) then 
   else 
      call execute_command_line("mkdir bon")
   end if

prev_structures=structurecounter("bon")

write(tmp,'(A8,I0.3)') "bon/BON_",structures+prev_structures
inquire(file=tmp, exist=dir_e)
if(dir_e) then
else
   write(command,'(A14,I0.3)')"mkdir bon/BON_",structures+prev_structures
   call execute_command_line(command)
end if




i=1
tmpint=atomnumber 
do while (i.eq.1) 
   if((tmpint-stochio(i)).gt.0) then 
      tmpint=tmpint-stochio(i)
   else 
      i=0
   end if
end do
write(name, '(A8,I0.3,A1,A,A,I0.3)') "bon/BON_",structures+prev_structures,"/",&
     &trim(adjustl(array(structures,atomnumber)%name)),"_",tmpint
open(101, file=name)
m=1
do j=1, eltot  
   write(101,*) array(structures,m)%name 
   do i=1, stochio(j)
      do l=1, 27
         
         write(101,*) bondlength(repeatedarray(structures,m)%position,&
              &array(structures,atomnumber)%position) 
         m=m+1
      end do
   end do
end do
!write(tmp,'(A8,I0.3)') "don/DON_",structures+prev_structures
!inquire(file=tmp, exist=dir_e)
!if(dir_e) then
!else
!   write(command,'(A14,I0.3)')"mkdir don/DON_",structures+prev_structures
!   call execute_command_line(command)
!end if


end subroutine generatebondfiles

 subroutine addposcar()
    type(unitcell), dimension(:), allocatable :: formula
    character(1024) :: tmp, name
    character(1024), dimension(:), allocatable :: elnames 
    double precision :: cellmultiplier
    integer :: k,j,i,structno, ecount, eltot, leng,structures, prev_structures
    type (atom), dimension(:,:), allocatable :: atomlist, alistrep
    integer, dimension(:), allocatable :: elno
    !! Wipes the randomly generated formula
    structures=1
    structno=1
    prev_structures=structurecounter("pos")
    print*, prev_structures
    call touchpos()
    call touchposdir(structures,prev_structures)
    allocate(formula(structno))

    write(6,*) "Please enter the filename you wish to add to the database"
    !read(*, *) name
    write(name,"(A6)") "POSCAR"
    open(50, file=name)
    print*, "step 1"
    read(50, '(A)') tmp
    read(50, '(F16.0)') cellmultiplier

    do i=1, 3
       read(50, * ) formula(1)%cell(1,i),&
            &formula(1)%cell(2,i),formula(1)%cell(3,i)
    end do
    write(name,'(A11,I0.3,A7)')"pos/POSCAR_",structures+prev_structures,"/POSCAR"
    open(structures+10000, file=name,status="new")
    write(structures+10000,*) "Test"
    write(structures+10000,*) 1.0
    do i=1, 3
       write(structures+10000,*)  formula(1)%cell(1,i),formula(1)%cell(2,i),formula(1)%cell(3,i)
    end do
    read(50,'(A)') tmp
    
    ecount=0
    eltot=0
    allocate(atomlist(1,1000))
    allocate(elnames(1024))
    do i=1,len(tmp)
       if((scan(tmp(i:i+1)," ").eq.0).or.&
            &((scan(tmp(i:i+1)," ").eq.2).and.(scan(tmp(i-1:i)," ")&
            &.eq.1))) then
          eltot=eltot+1
          !print*, tmp(i:i+1)
          if(scan(tmp(i:i+1)," ").eq.0) then
             elnames(eltot)=tmp(i:i+1)
          end if
          if((scan(tmp(i:i+1)," ").eq.2).and.(scan(tmp(i-1:i)," ").eq.1)) then
             elnames(eltot)=tmp(i:i)
          end if
       else
       end if
    end do
    
    allocate(elno(eltot))
    read(50,'(4X)', advance='no')
    k=0
    read(50,*) elno
    print*, elno(:)
    read(50, *) tmp
    k=0
    do i=1, eltot 
       do j=1, elno(i)
          k=k+1
          atomlist(structures,k)%name=elnames(i)
       end do
    end do
    print*, k

    do i=1, eltot
       do j=1, elno(i)
          read(50, *)&
               &atomlist(structures,j)%position(1),atomlist(structures,j)%position(2),&
               &atomlist(structures,j)%position(3)
          print*, atomlist(structures,j)%position(:)
          !print*, atomlist(1,j)%position(1),atomlist(1,j)%position(2),atomlist(1,j)%position(3\
          
       end do
    end do
    
    call poswrite(formula(structures)%cell,atomlist,k,1,1,prev_structures)
    allocate(alistrep(1,k*27))
    do i=1, k
       call atomrepeater(structures,atomlist(structures,i)%position,alistrep,&
            &formula,i,k)
    end do
    do i=1, k
       call generatebondfiles(structures,atomlist,alistrep,eltot,elno,i)
    end do
    
  end subroutine addposcar

end module gen
