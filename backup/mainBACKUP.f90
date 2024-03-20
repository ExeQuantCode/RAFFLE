PROGRAM random
use help
use gen
use atomtype
implicit none 
type(unitcell), dimension(:), allocatable :: formula
double precision, dimension(3,3) :: box
double precision, dimension(3) :: test,spacelist
integer :: i,j,k,m,len, n, clock, reason, nbin, nbin2,x,y,z, nbinf, structno, structures, options, eltot, coordination
integer, dimension(:), allocatable :: stochio
type (atom), dimension(:,:), allocatable :: atomlist, alistrep, atomlistt
double precision :: r, meanvol, alpha, beta, gamma,normvol,posneg, pi,q,a,b,c, sigma, bond_test, returned_val 
double precision , dimension(:), allocatable :: cutoff
integer, dimension(:), allocatable :: seed 
double precision, dimension(3) :: angle, x1,x2,x3,x4
type (densitymatrix), dimension(:), allocatable :: density
character(3), dimension(:), allocatable :: elnames
double precision, dimension(:,:,:), allocatable :: elrad
character(1024) :: buffer, command 
!!For input 
character(1024), dimension(:), allocatable :: tmpels
integer, dimension(:), allocatable :: tmpdig
character(1024) :: dummy

!!Is Pie
pi=3.14159265358979323846


!!! Reads Input file !!! 

!CALL angledistribution("C  ")
 
!x1=0.0
!x2=0.0
!x3=0.0
! x4=0.0
! x1(1)=1.0
! x3(2)=1.0
! x4(3)=1.0
! print*, x1
! print*, x2
! print*, x3
! print*, fourbody(x1,x2,x3,x4)
! stop




!! This needs to be activated for a full run that reads in new samples
!call bond_evolution()

CALL invar(1,tmpdig,tmpels)
! GENERATE THIS MANY STRUCTURES
structno=tmpdig(1)
deallocate(tmpdig)
structno=structno!*7
! THINK THIS IS DEFUNCT
nbin=100
! RESOLUTION OF DON
nbin2=200
! RESOLUTOIN OF REPEATED DON
nbinf=2
! WIDTH OF GAUSSIAN FIT TO DON
sigma=0.2
! SPECIES NUMBER
CALL invar(2,tmpdig,tmpels) 
eltot=tmpdig(1)
deallocate(tmpdig)
! SPECIAL CASES
CALL invar(3,tmpdig,tmpels)
options=tmpdig(1)
deallocate(tmpdig)





!!!! 0) Run RSS
!!!! 1) Regenerate DIst Files (WIP)
!!!! 2) Run HOST_RSS
!!!! 3) Test
!!!! 4) Sphere_Overlap
!!!! 5) Bondangle_test 
!!!! 6) Run evo (Should be run after any set created)
!!!! 7) Add new poscar  
!!!! 8) Run evo, but don't regen energies or evolve distributions (only reformat gaussians) 
!!!! 9) Run evo, don't get energies but do evolve distributions
allocate(elnames(eltot))
  
!!!---------------------------------------------------------------------------------------!
!!!Assign the elements of each atom                                                       !
!!!--------------------------------------------------------------------------------------!
  

  CALL invar(4,tmpdig,tmpels)
  do i=1, eltot 
     elnames(i)=tmpels(i)
  end do
  deallocate(tmpels)
  allocate(stochio(eltot))
  !! How many of each would you like
  CALL invar(5,tmpdig,tmpels)
  do i=1, eltot
     stochio(i)=tmpdig(i)
  end do
  deallocate(tmpdig)
  len=0
  !! Total number of atoms 
  do i=1, eltot
     len=len+stochio(i)
  end do
  


!!FIX here
if(options.eq.7) then 
   call addposcar(0,dummy,1,0)
   !call addxyzfile()
   stop
end if


if(options.eq.1) then
   call regenerate_distribution_files (800)
stop   
end if

if (options.eq.3) then
   bond_test=1.430
   print*, bond_test
   CALL evaluate_contribution ("C  ","C  ",bond_test,returned_val)
   stop
end if



if(options.eq.4) then
   CALL chemread(elnames,eltot,elrad)
   print*, sphereoverlap(2*elrad(1,1,1),elrad(1,1,1),elrad(1,1,1),pi)
   print*, (4.0/3.0)*pi*elrad(1,1,1)**3
   stop
end if
if(options.eq.5) then 
   x1=0.0
   x1(1)=-1.0
   x2=0.0
   x2(2)=0.0
   x3=0.0
   x3(3)=1.0
   print*, bondangle(x1,x2,x3)
   stop
end if
if (options.eq.6) then 
   call bond_evolution(1)
   stop 
end if
if (options.eq.8) then 
   call bond_evolution(0)
   stop
end if
if (options.eq.9) then 
   call bond_evolution(2) 
   stop 
end if

!CALL chemread(elnames,eltot,elrad)
allocate(formula(structno))

!!!-----------------------------!!!
!!!Random number initialisation !!!
!!!-----------------------------!!!
CALL RANDOM_SEED(size=n)
ALLOCATE(seed(n))
CALL SYSTEM_CLOCK(COUNT=clock)
seed = clock + 37 * (/ (i - 1, i = 1, n) /)
CALL RANDOM_SEED(PUT = seed)
DEALLOCATE(seed)

!!!--------------------------------------------------!!!
!!!Set the number of atoms and generate the unit cell!!!
!!!--------------------------------------------------!!!



call generation(len, atomlist, alistrep, spacelist, formula, structno,options, eltot, elnames, stochio, elrad)
print*, "The structures requested have been successfully generated and saved"

end PROGRAM random

