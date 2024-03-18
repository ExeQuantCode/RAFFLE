program raffle
  use constants, only: real12, pi
  use geom, only: get_sphere_overlap
  use inputs
  use help
  use gen
  use atomtype
  implicit none


  integer :: i, len, nbin, nbin2, nbinf
  real(real12), dimension(3) :: x1, x2, x3

  type(unitcell), dimension(:), allocatable :: formula
  real(real12), dimension(3) :: spacelist
  type (atom), dimension(:,:), allocatable :: atomlist, alistrep
  real(real12) :: sigma, bond_test, returned_val
  real(real12), dimension(:,:,:), allocatable :: elrad
  !!For input
  character(1024) :: dummy

!!!Is Pie
  !pi=3.14159265358979323846


!!! Reads Input file !!! 

  !call angledistribution("C  ")

  !x1=0.0
  !x2=0.0
  !x3=0.0
  ! x4=0.0
  ! x1(1)=1.0
  ! x3(2)=1.0
  ! x4(3)=1.0
  ! write(*,*) x1
  ! write(*,*) x2
  ! write(*,*) x3
  ! write(*,*) get_dihedral_angle(x1,x2,x3,x4)
  ! stop


  call set_global_vars()
  
  !! This needs to be activated for a full run that reads in new samples
  !call bond_evolution()
  
  nbin=100
  ! RESOLUTION OF DON
  nbin2=200
  ! RESOLUTOIN OF REPEATED DON
  nbinf=2
  ! WIDTH OF GAUSSIAN FIT TO DON
  sigma=sigma_don





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
!!! Assign the elements of each atom                                                       !
!!!--------------------------------------------------------------------------------------!


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
     call regenerate_distribution_files (1)
     stop   
  end if

  if (options.eq.3) then
     bond_test=1.430
     write(*,*) bond_test
     call evaluate_contribution ("C  ","C  ",bond_test,returned_val)
     stop
  end if



  if(options.eq.4) then
     call chemread(elnames,eltot,elrad)
     write(*,*) get_sphere_overlap(2*elrad(1,1,1),elrad(1,1,1),elrad(1,1,1))
     write(*,*) (4.0/3.0)*pi*elrad(1,1,1)**3
     stop
  end if
  if(options.eq.5) then 
     x1=0.0
     x1(1)=-1.0
     x2=0.0
     x2(2)=0.0
     x3=0.0
     x3(3)=1.0
     write(*,*) get_bondangle(x1,x2,x3)
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

  !call chemread(elnames,eltot,elrad)
  allocate(formula(structno))


!!!--------------------------------------------------!!!
!!!Set the number of atoms and generate the unit cell!!!
!!!--------------------------------------------------!!!



  call generation(len, atomlist, alistrep, spacelist, formula, structno,options, eltot, elnames, stochio, elrad, c_cut, c_min)
  write(*,*) "The structures requested have been successfully generated and saved"

end program raffle