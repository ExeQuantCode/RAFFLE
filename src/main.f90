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


!!! Reads Input file !!! 

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
!!!! 5) Bondangle_test !!! THIS LITERALLY JUST TESTS THAT THE BONDANGLE METHOD WORKS! DO NOT USE! !!!
!!!! 6) Run evo (Should be run after any set created)
!!!! 7) Add new poscar  
!!!! 8) Run evo, but don't regen energies or evolve distributions (only reformat gaussians) 
!!!! 9) Run evo, don't get energies but do evolve distributions

!!!---------------------------------------------------------------------------------------!
!!! Assign the elements of each atom                                                       !
!!!--------------------------------------------------------------------------------------!

  len = 0
  allocate(elnames(eltot))
  !! Total number of atoms 
  do i=1, eltot
     len = len + stochio(i)
  end do


  select case(options)
  case(0)
     write(*,*) "NOTHING WAS EVER SET UP FOR CASE 0"
  case(1)
     write(*,*) "Regenerating Distribution Files"
     call regenerate_distribution_files (1)
     stop 0
  case(2)
     write(*,*) "Running HOST_RSS"
     write(*,*) "NOTHING WAS EVER SET UP FOR CASE 2"
  case(3)
     write(*,*) "Testing"
     bond_test = 1.428_real12
     write(*,*) bond_test
     call evaluate_contribution ("C  ","C  ",bond_test,returned_val)
     write(*,*) "Returned value: ",returned_val
     stop 0
  case(4)
     call chemread(elnames,eltot,elrad)
     write(*,*) get_sphere_overlap(2*elrad(1,1,1),elrad(1,1,1),elrad(1,1,1))
     write(*,*) (4.0/3.0)*pi*elrad(1,1,1)**3
     stop 0
  case(5)
     write(*,*) "Bondangle test"
     write(*,*) "DEPRECATED"
     x1 = 0._real12
     x1(1) = 1._real12
     x2 = 0._real12
     x2(2) = 0._real12
     x3 = 0._real12
     x3(3) = 1._real12
     write(*,*) get_bondangle(x1,x2,x3)
     stop 0
  case(6)
      call bond_evolution(1)
      stop 0
  case(7)
     call addposcar(0,dummy,1,0)
     !call addxyzfile()
     stop 0
  case(8)
      call bond_evolution(0)
      stop 0
  case(9)
      call bond_evolution(2) 
      stop 0
  case default
     write(*,*) "Invalid option"
     stop 1
  end select
 

  !call chemread(elnames,eltot,elrad)
  allocate(formula(structno))


!!!--------------------------------------------------!!!
!!!Set the number of atoms and generate the unit cell!!!
!!!--------------------------------------------------!!!

  call generation(len, atomlist, alistrep, spacelist, formula, structno,options, eltot, elnames, stochio, elrad, c_cut, c_min)
  write(*,*) "The structures requested have been successfully generated and saved"

end program raffle