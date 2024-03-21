program raffle
  use constants, only: real12, pi
  use geom, only: get_sphere_overlap
  use inputs
  use help
  use gen
  use atomtype
  use rw_geom, only: bas_type

  use read_chem, only: get_element_radius
  implicit none


  integer :: i, nbin, nbin2, nbinf
  real(real12), dimension(3) :: x1, x2, x3

  type(unitcell), dimension(:), allocatable :: formula
  real(real12), dimension(3) :: spacelist
  type(bas_type) :: bas
  type(atom), dimension(:,:), allocatable :: atomlist, alistrep
  real(real12) :: sigma, bond_test, returned_val
  real(real12), dimension(:,:,:), allocatable :: radius_arr
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
     write(*,*) "DEPRECATED"
     stop 0
  case(4)
     write(*,*) "Sphere Overlap Test"
     write(*,*) "DEPRECATED"
     stop 0
  case(5)
     write(*,*) "Bondangle test"
     write(*,*) "DEPRECATED"
     stop 0
  case(6)
      call bond_evolution(1)
      stop 0
  case(7)
     write(*,*) "Add individual POSCAR to the database to be run and learned from"
     write(*,*) "DEPRECATED"
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


!!!--------------------------------------------------!!!
!!!Set the number of atoms and generate the unit cell!!!
!!!--------------------------------------------------!!!

  call generation(alistrep, structno, options, element_list, stoichiometry_list)
  write(*,*) "The structures requested have been successfully generated and saved"

end program raffle