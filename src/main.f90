program raffle
  use constants, only: real12, pi
  use inputs
  use evolve, only: bond_evolution
  use gen, only: generation
  use rw_geom, only: bas_type
  implicit none

  integer :: i, nbin, nbin2, nbinf
  type(bas_type) :: bas


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
  ! sigma=sigma_don



    !! task=1 is a special task allowing a new poscar to be added in at user specification




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
     write(*,*) "DEPRECATED"
     stop 0
  case(2)
     write(*,*) "Running HOST_RSS"
  case(6)
      call bond_evolution("database")
      stop 0
  case default
     write(*,*) "Invalid option"
     stop 1
  end select


!!!--------------------------------------------------!!!
!!!Set the number of atoms and generate the unit cell!!!
!!!--------------------------------------------------!!!

  call generation(structno, options, element_list, stoichiometry_list)
  write(*,*) "The structures requested have been successfully generated and saved"

end program raffle