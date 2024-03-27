program raffle
  use constants, only: real12, pi
  use inputs
  use read_structures, only: bond_evolution
  use gen, only: generation
  use rw_geom, only: bas_type
  use evolver, only: gvector_container_type
  implicit none

  integer :: i, nbin, nbin2, nbinf
  type(bas_type) :: bas
  type(gvector_container_type) :: gvector_container
  real(real12), dimension(3) :: method_probab


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


  select case(task)
  case(0)
     write(*,*) "NOTHING WAS EVER SET UP FOR CASE 0"
  case(1)
     write(*,*) "Regenerating Distribution Files"
     write(*,*) "DEPRECATED"
     stop 0
  case(2)
     write(*,*) "Running HOST_RSS"
  case default
     write(*,*) "Invalid option"
     stop 1
  end select
  gvector_container = bond_evolution(input_dir="database/", element_file="elements.dat", element_list = element_list)

  call gvector_container%write_2body(file="2body.txt")
  call gvector_container%write_3body(file="3body.txt")
  call gvector_container%write_4body(file="4body.txt")


  !!--------------------------------------------------------------------------
  !! calculate the probability of each placement method
  !!--------------------------------------------------------------------------
  method_probab(1) = vps_ratio(1)/real(sum(vps_ratio),real12)
  method_probab(2) = method_probab(1) + &
       vps_ratio(2)/real(sum(vps_ratio),real12)
  method_probab(3) = method_probab(2) + &
       vps_ratio(3)/real(sum(vps_ratio),real12)
  write(*,*) method_probab


!!! READ THE GVECTORS IN USING bond_evolution FUNCTION.
!!! CALL generation AND PROVIDE THE GVECTORS TO GENERATE THE STRUCTURES
!!! change generation to stop using isolated calculation setup and just use the elements_database


  write(*,*) "Generating the structures requested"
  !!!--------------------------------------------------!!!
  !!!Set the number of atoms and generate the unit cell!!!
  !!!--------------------------------------------------!!!
  call generation( gvector_container, num_structures, task, &
       element_list, stoichiometry_list, &
       method_probab )
  write(*,*) "The structures requested have been successfully generated and saved"

end program raffle