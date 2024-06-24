module raffle
  use constants, only: real12, pi
  use rw_geom, only: bas_type
  use read_structures, only: get_evolved_gvectors_from_data
  use gen, only: generation
  use evolver, only: gvector_container_type
  implicit none

 ! type(gvector_container_type) :: global_gvector_container
 ! real(real12), dimension(3) :: method_probab

  private
  public :: get_gvector_evolved



 contains

  function get_gvector_evolved( &
       input_dir, &
       element_file, bond_file, element_list, file_format, &
       width, sigma, cutoff_min, cutoff_max ) result(tmp)
    implicit none
    character(len=*), intent(in) :: input_dir
    character(len=*), intent(in), optional :: element_file
    character(len=*), intent(in), optional :: bond_file
    character(3), allocatable, dimension(:), intent(in), optional :: element_list
    character(len=*), intent(in), optional :: file_format
    real(real12), dimension(:), intent(in), optional :: width
    real(real12), dimension(:), intent(in), optional :: sigma
    real(real12), dimension(:), intent(in), optional :: cutoff_min
    real(real12), dimension(:), intent(in), optional :: cutoff_max
    type(gvector_container_type) :: gvector_container

    integer :: tmp
    character(len=256) :: element_file_, bond_file_
    real(real12), dimension(3) :: width_, sigma_, cutoff_min_, cutoff_max_

    if( .not. present(element_file) ) then
      element_file_ = "elements.dat"
    end if
    if( .not. present(bond_file) ) then
      bond_file_ = "chem.in"
    end if
    if( .not. present(width) ) then
      width_ = [ 0.025_real12, pi/24._real12, pi/32._real12 ]
    end if
    if( .not. present(sigma) ) then
      sigma_ = [ 0.1_real12, 0.05_real12, 0.05_real12 ]
    end if
    if( .not. present(cutoff_min) ) then
      cutoff_min_ = [ 0.5_real12, 0._real12, 0._real12 ]
    end if
    if( .not. present(cutoff_max) ) then
      cutoff_max_ = [ 6._real12, pi, pi/2._real12 ]
    end if

    if( present(element_list))then
       gvector_container = get_evolved_gvectors_from_data( &
           input_dir    = input_dir, &
           element_file = element_file, &
           bond_file    = bond_file, &
           element_list = element_list, &
           file_format  = file_format, &
           gvector_container_template = gvector_container_type( &
                width = width_, &
                sigma = sigma_, &
                cutoff_min = cutoff_min_, &
                cutoff_max = cutoff_max_ ) &
           )
    else
       gvector_container = get_evolved_gvectors_from_data( &
           input_dir    = input_dir, &
           element_file = element_file, &
           bond_file    = bond_file, &
           file_format  = file_format, &
           gvector_container_template = gvector_container_type( &
                width = width_, &
                sigma = sigma_, &
                cutoff_min = cutoff_min_, &
                cutoff_max = cutoff_max_ ) &
           )
    end if

  end function get_gvector_evolved


  function get_predicted_structures( &
       gvector_container, &
       lattice_host, basis_host, &
       num_structures, element_list, stoichiometry_list, method_probab, &
       task ) &
       result(bases)
    implicit none
    type(gvector_container_type), intent(in) :: gvector_container
    integer, intent(in) :: num_structures
    integer, intent(in) :: task
    character(3), dimension(:), intent(in) :: element_list
    integer, dimension(:), intent(in) :: stoichiometry_list
    real(real12), dimension(3,3), intent(in) :: lattice_host
    type(bas_type) :: basis_host
    real(real12), dimension(:), intent(in) :: method_probab

    type(bas_type), dimension(num_structures) :: bases

    call generation( gvector_container, num_structures, task, &
        element_list, stoichiometry_list, &
        method_probab )
  end function get_predicted_structures

  subroutine get_energies( &
       lattice, bases, &
       energies, energies_err )
    implicit none
    type(bas_type), dimension(:), intent(in) :: bases
    real(real12), dimension(3,3), intent(in) :: lattice
    real(real12), dimension(:), intent(out) :: energies
    real(real12), dimension(:), intent(out) :: energies_err

  end subroutine get_energies

! !!!-----------------------------------------------------------------------------
! !!! read input file
! !!!-----------------------------------------------------------------------------
!   call set_global_vars()


! !!!-----------------------------------------------------------------------------
! !!! check the task and run the appropriate case
! !!!-----------------------------------------------------------------------------
! !!! OLD TASKS !!!
! !!! 0) Run RSS
! !!! 1) Regenerate DIst Files (WIP)
! !!! 2) Run HOST_RSS
! !!! 3) Test
! !!! 4) Sphere_Overlap
! !!! 5) Bondangle_test !!! THIS LITERALLY JUST TESTS THAT THE BONDANGLE METHOD WORKS! DO NOT USE! !!!
! !!! 6) Run evo (Should be run after any set created)
! !!! 7) Add new poscar  
! !!! 8) Run evo, but don't regen energies or evolve distributions (only reformat gaussians) 
! !!! 9) Run evo, don't get energies but do evolve distributions
!   select case(task)
!   case(0)
!      write(*,*) "NOTHING WAS EVER SET UP FOR CASE 0"
!   case(1)
!      write(*,*) "Regenerating Distribution Files"
!      write(*,*) "DEPRECATED"
!      stop 0
!   case(2)
!      write(*,*) "Running HOST_RSS"
!   case default
!      write(*,*) "Invalid option"
!      stop 1
!   end select


! !!!-----------------------------------------------------------------------------
! !!! read structures from the database and generate gvectors
! !!!-----------------------------------------------------------------------------
!   gvector_container = get_evolved_gvectors_from_data( &
!        input_dir    = database_list, &
!        element_file = "elements.dat", &
!        bond_file    = "chem.in", &
!        element_list = element_list, &
!        file_format  = database_format, &
!        gvector_container_template = gvector_container_type(&
!             width = width_list, &
!             sigma = sigma_list, &
!             cutoff_min = cutoff_min_list, &
!             cutoff_max = cutoff_max_list ) )

!   call gvector_container%write_2body(file="2body.txt")
!   call gvector_container%write_3body(file="3body.txt")
!   call gvector_container%write_4body(file="4body.txt")


! !!!-----------------------------------------------------------------------------
! !!! calculate the probability of each placement method
! !!!-----------------------------------------------------------------------------
!   method_probab(1) = vps_ratio(1)/real(sum(vps_ratio),real12)
!   method_probab(2) = method_probab(1) + &
!        vps_ratio(2)/real(sum(vps_ratio),real12)
!   method_probab(3) = method_probab(2) + &
!        vps_ratio(3)/real(sum(vps_ratio),real12)
!   write(*,*) "Method probabilities (void, scan, pseudorandom-walk): ", &
!        method_probab


! !!!-----------------------------------------------------------------------------
! !!! generate random structures
! !!!-----------------------------------------------------------------------------
!   write(*,*) "Generating structures"
!   call generation( gvector_container, num_structures, task, &
!        element_list, stoichiometry_list, &
!        method_probab )
!   write(*,*) "Structures have been successfully generated and saved"

end module raffle