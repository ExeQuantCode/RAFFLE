program test_evolver
  use evolver, only: &
       gvector_container_type
  use constants, only: real12, pi
  use rw_geom, only: basis_type
  implicit none

  logical :: success = .true.
  type(basis_type) :: basis_diamond, basis_graphite

  ! diamond cell
  basis_diamond%nspec = 1
  basis_diamond%natom = 8
  allocate(basis_diamond%spec(basis_diamond%nspec))
  basis_diamond%spec(1)%num = 8
  basis_diamond%spec(1)%name = 'C'
  allocate(basis_diamond%spec(1)%atom(basis_diamond%spec(1)%num, 3))
  basis_diamond%spec(1)%atom(1, :3) = [0.0, 0.0, 0.0]
  basis_diamond%spec(1)%atom(2, :3) = [0.5, 0.5, 0.0]
  basis_diamond%spec(1)%atom(3, :3) = [0.5, 0.0, 0.5]
  basis_diamond%spec(1)%atom(4, :3) = [0.0, 0.5, 0.5]
  basis_diamond%spec(1)%atom(5, :3) = [0.25, 0.25, 0.25]
  basis_diamond%spec(1)%atom(6, :3) = [0.75, 0.75, 0.25]
  basis_diamond%spec(1)%atom(7, :3) = [0.75, 0.25, 0.75]
  basis_diamond%spec(1)%atom(8, :3) = [0.25, 0.75, 0.75]

  basis_diamond%lat(1,:) = [3.5607451090903233, 0.0, 0.0]
  basis_diamond%lat(2,:) = [0.0, 3.5607451090903233, 0.0]
  basis_diamond%lat(3,:) = [0.0, 0.0, 3.5607451090903233]
  basis_diamond%energy = -72.213492

  ! graphite cell
  basis_graphite%nspec = 1
  basis_graphite%natom = 4
  allocate(basis_graphite%spec(basis_graphite%nspec))
  basis_graphite%spec(1)%num = 4
  basis_graphite%spec(1)%name = 'C'
  allocate(basis_graphite%spec(1)%atom(basis_graphite%spec(1)%num, 3))
  basis_graphite%spec(1)%atom(1, :3) = [0.0, 0.0, 0.25]
  basis_graphite%spec(1)%atom(2, :3) = [0.0, 0.0, 0.75]
  basis_graphite%spec(1)%atom(3, :3) = [1.0/3.0, 2.0/3.0, 0.25]
  basis_graphite%spec(1)%atom(4, :3) = [2.0/3.0, 1.0/3.0, 0.75]

  basis_graphite%lat(1,:) = [1.2336456308015413, -2.1367369110836267, 0.0]
  basis_graphite%lat(2,:) = [1.2336456308015413,  2.1367369110836267, 0.0]
  basis_graphite%lat(3,:) = [0.0, 0.0, 7.8030730000000004]
  basis_graphite%energy = -36.86795585


  call test_init_gvector_container(success)
  call test_set_width(success)
  call test_set_sigma(success)
  call test_set_cutoff_min(success)
  call test_set_cutoff_max(success)
  call test_set_radius_distance_tol(success)
  call test_create([basis_graphite], success)
  call test_update([basis_graphite, basis_diamond], success)
  !  call test_write_read(success)
  !  call test_write_2body(success)
  !  call test_write_3body(success)
  !  call test_write_4body(success)
  call test_add(basis_diamond, success)
  call test_get_element_energies(basis_diamond, success)

  !-----------------------------------------------------------------------------
  ! check for any failed tests
  !-----------------------------------------------------------------------------
  write(*,*) "----------------------------------------"
  if(success)then
     write(*,*) 'test_evolver passed all tests'
  else
     write(0,*) 'test_evolver failed one or more tests'
     stop 1
  end if


contains

  subroutine test_init_gvector_container(success)
    implicit none
    logical, intent(inout) :: success

    integer :: i
    class(gvector_container_type), allocatable :: gvector_container
    character(len=10) :: test_name

    integer, dimension(3) :: nbins
    real(real12), dimension(3) :: width, sigma, cutoff_min, cutoff_max

    ! Test case 1: Default initialisation
    gvector_container = gvector_container_type()
    nbins = [-1, -1, -1]
    width = [0.025_real12, pi/64._real12, pi/64._real12]
    sigma = [0.1_real12, 0.1_real12, 0.1_real12]
    cutoff_min = [0.5_real12, 0._real12, 0._real12]
    cutoff_max = [6._real12, pi, pi]
    test_name = "Default"

    do i = 1, 2
       call assert( &
            all( gvector_container%nbins .eq. nbins ), &
            trim(test_name)//" nbins initialisation failed", &
            success &
       )
       call assert( &
            all( abs( gvector_container%width - width ) .lt. 1.E-6_real12 ), &
            trim(test_name)//" width initialisation failed", &
            success &
       )
       call assert( &
            all( abs( gvector_container%sigma - sigma ) .lt. 1.E-6_real12 ), &
            trim(test_name)//" sigma initialisation failed", &
            success &
       )
       call assert( &
            all( abs( gvector_container%cutoff_min - cutoff_min ) .lt. &
                 1.E-6_real12 &
            ), &
            trim(test_name)//" cutoff_min initialisation failed", &
            success &
       )
       call assert( &
            all( abs( gvector_container%cutoff_max - cutoff_max ) .lt. &
                 1.E-6_real12 &
            ), &
            trim(test_name)//" cutoff_max initialisation failed", &
            success &
       )

       if(i.eq.2) exit
       ! Test case 2: Custom initialisation
       nbins = [10, 20, 30]
       width = [0.05_real12, 0.1_real12, 0.15_real12]
       sigma = [0.2_real12, 0.3_real12, 0.4_real12]
       cutoff_min = [1.0_real12, 2.0_real12, 3.0_real12]
       cutoff_max = [5.0_real12, 6.0_real12, 7.0_real12]
       gvector_container = gvector_container_type( &
            nbins=nbins,  &
            width=width,  &
            sigma=sigma,  &
            cutoff_min=cutoff_min,  &
            cutoff_max=cutoff_max &
       )
       test_name = "Custom"
    end do
  
  end subroutine test_init_gvector_container

  subroutine test_set_width(success)
    implicit none
    logical, intent(inout) :: success

    type(gvector_container_type) :: gvector_container
    real(real12), dimension(3) :: width

    ! Initialise test data
    width = [0.05_real12, 0.1_real12, 0.15_real12]

    ! Call the subroutine to set the width
    call gvector_container%set_width(width)

    ! Check if the width was set correctly
    call assert( &
         all( abs( gvector_container%width - width ) .lt. 1.E-6_real12 ), &
         "Width was not set correctly", &
         success &
    )

  end subroutine test_set_width

  subroutine test_set_sigma(success)
    implicit none
    logical, intent(inout) :: success

    type(gvector_container_type) :: gvector_container
    real(real12), dimension(3) :: sigma

    ! Initialise test data
    sigma = [0.05_real12, 0.1_real12, 0.15_real12]

    ! Call the subroutine to set the width
    call gvector_container%set_sigma(sigma)

    ! Check if the width was set correctly
    call assert( &
         all( abs( gvector_container%sigma - sigma ) .lt. 1.E-6_real12 ), &
         "Sigma was not set correctly", &
         success &
    )

  end subroutine test_set_sigma

  subroutine test_set_cutoff_min(success)
    implicit none
    logical, intent(inout) :: success

    type(gvector_container_type) :: gvector_container
    real(real12), dimension(3) :: cutoff_min

    ! Initialise test data
    cutoff_min = [0.5_real12, 0.5_real12, 0.5_real12]

    ! Call the subroutine to set the cutoff_min
    call gvector_container%set_cutoff_min(cutoff_min)

    ! Check if the cutoff_min was set correctly
    call assert( &
         all( abs( gvector_container%cutoff_min - cutoff_min ) .lt. &
              1.E-6_real12 &
         ), &
         "Cutoff_min was not set correctly", &
         success &
    )

  end subroutine test_set_cutoff_min

  subroutine test_set_cutoff_max(success)
    implicit none
    logical, intent(inout) :: success

    type(gvector_container_type) :: gvector_container
    real(real12), dimension(3) :: cutoff_max

    ! Initialise test data
    cutoff_max = [6.0_real12, 6.0_real12, 6.0_real12]

    ! Call the subroutine to set the cutoff_max
    call gvector_container%set_cutoff_max(cutoff_max)

    ! Check if the cutoff_max was set correctly
    call assert( &
         all( abs( gvector_container%cutoff_max - cutoff_max ) .lt. &
              1.E-6_real12 &
         ), &
         "Cutoff_max was not set correctly", &
         success &
    )

  end subroutine test_set_cutoff_max

  subroutine test_set_radius_distance_tol(success)
    implicit none
    logical, intent(inout) :: success

    type(gvector_container_type) :: gvector_container
    real(real12), dimension(4) :: radius_distance_tol

    ! Initialise test data
    radius_distance_tol = [1.5_real12, 2.5_real12, 3.0_real12, 6.0_real12]

    ! Call the subroutine to set the radius_distance_tol
    call gvector_container%set_radius_distance_tol(radius_distance_tol)

    ! Check if the radius_distance_tol was set correctly
    call assert( &
         all( &
              abs( &
                   gvector_container%radius_distance_tol - radius_distance_tol &
              ) .lt. 1.E-6_real12 &
         ), &
         "Radius_distance_tol was not set correctly", &
         success &
    )

  end subroutine test_set_radius_distance_tol



  subroutine test_create(basis, success)
    !! Test the create subroutine of gvector_container_type
    implicit none
    logical, intent(inout) :: success
    type(basis_type), dimension(:), intent(in) :: basis

    integer :: i
    type(gvector_container_type) :: gvector_container
    type(basis_type), dimension(size(basis,1)) :: basis_list

    ! Initialise basis_list
    do i = 1, size(basis,1)
       call basis_list(i)%copy(basis(i))
    end do

    call gvector_container%set_element_energies(['C  '], [-9.027_real12])

    ! Call the create subroutine
    call gvector_container%create(basis_list, deallocate_systems=.false.)

    ! Check if the system is allocated
    call assert( &
         allocated(gvector_container%system),  &
         "system not allocated",  &
         success &
    )

    ! Check number of elements in element_info is correct
    call assert( &
         size(gvector_container%element_info, dim=1) .eq. 1,  &
         "Number of elements in element_info is incorrect",  &
         success &
    )

    ! Check symbol of the element is correct
    call assert( &
         trim(gvector_container%element_info(1)%name) .eq. 'C',  &
         "Symbol of the element is incorrect",  &
         success &
    )

    ! Check element energies are set correctly
    call assert( &
         abs( gvector_container%element_info(1)%energy + 9.027_real12 ) .lt. &
         1.E-6_real12 , &
         "element energies not set correctly",  &
         success &
    )

    ! Check if the 2-/3-/4-body distribution functions are not allocated
    call assert( &
         ( &
              allocated(gvector_container%total%df_2body) .or. &
              allocated(gvector_container%total%df_3body) .or. &
              allocated(gvector_container%total%df_4body) &
         ),  &
          "2-/3-/4-body distribution functions are allocated",  &
          success &
    )

    ! Check number of species and species pairs are correct
    call assert( &
         size(gvector_container%total%df_2body, dim=2) .eq. 1,  &
         "Number of species pairs in 2-body distribution function &
         &is incorrect",  &
         success &
    )
    call assert( &
         size(gvector_container%total%df_3body, dim=2) .eq. 1,  &
         "Number of species in 3-body distribution function &
         &is incorrect",  &
         success &
    )
    call assert( &
         size(gvector_container%total%df_4body, dim=2) .eq. 1,  &
         "Number of species in 4-body distribution function &
         &is incorrect",  &
         success &
    )

    ! Check if the 2-/3-/4-body distribution functions are not zero
    call assert( &
         any( abs( gvector_container%total%df_2body ) .gt. 1.E-6_real12 ),  &
         "2-body distribution functions are zero",  &
         success &
    )
    call assert( &
         any( abs( gvector_container%total%df_3body ) .gt. 1.E-6_real12 ),  &
         "3-body distribution functions are zero",  &
         success &
    )
    call assert( &
         any( abs( gvector_container%total%df_4body ) .gt. 1.E-6_real12 ),  &
         "4-body distribution functions are zero",  &
         success &
    )

    ! Check if the 2-/3-/4-body distribution functions are not NaN
    call assert( &
         all( .not. isnan( gvector_container%total%df_2body ) ),  &
         "2-body distribution functions are NaN",  &
         success &
    )
    call assert( &
         all( .not. isnan( gvector_container%total%df_3body ) ),  &
         "3-body distribution functions are NaN",  &
         success &
    )
    call assert( &
         all( .not. isnan( gvector_container%total%df_4body ) ),  &
         "4-body distribution functions are NaN",  &
         success &
    )

    ! Check that the maximum value of 2-/3-/4-body distribution functions is 1
    do i = 1, size(gvector_container%total%df_2body, dim=2)
       call assert( &
            abs( &
                 maxval(gvector_container%total%df_2body(:,i)) - &
                 1._real12 &
            ) .lt. 1.E-6_real12, &
            "Maximum value of 2-body distribution functions is not 1",  &
            success &
       )
    end do
    do i = 1, size(gvector_container%total%df_3body, dim=2)
       call assert( &
            abs( &
                 maxval(gvector_container%total%df_3body(:,i)) - &
                 1._real12 &
            ) .lt. 1.E-6_real12, &
            "Maximum value of 3-body distribution functions is not 1",  &
            success &
       )
       call assert( &
            abs( &
                 maxval(gvector_container%total%df_4body(:,i)) - &
                 1._real12 &
            ) .lt. 1.E-6_real12, &
            "Maximum value of 4-body distribution functions is not 1",  &
            success &
       )
    end do

    ! Check if norm is allocated and not zero
    call assert( &
         allocated(gvector_container%norm_2body) .and. &
         all( abs( gvector_container%norm_2body ) .gt. 1.E-6_real12 ),  &
         "2-body norm is not allocated or zero",  &
         success &
    )
    call assert( &
         allocated(gvector_container%norm_3body) .and. &
         all( abs( gvector_container%norm_3body ) .gt. 1.E-6_real12 ),  &
         "3-body norm is not allocated or zero",  &
         success &
    )
    call assert( &
         allocated(gvector_container%norm_4body) .and. &
         all( abs( gvector_container%norm_4body ) .gt. 1.E-6_real12 ),  &
         "4-body norm is not allocated or zero",  &
         success &
    )

    ! Call the create subroutine again
    call gvector_container%create(basis_list, deallocate_systems=.true.)

    ! Check if the system is deallocated
    call assert( &
         .not.allocated(gvector_container%system),  &
         "system not correctly deallocated",  &
         success &
    )
    
  end subroutine test_create

  subroutine test_update(basis, success)
    implicit none
    logical, intent(inout) :: success
    type(basis_type), dimension(:), intent(in) :: basis

    integer :: i
    type(gvector_container_type) :: gvector_container
    type(basis_type), dimension(size(basis,1)) :: basis_list

    ! Initialise basis_list
    do i = 1, size(basis,1)
       call basis_list(i)%copy(basis(i))
    end do

    call gvector_container%set_element_energies(['C  '], [-9.027_real12])

    ! Call the create subroutine
    call gvector_container%create([basis_list(1)], deallocate_systems=.false.)

    ! Call the update subroutine
    call gvector_container%update([basis_list(2)], deallocate_systems=.false.)

    ! Check if the system is allocated
    call assert( &
         allocated(gvector_container%system),  &
         "system not allocated",  &
         success &
    )

    ! Call the create subroutine again
    call gvector_container%update([basis_list(2)], deallocate_systems=.true.)

    ! Check if the system is deallocated
    call assert( &
         .not.allocated(gvector_container%system),  &
         "system not correctly deallocated",  &
         success &
    )
    
  end subroutine test_update

  subroutine test_add(basis, success)
    implicit none
    logical, intent(inout) :: success
    type(basis_type), intent(in) :: basis

    type(gvector_container_type) :: gvector_container

    ! Call the add subroutine
    call gvector_container%add(basis)

    ! Check number of systems is correct
    call assert( &
         size(gvector_container%system, dim=1) .eq. 1,  &
         "Number of systems is incorrect",  &
         success &
    )

    ! Check if the system information is correct
    call assert( &
         gvector_container%system(1)%energy .eq. basis%energy,  &
         "System energy is incorrect",  &
         success &
    )
    call assert( &
         gvector_container%system(1)%num_atoms .eq. basis%natom,  &
         "Number of atoms is incorrect",  &
         success &
    )

    ! Call the add subroutine
    call gvector_container%add([basis])

    ! Check number of systems is correct
    call assert( &
         size(gvector_container%system, dim=1) .eq. 2,  &
         "Number of systems is incorrect",  &
         success &
    )

    ! Check the add subroutine
    call gvector_container%add(gvector_container%system(1))

    ! Check number of systems is correct
    call assert( &
         size(gvector_container%system, dim=1) .eq. 3,  &
         "Number of systems is incorrect",  &
         success &
    )

    ! Check the add subroutine
    call gvector_container%add(gvector_container%system)

    ! Check number of systems is correct
    call assert( &
         size(gvector_container%system, dim=1) .eq. 6,  &
         "Number of systems is incorrect",  &
         success &
    )

  end subroutine test_add

  subroutine test_get_element_energies(basis, success)
    implicit none
    logical, intent(inout) :: success
    type(basis_type), intent(in) :: basis

    type(gvector_container_type) :: gvector_container
    character(len=3), dimension(:), allocatable :: elements
    real(real12), dimension(:), allocatable :: energies

    call gvector_container%set_element_energies(['C  '], [-9.027_real12])
    call gvector_container%add(basis)

    ! Call the get_element_energies subroutine
    call gvector_container%get_element_energies(elements, energies)

    ! Check if the element energies are retrieved correctly
    call assert( &
         size(elements, dim=1) .eq. 1,  &
         "Number of elements is incorrect",  &
         success &
    )
    call assert( &
         size(energies, dim=1) .eq. 1,  &
         "Number of energies is incorrect",  &
         success &
    )
    call assert( &
         trim(elements(1)) .eq. 'C',  &
         "Element symbol is incorrect",  &
         success &
    )
    call assert( &
         abs(energies(1) + 9.027_real12) .lt. 1.E-6_real12,  &
         "Element energy is incorrect",  &
         success &
    )
    
  end subroutine test_get_element_energies

!###############################################################################

  subroutine assert(condition, message, success)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    logical, intent(inout) :: success
    if (.not. condition) then
      write(0,*) "Test failed: ", message
      success = .false.
    end if
  end subroutine assert

end program test_evolver