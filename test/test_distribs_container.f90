program test_distribs_container
  use raffle__io_utils, only: test_error_handling
  use raffle__distribs_container, only: &
       distribs_container_type
  use raffle__constants, only: real32, pi
  use raffle__geom_rw, only: basis_type
  implicit none

  logical :: success = .true.
  type(basis_type) :: basis_diamond, basis_graphite, basis_mgo

  test_error_handling = .true.

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

  ! diamond cell
  basis_mgo%nspec = 2
  basis_mgo%natom = 8
  allocate(basis_mgo%spec(basis_mgo%nspec))
  basis_mgo%spec(1)%num = 4
  basis_mgo%spec(1)%name = 'Mg'
  allocate(basis_mgo%spec(1)%atom(basis_mgo%spec(1)%num, 3))
  basis_mgo%spec(1)%atom(1, :3) = [0.0, 0.0, 0.0]
  basis_mgo%spec(1)%atom(2, :3) = [0.5, 0.5, 0.0]
  basis_mgo%spec(1)%atom(3, :3) = [0.5, 0.0, 0.5]
  basis_mgo%spec(1)%atom(4, :3) = [0.0, 0.5, 0.5]
  basis_mgo%spec(2)%num = 4
  basis_mgo%spec(2)%name = 'O'
  allocate(basis_mgo%spec(2)%atom(basis_mgo%spec(2)%num, 3))
  basis_mgo%spec(2)%atom(1, :3) = [0.5, 0.0, 0.0]
  basis_mgo%spec(2)%atom(2, :3) = [0.0, 0.5, 0.0]
  basis_mgo%spec(2)%atom(3, :3) = [0.0, 0.0, 0.5]
  basis_mgo%spec(2)%atom(4, :3) = [0.5, 0.5, 0.5]

  basis_mgo%lat(1,:) = [4.19, 0.0, 0.0]
  basis_mgo%lat(2,:) = [0.0, 4.19, 0.0]
  basis_mgo%lat(3,:) = [0.0, 0.0, 4.19]
  basis_mgo%energy = -20.0

  call test_init_distribs_container(success)
  call test_set_width(success)
  call test_set_sigma(success)
  call test_set_cutoff_min(success)
  call test_set_cutoff_max(success)
  call test_set_radius_distance_tol(success)
  call test_create([basis_graphite], success)
  call test_update([basis_graphite, basis_diamond], success)

  call test_create([basis_diamond, basis_mgo], success)
  call test_update([basis_diamond, basis_mgo], success)
  !  call test_write_read(success)
  !  call test_write_2body(success)
  !  call test_write_3body(success)
  !  call test_write_4body(success)
  call test_add(basis_diamond, success)
  call test_get_element_energies(basis_diamond, success)
  call test_get_element_energies_staticmem(basis_diamond, success)
  call test_set_bond_radii(basis_diamond, success)
  call test_get_bin(success)

  !-----------------------------------------------------------------------------
  ! check for any failed tests
  !-----------------------------------------------------------------------------
  write(*,*) "----------------------------------------"
  if(success)then
     write(*,*) 'test_distribs_container passed all tests'
  else
     write(0,*) 'test_distribs_container failed one or more tests'
     stop 1
  end if

contains

  subroutine test_init_distribs_container(success)
    implicit none
    logical, intent(inout) :: success

    integer :: i
    class(distribs_container_type), allocatable :: distribs_container
    character(len=10) :: test_name

    integer, dimension(3) :: nbins
    real(real32), dimension(3) :: width, sigma, cutoff_min, cutoff_max

    ! Test case 1: Default initialisation
    distribs_container = distribs_container_type()
    nbins = [-1, -1, -1]
    width = [0.025_real32, pi/64._real32, pi/64._real32]
    sigma = [0.1_real32, 0.1_real32, 0.1_real32]
    cutoff_min = [0.5_real32, 0._real32, 0._real32]
    cutoff_max = [6._real32, pi, pi]
    test_name = "Default"

    do i = 1, 2
       call assert( &
            all( distribs_container%nbins .eq. nbins ), &
            trim(test_name)//" nbins initialisation failed", &
            success &
       )
       call assert( &
            all( abs( distribs_container%width - width ) .lt. 1.E-6_real32 ), &
            trim(test_name)//" width initialisation failed", &
            success &
       )
       call assert( &
            all( abs( distribs_container%sigma - sigma ) .lt. 1.E-6_real32 ), &
            trim(test_name)//" sigma initialisation failed", &
            success &
       )
       call assert( &
            all( abs( distribs_container%cutoff_min - cutoff_min ) .lt. &
                 1.E-6_real32 &
            ), &
            trim(test_name)//" cutoff_min initialisation failed", &
            success &
       )
       call assert( &
            all( abs( distribs_container%cutoff_max - cutoff_max ) .lt. &
                 1.E-6_real32 &
            ), &
            trim(test_name)//" cutoff_max initialisation failed", &
            success &
       )

       if(i.eq.2) exit
       ! Test case 2: Custom initialisation
       nbins = [10, 20, 30]
       width = [0.05_real32, 0.1_real32, 0.15_real32]
       sigma = [0.2_real32, 0.3_real32, 0.4_real32]
       cutoff_min = [1.0_real32, 2.0_real32, 3.0_real32]
       cutoff_max = [5.0_real32, 6.0_real32, 7.0_real32]
       distribs_container = distribs_container_type( &
            nbins=nbins,  &
            width=width,  &
            sigma=sigma,  &
            cutoff_min=cutoff_min,  &
            cutoff_max=cutoff_max &
       )
       test_name = "Custom"
    end do

    write(*,*) "Testing distribs_container_type initialisation error handling"
    cutoff_min = [6.0_real32, 6.0_real32, 6.0_real32]
    cutoff_max = [1.0_real32, 1.0_real32, 1.0_real32]
    distribs_container = distribs_container_type( &
         cutoff_min=cutoff_min,  &
         cutoff_max=cutoff_max &
    )
    write(*,*) "Handled error: cutoff_min > cutoff_max"
  
  end subroutine test_init_distribs_container

  subroutine test_set_width(success)
    implicit none
    logical, intent(inout) :: success

    type(distribs_container_type) :: distribs_container
    real(real32), dimension(3) :: width

    ! Initialise test data
    width = [0.05_real32, 0.1_real32, 0.15_real32]

    ! Call the subroutine to set the width
    call distribs_container%set_width(width)

    ! Check if the width was set correctly
    call assert( &
         all( abs( distribs_container%width - width ) .lt. 1.E-6_real32 ), &
         "Width was not set correctly", &
         success &
    )

  end subroutine test_set_width

  subroutine test_set_sigma(success)
    implicit none
    logical, intent(inout) :: success

    type(distribs_container_type) :: distribs_container
    real(real32), dimension(3) :: sigma

    ! Initialise test data
    sigma = [0.05_real32, 0.1_real32, 0.15_real32]

    ! Call the subroutine to set the width
    call distribs_container%set_sigma(sigma)

    ! Check if the width was set correctly
    call assert( &
         all( abs( distribs_container%sigma - sigma ) .lt. 1.E-6_real32 ), &
         "Sigma was not set correctly", &
         success &
    )

  end subroutine test_set_sigma

  subroutine test_set_cutoff_min(success)
    implicit none
    logical, intent(inout) :: success

    type(distribs_container_type) :: distribs_container
    real(real32), dimension(3) :: cutoff_min

    ! Initialise test data
    cutoff_min = [0.5_real32, 0.5_real32, 0.5_real32]

    ! Call the subroutine to set the cutoff_min
    call distribs_container%set_cutoff_min(cutoff_min)

    ! Check if the cutoff_min was set correctly
    call assert( &
         all( abs( distribs_container%cutoff_min - cutoff_min ) .lt. &
              1.E-6_real32 &
         ), &
         "Cutoff_min was not set correctly", &
         success &
    )

  end subroutine test_set_cutoff_min

  subroutine test_set_cutoff_max(success)
    implicit none
    logical, intent(inout) :: success

    type(distribs_container_type) :: distribs_container
    real(real32), dimension(3) :: cutoff_max

    ! Initialise test data
    cutoff_max = [6.0_real32, 6.0_real32, 6.0_real32]

    ! Call the subroutine to set the cutoff_max
    call distribs_container%set_cutoff_max(cutoff_max)

    ! Check if the cutoff_max was set correctly
    call assert( &
         all( abs( distribs_container%cutoff_max - cutoff_max ) .lt. &
              1.E-6_real32 &
         ), &
         "Cutoff_max was not set correctly", &
         success &
    )

  end subroutine test_set_cutoff_max

  subroutine test_set_radius_distance_tol(success)
    implicit none
    logical, intent(inout) :: success

    type(distribs_container_type) :: distribs_container
    real(real32), dimension(4) :: radius_distance_tol

    ! Initialise test data
    radius_distance_tol = [1.5_real32, 2.5_real32, 3.0_real32, 6.0_real32]

    ! Call the subroutine to set the radius_distance_tol
    call distribs_container%set_radius_distance_tol(radius_distance_tol)

    ! Check if the radius_distance_tol was set correctly
    call assert( &
         all( &
              abs( &
                   distribs_container%radius_distance_tol - radius_distance_tol &
              ) .lt. 1.E-6_real32 &
         ), &
         "Radius_distance_tol was not set correctly", &
         success &
    )

  end subroutine test_set_radius_distance_tol



  subroutine test_create(basis, success)
    !! Test the create subroutine of distribs_container_type
    implicit none
    logical, intent(inout) :: success
    type(basis_type), dimension(:), intent(in) :: basis

    integer :: i, j, k
    integer :: num_pairs
    character(len=3) :: species_tmp
    type(distribs_container_type) :: distribs_container
    type(basis_type), dimension(size(basis,1)) :: basis_list
    character(len=3), dimension(:), allocatable :: elements

    ! Initialise basis_list
    do i = 1, size(basis,1)
       call basis_list(i)%copy(basis(i))
    end do

    ! Test element_database uninitiaised error handling
    write(*,*) "Testing distribs_container_type create error handling"
    call distribs_container%create(basis_list, deallocate_systems=.false.)
    write(*,*) "Handled error: element_database not initialised"

    ! Set element energies
    allocate(elements(0))
    do i = 1, size(basis_list)
       species_loop: do j = 1, basis_list(i)%nspec
          call distribs_container%set_element_energies( &
               [ basis_list(i)%spec(1)%name ], [ -9.027_real32 ] &
          )
          species_tmp = basis_list(i)%spec(j)%name(1:3)
          if(.not.allocated(elements).or.size(elements,1).eq.0) then
             elements = [ species_tmp ]
             cycle species_loop
          end if
          do k = 1, size(elements), 1
             if ( trim(elements(k)) .eq. trim(species_tmp) ) &
                  cycle species_loop
          end do
          elements = [ elements, species_tmp ]
       end do species_loop
    end do
    num_pairs = nint( &
         gamma(real(size(elements) + 2, real32) ) / &
         ( gamma(real(size(elements), real32)) * gamma( 3._real32 ) ) &
    )

    ! Call the create subroutine
    call distribs_container%create(basis_list, deallocate_systems=.false.)

    ! Check if the system is allocated
    call assert( &
         allocated(distribs_container%system),  &
         "system not allocated",  &
         success &
    )

    ! Check number of elements in element_info is correct
    call assert( &
         size(distribs_container%element_info, dim=1) .eq. size(elements),  &
         "Number of elements in element_info is incorrect",  &
         success &
    )

    ! Check symbol of the element is correct
    call assert( &
         any( elements .eq. distribs_container%element_info(1)%name ),  &
         "Symbol of the element is incorrect",  &
         success &
    )

    ! Check element energies are set correctly
    call assert( &
         abs( distribs_container%element_info(1)%energy + 9.027_real32 ) .lt. &
         1.E-6_real32 , &
         "element energies not set correctly",  &
         success &
    )

    ! Check if the 2-/3-/4-body distribution functions are not allocated
    call assert( &
         ( &
              allocated(distribs_container%gdf%df_2body) .or. &
              allocated(distribs_container%gdf%df_3body) .or. &
              allocated(distribs_container%gdf%df_4body) &
         ),  &
         "2-/3-/4-body distribution functions are allocated",  &
         success &
    )

    ! Check number of species and species pairs are correct
    call assert( &
         size(distribs_container%gdf%df_2body, dim=2) .eq. num_pairs,  &
         "Number of species pairs in 2-body distribution function &
         &is incorrect",  &
         success &
    )
    call assert( &
         size(distribs_container%gdf%df_3body, dim=2) .eq. size(elements),  &
         "Number of species in 3-body distribution function &
         &is incorrect",  &
         success &
    )
    call assert( &
         size(distribs_container%gdf%df_4body, dim=2) .eq. size(elements),  &
         "Number of species in 4-body distribution function &
         &is incorrect",  &
         success &
    )

    ! Check if the 2-/3-/4-body distribution functions are not zero
    call assert( &
         any( abs( distribs_container%gdf%df_2body ) .gt. 1.E-6_real32 ),  &
         "2-body distribution functions are zero",  &
         success &
    )
    call assert( &
         any( abs( distribs_container%gdf%df_3body ) .gt. 1.E-6_real32 ),  &
         "3-body distribution functions are zero",  &
         success &
    )
    call assert( &
         any( abs( distribs_container%gdf%df_4body ) .gt. 1.E-6_real32 ),  &
         "4-body distribution functions are zero",  &
         success &
    )

    ! Check if the 2-/3-/4-body distribution functions are not NaN
    call assert( &
         all( .not. isnan( distribs_container%gdf%df_2body ) ),  &
         "2-body distribution functions are NaN",  &
         success &
    )
    call assert( &
         all( .not. isnan( distribs_container%gdf%df_3body ) ),  &
         "3-body distribution functions are NaN",  &
         success &
    )
    call assert( &
         all( .not. isnan( distribs_container%gdf%df_4body ) ),  &
         "4-body distribution functions are NaN",  &
         success &
    )

    ! Check that the maximum value of 2-/3-/4-body distribution functions is 1
    do i = 1, size(distribs_container%gdf%df_2body, dim=2)
       call assert( &
            abs( &
                 maxval(distribs_container%gdf%df_2body(:,i)) - &
                 1._real32 &
            ) .lt. 1.E-6_real32, &
            "Maximum value of 2-body distribution functions is not 1",  &
            success &
       )
    end do
    do i = 1, size(distribs_container%gdf%df_3body, dim=2)
       call assert( &
            abs( &
                 maxval(distribs_container%gdf%df_3body(:,i)) - &
                 1._real32 &
            ) .lt. 1.E-6_real32, &
            "Maximum value of 3-body distribution functions is not 1",  &
            success &
       )
       call assert( &
            abs( &
                 maxval(distribs_container%gdf%df_4body(:,i)) - &
                 1._real32 &
            ) .lt. 1.E-6_real32, &
            "Maximum value of 4-body distribution functions is not 1",  &
            success &
       )
    end do

    ! Check if norm is allocated and not zero
    call assert( &
         allocated(distribs_container%norm_2body) .and. &
         all( abs( distribs_container%norm_2body ) .gt. 1.E-6_real32 ),  &
         "2-body norm is not allocated or zero",  &
         success &
    )
    call assert( &
         allocated(distribs_container%norm_3body) .and. &
         all( abs( distribs_container%norm_3body ) .gt. 1.E-6_real32 ),  &
         "3-body norm is not allocated or zero",  &
         success &
    )
    call assert( &
         allocated(distribs_container%norm_4body) .and. &
         all( abs( distribs_container%norm_4body ) .gt. 1.E-6_real32 ),  &
         "4-body norm is not allocated or zero",  &
         success &
    )

    ! Call the create subroutine again
    call distribs_container%create(basis_list, deallocate_systems=.true.)

    ! Check if the system is deallocated
    call assert( &
         .not.allocated(distribs_container%system),  &
         "system not correctly deallocated",  &
         success &
    )
    
  end subroutine test_create

  subroutine test_update(basis, success)
    implicit none
    logical, intent(inout) :: success
    type(basis_type), dimension(:), intent(in) :: basis

    integer :: i, j, k
    integer :: num_pairs
    character(len=3) :: species_tmp
    type(distribs_container_type) :: distribs_container
    type(basis_type), dimension(size(basis,1)) :: basis_list
    character(len=3), dimension(:), allocatable :: elements

    ! Initialise basis_list
    do i = 1, size(basis,1)
       call basis_list(i)%copy(basis(i))
    end do

    ! Set host system
    call distribs_container%host_system%set(basis(1))

    ! Set element energies
    do i = 1, size(basis_list)
       species_loop: do j = 1, basis_list(i)%nspec
          call distribs_container%set_element_energies( &
               [ basis_list(i)%spec(1)%name ], [ -9.027_real32 ] &
          )
          species_tmp = basis_list(i)%spec(j)%name(1:3)
          if(.not.allocated(elements)) then
             elements = [ species_tmp ]
             cycle species_loop
          end if
          do k = 1, size(elements), 1
             if ( trim(elements(k)) .eq. trim(species_tmp) ) &
                  cycle species_loop
          end do
          elements = [ elements, species_tmp ]
       end do species_loop
    end do
    num_pairs = nint( &
         gamma(real(size(elements) + 2, real32)) / &
         ( gamma(real(size(elements), real32)) * gamma( 3._real32 ) ) &
    )

    ! Call the create subroutine
    call distribs_container%create([basis_list(1)], deallocate_systems=.false.)

    ! Call the update subroutine
    call distribs_container%update([basis_list(2)], deallocate_systems=.false.)

    ! Check if the system is allocated
    call assert( &
         allocated(distribs_container%system),  &
         "system not allocated",  &
         success &
    )

    ! Call the create subroutine again
    call distribs_container%update([basis_list(2)], deallocate_systems=.true.)

    ! Check if the system is deallocated
    call assert( &
         .not.allocated(distribs_container%system),  &
         "system not correctly deallocated",  &
         success &
    )
    
  end subroutine test_update

  subroutine test_add(basis, success)
    implicit none
    logical, intent(inout) :: success
    type(basis_type), intent(in) :: basis

    type(distribs_container_type) :: distribs_container
    integer, dimension(1,1,1,1) :: test_array = 1

    ! Call the add subroutine
    call distribs_container%add(basis)

    ! Check number of systems is correct
    call assert( &
         size(distribs_container%system, dim=1) .eq. 1,  &
         "Number of systems is incorrect",  &
         success &
    )

    ! Check if the system information is correct
    call assert( &
         abs( &
              distribs_container%system(1)%energy - basis%energy &
         ) .lt. 1.E-6,  &
         "System energy is incorrect",  &
         success &
    )
    call assert( &
         distribs_container%system(1)%num_atoms .eq. basis%natom,  &
         "Number of atoms is incorrect",  &
         success &
    )

    ! Call the add subroutine
    call distribs_container%add([basis])

    ! Check number of systems is correct
    call assert( &
         size(distribs_container%system, dim=1) .eq. 2,  &
         "Number of systems is incorrect",  &
         success &
    )

    ! Check the add subroutine
    call distribs_container%add(distribs_container%system(1))

    ! Check number of systems is correct
    call assert( &
         size(distribs_container%system, dim=1) .eq. 3,  &
         "Number of systems is incorrect",  &
         success &
    )

    ! Check the add subroutine
    call distribs_container%add(distribs_container%system)

    ! Check number of systems is correct
    call assert( &
         size(distribs_container%system, dim=1) .eq. 6,  &
         "Number of systems is incorrect",  &
         success &
    )

    ! Test unknown type and rank error handling
    write(*,*) "Testing distribs_container_type add error handling"
    call distribs_container%add(1)
    call distribs_container%add([1])
    write(*,*) "Handled error: system default type"
    call distribs_container%add(test_array)
    write(*,*) "Handled error: system default rank"



  end subroutine test_add

  subroutine test_get_element_energies(basis, success)
    implicit none
    logical, intent(inout) :: success
    type(basis_type), intent(in) :: basis

    type(distribs_container_type) :: distribs_container
    character(len=3), dimension(:), allocatable :: elements
    real(real32), dimension(:), allocatable :: energies

    call distribs_container%set_element_energies(['C  '], [-9.027_real32])
    call distribs_container%add(basis)

    ! Call the get_element_energies subroutine
    call distribs_container%get_element_energies(elements, energies)

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
         abs(energies(1) + 9.027_real32) .lt. 1.E-6_real32,  &
         "Element energy is incorrect",  &
         success &
    )
    
  end subroutine test_get_element_energies

  subroutine test_get_element_energies_staticmem(basis, success)
    implicit none
    logical, intent(inout) :: success
    type(basis_type), intent(in) :: basis

    type(distribs_container_type) :: distribs_container
    character(len=3), dimension(1) :: elements
    real(real32), dimension(1) :: energies

    call distribs_container%set_element_energies(['C  '], [-9.027_real32])
    call distribs_container%add(basis)

    ! Call the get_element_energies_staticmem subroutine
    call distribs_container%get_element_energies_staticmem(elements, energies)

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
         abs(energies(1) + 9.027_real32) .lt. 1.E-6_real32,  &
         "Element energy is incorrect",  &
         success &
    )
    
  end subroutine test_get_element_energies_staticmem

  subroutine test_set_bond_radii(basis, success)
    implicit none
    logical, intent(inout) :: success
    type(basis_type), intent(in) :: basis

    integer :: i
    type(distribs_container_type) :: distribs_container
    real(real32), dimension(1) :: radii
    character(len=3), dimension(1,2) :: elements
    real(real32), dimension(:), allocatable :: radii_get
    character(len=3), dimension(:,:), allocatable :: elements_get

    ! Initialise test data
    radii(1) = 12.5_real32
    elements(1,:) = ['C  ', 'C  ']

    ! Call the subroutine to set the bond radii
    call distribs_container%set_bond_radii(elements, radii)
    call distribs_container%add(basis)

    ! Get the bond radii
    call distribs_container%get_bond_radii(elements_get, radii_get)

    ! Check if the number of bond elements is correct
    call assert( &
         size(elements_get, dim=1) .eq. 1,  &
         "Number of bond elements is incorrect",  &
         success &
    )

    do i = 1, 2
       ! Check if the bond elements were set correctly
       call assert( &
            all( elements_get .eq. elements ), &
            "Bond elements were not set correctly", &
            success &
       )
       ! Check if the bond radius was set correctly
       call assert( &
            all( abs( radii_get - radii ) .lt. 1.E-6 ), &
            "Bond radius was not set correctly", &
            success &
       )

       elements_get = '   '
       radii_get = 0.0_real32
       ! Get the bond radii from staticmem
       call distribs_container%get_bond_radii_staticmem(elements_get, radii_get)
    end do

    radii(1) = 14.2_real32
    call distribs_container%set_bond_radii(elements, radii)
    call distribs_container%get_bond_radii(elements_get, radii_get)

    ! Check if the bond radius was set correctly
    call assert( &
         all( abs( radii_get - radii ) .lt. 1.E-6 ), &
         "Bond radius was not set correctly", &
         success &
    )
    
    

  end subroutine test_set_bond_radii


  subroutine test_get_bin(success)
    implicit none
    logical, intent(inout) :: success

    integer :: bin
    type(distribs_container_type) :: distribs_container

    distribs_container%nbins(1) = 10

    ! Check lower bound correct handling
    bin = distribs_container%get_bin(0._real32, 1)
    call assert( &
         bin .eq. 0,  &
         "Bin is incorrect",  &
         success &
    )

    ! Check upper bound correct handling
    bin = distribs_container%get_bin(100._real32, 1)
    call assert( &
         bin .eq. 0,  &
         "Bin is incorrect",  &
         success &
    )

    ! Check middle value correct handling
    bin = distribs_container%get_bin(3._real32, 1)
    call assert( &
         bin .eq. 1 + nint( (distribs_container%nbins(1) - 1) * 2.5 / 5.5 ),  &
         "Bin is incorrect",  &
         success &
    )


  end subroutine test_get_bin

!###############################################################################

  subroutine assert(condition, message, success)
    implicit none
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    logical, intent(inout) :: success
    if (.not. condition) then
       write(0,*) "Test failed: ", message
       success = .false.
    end if
  end subroutine assert

end program test_distribs_container