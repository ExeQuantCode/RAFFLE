module raffle__place_methods
  !! Module containing the placement methods available within RAFFLE.
  !!
  !! This module contains procedures to query points for atom placement.
  !! The available placement methods are:
  !! - void: place the atom in the gridpoint with the largest void space
  !! - rand: place the atom at a random gridpoint
  !! - walk: place the atom using a random walk method
  !! - growth: place the atom using a random walk, with last placement point
  !!        as the starting point
  !! - min:  place the atom at the gridpoint with the highest viability
  use raffle__constants, only: real32, pi
  use raffle__misc_linalg, only: modu, inverse_3x3
  use raffle__geom_extd, only: extended_basis_type
  use raffle__dist_calcs, only: &
       get_min_dist, &
       get_min_dist_between_point_and_species
  use evaluator, only: evaluate_point
  use raffle__distribs_container, only: distribs_container_type
  implicit none


  private

  public :: &
       place_method_void, place_method_rand, &
       place_method_walk, place_method_growth, &
       place_method_min


contains

!###############################################################################
  function place_method_void( &
       grid, grid_offset, basis, atom_ignore_list, viable &
   ) result(point)
    !! VOID placement method.
    !!
    !! This method returns the gridpoint with the lowest neighbour density.
    !! i.e. the point with the lowest density in the cell.
    implicit none

    ! Arguments
    type(extended_basis_type), intent(inout) :: basis
    !! Structure to add atom to.
    integer, dimension(3), intent(in) :: grid
    !! Number of gridpoints in each direction.
    real(real32), dimension(3), intent(in) :: grid_offset
    !! Offset for gridpoints.
    integer, dimension(:,:), intent(in) :: atom_ignore_list
    !! List of atoms to ignore (i.e. indices of atoms not yet placed).
    logical, intent(out) :: viable
    !! Boolean to indicate if point is viable.
    real(real32), dimension(3) :: point
    !! Point to add atom to.
    
    ! Local variables
    integer :: i, j, k
    !! Loop indices.
    real(real32), dimension(3) :: best_location
    !! Index of best location to place atom.
    real(real32) :: best_location_bond, smallest_bond
    !! Bond lengths.
    real(real32), dimension(3) :: tmpvector
    !! Temporary vector for gridpoint.

   
    !---------------------------------------------------------------------------
    ! loop over all gridpoints in the unit cell and find the one with the ...
    ! ... largest void space
    !---------------------------------------------------------------------------
    viable = .false.
    best_location = 0._real32
    best_location_bond = -huge(1._real32)
    do i = 0, grid(1) - 1, 1
       do j = 0, grid(2) - 1, 1
          do k = 0, grid(3) - 1, 1
             tmpvector = [ &
                  i + grid_offset(1), &
                  j + grid_offset(2), &
                  k + grid_offset(3) &
             ] / real(grid,real32)
             smallest_bond = modu(get_min_dist(&
                  basis, tmpvector, .false., &
                  ignore_list = atom_ignore_list))
             if( smallest_bond .gt. best_location_bond ) then
                best_location_bond = smallest_bond
                best_location = tmpvector
             end if
          end do
       end do
    end do

    ! return the gridpoint with the largest void space
    point = best_location
    viable = .true.

  end function place_method_void
!###############################################################################


!###############################################################################
  function place_method_rand( &
       basis, atom_ignore_list, radius_list, max_attempts, viable &
  ) result(point)
    !! Random placement method.
    !!
    !! This method places the atom at a random gridpoint.
    implicit none

    ! Arguments
    type(extended_basis_type), intent(inout) :: basis
    !! Structure to add atom to.
    integer, dimension(:,:), intent(in) :: atom_ignore_list
    !! List of atoms to ignore (i.e. indices of atoms not yet placed).
    real(real32), dimension(:), intent(in) :: radius_list
    !! List of radii for each pair of elements.
    integer, intent(in) :: max_attempts
    !! Limit on number of attempts.
    logical, intent(out) :: viable
    !! Boolean to indicate if point is viable.
    real(real32), dimension(3) :: point
    !! Point to add atom to.

    ! Local variables
    integer :: i, j, is, js
    !! Loop indices.
    real(real32) :: rtmp1
    !! random number.
    logical :: ltmp1
    !! logical variable.
    integer, dimension(basis%nspec,basis%nspec) :: pair_index


    viable = .false.

    i = 0
    do is = 1, basis%nspec
       do js = is, basis%nspec, 1
          i = i + 1
          pair_index(js,is) = i
          pair_index(is,js) = i
       end do
   end do

    ! find a random gridpoint
    atom_loop: do i = 1, max_attempts
       do j = 1, 3
          call random_number(rtmp1)
          point(j) = rtmp1
       end do
       do js = 1, basis%nspec
          if( &
               get_min_dist_between_point_and_species( &
                    basis, point, &
                    species = atom_ignore_list(1,1), &
                    ignore_list = atom_ignore_list &
               ) .lt. radius_list(pair_index(atom_ignore_list(1,1),js)) &
          )then
             cycle atom_loop
          end if
       end do
       exit atom_loop
    end do atom_loop

  end function place_method_rand
!###############################################################################


!###############################################################################
  function place_method_walk( distribs_container, &
       basis, atom_ignore_list, &
       radius_list, max_attempts, &
       step_size_coarse, step_size_fine, &
       viable &
  ) result(point)
    !! Random walk placement method.
    !!
    !! This method places the atom using a random walk method.
    !! An initial point is chosen at random, and then points nearby are tested
    !! to see if they are more suitable than the current point. If they are,
    !! the query point is moved to the new point and the process is repeated.
    !! The process is repeated, with each point being tested against a random
    !! number. If the random number is less than the suitability of the point,
    !! the atom is placed at that point.
    implicit none

    ! Arguments
    type(distribs_container_type), intent(in) :: distribs_container
    !! Distribution function (gvector) container.
    type(extended_basis_type), intent(inout) :: basis
    !! Structure to add atom to.
    integer, intent(in) :: max_attempts
    !! Limit on number of attempts.
    real(real32), intent(in) :: step_size_coarse, step_size_fine
    !! Step sizes for random walk.
    logical, intent(out) :: viable
    !! Boolean to indicate if point is viable.
    integer, dimension(:,:), intent(in) :: atom_ignore_list
    !! List of atoms to ignore (i.e. indices of atoms not yet placed).
    real(real32), dimension(:), intent(in) :: radius_list
    !! List of radii for each pair of elements.
    real(real32), dimension(3) :: point
    !! Point to add atom to.

    ! Local variables
    integer :: i
    !! Loop indices.
    integer :: nattempt, nstuck
    !! Number of attempts and number of times stuck at same site
    real(real32) :: rtmp1
    !! Random number.
    real(real32), dimension(3) :: rvec1, abc
    !! Random vector and lattice constants.
    real(real32) :: crude_norm
    !! Crude normalisation.
    real(real32) :: site_value, test_value
    !! Viability values.
    real(real32), dimension(3) :: site_vector, test_vector
    !! Vectors for gridpoints.
   

    viable = .false.

    !---------------------------------------------------------------------------
    ! test a random point in the unit cell
    !---------------------------------------------------------------------------
    do i = 1, 3
      abc(i) = modu(basis%lat(i,:))
    end do
    i = 0
    random_loop : do 
       i = i + 1      
       call random_number(site_vector)

       site_value = evaluate_point( distribs_container, &
            site_vector, atom_ignore_list(1,1), basis, &
            atom_ignore_list, radius_list &
       )
       call random_number(rtmp1)
       if(rtmp1.lt.site_value) exit random_loop
 
       if(i.ge.max_attempts) return
    end do random_loop


    !---------------------------------------------------------------------------
    ! now do a random walk to find a suitable point to place the atom
    !---------------------------------------------------------------------------
    nattempt = 0
    nstuck = 0
    crude_norm = 0.5_real32
    walk_loop : do
       call random_number(rtmp1)
       if(nattempt.ge.10) then 
          test_vector = site_vector + &
               ( rvec1 * 2._real32 - 1._real32 ) * step_size_fine / abc
       else
          test_vector = site_vector + &
          ( rvec1 * 2._real32 - 1._real32 ) * step_size_coarse / abc
       end if
       test_vector = test_vector - floor(test_vector)

       test_value = evaluate_point( distribs_container, &
            test_vector, atom_ignore_list(1,1), basis, &
            atom_ignore_list, radius_list &
       )
       if(test_value.lt.site_value) then 
          nstuck = nstuck + 1
          if(nstuck.ge.10) then
             nattempt = nattempt + 1
             if(crude_norm.lt.site_value) &
                  crude_norm = &
                       ( crude_norm + site_value/real(nattempt) ) / 2._real32

             ! if we have tried 10 times, and still no luck, then we need to
             ! reduce the tolerance
             if(nattempt.ge.10) site_value = site_value / crude_norm
             call random_number(rtmp1)
             if(rtmp1.lt.site_value) exit walk_loop   
          end if   
       else
          nstuck = 0
          site_vector = test_vector
          site_value  = test_value
 
          if(nattempt.ge.10) test_value = test_value / crude_norm
          call random_number(rtmp1)
          if(rtmp1.lt.test_value) exit walk_loop
       end if
 
    end do walk_loop
 
    point = site_vector
    viable=.true.
   
  end function place_method_walk
!###############################################################################


!###############################################################################
  function place_method_growth( distribs_container, &
       prior_point, prior_species, &
       basis, atom_ignore_list, &
       radius_list, max_attempts, &
       step_size_coarse, step_size_fine, &
       viable &
  ) result(point)
    !! Random walk placement method.
    !!
    !! This method places the atom using a random walk method.
    !! An initial point is chosen at random, and then points nearby are tested
    !! to see if they are more suitable than the current point. If they are,
    !! the query point is moved to the new point and the process is repeated.
    !! The process is repeated, with each point being tested against a random
    !! number. If the random number is less than the suitability of the point,
    !! the atom is placed at that point.
    implicit none

    ! Arguments
    type(distribs_container_type), intent(in) :: distribs_container
    !! Distribution function (gvector) container.
    real(real32), dimension(3), intent(in) :: prior_point
    !! Point to start walk from.
    integer, intent(in) :: prior_species
    !! Species of last atom placed.
    type(extended_basis_type), intent(inout) :: basis
    !! Structure to add atom to.
    integer, intent(in) :: max_attempts
    !! Limit on number of attempts.
    real(real32), intent(in) :: step_size_coarse, step_size_fine
    !! Step sizes for random walk.
    logical, intent(out) :: viable
    !! Boolean to indicate if point is viable.
    integer, dimension(:,:), intent(in) :: atom_ignore_list
    !! List of atoms to ignore (i.e. indices of atoms not yet placed).
    real(real32), dimension(:), intent(in) :: radius_list
    !! List of radii for each pair of elements.
    real(real32), dimension(3) :: point
    !! Point to add atom to.

    ! Local variables
    integer :: i, idx
    !! Loop indices.
    integer :: nattempt, nstuck
    !! Number of attempts and number of times stuck at same site
    real(real32) :: rtmp1, min_radius
    !! Random number and minimum radius.
    real(real32), dimension(3) :: rvec1, abc
    !! Random vector and lattice constants.
    real(real32) :: crude_norm
    !! Crude normalisation.
    real(real32) :: site_value, test_value
    !! Viability values.
    real(real32), dimension(3) :: site_vector, test_vector
    !! Vectors for gridpoints.
    real(real32), dimension(3,3) :: inverse_lattice
   

    viable = .false.

    !---------------------------------------------------------------------------
    ! get the lattice constants and the inverse lattice
    !---------------------------------------------------------------------------
    do i = 1, 3
      abc(i) = modu(basis%lat(i,:))
    end do
    inverse_lattice = inverse_3x3(basis%lat)


    !---------------------------------------------------------------------------
    ! get the index of the pair of species
    !---------------------------------------------------------------------------
    idx = nint( &
         ( &
              basis%nspec - &
              min( prior_species, atom_ignore_list(1,1) &
         ) / 2._real32 ) * &
         (  &
              min( prior_species, atom_ignore_list(1,1) ) - &
              1._real32 &
         ) +  max( prior_species, atom_ignore_list(1,1) ) &
    )
    min_radius = radius_list(idx) * distribs_container%radius_distance_tol(1)


    !---------------------------------------------------------------------------
    ! test a random point within a spherical shell around the prior point
    !---------------------------------------------------------------------------
    i = 0
    shell_loop: do
       i = i + 1
       call random_number(rvec1)
       ! map rvec1(1) to ring between min_radius and min_radius + 1.0
       rvec1(1) = rvec1(1) + min_radius ! r
       rvec1(2) = rvec1(2) * 2._real32 * pi ! theta
       rvec1(3) = rvec1(3) * pi ! phi
       ! convert from spherical to cartesian
       rvec1 = [ &
             rvec1(1) * cos(rvec1(2)) * sin(rvec1(3)), &
             rvec1(1) * sin(rvec1(2)) * sin(rvec1(3)), &
             rvec1(1) * cos(rvec1(3)) &
       ]
       ! convert from cartesian to direct
       rvec1 = matmul(rvec1, inverse_lattice)
       site_vector = prior_point + rvec1
       ! now evaluate the point and check if it passes the initial criteria
       site_value = evaluate_point( distribs_container, &
            site_vector, atom_ignore_list(1,1), basis, &
            atom_ignore_list, radius_list &
       )
       call random_number(rtmp1)
       if(rtmp1.lt.site_value) exit shell_loop

       if(i.ge.max_attempts) return
    end do shell_loop


    !---------------------------------------------------------------------------
    ! now do a random walk to find a suitable point to place the atom
    !---------------------------------------------------------------------------
    nattempt = 0
    nstuck = 0
    crude_norm = 0.5_real32
    walk_loop : do
       call random_number(rtmp1)
       if(nattempt.ge.10) then 
          test_vector = site_vector + &
               ( rvec1 * 2._real32 - 1._real32 ) * step_size_fine / abc
       else
          test_vector = site_vector + &
          ( rvec1 * 2._real32 - 1._real32 ) * step_size_coarse / abc
       end if 
       test_vector = test_vector - floor(test_vector)

       test_value = evaluate_point( distribs_container, &
            test_vector, atom_ignore_list(1,1), basis, &
            atom_ignore_list, radius_list &
       )
       if(test_value.lt.site_value) then 
          nstuck = nstuck + 1
          if(nstuck.ge.10) then
             nattempt = nattempt + 1
             if(crude_norm.lt.site_value) &
                  crude_norm = &
                       ( crude_norm + site_value/real(nattempt) ) / 2._real32

             ! if we have tried 10 times, and still no luck, then we need to
             ! reduce the tolerance
             if(nattempt.ge.10) site_value = site_value / crude_norm
             call random_number(rtmp1)
             if(rtmp1.lt.site_value) exit walk_loop   
          end if   
       else
          nstuck = 0
          site_vector = test_vector
          site_value  = test_value
 
          if(nattempt.ge.10) test_value = test_value / crude_norm
          call random_number(rtmp1)
          if(rtmp1.lt.test_value) exit walk_loop
       end if
 
    end do walk_loop
 
    point = site_vector
    viable=.true.
   
  end function place_method_growth
!###############################################################################


!###############################################################################
  function place_method_min( &
       points, species, &
       species_index_list, viable &
  ) result(point)
    !! Global minimum placement method.
    !!
    !! This method places the atom at the gridpoint with the highest
    !! suitability.
    implicit none

    ! Arguments
    logical, intent(out) :: viable
    !! Boolean to indicate if point is viable.
    integer, intent(in) :: species
    !! Species index to add atom to.
    integer, dimension(:), intent(in) :: species_index_list
    !! List of species indices to add atoms to.
    real(real32), dimension(:,:), intent(in) :: points
    !! List of gridpoints to consider.
    real(real32), dimension(3) :: point
    !! Point to add atom to.

    ! Local variables
    integer :: species_index
    !! Index of species in list.
    integer :: best_gridpoint
    !! Index of best gridpoint.


    viable = .false.

    ! find the gridpoint with the highest suitability
    species_index = findloc(species_index_list, species, 1)
    best_gridpoint = maxloc(points(3+species_index,:), dim=1)
    if(best_gridpoint.eq.0)then
       return
    elseif(points(3+species,best_gridpoint).lt.1.E-6)then
       return
    end if

    ! return the gridpoint with the highest suitability
    point = points(1:3,best_gridpoint)
    viable = .true.
   
  end function place_method_min
!###############################################################################

end module raffle__place_methods