module add_atom
  !! Module to add atoms to a cell.
  !!
  !! This module contains subroutines to add atoms to a cell using different
  !! placement methods. The methods are:
  !! - min:  place the atom at the gridpoint with the highest suitability
  !! - void: place the atom in the gridpoint with the largest void space
  !! - walk: place the atom using a pseudo-random walk method
  use constants, only: real12
  use misc_linalg, only: modu
  use rw_geom, only: basis_type
  use extended_geom, only: extended_basis_type
  use edit_geom, only: get_min_dist, get_min_dist_between_point_and_atom
  use evaluator, only: evaluate_point
  use evolver, only: gvector_container_type
  implicit none


  private

  public :: add_atom_min, add_atom_void, add_atom_walk
  public :: get_gridpoints_and_viability, update_gridpoints_and_viability


contains

!###############################################################################
  function add_atom_min(points, species, &
       species_index_list, viable) result(point)
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
    real(real12), dimension(:,:), intent(in) :: points
    !! List of gridpoints to consider.
    real(real12), dimension(3) :: point
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
   
  end function add_atom_min
!###############################################################################


!###############################################################################
  function add_atom_void(grid, grid_offset, basis, atom_ignore_list, viable) &
       result(point)
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
    real(real12), dimension(3), intent(in) :: grid_offset
    !! Offset for gridpoints.
    integer, dimension(:,:), intent(in) :: atom_ignore_list
    !! List of atoms to ignore (i.e. indices of atoms not yet placed).
    logical, intent(out) :: viable
    !! Boolean to indicate if point is viable.
    real(real12), dimension(3) :: point
    !! Point to add atom to.
    
    ! Local variables
    integer :: i, j, k
    !! Loop indices.
    real(real12), dimension(3) :: best_location
    !! Index of best location to place atom.
    real(real12) :: best_location_bond, smallest_bond
    !! Bond lengths.
    real(real12), dimension(3) :: tmpvector
    !! Temporary vector for gridpoint.

   
    !---------------------------------------------------------------------------
    ! loop over all gridpoints in the unit cell and find the one with the ...
    ! ... largest void space
    !---------------------------------------------------------------------------
    viable = .false.
    best_location_bond = -huge(1._real12)
    do i = 0, grid(1) - 1, 1
       do j = 0, grid(2) - 1, 1
          do k = 0, grid(3) - 1, 1
             tmpvector = [ &
                  i + grid_offset(1), &
                  j + grid_offset(2), &
                  k + grid_offset(3) &
             ] / real(grid,real12)
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

  end function add_atom_void
!###############################################################################


!###############################################################################
  function add_atom_walk ( gvector_container, &
       basis, atom_ignore_list, &
       radius_list, viable) result(point)
    !! Pseudo-random walk placement method.
    !!
    !! This method places the atom using a pseudo-random walk method.
    !! An initial point is chosen at random, and then points nearby are tested
    !! to see if they are more suitable than the current point. If they are,
    !! the query point is moved to the new point and the process is repeated.
    !! The process is repeated, with each point being tested against a random
    !! number. If the random number is less than the suitability of the point,
    !! the atom is placed at that point.
    implicit none

    ! Arguments
    type(gvector_container_type), intent(in) :: gvector_container
    !! Distribution function (gvector) container.
    type(extended_basis_type), intent(inout) :: basis
    !! Structure to add atom to.
    logical, intent(out) :: viable
    !! Boolean to indicate if point is viable.
    integer, dimension(:,:), intent(in) :: atom_ignore_list
    !! List of atoms to ignore (i.e. indices of atoms not yet placed).
    real(real12), dimension(:), intent(in) :: radius_list
    !! List of radii for each pair of elements.
    real(real12), dimension(3) :: point
    !! Point to add atom to.

    ! Local variables
    integer :: i, j, k, l
    !! Loop indices.
    real(real12) :: rtmp1
    !! Random number.
    real(real12) :: crude_norm
    !! Crude normalisation.
    real(real12) :: calculated_value, calculated_test
    !! Viability values.
    real(real12), dimension(3) :: tmpvector, testvector
    !! Vectors for gridpoints.
   

    !---------------------------------------------------------------------------
    ! test a random point in the unit cell
    !---------------------------------------------------------------------------
    i = 0
    viable = .false.
    crude_norm = 0._real12
    random_loop : do 
       i = i + 1
      
       !call random_number(rtmp1)
       !tmpvector = gridpoints(:,floor(rtmp1*size(gridpoints,dim=2))+1)
       do j = 1, 3
          call random_number(rtmp1)
          tmpvector(j) = rtmp1
       end do

       calculated_value = evaluate_point( gvector_container, &
            tmpvector, atom_ignore_list(1,1), basis, &
            atom_ignore_list, radius_list &
       )
       call random_number(rtmp1)
       if (rtmp1.lt.calculated_value) exit random_loop
 
       !! NOTE: HARDCODED LIMIT ON NUMBER OF TRIES. SET IN INPUT FILE
       if(i.ge.10000) return
    end do random_loop


    !---------------------------------------------------------------------------
    ! now do a pseudo-random walk to find a suitable point to place the atom
    !---------------------------------------------------------------------------
    k = 0
    l = 0
    walk_loop : do
       do j=1, 3
          call random_number(rtmp1)
          if(k.gt.10) then 
             testvector(j) = tmpvector(j) + &
                  ( rtmp1 * 2._real12 - 1._real12 ) * &
                  0.01_real12
          else
             testvector(j) = tmpvector(j) + &
             ( rtmp1 * 2._real12 - 1._real12 ) * &
             0.1_real12
          end if 
       end do
       testvector = testvector - floor(testvector)

       calculated_test = evaluate_point( gvector_container, &
            testvector, atom_ignore_list(1,1), basis, &
            atom_ignore_list, radius_list &
       )     
       if(calculated_test.lt.calculated_value) then 
          l = l + 1
          if(l.ge.10) then
             call random_number(rtmp1)
             if(crude_norm.lt.calculated_value) crude_norm = calculated_value

             if (rtmp1.lt.calculated_value) exit walk_loop
             if(k.ge.10) then 
                calculated_value = calculated_value / crude_norm
                if (rtmp1.lt.calculated_value) exit walk_loop
             end if
   
             !! if we have tried 10 times, and still no luck, then we need to
             !! reduce the tolerance
             k = k + 1
          end if   
          cycle walk_loop
       end if

       l=0
       tmpvector = testvector
       calculated_value = calculated_test
       call random_number(rtmp1)
 
       if(k.gt.10) calculated_test = calculated_test / crude_norm
       if (rtmp1.lt.calculated_test) then
          tmpvector = testvector
          exit walk_loop
       end if
 
    end do walk_loop
 
    point = tmpvector
    viable=.true.
   
  end function add_atom_walk
!###############################################################################


!###############################################################################
  function get_gridpoints_and_viability(gvector_container, grid, basis, &
       species_index_list, &
       radius_list, atom_ignore_list, grid_offset) result(points)
    !! Return a list of viable gridpoints and their viability for each species.
    !!
    !! This function returns the viability of all viable gridpoints.
    implicit none

    ! Arguments
    type(gvector_container_type), intent(in) :: gvector_container
    !! Distribution function (gvector) container.
    type(extended_basis_type), intent(in) :: basis
    !! Structure to add atom to.
    integer, dimension(3), intent(in) :: grid
    !! Number of gridpoints in each direction.
    real(real12), dimension(:), intent(in) :: radius_list
    !! List of radii for each pair of elements.
    integer, dimension(:), intent(in) :: species_index_list
    !! List of species indices to add atoms to.
    integer, dimension(:,:), intent(in) :: atom_ignore_list
    !! List of atoms to ignore (i.e. indices of atoms not yet placed).
    real(real12), dimension(3), intent(in) :: grid_offset
    !! Offset for gridpoints.
    real(real12), dimension(:,:), allocatable :: points
    !! List of gridpoints.

    ! Local variables
    integer :: i, j, k, l, is, ia
    !! Loop indices.
    integer :: num_points
    !! Number of gridpoints.
    real(real12) :: min_radius
    !! Minimum radius.
    real(real12), dimension(:,:), allocatable :: points_tmp
    !! Temporary list of gridpoints.


    !---------------------------------------------------------------------------
    ! loop over all gridpoints in the unit cell and check if they are too ...
    ! ... close to an existing atom. If they are, remove them from the list ...
    ! ... of viable gridpoints
    !---------------------------------------------------------------------------
    min_radius = minval(radius_list) * gvector_container%radius_distance_tol(1)
    allocate(points_tmp(3,product(grid)))
    num_points = 0
    grid_loop1: do i = 0, grid(1) - 1, 1
       grid_loop2: do j = 0, grid(2) - 1, 1
          grid_loop3: do k = 0, grid(3) - 1, 1
             do is = 1, basis%nspec
                atom_loop: do ia = 1, basis%spec(is)%num
                   do l = 1, size(atom_ignore_list,dim=1), 1
                      if(all(atom_ignore_list(l,:).eq.[is,ia])) cycle atom_loop
                   end do
                   if( get_min_dist_between_point_and_atom( &
                             basis, &
                             [ &
                                  i + grid_offset(1), &
                                  j + grid_offset(2), &
                                  k + grid_offset(3) &
                             ] / &
                                  real(grid,real12), &
                             [is,ia] &
                        ) .lt. &
                        min_radius &
                   ) cycle grid_loop3
                end do atom_loop
             end do
             num_points = num_points + 1
             points_tmp(:,num_points) = [ &
                       i + grid_offset(1), &
                       j + grid_offset(2), &
                       k + grid_offset(3) &
                ] / real(grid,real12)
          end do grid_loop3
       end do grid_loop2
    end do grid_loop1
    allocate(points( 3 + basis%nspec, num_points), source = 0._real12)
    points(1:3,:) = points_tmp(1:3,1:num_points)

    deallocate(points_tmp)


    !---------------------------------------------------------------------------
    ! run evaluate_point for a set of points in the unit cell
    !---------------------------------------------------------------------------
    !do concurrent( i = 1:size(gridpoints,dim=2) )
    do i = 1, size(points,dim=2)
       do is = 1, size(species_index_list,1)
          points(3+is,i) = &
               evaluate_point( gvector_container, &
                    points(1:3,i), species_index_list(is), basis, &
                    atom_ignore_list, radius_list &
               )
       end do
    end do

  end function get_gridpoints_and_viability
!###############################################################################


!###############################################################################
  subroutine update_gridpoints_and_viability(points, gvector_container, basis, &
       species_index_list, &
       atom, radius_list, atom_ignore_list)
    !! Update the list of viable gridpoints and their viability for each 
    !! species.
    !!
    !! This subroutine updates the viability of all viable gridpoints.
    implicit none

    ! Arguments
    real(real12), dimension(:,:), allocatable, intent(inout) :: points
    !! List of gridpoints.
    type(gvector_container_type), intent(in) :: gvector_container
    !! Distribution function (gvector) container.
    type(extended_basis_type), intent(in) :: basis
    !! Structure to add atom to.
    integer, dimension(2), intent(in) :: atom
    !! Index of atom to add.
    real(real12), dimension(:), intent(in) :: radius_list
    !! List of radii for each pair of elements.
    integer, dimension(:), intent(in) :: species_index_list
    !! List of species indices to add atoms to.
    integer, dimension(:,:), intent(in) :: atom_ignore_list
    !! List of atoms to ignore (i.e. indices of atoms not yet placed).

    ! Local variables
    integer :: i, j, is
    !! Loop indices.
    integer :: num_points
    !! Number of gridpoints.
    real(real12) :: min_radius
    !! Minimum radius.
    real(real12) :: distance
    !! Distance between atom and gridpoint.
    logical, dimension(size(points,dim=2)) :: viable
    !! Temporary list of gridpoints.
    real(real12), dimension(:,:), allocatable :: points_tmp
    !! Temporary list of gridpoints.


    !---------------------------------------------------------------------------
    ! run evaluate_point for a set of points in the unit cell
    !---------------------------------------------------------------------------
    if(.not.allocated(points)) return
    num_points = size(points,dim=2)
    viable = .true.
    !do concurrent( i = 1:size(gridpoints,dim=2) )
    min_radius = minval(radius_list) * gvector_container%radius_distance_tol(1)
    associate( atom_pos => [ basis%spec(atom(1))%atom(atom(2),1:3) ] )
       do i = 1, num_points
          distance = modu( matmul( atom_pos - points(1:3,i), basis%lat ) )
          if( distance .lt. min_radius )then
             points(4:,i) = 0._real12
             viable(i) = .false.
             cycle
          elseif( distance .gt. gvector_container%cutoff_max(1) )then
             points(4:,i) = 0._real12
             cycle
          end if
          do is = 1, size(species_index_list,1)
             points(3+is,i) = &
                  evaluate_point( gvector_container, &
                       points(1:3,i), species_index_list(is), basis, &
                       atom_ignore_list, radius_list &
                  )
          end do
       end do
    end associate

    num_points = count(viable)
    if(num_points.lt.1)then
       deallocate(points)
       return
    end if
    allocate(points_tmp(3+basis%nspec,num_points))

    i = 0
    j = 0
    do while (i .lt. num_points)
       j = j + 1
       if(.not.viable(j)) cycle
       i = i + 1
       points_tmp(:,i) = points(:,j)
    end do
    deallocate(points)
    allocate(points, source = points_tmp)

  end subroutine update_gridpoints_and_viability
!###############################################################################

end module add_atom