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
  use rw_geom, only: bas_type
  use edit_geom, only: get_min_dist, get_min_dist_between_point_and_atom
  use buildmap, only: buildmap_POINT
  use evolver, only: gvector_container_type
  implicit none


  private

  public :: add_atom_min, add_atom_void, add_atom_walk
  public :: get_viable_gridpoints, update_viable_gridpoints


contains

!###############################################################################
  function add_atom_min(gridpoints, gvector_container, &
       basis, atom_ignore_list, &
       radius_list, viable) result(point)
    !! Global minimum placement method.
    !!
    !! This method places the atom at the gridpoint with the highest
    !! suitability.
    implicit none

    ! Arguments
    type(gvector_container_type), intent(in) :: gvector_container
    !! Distribution function (gvector) container.
    type(bas_type), intent(inout) :: basis
    !! Structure to add atom to.
    logical, intent(out) :: viable
    !! Boolean to indicate if point is viable.
    integer, dimension(:,:), intent(in) :: atom_ignore_list
    !! List of atoms to ignore (i.e. indices of atoms not yet placed).
    real(real12), dimension(:,:), intent(in) :: gridpoints
    !! List of gridpoints to consider.
    real(real12), dimension(:) :: radius_list
    !! List of radii for each pair of elements.
    real(real12), dimension(3) :: point
    !! Point to add atom to.

    ! Local variables
    integer :: i
    !! Loop indices.
    integer :: best_gridpoint
    !! Index of best gridpoint.
    real(real12), dimension(:), allocatable :: suitability_grid
    !! Suitability of each gridpoint.


    !---------------------------------------------------------------------------
    ! run buildmap_point for a set of points in the unit cell
    !---------------------------------------------------------------------------
    viable = .false.
    allocate(suitability_grid(size(gridpoints,dim=2)))
    do concurrent( i = 1:size(gridpoints,dim=2) )
       suitability_grid(i) = buildmap_POINT( gvector_container, &
            gridpoints(:,i), basis, &
            atom_ignore_list, radius_list, &
            1.1_real12, 0.95_real12)
    end do
    if(abs(maxval(suitability_grid)).lt.1.E-6) then
      deallocate(suitability_grid)
      return
    end if

    best_gridpoint = maxloc(suitability_grid, dim=1)
    deallocate(suitability_grid)

    point = gridpoints(:,best_gridpoint)
    viable = .true.
   
  end function add_atom_min
!###############################################################################


!###############################################################################
  function add_atom_void(bin_size, basis, atom_ignore_list, viable) &
       result(point)
    !! VOID placement method.
    !!
    !! This method returns the gridpoint with the lowest neighbour density.
    !! i.e. the point with the lowest density in the cell.
    implicit none

    ! Arguments
    type(bas_type), intent(inout) :: basis
    !! Structure to add atom to.
    integer, dimension(3), intent(in) :: bin_size
    !! Number of gridpoints in each direction.
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
    do i = 0, bin_size(1) - 1, 1
       do j = 0, bin_size(2) - 1, 1
          do k = 0, bin_size(3) - 1, 1
             tmpvector = [i, j, k] / real(bin_size,real12)
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
    type(bas_type), intent(inout) :: basis
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
    integer :: crude_norm
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
100 random_loop : do 
       i = i + 1
      
       !call random_number(rtmp1)
       !tmpvector = gridpoints(:,floor(rtmp1*size(gridpoints,dim=2))+1)
       do j = 1, 3
          call random_number(rtmp1)
          tmpvector(j) = rtmp1
       end do

       calculated_value = buildmap_POINT( gvector_container, &
            tmpvector, basis, &
            atom_ignore_list, radius_list, &
            1.1_real12, 0.95_real12)

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

       calculated_test = buildmap_POINT( gvector_container, &
            testvector, basis, &
            atom_ignore_list, radius_list, &
            1.1_real12, 0.95_real12)
     
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
   
             !! if we have tried 10 times, and still no luck, then we need to ...
             !! ... reduce the tolerance
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
  function get_viable_gridpoints(bin_size, basis, &
       radius_list, atom_ignore_list) result(points)
    !! Get the viable gridpoints for adding an atom.
    !!
    !! This function returns a list of gridpoints that are not too close to an
    !! existing atom.
    implicit none

    ! Arguments
    type(bas_type), intent(in) :: basis
    !! Structure to add atom to.
    integer, dimension(3), intent(in) :: bin_size
    !! Number of gridpoints in each direction.
    integer, dimension(:,:), intent(in) :: atom_ignore_list
    !! List of atoms to ignore (i.e. indices of atoms not yet placed).
    real(real12), dimension(:), intent(in) :: radius_list
    !! List of radii for each pair of elements.

    ! Local variables
    integer, dimension(:), allocatable :: pair_index
    !! List of element pair indices.
    real(real12), dimension(:,:), allocatable :: points_tmp, points
    !! List of gridpoints.

    ! Local variables
    integer :: i, j, k, l, is, ia
    !! Loop indices.
    integer :: num_points
    !! Number of gridpoints.
   

    !---------------------------------------------------------------------------
    ! get list of element pair indices
    !---------------------------------------------------------------------------
    allocate(pair_index(basis%nspec), source = 0)
    do is = 1, basis%nspec
       pair_index(is) = nint( ( basis%nspec - is/2._real12 ) * ( is - 1 ) + is )
    end do


    !---------------------------------------------------------------------------
    ! loop over all gridpoints in the unit cell and check if they are too ...
    ! ... close to an existing atom. If they are, remove them from the list ...
    ! ... of viable gridpoints
    !---------------------------------------------------------------------------
    allocate(points_tmp(3,product(bin_size)))
    num_points = 0
    grid_loop1: do i = 0, bin_size(1) - 1, 1
       grid_loop2: do j = 0, bin_size(2) - 1, 1
          grid_loop3: do k = 0, bin_size(3) - 1, 1
             do is = 1, basis%nspec
                do ia = 1, basis%spec(is)%num
                   do l = 1, size(atom_ignore_list,dim=1), 1
                      if(all(atom_ignore_list(l,:).eq.[is,ia])) cycle
                   end do
                   if( get_min_dist_between_point_and_atom( &
                        basis, &
                        [i, j, k] / real(bin_size,real12), [is,ia] ) .lt. &
                        radius_list(pair_index(is)) * 0.95_real12 ) &
                        cycle grid_loop3
                end do
             end do
             num_points = num_points + 1
             points_tmp(:,num_points) = [i, j, k] / real(bin_size,real12)
          end do grid_loop3
       end do grid_loop2
    end do grid_loop1
    allocate(points, source = points_tmp(:,:num_points))

    deallocate(points_tmp, pair_index)

  end function get_viable_gridpoints
!###############################################################################


!###############################################################################
  subroutine update_viable_gridpoints(points, basis, atom, radius)
    !! Update the viable gridpoints after a new atom has been added.
    implicit none

    ! Arguments
    type(bas_type), intent(in) :: basis
    !! Structure to add atom to.
    integer, dimension(2) :: atom
    !! Index of atom to add.
    real(real12), dimension(:,:), allocatable, intent(inout) :: points
    !! List of gridpoints.
    real(real12), intent(in) :: radius
    !! Radius of added atom.

    ! Local variables
    integer :: i
    !! Loop indices.
    integer :: num_points
    !! Number of gridpoints.
    real(real12), dimension(:,:), allocatable :: points_tmp
    !! Temporary list of gridpoints.


    !---------------------------------------------------------------------------
    ! loop over all gridpoints in the unit cell and check if they are too ...
    ! ... close to the new atom. If they are, remove them from the list of ...
    ! ... viable gridpoints
    !---------------------------------------------------------------------------
    if(.not.allocated(points)) return
    num_points = size(points,dim=2)
    i = 0
    points_tmp = points
    do while (i .le. num_points)
       i = i + 1
       if( get_min_dist_between_point_and_atom( &
             basis, points(:,i), atom ) .lt. &
             radius * 0.95_real12 ) then
          num_points = num_points - 1
          points_tmp(:,i:num_points) = points_tmp(:,i+1:num_points+1)
          i = i - 1
       end if
    end do
    if(num_points.lt.1)then
       deallocate(points)
    else
       points = points_tmp(:,:num_points)
    end if

    deallocate(points_tmp)

  end subroutine update_viable_gridpoints
!###############################################################################

end module add_atom