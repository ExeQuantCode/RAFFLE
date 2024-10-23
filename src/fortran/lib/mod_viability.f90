module raffle__viability
  !! Module to determine the viability of a set of gridpoints
  !!
  !! This module contains procedures to determine the viability of a set of
  !! points and update the viability based on new atoms being added to the cell.
  use raffle__constants, only: real32
  use raffle__misc_linalg, only: modu
  use raffle__geom_extd, only: extended_basis_type
  use raffle__dist_calcs, only: get_min_dist_between_point_and_atom
  use evaluator, only: evaluate_point
  use raffle__distribs_container, only: distribs_container_type
  implicit none


  private

  public :: get_gridpoints_and_viability, update_gridpoints_and_viability


contains

!###############################################################################
  function get_gridpoints_and_viability(distribs_container, grid, basis, &
       species_index_list, &
       radius_list, atom_ignore_list, grid_offset) result(points)
    !! Return a list of viable gridpoints and their viability for each species.
    !!
    !! This function returns the viability of all viable gridpoints.
    implicit none

    ! Arguments
    type(distribs_container_type), intent(in) :: distribs_container
    !! Distribution function (gvector) container.
    type(extended_basis_type), intent(in) :: basis
    !! Structure to add atom to.
    integer, dimension(3), intent(in) :: grid
    !! Number of gridpoints in each direction.
    real(real32), dimension(:), intent(in) :: radius_list
    !! List of radii for each pair of elements.
    integer, dimension(:), intent(in) :: species_index_list
    !! List of species indices to add atoms to.
    integer, dimension(:,:), intent(in) :: atom_ignore_list
    !! List of atoms to ignore (i.e. indices of atoms not yet placed).
    real(real32), dimension(3), intent(in) :: grid_offset
    !! Offset for gridpoints.
    real(real32), dimension(:,:), allocatable :: points
    !! List of gridpoints.

    ! Local variables
    integer :: i, j, k, l, is, ia
    !! Loop indices.
    integer :: num_points
    !! Number of gridpoints.
    real(real32) :: min_radius
    !! Minimum radius.
    real(real32), dimension(:,:), allocatable :: points_tmp
    !! Temporary list of gridpoints.


    !---------------------------------------------------------------------------
    ! loop over all gridpoints in the unit cell and check if they are too ...
    ! ... close to an existing atom. If they are, remove them from the list ...
    ! ... of viable gridpoints
    !---------------------------------------------------------------------------
    min_radius = minval(radius_list) * distribs_container%radius_distance_tol(1)
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
                                  real(grid,real32), &
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
                ] / real(grid,real32)
          end do grid_loop3
       end do grid_loop2
    end do grid_loop1
    allocate(points( 3 + basis%nspec, num_points), source = 0._real32)
    points(1:3,:) = points_tmp(1:3,1:num_points)

    deallocate(points_tmp)


    !---------------------------------------------------------------------------
    ! run evaluate_point for a set of points in the unit cell
    !---------------------------------------------------------------------------
    !do concurrent( i = 1:size(gridpoints,dim=2) )
    do i = 1, size(points,dim=2)
       do is = 1, size(species_index_list,1)
          points(3+is,i) = &
               evaluate_point( distribs_container, &
                    points(1:3,i), species_index_list(is), basis, &
                    atom_ignore_list, radius_list &
               )
       end do
    end do

  end function get_gridpoints_and_viability
!###############################################################################


!###############################################################################
  subroutine update_gridpoints_and_viability(points, distribs_container, basis, &
       species_index_list, &
       atom, radius_list, atom_ignore_list)
    !! Update the list of viable gridpoints and their viability for each 
    !! species.
    !!
    !! This subroutine updates the viability of all viable gridpoints.
    implicit none

    ! Arguments
    real(real32), dimension(:,:), allocatable, intent(inout) :: points
    !! List of gridpoints.
    type(distribs_container_type), intent(in) :: distribs_container
    !! Distribution function (gvector) container.
    type(extended_basis_type), intent(in) :: basis
    !! Structure to add atom to.
    integer, dimension(2), intent(in) :: atom
    !! Index of atom to add.
    real(real32), dimension(:), intent(in) :: radius_list
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
    real(real32) :: min_radius
    !! Minimum radius.
    real(real32) :: distance
    !! Distance between atom and gridpoint.
    logical, dimension(size(points,dim=2)) :: viable
    !! Temporary list of gridpoints.
    real(real32), dimension(:,:), allocatable :: points_tmp
    !! Temporary list of gridpoints.


    !---------------------------------------------------------------------------
    ! run evaluate_point for a set of points in the unit cell
    !---------------------------------------------------------------------------
    if(.not.allocated(points)) return
    num_points = size(points,dim=2)
    viable = .true.
    !do concurrent( i = 1:size(gridpoints,dim=2) )
    min_radius = minval(radius_list) * distribs_container%radius_distance_tol(1)
    associate( atom_pos => [ basis%spec(atom(1))%atom(atom(2),1:3) ] )
       do i = 1, num_points
          distance = modu( matmul( atom_pos - points(1:3,i), basis%lat ) )
          if( distance .lt. min_radius )then
             points(4:,i) = 0._real32
             viable(i) = .false.
             cycle
          elseif( distance .gt. distribs_container%cutoff_max(1) )then
             points(4:,i) = 0._real32
             cycle
          end if
          do is = 1, size(species_index_list,1)
             points(3+is,i) = &
                  evaluate_point( distribs_container, &
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

end module raffle__viability