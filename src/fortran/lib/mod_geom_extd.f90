module raffle__geom_extd
  !! Module to extend the basis set to include images of atoms.
  !!
  !! This module is designed to extend the basis set to include images of atoms
  !! within a specified distance of the unit cell. This is useful for
  !! calculating interactions between atoms that are not within the unit cell.
  use raffle__constants, only: real32, pi
  use raffle__misc_linalg, only: modu, cross, inverse_3x3
  use raffle__geom_rw, only: basis_type, species_type
  implicit none


  private

  public :: extended_basis_type


  type, extends(basis_type) :: extended_basis_type
     !! Extended basis set type
     real(real32) :: max_extension
     !! Maximum distance to extend the basis set
     integer :: num_images
     !! Number of images in the extended basis set
     type(species_type), dimension(:), allocatable :: image_spec
     !! Species type for the images
   contains
     procedure, pass(this) :: create_images
     !! Create the images for the basis set
     procedure, pass(this) :: update_images
     !! Update the images for a specific atom
  end type extended_basis_type

contains

!###############################################################################
  subroutine create_images(this, max_bondlength, atom_ignore_list)
    !! Create the images for the basis set.
    implicit none

    ! Arguments
    class(extended_basis_type), intent(inout) :: this
    !! Parent of the procedure. Instance of the extended basis.
    real(real32), intent(in) :: max_bondlength
    !! Maximum distance to extend the basis set.
    integer, dimension(:,:), intent(in), optional :: atom_ignore_list
    !! List of atoms to ignore when creating images.

    ! Local variables
    integer :: is, ia, i, j, k
    !! Loop indices.
    integer :: amax, bmax, cmax
    !! Maximum number of lattice vectors to consider.
    real(real32), dimension(3) :: vtmp1
    !! Temporary vector for storing atom positions.
    type(species_type), dimension(this%nspec) :: image_species
    !! Temporary store for the images.
    integer, dimension(:,:), allocatable :: atom_ignore_list_
    !! List of atoms to ignore when creating images.


    !---------------------------------------------------------------------------
    ! initialise the atom_ignore_list
    !---------------------------------------------------------------------------
    if(present(atom_ignore_list)) then
       atom_ignore_list_ = atom_ignore_list(:,:)
    else
       allocate(atom_ignore_list_(0,0))
    end if


    !---------------------------------------------------------------------------
    ! get the maximum number of lattice vectors to consider
    ! NOTE: this is not perfect
    !       won't work for extremely acute/obtuse angle cells
    !       (due to diagonal path being shorter than individual lattice vectors)
    !---------------------------------------------------------------------------
    amax = ceiling(max_bondlength/modu(this%lat(1,:)))
    bmax = ceiling(max_bondlength/modu(this%lat(2,:)))
    cmax = ceiling(max_bondlength/modu(this%lat(3,:)))


    spec_loop: do is = 1, this%nspec
       allocate( &
            image_species(is)%atom( &
                 this%spec(is)%num*(2*amax+2)*(2*bmax+2)*(2*cmax+2), &
                 size(this%spec(is)%atom,2) &
            ) &
       )
       image_species(is)%num = 0
       image_species(is)%mass = this%spec(is)%mass
       image_species(is)%charge = this%spec(is)%charge
       image_species(is)%radius = this%spec(is)%radius
       image_species(is)%name = this%spec(is)%name
       atom_loop: do ia = 1, this%spec(is)%num
          do i = 1, size(atom_ignore_list_,2), 1
             if(all(atom_ignore_list_(:,i).eq.[is,ia])) cycle atom_loop
          end do
          do i=-amax,amax+1,1
             vtmp1(1) = this%spec(is)%atom(ia,1) + real(i, real32)
             do j=-bmax,bmax+1,1
                vtmp1(2) = this%spec(is)%atom(ia,2) + real(j, real32)
                do k=-cmax,cmax+1,1
                   if( i .eq. 0 .and. j .eq. 0 .and. k .eq. 0 ) cycle
                   vtmp1(3) = this%spec(is)%atom(ia,3) + real(k, real32)
                   if( get_distance_from_unit_cell(vtmp1, this%lat) .le. &
                        max_bondlength ) then
                      ! add the image to the list
                      image_species(is)%num = image_species(is)%num + 1
                      image_species(is)%atom(image_species(is)%num,:3) = vtmp1
                   end if
                end do
             end do
          end do
       end do atom_loop
    end do spec_loop


    allocate(this%image_spec(this%nspec))
    do is = 1, this%nspec
       this%image_spec(is)%num = image_species(is)%num
       this%image_spec(is)%mass = image_species(is)%mass
       this%image_spec(is)%charge = image_species(is)%charge
       this%image_spec(is)%radius = image_species(is)%radius
       this%image_spec(is)%name = image_species(is)%name
       if(image_species(is)%num .eq. 0) cycle
       allocate(this%image_spec(is)%atom( &
            image_species(is)%num, &
            size(image_species(is)%atom,2) &
       ) )
       this%image_spec(is)%atom(:,:) = &
            image_species(is)%atom(:image_species(is)%num,:)
    end do
    this%num_images = sum( this%image_spec(:)%num )

  end subroutine create_images
!###############################################################################


!###############################################################################
  subroutine update_images(this, max_bondlength, is, ia)
    !! Update the images for a specific atom
    implicit none

    ! Arguments
    class(extended_basis_type), intent(inout) :: this
    !! Parent of the procedure. Instance of the extended basis.
    real(real32), intent(in) :: max_bondlength
    !! Maximum distance to extend the basis set.
    integer, intent(in) :: is, ia
    !! Species and atom index to update.


    ! Local variables
    integer :: i, j, k, num_images, dim
    !! Loop indices.
    integer :: amax, bmax, cmax
    !! Maximum number of lattice vectors to consider.
    type(species_type) :: image_species
    !! Temporary store for the images.
    real(real32), dimension(3) :: vtmp1
    !! Temporary vector for storing atom positions.


    !---------------------------------------------------------------------------
    ! get the maximum number of lattice vectors to consider
    ! NOTE: this is not perfect
    !       won't work for extremely acute/obtuse angle cells
    !       (due to diagonal path being shorter than individual lattice vectors)
    !---------------------------------------------------------------------------
    num_images = this%image_spec(is)%num
    amax = ceiling(max_bondlength/modu(this%lat(1,:)))
    bmax = ceiling(max_bondlength/modu(this%lat(2,:)))
    cmax = ceiling(max_bondlength/modu(this%lat(3,:)))
    dim = 3
    do i = 1, this%nspec
       if ( size(this%spec(i)%atom,2) .gt. dim) dim =  size(this%spec(i)%atom,2)
    end do
    allocate( &
         image_species%atom( &
              num_images + &
              (2*amax+2)*(2*bmax+2)*(2*cmax+2), &
              dim &
         ) &
    )
    if( num_images .ne. 0 ) then
       image_species%atom(:num_images,:3) = &
            this%image_spec(is)%atom(:num_images,:3)
    end if


    do i=-amax,amax+1,1
       vtmp1(1) = this%spec(is)%atom(ia,1) + real(i, real32)
       do j=-bmax,bmax+1,1
          vtmp1(2) = this%spec(is)%atom(ia,2) + real(j, real32)
          do k=-cmax,cmax+1,1
             if( i .eq. 0 .and. j .eq. 0 .and. k .eq. 0 ) cycle
             vtmp1(3) = this%spec(is)%atom(ia,3) + real(k, real32)
             if( get_distance_from_unit_cell(vtmp1, this%lat) .le. &
                  max_bondlength ) then
                ! add the image to the list
                num_images = num_images + 1
                image_species%atom(num_images,:3) = vtmp1
             end if
          end do
       end do
    end do
    if( num_images .eq. this%image_spec(is)%num ) return


    this%image_spec(is)%num = num_images
    if(allocated(this%image_spec(is)%atom)) deallocate(this%image_spec(is)%atom)
    allocate(this%image_spec(is)%atom( &
         num_images, &
         size(image_species%atom,2) &
    ) )
    this%image_spec(is)%atom(:,:) = &
         image_species%atom(:num_images,:)
    deallocate(image_species%atom)
    this%num_images = sum( this%image_spec(:)%num )

  end subroutine update_images
!###############################################################################


!###############################################################################
  function get_distance_from_unit_cell( &
       point, lattice, closest_point, is_cartesian&
  ) result(distance)
    !! Get the distance of a point from the unit cell.
    implicit none

    ! Arguments
    real(real32), intent(in) :: point(3)
    !! Query point.
    real(real32), intent(in) :: lattice(3,3)
    !! Lattice vectors.
    real(real32), intent(out), optional :: closest_point(3)
    !! Closest point on the unit cell surface.
    logical, optional, intent(in) :: is_cartesian
    !! Boolean whether the point is in cartesian coordinates.
    real(real32) :: distance
    !! Distance of the point from the unit cell.

    ! Local variables
    integer :: i, j, k
    !! Loop indices.
    real(real32), dimension(3) :: point_
    !! Point in cartesian coordinates.
    real(real32), dimension(3,3) :: inverse_lattice
    !! Inverse of the lattice vectors.
    real(real32), dimension(3) :: normal
    !! Normal vector to the plane.
    real(real32), dimension(3) :: plane_point
    !! Point on the plane.
    real(real32), dimension(3) :: projection, closest_point_
    !! Projection of the point onto the plane.
    real(real32), dimension(3) :: inverse_projection
    !! Inverse projection of the point onto the plane.
    real(real32) :: min_distance
    !! Minimum distance to the unit cell.
    logical :: is_outside
    !! Boolean whether the point is outside the unit cell.
    integer, dimension(3) :: index_list
    !! List of indices for the lattice vectors.
    logical :: is_cartesian_
    !! Boolean whether the point is in cartesian coordinates.
        

    !---------------------------------------------------------------------------
    ! check if the point is in cartesian coordinates
    !---------------------------------------------------------------------------
    is_outside = .false.
    is_cartesian_ = .false.
    index_list = [1, 2, 3]
    if(present(is_cartesian)) is_cartesian_ = is_cartesian
    inverse_lattice = inverse_3x3( lattice )
    if(is_cartesian_) then
       ! Convert point to fractional coordinates
       ! point_ = matmul(LUinv(lattice), point)
       point_ = point
    else
       point_ = matmul(point, lattice)
    end if

    min_distance = huge(1._real32)

    ! get projection of point onto each face of the lattice
    ! get the length of the projection vector
    ! if negative, then the point is inside the unit cell
    ! if positive, then the point is outside the unit cell
    ! if the projection falls outside of the cell edges, use edge or corner 
    ! distances
    face_loop: do i = 1, 3
       index_list = cshift(index_list, 1)
       plane_point = 0._real32
       direction_loop: do j = 1, 2
          normal = (-1._real32)**j * cross( &
               [ lattice(index_list(2),:3) ], &
               [ lattice(index_list(3),:3) ] &
          )
          normal = normal / norm2(normal)
          projection = project_point_onto_plane(point_, plane_point, normal)

          ! check if point minus projection is negative
          ! if so, it is on the wrong side of the plane and should be ignored
          if( dot_product(point_ - projection, normal) .lt. 0._real32 )then
             plane_point = plane_point + lattice(index_list(1),:)
             cycle direction_loop
          end if
          is_outside = .true.

          ! check if projection is outside the surface
          
          inverse_projection = matmul(projection, inverse_lattice)
          if( &
               any( inverse_projection .lt. 0._real32 ) .or. &
               any( inverse_projection .gt. 1._real32 ) &
          ) then
             ! projection is outside the surface
             ! check if the projection is outside the edges
             ! if it is, then the closest point is the edge or corner
             ! if it is not, then the closest point is the projection
             do k = 1, 3
                if( inverse_projection(k) .lt. 0._real32 ) then
                   inverse_projection(k) = 0._real32
                else if( inverse_projection(k) .gt. 1._real32 ) then
                   inverse_projection(k) = 1._real32
                end if
             end do
          end if
          projection = matmul(inverse_projection, lattice)
          distance = norm2(point_ - projection)
          if( distance .lt. min_distance ) then
             min_distance = distance
             closest_point_ = projection
          end if

          !! makes it apply to the next iteration
          plane_point = plane_point + lattice(index_list(1),:)
       end do direction_loop
    end do face_loop

    if( is_outside ) then
       distance = min_distance
    else
       distance = 0._real32
    end if

    if( present(closest_point) ) then
       if(is_cartesian_) then
          closest_point = closest_point_
       else
          closest_point = matmul(closest_point_, inverse_lattice)
       end if
    end if

  end function get_distance_from_unit_cell
!###############################################################################


!###############################################################################
  function project_point_onto_plane(point, plane_point, normal) result(output)
    !! Project a point onto a plane.
    implicit none

    ! Arguments
    real(real32), dimension(3), intent(in) :: point
    !! Point to project.
    real(real32), dimension(3), intent(in) :: plane_point
    !! Point on the plane.
    real(real32), dimension(3), intent(in) :: normal
    !! Normal vector to the plane.
    real(real32), dimension(3) :: output
    !! Projected point.

    ! Local variables
    real(real32) :: distance
    !! Distance of the point from the plane.
    real(real32), dimension(3) :: vector_to_plane
    !! Vector from the point to the plane.
    
    vector_to_plane = point - plane_point

    distance = &
         dot_product(vector_to_plane, normal) / dot_product(normal, normal)

    output = point - distance * normal

  end function project_point_onto_plane
!###############################################################################

end module raffle__geom_extd