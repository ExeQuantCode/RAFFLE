module extended_geom
  use constants, only: real12, pi
  use misc_linalg, only: modu, LUinv, cross
  use rw_geom, only: basis_type, species_type
  implicit none


  private

  public :: extended_basis_type


  type, extends(basis_type) :: extended_basis_type
     real(real12) :: max_extension
     integer :: num_images
     type(species_type), dimension(:), allocatable :: image_spec
   contains
     procedure, pass(this) :: create_images
     procedure, pass(this) :: update_images
  end type extended_basis_type

contains

  subroutine create_images(this, max_bondlength)
    class(extended_basis_type), intent(inout) :: this
    real(real12), intent(in) :: max_bondlength

    type(species_type), dimension(this%nspec) :: image_species

    integer :: is, ia, i, j, k
    integer :: amax, bmax, cmax
    real(real12), dimension(3) :: vtmp1


    !---------------------------------------------------------------------------
    ! get the maximum number of lattice vectors to consider
    ! NOTE: this is not perfect
    !       won't work for extremely acute/obtuse angle cells
    !       (due to diagonal path being shorter than individual lattice vectors)
    !---------------------------------------------------------------------------
    amax = ceiling(max_bondlength/modu(this%lat(1,:)))
    bmax = ceiling(max_bondlength/modu(this%lat(2,:)))
    cmax = ceiling(max_bondlength/modu(this%lat(3,:)))


    !!! WARNING, NEED IGNORE LIST IN HERE TO ONLY APPLY TO ATOMS WE WANT TO EXTEND !!!
    !!! needs and update_images subroutine !!!
    spec_loop1: do is=1,this%nspec
       allocate( image_species(is)%atom( &
            this%spec(is)%num*(2*amax+1)*(2*bmax+1)*(2*cmax+1), &
            size(this%spec(is)%atom,2) ) &
       )
       image_species(is)%num = 0
       image_species(is)%mass = this%spec(is)%mass
       image_species(is)%charge = this%spec(is)%charge
       image_species(is)%radius = this%spec(is)%radius
       image_species(is)%name = this%spec(is)%name
       atom_loop1: do ia=1,this%spec(is)%num
          do i=-amax,amax+1,1
             vtmp1(1) = this%spec(is)%atom(ia,1) + real(i, real12)
             do j=-bmax,bmax+1,1
                vtmp1(2) = this%spec(is)%atom(ia,2) + real(j, real12)
                do k=-cmax,cmax+1,1
                   vtmp1(3) = this%spec(is)%atom(ia,3) + real(k, real12)
                   if( get_distance_from_unit_cell(vtmp1, this%lat) .le. max_bondlength ) then
                      ! add the image to the list
                      image_species(is)%num = image_species(is)%num + 1
                      image_species(is)%atom(image_species(is)%num,:3) = vtmp1
                   end if
                end do
             end do
          end do
       end do atom_loop1
    end do spec_loop1


   !  this%nspec_images = count(image_species%num.gt.0)
    allocate(this%image_spec(this%nspec))
    !  js = 0
    do is = 1, this%nspec
       ! js = js + 1
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

  end subroutine create_images


  subroutine update_images(this, max_bondlength, is, ia)
    !! Update the images for a specific atom
    implicit none

    ! Arguments
    class(extended_basis_type), intent(inout) :: this
    real(real12), intent(in) :: max_bondlength
    integer, intent(in) :: is, ia


    type(species_type) :: image_species

    integer :: i, j, k, num_images
    integer :: amax, bmax, cmax
    real(real12), dimension(3) :: vtmp1


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
    allocate( image_species%atom( &
         num_images + &
         (2*amax+1)*(2*bmax+1)*(2*cmax+1), &
         size(this%spec(is)%atom,2) ) &
    )
    image_species%atom(:num_images,:) = &
         this%image_spec(is)%atom(:num_images,:)


    !!! WARNING, NEED IGNORE LIST IN HERE TO ONLY APPLY TO ATOMS WE WANT TO EXTEND !!!
    !!! needs and update_images subroutine !!!
    do i=-amax,amax+1,1
       vtmp1(1) = this%spec(is)%atom(ia,1) + real(i, real12)
       do j=-bmax,bmax+1,1
          vtmp1(2) = this%spec(is)%atom(ia,2) + real(j, real12)
          do k=-cmax,cmax+1,1
             vtmp1(3) = this%spec(is)%atom(ia,3) + real(k, real12)
             if( get_distance_from_unit_cell(vtmp1, this%lat) .le. max_bondlength ) then
                ! add the image to the list
                num_images = num_images + 1
                image_species%atom(num_images,:3) = vtmp1
             end if
          end do
       end do
    end do
    if( num_images .eq. this%image_spec(is)%num ) return


    this%image_spec(is)%num = num_images
    deallocate(this%image_spec(is)%atom)
    allocate(this%image_spec(is)%atom( &
         num_images, &
         size(image_species%atom,2) &
    ) )
    this%image_spec(is)%atom(:,:) = &
         image_species%atom(:num_images,:)

  end subroutine update_images




  function get_distance_from_unit_cell(point, lattice, closest_point, is_cartesian) result(distance)
    implicit none
    ! Input
    real(real12), intent(in) :: point(3)        ! Point in 3D space
    real(real12), intent(in) :: lattice(3,3) ! 3x3 array representing lattice vectors as columns
    ! Output
    real(real12), intent(out), optional :: closest_point(3)  ! Closest point on the unit cell surface
    logical, optional, intent(in) :: is_cartesian ! Flag indicating whether the point is in Cartesian coordinates
    real(real12) :: distance           ! Shortest distance to the unit cell surface
    ! Local variables
    real(real12), dimension(3) :: point_
    real(real12), dimension(3,3) :: inverse_lattice
    real(real12), dimension(3) :: normal
    real(real12), dimension(3) :: plane_point
    real(real12), dimension(3) :: projection, closest_point_
    real(real12), dimension(3) :: inverse_projection
    real(real12) :: min_distance
    logical :: is_outside = .false.
    integer :: i, j, k
    integer, dimension(3) :: index_list = [1, 2, 3]
    logical :: is_cartesian_ = .false.
        


    if(present(is_cartesian)) is_cartesian_ = is_cartesian
      
    inverse_lattice = transpose( LUinv( lattice ) ) * 2._real12 * pi
    if(is_cartesian_) then
        ! Convert point to fractional coordinates
        ! point_ = matmul(LUinv(lattice), point)
        point_ = point
    else
        point_ = matmul(point, lattice)
    end if

    min_distance = huge(1._real12)

  ! get projection of point onto each face of the lattice
  ! get the length of the projection vector
  ! if negative, then the point is inside the unit cell
  ! if positive, then the point is outside the unit cell
  ! if the projection falls outside of the cell edges, use edge or corner distances
    
    do i = 1, 3
       index_list = cshift(index_list, 1)
       plane_point = 0._real12
       do j = 1, 2
          normal = (-1._real12)**j * cross(lattice(index_list(2),:), lattice(index_list(3),:))
          projection = project_point_onto_plane(point_, plane_point, normal)

          ! check if point minus projection is negative
          ! if so, it is on the wrong side of the plane and should be ignored
          if( dot_product(point_ - projection, normal) .lt. 0._real12 ) cycle
          is_outside = .true.

          ! check if projection is outside the surface
          
          inverse_projection = matmul(projection, inverse_lattice)
          if( any( inverse_projection .lt. 0._real12 ) .or. &
              any( inverse_projection .gt. 1._real12 ) ) then
             ! projection is outside the surface
             ! check if the projection is outside the edges
             ! if it is, then the closest point is the edge or corner
             ! if it is not, then the closest point is the projection
             do k = 1, 3
                if( inverse_projection(k) .lt. 0._real12 ) then
                   inverse_projection(k) = 0._real12
                else if( inverse_projection(k) .gt. 1._real12 ) then
                   inverse_projection(k) = 1._real12
                end if
             end do
             projection = matmul(inverse_projection, lattice)
          end if
          distance = norm2(point_ - projection)
          if( distance .lt. min_distance ) then
             min_distance = distance
             closest_point_ = projection
          end if

          !! makes it apply to the next iteration
          plane_point = plane_point + lattice(index_list(1),:)
       end do
    end do

    if( is_outside ) then
       distance = min_distance
    else
       distance = 0._real12
    end if

    if( present(closest_point) ) then
       if(is_cartesian_) then
          closest_point = closest_point_
       else
          closest_point = matmul(closest_point_, inverse_lattice)
       end if
    end if

  end function get_distance_from_unit_cell




  function project_point_onto_plane(point, plane_point, normal) result(output)
    implicit none
    real(real12), dimension(3), intent(in) :: point
    real(real12), dimension(3), intent(in) :: plane_point
    real(real12), dimension(3), intent(in) :: normal
    real(real12), dimension(3) :: output

    real(real12) :: distance
    real(real12), dimension(3) :: vector_to_plane
    
    vector_to_plane = point - plane_point

    distance = dot_product(vector_to_plane, normal) / dot_product(normal, normal)

    output = point - distance * normal

  end function project_point_onto_plane


end module extended_geom