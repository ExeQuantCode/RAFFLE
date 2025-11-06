module raffle__bounds
  use raffle__constants, only: real32, pi
  use raffle__io_utils, only: stop_program
  use raffle__misc, only: to_upper
  implicit none

  private

  public :: abstract_bounds_type, box_bounds_type, sphere_bounds_type
  public :: bounds_container_type


  type, abstract :: abstract_bounds_type
     character(len=20) :: name
     !! Name of bound shape.
     real(real32), dimension(3) :: origin
     !! Origin of bound shape.
     logical :: is_fractional_coordinates = .false.
     !! Are the coordinates fractional?
   contains
     procedure(abstract_is_within_bounds), deferred, pass(this) :: is_within_bounds
     !! Function to check if a point is within the bounds.
     procedure(abstract_get_random_point_within_bounds), deferred, pass(this) :: &
          get_random_point_within_bounds
     !! Function to get a random point within the bounds.
     procedure(abstract_get_extent), deferred, pass(this) :: get_extent
     !! Function to get the extent of the bounds.
  end type abstract_bounds_type

  abstract interface
     module function abstract_is_within_bounds(this, point, lattice) result(output)
       class(abstract_bounds_type), intent(in) :: this
       real(real32), dimension(3), intent(in) :: point
       real(real32), dimension(3,3), intent(in) :: lattice
       logical :: output
     end function abstract_is_within_bounds

     module function abstract_get_random_point_within_bounds(this) result(point)
       class(abstract_bounds_type), intent(in) :: this
       real(real32), dimension(3) :: point
     end function abstract_get_random_point_within_bounds

     module function abstract_get_extent(this) result(extent)
       class(abstract_bounds_type), intent(in) :: this
       real(real32), dimension(2,3) :: extent
     end function abstract_get_extent
  end interface


  type, extends(abstract_bounds_type) :: box_bounds_type
     real(real32), dimension(3) :: lengths
     !! Lengths of the box in each direction.
   contains
     procedure :: is_within_bounds => box_is_within_bounds
     procedure :: get_random_point_within_bounds => box_get_random_point_within_bounds
     procedure :: get_extent => box_get_extent
  end type box_bounds_type

  interface box_bounds_type
     module procedure :: init_box_bounds
  end interface box_bounds_type

  ! type, extends(abstract_bounds_type) :: parallelepiped_bounds_type
  !    real(real32), dimension(3,3) :: vectors
  !    !! Vectors defining the parallelepiped.
  !  contains
  !    procedure :: is_within_bounds => parallelepiped_is_within_bounds
  !    procedure :: get_random_point_within_bounds => box_get_random_point_within_bounds
  !    procedure :: get_extent => box_get_extent
  ! end type parallelepiped_bounds_type

  ! interface parallelepiped_bounds_type
  !    module procedure :: init_parallelepiped_bounds
  ! end interface parallelepiped_bounds_type

  type, extends(abstract_bounds_type) :: sphere_bounds_type
     real(real32) :: radius
     !! Radius of the sphere.
   contains
     procedure :: is_within_bounds => sphere_is_within_bounds
     procedure :: &
          get_random_point_within_bounds => sphere_get_random_point_within_bounds
     procedure :: get_extent => sphere_get_extent
  end type sphere_bounds_type

  interface sphere_bounds_type
     module procedure :: init_sphere_bounds
  end interface sphere_bounds_type



  type :: bounds_container_type
     class(abstract_bounds_type), allocatable :: bounds
     !! Container for different types of bounds.
   contains
     procedure :: add_bounds
  end type bounds_container_type

contains

  function init_box_bounds(lengths, origin, is_fractional_coordinates) result(box)
    implicit none
    real(real32), dimension(:), intent(in) :: lengths
    real(real32), dimension(3), intent(in) :: origin
    logical, intent(in), optional :: is_fractional_coordinates
    type(box_bounds_type) :: box

    if(size(lengths) .ne. 3 .and. size(lengths) .ne. 1) then
       call stop_program('Box lengths must be a 1- or 3-element array.')
       return
    elseif(any(lengths .le. 0.0_real32)) then
       call stop_program('Box lengths must be positive.')
       return
    end if

    box%name = 'box'
    box%origin = origin
    box%lengths = lengths
    if(present(is_fractional_coordinates)) then
       box%is_fractional_coordinates = is_fractional_coordinates
    end if
  end function init_box_bounds

  ! function init_parallelepiped_bounds(vectors, origin, is_fractional_coordinates) result(parallelepiped)
  !   implicit none
  !   real(real32), dimension(3,3), intent(in) :: vectors
  !   real(real32), dimension(3), intent(in) :: origin
  !   logical, intent(in), optional :: is_fractional_coordinates
  !   type(parallelepiped_bounds_type) :: parallelepiped

  !   ! Check that the vectors are linearly independent
  !   if (abs(det(vectors)) .lt. 1.e-6_real32) then
  !      call stop_program('Parallelepiped vectors must be linearly independent.')
  !      return
  !   end if

  !   parallelepiped%name = 'parallelepiped'
  !   parallelepiped%origin = origin
  !   parallelepiped%vectors = vectors
  !   if(present(is_fractional_coordinates)) then
  !      parallelepiped%is_fractional_coordinates = is_fractional_coordinates
  !   end if
  ! end function init_parallelepiped_bounds

  function init_sphere_bounds(radius, origin, is_fractional_coordinates) result(sphere)
    implicit none
    real(real32), intent(in) :: radius
    real(real32), dimension(3), intent(in) :: origin
    logical, intent(in), optional :: is_fractional_coordinates
    type(sphere_bounds_type) :: sphere

    if(radius .le. 0.0_real32) then
       call stop_program('Sphere radius must be positive.')
       return
    end if

    sphere%name = 'sphere'
    sphere%origin = origin
    sphere%radius = radius
    if(present(is_fractional_coordinates)) then
       sphere%is_fractional_coordinates = is_fractional_coordinates
    end if
  end function init_sphere_bounds

  function box_is_within_bounds(this, point, lattice) result(output)
    implicit none
    class(box_bounds_type), intent(in) :: this
    real(real32), dimension(3), intent(in) :: point
    real(real32), dimension(3,3), intent(in) :: lattice
    logical :: output
    real(real32), dimension(3) :: diff
    diff = point - this%origin
    output = all(diff >= 0.0_real32) .and. all(diff <= this%lengths)
  end function box_is_within_bounds

  function sphere_is_within_bounds(this, point, lattice) result(output)
    implicit none
    class(sphere_bounds_type), intent(in) :: this
    real(real32), dimension(3), intent(in) :: point
    real(real32), dimension(3,3), intent(in) :: lattice
    logical :: output
    real(real32), dimension(3) :: diff
    diff = point - this%origin
    ! map into shortest vector in periodic box if necessary
    diff = diff - ceiling(diff - 0.5_real32)
    output = norm2(matmul(diff, lattice)) <= this%radius**2
  end function sphere_is_within_bounds

  function box_get_random_point_within_bounds(this) result(point)
    implicit none
    class(box_bounds_type), intent(in) :: this
    real(real32), dimension(3) :: point
    integer :: i
    do i = 1, 3
       call random_number(point(i))
       point(i) = &
            this%origin(i) - this%lengths(i)/2.0_real32 + point(i)*this%lengths(i)
    end do
  end function box_get_random_point_within_bounds

  function sphere_get_random_point_within_bounds(this) result(point)
    implicit none
    class(sphere_bounds_type), intent(in) :: this
    real(real32), dimension(3) :: point
    real(real32) :: u, v, theta, phi, r
    call random_number(r)
    call random_number(u)
    call random_number(v)
    theta = acos(1.0_real32 - 2.0_real32 * u)
    phi = 2.0_real32 * pi * v
    r = this%radius * (r**(1.0_real32/3.0_real32))
    point(1) = this%origin(1) + r * sin(theta) * cos(phi)
    point(2) = this%origin(2) + r * sin(theta) * sin(phi)
    point(3) = this%origin(3) + r * cos(theta)
  end function sphere_get_random_point_within_bounds

  function box_get_extent(this) result(extent)
    implicit none
    class(box_bounds_type), intent(in) :: this
    real(real32), dimension(2,3) :: extent
    extent(1,:) = this%origin
    extent(2,:) = this%origin + this%lengths
  end function box_get_extent

  function sphere_get_extent(this) result(extent)
    implicit none
    class(sphere_bounds_type), intent(in) :: this
    real(real32), dimension(2,3) :: extent
    extent(1,:) = this%origin - this%radius
    extent(2,:) = this%origin + this%radius
  end function sphere_get_extent

  subroutine add_bounds( &
       this, shape, origin, lengths, vectors, is_fractional_coordinates, exit_code &
  )
    class(bounds_container_type), intent(inout) :: this
    character(len=*), intent(in) :: shape
    real(real32), dimension(3), intent(in) :: origin
    real(real32), dimension(:), intent(in), optional :: lengths
    real(real32), dimension(3,3), intent(in), optional :: vectors
    logical, intent(in), optional :: is_fractional_coordinates
    integer, intent(out), optional :: exit_code

    integer :: exit_code_

    ! Local variables
    character(len=20) :: shape_u
    !! Uppercase version of shape.

    exit_code_ = 0
    if (allocated(this%bounds)) deallocate(this%bounds)

    shape_u = to_upper(trim(shape))
    select case(shape_u)
    case('BOX','CUBOID','CUBE')
       if(.not.present(lengths))then
          call stop_program("Lengths must be provided for box bounds")
          exit_code_ = 1
       else
          this%bounds = box_bounds_type(lengths=lengths, origin=origin)
       end if
    case('SPHERE')
       if(.not.present(lengths))then
          call stop_program("Lengths must be provided for sphere bounds")
          exit_code_ = 1
       elseif(lengths(1).le.0.0_real32)then
          call stop_program("Radius must be positive for sphere bounds")
          exit_code_ = 1
       else
          this%bounds = sphere_bounds_type(radius=lengths(1), origin=origin)
       end if
    case('PARALLELEPIPED')
       if(.not.present(vectors))then
          call stop_program("Vectors must be provided for parallelepiped bounds")
          exit_code_ = 1
       else
          write(*,*) "Parallelepiped bounds are not currently supported."
          ! this%bounds = parallelepiped_bounds_type(vectors=vectors, origin=origin)
       end if
    case default
       call stop_program("Invalid shape: "//trim(shape))
       exit_code_ = 1
    end select

    if(present(exit_code)) exit_code = exit_code_
    if(exit_code_ .ne. 0)then
       return
    end if
    if(present(is_fractional_coordinates)) &
         this%bounds%is_fractional_coordinates = is_fractional_coordinates

  end subroutine add_bounds

end module raffle__bounds
