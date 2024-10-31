module raffle__misc_linalg
  !! Module contains various linear algebra functions and subroutines.
  use raffle__constants, only: real32, pi
  implicit none


  private

  public :: modu, cross
  public :: get_distance, get_angle, get_dihedral_angle
  public :: get_improper_dihedral_angle
  public :: inverse_3x3


  interface get_angle
     procedure get_angle_from_points, get_angle_from_vectors
  end interface get_angle

  interface get_dihedral_angle
     procedure get_dihedral_angle_from_points, get_dihedral_angle_from_vectors
  end interface get_dihedral_angle

  interface get_improper_dihedral_angle
     procedure get_improper_dihedral_angle_from_points, &
          get_improper_dihedral_angle_from_vectors
  end interface get_improper_dihedral_angle

contains

!###############################################################################
  pure function modu(vector)
    !! Return the magnitude of a vector of any size.
    implicit none

    ! Arguments
    real(real32),dimension(:), intent(in) :: vector
    !! Input vector.
    real(real32) :: modu
    !! Output magnitude.

    modu = abs( sqrt( sum(vector(:)**2) ) )
  end function modu
!###############################################################################


!###############################################################################
  pure function cross(a,b)
    !! Return the cross product of two vectors.
    implicit none

    ! Arguments
    real(real32), dimension(3), intent(in) :: a, b
    !! Input vectors.
    real(real32), dimension(3) :: cross
    !! Output cross product.

    cross(1) = a(2) * b(3) - a(3) * b(2)
    cross(2) = a(3) * b(1) - a(1) * b(3)
    cross(3) = a(1) * b(2) - a(2) * b(1)

  end function cross
!###############################################################################


!###############################################################################
  pure function get_distance(point1, point2) result(distance)
    !! Return the distance between two points.
    implicit none

    ! Arguments
    real(real32), dimension(3), intent(in) :: point1, point2
    !! Input points.
    real(real32) :: distance
    !! Output distance.

    distance = modu( point1 - point2 )

    return
  end function get_distance
!###############################################################################


!###############################################################################
  pure function get_angle_from_vectors(vector1, vector2) result(angle)
    !! Return the angle between two vectors.
    implicit none

    ! Arguments
    real(real32), dimension(3), intent(in) :: vector1, vector2
    !! Input vectors.
    real(real32) :: angle
    !! Output angle.

    angle =  dot_product(vector1,vector2) / &
         ( modu(vector1) * modu(vector2) )
    if(angle .ge. 1._real32)then
       angle = 0._real32
    elseif(angle .le. -1._real32)then
       angle = pi
    else
       angle = acos(angle)
    end if
  end function get_angle_from_vectors
!###############################################################################


!###############################################################################
  pure function get_angle_from_points(point1, point2, point3) result(angle)
    !! Return the angle formed by three points.
    !!
    !! The angle is formed by the path point1 -> point2 -> point3.
    implicit none

    ! Arguments
    real(real32), dimension(3), intent(in) :: point1, point2, point3
    !! Input points.
    real(real32) :: angle
    !! Output angle.

    angle = dot_product( point1 - point2, point3 - point2 ) / &
         ( modu( point1 - point2 ) * modu( point3 - point2 ) )
    if(angle .ge. 1._real32)then
       angle = 0._real32
    elseif(angle .le. -1._real32)then
       angle = pi
    else
       angle = acos(angle)
    end if
  end function get_angle_from_points
!###############################################################################


!###############################################################################
  pure function get_dihedral_angle_from_vectors( &
       vector1, vector2, vector3) result(angle)
    !! Return the dihedral angle between two planes.
    !!
    !! The dihedral angle is the angle between the plane defined by the cross
    !! product of two vectors and a third vector.
    !! i.e. ( vector1 x vector2 ) . vector3
    implicit none

    ! Arguments
    real(real32), dimension(3), intent(in) :: vector1, vector2, vector3
    !! Input vectors.
    real(real32) :: angle
    !! Output angle.

    angle = get_angle(cross(vector1, vector2), vector3)

  end function get_dihedral_angle_from_vectors
!###############################################################################


!###############################################################################
  pure function get_dihedral_angle_from_points( &
       point1, point2, point3, point4 &
  ) result(angle)
    !! Return the dihedral angle between two planes.
    !!
    !! The dihedral angle is the angle between the plane defined by four points.
    !! i.e. ( point2 - point1 ) x ( point3 - point2 ) . ( point4 - point2 )
    !! alt. angle between plane point1point2point3 and vector point2point4
    implicit none
    real(real32), dimension(3), intent(in) :: point1, point2, point3, point4
    real(real32) :: angle
  
    angle = get_angle(cross(point2 - point1, point3 - point2), point4 - point2)
  
  end function get_dihedral_angle_from_points
!###############################################################################


!###############################################################################
  pure function get_improper_dihedral_angle_from_vectors( &
       vector1, vector2, vector3 &
  ) result(angle)
    !! Return the improper dihedral angle between two planes.
    !!
    !! The improper dihedral angle is the angle between two planes made by
    !! three vectors.
    !! i.e. ( vector1 x vector2 ) . ( vector2 x vector3 )
    !! alt. angle between plane vector1vector2 and vector2vector3
    implicit none
    real(real32), dimension(3), intent(in) :: vector1, vector2, vector3
    real(real32) :: angle

    angle = get_angle( &
         cross(vector1, vector2), &
         cross(vector2, vector3) &
    )
    !! map angle back into the range [0, pi]
    if(angle .gt. pi) angle = 2._real32 * pi - angle


  end function get_improper_dihedral_angle_from_vectors
!###############################################################################


!###############################################################################
  pure function get_improper_dihedral_angle_from_points( &
       point1, point2, point3, point4 &
  ) result(angle)
    !! Return the improper dihedral angle between two planes.
    !!
    !! The dihedral angle is the angle between the plane defined by four points.
    !! i.e. ( point2 - point1 ) x ( point3 - point1 ) . 
    !! ( point4 - point2 ) x ( point3 - point1 )
    !! alt. angle between plane point1point2point3 and point1point3point4
    implicit none
    real(real32), dimension(3), intent(in) :: point1, point2, point3, point4
    real(real32) :: angle

    angle = get_angle( &
         cross(point2 - point1, point3 - point1), &
         cross(point3 - point1, point4 - point1) &
    )

  end function get_improper_dihedral_angle_from_points
!###############################################################################


!###############################################################################
  pure function inverse_3x3(mat) result(output)
    implicit none
    real(real32) :: det
    real(real32), dimension(3,3) :: output
    real(real32), dimension(3,3), intent(in) :: mat

    det = &
         mat(1,1) * mat(2,2) * mat(3,3) - mat(1,1) * mat(2,3) * mat(3,2) - &
         mat(1,2) * mat(2,1) * mat(3,3) + mat(1,2) * mat(2,3) * mat(3,1) + &
         mat(1,3) * mat(2,1) * mat(3,2) - mat(1,3) * mat(2,2) * mat(3,1)

    output(1,1) = +1._real32 / det * (mat(2,2) * mat(3,3) - mat(2,3) * mat(3,2))
    output(2,1) = -1._real32 / det * (mat(2,1) * mat(3,3) - mat(2,3) * mat(3,1))
    output(3,1) = +1._real32 / det * (mat(2,1) * mat(3,2) - mat(2,2) * mat(3,1))
    output(1,2) = -1._real32 / det * (mat(1,2) * mat(3,3) - mat(1,3) * mat(3,2))
    output(2,2) = +1._real32 / det * (mat(1,1) * mat(3,3) - mat(1,3) * mat(3,1))
    output(3,2) = -1._real32 / det * (mat(1,1) * mat(3,2) - mat(1,2) * mat(3,1))
    output(1,3) = +1._real32 / det * (mat(1,2) * mat(2,3) - mat(1,3) * mat(2,2))
    output(2,3) = -1._real32 / det * (mat(1,1) * mat(2,3) - mat(1,3) * mat(2,1))
    output(3,3) = +1._real32 / det * (mat(1,1) * mat(2,2) - mat(1,2) * mat(2,1))

  end function inverse_3x3
!###############################################################################

end module raffle__misc_linalg
