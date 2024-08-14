module misc_linalg
  !! Module contains various linear algebra functions and subroutines.
  use constants, only: real12, pi
  implicit none


  private

  public :: uvec, modu, cross
  public :: get_distance, get_angle, get_dihedral_angle, get_area, get_vol
  public :: get_improper_dihedral_angle
  public :: LUinv, LUdecompose


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
pure function uvec(vector)
    !! Return the unit vector of a vector of any size.
    implicit none

    ! Arguments
    real(real12),dimension(:), intent(in)::vector
    !! Input vector.
    real(real12),allocatable,dimension(:)::uvec
    !! Output unit vector.

    allocate(uvec(size(vector)))
    uvec = vector/modu(vector)
  end function uvec
!###############################################################################


!###############################################################################
  pure function modu(vector)
    !! Return the magnitude of a vector of any size.
    implicit none

    ! Arguments
    real(real12),dimension(:), intent(in)::vector
    !! Input vector.
    real(real12)::modu
    !! Output magnitude.

    modu = abs(sqrt(sum(vector(:)**2)))
  end function modu
!###############################################################################


!###############################################################################
  pure function cross(a,b)
    !! Return the cross product of two vectors.
    implicit none

    ! Arguments
    real(real12), dimension(3), intent(in) :: a,b
    !! Input vectors.
    real(real12), dimension(3) :: cross
    !! Output cross product.

    cross(1) = a(2)*b(3) - a(3)*b(2)
    cross(2) = a(3)*b(1) - a(1)*b(3)
    cross(3) = a(1)*b(2) - a(2)*b(1)

  end function cross
!###############################################################################


!###############################################################################
  pure function get_distance(point1, point2) result(distance)
    !! Return the distance between two points.
    implicit none

    ! Arguments
    real(real12), dimension(3), intent(in) :: point1,point2
    !! Input points.
    real(real12) :: distance
    !! Output distance.

    distance = modu(point1-point2)

    return
  end function get_distance
!###############################################################################


!###############################################################################
  pure function get_angle_from_vectors(vector1, vector2) result(angle)
    !! Return the angle between two vectors.
    implicit none

    ! Arguments
    real(real12), dimension(3), intent(in) :: vector1,vector2
    !! Input vectors.
    real(real12) :: angle
    !! Output angle.

    angle = acos( dot_product(vector1,vector2)/&
         ( modu(vector1) * modu(vector2) ))
    if (isnan(angle)) angle = 0._real12

  end function get_angle_from_vectors
!###############################################################################


!###############################################################################
  pure function get_angle_from_points(point1, point2, point3) result(angle)
    !! Return the angle formed by three points.
    !!
    !! The angle is formed by the path point1 -> point2 -> point3.
    implicit none

    ! Arguments
    real(real12), dimension(3), intent(in) :: point1, point2, point3
    !! Input points.
    real(real12) :: angle
    !! Output angle.

    angle = acos( ( dot_product( point2 - point1, point3 - point2 ) ) / &
         ( modu( point2 - point1 ) * modu( point3 - point2 ) ) )
    if(isnan(angle)) angle = 0._real12
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
    real(real12), dimension(3), intent(in) :: vector1,vector2,vector3
    !! Input vectors.
    real(real12) :: angle
    !! Output angle.

    angle = get_angle(cross(vector1, vector2), vector3)

  end function get_dihedral_angle_from_vectors
!###############################################################################


!###############################################################################
  pure function get_dihedral_angle_from_points(point1, point2, point3, point4) &
         result(angle)
    !! Return the dihedral angle between two planes.
    !!
    !! The dihedral angle is the angle between the plane defined by four points.
    !! i.e. ( point2 - point1 ) x ( point3 - point2 ) . ( point4 - point2 )
    !! alt. angle between plane point1point2point3 and vector point2point4
    implicit none
    real(real12), dimension(3), intent(in) :: point1, point2, point3, point4
    real(real12) :: angle
  
    angle = get_angle(cross(point2 - point1, point3 - point2), point4 - point2)
  
  end function get_dihedral_angle_from_points
!###############################################################################


!###############################################################################
  pure function get_improper_dihedral_angle_from_vectors( &
       vector1, vector2, vector3 ) &
       result(angle)
    !! Return the improper dihedral angle between two planes.
    !!
    !! The improper dihedral angle is the angle between two planes made by
    !! three vectors.
    !! i.e. ( vector1 x vector2 ) . ( vector2 x vector3 )
       !! alt. angle between plane vector1vector2 and vector2vector3
       implicit none
    real(real12), dimension(3), intent(in) :: vector1, vector2, vector3
    real(real12) :: angle

    angle = get_angle( &
         cross(vector1, vector2), &
          cross(vector2, vector3) &
    )
    !! map angle back into the range [0, pi]
    if(angle .gt. pi) angle = 2*pi - angle


  end function get_improper_dihedral_angle_from_vectors
!###############################################################################


!###############################################################################
  pure function get_improper_dihedral_angle_from_points( &
       point1, point2, point3, point4 ) &
       result(angle)
    !! Return the improper dihedral angle between two planes.
    !!
    !! The dihedral angle is the angle between the plane defined by four points.
    !! i.e. ( point2 - point1 ) x ( point3 - point1 ) . 
    !! ( point4 - point2 ) x ( point3 - point1 )
    !! alt. angle between plane point1point2point3 and point1point3point4
    implicit none
    real(real12), dimension(3), intent(in) :: point1, point2, point3, point4
    real(real12) :: angle

    angle = get_angle( &
         cross(point2 - point1, point3 - point1), &
          cross(point3 - point1, point4 - point1) &
    )

  end function get_improper_dihedral_angle_from_points
!###############################################################################


!###############################################################################
  pure function get_area(a,b) result(area)
    !! Return the area made by two vectors.
    implicit none

    ! Arguments
    real(real12), dimension(3), intent(in) :: a,b
    !! Input vectors.
    real(real12) :: area
    !! Output area.
    real(real12), dimension(3) :: vec
    !! Cross product of a and b.

    vec = cross(a,b)
    area = sqrt(dot_product(vec,vec))

  end function get_area
!###############################################################################


!###############################################################################
  function get_vol(matrix) result(vol)
    !! Return the volume of a matrix.
    implicit none

    ! Arguments
    real(real12), dimension(3,3), intent(in) :: matrix
    !! Input matrix.

    ! Local variables
    integer :: n,i,j,k,l
    !! Loop indices.
    real(real12) :: vol,scale
    !! Volume and scale factor.
    real(real12), dimension(3) :: a,b,c
    !! Vectors of the matrix.


    a=matrix(1,:)
    b=matrix(2,:)
    c=matrix(3,:)
    vol = 0._real12;scale = 1._real12
    i=1;j=2;k=3
1   do n=1,3
       vol = vol+scale*a(i)*b(j)*c(k)
       l=i;i=j;j=k;k=l
    end do
    i=2;j=1;k=3;scale=-scale
    if(scale<0._real12) goto 1

  end function get_vol
!###############################################################################


!###############################################################################
  function LUinv(matrix)
    !! Inverse of a matrix of any size using LU decomposition.
    !!
    !! The function uses LU decomposition to solve the equation Ax=b for x.
    !! The function returns the inverse of the input matrix.
    !! L = lower matrix, U = upper matrix.
    !! NOTE: The function does not work if a diagonal element = 0.
    implicit none

    ! Arguments
    real(real12), dimension(:,:), intent(in) :: matrix
    !! Input matrix.

    ! Local variables
    integer :: i,m,N
    !! Loop indices.
    real(real12), dimension(size(matrix,1),size(matrix,1)) :: LUinv
    !! Inverse of the input matrix.
    real(real12), dimension(size(matrix,1),size(matrix,1)) :: L,U
    !! Lower and upper matrices.
    real(real12), dimension(size(matrix,1)) :: c,z,x
    !! Column vectors of the identity matrix.


    L=0._real12
    U=0._real12
    N=size(matrix,1)
    call LUdecompose(matrix,L,U)

    !---------------------------------------------------------------------------
    ! Lz=c
    ! c are column vectors of the identity matrix
    ! use forward substitution to solve
    !---------------------------------------------------------------------------
    do m=1,N
       c=0._real12
       c(m)=1._real12

       z(1)=c(1)
       do i=2,N
          z(i)=c(i)-dot_product(L(i,1:i-1),z(1:i-1))
       end do


       !------------------------------------------------------------------------
       ! Ux=z
       ! x are the rows of the inversion matrix
       ! use backwards substitution to solve
       !------------------------------------------------------------------------
       x(N)=z(N)/U(N,N)
       do i=N-1,1,-1
          x(i)=z(i)-dot_product(U(i,i+1:N),x(i+1:N))
          x(i)= x(i)/U(i,i)
       end do

       LUinv(:,m)=x(:)
    end do

  end function LUinv
!###############################################################################


!###############################################################################
  subroutine LUdecompose(matrix,L,U)
    !! Decompose a matrix into upper and lower matrices. A=LU
    !!
    !! The subroutine uses LU decomposition to solve the equation Ax=b for x.
    !! NOTE: The subroutine does not work if a diagonal element = 0.
    implicit none

    ! Arguments
    integer :: i,j,N
    !! Loop indices.
    real(real12), dimension(:,:), intent(in) :: matrix
    !! Input matrix, lower and upper matrices.
    real(real12), dimension(size(matrix,1),size(matrix,1)) :: L,U
    !! Lower and upper matrices.

    ! Local variables
    real(real12), dimension(size(matrix,1),size(matrix,1)) :: mat
    !! Temporary matrix.


    N=size(matrix,1)
    mat=matrix
    L=0._real12
    U=0._real12

    do j=1,N
       L(j,j)=1._real12
    end do


    !---------------------------------------------------------------------------
    ! Solve the lower matrix
    !---------------------------------------------------------------------------
    do j=1,N-1
       do i=j+1,N
          L(i,j)=mat(i,j)/mat(j,j)
          mat(i,j+1:N)=mat(i,j+1:N)-L(i,j)*mat(j,j+1:N)
       end do
    end do


    !---------------------------------------------------------------------------
    ! Equate upper half of remaining mat to upper matrix
    !---------------------------------------------------------------------------
    do j=1,N
       do i=1,j
          U(i,j)=mat(i,j)
       end do
    end do

  end subroutine LUdecompose
!###############################################################################

end module misc_linalg
