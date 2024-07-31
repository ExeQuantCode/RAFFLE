!!!#############################################################################
!!! Code written by Ned Thaddeus Taylor and Francis Huw Davies
!!! Code part of the ARTEMIS group (Hepplestone research group).
!!! Think Hepplestone, think HRG.
!!!#############################################################################
!!! module contains various linear algebra functions and subroutines.
!!! module includes the following functions and subroutines:
!!! uvec             (unit vector of vector of any size)
!!! modu             (magnitude of vector of any size)
!!! cross            (cross product of two vectors)
!!!##################
!!! get_distance     (get the distance between two points)
!!! get_angle        (get the angle between two vectors)
!!! get_dihedral_angle (get the dihedral angle between two planes)
!!! get_area         (get the area made by two vectors)
!!! get_vol          (get the volume of a matrix)
!!! LUinv            (inverse of a matrix of any size using LUdecomposition)
!!! LUdecompose      (decompose a matrix into upper and lower matrices. A=LU)
!!!#############################################################################
module misc_linalg
  use constants, only: real12, pi
  implicit none


  private

  public :: uvec, modu, cross
  public :: get_distance, get_angle, get_dihedral_angle, get_area, get_vol
  public :: LUinv, LUdecompose


  interface get_angle
     procedure get_angle_from_points, get_angle_from_vectors
  end interface get_angle

  interface get_dihedral_angle
     procedure get_dihedral_angle_from_points, get_dihedral_angle_from_vectors
  end interface get_dihedral_angle


!!!updated 2021/12/09


contains
!!!#####################################################
!!! finds unit vector of an arbitrary vector
!!!#####################################################
  pure function uvec(vec)
    implicit none
    real(real12),dimension(:), intent(in)::vec
    real(real12),allocatable,dimension(:)::uvec
    allocate(uvec(size(vec)))
    uvec=vec/modu(vec)
  end function uvec
!!!#####################################################


!!!#####################################################
!!! finds modulus of an arbitrary length vector
!!!#####################################################
  pure function modu(vec)
    implicit none
    real(real12),dimension(:), intent(in)::vec
    real(real12)::modu
    modu=abs(sqrt(sum(vec(:)**2)))
  end function modu
!!!#####################################################


!!!#####################################################
!!! cross product
!!!#####################################################
  pure function cross(a,b)
    implicit none
    real(real12), dimension(3) :: cross
    real(real12), dimension(3), intent(in) :: a,b

    cross(1) = a(2)*b(3) - a(3)*b(2)
    cross(2) = a(3)*b(1) - a(1)*b(3)
    cross(3) = a(1)*b(2) - a(2)*b(1)

    return
  end function cross
!!!#####################################################



!!!#############################################################################
!!!#############################################################################
!!!  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
!!!#############################################################################
!!!#############################################################################



!!!#####################################################
!!! returns distance between two points
!!!#####################################################
  pure function get_distance(point1,point2) result(distance)
    implicit none
    real(real12) :: distance
    real(real12), dimension(3), intent(in) :: point1,point2

    distance = modu(point1-point2)

    return
  end function get_distance
!!!#####################################################


!!!#####################################################
!!! returns angle between two vectors
!!!#####################################################
  pure function get_angle_from_vectors(vec1,vec2) result(angle)
    implicit none
    real(real12), dimension(3), intent(in) :: vec1,vec2
    real(real12) :: angle

    angle = acos( dot_product(vec1,vec2)/&
         ( modu(vec1) * modu(vec2) ))
    if (isnan(angle)) angle = 0._real12

    return
  end function get_angle_from_vectors
!!!-----------------------------------------------------
!!! get the angle between vectors point1point2 and point2point3
!!! i.e. follow the path of point1 -> point2 -> point3
!!!-----------------------------------------------------
  pure function get_angle_from_points(point1, point2, point3) result(angle)
    implicit none
    real(real12), dimension(3), intent(in) :: point1, point2, point3
    real(real12) :: angle

    angle = acos( ( dot_product( point2 - point1, point3 - point2 ) ) / &
         ( modu( point2 - point1 ) * modu( point3 - point2 ) ) )
    if(isnan(angle)) angle = 0._real12
  end function get_angle_from_points
!!!#####################################################


!!!#####################################################
!!! returns the dihedral angle between the plane defined by the vectors ...
!!! vec1 x vec2 and the vector vec3
!!!#####################################################
  pure function get_dihedral_angle_from_vectors(vec1,vec2,vec3) result(angle)
    implicit none
    real(real12), dimension(3), intent(in) :: vec1,vec2,vec3
    real(real12) :: angle

    angle = get_angle(cross(vec1, vec2), vec3)

  end function get_dihedral_angle_from_vectors
!!!-----------------------------------------------------
!!! get the angle between the plane defined by ...
!!! ... point1point2point3 and the vector point2point4
!!!-----------------------------------------------------
  pure function get_dihedral_angle_from_points(point1, point2, point3, point4) &
         result(angle)
     implicit none
     real(real12), dimension(3), intent(in) :: point1, point2, point3, point4
     real(real12) :: angle
  
     angle = get_angle(cross(point2 - point1, point3 - point2), point4 - point2)
  
  end function get_dihedral_angle_from_points
!!!#####################################################


!!!#####################################################
!!! returns area made by two vectors
!!!#####################################################
  pure function get_area(a,b) result(area)
    implicit none
    real(real12), dimension(3), intent(in) :: a,b
    real(real12) :: area
    real(real12), dimension(3) :: vec

    vec = cross(a,b)
    area = sqrt(dot_product(vec,vec))

    return
  end function get_area
!!!#####################################################


!!!#####################################################
!!! returns volume of a lattice
!!!#####################################################
  function get_vol(lat) result(vol)
    implicit none
    integer :: n,i,j,k,l
    real(real12) :: vol,scale
    real(real12), dimension(3,3) :: lat
    real(real12), dimension(3) :: a,b,c

    a=lat(1,:)
    b=lat(2,:)
    c=lat(3,:)
    vol = 0._real12;scale = 1._real12
    i=1;j=2;k=3
1   do n=1,3
       vol = vol+scale*a(i)*b(j)*c(k)
       l=i;i=j;j=k;k=l
    end do
    i=2;j=1;k=3;scale=-scale
    if(scale<0._real12) goto 1

    return
  end function get_vol
!!!#####################################################


!!!#####################################################
!!! inverse of n x n matrix
!!!#####################################################
!!! doesn't work if a diagonal element = 0
!!! L = lower
!!! U = upper
!!! inmat = input nxn matrix
!!! LUinv = output nxn inverse of matrix
!!! Lz=b
!!! Ux=z
!!! x=column vectors of the inverse matrix
  function LUinv(inmat)
    implicit none
    integer :: i,m,N
    real(real12), dimension(:,:) :: inmat
    real(real12), dimension(size(inmat,1),size(inmat,1)) :: LUinv
    real(real12), dimension(size(inmat,1),size(inmat,1)) :: L,U
    real(real12), dimension(size(inmat,1)) :: c,z,x

    L=0._real12
    U=0._real12
    N=size(inmat,1)
    call LUdecompose(inmat,L,U)

!!! Lz=c
!!! c are column vectors of the identity matrix
!!! uses forward substitution to solve
    do m=1,N
       c=0._real12
       c(m)=1._real12

       z(1)=c(1)
       do i=2,N
          z(i)=c(i)-dot_product(L(i,1:i-1),z(1:i-1))
       end do


!!! Ux=z
!!! x are the rows of the inversion matrix
!!! uses backwards substitution to solve
       x(N)=z(N)/U(N,N)
       do i=N-1,1,-1
          x(i)=z(i)-dot_product(U(i,i+1:N),x(i+1:N))
          x(i)= x(i)/U(i,i)
       end do

       LUinv(:,m)=x(:)
    end do

    return
  end function LUinv
!!!#####################################################


!!!#####################################################
!!! A=LU matrix decomposer
!!!#####################################################
!!! Method: Based on Doolittle LU factorization for Ax=b
!!! doesn't work if a diagonal element = 0
!!! L = lower
!!! U = upper
!!! inmat = input nxn matrix
  subroutine LUdecompose(inmat,L,U)
    implicit none
    integer :: i,j,N
    real(real12), dimension(:,:) :: inmat,L,U
    real(real12), dimension(size(inmat,1),size(inmat,1)) :: mat

    N=size(inmat,1)
    mat=inmat
    L=0._real12
    U=0._real12

    do j=1,N
       L(j,j)=1._real12
    end do
!!! Solves the lower matrix
    do j=1,N-1
       do i=j+1,N
          L(i,j)=mat(i,j)/mat(j,j)
          mat(i,j+1:N)=mat(i,j+1:N)-L(i,j)*mat(j,j+1:N)
       end do
    end do

!!! Equates upper half of remaining mat to upper matrix
    do j=1,N
       do i=1,j
          U(i,j)=mat(i,j)
       end do
    end do

    return
  end subroutine LUdecompose
!!!#####################################################

end module misc_linalg
