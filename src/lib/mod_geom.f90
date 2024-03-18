module geom
  use constants, only: real12, pi
  use misc_linalg, only: cross, get_angle
  use vasp_file_handler, only: unitcell
  use atomtype
  implicit none


contains
  
!!!#############################################################################
!!! return the distance between two points
!!!#############################################################################
  pure function get_bondlength(A,B) result(C)
    implicit none
    real(real12), dimension(3), intent(in) :: A,B
    real(real12) :: C

    C = sqrt(dot_product(A - B, A - B))
  end function get_bondlength
!!!#############################################################################

  
!!!#############################################################################
!!! return the angle between the vector AB and the vector BC
!!!#############################################################################
  function get_bondangle(A,B,C) result(theta)
    implicit none
    real(real12), dimension(3), intent(in) :: A,B,C
    real(real12) :: theta, x
    real(real12), dimension(3) :: bond1, bond2
    integer :: i

    theta=Acos((dot_product(A-B,C-B))/(norm2(A-B)*norm2(C-B)))
    if(isnan(theta)) theta=0
  end function get_bondangle
!!!#############################################################################


!!!#############################################################################
!!! return the dihedral angle between the planes defined by the vectors ABC and BCD
!!!#############################################################################
  function get_dihedral_angle(A,B,C,D) result(theta)
    implicit none
    real(real12), dimension(3), intent(in) :: A,B,C,D
    real(real12) :: theta, x
    real(real12), dimension(3) :: bond1, bond2, AB, AC, AD
    !write(*,*) ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
    !write(*,*) A, B, C, D
    AB(:)=B(:)-A(:)
    AC(:)=C(:)-B(:)
    AD(:)=D(:)-B(:)
    !write(*,*) AB, AC, AD
    !write(*,*) "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
    
    theta=get_angle(cross(AB,AC),AD)
    !write(*,*) theta*180/3.141592
    
    !write(*,*) cross(AB,AC)!
    !Theta set to 1000 so contribution function can recognise the special case
    if(isnan(theta)) theta=1000
    !if(theta.ne.0) write(*,*) theta
  end function get_dihedral_angle
!!!#############################################################################



!!!#############################################################################
!!! return the volume of the unit cell
!!!#############################################################################
  pure function get_volume(lat) result(vol)
    implicit none
    real(real12), dimension(3,3), intent(in) :: lat
    real(real12), dimension(3) :: a,b,c, tmp
    integer :: i
    real(real12) :: vol
    a=lat(:,1)
    b=lat(:,2)
    c=lat(:,3)
    tmp=cross(b,c)
    vol = 0._real12
    do i=1,3
       vol = vol + a(i) * tmp(i)
    end do
  end function get_volume
!!!#############################################################################


!!!#############################################################################
!!! return the spherical overlap volume
!!!#############################################################################
  pure function get_sphere_overlap(r,rp,b) result(volume) 
    implicit none
    real(real12), intent(in) :: r,rp,b
    real(real12) :: volume
    real(real12) :: step1, step2, step3

    step1=pi*(r+rp-b)**2
    step2=(b**2)+(2*b*rp)-(3*(rp**2))+(2*b*r)+(6*rp*r)-3*(r**2)
    step3=1.0_real12/(12.0*b)
    volume = step1 * step2 * step3
    if(volume.lt.0._real12) volume = 0._real12
  end function get_sphere_overlap
!!!#############################################################################


!!!#############################################################################
!!! return the projection of the atoms in the unit cell
!!!#############################################################################
  pure subroutine atomprojector(position,array,unit,atomnumber,structures)
    use atomtype
    implicit none
    type (atom), dimension(:,:), intent(out) :: array
    real(real12), dimension(3), intent(in) :: position
    integer, intent(in) :: atomnumber, structures
    type(unitcell), dimension(:), intent(in) :: unit
    integer :: j,length,x,y,z,m
    
    m=0
    do x=-1,1
       do y=-1,1
          do z=-1,1
             if((x.eq.0).and.(y.eq.0).and.(z.eq.0)) cycle
             m = m + 1
             do j=1, 3
                   array(1,m)%position(j)=position(j)+&
                     &(x*unit(structures)%cell(j,1))+&
                     &(y*unit(structures)%cell(j,2))+&
                     &(z*unit(structures)%cell(j,3))
             end do
          end do
       end do
    end do

  end subroutine atomprojector
!!!#############################################################################
  
end module geom