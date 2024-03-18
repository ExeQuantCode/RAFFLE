module geom
  use constants, only: real12, pi
  use misc_linalg, only: cross
  implicit none


contains
  
  pure function get_bondlength(A,B) result(C)
    implicit none
    real(real12), dimension(3), intent(in) :: A,B
    real(real12) :: C

    C = sqrt(dot_product(A - B, A - B))
  end function get_bondlength

  
  function get_bondangle(A,B,C) result(theta)
    implicit none
    real(real12) :: theta, x
    real(real12), dimension(3) :: A,B,C, bond1, bond2
    integer :: i

    theta=Acos((dot_product(A-B,C-B))/(norm2(A-B)*norm2(C-B)))
    if(isnan(theta)) then 
       theta=0
    end if
  end function get_bondangle



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

  
end module geom
