module geom
  use constants, only: real12, pi
  use misc_linalg, only: cross, get_angle
  use vasp_file_handler, only: unitcell
  use atomtype
  implicit none

  private

  public :: get_bondlength, get_bondangle, get_dihedral_angle
  public :: get_volume, get_sphere_overlap
  public :: atomprojector
  public :: get_random_unit_cell


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
  pure subroutine atomprojector(position,array,lattice,atomnumber,structures)
    use atomtype
    implicit none
    type (atom), dimension(:,:), intent(out) :: array
    real(real12), dimension(3), intent(in) :: position
    integer, intent(in) :: atomnumber, structures
    real(real12), dimension(3,3), intent(in) :: lattice
    integer :: j,length,x,y,z,m
    
    m=0
    do x=-1,1
       do y=-1,1
          do z=-1,1
             if((x.eq.0).and.(y.eq.0).and.(z.eq.0)) cycle
             m = m + 1
             do j=1, 3
                   array(1,m)%position(j)=position(j)+&
                     &(x*lattice(j,1))+&
                     &(y*lattice(j,2))+&
                     &(z*lattice(j,3))
             end do
          end do
       end do
    end do

  end subroutine atomprojector
!!!#############################################################################
  

!!!#############################################################################
!!! return a random unit cell
!!!#############################################################################
  function get_random_unit_cell(bravais_type, angle, volume) result(lattice)
    implicit none
    integer, intent(in) :: bravais_type
    real(real12), intent(out) :: volume
    real(real12), dimension(3), intent(out) :: angle
    real(real12), dimension(3,3) :: lattice

    integer :: i
    real(real12) :: rtmp1


    !! Initialises lattice, which is a cubic unit serving as the basis for the random unit cel
    lattice = 0._real12
    !! volume keeps a running total of the "volume" in the loosest sense of the word  
    volume = 1._real12

    do i = 1, 3
       call random_number(rtmp1)
       lattice(i,i) = 0.75_real12 + rtmp1 * 2.25_real12
       volume = volume * rtmp1
    end do
    !! to convert from nanometres to angstroms???
    volume = volume * 1000._real12



!!!--------------------------------------------------------------!!!
!!!Sets the random angles between the unit vectors between 60-120!!!
!!!--------------------------------------------------------------!!!
    angle(:)=0
!!! BRAVAIS LATTICES 
    select case(bravais_type)
    case(1) !!! Triclinic
       do i=1, 3
          call random_number(rtmp1)
          rtmp1 = ( rtmp1 * 60._real12 ) + 60._real12  
          rtmp1 = ( rtmp1 * pi ) / ( 180._real12 )                                                         !     
          angle(i) = rtmp1 
       end do
    case(2) !!! Cubic
       volume = 1._real12
       call random_number(rtmp1) 
       lattice(1,1) = 0.75_real12 + rtmp1 * 2.25_real12
       volume = volume * ( rtmp1 ** 3._real12 ) * 1000._real12
       lattice(2,2) = lattice(1,1) 
       lattice(3,3) = lattice(1,1)
       angle(:) = pi / 2._real12
    case(3) !!! Monoclinic
       angle(3) = pi / 2._real12
       angle(1) = pi / 2._real12
       call random_number(rtmp1) 
       rtmp1 = ( rtmp1 * 60._real12 ) + 60._real12  
       rtmp1 = ( rtmp1 * pi ) / ( 180.0_real12 )
       angle(2) = rtmp1
    case(4) !!! Orthorhombic
       angle(:) = pi / 2._real12
    case(5) !!! Tetragonal B (WHAT DOES THE B MEAN??)
       volume = 1._real12 
       call random_number(rtmp1) 
       lattice(1,1) = 0.75_real12 + rtmp1 * 2.25_real12
       lattice(2,2) = lattice(1,1) 
       volume = volume * ( rtmp1 ** 2._real12 )
       call random_number(rtmp1) 
       lattice(3,3) = 0.75_real12 + rtmp1 * 2.25_real12 
       volume = volume * rtmp1 * 1000._real12
       angle(:) = pi / 2._real12
    case(6) !!! Rhombohedral very broken/ Trigonal (WHY???)
       volume = 1._real12
       call random_number(rtmp1)
       lattice(1,1) = 0.75_real12 + rtmp1 * 2.25_real12
       ! lattice(1,1)=5.0
       lattice(2,2) = lattice(1,1) 
       lattice(3,3) = lattice(1,1) 
       volume = volume * ( rtmp1 ** 3._real12 ) * 1000._real12 
       !write(*,*) lattice(1,1)
       call random_number(rtmp1)
       rtmp1 = ( rtmp1 * 60._real12 ) + 60._real12 
       rtmp1 = ( rtmp1 * pi ) / ( 180._real12 )
       angle(:) = rtmp1
       ! angle(:) = 1.75_real12 * pi / 3._real12
    case(7) !!! Hexagonal
       angle(1) = pi / 2._real12
       angle(2) = pi / 2._real12
       angle(3) = 2._real12 * pi / 3._real12
       call random_number(rtmp1) 
       lattice(1,1) = 0.75_real12 + rtmp1 * 2.25_real12
       lattice(2,2) = lattice(1,1) 
       volume = rtmp1 ** 2._real12
       call random_number(rtmp1) 
       lattice(3,3) = 0.75_real12 + rtmp1 * 2.25_real12
       volume = volume * rtmp1 * 1000._real12
    case default
       write(*,'("Bravais type ",I0," not recognised")') bravais_type
       stop 1
    end select

  end function get_random_unit_cell
!!!#############################################################################

end module geom