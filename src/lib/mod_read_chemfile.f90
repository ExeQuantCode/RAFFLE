module read_chem
  use constants, only: real12
  implicit none

  private

  public :: get_element_radius

contains

  function get_element_radius(element) result(radius)
    implicit none
    character(3), dimension(:), intent(in) :: element
    real(real12), dimension(:,:,:), allocatable :: radius

    integer :: unit
    integer :: num_elements
    character(3) :: element_1, element_2
    real(real12) :: r_vdw, r_cov, c1, c2
    integer :: increment, i, j, ierror
  
    
    num_elements = size(element)
    allocate(radius(4,num_elements,num_elements), source = 0._real12)

    open(newunit=unit, file="chem.in", status="old")
    read(unit, *) buffer
    if(  index(trim(adjustl(buffer)),"#").ne.1 .and. &
         index(trim(adjustl(buffer)),"element_1").eq.0)then
       write(0,*) 'Invalid elements file'
       write(0,*) 'Expected "element_1" in header, found "', trim(buffer), '"'
       stop 1
    end if
    do 
       read(unit, *, iostat=ierror) element_1, element_2, r_cov, r_vdw, c1, c2
       if(is_iostat_end(ierror))then
          exit
       elseif(ierror.ne.0) then
          stop 1
       end if
       write(*,'(A," ",A," ",F0.3," ",F0.3)') &
            trim(adjustl(element_1)), trim(adjustl(element_2)), r_cov, r_vdw
       do i=1, num_elements 
          do j=1, num_elements
             if(element(i).eq.trim(adjustl(element_1))) then;
                if(element(j).eq.trim(adjustl(element_2))) then; 
                   radius(1,i,j) = r_cov
                   radius(2,i,j) = r_vdw
                   radius(3,i,j) = c1
                   radius(1,j,i) = r_cov
                   radius(2,j,i) = r_vdw
                   radius(3,j,i) = c1
                   if(trim(element(i)).ne.trim(element(j))) then; 
                      radius(3,i,j) = c1
                      radius(4,i,j) = c2
                      radius(4,j,i) = c1
                      radius(3,j,i) = c2 
                      write(*,*) element(i),c1,",",element(j),c2
                   end if
                   continue
                end if
             end if
             if(element(j).eq.trim(adjustl(element_1))) then;
                if(element(i).eq.trim(adjustl(element_2))) then;
                   radius(1,j,i) = r_cov
                   radius(2,j,i) = r_vdw
                   radius(3,j,i) = c1
                   if(element(i).ne.element(j)) then;
                      radius(3,j,i) = c1
                      radius(4,j,i) = c2
                      radius(4,j,i) = c1
                      radius(3,j,i) = c2 
                   end if
                   continue 
                end if
             end if
          end do 
       end do
    end do
  
    close(unit)
  end function get_element_radius

end module read_chem