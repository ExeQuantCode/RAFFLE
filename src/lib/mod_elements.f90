module elements
  use constants, only: real12
  implicit none

  private

  public :: element_type, elements_database, load_elements


  type :: element_type
     character(len=3) :: name
     real(real12) :: mass
     real(real12) :: charge
     real(real12) :: energy
   contains
     procedure, pass(this) :: set
  end type element_type
  type(element_type), dimension(:), allocatable :: elements_database


contains

!!!#############################################################################
!!! set element properties from database
!!!#############################################################################
  subroutine set(this, name)
    implicit none
    class(element_type), intent(inout) :: this
    character(len=3), intent(in) :: name

    integer :: i

    do i = 1, size(elements_database)
       if(elements_database(i)%name .eq. name)then
          this%name = elements_database(i)%name
          this%mass = elements_database(i)%mass
          this%charge = elements_database(i)%charge
          this%energy = elements_database(i)%energy
          return
       end if
    end do

    write(0,*) 'Element ', name, ' not found in database'
    stop 1

  end subroutine set
!!!#############################################################################


!!!#############################################################################
!!! load elements data from file
!!!#############################################################################
  subroutine load_elements(file)
    implicit none
    character(*), intent(in), optional :: file

    integer :: i, unit, ierror
    type(element_type) :: element
    character(1024) :: buffer, file_ = "elements.dat"


    if (present(file)) file_ = file

    if(allocated(elements_database)) deallocate(elements_database)
    allocate(elements_database(0))
    open(newunit=unit, file=file_, status='old', action='read')
    read(unit, *) buffer
    if(index(trim(adjustl(buffer)),"element").ne.1)then
       write(0,*) 'Invalid elements file'
       write(0,*) 'Expected "elements" in header, found "', trim(buffer), '"'
       stop 1
    end if
    do
       read(unit, '(A)', iostat=ierror) buffer
       if(is_iostat_end(ierror))then
          exit
       elseif(ierror .ne. 0)then
          write(0,*) 'Error reading elements file'
          stop 1
       end if
       buffer = trim(adjustl(buffer))
       if(trim(buffer) .eq. "" .or. buffer(1:1) .eq. "!") cycle
       read(buffer, *) element%name, element%energy, element%mass, element%charge
       elements_database = [elements_database, element]
    end do
    close(unit)

  end subroutine load_elements
!!!#############################################################################

end module elements