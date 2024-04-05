module elements
  use constants, only: real12
  implicit none

  private

  public :: element_type, element_bond_type
  public :: element_database, element_bond_database
  public :: load_elements, load_element_bonds


  type :: element_type
     character(len=3) :: name
     real(real12) :: mass
     real(real12) :: charge
     real(real12) :: energy
   contains
     procedure, pass(this) :: set
  end type element_type
  type(element_type), dimension(:), allocatable :: element_database


  type :: element_bond_type
     real(real12) :: radius_covalent
     real(real12) :: radius_vdw
     integer, dimension(2) :: coordination
     character(3), dimension(2) :: element
   end type element_bond_type
   type(element_bond_type), dimension(:), allocatable :: element_bond_database
  

contains

!!!#############################################################################
!!! set element properties from database
!!!#############################################################################
  subroutine set(this, name)
    implicit none
    class(element_type), intent(inout) :: this
    character(len=3), intent(in) :: name

    integer :: i

    do i = 1, size(element_database)
       if(element_database(i)%name .eq. name)then
          this%name = element_database(i)%name
          this%mass = element_database(i)%mass
          this%charge = element_database(i)%charge
          this%energy = element_database(i)%energy
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

    if(allocated(element_database)) deallocate(element_database)
    allocate(element_database(0))
    open(newunit=unit, file=file_, status='old', action='read')
    read(unit, *) buffer
    if(  index(trim(adjustl(buffer)),"#").ne.1 .and. &
         index(trim(adjustl(buffer)),"element").eq.0)then
       write(0,*) 'Invalid elements file'
       write(0,*) 'Expected "element" in header, found "', trim(buffer), '"'
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
       element_database = [element_database, element]
    end do
    close(unit)

  end subroutine load_elements
!!!#############################################################################


!!!#############################################################################
!!! get element bond data from file
!!!#############################################################################
  subroutine load_element_bonds(file)
    implicit none
    character(*), intent(in), optional :: file

    integer :: unit
    integer :: i, ierror
    type(element_bond_type) :: bond
    character(1024) :: buffer, file_ = "chem.in"


    if (present(file)) file_ = file
    if(allocated(element_bond_database)) deallocate(element_bond_database)
    allocate(element_bond_database(0))


    !! open file containing element bond data
    open(newunit=unit, file=file, status="old")
    read(unit, *) buffer
    if(  index(trim(adjustl(buffer)),"#").ne.1 .and. &
         index(trim(adjustl(buffer)),"element_1").eq.0)then
       write(0,*) 'Invalid elements file'
       write(0,*) 'Expected "element_1" in header, found "', trim(buffer), '"'
       stop 1
    end if


    !! read element bond data
    do 
       read(unit, '(A)', iostat=ierror) buffer
       if(is_iostat_end(ierror))then
          exit
       elseif(ierror.ne.0) then
          stop 1
       end if
       buffer = trim(adjustl(buffer))
       if(trim(buffer) .eq. "" .or. buffer(1:1) .eq. "!") cycle
       read(buffer, *) &
            bond%element(:), &
            bond%radius_covalent, bond%radius_vdw, &
            bond%coordination(:)
       check_loop: do i = 1, size(element_bond_database)
          if( ( element_bond_database(i)%element(1) .eq. bond%element(1) .and. &
                element_bond_database(i)%element(2) .eq. bond%element(2) ) .or. &
              ( element_bond_database(i)%element(1) .eq. bond%element(2) .and. &
                element_bond_database(i)%element(2) .eq. bond%element(1) ) ) then
             write(0,*) 'Error reading element bond data'
             write(0,*) 'Duplicate entry for ', bond%element(1), bond%element(2)
             stop 1
          end if
       end do check_loop
       element_bond_database = [element_bond_database, bond]
    end do
    close(unit)

  end subroutine load_element_bonds
!!!#############################################################################

end module elements