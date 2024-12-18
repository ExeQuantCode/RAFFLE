module raffle__element_utils
  !! Module for storing and handling element and bond data.
  !!
  !! This module contains the element and bond types, and the element and bond
  !! databases. The element and bond databases are used to store the properties
  !! of the elements and bonds in the system, respectively.
  !! The element and bond types are used by other modules to store the
  !! properties relevant to an individual system.
  use raffle__constants, only: real32
  use raffle__io_utils, only: print_warning
  implicit none

  private

  public :: element_type, element_bond_type
  public :: element_database, element_bond_database


  type :: element_type
     !! Type for storing the properties of an element.
     character(len=3) :: name
     real(real32) :: mass = 0._real32
     real(real32) :: charge = 0._real32
     real(real32) :: radius = 0._real32
     real(real32) :: energy = 0._real32
   contains
     procedure, pass(this) :: set => set_element
  end type element_type
  type(element_type), dimension(:), allocatable :: element_database


  type :: element_bond_type
     !! Type for storing the properties of a bond between two elements.    
     real(real32) :: radius_covalent = 0._real32
     character(3), dimension(2) :: element
   contains
     procedure, pass(this) :: set => set_bond
  end type element_bond_type
  type(element_bond_type), dimension(:), allocatable :: element_bond_database
  
   
  interface element_type
     !! Constructor for the element type.
     module function init_element_type( &
          name, mass, charge, energy) result(element)
       character(len=3), intent(in) :: name
       real(real32), intent(in), optional :: mass, charge, energy
       type(element_type) :: element
     end function init_element_type
  end interface element_type

   
  interface element_bond_type
     !! Constructor for the element bond type.
     module function init_element_bond_type( &
          elements, radius) result(bond)
       character(len=3), dimension(2), intent(in) :: elements
       real(real32), intent(in), optional :: radius
       type(element_bond_type) :: bond
     end function init_element_bond_type
  end interface element_bond_type


contains

!###############################################################################
  module function init_element_type(name, mass, charge, energy) result(element)
    !! Initialise an instance of the element_type.
    !!
    !! This function initialises an instance of the element_type with the
    !! provided properties.
    implicit none

    ! Arguments
    character(len=3), intent(in) :: name
    !! Element name.
    real(real32), intent(in), optional :: mass, charge, energy
    !! Element mass, charge, and energy.

    type(element_type) :: element
    !! Instance of element_type.

    element%name = name
    if(present(mass)) element%mass= mass
    if(present(charge)) element%charge = charge
    if(present(energy)) element%energy = energy

  end function init_element_type
!###############################################################################


!###############################################################################
  module function init_element_bond_type(elements, radius) result(bond)
    !! Initialise an instance of the element_bond_type.
    !!
    !! This function initialises an instance of the element_bond_type with the
    !! provided properties.
    implicit none

    ! Arguments
    character(len=3), dimension(2), intent(in) :: elements
    !! Element names.
    real(real32), intent(in), optional :: radius
    !! Element radius.

    type(element_bond_type) :: bond
    !! Instance of bond_type.

    bond%element = elements
    if(present(radius)) bond%radius_covalent = radius

  end function init_element_bond_type
!###############################################################################


!###############################################################################
  subroutine set_element(this, name, in_database)
    !! Set the element properties.
    !!
    !! This subroutine sets the properties of an element instance with data from
    !! the element database.
    !! Element properties include the mass, charge, radius, and reference energy
    !! of the element.
    implicit none

    ! Arguments
    class(element_type), intent(inout) :: this
    !! Parent. Instance of element_type.
    character(len=3), intent(in) :: name
    !! Element name.
    logical, intent(out), optional :: in_database
    !! Boolean whether pair is in database.

    ! Local variables
    integer :: i
    !! Loop index.
    character(256) :: warn_msg
    !! Warning message.


    if(present(in_database)) in_database = .false.
    if(allocated(element_database))then
       do i = 1, size(element_database)
          if(trim(element_database(i)%name) .eq. trim(name))then
             this%name = element_database(i)%name
             this%mass = element_database(i)%mass
             this%charge = element_database(i)%charge
             this%radius = element_database(i)%radius
             this%energy = element_database(i)%energy
             if(present(in_database)) in_database = .true.
             return
          end if
       end do
    end if

    write(warn_msg,'("Element ",A," not found in element database")') trim(name)
    call print_warning(warn_msg)

  end subroutine set_element
!###############################################################################


!###############################################################################
  subroutine set_bond(this, element_1, element_2, in_database)
    !! Set the bond properties for a pair of elements.
    !!
    !! This subroutine sets the properties of a bond instance with data from
    !! the bond database.
    !! Bond properties include the covalent radius of the bond.
    implicit none

    ! Arguments
    class(element_bond_type), intent(inout) :: this
    !! Parent. Instance of element_bond_type.
    character(len=3), intent(in) :: element_1, element_2
    !! Element names.
    logical, intent(out), optional :: in_database
    !! Boolean whether pair is in database.

    ! Local variables
    integer :: i
    !! Loop index.
    character(256) :: warn_msg
    !! Warning message.


    if(present(in_database)) in_database = .false.
    if(allocated(element_bond_database))then
       do i = 1, size(element_bond_database)
          if( &
               ( &
                    trim(element_bond_database(i)%element(1)) .eq. &
                    trim(element_1) .and. &
                    ( &
                         trim(element_bond_database(i)%element(2)) .eq. &
                         trim(element_2) &
                    ) &
               ) .or. ( &
                    trim(element_bond_database(i)%element(1)) .eq. &
                    trim(element_2) .and. &
                    ( &
                         trim(element_bond_database(i)%element(2)) .eq. &
                         trim(element_1) &
                    ) &
               ) &
          )then
             this%element = element_bond_database(i)%element
             this%radius_covalent = element_bond_database(i)%radius_covalent
             if(present(in_database)) in_database = .true.
             return
          end if
       end do
    end if

    write(warn_msg, &
         '("Bond between ",A," and ",A," not found in bond database")' &
    ) &
         trim(element_1), trim(element_2)
    call print_warning(warn_msg)

  end subroutine set_bond
!###############################################################################

end module raffle__element_utils