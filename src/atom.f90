module atomtype
  use constants, only: real12
  implicit none

  type atom
    character(3) :: name
    integer :: register, element_index
    real(real12), dimension(3) :: position
  end type atom

end module atomtype
