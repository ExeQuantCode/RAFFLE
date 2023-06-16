module atomtype 
implicit none 
 TYPE atom
 character(3) :: name
 integer :: register, element_index
 double precision, dimension(3) :: position
 
end type atom
end module atomtype
