program test_mod_element_utils
  use raffle__element_utils
  use raffle__geom_rw, only: get_element_properties
  use raffle__constants, only: real32
  implicit none

  ! Local variables
  class(element_type), allocatable :: element
  class(element_bond_type), allocatable :: bond

  integer :: i
  logical :: in_database
  character(len=3), dimension(:), allocatable :: element_list

  logical :: success = .true.

  ! Initialise the element database
  allocate(element_list(118))
  element_list = [ &
       'H  ', 'He ', 'Li ', 'Be ', 'B  ', 'C  ', 'N  ', 'O  ', 'F  ', 'Ne ', &
       'Na ', 'Mg ', 'Al ', 'Si ', 'P  ', 'S  ', 'Cl ', 'Ar ', 'K  ', 'Ca ', &
       'Sc ', 'Ti ', 'V  ', 'Cr ', 'Mn ', 'Fe ', 'Co ', 'Ni ', 'Cu ', 'Zn ', &
       'Ga ', 'Ge ', 'As ', 'Se ', 'Br ', 'Kr ', 'Rb ', 'Sr ', 'Y  ', 'Zr ', &
       'Nb ', 'Mo ', 'Tc ', 'Ru ', 'Rh ', 'Pd ', 'Ag ', 'Cd ', 'In ', 'Sn ', &
       'Sb ', 'Te ', 'I  ', 'Xe ', 'Cs ', 'Ba ', 'La ', 'Ce ', 'Pr ', 'Nd ', &
       'Pm ', 'Sm ', 'Eu ', 'Gd ', 'Tb ', 'Dy ', 'Ho ', 'Er ', 'Tm ', 'Yb ', &
       'Lu ', 'Hf ', 'Ta ', 'W  ', 'Re ', 'Os ', 'Ir ', 'Pt ', 'Au ', 'Hg ', &
       'Tl ', 'Pb ', 'Bi ', 'Po ', 'At ', 'Rn ', 'Fr ', 'Ra ', 'Ac ', 'Th ', &
       'Pa ', 'U  ', 'Np ', 'Pu ', 'Am ', 'Cm ', 'Bk ', 'Cf ', 'Es ', 'Fm ', &
       'Md ', 'No ', 'Lr ', 'Rf ', 'Db ', 'Sg ', 'Bh ', 'Hs ', 'Mt ', 'Ds ', &
       'Rg ', 'Cn ', 'Nh ', 'Fl ', 'Mc ', 'Lv ', 'Ts ', 'Og ' &
  ]

  allocate(element_database(size(element_list)))
  do i = 1, size(element_list)
     element_database(i)%name = element_list(i)
     call get_element_properties( &
          element_list(i), &
          mass = element_database(i)%mass, &
          charge = element_database(i)%charge, &
          radius = element_database(i)%radius &
     )
  end do


  ! Test element initialization
  element = element_type( &
       name='H  ', &
       mass=1.008_real32, &
       charge=0.0_real32, &
       energy=13.6_real32 &
  )
  call assert( &
       trim(element%name) .eq. "H", &
       "Element name initialisation failed", &
       success &
  )
  call assert( &
       abs(element%mass - 1.008_real32) .lt. 1.E-6, &
       "Element mass initialisation failed", &
       success &
  )
  call assert( &
       abs(element%charge) .lt. 1.E-6, &
       "Element charge initialisation failed", &
       success &
  )
  call assert( &
       abs(element%energy - 13.6_real32) .lt. 1.E-6, &
       "Element energy initialisation failed", &
       success &
  )

  ! Test setting element properties
  call element%set('C  ', in_database)
  call assert( &
       trim(element%name) .eq. "C", &
       "Element name setting failed", &
       success &
  )
  call assert( &
       abs(element%mass - 12.011_real32) .lt. 1.E-6, &
       "Element mass setting failed", &
       success &
  )
  call assert( &
       abs(element%charge - 6._real32) .lt. 1.E-6, &
       "Element charge setting failed", &
       success &
  )
  call assert( &
       abs(element%radius - 0.76_real32) .lt. 1.E-6, &
       "Element radius setting failed", &
       success &
  )
  call assert( &
       in_database, &
       "Element in_database setting failed", &
       success &
  )

  ! Test bond initialisation
  bond = element_bond_type(elements=['H  ', 'O  '], radius=0.96_real32)

  ! Test setting bond properties
  call bond%set('C  ', 'O  ', in_database)


  !-----------------------------------------------------------------------------
  ! check for any failed tests
  !-----------------------------------------------------------------------------
  write(*,*) "----------------------------------------"
  if(success)then
     write(*,*) 'test_elements passed all tests'
  else
     write(0,*) 'test_elements failed one or more tests'
     stop 1
  end if

contains

!###############################################################################

  subroutine assert(condition, message, success)
    implicit none
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    logical, intent(inout) :: success
    if (.not. condition) then
       write(0,*) "Test failed: ", message
       success = .false.
    end if
  end subroutine assert

end program test_mod_element_utils