program test_edit_geom
  !! Test program for the module edit_geom.
  use raffle__constants, only: real32
  use raffle__geom_rw, only: basis_type
  use raffle__misc_linalg, only: modu
  use raffle__dist_calcs, only: &
       get_min_dist, &
       get_min_dist_between_point_and_atom, &
       get_min_dist_between_point_and_species, &
       get_dist_between_point_and_atom

  implicit none

  type(basis_type) :: bas, bas2
  real(real32) :: rtmp1, rtmp2
  real(real32), dimension(3) :: loc

  logical :: success = .true.


  ! Initialise silicon basis
  bas%sysname = "Silicon"
  bas%nspec = 1
  bas%natom = 2
  allocate(bas%spec(bas%nspec))
  bas%spec(1)%num = 2
  bas%spec(1)%name = 'Si'
  allocate(bas%spec(1)%atom(bas%spec(1)%num, 3))
  bas%spec(1)%atom(1, :) = [0.0, 0.0, 0.0]
  bas%spec(1)%atom(2, :) = [0.25, 0.25, 0.25]

  ! Initialise silicon lattice
  bas%lat(1,:) = [0.0, 2.14, 2.14]
  bas%lat(2,:) = [2.14, 0.0, 2.14]
  bas%lat(3,:) = [2.14, 2.14, 0.0]


  !-----------------------------------------------------------------------------
  ! Test get_min_dist
  !-----------------------------------------------------------------------------
  rtmp1 = modu(get_min_dist(bas, loc=[0.9, 0.9, 0.9], lignore_close = .true.))

  loc = [1.0, 1.0, 1.0] - [0.9, 0.9, 0.9]
  loc = loc - ceiling(loc - 0.5)
  loc = matmul(loc, bas%lat)
  rtmp2 = modu(loc)

  if ( abs(rtmp1 - rtmp2) .gt. 1.E-6 ) then
     write(0,*) 'get_min_dist failed'
     success = .false.
  end if


  !-----------------------------------------------------------------------------
  ! Test get_min_dist_between_point_and_atom
  !-----------------------------------------------------------------------------
  rtmp1 = get_min_dist_between_point_and_atom( &
       bas,  &
       loc=[0.9, 0.9, 0.9],  &
       atom=[1, 1] &
  )

  loc = [1.0, 1.0, 1.0] - [0.9, 0.9, 0.9]
  loc = loc - ceiling(loc - 0.5)
  loc = matmul(loc, bas%lat)
  rtmp2 = modu(loc)

  if ( abs(rtmp1 - rtmp2) .gt. 1.E-6 ) then
     write(0,*) 'get_min_dist_between_point_and_atom failed'
     success = .false.
  end if


  !-----------------------------------------------------------------------------
  ! Test get_min_dist_between_point_and_species
  !-----------------------------------------------------------------------------
  bas2%sysname = "Si|Ge"
  bas2%nspec = 2
  bas2%natom = 2
  allocate(bas2%spec(bas2%nspec))
  bas2%spec(1)%num = 1
  bas2%spec(1)%name = 'Si'
  allocate(bas2%spec(1)%atom(bas2%spec(1)%num, 3))
  bas2%spec(1)%atom(1, :) = [0.0, 0.0, 0.0]
  bas2%spec(2)%num = 1
  bas2%spec(2)%name = 'Ge'
  allocate(bas2%spec(2)%atom(bas2%spec(2)%num, 3))
  bas2%spec(2)%atom(1, :) = [0.25, 0.25, 0.25]

  rtmp1 = get_min_dist_between_point_and_species( &
       bas2,  &
       loc=[0.9, 0.9, 0.9],  &
       species=1 &
  )
  loc = bas2%spec(1)%atom(1,:3) - [0.9, 0.9, 0.9]
  loc = loc - ceiling(loc - 0.5)
  loc = matmul(loc, bas2%lat)
  rtmp2 = modu(loc)

  if ( abs(rtmp1 - rtmp2) .gt. 1.E-6 ) then
      write(0,*) 'get_min_dist_between_point_and_species failed'
      success = .false.
  end if

  rtmp1 = get_min_dist_between_point_and_species( &
       bas2,  &
       loc=[0.9, 0.9, 0.9],  &
       species=2 &
  )
  loc = bas2%spec(2)%atom(1,:3) - [0.9, 0.9, 0.9]
  loc = loc - ceiling(loc - 0.5)
  loc = matmul(loc, bas2%lat)
  rtmp2 = modu(loc)
 
  if ( abs(rtmp1 - rtmp2) .gt. 1.E-6 ) then
      write(0,*) 'get_min_dist_between_point_and_species failed'
      success = .false.
  end if


  !-----------------------------------------------------------------------------
  ! Test get_dist_between_point_and_atom
  !-----------------------------------------------------------------------------
  rtmp1 = get_dist_between_point_and_atom( &
         bas,  &
         loc=[0.9, 0.9, 0.9],  &
         atom=[1, 1] &
  )
  loc = bas%spec(1)%atom(1,:3) - [0.9, 0.9, 0.9]
  loc = matmul(loc, bas%lat)
  rtmp2 = modu(loc)
  if ( abs(rtmp1 - rtmp2) .gt. 1.E-6 ) then
     write(*,*) rtmp1, rtmp2
     write(0,*) 'get_dist_between_point_and_atom failed'
     success = .false.
  end if
  rtmp1 = get_dist_between_point_and_atom( &
         bas,  &
         loc=[0.9, 0.9, 0.9],  &
         atom=[1, 2] &
  )
  loc = bas%spec(1)%atom(2,:3) - [0.9, 0.9, 0.9]
  loc = matmul(loc, bas%lat)
  rtmp2 = modu(loc)
  if ( abs(rtmp1 - rtmp2) .gt. 1.E-6 ) then
   write(*,*) rtmp1, rtmp2
     write(0,*) 'get_dist_between_point_and_atom failed'
     success = .false.
  end if


  !-----------------------------------------------------------------------------
  ! check for any failed tests
  !-----------------------------------------------------------------------------
  write(*,*) "----------------------------------------"
  if(success)then
     write(*,*) 'test_dist_calcs passed all tests'
  else
     write(0,*) 'test_dist_calcs failed one or more tests'
     stop 1
  end if

end program test_edit_geom