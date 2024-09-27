program test_edit_geom
  !! Test program for the module edit_geom.
  use constants, only: real12
  use rw_geom, only: basis_type
  use misc_linalg, only: modu
  use edit_geom, only: &
       get_min_dist, &
       get_min_dist_between_point_and_atom, &
       basis_merge

  implicit none

  type(basis_type) :: bas, bas2, basis_merged
  real(real12) :: rtmp1, rtmp2
  real(real12), dimension(3) :: loc

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
  rtmp1 = get_min_dist_between_point_and_atom(bas, loc=[0.9, 0.9, 0.9], atom=[1, 1])

  loc = [1.0, 1.0, 1.0] - [0.9, 0.9, 0.9]
  loc = loc - ceiling(loc - 0.5)
  loc = matmul(loc, bas%lat)
  rtmp2 = modu(loc)

  if ( abs(rtmp1 - rtmp2) .gt. 1.E-6 ) then
    write(0,*) 'get_min_dist_between_point_and_atom failed'
    success = .false.
  end if


  !-----------------------------------------------------------------------------
  ! Test basis_merge
  !-----------------------------------------------------------------------------

  ! Initialise second basis
  bas2%sysname = "SiO2"
  bas2%nspec = 2
  bas2%natom = 3
  allocate(bas2%spec(bas2%nspec))
  bas2%spec(1)%num = 2
  bas2%spec(1)%name = 'Si'
  allocate(bas2%spec(1)%atom(bas2%spec(1)%num, 3))
  bas2%spec(1)%atom(1, :) = [0.5, 0.5, 0.5]
  bas2%spec(1)%atom(2, :) = [0.5, 0.0, 0.0]
  bas2%spec(2)%num = 1
  bas2%spec(2)%name = 'O'
  allocate(bas2%spec(2)%atom(bas2%spec(2)%num, 3))
  bas2%spec(2)%atom(1, :) = [0.0, 0.7, 0.0]

  ! Initialise second lattice
  bas2%lat(1,:) = [0.0, 2.14, 2.14]
  bas2%lat(2,:) = [2.14, 0.0, 2.14]
  bas2%lat(3,:) = [2.14, 2.14, 0.0]

  
  basis_merged = basis_merge(bas, bas2)

  if ( basis_merged%nspec .ne. 2 ) then
    write(0,*) 'basis_merge failed, number of species not equal to 2: ', basis_merged%nspec
    success = .false.
  end if
  if ( basis_merged%natom .ne. 5 ) then
    write(0,*) 'basis_merge failed, number of atoms not equal to 5: ', basis_merged%natom
    success = .false.
  end if
  if ( basis_merged%spec(1)%num .ne. 4 ) then
    write(0,*) 'basis_merge failed, number of atoms for species 1 not equal to 4: ', basis_merged%spec(1)%num
    success = .false.
  end if
  if ( basis_merged%spec(2)%num .ne. 1 ) then
    write(0,*) 'basis_merge failed, number of atoms for species 2 not equal to 1: ', basis_merged%spec(2)%num
    success = .false.
  end if
  if( basis_merged%spec(1)%name .ne. 'Si' ) then
    write(0,*) 'basis_merge failed, name of species 1 not equal to Si: ', basis_merged%spec(1)%name
    success = .false.
  end if
  if( basis_merged%spec(2)%name .ne. 'O' ) then
    write(0,*) 'basis_merge failed, name of species 2 not equal to O: ', basis_merged%spec(2)%name
    success = .false.
  end if


  !-----------------------------------------------------------------------------
  ! check for any failed tests
  !-----------------------------------------------------------------------------
  write(*,*) "----------------------------------------"
  if(success)then
     write(*,*) 'test_edit_geom passed all tests'
  else
     write(0,*) 'test_edit_geom failed one or more tests'
     stop 1
  end if

end program test_edit_geom