program test_geom_utils
  !! Test program for the module edit_geom.
  use raffle__constants, only: real32
  use raffle__geom_rw, only: basis_type
  use raffle__misc_linalg, only: modu
  use raffle__geom_utils, only: basis_merge

  implicit none

  type(basis_type) :: bas, bas2, basis_merged

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
     write(0,*) &
          'basis_merge failed, number of species not equal to 2: ', &
          basis_merged%nspec
     success = .false.
  end if
  if ( basis_merged%natom .ne. 5 ) then
     write(0,*) &
          'basis_merge failed, number of atoms not equal to 5: ', &
          basis_merged%natom
     success = .false.
  end if
  if ( basis_merged%spec(1)%num .ne. 4 ) then
     write(0,*) &
          'basis_merge failed, number of atoms for species 1 not equal to 4: ', &
          basis_merged%spec(1)%num
     success = .false.
  end if
  if ( basis_merged%spec(2)%num .ne. 1 ) then
     write(0,*) &
          'basis_merge failed, number of atoms for species 2 not equal to 1: ', &
          basis_merged%spec(2)%num
     success = .false.
  end if
  if( basis_merged%spec(1)%name .ne. 'Si' ) then
     write(0,*) &
          'basis_merge failed, name of species 1 not equal to Si: ', &
          basis_merged%spec(1)%name
     success = .false.
  end if
  if( basis_merged%spec(2)%name .ne. 'O' ) then
     write(0,*) &
          'basis_merge failed, name of species 2 not equal to O: ', &
          basis_merged%spec(2)%name
     success = .false.
  end if


  !-----------------------------------------------------------------------------
  ! check for any failed tests
  !-----------------------------------------------------------------------------
  write(*,*) "----------------------------------------"
  if(success)then
     write(*,*) 'test_geom_utils passed all tests'
  else
     write(0,*) 'test_geom_utils failed one or more tests'
     stop 1
  end if

end program test_geom_utils