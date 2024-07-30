program test_rw_geom
  !! Test program for the module rw_geom.
  use constants, only: pi,real12
  use rw_geom, only: bas_type, geom_read, geom_write
  implicit none

  integer :: unit, iostat
  type(bas_type) :: bas1, bas2

  character(len=256) :: cwd, filename = 'test/data/POSCAR_Si'
  logical :: exist
  logical :: success = .true.


  ! Read the geometry
  call getcwd(cwd)
  filename = trim(cwd)//"/"//filename
  inquire(file=filename, exist=exist)
  if(iostat .ne. 0) then
     write(0,*) 'Geometry file not found'
     success = .false.
     stop 1
  end if
  open(newunit=unit, file=filename, status='old')
  call geom_read(unit, bas1, iostat=iostat)
  if(iostat .ne. 0) then
     write(0,*) 'Geometry read failed'
     success = .false.
  end if
  close(unit)

  ! Write the geometry
  open(newunit=unit, status='scratch')
  call geom_write(unit, bas1)
  rewind(unit)
  call geom_read(unit, bas2, iostat=iostat)
  if(iostat .ne. 0) then
     write(0,*) 'Geometry write failed'
     success = .false.
  end if
  close(unit)

  ! Compare the geometries
  if (any(abs(bas1%lat - bas2%lat).gt.1.E-6)) then
     write(0,*) 'Geometry read/write failed'
     success = .false.
  end if

  !-----------------------------------------------------------------------------
  ! check for any failed tests
  !-----------------------------------------------------------------------------
  write(*,*) "----------------------------------------"
  if(success)then
     write(*,*) 'test_accuracy passed all tests'
  else
     write(0,*) 'test_accuracy failed one or more tests'
     stop 1
  end if


end program test_rw_geom