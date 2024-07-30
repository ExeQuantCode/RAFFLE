program test_rw_geom
  !! Test program for the module rw_geom.
  use constants, only: pi,real12
  use rw_geom, only: bas_type, geom_read, geom_write, igeom_input, igeom_output
  implicit none

  integer :: unit, iostat
  type(bas_type) :: bas1, bas2

  character(len=256) :: cwd, filename = 'test/data/POSCAR_Si'
  logical :: exist, check
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


  !-----------------------------------------------------------------------------
  ! test VASP geometry read/write
  !-----------------------------------------------------------------------------
  igeom_input = 1
  igeom_output = 1
  ! Write the geometry
  open(newunit=unit, status='scratch')
  call geom_write(unit, bas1)
  rewind(unit)
  call geom_read(unit, bas2, iostat=iostat)
  if(iostat .ne. 0) then
     write(0,*) 'Geometry read or write failed'
     success = .false.
  end if
  close(unit)
  check = compare_bas(bas1, bas2)
  if(.not.check) success = .false.



  !-----------------------------------------------------------------------------
  ! test extXYZ geometry read/write
  !-----------------------------------------------------------------------------
  bas1%energy = 12.E0
  bas1%sysname = ""
  bas2%sysname = ""
  bas2%natom = 0
  bas2%lat = 0.E0
  deallocate(bas2%spec)
  igeom_input = 6
  igeom_output = 6
  write(*,*) "LOOK", bas1%natom
  ! Write the geometry
  open(newunit=unit, status='scratch')
  call geom_write(unit, bas1)
  rewind(unit)
  call geom_read(unit, bas2, iostat=iostat)
  if(iostat .ne. 0) then
     write(0,*) 'Geometry read or write failed'
     success = .false.
  end if
  close(unit)
  check = compare_bas(bas1, bas2)
  if(.not.check) success = .false.



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


contains

  function compare_bas(bas1, bas2) result(output)
    type(bas_type), intent(in) :: bas1, bas2
    logical :: output
    output = .true.

  ! Compare the geometries
    if(any(abs(bas1%lat - bas2%lat).gt.1.E-6)) then
      write(0,*) 'Geometry read/write failed, lattice check failed'
      output = .false.
   end if
   if(bas1%sysname .ne. bas2%sysname) then
      write(0,*) 'Geometry read/write failed, system name check failed'
      output = .false.
   end if
   if(bas1%natom .ne. bas2%natom) then
      write(0,*) 'Geometry read/write failed, number of atoms check failed'
      output = .false.
   end if
   if(abs(bas1%energy - bas2%energy).gt.1.E-6) then
      write(0,*) 'Geometry read/write failed, energy check failed'
      output = .false.
   end if

  end function compare_bas

end program test_rw_geom