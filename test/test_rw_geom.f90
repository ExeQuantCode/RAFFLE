program test_rw_geom
  !! Test program for the module rw_geom.
  use constants, only: pi,real12
  use rw_geom, only: &
       basis_type, &
       geom_read, geom_write, &
       igeom_input, igeom_output, &
       get_element_properties
  implicit none

  integer :: unit, iostat, i
  real :: mass, charge, radius
  type(basis_type) :: bas1, bas2

  character(len=256) :: cwd, filename = 'test/data/POSCAR_Si'
  logical :: exist, check
  logical :: success = .true.
  character(len=3), dimension(118) :: element_list
  real(real12), dimension(:,:), allocatable :: positions


  ! Read the geometry
  call getcwd(cwd)
  filename = trim(cwd)//"/"//filename
  inquire(file=filename, exist=exist)
  if(iostat .ne. 0) then
     write(0,*) "Filepath: ", filename
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
  if(.not.check)then
     write(0,*) 'VASP geometry read/write failed'
     success = .false.
  end if


  !-----------------------------------------------------------------------------
  ! test extXYZ geometry read/write
  !-----------------------------------------------------------------------------
  bas1%energy = 12.E0
  bas1%sysname = ""
  call uninitialise_bas(bas2)
  igeom_input = 6
  igeom_output = 6
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
  if(.not.check)then
     write(0,*) 'extXYZ geometry read/write failed'
     success = .false.
  end if

  
  !-----------------------------------------------------------------------------
  ! test copy geometry
  !-----------------------------------------------------------------------------
  call uninitialise_bas(bas2)
  call bas2%copy(bas1)
  check = compare_bas(bas1, bas2)
  if(.not.check) success = .false.

  
  !-----------------------------------------------------------------------------
  ! test coordinate system conversion
  !-----------------------------------------------------------------------------
  allocate(positions(bas1%spec(1)%num, 3))
  do i = 1, bas1%spec(1)%num
     positions(i,:) = matmul(bas1%lat, bas1%spec(1)%atom(i,:3))
  end do
  call bas1%convert()
  if(.not.bas1%lcart) then
     write(0,*) 'Coordinate system conversion failed, lcart check failed'
     write(*,*) bas1%lcart
     success = .false.
  end if
  if(any(abs(positions - bas1%spec(1)%atom(:,:3)).gt.1.E-6)) then
     write(0,*) 'Coordinate system conversion failed, atom positions &
          &check failed'
     success = .false.
  end if

  
  !-----------------------------------------------------------------------------
  ! test element parameters
  !-----------------------------------------------------------------------------
  if(abs(bas1%spec(1)%mass - 28.085E0).gt.1.E-6) then
     write(0,*) 'Element parameters failed, mass check failed'
     write(0,*) bas1%spec(1)%mass
     success = .false.
  end if
  if(abs(bas1%spec(1)%charge - 14.0E0).gt.1.E-6) then
     write(0,*) 'Element parameters failed, charge check failed'
     write(0,*) bas1%spec(1)%charge
     success = .false.
  end if
  if(abs(bas1%spec(1)%radius - 1.11E0).gt.1.E-6) then
     write(0,*) 'Element parameters failed, radius check failed'
     write(0,*) bas1%spec(1)%radius
     success = .false.
  end if


  !-----------------------------------------------------------------------------
  ! test element properties
  !-----------------------------------------------------------------------------
  element_list = [ &
       'H  ', 'He ', &
       'Li ', 'Be ', 'B  ', 'C  ', 'N  ', 'O  ', 'F  ', 'Ne ', &
       'Na ', 'Mg ', 'Al ', 'Si ', 'P  ', 'S  ', 'Cl ', 'Ar ', &
       'K  ', 'Ca ', &
       'Sc ', 'Ti ', 'V  ', 'Cr ', 'Mn ', 'Fe ', 'Co ', 'Ni ', 'Cu ', 'Zn ', &
       'Ga ', 'Ge ', 'As ', 'Se ', 'Br ', 'Kr ', &
       'Rb ', 'Sr ', 'Y  ', &
       'Zr ', 'Nb ', 'Mo ', 'Tc ', 'Ru ', 'Rh ', 'Pd ', 'Ag ', 'Cd ', &
       'In ', 'Sn ', 'Sb ', 'Te ', 'I  ', 'Xe ', &
       'Cs ', 'Ba ', 'La ', &
       'Ce ', 'Pr ', 'Nd ', 'Pm ', 'Sm ', 'Eu ', 'Gd ', 'Tb ', 'Dy ', &
       'Ho ', 'Er ', 'Tm ', 'Yb ', 'Lu ', &
       'Hf ', 'Ta ', 'W  ', 'Re ', 'Os ', 'Ir ', 'Pt ', 'Au ', 'Hg ', &
       'Tl ', 'Pb ', 'Bi ', 'Po ', 'At ', 'Rn ', &
       'Fr ', 'Ra ', 'Ac ', &
       'Th ', 'Pa ', 'U  ', 'Np ', 'Pu ', 'Am ', 'Cm ', 'Bk ', 'Cf ', &
       'Es ', 'Fm ', 'Md ', 'No ', 'Lr ', &
       'Rf ', 'Db ', 'Sg ', 'Bh ', 'Hs ', 'Mt ', 'Ds ', &
       'Rg ', 'Cn ', 'Nh ', 'Fl ', 'Mc ', 'Lv ', 'Ts ', 'Og ' &
  ]

  do i = 1, size(element_list)
     call get_element_properties(element_list(i), &
          mass = mass, &
          charge = charge, &
          radius = radius &
     )
     if(mass.lt.1.E-6) then
        write(0,*) 'Element properties failed, mass check failed'
         write(*,*) element_list(i), mass
        success = .false.
     end if
     if(charge.lt.1.E-6) then
         write(0,*) 'Element properties failed, charge check failed'
         write(*,*) element_list(i), charge
         success = .false.
     end if
     if(radius.lt.1.E-6) then
         write(0,*) 'Element properties failed, radius check failed'
         write(*,*) element_list(i), radius
         success = .false.
     end if
  end do

  call get_element_properties('X  ', &
       mass = mass, &
       charge = charge, &
       radius = radius &
  )
   if(mass.gt.1.E-6) then
       write(0,*) 'Element properties failed, mass check failed'
       write(*,*) 'X', mass
       success = .false.
   end if
   if(charge.gt.1.E-6) then
       write(0,*) 'Element properties failed, charge check failed'
       write(*,*) 'X', charge
       success = .false.
   end if
   if(radius.gt.1.E-6) then
       write(0,*) 'Element properties failed, radius check failed'
       write(*,*) 'X', radius
       success = .false.
   end if


  !-----------------------------------------------------------------------------
  ! check for any failed tests
  !-----------------------------------------------------------------------------
  write(*,*) "----------------------------------------"
  if(success)then
     write(*,*) 'test_rw_geom passed all tests'
  else
     write(0,*) 'test_rw_geom failed one or more tests'
     stop 1
  end if


contains

  function compare_bas(bas1, bas2) result(output)
    type(basis_type), intent(in) :: bas1, bas2
    logical :: output
    output = .true.

  ! Compare the geometries
    if(any(abs(bas1%lat - bas2%lat).gt.1.E-6)) then
      write(0,*) 'Geometry read/write failed, lattice check failed'
      output = .false.
   end if
   if(bas1%sysname .ne. bas2%sysname) then
      write(0,*) 'Geometry read/write failed, system name check failed'
      write(0,*) bas1%sysname, bas2%sysname
      output = .false.
   end if
   if(bas1%natom .ne. bas2%natom) then
      write(0,*) 'Geometry read/write failed, number of atoms check failed'
      write(0,*) bas1%natom, bas2%natom
      output = .false.
   end if
   if(abs(bas1%energy - bas2%energy).gt.1.E-6) then
      write(0,*) 'Geometry read/write failed, energy check failed'
      write(0,*) bas1%energy, bas2%energy
      output = .false.
   end if

  end function compare_bas

  subroutine uninitialise_bas(bas)
    implicit none
    type(basis_type), intent(inout) :: bas

    bas%natom = 0
    bas%nspec = 0
    bas%lat = 0.E0
    bas%energy = 0.E0
    bas%sysname = ""
    bas%lcart = .false.
    bas%pbc = .true.
    deallocate(bas%spec)
    
  end subroutine uninitialise_bas

end program test_rw_geom