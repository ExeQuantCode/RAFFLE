program test_geom_rw
  !! Test program for the module geom_rw.
  use raffle__constants, only: pi,real32
  use raffle__geom_rw, only: &
       basis_type, &
       geom_read, geom_write, &
       igeom_input, igeom_output, &
       get_element_properties
  implicit none

  integer :: unit, iostat, i
  real(real32) :: mass, charge, radius
  type(basis_type) :: bas1, bas2
  class(basis_type), allocatable :: bas
  real(real32), dimension(3,3) :: atoms

  character(len=256) :: cwd, filename = 'test/data/POSCAR_Si'
  logical :: exist, check
  logical :: success = .true.
  character(len=3), dimension(118) :: element_list
  real(real32), dimension(:,:), allocatable :: positions


  ! Read the geometry
  call getcwd(cwd)
  filename = trim(cwd)//"/"//filename
  inquire(file=trim(filename), exist=exist)
  if(.not.exist) then
     write(0,*) "Filepath: ", filename
     write(0,*) 'Geometry file not found'
     success = .false.
     stop 1
  end if
  open(newunit=unit, file=trim(filename), status='old', action='read')
  call geom_read(unit, bas1, iostat=iostat)
  if(iostat .ne. 0) then
     write(0,*) 'Geometry read failed'
     success = .false.
  end if
  close(unit)


  !-----------------------------------------------------------------------------
  ! test basis_type initialisation
  !-----------------------------------------------------------------------------
  bas = basis_type(bas1)
  check = compare_bas(bas, bas1)
  if(.not.check)then
     write(0,*) 'basis_type initialisation failed'
     success = .false.
  end if


  !-----------------------------------------------------------------------------
  ! test allocate_species
  !-----------------------------------------------------------------------------
  call bas%allocate_species()
  if(bas%nspec .ne. 1) then
     write(0,*) 'allocate_species failed, nspec changed'
     write(0,*) bas%nspec
     success = .false.
  end if
  call bas%allocate_species()
  if(size(bas%spec,1) .ne. bas%nspec) then
     write(0,*) 'allocate_species failed, spec array size not equal to nspec'
     write(0,*) bas%nspec
     success = .false.
  end if
  call bas%allocate_species(2)
  if(bas%nspec .ne. 2 .or. size(bas%spec,1) .ne. 2) then
     write(0,*) 'allocate_species failed, &
          &nspec or spec array size not equal to 2'
     write(0,*) bas%nspec
     success = .false.
  end if
  atoms(1,:3) = [0.0, 0.0, 0.0]
  atoms(2,:3) = [0.5, 0.5, 0.5]
  atoms(3,:3) = [0.0, 0.0, 0.5]
  call bas%allocate_species( &
       num_species = 2, &
       species_symbols = ['Si ', 'O  '], &
       species_count = [1, 2], &
       atoms = atoms &
  )
  if(bas%nspec .ne. 2 .or. size(bas%spec,1) .ne. 2) then
     write(0,*) 'allocate_species failed, &
          &nspec or spec array size not equal to 2'
     write(0,*) bas%nspec
     success = .false.
  end if
  if(bas%spec(1)%num .ne. 1 .or. bas%spec(2)%num .ne. 2) then
     write(0,*) 'allocate_species failed, &
          &num not equal to species_count'
     write(0,*) bas%spec(1)%num, bas%spec(2)%num
     success = .false.
  end if
  if(trim(bas%spec(1)%name) .ne. 'Si' .or. trim(bas%spec(2)%name) .ne. 'O') then
     write(0,*) 'allocate_species failed, &
          &name not equal to species_symbols'
     write(0,*) bas%spec(1)%name, bas%spec(2)%name
     success = .false.
  end if
  if(any(abs(bas%spec(2)%atom(1,1:3) - 0.5).gt.1.E-6))then
     write(0,*) 'allocate_species failed, atom positions not equal to atoms'
     write(0,*) bas%spec(2)%atom(1,1:3)
     success = .false.
  end if



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
     write(*,*) 'test_geom_rw passed all tests'
  else
     write(0,*) 'test_geom_rw failed one or more tests'
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

end program test_geom_rw