module raffle__geom_rw
  !! Module to store, read and write geometry files
  !!
  !! This module contains the procedures to read and write geometry files.
  !! It also contains the derived types used to store the geometry data.
  use raffle__constants, only: pi,real32
  use raffle__io_utils, only: stop_program, print_warning
  use raffle__misc, only: to_upper, to_lower, jump, icount, strip_null
  use raffle__misc_linalg, only: modu, inverse_3x3
  implicit none


  private

  public :: igeom_input, igeom_output
  public :: basis_type, species_type
  public :: geom_read, geom_write
  public :: get_element_properties


  integer :: igeom_input = 1
  !! geometry input file format
  !! 1 = VASP
  ! 2 = CASTEP
  !! 3 = Quantum Espresso
  !! 4 = CRYSTAL
  !! 5 = XYZ
  !! 6 = extended XYZ
  integer :: igeom_output = 1
  !! geometry output file format

  type :: species_type
     !! Derived type to store information about a species/element.
     real(real32), allocatable ,dimension(:,:) :: atom
     !! The atomic positions of the species.
     real(real32) :: mass
     !! The mass of the species.
     real(real32) :: charge
     !! The charge of the species.
     real(real32) :: radius
     !! The radius of the species.
     character(len=3) :: name
     !! The name of the species.
     integer :: num
     !! The number of atoms of this species.
  end type species_type
  type :: basis_type
     !! Derived type to store information about a basis.
     type(species_type), allocatable, dimension(:) :: spec
     !! Information about each species in the basis.
     integer :: nspec = 0
     !! The number of species in the basis.
     integer :: natom = 0
     !! The number of atoms in the basis.
     real(real32) :: energy = 0._real32
     !! The energy of the basis.
     real(real32) :: lat(3,3) = 0._real32
     !! The lattice vectors of the basis.
     logical :: lcart = .false.
     !! Boolean whether the basis is in cartesian coordinates.
     logical, dimension(3) :: pbc = .true.
     !! Boolean whether the basis has periodic boundary conditions.
     character(len=128) :: sysname = "default"
     !! The name of the system.
   contains
     procedure, pass(this) :: allocate_species
     !! Procedure to allocate the species in the basis.
     procedure, pass(this) :: convert
     !! Procedure to convert the basis to cartesian coordinates.
     procedure, pass(this) :: copy
     !! Procedure to copy the basis.
     procedure, pass(this) :: get_lattice_constants
     !! Procedure to get the lattice constants of the basis.
     procedure, pass(this) :: remove_atom
     !! Procedure to remove an atom from the basis.
     procedure, pass(this) :: remove_atoms
     !! Procedure to remove atoms from the basis.
  end type basis_type


  interface basis_type
     module function init_basis_type(basis) result(output)
       !! Initialise the basis type.
       type(basis_type), intent(in), optional :: basis
       !! Optional. Basis to copy.
       type(basis_type) :: output
       !! The basis to initialise.
     end function init_basis_type
  end interface basis_type



contains

!###############################################################################
  module function init_basis_type(basis) result(output)
    !! Initialise the basis type.
    implicit none

    ! Arguments
    type(basis_type), intent(in), optional :: basis
    !! Optional. Basis to copy.
    type(basis_type) :: output
    !! The basis to initialise.

    if(present(basis)) call output%copy(basis)

  end function init_basis_type
!###############################################################################


!###############################################################################
  subroutine allocate_species( &
       this, num_species, &
       species_symbols, species_count, atoms )
    !! Allocate the species in the basis.
    implicit none

    ! Arguments
    class(basis_type), intent(inout) :: this
    !! Parent. The basis to allocate the species in.
    integer, intent(in), optional :: num_species
    !! Optional. The number of species in the basis.
    character(3), dimension(:), intent(in), optional :: species_symbols
    !! Optional. The symbols of the species.
    integer, dimension(:), intent(in), optional :: species_count
    !! Optional. The number of atoms of each species.
    real(real32), dimension(:,:), intent(in), optional :: atoms
    !! Optional. The atomic positions of the species.
    
    ! Local variables
    integer :: i, istart, iend
    !! Loop index.

    if(present(num_species)) this%nspec = num_species

    if(allocated(this%spec)) deallocate(this%spec)
    allocate(this%spec(this%nspec))

    species_check: if(present(species_symbols))then
       if(size(species_symbols).ne.this%nspec) exit species_check
       this%spec(:)%name = species_symbols
    end if species_check

    natom_check: if(present(species_count))then
       if(size(species_count).ne.this%nspec) exit natom_check
       this%spec(:)%num = species_count
       istart = 1
       do i = 1, this%nspec
          iend = istart + this%spec(i)%num - 1
          allocate(this%spec(i)%atom(this%spec(i)%num,3))
          if(present(atoms))then
             this%spec(i)%atom = atoms(istart:iend,:3)
          end if
          istart = iend + 1
       end do
    end if natom_check
    
    do i = 1, this%nspec
       call get_element_properties( &
            this%spec(i)%name, &
            mass = this%spec(i)%mass, &
            charge = this%spec(i)%charge, &
            radius = this%spec(i)%radius )
    end do
       
  end subroutine allocate_species
!###############################################################################

  
!###############################################################################
  subroutine geom_read(UNIT, basis, length, iostat)
    !! Read geometry from a file.
    implicit none

    ! Arguments
    integer, intent(in) :: UNIT
    !! The unit number of the file to read from.
    type(basis_type), intent(out) :: basis
    !! The basis to read the geometry into.
    integer, optional, intent(in) :: length
    !! Optional. The dimension of the basis atom positions.
    integer, optional, intent(out) :: iostat
    !! Optional. The I/O status of the read.

    ! Local variables
    integer :: i
    !! Loop index.
    integer :: length_
    !! The dimension of the basis atom positions.
    integer :: iostat_
    !! The I/O status of the read.


    length_ = 3
    iostat_ = 0
    if(present(length)) length_=length

    select case(igeom_input)
    case(1)
       call VASP_geom_read(UNIT, basis, length_, iostat_)
    case(2)
       call CASTEP_geom_read(UNIT, basis, length_)
    case(3)
       call QE_geom_read(UNIT, basis, length_)
    case(4)
       call stop_program("Not yet set up for CRYSTAL")
       return
    case(5)
       call XYZ_geom_read(UNIT, basis, length_, iostat_)
       call print_warning("XYZ file format does not contain lattice data")
    case(6)
       call extXYZ_geom_read(UNIT, basis, length_, iostat_)
    end select
    if(iostat_.ne.0) then
       if(present(iostat)) iostat = iostat_
       return
    else
       if(present(iostat)) iostat = 0
    end if
    if(length_.eq.4)then
       do i = 1, basis%nspec
          basis%spec(i)%atom(:,4)=1._real32
       end do
    end if
    do i = 1, basis%nspec
       call get_element_properties( &
            basis%spec(i)%name, &
            mass = basis%spec(i)%mass, &
            charge = basis%spec(i)%charge, &
            radius = basis%spec(i)%radius )
    end do

  end subroutine geom_read
!###############################################################################


!###############################################################################
  subroutine geom_write(UNIT, basis)
    !! Write geometry to a file.
    implicit none

    ! Arguments
    integer, intent(in) :: UNIT
    !! The unit number of the file to write to.
    class(basis_type), intent(in) :: basis
    !! The basis to write the geometry from.

    ! MAKE IT CHANGE HERE IF USER SPECIFIES LCART OR NOT
    ! AND GIVE IT THE CASTEP AND QE OPTION OF LABC !

    select case(igeom_output)
    case(1)
       call VASP_geom_write(UNIT,basis)
    case(2)
       call CASTEP_geom_write(UNIT,basis)
    case(3)
       call QE_geom_write(UNIT,basis)
    case(4)
       call stop_program("ERROR: Not yet set up for CRYSTAL")
       return
    case(5)   
       call XYZ_geom_write(UNIT,basis)
    case(6)
       call extXYZ_geom_write(UNIT,basis)
    end select

  end subroutine geom_write
!###############################################################################


!###############################################################################
  subroutine VASP_geom_read(UNIT, basis, length, iostat)
    !! Read the structure in vasp poscar style format.
    implicit none

    ! Arguments
    integer, intent(in) :: UNIT
    !! The unit number of the file to read from.
    type(basis_type), intent(out) :: basis
    !! The basis to read the geometry into.
    integer, intent(in), optional :: length
    !! Optional. The dimension of the basis atom positions.
    integer, intent(out), optional :: iostat
    !! Optional. The I/O status of the read.

    integer :: Reason
    !! The I/O status of the read.
    integer :: pos, count
    !! Temporary integer variables.
    real(real32) :: scal
    !! The scaling factor of the lattice.
    character(len=100) :: lspec
    !! The species names and number of each atomic species.
    character(len=1024) :: buffer
    !! Temporary character variable.
    integer :: i, j, k
    !! Loop index.
    integer :: length_
    !! The dimension of the basis atom positions.
    integer :: iostat_
    !! The I/O status of the read.


    length_ = 3
    iostat_ = 0
    !---------------------------------------------------------------------------
    ! determine dimension of basis (include translation dimension for symmetry?)
    !---------------------------------------------------------------------------
    if(present(length)) length_ = length


    !---------------------------------------------------------------------------
    ! read system name
    !---------------------------------------------------------------------------
    read(UNIT,'(A)',iostat=Reason) basis%sysname
    if(Reason.ne.0)then
       write(0,'("ERROR: The file is not in POSCAR format.")')
       write(0,*) "Expected system name, got: ",trim(basis%sysname)
       iostat_ = 1
       if(present(iostat)) iostat = iostat_
       return
    end if
    read(UNIT,*) scal


    !---------------------------------------------------------------------------
    ! read lattice
    !---------------------------------------------------------------------------
    do i = 1, 3
       read(UNIT,*) (basis%lat(i,j),j=1,3)
    end do
    basis%lat=scal*basis%lat


    !---------------------------------------------------------------------------
    ! read species names and number of each atomic species
    !---------------------------------------------------------------------------
    read(UNIT,'(A)') lspec
    basis%nspec = icount(lspec)
    allocate(basis%spec(basis%nspec))
    if(verify(lspec,' 0123456789').ne.0) then
       count=0;pos=1
       speccount: do
          i=verify(lspec(pos:), ' ')
          if (i.eq.0) exit speccount
          count=count+1
          pos=i+pos-1
          i=scan(lspec(pos:), ' ')
          if (i.eq.0) exit speccount
          basis%spec(count)%name=lspec(pos:pos+i-1)
          pos=i+pos-1
       end do speccount

       read(UNIT,*) (basis%spec(j)%num,j=1,basis%nspec)
    else !only numbers
       do count = 1, basis%nspec
          write(basis%spec(count)%name,'(I0)') count
       end do
       read(lspec,*) (basis%spec(j)%num,j=1,basis%nspec)
    end if


    !---------------------------------------------------------------------------
    ! determines whether input basis is in direct or cartesian coordinates
    !---------------------------------------------------------------------------
    basis%lcart=.false.
    read(UNIT,'(A)') buffer
    buffer = to_lower(buffer)
    if(verify(trim(buffer),'direct').eq.0) basis%lcart=.false.
    if(verify(trim(buffer),'cartesian').eq.0) basis%lcart=.true.


    !---------------------------------------------------------------------------
    ! read basis
    !---------------------------------------------------------------------------
    do i = 1, basis%nspec
       allocate(basis%spec(i)%atom(basis%spec(i)%num,length_))
       basis%spec(i)%atom(:,:)=0._real32
       do j = 1, basis%spec(i)%num
          read(UNIT,*) (basis%spec(i)%atom(j,k),k=1,3)
       end do
    end do


    !---------------------------------------------------------------------------
    ! convert basis if in cartesian coordinates
    !---------------------------------------------------------------------------
    if(basis%lcart) call basis%convert()


    !---------------------------------------------------------------------------
    ! normalise basis to between 0 and 1 in direct coordinates
    !---------------------------------------------------------------------------
    do i = 1, basis%nspec
       do j = 1, basis%spec(i)%num
          do k = 1, 3
             basis%spec(i)%atom(j,k)=&
                  basis%spec(i)%atom(j,k)-floor(basis%spec(i)%atom(j,k))
          end do
       end do
    end do
    basis%natom=sum(basis%spec(:)%num)
 
    if(present(iostat)) iostat = iostat_
 
  end subroutine VASP_geom_read
!###############################################################################


!###############################################################################
  subroutine VASP_geom_write(UNIT, basis, cartesian)
    !! Write the structure in vasp poscar style format.
    implicit none

    ! Arguments
    integer, intent(in) :: UNIT
    !! The unit number of the file to write to.
    class(basis_type), intent(in) :: basis
    !! The basis to write the geometry from.
    logical, intent(in), optional :: cartesian
    !! Optional. Whether to write the basis in cartesian coordinates.

    ! Local variables
    integer :: i,j
    !! Loop index.
    character(100) :: fmt
    !! Format string.
    character(10) :: string
    !! String to determine whether to write in direct or cartesian coordinates.


    string="Direct"
    if(present(cartesian))then
       if(cartesian) string="Cartesian"
    end if

    write(UNIT,'(A)') trim(adjustl(basis%sysname))
    write(UNIT,'(F15.9)') 1._real32
    do i = 1, 3
       write(UNIT,'(3(F15.9))') basis%lat(i,:)
    end do
    write(fmt,'("(",I0,"(A,1X))")') basis%nspec
    write(UNIT,trim(adjustl(fmt))) (adjustl(basis%spec(j)%name),j=1,basis%nspec)
    write(fmt,'("(",I0,"(I0,5X))")') basis%nspec
    write(UNIT,trim(adjustl(fmt))) (basis%spec(j)%num,j=1,basis%nspec)
    write(UNIT,'(A)') trim(adjustl(string))
    do i = 1, basis%nspec
       do j = 1, basis%spec(i)%num
          write(UNIT,'(3(F15.9))') basis%spec(i)%atom(j,1:3)
       end do
    end do

  end subroutine VASP_geom_write
!###############################################################################


!###############################################################################
  subroutine QE_geom_read(UNIT,basis,length)
    !! Read the structure in Quantum Espresso style format.
    implicit none

    ! Arguments
    integer, intent(in) :: UNIT
    !! The unit number of the file to read from.
    type(basis_type), intent(out) :: basis
    !! The basis to read the geometry into.
    integer, intent(in), optional :: length
    !! Optional. The dimension of the basis atom positions.

    ! Local variables
    integer :: Reason
    !! The I/O status of the read.
    integer :: i, j, k, iline
    !! Loop index.
    integer :: length_ = 3
    !! The dimension of the basis atom positions.
    integer, dimension(1000) :: tmp_natom
    !! Temporary array to store the number of atoms of each species.
    real(real32), dimension(3) :: tmpvec
    !! Temporary array to store the atomic positions.
    character(len=3) :: ctmp
    !! Temporary character variable.
    character(256) :: stop_msg
    !! Error message.
    character(len=3), dimension(1000) :: tmp_spec
    !! Temporary array to store the species names.
    character(len=1024) :: buffer, buffer2
    !! Temporary character variables.


    !---------------------------------------------------------------------------
    ! determine dimension of basis (include translation dimension for symmetry?)
    !---------------------------------------------------------------------------
    if(present(length)) length_ = length
    basis%lcart = .false.
    basis%sysname = "Converted_from_geom_file"


    !---------------------------------------------------------------------------
    ! read lattice
    !---------------------------------------------------------------------------
    rewind UNIT
    cellparam: do
       read(UNIT,'(A)',iostat=Reason) buffer
       if(Reason.ne.0)then
          call stop_program( &
               "An issue with the QE input file format has been encountered." &
          )
          return
       end if
       if(index(trim(buffer),"ibrav").ne.0)then
          write(stop_msg,*) &
               "Internal error in QE_geom_read" // &
               achar(13) // achar(10) // &
               "  Subroutine not yet set up to read IBRAV lattices"
          call stop_program(stop_msg)
          return
       end if
       if(verify("CELL_PARAMETERS",buffer).eq.0) then
          exit cellparam
       end if
    end do cellparam
    do i = 1, 3
       read(UNIT,*) (basis%lat(i,j),j=1,3)
    end do


    !---------------------------------------------------------------------------
    ! determines whether input basis is in direct or cartesian coordinates
    !---------------------------------------------------------------------------
    iline=0
    rewind UNIT
    basfind: do
       read(UNIT,'(A)',iostat=Reason) buffer
       iline=iline+1
       if(verify("ATOMIC_POSITIONS",buffer).eq.0)then
          backspace(UNIT)
          read(UNIT,*) buffer,buffer2
          if(verify("crystal",buffer2).eq.0) basis%lcart = .false.
          if(verify("angstrom",buffer2).eq.0) basis%lcart = .true.
          exit basfind
       end if
    end do basfind
 
 
    !---------------------------------------------------------------------------
    ! read basis
    !---------------------------------------------------------------------------
    basis%natom = 0
    basis%nspec = 0
    tmp_natom   = 1
    basread: do
       read(UNIT,'(A)',iostat=Reason) buffer
       read(buffer,*) ctmp
       if(Reason.ne.0) exit
       if(trim(ctmp).eq.'') exit
       if(verify(buffer,' 0123456789').eq.0) exit
       basis%natom = basis%natom + 1
       if(.not.any(tmp_spec(1:basis%nspec).eq.ctmp))then
          basis%nspec = basis%nspec + 1
          tmp_spec(basis%nspec) = ctmp
       else
          where(tmp_spec(1:basis%nspec).eq.ctmp)
             tmp_natom(1:basis%nspec) = tmp_natom(1:basis%nspec) + 1
          end where
       end if
    end do basread

    allocate(basis%spec(basis%nspec))
    basis%spec(1:basis%nspec)%name = tmp_spec(1:basis%nspec)
    do i = 1, basis%nspec
       basis%spec(i)%num = 0
       allocate(basis%spec(i)%atom(tmp_natom(i),length_))
    end do

    call jump(UNIT,iline)
    basread2: do i = 1, basis%natom
       read(UNIT,*,iostat=Reason) ctmp,tmpvec(1:3)
       do j = 1, basis%nspec
          if(basis%spec(j)%name.eq.ctmp)then
             basis%spec(j)%num = basis%spec(j)%num + 1
             basis%spec(j)%atom(basis%spec(j)%num,1:3) = tmpvec(1:3)
             exit
          end if
       end do
    end do basread2


    !---------------------------------------------------------------------------
    ! convert basis if in cartesian coordinates
    !---------------------------------------------------------------------------
    if(basis%lcart) call basis%convert()


    !---------------------------------------------------------------------------
    ! normalise basis to between 0 and 1 in direct coordinates
    !---------------------------------------------------------------------------
    do i = 1, basis%nspec
       do j = 1, basis%spec(i)%num
          do k = 1, 3
             basis%spec(i)%atom(j,k) = &
                  basis%spec(i)%atom(j,k) - floor( basis%spec(i)%atom(j,k) )
          end do
       end do
    end do
    basis%natom=sum(basis%spec(:)%num)

  end subroutine QE_geom_read
!###############################################################################


!###############################################################################
  subroutine QE_geom_write(UNIT, basis, cartesian)
    !! Write the structure in Quantum Espresso style format.
    implicit none

    ! Arguments
    integer, intent(in) :: UNIT
    !! The unit number of the file to write to.
    class(basis_type), intent(in) :: basis
    !! The basis to write the geometry from.
    logical, intent(in), optional :: cartesian
    !! Optional. Whether to write the basis in cartesian coordinates.

    ! Local variables
    integer :: i,j
    !! Loop index.
    character(10) :: string
    !! String to determine whether to write in crystal or angstrom coordinates.


    string="crystal"
    if(present(cartesian))then
       if(cartesian) string="angstrom"
    end if


    write(UNIT,'("CELL_PARAMETERS angstrom")')
    do i = 1, 3
       write(UNIT,'(3(F15.9))') basis%lat(i,:)
    end do
    write(UNIT,'("ATOMIC_SPECIES")')
    do i = 1, basis%nspec
       write(UNIT,'(A)') trim(adjustl(basis%spec(i)%name))
    end do
    write(UNIT,'("ATOMIC_POSITIONS",1X,A)') trim(adjustl(string))
    do i = 1, basis%nspec
       do j = 1, basis%spec(i)%num
          write(UNIT,'(A5,1X,3(F15.9))') &
               basis%spec(i)%name,basis%spec(i)%atom(j,1:3)
       end do
    end do

  end subroutine QE_geom_write
!###############################################################################


!###############################################################################
  subroutine CASTEP_geom_read(UNIT, basis, length)
    !! Read the structure in CASTEP style format.
    implicit none

    ! Arguments
    integer, intent(in) :: UNIT
    !! The unit number of the file to read from.
    type(basis_type), intent(out) :: basis
    !! The basis to read the geometry into.
    integer, intent(in), optional :: length
    !! Optional. The dimension of the basis atom positions.

    ! Local variables
    integer :: Reason
    !! The I/O status of the read.
    integer :: i, j, k, iline
    !! Loop index.
    integer :: length_ = 3
    !! The dimension of the basis atom positions.
    integer :: itmp1
    !! Temporary integer variable.
    character(len=3) :: ctmp
    !! Temporary character variable.
    character(len=20) :: units
    !! Units of the lattice vectors.
    character(len=200) :: buffer, store
    !! Temporary character variables.
    logical :: labc
    !! Logical variable to determine whether the lattice is in abc or 
    !! cartesian coordinates.
    integer, dimension(1000) :: tmp_natom
    !! Temporary array to store the number of atoms of each species.
    real(real32), dimension(3) :: abc, angle, dvtmp1
    !! Temporary arrays to store the lattice vectors.
    character(len=3), dimension(1000) :: tmp_spec
    !! Temporary array to store the species names.

   
    !---------------------------------------------------------------------------
    ! determine dimension of basis (include translation dimension for symmetry?)
    !---------------------------------------------------------------------------
    if(present(length)) length_ = length


    !---------------------------------------------------------------------------
    ! reading loop of file
    !---------------------------------------------------------------------------
    tmp_spec = ""
    tmp_natom = 0
    iline = 0
    labc = .true.
    basis%sysname = "from CASTEP"
    rewind(UNIT)
    readloop: do
       iline=iline+1
       read(UNIT,'(A)',iostat=Reason) buffer
       if(Reason.ne.0) exit
       buffer=to_upper(buffer)
       if(scan(trim(adjustl(buffer)),'%').ne.1) cycle readloop
       if(index(trim(adjustl(buffer)),'%END').eq.1) cycle readloop
       read(buffer,*) store, buffer
       if(trim(buffer).eq.'') cycle readloop
       !------------------------------------------------------------------------
       ! read lattice
       !------------------------------------------------------------------------
       lattice_if: if(index(trim(buffer),"LATTICE").eq.1)then
          if(index(trim(buffer),"ABC").ne.0) labc = .true.
          if(index(trim(buffer),"CART").ne.0) labc = .false.
          store = ""
          itmp1 = 0
          lattice_loop: do 
             itmp1 = itmp1 + 1
             read(UNIT,'(A)',iostat=Reason) buffer
             if(Reason.ne.0) exit lattice_loop
             if(scan(trim(adjustl(buffer)),'%').eq.1) exit lattice_loop
             if(itmp1.eq.5)then
                call stop_program( &
                     "Too many lines in LATTICE block of structure file" &
                )
                return
             end if
             store=trim(store)//" "//trim(buffer)
          end do lattice_loop
          iline=iline+itmp1

          if(labc)then
             read(store,*) units,(abc(i),i=1,3), (angle(j),j=1,3)
             basis%lat = convert_abc_to_lat(abc,angle,.false.)
          else
             read(store,*) units,(basis%lat(i,:),i=1,3)
          end if
          cycle readloop
       end if lattice_if

       !------------------------------------------------------------------------
       ! read basis
       !------------------------------------------------------------------------
       basis_if: if(index(trim(buffer),"POSITIONS").eq.1) then
          if(index(trim(buffer),"ABS").ne.0) basis%lcart=.true.
          if(index(trim(buffer),"FRAC").ne.0) basis%lcart=.false.
          itmp1 = 0
          basis_loop1: do
             read(UNIT,'(A)',iostat=Reason) buffer
             if(Reason.ne.0) exit basis_loop1
             if(scan(trim(adjustl(buffer)),'%').eq.1) exit basis_loop1
             read(buffer,*) ctmp
             if(trim(ctmp).eq.'') exit
             if(verify(buffer,' 0123456789').eq.0) exit
             basis%natom = basis%natom + 1
             if(.not.any(tmp_spec(1:basis%nspec).eq.ctmp))then
                basis%nspec = basis%nspec+1
                tmp_natom(basis%nspec) = 1
                tmp_spec(basis%nspec)  = ctmp
             else
                where(tmp_spec(1:basis%nspec).eq.ctmp)
                   tmp_natom(1:basis%nspec) = tmp_natom(1:basis%nspec) + 1
                end where
             end if
          end do basis_loop1

          allocate(basis%spec(basis%nspec))
          basis%spec(1:basis%nspec)%name = tmp_spec(1:basis%nspec)
          do i = 1, basis%nspec
             basis%spec(i)%num = 0
             allocate(basis%spec(i)%atom(tmp_natom(i),length_))
          end do

          call jump(UNIT,iline)
          basis_loop2: do i = 1, basis%natom
             read(UNIT,'(A)',iostat=Reason) buffer
             if(Reason.ne.0)then
                call stop_program("Internal error in assigning the basis")
                return
             end if
             read(buffer,*) ctmp,dvtmp1(1:3)
             species_loop: do j = 1, basis%nspec
                if(basis%spec(j)%name.eq.ctmp)then
                   basis%spec(j)%num = basis%spec(j)%num + 1
                   basis%spec(j)%atom(basis%spec(j)%num,1:3) = dvtmp1(1:3)
                   exit species_loop
                end if
             end do species_loop
          end do basis_loop2

       end if basis_if
    end do readloop


    !---------------------------------------------------------------------------
    ! convert basis if in cartesian coordinates
    !---------------------------------------------------------------------------
    if(basis%lcart) call basis%convert()


    !---------------------------------------------------------------------------
    ! normalise basis to between 0 and 1 in direct coordinates
    !---------------------------------------------------------------------------
    do i = 1, basis%nspec
       do j = 1, basis%spec(i)%num
          do k = 1, 3
             basis%spec(i)%atom(j,k) = &
                  basis%spec(i)%atom(j,k) - floor( basis%spec(i)%atom(j,k) )
          end do
       end do
    end do
    basis%natom=sum(basis%spec(:)%num)

  end subroutine CASTEP_geom_read
!###############################################################################


!###############################################################################
  subroutine CASTEP_geom_write(UNIT, basis, labc, cartesian)
    !! Write the structure in CASTEP style format.
    implicit none

    ! Arguments
    integer :: UNIT
    !! The unit number of the file to write to.
    class(basis_type), intent(in) :: basis
    !! The basis to write the geometry from.
    logical, intent(in), optional :: labc
    !! Optional. Boolean whether to write the lattice in abc format.
    logical, intent(in), optional :: cartesian
    !! Optional. Boolean whether to write basis in cartesian coordinates.

    ! Local variables
    integer :: i, j
    !! Loop index.
    real(real32), dimension(2,3) :: abc_angle
    !! Temporary arrays to store the lattice vectors.
    character(4) :: string_lat, string_bas
    !! Strings specifying lattice and basis format
    character(len=256) :: stop_msg
    !! Error message.


    string_lat="CART"
    if(present(labc))then
       if(labc) string_lat="ABC"
    end if

    string_bas="FRAC"
    if(present(cartesian))then
       if(cartesian)then
          string_bas="ABS"
          write(stop_msg,*) &
               "Internal error in CASTEP_geom_write" // &
               achar(13) // achar(10) // &
               "  Subroutine not yet set up to output cartesian coordinates"
          call stop_program(stop_msg)
          return
       end if
    end if

    write(UNIT,'("%block LATTICE_",A)') trim(string_lat)
    write(UNIT,'("ang")')
    if(present(labc))then
       if(labc)then
          abc_angle = convert_lat_to_abc(basis%lat)
          write(UNIT,'(3(F15.9))') abc_angle(1,:)
          write(UNIT,'(3(F15.9))') abc_angle(2,:)
          goto 10
       end if
    end if
    do i = 1, 3
       write(UNIT,'(3(F15.9))') basis%lat(i,:)
    end do

10  write(UNIT,'("%endblock LATTICE_",A)') trim(string_lat)

    write(UNIT,*)
    write(UNIT,'("%block POSITIONS_",A)') trim(string_bas)
    do i = 1, basis%nspec
       do j = 1, basis%spec(i)%num
          write(UNIT,'(A5,1X,3(F15.9))') &
               basis%spec(i)%name,basis%spec(i)%atom(j,1:3)
       end do
    end do
    write(UNIT,'("%endblock POSITIONS_",A)') trim(string_bas)

  end subroutine CASTEP_geom_write
!###############################################################################


!###############################################################################
  subroutine XYZ_geom_read(UNIT, basis, length, iostat)
    !! Read the structure in xyz style format.
    implicit none

    ! Arguments
    integer, intent(in) :: UNIT
    !! The unit number of the file to read from.
    type(basis_type), intent(out) :: basis
    !! The basis to read the geometry into.
    integer, intent(in), optional :: length
    !! Optional. The dimension of the basis atom positions.
    integer, intent(out), optional :: iostat
    !! Optional. The I/O status of the read.

    ! Local variables
    integer :: Reason
    !! The I/O status of the read.
    integer :: i, j
    !! Loop index.
    integer, allocatable, dimension(:) :: tmp_num
    !! Temporary array to store the number of atoms of each species.
    real(real32), dimension(3) :: vec
    !! Temporary array to store the atomic positions.
    real(real32), allocatable, dimension(:,:,:) :: tmp_bas
    !! Temporary array to store the atomic positions.
    character(len=3) :: ctmp
    !! Temporary character variable.
    character(len=3), allocatable, dimension(:) :: tmp_spec
    !! Temporary array to store the species names.
    integer :: length_
    !! The dimension of the basis atom positions.
    integer :: iostat_
    !! The I/O status of the read.


    length_ = 3
    iostat_ = 0
    if(present(length)) length_ = length


    read(UNIT,*,iostat=Reason) basis%natom
    if(Reason.ne.0)then
       write(0,'("ERROR: The file is not in xyz format.")')
       iostat_ = 1
       if(present(iostat)) iostat = iostat_
       return
    end if
    read(UNIT,'(A)',iostat=Reason) basis%sysname


    !---------------------------------------------------------------------------
    ! read basis
    !---------------------------------------------------------------------------
    allocate(tmp_spec(basis%natom))
    allocate(tmp_num(basis%natom))
    allocate(tmp_bas(basis%natom,basis%natom,length_))
    tmp_num(:) = 0
    tmp_spec = ""
    tmp_bas = 0
    basis%nspec = 0
    do i = 1, basis%natom
       read(UNIT,*,iostat=Reason) ctmp,vec(1:3)
       if(.not.any(tmp_spec(1:basis%nspec).eq.ctmp))then
          basis%nspec = basis%nspec + 1
          tmp_spec(basis%nspec) = ctmp
          tmp_bas(basis%nspec,1,1:3) = vec(1:3)
          tmp_num(basis%nspec) = 1
       else
          checkspec: do j = 1, basis%nspec
             if(tmp_spec(j).eq.ctmp)then
                tmp_num(j) = tmp_num(j)+1
                tmp_bas(j,tmp_num(j),1:3) = vec(1:3)
                exit checkspec
             end if
          end do checkspec
       end if
    end do
 
 
    !---------------------------------------------------------------------------
    ! move basis from temporary basis to main basis.
    ! done to allow for correct allocation of number of and per species
    !---------------------------------------------------------------------------
    allocate(basis%spec(basis%nspec))
    do i = 1, basis%nspec
       basis%spec(i)%name = tmp_spec(i)
       basis%spec(i)%num  = tmp_num(i)
       allocate(basis%spec(i)%atom(tmp_num(i),length_))
       basis%spec(i)%atom(:,:) = 0
       basis%spec(i)%atom(1:tmp_num(i),1:3) = tmp_bas(i,1:tmp_num(i),1:3)
    end do
 
    if(present(iostat)) iostat = iostat_
 
  end subroutine XYZ_geom_read
!###############################################################################


!###############################################################################
  subroutine XYZ_geom_write(UNIT,basis)
    !! Write the structure in xyz style format.
    implicit none

    ! Arguments
    integer, intent(in) :: UNIT
    !! The unit number of the file to write to.
    class(basis_type), intent(in) :: basis
    !! The basis to write the geometry from.

    ! Local variables
    integer :: i, j
    !! Loop index.


    write(UNIT,'("I0")') basis%natom
    write(UNIT,'("A")') basis%sysname
    do i = 1, basis%nspec
       do j = 1, basis%spec(i)%num
          write(UNIT,'(A5,1X,3(F15.9))') &
               basis%spec(i)%name,basis%spec(i)%atom(j,1:3)
       end do
    end do

  end subroutine XYZ_geom_write
!###############################################################################


!###############################################################################
  subroutine extXYZ_geom_read(UNIT, basis, length, iostat)
    !! Read the structure in extended xyz style format.
    implicit none

    ! Arguments
    integer, intent(in) :: UNIT
    !! The unit number of the file to read from.
    type(basis_type), intent(out) :: basis
    !! The basis to read the geometry into.
    integer, intent(in), optional :: length
    !! Optional. The dimension of the basis atom positions.
    integer, intent(out), optional :: iostat
    !! Optional. The I/O status of the read.

    ! Local variables
    integer :: Reason
    !! The I/O status of the read.
    integer :: i, j
    !! Loop index.
    integer :: index1, index2
    !! Index variables.
    integer, allocatable, dimension(:) :: tmp_num
    !! Temporary array to store the number of atoms of each species.
    real(real32), dimension(3) :: vec
    !! Temporary array to store the atomic positions.
    real(real32), allocatable, dimension(:,:,:) :: tmp_bas
    !! Temporary array to store the atomic positions.
    character(len=3) :: ctmp
    !! Temporary character variable.
    character(len=3), allocatable, dimension(:) :: tmp_spec
    !! Temporary array to store the species names.
    character(len=1024) :: buffer
    !! Temporary character variable.
    integer :: length_ = 3
    !! The dimension of the basis atom positions.
    integer :: iostat_ = 0
    !! The I/O status of the read.


    basis%lcart=.true.
    if(present(length)) length_ = length


    !---------------------------------------------------------------------------
    ! read system information
    !---------------------------------------------------------------------------
    read(UNIT,*,iostat=Reason) basis%natom
    if(Reason.ne.0)then
       write(0,'("ERROR: The file is not in xyz format.")')
       iostat_ = 1
       if(present(iostat)) iostat = iostat_
       return
    end if
    read(UNIT,'(A)',iostat=Reason) buffer
    if(Reason.ne.0)then
       write(0,'("ERROR: The file is not in xyz format.")')
       iostat_ = 1
       if(present(iostat)) iostat = iostat_
       return
    end if
    index1 = index(buffer,'Lattice="') + 9
    index2 = index(buffer(index1:),'"') + index1 - 2
    read(buffer(index1:index2),*) ( ( basis%lat(i,j), j = 1, 3), i = 1, 3)
 
    index1 = index(buffer,'free_energy=') + 12
    read(buffer(index1:),*) basis%energy
 
 
    !---------------------------------------------------------------------------
    ! read basis
    !---------------------------------------------------------------------------
    allocate(tmp_spec(basis%natom))
    allocate(tmp_num(basis%natom))
    allocate(tmp_bas(basis%natom,basis%natom,length_))
    tmp_num(:) = 0
    tmp_spec = ""
    tmp_bas = 0
    basis%nspec = 0
    do i = 1, basis%natom
       read(UNIT,*,iostat=Reason) ctmp, vec(1:3)
       if(.not.any(tmp_spec(1:basis%nspec).eq.ctmp))then
          basis%nspec=basis%nspec+1
          tmp_spec(basis%nspec) = trim(adjustl(ctmp))
          tmp_bas(basis%nspec,1,1:3) = vec(1:3)
          tmp_num(basis%nspec) = 1
       else
          checkspec: do j = 1, basis%nspec
             if(tmp_spec(j).eq.ctmp)then
                tmp_num(j) = tmp_num(j) + 1
                tmp_bas(j,tmp_num(j),1:3) = vec(1:3)
                exit checkspec
             end if
          end do checkspec
       end if
    end do
 
 
    !---------------------------------------------------------------------------
    ! move basis from temporary basis to main basis.
    ! done to allow for correct allocation of number of and per species
    !---------------------------------------------------------------------------
    allocate(basis%spec(basis%nspec))
    basis%sysname = ""
    do i = 1, basis%nspec
       basis%spec(i)%name = tmp_spec(i)
       basis%spec(i)%num = tmp_num(i)
       allocate(basis%spec(i)%atom(tmp_num(i),length_))
       basis%spec(i)%atom(:,:) = 0
       basis%spec(i)%atom(1:tmp_num(i),1:3) = tmp_bas(i,1:tmp_num(i),1:3)
       write(buffer,'(I0,A)') basis%spec(i)%num,trim(basis%spec(i)%name)
       basis%sysname = basis%sysname//trim(buffer)
       if(i.lt.basis%nspec) basis%sysname = trim(adjustl(basis%sysname))//"_"
    end do
 
    if(present(iostat)) iostat = iostat_
 
  end subroutine extXYZ_geom_read
!###############################################################################


!###############################################################################
  subroutine extXYZ_geom_write(UNIT, basis)
    !! Write the structure in extended xyz style format.
    implicit none

    ! Arguments
    integer, intent(in) :: UNIT
    !! The unit number of the file to write to.
    class(basis_type), intent(in) :: basis
    !! The basis to write the geometry from.

    ! Local variables
    integer :: i, j
    !! Loop index.


    write(UNIT,'(I0)') basis%natom
    write(UNIT,'(A,8(F0.8,1X),F0.8,A)', advance="no") &
         'Lattice="',((basis%lat(i,j),j=1,3),i=1,3),'"'
    write(UNIT,'(A,F0.8)', advance="no") ' free_energy=',basis%energy
    write(UNIT,'(A)', advance="no") ' pbc="T T T"'
    if(basis%lcart)then
       do i = 1, basis%nspec
          do j = 1, basis%spec(i)%num
             write(UNIT,'(A8,3(1X, F16.8))') &
                  basis%spec(i)%name,basis%spec(i)%atom(j,1:3)
          end do
       end do
    else
       do i = 1, basis%nspec
          do j = 1, basis%spec(i)%num
             write(UNIT,'(A8,3(1X, F16.8))') basis%spec(i)%name, &
                  matmul(basis%spec(i)%atom(j,1:3),basis%lat)
          end do
       end do
    end if

  end subroutine extXYZ_geom_write
!###############################################################################


!###############################################################################
  subroutine convert(this)
    !! Convert the basis between direct and cartesian coordinates.
    implicit none
   
    ! Arguments
    class(basis_type), intent(inout) :: this
    !! Parent. The basis to convert.

    ! Local variables
    integer :: is, ia
    !! Loop index.
    real(real32), dimension(3,3) :: lattice
    !! The reciprocal lattice vectors.

   
    if(this%lcart)then
       lattice = inverse_3x3( this%lat )
    else
       lattice = this%lat
    end if

    this%lcart = .not.this%lcart
    do is = 1, this%nspec
       do ia = 1, this%spec(is)%num
          this%spec(is)%atom(ia,1:3) = &
               matmul( this%spec(is)%atom(ia,1:3), lattice )
       end do
    end do
   
  end subroutine convert
!###############################################################################


!#############################################################################
  function convert_abc_to_lat(abc,angle,radians) result(lattice)
    !! Convert the lattice from abc and αβγ to lattice matrix.
    implicit none

    ! Arguments
    real(real32), dimension(3), intent(in) :: abc, angle
    !! lattice constants
    logical, intent(in), optional :: radians
    !! Optional. Boolean whether angles are in radians.
    real(real32), dimension(3,3) :: lattice
    !! The lattice matrix.

    ! Local variables
    real(real32), dimension(3) :: in_angle
    !! The lattice angles in radians.



    in_angle = angle
    if(present(radians))then
       if(.not.radians) in_angle = angle*pi/180._real32
    end if

    lattice=0._real32

    lattice(1,1)=abc(1)
    lattice(2,:2)=(/abc(2)*cos(in_angle(3)),abc(2)*sin(in_angle(3))/)

    lattice(3,1) = abc(3)*cos(in_angle(2))
    lattice(3,2) = abc(3)*(cos(in_angle(1)) - cos(in_angle(2))*&
         cos(in_angle(3)))/sin(in_angle(3))
    lattice(3,3) = sqrt(abc(3)**2._real32 - &
         lattice(3,1)**2._real32 - &
         lattice(3,2)**2._real32)

  end function convert_abc_to_lat
!###############################################################################


!###############################################################################
  function convert_lat_to_abc(lattice, radians) result(abc_angle)
    !! Convert the lattice from lattice matrix to abc and αβγ.
    implicit none

    ! Arguments
    real(real32), dimension(3,3), intent(in) :: lattice
    !! The lattice matrix.
    logical, intent(in), optional :: radians
    !! Optional. Boolean whether to return angles in radians.
    real(real32), dimension(2,3) :: abc_angle
    !! The lattice constants and angles.

    ! Local variables
    integer :: i
    !! Loop index.


    do i = 1, 3
       abc_angle(1,i)=modu(lattice(i,:))
    end do
    do i = 1, 3
    end do
    abc_angle(2,1)=acos(dot_product(lattice(2,:),lattice(3,:))/&
         (abc_angle(1,2)*abc_angle(1,3)))
    abc_angle(2,3)=acos(dot_product(lattice(1,:),lattice(3,:))/&
         (abc_angle(1,1)*abc_angle(1,3)))
    abc_angle(2,3)=acos(dot_product(lattice(1,:),lattice(2,:))/&
         (abc_angle(1,1)*abc_angle(1,2)))

    if(present(radians))then
       if(.not.radians) abc_angle(2,:)=abc_angle(2,:)*180._real32/pi
    end if

  end function convert_lat_to_abc
!###############################################################################


!###############################################################################
  function get_lattice_constants(this, radians) result(output)
    !! Convert the lattice from lattice matrix to abc and αβγ.
    implicit none

    ! Arguments
    class(basis_type), intent(in) :: this
    !! Parent. The basis.
    logical, intent(in), optional :: radians
    !! Optional. Boolean whether to return angles in radians.
    real(real32), dimension(2,3) :: output
    !! The lattice constants and angles.

    ! Local variables
    logical :: radians_
    !! Boolean whether to return angles in radians.


    radians_ = .true.
    if(present(radians)) radians_ = radians

    output = convert_lat_to_abc(this%lat, radians_)

  end function get_lattice_constants
!###############################################################################


!###############################################################################
  subroutine copy(this, basis, length)
    !! Copy the basis.
    implicit none

    ! Arguments
    class(basis_type), intent(inout) :: this
    !! Parent. The basis to copy into.
    class(basis_type), intent(in) :: basis
    !! The basis to copy from.
    integer, intent(in), optional :: length
    !! The dimension of the basis atom positions.


    ! Local variables
    integer :: i
    !! Loop index.
    integer :: length_, length_input
    !! The dimension of the basis atom positions.


    !---------------------------------------------------------------------------
    ! determines whether user wants output basis extra translational dimension
    !---------------------------------------------------------------------------
    length_input = size(basis%spec(1)%atom(1,:),dim=1)
    if(present(length))then
       length_ = length
    else
       length_ = length_input
    end if


    !---------------------------------------------------------------------------
    ! if already allocated, deallocates output basis
    !---------------------------------------------------------------------------
    if(allocated(this%spec))then
       do i = 1, this%nspec
          if(allocated(this%spec(i)%atom)) deallocate(this%spec(i)%atom)
       end do
       deallocate(this%spec)
    end if


    !---------------------------------------------------------------------------
    ! allocates output basis and clones data from input basis to output basis
    !---------------------------------------------------------------------------
    allocate(this%spec(basis%nspec))
    do i = 1, basis%nspec
       allocate(this%spec(i)%atom(&
            basis%spec(i)%num,length_))

       this%spec(i)%atom(:,:3) = basis%spec(i)%atom(:,:3)
       if(length_input.eq.length_)then
          this%spec(i)%atom(:,:length_) = basis%spec(i)%atom(:,:length_)
       elseif(length_input.gt.length_)then
          this%spec(i)%atom(:,:3) = basis%spec(i)%atom(:,:3)
       elseif(length_input.lt.length_)then
          this%spec(i)%atom(:,:3) = basis%spec(i)%atom(:,:3)
          this%spec(i)%atom(:,4) = 1._real32
       end if
       this%spec(i)%num = basis%spec(i)%num
       this%spec(i)%name = strip_null(basis%spec(i)%name)

       this%spec(i)%mass = basis%spec(i)%mass
       this%spec(i)%charge = basis%spec(i)%charge
       this%spec(i)%radius = basis%spec(i)%radius
    end do
    this%nspec = basis%nspec
    this%natom = basis%natom
    this%lcart = basis%lcart
    this%sysname = basis%sysname
    this%energy = basis%energy
    this%lat = basis%lat
    this%pbc = basis%pbc

  end subroutine copy
!###############################################################################


!###############################################################################
  subroutine remove_atom(this, ispec, iatom)
    !! Remove an atom from the basis.
    implicit none

    ! Arguments
    class(basis_type), intent(inout) :: this
    !! Parent. The basis.
    integer, intent(in) :: ispec, iatom
    !! The species and atom to remove.

    ! Local variables
    integer :: i
    !! Loop index.
    real(real32), dimension(:,:), allocatable :: atom
    !! Temporary array to store the atomic positions.


    !---------------------------------------------------------------------------
    ! remove atom from basis
    !---------------------------------------------------------------------------
    do i = 1, this%nspec
       if(i.eq.ispec)then
          if(iatom.gt.this%spec(i)%num)then
             call stop_program("Atom to remove does not exist")
             return
          end if
          allocate(atom(this%spec(i)%num-1,size(this%spec(i)%atom,2)))
          atom(1:iatom-1:1,:) = this%spec(i)%atom(1:iatom-1:1,:)
          atom(iatom:this%spec(i)%num-1:1,:) = &
               this%spec(i)%atom(iatom+1:this%spec(i)%num:1,:)
          this%spec(i)%atom = atom
          deallocate(atom)
          this%spec(i)%num = this%spec(i)%num - 1
          this%natom = this%natom - 1
          if(this%spec(i)%num.eq.0)then
             deallocate(this%spec(i)%atom)
             if(this%nspec.eq.0)then
                deallocate(this%spec)
                this%lcart = .true.
                this%sysname = ""
                this%energy = 0._real32
                this%lat = 0._real32
                this%pbc = .true.
             end if
          end if
       end if
    end do

  end subroutine remove_atom
!###############################################################################


!###############################################################################
  subroutine remove_atoms(this, atoms)
    !! Remove atoms from the basis.
    use raffle__misc, only: swap
    implicit none

    ! Arguments
    class(basis_type), intent(inout) :: this
    !! Parent. The basis.
    integer, dimension(:,:), intent(in) :: atoms
    !! The atoms to remove (2, number of atoms to remove)
    !! 1st value of 1st dimension is the species number
    !! 2nd value of 1st dimension is the atom number
    !! 2nd dimension is the number of atoms to remove

    ! Local variables
    integer :: is, ia, i
    !! Loop index.
    integer :: n, m, start_idx, end_idx, loc
    !! Index variables.
    integer :: num_species
    !! The number of species.
    integer, dimension(:,:), allocatable :: atoms_ordered
    !! The atoms to remove ordered by species and atom
    real(real32), dimension(:,:), allocatable :: atom
    !! Temporary array to store the atomic positions.


    !---------------------------------------------------------------------------
    ! reorder atoms to remove
    !---------------------------------------------------------------------------
    allocate(atoms_ordered, source=atoms)
    n = size(atoms_ordered, 1)
    m = size(atoms_ordered, 2)

    do i = 1, m
       loc = maxloc(atoms_ordered(1, i:n), dim=1) + i - 1
       if (loc .ne. i) then
          call swap(atoms_ordered(1, i), atoms_ordered(1, loc))
          call swap(atoms_ordered(2, i), atoms_ordered(2, loc))
       end if
    end do
    num_species = this%nspec
    do is = 1, num_species
       start_idx = findloc(atoms_ordered(1, :), is, dim=1)
       end_idx   = findloc(atoms_ordered(1, :), is, dim=1, back=.true.)
       if (start_idx .eq. 0) cycle
       do ia = start_idx, end_idx, 1
          loc = maxloc( &
               atoms_ordered(2, ia:end_idx), &
               dim=1 &
          ) + ia - 1
          if (loc .ne. ia) then
             call swap(atoms_ordered(1, ia), atoms_ordered(1, loc))
             call swap(atoms_ordered(2, ia), atoms_ordered(2, loc))
          end if
       end do
    end do


    !---------------------------------------------------------------------------
    ! remove atoms from basis
    !---------------------------------------------------------------------------
    do i = 1, size(atoms_ordered, 2)
       call this%remove_atom(atoms_ordered(1, i), atoms_ordered(2, i))
    end do

    do is = 1, this%nspec
       if (this%spec(is)%num .eq. 0) then
          this%spec = [ this%spec(1:is-1), this%spec(is+1:) ]
          this%nspec = this%nspec - 1
       end if
    end do

  end subroutine remove_atoms
!###############################################################################


!###############################################################################
  subroutine get_element_properties(element, charge, mass, radius)
    !! Set the mass and charge of the element
    implicit none

    ! Arguments
    character(len=3), intent(in) :: element
    !! Element name.
    real(real32), intent(out), optional :: charge
    !! Charge of the element.
    real(real32), intent(out), optional :: mass
    !! Mass of the element.
    real(real32), intent(out), optional :: radius
    !! Radius of the element.

    ! Local variables
    real(real32) :: mass_, charge_, radius_
    !! Mass, charge and radius of the element.

    select case(element)
    case('H')
       mass_ = 1.00784_real32
       charge_ = 1.0_real32
       radius_ = 0.31_real32
    case('He')
       mass_ = 4.0026_real32
       charge_ = 2.0_real32
       radius_ = 0.28_real32
    case('Li')
       mass_ = 6.94_real32
       charge_ = 3.0_real32
       radius_ = 1.28_real32
    case('Be')
       mass_ = 9.0122_real32
       charge_ = 4.0_real32
       radius_ = 0.96_real32
    case('B')
       mass_ = 10.81_real32
       charge_ = 5.0_real32
       radius_ = 0.84_real32
    case('C')
       mass_ = 12.011_real32
       charge_ = 6.0_real32
       radius_ = 0.76_real32
    case('N')
       mass_ = 14.007_real32
       charge_ = 7.0_real32
       radius_ = 0.71_real32
    case('O')
       mass_ = 15.999_real32
       charge_ = 8.0_real32
       radius_ = 0.66_real32
    case('F')
       mass_ = 18.998_real32
       charge_ = 9.0_real32
       radius_ = 0.57_real32
    case('Ne')
       mass_ = 20.180_real32
       charge_ = 10.0_real32
       radius_ = 0.58_real32
    case('Na')
       mass_ = 22.989_real32
       charge_ = 11.0_real32
       radius_ = 1.66_real32
    case('Mg')
       mass_ = 24.305_real32
       charge_ = 12.0_real32
       radius_ = 1.41_real32
    case('Al')
       mass_ = 26.982_real32
       charge_ = 13.0_real32
       radius_ = 1.21_real32
    case('Si')
       mass_ = 28.085_real32
       charge_ = 14.0_real32
       radius_ = 1.11_real32
    case('P')
       mass_ = 30.974_real32
       charge_ = 15.0_real32
       radius_ = 1.07_real32
    case('S')  
       mass_ = 32.06_real32
       charge_ = 16.0_real32
       radius_ = 1.05_real32
    case('Cl')
       mass_ = 35.453_real32
       charge_ = 17.0_real32
       radius_ = 1.02_real32
    case('Ar')
       mass_ = 39.948_real32
       charge_ = 18.0_real32
       radius_ = 1.06_real32
    case('K')
       mass_ = 39.098_real32
       charge_ = 19.0_real32
       radius_ = 2.03_real32
    case('Ca')
       mass_ = 40.078_real32
       charge_ = 20.0_real32
       radius_ = 1.74_real32
    case('Sc')
       mass_ = 44.956_real32
       charge_ = 21.0_real32
       radius_ = 1.44_real32
    case('Ti')
       mass_ = 47.867_real32
       charge_ = 22.0_real32
       radius_ = 1.32_real32
    case('V')
       mass_ = 50.942_real32
       charge_ = 23.0_real32
       radius_ = 1.22_real32
    case('Cr')
       mass_ = 51.996_real32
       charge_ = 24.0_real32
       radius_ = 1.18_real32
    case('Mn')
       mass_ = 54.938_real32
       charge_ = 25.0_real32
       radius_ = 1.17_real32
    case('Fe')
       mass_ = 55.845_real32
       charge_ = 26.0_real32
       radius_ = 1.17_real32
    case('Co')
       mass_ = 58.933_real32
       charge_ = 27.0_real32
       radius_ = 1.16_real32
    case('Ni')
       mass_ = 58.693_real32
       charge_ = 28.0_real32
       radius_ = 1.15_real32
    case('Cu')
       mass_ = 63.546_real32
       charge_ = 29.0_real32
       radius_ = 1.17_real32
    case('Zn')
       mass_ = 65.38_real32
       charge_ = 30.0_real32
       radius_ = 1.25_real32
    case('Ga')
       mass_ = 69.723_real32
       charge_ = 31.0_real32
       radius_ = 1.26_real32
    case('Ge')
       mass_ = 72.63_real32
       charge_ = 32.0_real32
       radius_ = 1.22_real32
    case('As')
       mass_ = 74.922_real32
       charge_ = 33.0_real32
       radius_ = 1.19_real32
    case('Se')
       mass_ = 78.971_real32
       charge_ = 34.0_real32
       radius_ = 1.16_real32
    case('Br')
       mass_ = 79.904_real32
       charge_ = 35.0_real32
       radius_ = 1.14_real32
    case('Kr')
       mass_ = 83.798_real32
       charge_ = 36.0_real32
       radius_ = 1.12_real32
    case('Rb')
       mass_ = 85.468_real32
       charge_ = 37.0_real32
       radius_ = 2.16_real32
    case('Sr')
       mass_ = 87.62_real32
       charge_ = 38.0_real32
       radius_ = 1.91_real32
    case('Y')
       mass_ = 88.906_real32
       charge_ = 39.0_real32
       radius_ = 1.62_real32
    case('Zr')
       mass_ = 91.224_real32
       charge_ = 40.0_real32
       radius_ = 1.45_real32
    case('Nb')
       mass_ = 92.906_real32
       charge_ = 41.0_real32
       radius_ = 1.34_real32
    case('Mo')
       mass_ = 95.95_real32
       charge_ = 42.0_real32
       radius_ = 1.3_real32
    case('Tc')
       mass_ = 98.0_real32
       charge_ = 43.0_real32
       radius_ = 1.27_real32
    case('Ru')
       mass_ = 101.07_real32
       charge_ = 44.0_real32
       radius_ = 1.25_real32
    case('Rh')
       mass_ = 102.91_real32
       charge_ = 45.0_real32
       radius_ = 1.25_real32
    case('Pd')
       mass_ = 106.42_real32
       charge_ = 46.0_real32
       radius_ = 1.28_real32
    case('Ag')
       mass_ = 107.87_real32
       charge_ = 47.0_real32
       radius_ = 1.34_real32
    case('Cd')
       mass_ = 112.41_real32
       charge_ = 48.0_real32
       radius_ = 1.48_real32
    case('In')
       mass_ = 114.82_real32
       charge_ = 49.0_real32
       radius_ = 1.44_real32
    case('Sn')
       mass_ = 118.71_real32
       charge_ = 50.0_real32
       radius_ = 1.41_real32
    case('Sb')
       mass_ = 121.76_real32
       charge_ = 51.0_real32
       radius_ = 1.38_real32
    case('Te')
       mass_ = 127.6_real32
       charge_ = 52.0_real32
       radius_ = 1.35_real32
    case('I')
       mass_ = 126.9_real32
       charge_ = 53.0_real32
       radius_ = 1.33_real32
    case('Xe')
       mass_ = 131.29_real32
       charge_ = 54.0_real32
       radius_ = 1.31_real32
    case('Cs')
       mass_ = 132.91_real32
       charge_ = 55.0_real32
       radius_ = 2.35_real32
    case('Ba')
       mass_ = 137.33_real32
       charge_ = 56.0_real32
       radius_ = 1.98_real32
    case('La')
       mass_ = 138.91_real32
       charge_ = 57.0_real32
       radius_ = 1.69_real32
    case('Ce')
       mass_ = 140.12_real32
       charge_ = 58.0_real32
       radius_ = 1.65_real32
    case('Pr')
       mass_ = 140.91_real32
       charge_ = 59.0_real32
       radius_ = 1.65_real32
    case('Nd')
       mass_ = 144.24_real32
       charge_ = 60.0_real32
       radius_ = 1.64_real32
    case('Pm')
       mass_ = 145.0_real32
       charge_ = 61.0_real32
       radius_ = 1.63_real32
    case('Sm')
       mass_ = 150.36_real32
       charge_ = 62.0_real32
       radius_ = 1.62_real32
    case('Eu')
       mass_ = 152.0_real32
       charge_ = 63.0_real32
       radius_ = 1.85_real32
    case('Gd')
       mass_ = 157.25_real32
       charge_ = 64.0_real32
       radius_ = 1.61_real32
    case('Tb')
       mass_ = 158.93_real32
       charge_ = 65.0_real32
       radius_ = 1.59_real32
    case('Dy')
       mass_ = 162.5_real32
       charge_ = 66.0_real32
       radius_ = 1.59_real32
    case('Ho')
       mass_ = 164.93_real32
       charge_ = 67.0_real32
       radius_ = 1.58_real32
    case('Er')
       mass_ = 167.26_real32
       charge_ = 68.0_real32
       radius_ = 1.57_real32
    case('Tm')
       mass_ = 168.93_real32
       charge_ = 69.0_real32
       radius_ = 1.56_real32
    case('Yb')
       mass_ = 173.05_real32
       charge_ = 70.0_real32
       radius_ = 1.74_real32
    case('Lu')
       mass_ = 174.97_real32
       charge_ = 71.0_real32
       radius_ = 1.56_real32
    case('Hf')
       mass_ = 178.49_real32
       charge_ = 72.0_real32
       radius_ = 1.44_real32
    case('Ta')
       mass_ = 180.95_real32
       charge_ = 73.0_real32
       radius_ = 1.34_real32
    case('W')
       mass_ = 183.84_real32
       charge_ = 74.0_real32
       radius_ = 1.3_real32
    case('Re')
       mass_ = 186.21_real32
       charge_ = 75.0_real32
       radius_ = 1.28_real32
    case('Os')
       mass_ = 190.23_real32
       charge_ = 76.0_real32
       radius_ = 1.26_real32
    case('Ir')
       mass_ = 192.22_real32
       charge_ = 77.0_real32
       radius_ = 1.27_real32
    case('Pt')
       mass_ = 195.08_real32
       charge_ = 78.0_real32
       radius_ = 1.3_real32
    case('Au')
       mass_ = 196.97_real32
       charge_ = 79.0_real32
       radius_ = 1.34_real32
    case('Hg')
       mass_ = 200.59_real32
       charge_ = 80.0_real32
       radius_ = 1.49_real32
    case('Tl')
       mass_ = 204.38_real32
       charge_ = 81.0_real32
       radius_ = 1.48_real32
    case('Pb')
       mass_ = 207.2_real32
       charge_ = 82.0_real32
       radius_ = 1.47_real32
    case('Bi')
       mass_ = 208.98_real32
       charge_ = 83.0_real32
       radius_ = 1.46_real32
    case('Po')
       mass_ = 209.0_real32
       charge_ = 84.0_real32
       radius_ = 1.45_real32
    case('At')
       mass_ = 210.0_real32
       charge_ = 85.0_real32
       radius_ = 1.44_real32
    case('Rn')
       mass_ = 222.0_real32
       charge_ = 86.0_real32
       radius_ = 1.43_real32
    case('Fr')
       mass_ = 223.0_real32
       charge_ = 87.0_real32
       radius_ = 2.6_real32
    case('Ra')
       mass_ = 226.0_real32
       charge_ = 88.0_real32
       radius_ = 2.21_real32
    case('Ac')
       mass_ = 227.0_real32
       charge_ = 89.0_real32
       radius_ = 1.86_real32
    case('Th')
       mass_ = 232.04_real32
       charge_ = 90.0_real32
       radius_ = 1.75_real32
    case('Pa')
       mass_ = 231.04_real32
       charge_ = 91.0_real32
       radius_ = 1.61_real32
    case('U')
       mass_ = 238.03_real32
       charge_ = 92.0_real32
       radius_ = 1.58_real32
    case('Np')
       mass_ = 237.0_real32
       charge_ = 93.0_real32
       radius_ = 1.55_real32
    case('Pu')
       mass_ = 244.0_real32
       charge_ = 94.0_real32
       radius_ = 1.53_real32
    case('Am')
       mass_ = 243.0_real32
       charge_ = 95.0_real32
       radius_ = 1.51_real32
    case('Cm')
       mass_ = 247.0_real32
       charge_ = 96.0_real32
       radius_ = 1.69_real32
    case('Bk')
       mass_ = 247.0_real32
       charge_ = 97.0_real32
       radius_ = 1.48_real32
    case('Cf')
       mass_ = 251.0_real32
       charge_ = 98.0_real32
       radius_ = 1.47_real32
    case('Es')
       mass_ = 252.0_real32
       charge_ = 99.0_real32
       radius_ = 1.46_real32
    case('Fm')
       mass_ = 257.0_real32
       charge_ = 100.0_real32
       radius_ = 1.45_real32
    case('Md')
       mass_ = 258.0_real32
       charge_ = 101.0_real32
       radius_ = 1.44_real32
    case('No')
       mass_ = 259.0_real32
       charge_ = 102.0_real32
       radius_ = 1.43_real32
    case('Lr')
       mass_ = 262.0_real32
       charge_ = 103.0_real32
       radius_ = 1.62_real32
    case('Rf')
       mass_ = 267.0_real32
       charge_ = 104.0_real32
       radius_ = 1.57_real32
    case('Db')
       mass_ = 270.0_real32
       charge_ = 105.0_real32
       radius_ = 1.49_real32
    case('Sg')
       mass_ = 271.0_real32
       charge_ = 106.0_real32
       radius_ = 1.43_real32
    case('Bh')
       mass_ = 270.0_real32
       charge_ = 107.0_real32
       radius_ = 1.41_real32
    case('Hs')
       mass_ = 277.0_real32
       charge_ = 108.0_real32
       radius_ = 1.34_real32
    case('Mt')
       mass_ = 276.0_real32
       charge_ = 109.0_real32
       radius_ = 1.29_real32
    case('Ds')
       mass_ = 281.0_real32
       charge_ = 110.0_real32
       radius_ = 1.28_real32
    case('Rg')
       mass_ = 280.0_real32
       charge_ = 111.0_real32
       radius_ = 1.21_real32
    case('Cn')
       mass_ = 285.0_real32
       charge_ = 112.0_real32
       radius_ = 1.22_real32
    case('Nh')
       mass_ = 284.0_real32
       charge_ = 113.0_real32
       radius_ = 1.21_real32
    case('Fl')
       mass_ = 289.0_real32
       charge_ = 114.0_real32
       radius_ = 1.21_real32
    case('Mc')
       mass_ = 288.0_real32
       charge_ = 115.0_real32
       radius_ = 1.21_real32
    case('Lv')
       mass_ = 293.0_real32
       charge_ = 116.0_real32
       radius_ = 1.21_real32
    case('Ts')
       mass_ = 294.0_real32
       charge_ = 117.0_real32
       radius_ = 1.21_real32
    case('Og')
       mass_ = 294.0_real32
       charge_ = 118.0_real32
       radius_ = 1.21_real32
    case default
       ! handle unknown element
       mass_ = 0.0_real32
       charge_ = 0.0_real32
       radius_ = 0.0_real32
    end select

    !---------------------------------------------------------------------------
    ! Return the values
    !---------------------------------------------------------------------------
    if(present(mass)) mass = mass_
    if(present(charge)) charge = charge_
    if(present(radius)) radius = radius_

  end subroutine get_element_properties
!###############################################################################

end module raffle__geom_rw