! Module rw_geom defined in file ../src/lib/mod_rw_geom.f90

subroutine f90wrap_spec_type__array__atom(this, nd, dtype, dshape, dloc)
    use rw_geom, only: spec_type
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type spec_type_ptr_type
        type(spec_type), pointer :: p => NULL()
    end type spec_type_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(spec_type_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 11
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%atom)) then
        dshape(1:2) = shape(this_ptr%p%atom)
        dloc = loc(this_ptr%p%atom)
    else
        dloc = 0
    end if
end subroutine f90wrap_spec_type__array__atom

subroutine f90wrap_spec_type__get__mass(this, f90wrap_mass)
    use rw_geom, only: spec_type
    implicit none
    type spec_type_ptr_type
        type(spec_type), pointer :: p => NULL()
    end type spec_type_ptr_type
    integer, intent(in)   :: this(2)
    type(spec_type_ptr_type) :: this_ptr
    real(4), intent(out) :: f90wrap_mass
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_mass = this_ptr%p%mass
end subroutine f90wrap_spec_type__get__mass

subroutine f90wrap_spec_type__set__mass(this, f90wrap_mass)
    use rw_geom, only: spec_type
    implicit none
    type spec_type_ptr_type
        type(spec_type), pointer :: p => NULL()
    end type spec_type_ptr_type
    integer, intent(in)   :: this(2)
    type(spec_type_ptr_type) :: this_ptr
    real(4), intent(in) :: f90wrap_mass
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%mass = f90wrap_mass
end subroutine f90wrap_spec_type__set__mass

subroutine f90wrap_spec_type__get__charge(this, f90wrap_charge)
    use rw_geom, only: spec_type
    implicit none
    type spec_type_ptr_type
        type(spec_type), pointer :: p => NULL()
    end type spec_type_ptr_type
    integer, intent(in)   :: this(2)
    type(spec_type_ptr_type) :: this_ptr
    real(4), intent(out) :: f90wrap_charge
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_charge = this_ptr%p%charge
end subroutine f90wrap_spec_type__get__charge

subroutine f90wrap_spec_type__set__charge(this, f90wrap_charge)
    use rw_geom, only: spec_type
    implicit none
    type spec_type_ptr_type
        type(spec_type), pointer :: p => NULL()
    end type spec_type_ptr_type
    integer, intent(in)   :: this(2)
    type(spec_type_ptr_type) :: this_ptr
    real(4), intent(in) :: f90wrap_charge
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%charge = f90wrap_charge
end subroutine f90wrap_spec_type__set__charge

subroutine f90wrap_spec_type__get__name(this, f90wrap_name)
    use rw_geom, only: spec_type
    implicit none
    type spec_type_ptr_type
        type(spec_type), pointer :: p => NULL()
    end type spec_type_ptr_type
    integer, intent(in)   :: this(2)
    type(spec_type_ptr_type) :: this_ptr
    character(3), intent(out) :: f90wrap_name
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_name = this_ptr%p%name
end subroutine f90wrap_spec_type__get__name

subroutine f90wrap_spec_type__set__name(this, f90wrap_name)
    use rw_geom, only: spec_type
    implicit none
    type spec_type_ptr_type
        type(spec_type), pointer :: p => NULL()
    end type spec_type_ptr_type
    integer, intent(in)   :: this(2)
    type(spec_type_ptr_type) :: this_ptr
    character(3), intent(in) :: f90wrap_name
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%name = f90wrap_name
end subroutine f90wrap_spec_type__set__name

subroutine f90wrap_spec_type__get__num(this, f90wrap_num)
    use rw_geom, only: spec_type
    implicit none
    type spec_type_ptr_type
        type(spec_type), pointer :: p => NULL()
    end type spec_type_ptr_type
    integer, intent(in)   :: this(2)
    type(spec_type_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_num
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_num = this_ptr%p%num
end subroutine f90wrap_spec_type__get__num

subroutine f90wrap_spec_type__set__num(this, f90wrap_num)
    use rw_geom, only: spec_type
    implicit none
    type spec_type_ptr_type
        type(spec_type), pointer :: p => NULL()
    end type spec_type_ptr_type
    integer, intent(in)   :: this(2)
    type(spec_type_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_num
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%num = f90wrap_num
end subroutine f90wrap_spec_type__set__num

subroutine f90wrap_spec_type__array__lat(this, nd, dtype, dshape, dloc)
    use rw_geom, only: spec_type
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type spec_type_ptr_type
        type(spec_type), pointer :: p => NULL()
    end type spec_type_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(spec_type_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 11
    this_ptr = transfer(this, this_ptr)
    dshape(1:2) = shape(this_ptr%p%lat)
    dloc = loc(this_ptr%p%lat)
end subroutine f90wrap_spec_type__array__lat

subroutine f90wrap_rw_geom__spec_type_initialise(this)
    use rw_geom, only: spec_type
    implicit none
    
    type spec_type_ptr_type
        type(spec_type), pointer :: p => NULL()
    end type spec_type_ptr_type
    type(spec_type_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_rw_geom__spec_type_initialise

subroutine f90wrap_rw_geom__spec_type_finalise(this)
    use rw_geom, only: spec_type
    implicit none
    
    type spec_type_ptr_type
        type(spec_type), pointer :: p => NULL()
    end type spec_type_ptr_type
    type(spec_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_rw_geom__spec_type_finalise

subroutine f90wrap_bas_type__array_getitem__spec(f90wrap_this, f90wrap_i, specitem)
    
    use rw_geom, only: bas_type, spec_type
    implicit none
    
    type bas_type_ptr_type
        type(bas_type), pointer :: p => NULL()
    end type bas_type_ptr_type
    type spec_type_ptr_type
        type(spec_type), pointer :: p => NULL()
    end type spec_type_ptr_type
    integer, intent(in) :: f90wrap_this(2)
    type(bas_type_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_i
    integer, intent(out) :: specitem(2)
    type(spec_type_ptr_type) :: spec_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%spec)) then
        if (f90wrap_i < 1 .or. f90wrap_i > size(this_ptr%p%spec)) then
            call f90wrap_abort("array index out of range")
        else
            spec_ptr%p => this_ptr%p%spec(f90wrap_i)
            specitem = transfer(spec_ptr,specitem)
        endif
    else
        call f90wrap_abort("derived type array not allocated")
    end if
end subroutine f90wrap_bas_type__array_getitem__spec

subroutine f90wrap_bas_type__array_setitem__spec(f90wrap_this, f90wrap_i, specitem)
    
    use rw_geom, only: bas_type, spec_type
    implicit none
    
    type bas_type_ptr_type
        type(bas_type), pointer :: p => NULL()
    end type bas_type_ptr_type
    type spec_type_ptr_type
        type(spec_type), pointer :: p => NULL()
    end type spec_type_ptr_type
    integer, intent(in) :: f90wrap_this(2)
    type(bas_type_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_i
    integer, intent(in) :: specitem(2)
    type(spec_type_ptr_type) :: spec_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%spec)) then
        if (f90wrap_i < 1 .or. f90wrap_i > size(this_ptr%p%spec)) then
            call f90wrap_abort("array index out of range")
        else
            spec_ptr = transfer(specitem,spec_ptr)
            this_ptr%p%spec(f90wrap_i) = spec_ptr%p
        endif
    else
        call f90wrap_abort("derived type array not allocated")
    end if
end subroutine f90wrap_bas_type__array_setitem__spec

subroutine f90wrap_bas_type__array_len__spec(f90wrap_this, f90wrap_n)
    
    use rw_geom, only: bas_type, spec_type
    implicit none
    
    type bas_type_ptr_type
        type(bas_type), pointer :: p => NULL()
    end type bas_type_ptr_type
    type spec_type_ptr_type
        type(spec_type), pointer :: p => NULL()
    end type spec_type_ptr_type
    integer, intent(out) :: f90wrap_n
    integer, intent(in) :: f90wrap_this(2)
    type(bas_type_ptr_type) :: this_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%spec)) then
        f90wrap_n = size(this_ptr%p%spec)
    else
        f90wrap_n = 0
    end if
end subroutine f90wrap_bas_type__array_len__spec

subroutine f90wrap_bas_type__get__nspec(this, f90wrap_nspec)
    use rw_geom, only: bas_type
    implicit none
    type bas_type_ptr_type
        type(bas_type), pointer :: p => NULL()
    end type bas_type_ptr_type
    integer, intent(in)   :: this(2)
    type(bas_type_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_nspec
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_nspec = this_ptr%p%nspec
end subroutine f90wrap_bas_type__get__nspec

subroutine f90wrap_bas_type__set__nspec(this, f90wrap_nspec)
    use rw_geom, only: bas_type
    implicit none
    type bas_type_ptr_type
        type(bas_type), pointer :: p => NULL()
    end type bas_type_ptr_type
    integer, intent(in)   :: this(2)
    type(bas_type_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_nspec
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%nspec = f90wrap_nspec
end subroutine f90wrap_bas_type__set__nspec

subroutine f90wrap_bas_type__get__natom(this, f90wrap_natom)
    use rw_geom, only: bas_type
    implicit none
    type bas_type_ptr_type
        type(bas_type), pointer :: p => NULL()
    end type bas_type_ptr_type
    integer, intent(in)   :: this(2)
    type(bas_type_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_natom
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_natom = this_ptr%p%natom
end subroutine f90wrap_bas_type__get__natom

subroutine f90wrap_bas_type__set__natom(this, f90wrap_natom)
    use rw_geom, only: bas_type
    implicit none
    type bas_type_ptr_type
        type(bas_type), pointer :: p => NULL()
    end type bas_type_ptr_type
    integer, intent(in)   :: this(2)
    type(bas_type_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_natom
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%natom = f90wrap_natom
end subroutine f90wrap_bas_type__set__natom

subroutine f90wrap_bas_type__get__energy(this, f90wrap_energy)
    use rw_geom, only: bas_type
    implicit none
    type bas_type_ptr_type
        type(bas_type), pointer :: p => NULL()
    end type bas_type_ptr_type
    integer, intent(in)   :: this(2)
    type(bas_type_ptr_type) :: this_ptr
    real(4), intent(out) :: f90wrap_energy
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_energy = this_ptr%p%energy
end subroutine f90wrap_bas_type__get__energy

subroutine f90wrap_bas_type__set__energy(this, f90wrap_energy)
    use rw_geom, only: bas_type
    implicit none
    type bas_type_ptr_type
        type(bas_type), pointer :: p => NULL()
    end type bas_type_ptr_type
    integer, intent(in)   :: this(2)
    type(bas_type_ptr_type) :: this_ptr
    real(4), intent(in) :: f90wrap_energy
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%energy = f90wrap_energy
end subroutine f90wrap_bas_type__set__energy

subroutine f90wrap_bas_type__get__lcart(this, f90wrap_lcart)
    use rw_geom, only: bas_type
    implicit none
    type bas_type_ptr_type
        type(bas_type), pointer :: p => NULL()
    end type bas_type_ptr_type
    integer, intent(in)   :: this(2)
    type(bas_type_ptr_type) :: this_ptr
    logical, intent(out) :: f90wrap_lcart
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_lcart = this_ptr%p%lcart
end subroutine f90wrap_bas_type__get__lcart

subroutine f90wrap_bas_type__set__lcart(this, f90wrap_lcart)
    use rw_geom, only: bas_type
    implicit none
    type bas_type_ptr_type
        type(bas_type), pointer :: p => NULL()
    end type bas_type_ptr_type
    integer, intent(in)   :: this(2)
    type(bas_type_ptr_type) :: this_ptr
    logical, intent(in) :: f90wrap_lcart
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%lcart = f90wrap_lcart
end subroutine f90wrap_bas_type__set__lcart

subroutine f90wrap_bas_type__get__sysname(this, f90wrap_sysname)
    use rw_geom, only: bas_type
    implicit none
    type bas_type_ptr_type
        type(bas_type), pointer :: p => NULL()
    end type bas_type_ptr_type
    integer, intent(in)   :: this(2)
    type(bas_type_ptr_type) :: this_ptr
    character(1024), intent(out) :: f90wrap_sysname
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_sysname = this_ptr%p%sysname
end subroutine f90wrap_bas_type__get__sysname

subroutine f90wrap_bas_type__set__sysname(this, f90wrap_sysname)
    use rw_geom, only: bas_type
    implicit none
    type bas_type_ptr_type
        type(bas_type), pointer :: p => NULL()
    end type bas_type_ptr_type
    integer, intent(in)   :: this(2)
    type(bas_type_ptr_type) :: this_ptr
    character(1024), intent(in) :: f90wrap_sysname
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%sysname = f90wrap_sysname
end subroutine f90wrap_bas_type__set__sysname

subroutine f90wrap_rw_geom__bas_type_initialise(this)
    use rw_geom, only: bas_type
    implicit none
    
    type bas_type_ptr_type
        type(bas_type), pointer :: p => NULL()
    end type bas_type_ptr_type
    type(bas_type_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_rw_geom__bas_type_initialise

subroutine f90wrap_rw_geom__bas_type_finalise(this)
    use rw_geom, only: bas_type
    implicit none
    
    type bas_type_ptr_type
        type(bas_type), pointer :: p => NULL()
    end type bas_type_ptr_type
    type(bas_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_rw_geom__bas_type_finalise





subroutine f90wrap_bas_type_xnum_array__array_getitem__items( &
    this, f90wrap_i, itemsitem)
    use rw_geom, only: bas_type
    implicit none

    type bas_type_xnum_array
        type(bas_type), dimension(:), allocatable :: items
    end type bas_type_xnum_array

    type bas_type_xnum_array_ptr_type
        type(bas_type_xnum_array), pointer :: p => NULL()
    end type bas_type_xnum_array_ptr_type
    type bas_type_ptr_type
        type(bas_type), pointer :: p => NULL()
    end type bas_type_ptr_type
    integer, intent(in), dimension(2) :: this
    type(bas_type_xnum_array_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_i
    integer, intent(out) :: itemsitem(2)
    type(bas_type_ptr_type) :: items_ptr
    
    this_ptr = transfer(this, this_ptr)
    if (f90wrap_i < 1 .or. f90wrap_i > size(this_ptr%p%items)) then
        call f90wrap_abort("array index out of range")
    else
        items_ptr%p => this_ptr%p%items(f90wrap_i)
        itemsitem = transfer(items_ptr,itemsitem)
    endif
end subroutine f90wrap_bas_type_xnum_array__array_getitem__items

subroutine f90wrap_bas_type_xnum_array__array_setitem__items(this, f90wrap_i, itemsitem)
    use rw_geom, only: bas_type
    implicit none

    type bas_type_xnum_array
        type(bas_type), dimension(:), allocatable :: items
    end type bas_type_xnum_array

    type bas_type_xnum_array_ptr_type
        type(bas_type_xnum_array), pointer :: p => NULL()
    end type bas_type_xnum_array_ptr_type
    type bas_type_ptr_type
        type(bas_type), pointer :: p => NULL()
    end type bas_type_ptr_type
    integer, intent(in), dimension(2) :: this
    type(bas_type_xnum_array_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_i
    integer, intent(out) :: itemsitem(2)
    type(bas_type_ptr_type) :: items_ptr
    
    this_ptr = transfer(this, this_ptr)
    if (f90wrap_i < 1 .or. f90wrap_i > size(this_ptr%p%items)) then
        call f90wrap_abort("array index out of range")
    else
        items_ptr = transfer(itemsitem,items_ptr)
        this_ptr%p%items(f90wrap_i) = items_ptr%p
    endif
end subroutine f90wrap_bas_type_xnum_array__array_setitem__items

subroutine f90wrap_bas_type_xnum_array__array_len__items(this, f90wrap_n)
    use rw_geom, only: bas_type
    implicit none

    type bas_type_xnum_array
        type(bas_type), dimension(:), allocatable :: items
    end type bas_type_xnum_array

    type bas_type_xnum_array_ptr_type
        type(bas_type_xnum_array), pointer :: p => NULL()
    end type bas_type_xnum_array_ptr_type
    integer, intent(in), dimension(2) :: this
    type(bas_type_xnum_array_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_n
    this_ptr = transfer(this, this_ptr)
    f90wrap_n = size(this_ptr%p%items)
end subroutine f90wrap_bas_type_xnum_array__array_len__items

subroutine f90wrap_bas_type_xnum_array__array_alloc__items(this, num)
    use rw_geom, only: bas_type
    implicit none

    type bas_type_xnum_array
        type(bas_type), dimension(:), allocatable :: items
    end type bas_type_xnum_array

    type bas_type_xnum_array_ptr_type
        type(bas_type_xnum_array), pointer :: p => NULL()
    end type bas_type_xnum_array_ptr_type
    type(bas_type_xnum_array_ptr_type) :: this_ptr
    integer, intent(in) :: num
    integer, intent(inout), dimension(2) :: this

    this_ptr = transfer(this, this_ptr)
    allocate(this_ptr%p%items(num))
    this = transfer(this_ptr, this)
end subroutine f90wrap_bas_type_xnum_array__array_alloc__items

subroutine f90wrap_bas_type_xnum_array__array_dealloc__items(this)
    use rw_geom, only: bas_type
    implicit none

    type bas_type_xnum_array
        type(bas_type), dimension(:), allocatable :: items
    end type bas_type_xnum_array

    type bas_type_xnum_array_ptr_type
        type(bas_type_xnum_array), pointer :: p => NULL()
    end type bas_type_xnum_array_ptr_type
    type(bas_type_xnum_array_ptr_type) :: this_ptr
    integer, intent(inout), dimension(2) :: this

    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p%items)
    this = transfer(this_ptr, this)
end subroutine f90wrap_bas_type_xnum_array__array_dealloc__items

subroutine f90wrap_rw_geom__bas_type_xnum_array_initialise(this)
    use rw_geom, only: bas_type
    implicit none

    type bas_type_xnum_array
        type(bas_type), dimension(:), allocatable :: items
    end type bas_type_xnum_array

    type bas_type_xnum_array_ptr_type
        type(bas_type_xnum_array), pointer :: p => NULL()
    end type bas_type_xnum_array_ptr_type
    type(bas_type_xnum_array_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_rw_geom__bas_type_xnum_array_initialise

subroutine f90wrap_rw_geom__bas_type_xnum_array_finalise(this)
    use rw_geom, only: bas_type
    implicit none

    type bas_type_xnum_array
        type(bas_type), dimension(:), allocatable :: items
    end type bas_type_xnum_array

    type bas_type_xnum_array_ptr_type
        type(bas_type_xnum_array), pointer :: p => NULL()
    end type bas_type_xnum_array_ptr_type
    type(bas_type_xnum_array_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_rw_geom__bas_type_xnum_array_finalise






subroutine f90wrap_rw_geom__allocate_species__binding__bas_type( &
       this, num_species, species_symbols, species_count, atoms, n0, &
       n1, n2, n3)
    use rw_geom, only: bas_type
    implicit none
    
    type bas_type_ptr_type
        type(bas_type), pointer :: p => NULL()
    end type bas_type_ptr_type
    type(bas_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    integer, intent(in), optional :: num_species
    character(3), intent(in), optional, dimension(n0) :: species_symbols
    integer, intent(in), optional, dimension(n1) :: species_count
    real(4), intent(in), optional, dimension(n2,n3) :: atoms
    integer :: n0
    !f2py intent(hide), depend(species_symbols) :: n0 = shape(species_symbols,0)
    integer :: n1
    !f2py intent(hide), depend(species_count) :: n1 = shape(species_count,0)
    integer :: n2
    !f2py intent(hide), depend(atoms) :: n2 = shape(atoms,0)
    integer :: n3
    !f2py intent(hide), depend(atoms) :: n3 = shape(atoms,1)
    this_ptr = transfer(this, this_ptr)
    call this_ptr%p%allocate_species( &
         num_species=num_species, &
         species_symbols=species_symbols, &
         species_count=species_count, &
         atoms=atoms &
    )
end subroutine f90wrap_rw_geom__allocate_species__binding__bas_type

! subroutine f90wrap_rw_geom__geom_read(unit, lat, bas, length)
!     use rw_geom, only: geom_read, bas_type
!     implicit none
    
!     type bas_type_ptr_type
!         type(bas_type), pointer :: p => NULL()
!     end type bas_type_ptr_type
!     integer :: unit
!     !f2py intent(inout) unit
!     real(4), dimension(3,3) :: lat
!     !f2py intent(inout) lat
!     type(bas_type_ptr_type) :: bas_ptr
!     integer, intent(in), dimension(2) :: bas
!     integer, optional, intent(in) :: length
!     bas_ptr = transfer(bas, bas_ptr)
!     call geom_read(UNIT=unit, lat=lat, bas=bas_ptr%p, length=length)
! end subroutine f90wrap_rw_geom__geom_read

! subroutine f90wrap_rw_geom__geom_write(unit, lat, bas)
!     use rw_geom, only: bas_type, geom_write
!     implicit none
    
!     type bas_type_ptr_type
!         type(bas_type), pointer :: p => NULL()
!     end type bas_type_ptr_type
!     integer :: unit
!     !f2py intent(inout) unit
!     real(4), dimension(3,3) :: lat
!     !f2py intent(inout) lat
!     type(bas_type_ptr_type) :: bas_ptr
!     integer, intent(in), dimension(2) :: bas
!     bas_ptr = transfer(bas, bas_ptr)
!     call geom_write(UNIT=unit, lat=lat, bas=bas_ptr%p)
! end subroutine f90wrap_rw_geom__geom_write

! subroutine f90wrap_rw_geom__convert_bas(inbas, ret_outbas, latconv)
!     use rw_geom, only: convert_bas, bas_type
!     implicit none
    
!     type bas_type_ptr_type
!         type(bas_type), pointer :: p => NULL()
!     end type bas_type_ptr_type
!     type(bas_type_ptr_type) :: inbas_ptr
!     integer, intent(in), dimension(2) :: inbas
!     type(bas_type_ptr_type) :: ret_outbas_ptr
!     integer, intent(out), dimension(2) :: ret_outbas
!     real(4), dimension(3,3), intent(in) :: latconv
!     inbas_ptr = transfer(inbas, inbas_ptr)
!     allocate(ret_outbas_ptr%p)
!     ret_outbas_ptr%p = convert_bas(inbas=inbas_ptr%p, latconv=latconv)
!     ret_outbas = transfer(ret_outbas_ptr, ret_outbas)
! end subroutine f90wrap_rw_geom__convert_bas

! subroutine f90wrap_rw_geom__clone_bas(inbas, outbas, inlat, outlat, trans_dim)
!     use rw_geom, only: clone_bas, bas_type
!     implicit none
    
!     type bas_type_ptr_type
!         type(bas_type), pointer :: p => NULL()
!     end type bas_type_ptr_type
!     type(bas_type_ptr_type) :: inbas_ptr
!     integer, intent(in), dimension(2) :: inbas
!     type(bas_type_ptr_type) :: outbas_ptr
!     integer, intent(in), dimension(2) :: outbas
!     real(4), dimension(3,3), optional :: inlat
!     !f2py intent(inout) inlat
!     real(4), dimension(3,3), optional :: outlat
!     !f2py intent(inout) outlat
!     logical, optional, intent(in) :: trans_dim
!     inbas_ptr = transfer(inbas, inbas_ptr)
!     outbas_ptr = transfer(outbas, outbas_ptr)
!     call clone_bas(inbas=inbas_ptr%p, outbas=outbas_ptr%p, inlat=inlat, outlat=outlat, trans_dim=trans_dim)
! end subroutine f90wrap_rw_geom__clone_bas

! subroutine f90wrap_rw_geom__get__igeom_input(f90wrap_igeom_input)
!     use rw_geom, only: rw_geom_igeom_input => igeom_input
!     implicit none
!     integer, intent(out) :: f90wrap_igeom_input
    
!     f90wrap_igeom_input = rw_geom_igeom_input
! end subroutine f90wrap_rw_geom__get__igeom_input

! subroutine f90wrap_rw_geom__set__igeom_input(f90wrap_igeom_input)
!     use rw_geom, only: rw_geom_igeom_input => igeom_input
!     implicit none
!     integer, intent(in) :: f90wrap_igeom_input
    
!     rw_geom_igeom_input = f90wrap_igeom_input
! end subroutine f90wrap_rw_geom__set__igeom_input

! subroutine f90wrap_rw_geom__get__igeom_output(f90wrap_igeom_output)
!     use rw_geom, only: rw_geom_igeom_output => igeom_output
!     implicit none
!     integer, intent(out) :: f90wrap_igeom_output
    
!     f90wrap_igeom_output = rw_geom_igeom_output
! end subroutine f90wrap_rw_geom__get__igeom_output

! subroutine f90wrap_rw_geom__set__igeom_output(f90wrap_igeom_output)
!     use rw_geom, only: rw_geom_igeom_output => igeom_output
!     implicit none
!     integer, intent(in) :: f90wrap_igeom_output
    
!     rw_geom_igeom_output = f90wrap_igeom_output
! end subroutine f90wrap_rw_geom__set__igeom_output

! End of module rw_geom defined in file ../src/lib/mod_rw_geom.f90

