! Module generator defined in file ../src/lib/mod_generator.f90

subroutine f90wrap_stoichiometry_type__get__element(this, f90wrap_element)
    use generator, only: stoichiometry_type
    implicit none
    type stoichiometry_type_ptr_type
        type(stoichiometry_type), pointer :: p => NULL()
    end type stoichiometry_type_ptr_type
    integer, intent(in)   :: this(2)
    type(stoichiometry_type_ptr_type) :: this_ptr
    character(3), intent(out) :: f90wrap_element
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_element = this_ptr%p%element
end subroutine f90wrap_stoichiometry_type__get__element

subroutine f90wrap_stoichiometry_type__set__element(this, f90wrap_element)
    use generator, only: stoichiometry_type
    implicit none
    type stoichiometry_type_ptr_type
        type(stoichiometry_type), pointer :: p => NULL()
    end type stoichiometry_type_ptr_type
    integer, intent(in)   :: this(2)
    type(stoichiometry_type_ptr_type) :: this_ptr
    character(3), intent(in) :: f90wrap_element
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%element = f90wrap_element
end subroutine f90wrap_stoichiometry_type__set__element

subroutine f90wrap_stoichiometry_type__get__num(this, f90wrap_num)
    use generator, only: stoichiometry_type
    implicit none
    type stoichiometry_type_ptr_type
        type(stoichiometry_type), pointer :: p => NULL()
    end type stoichiometry_type_ptr_type
    integer, intent(in)   :: this(2)
    type(stoichiometry_type_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_num
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_num = this_ptr%p%num
end subroutine f90wrap_stoichiometry_type__get__num

subroutine f90wrap_stoichiometry_type__set__num(this, f90wrap_num)
    use generator, only: stoichiometry_type
    implicit none
    type stoichiometry_type_ptr_type
        type(stoichiometry_type), pointer :: p => NULL()
    end type stoichiometry_type_ptr_type
    integer, intent(in)   :: this(2)
    type(stoichiometry_type_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_num
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%num = f90wrap_num
end subroutine f90wrap_stoichiometry_type__set__num

subroutine f90wrap_stoichiometry_type_initialise(this)
    use generator, only: stoichiometry_type
    implicit none
    
    type stoichiometry_type_ptr_type
        type(stoichiometry_type), pointer :: p => NULL()
    end type stoichiometry_type_ptr_type
    type(stoichiometry_type_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_stoichiometry_type_initialise

subroutine f90wrap_stoichiometry_type_finalise(this)
    use generator, only: stoichiometry_type
    implicit none
    
    type stoichiometry_type_ptr_type
        type(stoichiometry_type), pointer :: p => NULL()
    end type stoichiometry_type_ptr_type
    type(stoichiometry_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_stoichiometry_type_finalise


subroutine f90wrap_stoich_type_xnum_array__array_getitem__items( &
       this, f90wrap_i, itemsitem)
    use generator, only: stoichiometry_type
    implicit none

    type stoichiometry_type_xnum_array
        type(stoichiometry_type), dimension(:), allocatable :: items
    end type stoichiometry_type_xnum_array

    type stoichiometry_type_xnum_array_ptr_type
        type(stoichiometry_type_xnum_array), pointer :: p => NULL()
    end type stoichiometry_type_xnum_array_ptr_type
    type stoichiometry_type_ptr_type
        type(stoichiometry_type), pointer :: p => NULL()
    end type stoichiometry_type_ptr_type
    integer, intent(in), dimension(2) :: this
    type(stoichiometry_type_xnum_array_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_i
    integer, intent(out) :: itemsitem(2)
    type(stoichiometry_type_ptr_type) :: items_ptr
    
    this_ptr = transfer(this, this_ptr)
    if (f90wrap_i < 1 .or. f90wrap_i > size(this_ptr%p%items)) then
        call f90wrap_abort("array index out of range")
    else
        items_ptr%p => this_ptr%p%items(f90wrap_i)
        itemsitem = transfer(items_ptr,itemsitem)
    endif
end subroutine f90wrap_stoich_type_xnum_array__array_getitem__items

subroutine f90wrap_stoich_type_xnum_array__array_setitem__items(this, f90wrap_i, itemsitem)
    use generator, only: stoichiometry_type
    implicit none

    type stoichiometry_type_xnum_array
        type(stoichiometry_type), dimension(:), allocatable :: items
    end type stoichiometry_type_xnum_array

    type stoichiometry_type_xnum_array_ptr_type
        type(stoichiometry_type_xnum_array), pointer :: p => NULL()
    end type stoichiometry_type_xnum_array_ptr_type
    type stoichiometry_type_ptr_type
        type(stoichiometry_type), pointer :: p => NULL()
    end type stoichiometry_type_ptr_type
    integer, intent(in), dimension(2) :: this
    type(stoichiometry_type_xnum_array_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_i
    integer, intent(out) :: itemsitem(2)
    type(stoichiometry_type_ptr_type) :: items_ptr
    
    this_ptr = transfer(this, this_ptr)
    if (f90wrap_i < 1 .or. f90wrap_i > size(this_ptr%p%items)) then
        call f90wrap_abort("array index out of range")
    else
        items_ptr = transfer(itemsitem,items_ptr)
        this_ptr%p%items(f90wrap_i) = items_ptr%p
    endif
end subroutine f90wrap_stoich_type_xnum_array__array_setitem__items

subroutine f90wrap_stoich_type_xnum_array__array_len__items(this, f90wrap_n)
    use generator, only: stoichiometry_type
    implicit none

    type stoichiometry_type_xnum_array
        type(stoichiometry_type), dimension(:), allocatable :: items
    end type stoichiometry_type_xnum_array

    type stoichiometry_type_xnum_array_ptr_type
        type(stoichiometry_type_xnum_array), pointer :: p => NULL()
    end type stoichiometry_type_xnum_array_ptr_type
    integer, intent(in), dimension(2) :: this
    type(stoichiometry_type_xnum_array_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_n
    this_ptr = transfer(this, this_ptr)
    f90wrap_n = size(this_ptr%p%items)
end subroutine f90wrap_stoich_type_xnum_array__array_len__items

subroutine f90wrap_stoich_type_xnum_array__array_alloc__items(this, num)
    use generator, only: stoichiometry_type
    implicit none

    type stoichiometry_type_xnum_array
        type(stoichiometry_type), dimension(:), allocatable :: items
    end type stoichiometry_type_xnum_array

    type stoichiometry_type_xnum_array_ptr_type
        type(stoichiometry_type_xnum_array), pointer :: p => NULL()
    end type stoichiometry_type_xnum_array_ptr_type
    type(stoichiometry_type_xnum_array_ptr_type) :: this_ptr
    integer, intent(in) :: num
    integer, intent(inout), dimension(2) :: this

    this_ptr = transfer(this, this_ptr)
    allocate(this_ptr%p%items(num))
    this = transfer(this_ptr, this)
end subroutine f90wrap_stoich_type_xnum_array__array_alloc__items

subroutine f90wrap_stoich_type_xnum_array__array_dealloc__items(this)
    use generator, only: stoichiometry_type
    implicit none

    type stoichiometry_type_xnum_array
        type(stoichiometry_type), dimension(:), allocatable :: items
    end type stoichiometry_type_xnum_array

    type stoichiometry_type_xnum_array_ptr_type
        type(stoichiometry_type_xnum_array), pointer :: p => NULL()
    end type stoichiometry_type_xnum_array_ptr_type
    type(stoichiometry_type_xnum_array_ptr_type) :: this_ptr
    integer, intent(inout), dimension(2) :: this

    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p%items)
    this = transfer(this_ptr, this)
end subroutine f90wrap_stoich_type_xnum_array__array_dealloc__items

subroutine f90wrap_generator__stoich_type_xnum_array_initialise(this)
    use generator, only: stoichiometry_type
    implicit none

    type stoichiometry_type_xnum_array
        type(stoichiometry_type), dimension(:), allocatable :: items
    end type stoichiometry_type_xnum_array

    type stoichiometry_type_xnum_array_ptr_type
        type(stoichiometry_type_xnum_array), pointer :: p => NULL()
    end type stoichiometry_type_xnum_array_ptr_type
    type(stoichiometry_type_xnum_array_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_generator__stoich_type_xnum_array_initialise

subroutine f90wrap_generator__stoich_type_xnum_array_finalise(this)
    use generator, only: stoichiometry_type
    implicit none

    type stoichiometry_type_xnum_array
        type(stoichiometry_type), dimension(:), allocatable :: items
    end type stoichiometry_type_xnum_array

    type stoichiometry_type_xnum_array_ptr_type
        type(stoichiometry_type_xnum_array), pointer :: p => NULL()
    end type stoichiometry_type_xnum_array_ptr_type
    type(stoichiometry_type_xnum_array_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_generator__stoich_type_xnum_array_finalise




subroutine f90wrap_raffle_generator_type__get__num_structures(this, f90wrap_num_structures)
    use generator, only: raffle_generator_type
    implicit none
    type raffle_generator_type_ptr_type
        type(raffle_generator_type), pointer :: p => NULL()
    end type raffle_generator_type_ptr_type
    integer, intent(in)   :: this(2)
    type(raffle_generator_type_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_num_structures
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_num_structures = this_ptr%p%num_structures
end subroutine f90wrap_raffle_generator_type__get__num_structures

subroutine f90wrap_raffle_generator_type__set__num_structures(this, f90wrap_num_structures)
    use generator, only: raffle_generator_type
    implicit none
    type raffle_generator_type_ptr_type
        type(raffle_generator_type), pointer :: p => NULL()
    end type raffle_generator_type_ptr_type
    integer, intent(in)   :: this(2)
    type(raffle_generator_type_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_num_structures
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%num_structures = f90wrap_num_structures
end subroutine f90wrap_raffle_generator_type__set__num_structures

subroutine f90wrap_raffle_generator_type__get__host(this, f90wrap_host)
    use generator, only: raffle_generator_type
    use rw_geom, only: basis_type
    implicit none
    type raffle_generator_type_ptr_type
        type(raffle_generator_type), pointer :: p => NULL()
    end type raffle_generator_type_ptr_type
    type basis_type_ptr_type
        type(basis_type), pointer :: p => NULL()
    end type basis_type_ptr_type
    integer, intent(in)   :: this(2)
    type(raffle_generator_type_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_host(2)
    type(basis_type_ptr_type) :: host_ptr

    this_ptr = transfer(this, this_ptr)
    host_ptr%p => this_ptr%p%host
    f90wrap_host = transfer(host_ptr, f90wrap_host)
end subroutine f90wrap_raffle_generator_type__get__host

subroutine f90wrap_raffle_generator_type__set__host(this, f90wrap_host)
    use generator, only: raffle_generator_type
    use rw_geom, only: basis_type
    implicit none
    type raffle_generator_type_ptr_type
        type(raffle_generator_type), pointer :: p => NULL()
    end type raffle_generator_type_ptr_type
    type basis_type_ptr_type
        type(basis_type), pointer :: p => NULL()
    end type basis_type_ptr_type
    integer, intent(in)   :: this(2)
    type(raffle_generator_type_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_host(2)
    type(basis_type_ptr_type) :: host_ptr

    this_ptr = transfer(this, this_ptr)
    host_ptr = transfer(f90wrap_host, host_ptr)
    this_ptr%p%host = host_ptr%p
end subroutine f90wrap_raffle_generator_type__set__host

subroutine f90wrap_raffle_generator_type__array__grid(this, nd, dtype, dshape, dloc)
    use generator, only: raffle_generator_type
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type raffle_generator_type_ptr_type
        type(raffle_generator_type), pointer :: p => NULL()
    end type raffle_generator_type_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(raffle_generator_type_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    dshape(1:1) = shape(this_ptr%p%grid)
    dloc = loc(this_ptr%p%grid)
end subroutine f90wrap_raffle_generator_type__array__grid

subroutine f90wrap_raffle_generator_type__array__grid_offset(this, nd, dtype, dshape, dloc)
    use generator, only: raffle_generator_type
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type raffle_generator_type_ptr_type
        type(raffle_generator_type), pointer :: p => NULL()
    end type raffle_generator_type_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(raffle_generator_type_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 11
    this_ptr = transfer(this, this_ptr)
    dshape(1:1) = shape(this_ptr%p%grid_offset)
    dloc = loc(this_ptr%p%grid_offset)
end subroutine f90wrap_raffle_generator_type__array__grid_offset

subroutine f90wrap_raffle_generator_type__get__grid_spacing(this, f90wrap_grid_spacing)
    use generator, only: raffle_generator_type
    implicit none
    type raffle_generator_type_ptr_type
        type(raffle_generator_type), pointer :: p => NULL()
    end type raffle_generator_type_ptr_type
    integer, intent(in)   :: this(2)
    type(raffle_generator_type_ptr_type) :: this_ptr
    real(4), intent(out) :: f90wrap_grid_spacing
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_grid_spacing = this_ptr%p%grid_spacing
end subroutine f90wrap_raffle_generator_type__get__grid_spacing

subroutine f90wrap_raffle_generator_type__set__grid_spacing(this, f90wrap_grid_spacing)
    use generator, only: raffle_generator_type
    implicit none
    type raffle_generator_type_ptr_type
        type(raffle_generator_type), pointer :: p => NULL()
    end type raffle_generator_type_ptr_type
    integer, intent(in)   :: this(2)
    type(raffle_generator_type_ptr_type) :: this_ptr
    real(4), intent(in) :: f90wrap_grid_spacing
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%grid_spacing = f90wrap_grid_spacing
end subroutine f90wrap_raffle_generator_type__set__grid_spacing

subroutine f90wrap_raffle_generator_type__get__distributions(this, f90wrap_distributions)
    use generator, only: raffle_generator_type
    use evolver, only: gvector_container_type
    implicit none
    type raffle_generator_type_ptr_type
        type(raffle_generator_type), pointer :: p => NULL()
    end type raffle_generator_type_ptr_type
    type gvector_container_type_ptr_type
        type(gvector_container_type), pointer :: p => NULL()
    end type gvector_container_type_ptr_type
    integer, intent(in)   :: this(2)
    type(raffle_generator_type_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_distributions(2)
    type(gvector_container_type_ptr_type) :: distributions_ptr
    
    this_ptr = transfer(this, this_ptr)
    distributions_ptr%p => this_ptr%p%distributions
    f90wrap_distributions = transfer(distributions_ptr,f90wrap_distributions)
end subroutine f90wrap_raffle_generator_type__get__distributions

subroutine f90wrap_raffle_generator_type__set__distributions(this, f90wrap_distributions)
    use generator, only: raffle_generator_type
    use evolver, only: gvector_container_type
    implicit none
    type raffle_generator_type_ptr_type
        type(raffle_generator_type), pointer :: p => NULL()
    end type raffle_generator_type_ptr_type
    type gvector_container_type_ptr_type
        type(gvector_container_type), pointer :: p => NULL()
    end type gvector_container_type_ptr_type
    integer, intent(in)   :: this(2)
    type(raffle_generator_type_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_distributions(2)
    type(gvector_container_type_ptr_type) :: distributions_ptr
    
    this_ptr = transfer(this, this_ptr)
    distributions_ptr = transfer(f90wrap_distributions,distributions_ptr)
    this_ptr%p%distributions = distributions_ptr%p
end subroutine f90wrap_raffle_generator_type__set__distributions

subroutine f90wrap_raffle_generator_type__array__method_probab(this, nd, dtype, dshape, dloc)
    use generator, only: raffle_generator_type
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type raffle_generator_type_ptr_type
        type(raffle_generator_type), pointer :: p => NULL()
    end type raffle_generator_type_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(raffle_generator_type_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 11
    this_ptr = transfer(this, this_ptr)
    dshape(1:1) = shape(this_ptr%p%method_probab)
    dloc = loc(this_ptr%p%method_probab)
end subroutine f90wrap_raffle_generator_type__array__method_probab

subroutine f90wrap_raffle_generator_type__array_getitem__structures(f90wrap_this, f90wrap_i, structuresitem)
    
    use generator, only: raffle_generator_type
    use rw_geom, only: basis_type
    implicit none
    
    type raffle_generator_type_ptr_type
        type(raffle_generator_type), pointer :: p => NULL()
    end type raffle_generator_type_ptr_type
    type basis_type_ptr_type
        type(basis_type), pointer :: p => NULL()
    end type basis_type_ptr_type
    integer, intent(in) :: f90wrap_this(2)
    type(raffle_generator_type_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_i
    integer, intent(out) :: structuresitem(2)
    type(basis_type_ptr_type) :: structures_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%structures)) then
        if (f90wrap_i < 1 .or. f90wrap_i > size(this_ptr%p%structures)) then
            call f90wrap_abort("array index out of range")
        else
            structures_ptr%p => this_ptr%p%structures(f90wrap_i)
            structuresitem = transfer(structures_ptr,structuresitem)
        endif
    else
        call f90wrap_abort("derived type array not allocated")
    end if
end subroutine f90wrap_raffle_generator_type__array_getitem__structures

subroutine f90wrap_raffle_generator_type__array_setitem__structures(f90wrap_this, f90wrap_i, structuresitem)
    
    use generator, only: raffle_generator_type
    use rw_geom, only: basis_type
    implicit none
    
    type raffle_generator_type_ptr_type
        type(raffle_generator_type), pointer :: p => NULL()
    end type raffle_generator_type_ptr_type
    type basis_type_ptr_type
        type(basis_type), pointer :: p => NULL()
    end type basis_type_ptr_type
    integer, intent(in) :: f90wrap_this(2)
    type(raffle_generator_type_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_i
    integer, intent(in) :: structuresitem(2)
    type(basis_type_ptr_type) :: structures_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%structures)) then
        if (f90wrap_i < 1 .or. f90wrap_i > size(this_ptr%p%structures)) then
            call f90wrap_abort("array index out of range")
        else
            structures_ptr = transfer(structuresitem,structures_ptr)
            this_ptr%p%structures(f90wrap_i) = structures_ptr%p
        endif
    else
        call f90wrap_abort("derived type array not allocated")
    end if
end subroutine f90wrap_raffle_generator_type__array_setitem__structures

subroutine f90wrap_raffle_generator_type__array_len__structures(f90wrap_this, f90wrap_n)
    
    use generator, only: raffle_generator_type
    use rw_geom, only: basis_type
    implicit none
    
    type raffle_generator_type_ptr_type
        type(raffle_generator_type), pointer :: p => NULL()
    end type raffle_generator_type_ptr_type
    integer, intent(out) :: f90wrap_n
    integer, intent(in) :: f90wrap_this(2)
    type(raffle_generator_type_ptr_type) :: this_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%structures)) then
        f90wrap_n = size(this_ptr%p%structures)
    else
        f90wrap_n = 0
    end if
end subroutine f90wrap_raffle_generator_type__array_len__structures

subroutine f90wrap_generator__raffle_generator_type_initialise(this)
    use generator, only: raffle_generator_type
    implicit none
    
    type raffle_generator_type_ptr_type
        type(raffle_generator_type), pointer :: p => NULL()
    end type raffle_generator_type_ptr_type
    type(raffle_generator_type_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_generator__raffle_generator_type_initialise

subroutine f90wrap_generator__raffle_generator_type_finalise(this)
    use generator, only: raffle_generator_type
    implicit none
    
    type raffle_generator_type_ptr_type
        type(raffle_generator_type), pointer :: p => NULL()
    end type raffle_generator_type_ptr_type
    type(raffle_generator_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_generator__raffle_generator_type_finalise

subroutine f90wrap_generator__set_host__binding__rgt(this, host)
    use rw_geom, only: basis_type
    use generator, only: raffle_generator_type
    implicit none
    
    type raffle_generator_type_ptr_type
        type(raffle_generator_type), pointer :: p => NULL()
    end type raffle_generator_type_ptr_type
    type basis_type_ptr_type
        type(basis_type), pointer :: p => NULL()
    end type basis_type_ptr_type
    type(raffle_generator_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    type(basis_type_ptr_type) :: host_ptr
    integer, intent(in), dimension(2) :: host
    this_ptr = transfer(this, this_ptr)
    host_ptr = transfer(host, host_ptr)
    call this_ptr%p%set_host(host=host_ptr%p)
end subroutine f90wrap_generator__set_host__binding__rgt

subroutine f90wrap_generator__set_grid__binding__raffle_generator_type(this, grid, grid_spacing, grid_offset)
    use generator, only: raffle_generator_type
    implicit none
    
    type raffle_generator_type_ptr_type
        type(raffle_generator_type), pointer :: p => NULL()
    end type raffle_generator_type_ptr_type
    type(raffle_generator_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    integer, dimension(3), intent(in), optional :: grid
    real(4), intent(in), optional :: grid_spacing
    real(4), dimension(3), intent(in), optional :: grid_offset
    this_ptr = transfer(this, this_ptr)
    call this_ptr%p%set_grid(grid=grid, grid_spacing=grid_spacing, grid_offset=grid_offset)
end subroutine f90wrap_generator__set_grid__binding__raffle_generator_type

subroutine f90wrap_generator__reset_grid__binding__raffle_generator_type(this)
    use generator, only: raffle_generator_type
    implicit none
    
    type raffle_generator_type_ptr_type
        type(raffle_generator_type), pointer :: p => NULL()
    end type raffle_generator_type_ptr_type
    type(raffle_generator_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    call this_ptr%p%reset_grid()
end subroutine f90wrap_generator__reset_grid__binding__raffle_generator_type

subroutine f90wrap_generator__generate__binding__rgt( &
       this, num_structures, stoichiometry, &
    method_probab, seed, verbose, n0)
    use generator, only: raffle_generator_type, stoichiometry_type
    implicit none
    
    type raffle_generator_type_ptr_type
        type(raffle_generator_type), pointer :: p => NULL()
    end type raffle_generator_type_ptr_type


    type stoichiometry_type_xnum_array
        type(stoichiometry_type), dimension(:), allocatable :: items
    end type stoichiometry_type_xnum_array

    type stoichiometry_type_xnum_array_ptr_type
        type(stoichiometry_type_xnum_array), pointer :: p => NULL()
    end type stoichiometry_type_xnum_array_ptr_type
    type(raffle_generator_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    integer, intent(in) :: num_structures
    type(stoichiometry_type_xnum_array_ptr_type) :: stoichiometry_ptr
    integer, intent(in), dimension(2) :: stoichiometry
    real(4), intent(in), optional, dimension(n0) :: method_probab
    integer, intent(in), optional :: seed
    integer, intent(in), optional :: verbose
    integer :: n0
    !f2py intent(hide), depend(method_probab) :: n0 = shape(method_probab,0)
    integer :: verbose_ = 0

    if(present(verbose)) verbose_ = verbose
    this_ptr = transfer(this, this_ptr)
    stoichiometry_ptr = transfer(stoichiometry, stoichiometry_ptr)
    if(present(method_probab)) then
        if(present(seed)) then
           call this_ptr%p%generate(num_structures=num_structures, &
                stoichiometry=stoichiometry_ptr%p%items, &
                method_probab=method_probab, seed=seed, verbose=verbose_)
        else
           call this_ptr%p%generate(num_structures=num_structures, &
                stoichiometry=stoichiometry_ptr%p%items, &
                method_probab=method_probab, verbose=verbose_)
        end if
    else
        if(present(seed)) then
           call this_ptr%p%generate(num_structures=num_structures, &
                stoichiometry=stoichiometry_ptr%p%items, &
                seed=seed, verbose=verbose_)
        else
           call this_ptr%p%generate(num_structures=num_structures, &
                stoichiometry=stoichiometry_ptr%p%items, verbose=verbose_)
        end if
    end if
end subroutine f90wrap_generator__generate__binding__rgt

subroutine f90wrap_generator__evaluate__binding__rgt(this, ret_viability, basis)
    use rw_geom, only: basis_type
    use generator, only: raffle_generator_type
    implicit none
    
    type raffle_generator_type_ptr_type
        type(raffle_generator_type), pointer :: p => NULL()
    end type raffle_generator_type_ptr_type
    type basis_type_ptr_type
        type(basis_type), pointer :: p => NULL()
    end type basis_type_ptr_type
    type(raffle_generator_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    real(4), intent(out) :: ret_viability
    type(basis_type_ptr_type) :: basis_ptr
    integer, intent(in), dimension(2) :: basis
    this_ptr = transfer(this, this_ptr)
    basis_ptr = transfer(basis, basis_ptr)
    ret_viability = this_ptr%p%evaluate(basis=basis_ptr%p)
end subroutine f90wrap_generator__evaluate__binding__rgt

subroutine f90wrap_generator__get_structures__binding__rgt(this, ret_structures)
    use rw_geom, only: basis_type
    use generator, only: raffle_generator_type
    implicit none

    type raffle_generator_type_ptr_type
        type(raffle_generator_type), pointer :: p => NULL()
    end type raffle_generator_type_ptr_type

    type basis_type_xnum_array
        type(basis_type), dimension(:), allocatable :: items
    end type basis_type_xnum_array

    type basis_type_xnum_array_ptr_type
        type(basis_type_xnum_array), pointer :: p => NULL()
    end type basis_type_xnum_array_ptr_type
    type(raffle_generator_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    integer, intent(out), dimension(2) :: ret_structures
    type(basis_type_xnum_array_ptr_type) :: ret_structures_ptr

    this_ptr = transfer(this, this_ptr)
    ret_structures_ptr%p%items = this_ptr%p%get_structures()
    ret_structures = transfer(ret_structures_ptr,ret_structures)
end subroutine f90wrap_generator__get_structures__binding__rgt

! End of module generator defined in file ../src/lib/mod_generator.f90

