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




subroutine f90wrap_raffle_generator_type__array__bins(this, nd, dtype, dshape, dloc)
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
    dshape(1:1) = shape(this_ptr%p%bins)
    dloc = loc(this_ptr%p%bins)
end subroutine f90wrap_raffle_generator_type__array__bins

subroutine f90wrap_raffle_generator_type__array__lattice_host(this, nd, dtype, dshape, dloc)
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
    
    nd = 2
    dtype = 11
    this_ptr = transfer(this, this_ptr)
    dshape(1:2) = shape(this_ptr%p%lattice_host)
    dloc = loc(this_ptr%p%lattice_host)
end subroutine f90wrap_raffle_generator_type__array__lattice_host

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

subroutine f90wrap_generator__generate__binding__rgt( &
       this, num_structures, stoichiometry, &
    method_probab, n0)
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
    integer :: n0
    !f2py intent(hide), depend(method_probab) :: n0 = shape(method_probab,0)
    write(*,*) "in generate"
    this_ptr = transfer(this, this_ptr)
    stoichiometry_ptr = transfer(stoichiometry, stoichiometry_ptr)
    if(present(method_probab)) then
        write(*,*) "method_probab present"
        call this_ptr%p%generate(num_structures=num_structures, stoichiometry=stoichiometry_ptr%p%items, &
            method_probab=method_probab)
    else
        write(*,*) "method_probab not present"
        call this_ptr%p%generate(num_structures=num_structures, stoichiometry=stoichiometry_ptr%p%items)
    end if
end subroutine f90wrap_generator__generate__binding__rgt

! End of module generator defined in file ../src/lib/mod_generator.f90

