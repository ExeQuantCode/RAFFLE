! Module generator defined in file ../src/lib/mod_generator.f90

!###############################################################################
! stoichiometry derived type
!###############################################################################
subroutine f90wrap_stoichiometry_type__get__element(this, f90wrap_element)
    use raffle__generator, only: stoichiometry_type
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
    use raffle__generator, only: stoichiometry_type
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
    use raffle__generator, only: stoichiometry_type
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
    use raffle__generator, only: stoichiometry_type
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
    use raffle__generator, only: stoichiometry_type
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
    use raffle__generator, only: stoichiometry_type
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
    use raffle__generator, only: stoichiometry_type
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
    use raffle__generator, only: stoichiometry_type
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
    use raffle__generator, only: stoichiometry_type
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
    use raffle__generator, only: stoichiometry_type
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
    use raffle__generator, only: stoichiometry_type
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
!###############################################################################


!###############################################################################
! generator contained stoichiometry
!###############################################################################
subroutine f90wrap_generator__stoich_type_xnum_array_initialise(this)
    use raffle__generator, only: stoichiometry_type
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
    use raffle__generator, only: stoichiometry_type
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
!###############################################################################


!###############################################################################
! number of generated structures
!###############################################################################
subroutine f90wrap_raffle_generator_type__get__num_structures( &
     this, f90wrap_num_structures &
)
    use raffle__generator, only: raffle_generator_type
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

subroutine f90wrap_raffle_generator_type__set__num_structures( &
     this, f90wrap_num_structures &
)
    use raffle__generator, only: raffle_generator_type
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
!###############################################################################


!###############################################################################
! host handling
!###############################################################################
subroutine f90wrap_raffle_generator_type__get__host(this, f90wrap_host)
    use raffle__generator, only: raffle_generator_type
    use raffle__geom_rw, only: basis_type
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
    use raffle__generator, only: raffle_generator_type
    use raffle__geom_rw, only: basis_type
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
!###############################################################################


!###############################################################################
! viability grid parameters
!###############################################################################
subroutine f90wrap_raffle_generator_type__array__grid( &
     this, nd, dtype, dshape, dloc &
)
    use raffle__generator, only: raffle_generator_type
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

subroutine f90wrap_raffle_generator_type__array__grid_offset( &
     this, nd, dtype, dshape, dloc &
)
    use raffle__generator, only: raffle_generator_type
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

subroutine f90wrap_raffle_generator_type__get__grid_spacing( &
     this, f90wrap_grid_spacing &
)
    use raffle__generator, only: raffle_generator_type
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

subroutine f90wrap_raffle_generator_type__set__grid_spacing( &
     this, f90wrap_grid_spacing &
)
    use raffle__generator, only: raffle_generator_type
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

subroutine f90wrap_raffle_generator_type__array__bounds( &
     this, nd, dtype, dshape, dloc &
)
    use raffle__generator, only: raffle_generator_type
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
    dshape(1:2) = shape(this_ptr%p%bounds)
    dloc = loc(this_ptr%p%bounds)
end subroutine f90wrap_raffle_generator_type__array__bounds
!###############################################################################


!###############################################################################
! distribution function handling
!###############################################################################
subroutine f90wrap_raffle_generator_type__get__distributions( &
     this, f90wrap_distributions &
)
    use raffle__generator, only: raffle_generator_type
    use raffle__distribs_container, only: distribs_container_type
    implicit none
    type raffle_generator_type_ptr_type
        type(raffle_generator_type), pointer :: p => NULL()
    end type raffle_generator_type_ptr_type
    type distribs_container_type_ptr_type
        type(distribs_container_type), pointer :: p => NULL()
    end type distribs_container_type_ptr_type
    integer, intent(in)   :: this(2)
    type(raffle_generator_type_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_distributions(2)
    type(distribs_container_type_ptr_type) :: distributions_ptr
    
    this_ptr = transfer(this, this_ptr)
    distributions_ptr%p => this_ptr%p%distributions
    f90wrap_distributions = transfer(distributions_ptr,f90wrap_distributions)
end subroutine f90wrap_raffle_generator_type__get__distributions

subroutine f90wrap_raffle_generator_type__set__distributions( &
     this, f90wrap_distributions &
)
    use raffle__generator, only: raffle_generator_type
    use raffle__distribs_container, only: distribs_container_type
    implicit none
    type raffle_generator_type_ptr_type
        type(raffle_generator_type), pointer :: p => NULL()
    end type raffle_generator_type_ptr_type
    type distribs_container_type_ptr_type
        type(distribs_container_type), pointer :: p => NULL()
    end type distribs_container_type_ptr_type
    integer, intent(in)   :: this(2)
    type(raffle_generator_type_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_distributions(2)
    type(distribs_container_type_ptr_type) :: distributions_ptr
    
    this_ptr = transfer(this, this_ptr)
    distributions_ptr = transfer(f90wrap_distributions,distributions_ptr)
    this_ptr%p%distributions = distributions_ptr%p
end subroutine f90wrap_raffle_generator_type__set__distributions
!###############################################################################


!###############################################################################
! random walk parameters
!###############################################################################
subroutine f90wrap_raffle_generator_type__get__max_attempts( &
     this, f90wrap_max_attempts &
)
    use raffle__generator, only: raffle_generator_type
    implicit none
    type raffle_generator_type_ptr_type
        type(raffle_generator_type), pointer :: p => NULL()
    end type raffle_generator_type_ptr_type
    integer, intent(in)   :: this(2)
    type(raffle_generator_type_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_max_attempts
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_max_attempts = this_ptr%p%max_attempts
end subroutine f90wrap_raffle_generator_type__get__max_attempts

subroutine f90wrap_raffle_generator_type__set__max_attempts( &
     this, f90wrap_max_attempts &
)
    use raffle__generator, only: raffle_generator_type
    implicit none
    type raffle_generator_type_ptr_type
        type(raffle_generator_type), pointer :: p => NULL()
    end type raffle_generator_type_ptr_type
    integer, intent(in)   :: this(2)
    type(raffle_generator_type_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_max_attempts
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%max_attempts = f90wrap_max_attempts
end subroutine f90wrap_raffle_generator_type__set__max_attempts

subroutine f90wrap_raffle_generator_type__get__walk_step_size_coarse( &
     this, f90wrap_walk_step_size_coarse &
)
    use raffle__generator, only: raffle_generator_type
    implicit none
    type raffle_generator_type_ptr_type
        type(raffle_generator_type), pointer :: p => NULL()
    end type raffle_generator_type_ptr_type
    integer, intent(in)   :: this(2)
    type(raffle_generator_type_ptr_type) :: this_ptr
    real(4), intent(out) :: f90wrap_walk_step_size_coarse
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_walk_step_size_coarse = this_ptr%p%walk_step_size_coarse
end subroutine f90wrap_raffle_generator_type__get__walk_step_size_coarse

subroutine f90wrap_raffle_generator_type__set__walk_step_size_coarse( &
     this, f90wrap_walk_step_size_coarse &
)
    use raffle__generator, only: raffle_generator_type
    implicit none
    type raffle_generator_type_ptr_type
        type(raffle_generator_type), pointer :: p => NULL()
    end type raffle_generator_type_ptr_type
    integer, intent(in)   :: this(2)
    type(raffle_generator_type_ptr_type) :: this_ptr
    real(4), intent(in) :: f90wrap_walk_step_size_coarse
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%walk_step_size_coarse = f90wrap_walk_step_size_coarse
end subroutine f90wrap_raffle_generator_type__set__walk_step_size_coarse

subroutine f90wrap_raffle_generator_type__get__walk_step_size_fine( &
     this, f90wrap_walk_step_size_fine &
)
    use raffle__generator, only: raffle_generator_type
    implicit none
    type raffle_generator_type_ptr_type
        type(raffle_generator_type), pointer :: p => NULL()
    end type raffle_generator_type_ptr_type
    integer, intent(in)   :: this(2)
    type(raffle_generator_type_ptr_type) :: this_ptr
    real(4), intent(out) :: f90wrap_walk_step_size_fine
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_walk_step_size_fine = this_ptr%p%walk_step_size_fine
end subroutine f90wrap_raffle_generator_type__get__walk_step_size_fine

subroutine f90wrap_raffle_generator_type__set__walk_step_size_fine( &
     this, f90wrap_walk_step_size_fine &
)
    use raffle__generator, only: raffle_generator_type
    implicit none
    type raffle_generator_type_ptr_type
        type(raffle_generator_type), pointer :: p => NULL()
    end type raffle_generator_type_ptr_type
    integer, intent(in)   :: this(2)
    type(raffle_generator_type_ptr_type) :: this_ptr
    real(4), intent(in) :: f90wrap_walk_step_size_fine
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%walk_step_size_fine = f90wrap_walk_step_size_fine
end subroutine f90wrap_raffle_generator_type__set__walk_step_size_fine
!###############################################################################


!###############################################################################
! placement method ratio
!###############################################################################
subroutine f90wrap_raffle_generator_type__array__method_ratio( &
     this, nd, dtype, dshape, dloc &
)
    use raffle__generator, only: raffle_generator_type
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
    dshape(1:1) = shape(this_ptr%p%method_ratio)
    dloc = loc(this_ptr%p%method_ratio)
end subroutine f90wrap_raffle_generator_type__array__method_ratio
!###############################################################################


!###############################################################################
! generated structures handling
!###############################################################################
subroutine f90wrap_raffle_generator_type__array_getitem__structures( &
     f90wrap_this, f90wrap_i, structuresitem &
)
    
    use raffle__generator, only: raffle_generator_type
    use raffle__geom_rw, only: basis_type
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

subroutine f90wrap_raffle_generator_type__array_setitem__structures( &
     f90wrap_this, f90wrap_i, structuresitem &
)
    
    use raffle__generator, only: raffle_generator_type
    use raffle__geom_rw, only: basis_type
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

subroutine f90wrap_raffle_generator_type__array_len__structures( &
     f90wrap_this, f90wrap_n &
)
    
    use raffle__generator, only: raffle_generator_type
    use raffle__geom_rw, only: basis_type
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
!###############################################################################


!###############################################################################
! generator derived type initialisation and finalisation
!###############################################################################
subroutine f90wrap_generator__raffle_generator_type_initialise(this)
    use raffle__generator, only: raffle_generator_type
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
    use raffle__generator, only: raffle_generator_type
    implicit none
    
    type raffle_generator_type_ptr_type
        type(raffle_generator_type), pointer :: p => NULL()
    end type raffle_generator_type_ptr_type
    type(raffle_generator_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_generator__raffle_generator_type_finalise
!###############################################################################


!###############################################################################
! generator type procedure bindings
!###############################################################################
subroutine f90wrap_generator__set_host__binding__rgt(this, host)
    use raffle__geom_rw, only: basis_type
    use raffle__generator, only: raffle_generator_type
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

subroutine f90wrap_generator__set_grid__binding__raffle_generator_type( &
     this, grid, grid_spacing, grid_offset &
)
    use raffle__generator, only: raffle_generator_type
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
    call this_ptr%p%set_grid( &
         grid = grid, &
         grid_spacing = grid_spacing, &
         grid_offset = grid_offset &
    )
end subroutine f90wrap_generator__set_grid__binding__raffle_generator_type

subroutine f90wrap_generator__reset_grid__binding__raffle_generator_type(this)
    use raffle__generator, only: raffle_generator_type
    implicit none
    
    type raffle_generator_type_ptr_type
        type(raffle_generator_type), pointer :: p => NULL()
    end type raffle_generator_type_ptr_type
    type(raffle_generator_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    call this_ptr%p%reset_grid()
end subroutine f90wrap_generator__reset_grid__binding__raffle_generator_type

subroutine f90wrap_generator__set_bounds__binding__rgt(this, bounds)
    use raffle__generator, only: raffle_generator_type
    implicit none
    
    type raffle_generator_type_ptr_type
        type(raffle_generator_type), pointer :: p => NULL()
    end type raffle_generator_type_ptr_type
    type(raffle_generator_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    real(4), dimension(2,3), intent(in) :: bounds
    this_ptr = transfer(this, this_ptr)
    call this_ptr%p%set_bounds(bounds=bounds)
end subroutine f90wrap_generator__set_bounds__binding__rgt

subroutine f90wrap_generator__reset_bounds__binding__rgt(this)
    use raffle__generator, only: raffle_generator_type
    implicit none
    
    type raffle_generator_type_ptr_type
        type(raffle_generator_type), pointer :: p => NULL()
    end type raffle_generator_type_ptr_type
    type(raffle_generator_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    call this_ptr%p%reset_bounds()
end subroutine f90wrap_generator__reset_bounds__binding__rgt

subroutine f90wrap_generator__generate__binding__rgt( &
       this, num_structures, stoichiometry, &
       method_ratio, seed, settings_out_file, verbose, exit_code &
)
    use raffle__geom_rw, only: basis_type
    use raffle__generator, only: raffle_generator_type, stoichiometry_type
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
    real(4), intent(in), optional, dimension(5) :: method_ratio
    character*(*), intent(in), optional :: settings_out_file
    integer, intent(in), optional :: seed
    integer, intent(in), optional :: verbose
    integer, intent(out), optional :: exit_code

    this_ptr = transfer(this, this_ptr)
    stoichiometry_ptr = transfer(stoichiometry, stoichiometry_ptr)
    call this_ptr%p%generate( &
         num_structures=num_structures, &
         stoichiometry=stoichiometry_ptr%p%items, &
         method_ratio=method_ratio, &
         seed=seed, &
         settings_out_file=settings_out_file, &
         verbose=verbose, &
         exit_code=exit_code &
    )
end subroutine f90wrap_generator__generate__binding__rgt

subroutine f90wrap_generator__get_structures__binding__rgt(this, ret_structures)
    use raffle__geom_rw, only: basis_type
    use raffle__generator, only: raffle_generator_type
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

subroutine f90wrap_generator__set_structures__binding__rgt(this, structures)
    use raffle__geom_rw, only: basis_type
    use raffle__generator, only: raffle_generator_type
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
    integer, intent(in), dimension(2) :: structures
    type(basis_type_xnum_array_ptr_type) :: structures_ptr

    this_ptr = transfer(this, this_ptr)
    structures_ptr = transfer(structures, structures_ptr)
    call this_ptr%p%set_structures(structures_ptr%p%items)
end subroutine f90wrap_generator__set_structures__binding__rgt

subroutine f90wrap_generator__remove_structure__binding__rgt(this, index_bn, n0)
    use raffle__generator, only: raffle_generator_type
    implicit none
    
    type raffle_generator_type_ptr_type
        type(raffle_generator_type), pointer :: p => NULL()
    end type raffle_generator_type_ptr_type
    type(raffle_generator_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    integer, dimension(n0), intent(in) :: index_bn
    integer :: n0
    !f2py intent(hide), depend(energy_above_hull_list) :: n0 = shape(energy_above_hull_list,0)
    this_ptr = transfer(this, this_ptr)
    call this_ptr%p%remove_structure(index=index_bn)
end subroutine f90wrap_generator__remove_structure__binding__rgt

subroutine f90wrap_generator__evaluate__binding__rgt(this, ret_viability, basis)
    use raffle__geom_rw, only: basis_type
    use raffle__generator, only: raffle_generator_type
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


subroutine f90wrap_generator__print_settings__binding__rgt( &
    this, file &
)
   use raffle__generator, only: raffle_generator_type
   implicit none
   
   type raffle_generator_type_ptr_type
       type(raffle_generator_type), pointer :: p => NULL()
   end type raffle_generator_type_ptr_type
   type(raffle_generator_type_ptr_type) :: this_ptr
   integer, intent(in), dimension(2) :: this
   character*(*), intent(in) :: file
   this_ptr = transfer(this, this_ptr)
   call this_ptr%p%print_settings(file=file)
end subroutine f90wrap_generator__print_settings__binding__rgt


subroutine f90wrap_generator__read_settings__binding__rgt( &
    this, file &
)
   use raffle__generator, only: raffle_generator_type
   implicit none
   
   type raffle_generator_type_ptr_type
       type(raffle_generator_type), pointer :: p => NULL()
   end type raffle_generator_type_ptr_type
   type(raffle_generator_type_ptr_type) :: this_ptr
   integer, intent(in), dimension(2) :: this
   character*(*), intent(in) :: file
   this_ptr = transfer(this, this_ptr)
   call this_ptr%p%read_settings(file=file)
end subroutine f90wrap_generator__read_settings__binding__rgt
!###############################################################################

! End of module generator defined in file ../src/lib/mod_generator.f90

