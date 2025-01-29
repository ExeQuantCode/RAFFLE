! Module raffle__distribs_container defined in file ../src/lib/mod_distribs_container.f90

!###############################################################################
! number of evaluated systems
!###############################################################################
subroutine f90wrap_distribs_container_type__get__num_evaluated( &
     this, f90wrap_num_evaluated &
)
    use raffle__distribs_container, only: distribs_container_type
    implicit none
    type distribs_container_type_ptr_type
        type(distribs_container_type), pointer :: p => NULL()
    end type distribs_container_type_ptr_type
    integer, intent(in)   :: this(2)
    type(distribs_container_type_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_num_evaluated
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_num_evaluated = this_ptr%p%num_evaluated
end subroutine f90wrap_distribs_container_type__get__num_evaluated

subroutine f90wrap_distribs_container_type__set__num_evaluated( &
     this, f90wrap_num_evaluated &
)
    use raffle__distribs_container, only: distribs_container_type
    implicit none
    type distribs_container_type_ptr_type
        type(distribs_container_type), pointer :: p => NULL()
    end type distribs_container_type_ptr_type
    integer, intent(in)   :: this(2)
    type(distribs_container_type_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_num_evaluated
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%num_evaluated = f90wrap_num_evaluated
end subroutine f90wrap_distribs_container_type__set__num_evaluated

subroutine f90wrap_distribs_container_type__get__num_evaluated_allocated( &
     this, f90wrap_num_evaluated_allocated &
)
    use raffle__distribs_container, only: distribs_container_type
    implicit none
    type distribs_container_type_ptr_type
        type(distribs_container_type), pointer :: p => NULL()
    end type distribs_container_type_ptr_type
    integer, intent(in)   :: this(2)
    type(distribs_container_type_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_num_evaluated_allocated
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_num_evaluated_allocated = this_ptr%p%num_evaluated_allocated
end subroutine f90wrap_distribs_container_type__get__num_evaluated_allocated

subroutine f90wrap_distribs_container_type__set__num_evaluated_allocated( &
     this, f90wrap_num_evaluated_allocated &
)
    use raffle__distribs_container, only: distribs_container_type
    implicit none
    type distribs_container_type_ptr_type
        type(distribs_container_type), pointer :: p => NULL()
    end type distribs_container_type_ptr_type
    integer, intent(in)   :: this(2)
    type(distribs_container_type_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_num_evaluated_allocated
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%num_evaluated_allocated = f90wrap_num_evaluated_allocated
end subroutine f90wrap_distribs_container_type__set__num_evaluated_allocated
!###############################################################################


!###############################################################################
! set energy scaling
!###############################################################################
subroutine f90wrap_distribs_container_type__get__kBT(this, f90wrap_kBT)
    use raffle__distribs_container, only: distribs_container_type
    implicit none
    type distribs_container_type_ptr_type
        type(distribs_container_type), pointer :: p => NULL()
    end type distribs_container_type_ptr_type
    integer, intent(in)   :: this(2)
    type(distribs_container_type_ptr_type) :: this_ptr
    real(4), intent(out) :: f90wrap_kBT
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_kBT = this_ptr%p%kBT
end subroutine f90wrap_distribs_container_type__get__kBT

subroutine f90wrap_distribs_container_type__set__kBT(this, f90wrap_kBT)
    use raffle__distribs_container, only: distribs_container_type
    implicit none
    type distribs_container_type_ptr_type
        type(distribs_container_type), pointer :: p => NULL()
    end type distribs_container_type_ptr_type
    integer, intent(in)   :: this(2)
    type(distribs_container_type_ptr_type) :: this_ptr
    real(4), intent(in) :: f90wrap_kBT
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%kBT = f90wrap_kBT
end subroutine f90wrap_distribs_container_type__set__kBT
!###############################################################################


!###############################################################################
! boolean for using hull weighting or empirical method
!###############################################################################
subroutine f90wrap_distribs_container_type__get__weight_by_hull( &
     this, f90wrap_weight_by_hull &
)
    use raffle__distribs_container, only: distribs_container_type
    implicit none
    type distribs_container_type_ptr_type
        type(distribs_container_type), pointer :: p => NULL()
    end type distribs_container_type_ptr_type
    integer, intent(in)   :: this(2)
    type(distribs_container_type_ptr_type) :: this_ptr
    logical, intent(out) :: f90wrap_weight_by_hull
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_weight_by_hull = this_ptr%p%weight_by_hull
end subroutine f90wrap_distribs_container_type__get__weight_by_hull

subroutine f90wrap_distribs_container_type__set__weight_by_hull( &
     this, f90wrap_weight_by_hull &
)
    use raffle__distribs_container, only: distribs_container_type
    implicit none
    type distribs_container_type_ptr_type
        type(distribs_container_type), pointer :: p => NULL()
    end type distribs_container_type_ptr_type
    integer, intent(in)   :: this(2)
    type(distribs_container_type_ptr_type) :: this_ptr
    logical, intent(in) :: f90wrap_weight_by_hull
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%weight_by_hull = f90wrap_weight_by_hull
end subroutine f90wrap_distribs_container_type__set__weight_by_hull
!###############################################################################


!###############################################################################
! viability default values
!###############################################################################
subroutine f90wrap_distribs_container_type__get__viability_3body_default( &
     this, f90wrap_viability_3body_default &
)
    use raffle__distribs_container, only: distribs_container_type
    implicit none
    type distribs_container_type_ptr_type
        type(distribs_container_type), pointer :: p => NULL()
    end type distribs_container_type_ptr_type
    integer, intent(in)   :: this(2)
    type(distribs_container_type_ptr_type) :: this_ptr
    real(4), intent(out) :: f90wrap_viability_3body_default
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_viability_3body_default = this_ptr%p%viability_3body_default
end subroutine f90wrap_distribs_container_type__get__viability_3body_default

subroutine f90wrap_distribs_container_type__set__viability_3body_default( &
     this, f90wrap_viability_3body_default &
)
    use raffle__distribs_container, only: distribs_container_type
    implicit none
    type distribs_container_type_ptr_type
        type(distribs_container_type), pointer :: p => NULL()
    end type distribs_container_type_ptr_type
    integer, intent(in)   :: this(2)
    type(distribs_container_type_ptr_type) :: this_ptr
    real(4), intent(in) :: f90wrap_viability_3body_default
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%viability_3body_default = f90wrap_viability_3body_default
end subroutine f90wrap_distribs_container_type__set__viability_3body_default

subroutine f90wrap_distribs_container_type__get__viability_4body_default( &
     this, f90wrap_viability_4body_default &
)
    use raffle__distribs_container, only: distribs_container_type
    implicit none
    type distribs_container_type_ptr_type
        type(distribs_container_type), pointer :: p => NULL()
    end type distribs_container_type_ptr_type
    integer, intent(in)   :: this(2)
    type(distribs_container_type_ptr_type) :: this_ptr
    real(4), intent(out) :: f90wrap_viability_4body_default
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_viability_4body_default = this_ptr%p%viability_4body_default
end subroutine f90wrap_distribs_container_type__get__viability_4body_default

subroutine f90wrap_distribs_container_type__set__viability_4body_default( &
     this, f90wrap_viability_4body_default &
)
    use raffle__distribs_container, only: distribs_container_type
    implicit none
    type distribs_container_type_ptr_type
        type(distribs_container_type), pointer :: p => NULL()
    end type distribs_container_type_ptr_type
    integer, intent(in)   :: this(2)
    type(distribs_container_type_ptr_type) :: this_ptr
    real(4), intent(in) :: f90wrap_viability_4body_default
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%viability_4body_default = f90wrap_viability_4body_default
end subroutine f90wrap_distribs_container_type__set__viability_4body_default
!###############################################################################


!###############################################################################
! parameters for the distribution functions
!###############################################################################
subroutine f90wrap_distribs_container_type__array__nbins( &
     this, nd, dtype, dshape, dloc &
)
    use raffle__distribs_container, only: distribs_container_type
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type distribs_container_type_ptr_type
        type(distribs_container_type), pointer :: p => NULL()
    end type distribs_container_type_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(distribs_container_type_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    dshape(1:1) = shape(this_ptr%p%nbins)
    dloc = loc(this_ptr%p%nbins)
end subroutine f90wrap_distribs_container_type__array__nbins

subroutine f90wrap_distribs_container_type__array__sigma( &
     this, nd, dtype, dshape, dloc &
)
    use raffle__distribs_container, only: distribs_container_type
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type distribs_container_type_ptr_type
        type(distribs_container_type), pointer :: p => NULL()
    end type distribs_container_type_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(distribs_container_type_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 11
    this_ptr = transfer(this, this_ptr)
    dshape(1:1) = shape(this_ptr%p%sigma)
    dloc = loc(this_ptr%p%sigma)
end subroutine f90wrap_distribs_container_type__array__sigma

subroutine f90wrap_distribs_container_type__array__width( &
     this, nd, dtype, dshape, dloc &
)
    use raffle__distribs_container, only: distribs_container_type
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type distribs_container_type_ptr_type
        type(distribs_container_type), pointer :: p => NULL()
    end type distribs_container_type_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(distribs_container_type_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 11
    this_ptr = transfer(this, this_ptr)
    dshape(1:1) = shape(this_ptr%p%width)
    dloc = loc(this_ptr%p%width)
end subroutine f90wrap_distribs_container_type__array__width

subroutine f90wrap_distribs_container_type__array__cutoff_min( &
     this, nd, dtype, dshape, dloc &
)
    use raffle__distribs_container, only: distribs_container_type
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type distribs_container_type_ptr_type
        type(distribs_container_type), pointer :: p => NULL()
    end type distribs_container_type_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(distribs_container_type_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 11
    this_ptr = transfer(this, this_ptr)
    dshape(1:1) = shape(this_ptr%p%cutoff_min)
    dloc = loc(this_ptr%p%cutoff_min)
end subroutine f90wrap_distribs_container_type__array__cutoff_min

subroutine f90wrap_distribs_container_type__array__cutoff_max( &
     this, nd, dtype, dshape, dloc &
)
    use raffle__distribs_container, only: distribs_container_type
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type distribs_container_type_ptr_type
        type(distribs_container_type), pointer :: p => NULL()
    end type distribs_container_type_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(distribs_container_type_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 11
    this_ptr = transfer(this, this_ptr)
    dshape(1:1) = shape(this_ptr%p%cutoff_max)
    dloc = loc(this_ptr%p%cutoff_max)
end subroutine f90wrap_distribs_container_type__array__cutoff_max

subroutine f90wrap_distribs_container_type__array__radius_distance_tol( &
     this, nd, dtype, dshape, dloc &
)
    use raffle__distribs_container, only: distribs_container_type
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type distribs_container_type_ptr_type
        type(distribs_container_type), pointer :: p => NULL()
    end type distribs_container_type_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(distribs_container_type_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 11
    this_ptr = transfer(this, this_ptr)
    dshape(1:1) = shape(this_ptr%p%radius_distance_tol)
    dloc = loc(this_ptr%p%radius_distance_tol)
end subroutine f90wrap_distribs_container_type__array__radius_distance_tol
!###############################################################################


!###############################################################################
! distributions container type initialiser and finaliser
!###############################################################################
subroutine f90wrap_raffle__dc__dc_type_initialise(this)
    use raffle__distribs_container, only: distribs_container_type
    implicit none
    
    type distribs_container_type_ptr_type
        type(distribs_container_type), pointer :: p => NULL()
    end type distribs_container_type_ptr_type
    type(distribs_container_type_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_raffle__dc__dc_type_initialise

subroutine f90wrap_raffle__dc__dc_type_finalise(this)
    use raffle__distribs_container, only: distribs_container_type
    implicit none
    
    type distribs_container_type_ptr_type
        type(distribs_container_type), pointer :: p => NULL()
    end type distribs_container_type_ptr_type
    type(distribs_container_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_raffle__dc__dc_type_finalise
!###############################################################################


!###############################################################################
! procedures to set distribution function parameters
!###############################################################################
subroutine f90wrap_raffle__dc__set_width__binding__dc_type( &
     this, width &
)
    use raffle__distribs_container, only: distribs_container_type
    implicit none
    
    type distribs_container_type_ptr_type
        type(distribs_container_type), pointer :: p => NULL()
    end type distribs_container_type_ptr_type
    type(distribs_container_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    real(4), dimension(3), intent(in) :: width
    this_ptr = transfer(this, this_ptr)
    call this_ptr%p%set_width(width=width)
end subroutine f90wrap_raffle__dc__set_width__binding__dc_type

subroutine f90wrap_raffle__dc__set_sigma__binding__dc_type( &
     this, sigma &
)
    use raffle__distribs_container, only: distribs_container_type
    implicit none
    
    type distribs_container_type_ptr_type
        type(distribs_container_type), pointer :: p => NULL()
    end type distribs_container_type_ptr_type
    type(distribs_container_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    real(4), dimension(3), intent(in) :: sigma
    this_ptr = transfer(this, this_ptr)
    call this_ptr%p%set_sigma(sigma=sigma)
end subroutine f90wrap_raffle__dc__set_sigma__binding__dc_type

subroutine f90wrap_raffle__dc__set_cutoff_min__binding__dc_type( &
     this, cutoff_min &
)
    use raffle__distribs_container, only: distribs_container_type
    implicit none
    
    type distribs_container_type_ptr_type
        type(distribs_container_type), pointer :: p => NULL()
    end type distribs_container_type_ptr_type
    type(distribs_container_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    real(4), dimension(3), intent(in) :: cutoff_min
    this_ptr = transfer(this, this_ptr)
    call this_ptr%p%set_cutoff_min(cutoff_min=cutoff_min)
end subroutine f90wrap_raffle__dc__set_cutoff_min__binding__dc_type

subroutine f90wrap_raffle__dc__set_cutoff_max__binding__dc_type( &
     this, cutoff_max &
)
    use raffle__distribs_container, only: distribs_container_type
    implicit none
    
    type distribs_container_type_ptr_type
        type(distribs_container_type), pointer :: p => NULL()
    end type distribs_container_type_ptr_type
    type(distribs_container_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    real(4), dimension(3), intent(in) :: cutoff_max
    this_ptr = transfer(this, this_ptr)
    call this_ptr%p%set_cutoff_max(cutoff_max=cutoff_max)
end subroutine f90wrap_raffle__dc__set_cutoff_max__binding__dc_type

subroutine f90wrap_raffle__dc__set_radius_distance_tol__binding__dc_type( &
     this, radius_distance_tol &
)
    use raffle__distribs_container, only: distribs_container_type
    implicit none
    
    type distribs_container_type_ptr_type
        type(distribs_container_type), pointer :: p => NULL()
    end type distribs_container_type_ptr_type
    type(distribs_container_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    real(4), dimension(4), intent(in) :: radius_distance_tol
    this_ptr = transfer(this, this_ptr)
    call this_ptr%p%set_radius_distance_tol( &
         radius_distance_tol=radius_distance_tol &
    )
end subroutine f90wrap_raffle__dc__set_radius_distance_tol__binding__dc_type
!###############################################################################


!###############################################################################
! create and update the generalised distribution functions
!###############################################################################
subroutine f90wrap_raffle__dc__create__binding__dc_type( &
     this, basis_list, deallocate_systems, verbose, energy_above_hull_list, n0 &
)
    use raffle__geom_rw, only: basis_type
    use raffle__distribs_container, only: distribs_container_type
    implicit none
    
    type distribs_container_type_ptr_type
        type(distribs_container_type), pointer :: p => NULL()
    end type distribs_container_type_ptr_type

    type basis_type_xnum_array
        type(basis_type), dimension(:), allocatable :: items
    end type basis_type_xnum_array

    type basis_type_xnum_array_ptr_type
        type(basis_type_xnum_array), pointer :: p => NULL()
    end type basis_type_xnum_array_ptr_type
    type(distribs_container_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    type(basis_type_xnum_array_ptr_type) :: basis_list_ptr
    integer, intent(in), dimension(2) :: basis_list
    logical, intent(in), optional :: deallocate_systems
    integer, intent(in), optional :: verbose
    real(4), dimension(n0), intent(in), optional :: energy_above_hull_list
    integer :: n0
    !f2py intent(hide), depend(energy_above_hull_list) :: n0 = shape(energy_above_hull_list,0)

    this_ptr = transfer(this, this_ptr)
    basis_list_ptr = transfer(basis_list, basis_list_ptr)
    call this_ptr%p%create( &
         basis_list=basis_list_ptr%p%items, &
         energy_above_hull_list=energy_above_hull_list, &
         deallocate_systems=deallocate_systems, &
         verbose=verbose &
    )
end subroutine f90wrap_raffle__dc__create__binding__dc_type

subroutine f90wrap_raffle__dc__update__binding__dc_type( &
       this, basis_list, &
         from_host, deallocate_systems, verbose, energy_above_hull_list, n0 &
)
    use raffle__geom_rw, only: basis_type
    use raffle__distribs_container, only: distribs_container_type
    implicit none
    
    type distribs_container_type_ptr_type
        type(distribs_container_type), pointer :: p => NULL()
    end type distribs_container_type_ptr_type

    type basis_type_xnum_array
        type(basis_type), dimension(:), allocatable :: items
    end type basis_type_xnum_array

    type basis_type_xnum_array_ptr_type
        type(basis_type_xnum_array), pointer :: p => NULL()
    end type basis_type_xnum_array_ptr_type
    type(distribs_container_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    type(basis_type_xnum_array_ptr_type) :: basis_list_ptr
    integer, intent(in), dimension(2) :: basis_list
    logical, intent(in), optional :: from_host
    logical, intent(in), optional :: deallocate_systems
    integer, intent(in), optional :: verbose
    real(4), dimension(n0), intent(in), optional :: energy_above_hull_list
    integer :: n0
    !f2py intent(hide), depend(energy_above_hull_list) :: n0 = shape(energy_above_hull_list,0)

    this_ptr = transfer(this, this_ptr)
    basis_list_ptr = transfer(basis_list, basis_list_ptr)
    call this_ptr%p%update(basis_list=basis_list_ptr%p%items, &
         energy_above_hull_list=energy_above_hull_list, &
         from_host=from_host, &
         deallocate_systems=deallocate_systems, &
         verbose=verbose &
    )
end subroutine f90wrap_raffle__dc__update__binding__dc_type
!###############################################################################


!###############################################################################
! deallocate systems procedure
!###############################################################################
subroutine f90wrap_raffle__dc__deallocate_systems__binding__dc_type( &
     this &
)
    use raffle__distribs_container, only: distribs_container_type
    implicit none
    
    type distribs_container_type_ptr_type
        type(distribs_container_type), pointer :: p => NULL()
    end type distribs_container_type_ptr_type
    type(distribs_container_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    call this_ptr%p%deallocate_systems()
end subroutine f90wrap_raffle__dc__deallocate_systems__binding__dc_type
!###############################################################################


!###############################################################################
! add an individual basis to the set of distribution functions
! this does not update the generalised distribution function
!###############################################################################
subroutine f90wrap_raffle__dc__add_basis__binding__dc_type( &
     this, basis &
)
    use raffle__geom_rw, only: basis_type
    use raffle__distribs_container, only: distribs_container_type
    implicit none
    
    type distribs_container_type_ptr_type
        type(distribs_container_type), pointer :: p => NULL()
    end type distribs_container_type_ptr_type
    type basis_type_ptr_type
        type(basis_type), pointer :: p => NULL()
    end type basis_type_ptr_type
    type(distribs_container_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    type(basis_type_ptr_type) :: basis_ptr
    integer, intent(in), dimension(2) :: basis
    this_ptr = transfer(this, this_ptr)
    basis_ptr = transfer(basis, basis_ptr)
    call this_ptr%p%add_basis(basis=basis_ptr%p)
end subroutine f90wrap_raffle__dc__add_basis__binding__dc_type
!###############################################################################


!###############################################################################
! get the number of elements in the distribution container
!###############################################################################
subroutine f90wrap_raffle__dc__get__num_elements( &
     this, ret_num_elements &
)
    use raffle__distribs_container, only: distribs_container_type
    implicit none
    type distribs_container_type_ptr_type
        type(distribs_container_type), pointer :: p => NULL()
    end type distribs_container_type_ptr_type
    integer, intent(in)   :: this(2)
    type(distribs_container_type_ptr_type) :: this_ptr
    integer, intent(out) :: ret_num_elements
    
    this_ptr = transfer(this, this_ptr)
    if(.not.allocated(this_ptr%p%element_info)) then
        ret_num_elements = 0
    else
        ret_num_elements = size(this_ptr%p%element_info,1)
    end if
end subroutine f90wrap_raffle__dc__get__num_elements
!###############################################################################

!###############################################################################
! handle element reference energies and element pair bond radii
!###############################################################################
subroutine f90wrap_raffle__dc__set_element_energy__binding__dc_type( &
     this, element, energy &
)
    use raffle__distribs_container, only: distribs_container_type
    implicit none
    
    type distribs_container_type_ptr_type
        type(distribs_container_type), pointer :: p => NULL()
    end type distribs_container_type_ptr_type
    type(distribs_container_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    character(3), intent(in) :: element
    real(4), intent(in) :: energy
    this_ptr = transfer(this, this_ptr)
    call this_ptr%p%set_element_energy(element=element, energy=energy)
end subroutine f90wrap_raffle__dc__set_element_energy__binding__dc_type

subroutine f90wrap_raffle__dc__set_element_energies__binding__dc_type( &
     this, elements, energies, n0 &
)
    use raffle__distribs_container, only: distribs_container_type
    implicit none
    
    type distribs_container_type_ptr_type
        type(distribs_container_type), pointer :: p => NULL()
    end type distribs_container_type_ptr_type
    type(distribs_container_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    character(3), dimension(n0), intent(in) :: elements
    real(4), dimension(n0), intent(in) :: energies
    integer :: n0
    !f2py intent(hide), depend(elements) :: n0 = shape(elements,0)
    this_ptr = transfer(this, this_ptr)
    call this_ptr%p%set_element_energies(elements=elements, energies=energies)
end subroutine f90wrap_raffle__dc__set_element_energies__binding__dc_type

subroutine f90wrap_raffle__dc__get_element_energies_sm__binding__dc_type( &
     this, elements, energies, n0 &
)
    use raffle__distribs_container, only: distribs_container_type
    implicit none
    
    type distribs_container_type_ptr_type
        type(distribs_container_type), pointer :: p => NULL()
    end type distribs_container_type_ptr_type
    type(distribs_container_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    character(3), intent(inout), dimension(n0) :: elements
    real(4), intent(inout), dimension(n0) :: energies
    integer :: n0
    !f2py intent(hide), depend(elements) :: n0 = shape(elements,0)
    this_ptr = transfer(this, this_ptr)
    call this_ptr%p%get_element_energies_staticmem(elements=elements, energies=energies)
end subroutine f90wrap_raffle__dc__get_element_energies_sm__binding__dc_type


subroutine f90wrap_raffle__dc__set_bond_radius__binding__dc_type( &
     this, elements, radius &
)
    use raffle__distribs_container, only: distribs_container_type
    implicit none
    
    type distribs_container_type_ptr_type
        type(distribs_container_type), pointer :: p => NULL()
    end type distribs_container_type_ptr_type
    type(distribs_container_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    character(3), dimension(2), intent(in) :: elements
    real(4), intent(in) :: radius
    this_ptr = transfer(this, this_ptr)
    call this_ptr%p%set_bond_radius(elements=elements, radius=radius)
end subroutine f90wrap_raffle__dc__set_bond_radius__binding__dc_type

subroutine f90wrap_raffle__dc__set_bond_radii__binding__dc_type( &
     this, elements, radii, n0, n1, n2 &
)
    use raffle__distribs_container, only: distribs_container_type
    implicit none
    
    type distribs_container_type_ptr_type
        type(distribs_container_type), pointer :: p => NULL()
    end type distribs_container_type_ptr_type
    type(distribs_container_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    character(3), intent(in), dimension(n0,n1) :: elements
    real(4), intent(in), dimension(n2) :: radii
    integer :: n0
    !f2py intent(hide), depend(elements) :: n0 = shape(elements,0)
    integer :: n1
    !f2py intent(hide), depend(elements) :: n1 = shape(elements,1)
    integer :: n2
    !f2py intent(hide), depend(radii) :: n2 = shape(radii,0)
    this_ptr = transfer(this, this_ptr)
    call this_ptr%p%set_bond_radii(elements=elements, radii=radii)
end subroutine f90wrap_raffle__dc__set_bond_radii__binding__dc_type

subroutine f90wrap_raffle__dc__get_bond_radii_staticmem__binding__dc_type( &
     this, elements, radii, n0 &
)
    use raffle__distribs_container, only: distribs_container_type
    implicit none
    
    type distribs_container_type_ptr_type
        type(distribs_container_type), pointer :: p => NULL()
    end type distribs_container_type_ptr_type
    type(distribs_container_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    character(3), intent(inout), dimension(n0,2) :: elements
    real(4), intent(inout), dimension(n0) :: radii
    integer :: n0
    !f2py intent(hide), depend(elements) :: n0 = shape(elements,0)
    this_ptr = transfer(this, this_ptr)
    call this_ptr%p%get_bond_radii_staticmem(elements=elements, radii=radii)
end subroutine f90wrap_raffle__dc__get_bond_radii_staticmem__binding__dc_type
!###############################################################################


!###############################################################################
! initialise generalised distribution functions
!###############################################################################
subroutine f90wrap_raffle__dc__initialise_gdfs__binding__dc_type(this)
    use raffle__distribs_container, only: distribs_container_type
    implicit none
    
    type distribs_container_type_ptr_type
        type(distribs_container_type), pointer :: p => NULL()
    end type distribs_container_type_ptr_type
    type(distribs_container_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    call this_ptr%p%initialise_gdfs()
end subroutine f90wrap_raffle__dc__initialise_gdfs__binding__dc_type
!###############################################################################


!###############################################################################
! evolve the generalised distribution functions
!###############################################################################
subroutine f90wrap_raffle__dc__evolve__binding__dc_type(this) !, system)
    use raffle__distribs_container, only: distribs_container_type
    implicit none
    
    ! type distribs_type_ptr_type
    !     type(distribs_type), pointer :: p => NULL()
    ! end type distribs_type_ptr_type
    type distribs_container_type_ptr_type
        type(distribs_container_type), pointer :: p => NULL()
    end type distribs_container_type_ptr_type
    type(distribs_container_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    ! type(distribs_type_ptr_type) :: system_ptr
    ! integer, optional, intent(in), dimension(2) :: system
    this_ptr = transfer(this, this_ptr)
    ! if (present(system)) then
    !     system_ptr = transfer(system, system_ptr)
    ! else
    !     system_ptr%p => null()
    ! end if
    call this_ptr%p%evolve() !system=system_ptr%p)
end subroutine f90wrap_raffle__dc__evolve__binding__dc_type
!###############################################################################


!###############################################################################
! read and write distribution functions to file
!###############################################################################
subroutine f90wrap_raffle__dc__read_gdfs__binding__dc_type( &
     this, file &
)
    use raffle__distribs_container, only: distribs_container_type
    implicit none
    
    type distribs_container_type_ptr_type
        type(distribs_container_type), pointer :: p => NULL()
    end type distribs_container_type_ptr_type
    type(distribs_container_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    character*(*), intent(in) :: file
    this_ptr = transfer(this, this_ptr)
    call this_ptr%p%read_gdfs(file=file)
end subroutine f90wrap_raffle__dc__read_gdfs__binding__dc_type

subroutine f90wrap_raffle__dc__write_gdfs__binding__dc_type( &
     this, file &
)
    use raffle__distribs_container, only: distribs_container_type
    implicit none
    
    type distribs_container_type_ptr_type
        type(distribs_container_type), pointer :: p => NULL()
    end type distribs_container_type_ptr_type
    type(distribs_container_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    character*(*), intent(in) :: file
    this_ptr = transfer(this, this_ptr)
    call this_ptr%p%write_gdfs(file=file)
end subroutine f90wrap_raffle__dc__write_gdfs__binding__dc_type

subroutine f90wrap_raffle__dc__read_dfs__binding__dc_type( &
     this, file &
)
    use raffle__distribs_container, only: distribs_container_type
    implicit none
    
    type distribs_container_type_ptr_type
        type(distribs_container_type), pointer :: p => NULL()
    end type distribs_container_type_ptr_type
    type(distribs_container_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    character*(*), intent(in) :: file
    this_ptr = transfer(this, this_ptr)
    call this_ptr%p%read_dfs(file=file)
end subroutine f90wrap_raffle__dc__read_dfs__binding__dc_type

subroutine f90wrap_raffle__dc__write_dfs__binding__dc_type( &
     this, file &
)
    use raffle__distribs_container, only: distribs_container_type
    implicit none
    
    type distribs_container_type_ptr_type
        type(distribs_container_type), pointer :: p => NULL()
    end type distribs_container_type_ptr_type
    type(distribs_container_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    character*(*), intent(in) :: file
    this_ptr = transfer(this, this_ptr)
    call this_ptr%p%write_dfs(file=file)
end subroutine f90wrap_raffle__dc__write_dfs__binding__dc_type

subroutine f90wrap_raffle__dc__write_2body__binding__dc_type( &
     this, file &
)
    use raffle__distribs_container, only: distribs_container_type
    implicit none
    
    type distribs_container_type_ptr_type
        type(distribs_container_type), pointer :: p => NULL()
    end type distribs_container_type_ptr_type
    type(distribs_container_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    character*(*), intent(in) :: file
    this_ptr = transfer(this, this_ptr)
    call this_ptr%p%write_2body(file=file)
end subroutine f90wrap_raffle__dc__write_2body__binding__dc_type

subroutine f90wrap_raffle__dc__write_3body__binding__dc_type( &
     this, file &
)
    use raffle__distribs_container, only: distribs_container_type
    implicit none
    
    type distribs_container_type_ptr_type
        type(distribs_container_type), pointer :: p => NULL()
    end type distribs_container_type_ptr_type
    type(distribs_container_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    character*(*), intent(in) :: file
    this_ptr = transfer(this, this_ptr)
    call this_ptr%p%write_3body(file=file)
end subroutine f90wrap_raffle__dc__write_3body__binding__dc_type

subroutine f90wrap_raffle__dc__write_4body__binding__dc_type( &
     this, file &
)
    use raffle__distribs_container, only: distribs_container_type
    implicit none
    
    type distribs_container_type_ptr_type
        type(distribs_container_type), pointer :: p => NULL()
    end type distribs_container_type_ptr_type
    type(distribs_container_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    character*(*), intent(in) :: file
    this_ptr = transfer(this, this_ptr)
    call this_ptr%p%write_4body(file=file)
end subroutine f90wrap_raffle__dc__write_4body__binding__dc_type
!###############################################################################


!###############################################################################
! bin and pair index handling
!###############################################################################
subroutine f90wrap_raffle__dc__get_pair_index__binding__dc_type( &
     this, species1, ret_idx, species2 &
)
    use raffle__distribs_container, only: distribs_container_type
    implicit none
    
    type distribs_container_type_ptr_type
        type(distribs_container_type), pointer :: p => NULL()
    end type distribs_container_type_ptr_type
    type(distribs_container_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    character(3), intent(in) :: species1
    integer, intent(out) :: ret_idx
    character(3), intent(in) :: species2
    this_ptr = transfer(this, this_ptr)
    ret_idx = this_ptr%p%get_pair_index(species1=species1, species2=species2)
end subroutine f90wrap_raffle__dc__get_pair_index__binding__dc_type

subroutine f90wrap_raffle__dc__get_bin__binding__dc_type( &
     this, value, ret_bin, dim &
)
    use raffle__distribs_container, only: distribs_container_type
    implicit none
    
    type distribs_container_type_ptr_type
        type(distribs_container_type), pointer :: p => NULL()
    end type distribs_container_type_ptr_type
    type(distribs_container_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    real(4), intent(in) :: value
    integer, intent(out) :: ret_bin
    integer, intent(in) :: dim
    this_ptr = transfer(this, this_ptr)
    ret_bin = this_ptr%p%get_bin(value=value, dim=dim)
end subroutine f90wrap_raffle__dc__get_bin__binding__dc_type
!###############################################################################

! End of module raffle__distribs_container defined in file ../src/lib/mod_distribs_container.f90

