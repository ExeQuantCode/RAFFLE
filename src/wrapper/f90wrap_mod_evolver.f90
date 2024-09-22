! Module evolver defined in file ../src/lib/mod_evolver.f90

subroutine f90wrap_gvector_base_type__array__df_2body(this, nd, dtype, dshape, dloc)
    use evolver, only: gvector_base_type
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type gvector_base_type_ptr_type
        type(gvector_base_type), pointer :: p => NULL()
    end type gvector_base_type_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(gvector_base_type_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 11
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%df_2body)) then
        dshape(1:2) = shape(this_ptr%p%df_2body)
        dloc = loc(this_ptr%p%df_2body)
    else
        dloc = 0
    end if
end subroutine f90wrap_gvector_base_type__array__df_2body

subroutine f90wrap_gvector_base_type__array__df_3body(this, nd, dtype, dshape, dloc)
    use evolver, only: gvector_base_type
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type gvector_base_type_ptr_type
        type(gvector_base_type), pointer :: p => NULL()
    end type gvector_base_type_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(gvector_base_type_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 11
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%df_3body)) then
        dshape(1:2) = shape(this_ptr%p%df_3body)
        dloc = loc(this_ptr%p%df_3body)
    else
        dloc = 0
    end if
end subroutine f90wrap_gvector_base_type__array__df_3body

subroutine f90wrap_gvector_base_type__array__df_4body(this, nd, dtype, dshape, dloc)
    use evolver, only: gvector_base_type
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type gvector_base_type_ptr_type
        type(gvector_base_type), pointer :: p => NULL()
    end type gvector_base_type_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(gvector_base_type_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 11
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%df_4body)) then
        dshape(1:2) = shape(this_ptr%p%df_4body)
        dloc = loc(this_ptr%p%df_4body)
    else
        dloc = 0
    end if
end subroutine f90wrap_gvector_base_type__array__df_4body

subroutine f90wrap_evolver__gvector_base_type_initialise(this)
    use evolver, only: gvector_base_type
    implicit none
    
    type gvector_base_type_ptr_type
        type(gvector_base_type), pointer :: p => NULL()
    end type gvector_base_type_ptr_type
    type(gvector_base_type_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_evolver__gvector_base_type_initialise

subroutine f90wrap_evolver__gvector_base_type_finalise(this)
    use evolver, only: gvector_base_type
    implicit none
    
    type gvector_base_type_ptr_type
        type(gvector_base_type), pointer :: p => NULL()
    end type gvector_base_type_ptr_type
    type(gvector_base_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_evolver__gvector_base_type_finalise

subroutine f90wrap_gvector_type__get__num_atoms(this, f90wrap_num_atoms)
    use evolver, only: gvector_type
    implicit none
    type gvector_type_ptr_type
        type(gvector_type), pointer :: p => NULL()
    end type gvector_type_ptr_type
    integer, intent(in)   :: this(2)
    type(gvector_type_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_num_atoms
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_num_atoms = this_ptr%p%num_atoms
end subroutine f90wrap_gvector_type__get__num_atoms

subroutine f90wrap_gvector_type__set__num_atoms(this, f90wrap_num_atoms)
    use evolver, only: gvector_type
    implicit none
    type gvector_type_ptr_type
        type(gvector_type), pointer :: p => NULL()
    end type gvector_type_ptr_type
    integer, intent(in)   :: this(2)
    type(gvector_type_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_num_atoms
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%num_atoms = f90wrap_num_atoms
end subroutine f90wrap_gvector_type__set__num_atoms

subroutine f90wrap_gvector_type__get__energy(this, f90wrap_energy)
    use evolver, only: gvector_type
    implicit none
    type gvector_type_ptr_type
        type(gvector_type), pointer :: p => NULL()
    end type gvector_type_ptr_type
    integer, intent(in)   :: this(2)
    type(gvector_type_ptr_type) :: this_ptr
    real(4), intent(out) :: f90wrap_energy
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_energy = this_ptr%p%energy
end subroutine f90wrap_gvector_type__get__energy

subroutine f90wrap_gvector_type__set__energy(this, f90wrap_energy)
    use evolver, only: gvector_type
    implicit none
    type gvector_type_ptr_type
        type(gvector_type), pointer :: p => NULL()
    end type gvector_type_ptr_type
    integer, intent(in)   :: this(2)
    type(gvector_type_ptr_type) :: this_ptr
    real(4), intent(in) :: f90wrap_energy
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%energy = f90wrap_energy
end subroutine f90wrap_gvector_type__set__energy

subroutine f90wrap_gvector_type__array__stoichiometry(this, nd, dtype, dshape, dloc)
    use evolver, only: gvector_type
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type gvector_type_ptr_type
        type(gvector_type), pointer :: p => NULL()
    end type gvector_type_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(gvector_type_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%stoichiometry)) then
        dshape(1:1) = shape(this_ptr%p%stoichiometry)
        dloc = loc(this_ptr%p%stoichiometry)
    else
        dloc = 0
    end if
end subroutine f90wrap_gvector_type__array__stoichiometry

subroutine f90wrap_gvector_type__array__element_symbols(this, nd, dtype, dshape, dloc)
    use evolver, only: gvector_type
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type gvector_type_ptr_type
        type(gvector_type), pointer :: p => NULL()
    end type gvector_type_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(gvector_type_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 2
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%element_symbols)) then
        dshape(1:2) = (/len(this_ptr%p%element_symbols(1)), &
             shape(this_ptr%p%element_symbols)/)
        dloc = loc(this_ptr%p%element_symbols)
    else
        dloc = 0
    end if
end subroutine f90wrap_gvector_type__array__element_symbols

subroutine f90wrap_evolver__gvector_type_initialise(this)
    use evolver, only: gvector_type
    implicit none
    
    type gvector_type_ptr_type
        type(gvector_type), pointer :: p => NULL()
    end type gvector_type_ptr_type
    type(gvector_type_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_evolver__gvector_type_initialise

subroutine f90wrap_evolver__gvector_type_finalise(this)
    use evolver, only: gvector_type
    implicit none
    
    type gvector_type_ptr_type
        type(gvector_type), pointer :: p => NULL()
    end type gvector_type_ptr_type
    type(gvector_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_evolver__gvector_type_finalise

subroutine f90wrap_evolver__calculate__binding__gvector_type(this, basis, nbins, width, sigma, cutoff_min, cutoff_max, &
    radius_distance_tol)
    use evolver, only: gvector_type
    use rw_geom, only: basis_type
    implicit none
    
    type basis_type_ptr_type
        type(basis_type), pointer :: p => NULL()
    end type basis_type_ptr_type
    type gvector_type_ptr_type
        type(gvector_type), pointer :: p => NULL()
    end type gvector_type_ptr_type
    type(gvector_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    type(basis_type_ptr_type) :: basis_ptr
    integer, intent(in), dimension(2) :: basis
    integer, dimension(3), intent(in), optional :: nbins
    real(4), dimension(3), intent(in), optional :: width
    real(4), dimension(3), intent(in), optional :: sigma
    real(4), dimension(3), intent(in), optional :: cutoff_min
    real(4), dimension(3), intent(in), optional :: cutoff_max
    real(4), dimension(4), intent(in), optional :: radius_distance_tol
    this_ptr = transfer(this, this_ptr)
    basis_ptr = transfer(basis, basis_ptr)
    call this_ptr%p%calculate(basis=basis_ptr%p, nbins=nbins, width=width, sigma=sigma, cutoff_min=cutoff_min, &
        cutoff_max=cutoff_max, radius_distance_tol=radius_distance_tol)
end subroutine f90wrap_evolver__calculate__binding__gvector_type

subroutine f90wrap_gvector_container_type__get__num_evaluated(this, f90wrap_num_evaluated)
    use evolver, only: gvector_container_type
    implicit none
    type gvector_container_type_ptr_type
        type(gvector_container_type), pointer :: p => NULL()
    end type gvector_container_type_ptr_type
    integer, intent(in)   :: this(2)
    type(gvector_container_type_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_num_evaluated
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_num_evaluated = this_ptr%p%num_evaluated
end subroutine f90wrap_gvector_container_type__get__num_evaluated

subroutine f90wrap_gvector_container_type__set__num_evaluated(this, f90wrap_num_evaluated)
    use evolver, only: gvector_container_type
    implicit none
    type gvector_container_type_ptr_type
        type(gvector_container_type), pointer :: p => NULL()
    end type gvector_container_type_ptr_type
    integer, intent(in)   :: this(2)
    type(gvector_container_type_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_num_evaluated
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%num_evaluated = f90wrap_num_evaluated
end subroutine f90wrap_gvector_container_type__set__num_evaluated

subroutine f90wrap_gvector_container_type__get__num_evaluated_allocated(this, f90wrap_num_evaluated_allocated)
    use evolver, only: gvector_container_type
    implicit none
    type gvector_container_type_ptr_type
        type(gvector_container_type), pointer :: p => NULL()
    end type gvector_container_type_ptr_type
    integer, intent(in)   :: this(2)
    type(gvector_container_type_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_num_evaluated_allocated
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_num_evaluated_allocated = this_ptr%p%num_evaluated_allocated
end subroutine f90wrap_gvector_container_type__get__num_evaluated_allocated

subroutine f90wrap_gvector_container_type__set__num_evaluated_allocated(this, f90wrap_num_evaluated_allocated)
    use evolver, only: gvector_container_type
    implicit none
    type gvector_container_type_ptr_type
        type(gvector_container_type), pointer :: p => NULL()
    end type gvector_container_type_ptr_type
    integer, intent(in)   :: this(2)
    type(gvector_container_type_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_num_evaluated_allocated
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%num_evaluated_allocated = f90wrap_num_evaluated_allocated
end subroutine f90wrap_gvector_container_type__set__num_evaluated_allocated

! subroutine f90wrap_gvector_container_type__get__best_system(this, f90wrap_best_system)
!     use evolver, only: gvector_container_type
!     implicit none
!     type gvector_container_type_ptr_type
!         type(gvector_container_type), pointer :: p => NULL()
!     end type gvector_container_type_ptr_type
!     integer, intent(in)   :: this(2)
!     type(gvector_container_type_ptr_type) :: this_ptr
!     integer, intent(out) :: f90wrap_best_system
    
!     this_ptr = transfer(this, this_ptr)
!     f90wrap_best_system = this_ptr%p%best_system
! end subroutine f90wrap_gvector_container_type__get__best_system

! subroutine f90wrap_gvector_container_type__set__best_system(this, f90wrap_best_system)
!     use evolver, only: gvector_container_type
!     implicit none
!     type gvector_container_type_ptr_type
!         type(gvector_container_type), pointer :: p => NULL()
!     end type gvector_container_type_ptr_type
!     integer, intent(in)   :: this(2)
!     type(gvector_container_type_ptr_type) :: this_ptr
!     integer, intent(in) :: f90wrap_best_system
    
!     this_ptr = transfer(this, this_ptr)
!     this_ptr%p%best_system = f90wrap_best_system
! end subroutine f90wrap_gvector_container_type__set__best_system

! subroutine f90wrap_gvector_container_type__get__best_energy(this, f90wrap_best_energy)
!     use evolver, only: gvector_container_type
!     implicit none
!     type gvector_container_type_ptr_type
!         type(gvector_container_type), pointer :: p => NULL()
!     end type gvector_container_type_ptr_type
!     integer, intent(in)   :: this(2)
!     type(gvector_container_type_ptr_type) :: this_ptr
!     real(4), intent(out) :: f90wrap_best_energy
    
!     this_ptr = transfer(this, this_ptr)
!     f90wrap_best_energy = this_ptr%p%best_energy
! end subroutine f90wrap_gvector_container_type__get__best_energy

! subroutine f90wrap_gvector_container_type__set__best_energy(this, f90wrap_best_energy)
!     use evolver, only: gvector_container_type
!     implicit none
!     type gvector_container_type_ptr_type
!         type(gvector_container_type), pointer :: p => NULL()
!     end type gvector_container_type_ptr_type
!     integer, intent(in)   :: this(2)
!     type(gvector_container_type_ptr_type) :: this_ptr
!     real(4), intent(in) :: f90wrap_best_energy
    
!     this_ptr = transfer(this, this_ptr)
!     this_ptr%p%best_energy = f90wrap_best_energy
! end subroutine f90wrap_gvector_container_type__set__best_energy

subroutine f90wrap_gvector_container_type__get__kbt(this, f90wrap_kbt)
    use evolver, only: gvector_container_type
    implicit none
    type gvector_container_type_ptr_type
        type(gvector_container_type), pointer :: p => NULL()
    end type gvector_container_type_ptr_type
    integer, intent(in)   :: this(2)
    type(gvector_container_type_ptr_type) :: this_ptr
    real(4), intent(out) :: f90wrap_kbt
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_kbt = this_ptr%p%kbt
end subroutine f90wrap_gvector_container_type__get__kbt

subroutine f90wrap_gvector_container_type__set__kbt(this, f90wrap_kbt)
    use evolver, only: gvector_container_type
    implicit none
    type gvector_container_type_ptr_type
        type(gvector_container_type), pointer :: p => NULL()
    end type gvector_container_type_ptr_type
    integer, intent(in)   :: this(2)
    type(gvector_container_type_ptr_type) :: this_ptr
    real(4), intent(in) :: f90wrap_kbt
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%kbt = f90wrap_kbt
end subroutine f90wrap_gvector_container_type__set__kbt

subroutine f90wrap_gvector_container_type__array__nbins(this, nd, dtype, dshape, dloc)
    use evolver, only: gvector_container_type
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type gvector_container_type_ptr_type
        type(gvector_container_type), pointer :: p => NULL()
    end type gvector_container_type_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(gvector_container_type_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    dshape(1:1) = shape(this_ptr%p%nbins)
    dloc = loc(this_ptr%p%nbins)
end subroutine f90wrap_gvector_container_type__array__nbins

subroutine f90wrap_gvector_container_type__array__sigma(this, nd, dtype, dshape, dloc)
    use evolver, only: gvector_container_type
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type gvector_container_type_ptr_type
        type(gvector_container_type), pointer :: p => NULL()
    end type gvector_container_type_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(gvector_container_type_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 11
    this_ptr = transfer(this, this_ptr)
    dshape(1:1) = shape(this_ptr%p%sigma)
    dloc = loc(this_ptr%p%sigma)
end subroutine f90wrap_gvector_container_type__array__sigma

subroutine f90wrap_gvector_container_type__array__width(this, nd, dtype, dshape, dloc)
    use evolver, only: gvector_container_type
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type gvector_container_type_ptr_type
        type(gvector_container_type), pointer :: p => NULL()
    end type gvector_container_type_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(gvector_container_type_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 11
    this_ptr = transfer(this, this_ptr)
    dshape(1:1) = shape(this_ptr%p%width)
    dloc = loc(this_ptr%p%width)
end subroutine f90wrap_gvector_container_type__array__width

subroutine f90wrap_gvector_container_type__array__cutoff_min(this, nd, dtype, dshape, dloc)
    use evolver, only: gvector_container_type
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type gvector_container_type_ptr_type
        type(gvector_container_type), pointer :: p => NULL()
    end type gvector_container_type_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(gvector_container_type_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 11
    this_ptr = transfer(this, this_ptr)
    dshape(1:1) = shape(this_ptr%p%cutoff_min)
    dloc = loc(this_ptr%p%cutoff_min)
end subroutine f90wrap_gvector_container_type__array__cutoff_min

subroutine f90wrap_gvector_container_type__array__cutoff_max(this, nd, dtype, dshape, dloc)
    use evolver, only: gvector_container_type
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type gvector_container_type_ptr_type
        type(gvector_container_type), pointer :: p => NULL()
    end type gvector_container_type_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(gvector_container_type_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 11
    this_ptr = transfer(this, this_ptr)
    dshape(1:1) = shape(this_ptr%p%cutoff_max)
    dloc = loc(this_ptr%p%cutoff_max)
end subroutine f90wrap_gvector_container_type__array__cutoff_max

subroutine f90wrap_gvector_container_type__get__total(this, f90wrap_total)
    use evolver, only: gvector_container_type, gvector_base_type
    implicit none
    type gvector_container_type_ptr_type
        type(gvector_container_type), pointer :: p => NULL()
    end type gvector_container_type_ptr_type
    type gvector_base_type_ptr_type
        type(gvector_base_type), pointer :: p => NULL()
    end type gvector_base_type_ptr_type
    integer, intent(in)   :: this(2)
    type(gvector_container_type_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_total(2)
    type(gvector_base_type_ptr_type) :: total_ptr
    
    this_ptr = transfer(this, this_ptr)
    total_ptr%p => this_ptr%p%total
    f90wrap_total = transfer(total_ptr,f90wrap_total)
end subroutine f90wrap_gvector_container_type__get__total

subroutine f90wrap_gvector_container_type__set__total(this, f90wrap_total)
    use evolver, only: gvector_container_type, gvector_base_type
    implicit none
    type gvector_container_type_ptr_type
        type(gvector_container_type), pointer :: p => NULL()
    end type gvector_container_type_ptr_type
    type gvector_base_type_ptr_type
        type(gvector_base_type), pointer :: p => NULL()
    end type gvector_base_type_ptr_type
    integer, intent(in)   :: this(2)
    type(gvector_container_type_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_total(2)
    type(gvector_base_type_ptr_type) :: total_ptr
    
    this_ptr = transfer(this, this_ptr)
    total_ptr = transfer(f90wrap_total,total_ptr)
    this_ptr%p%total = total_ptr%p
end subroutine f90wrap_gvector_container_type__set__total

subroutine f90wrap_gvector_container_type__array_getitem__system(f90wrap_this, f90wrap_i, systemitem)
    
    use evolver, only: gvector_type, gvector_container_type
    implicit none
    
    type gvector_container_type_ptr_type
        type(gvector_container_type), pointer :: p => NULL()
    end type gvector_container_type_ptr_type
    type gvector_type_ptr_type
        type(gvector_type), pointer :: p => NULL()
    end type gvector_type_ptr_type
    integer, intent(in) :: f90wrap_this(2)
    type(gvector_container_type_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_i
    integer, intent(out) :: systemitem(2)
    type(gvector_type_ptr_type) :: system_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%system)) then
        if (f90wrap_i < 1 .or. f90wrap_i > size(this_ptr%p%system)) then
            call f90wrap_abort("array index out of range")
        else
            system_ptr%p => this_ptr%p%system(f90wrap_i)
            systemitem = transfer(system_ptr,systemitem)
        end if
    else
        call f90wrap_abort("derived type array not allocated")
    end if
end subroutine f90wrap_gvector_container_type__array_getitem__system

subroutine f90wrap_gvector_container_type__array_setitem__system(f90wrap_this, f90wrap_i, systemitem)
    
    use evolver, only: gvector_type, gvector_container_type
    implicit none
    
    type gvector_container_type_ptr_type
        type(gvector_container_type), pointer :: p => NULL()
    end type gvector_container_type_ptr_type
    type gvector_type_ptr_type
        type(gvector_type), pointer :: p => NULL()
    end type gvector_type_ptr_type
    integer, intent(in) :: f90wrap_this(2)
    type(gvector_container_type_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_i
    integer, intent(in) :: systemitem(2)
    type(gvector_type_ptr_type) :: system_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%system)) then
        if (f90wrap_i < 1 .or. f90wrap_i > size(this_ptr%p%system)) then
            call f90wrap_abort("array index out of range")
        else
            system_ptr = transfer(systemitem,system_ptr)
            this_ptr%p%system(f90wrap_i) = system_ptr%p
        endif
    else
        call f90wrap_abort("derived type array not allocated")
    end if
end subroutine f90wrap_gvector_container_type__array_setitem__system

subroutine f90wrap_gvector_container_type__array_len__system(f90wrap_this, f90wrap_n)
    
    use evolver, only: gvector_type, gvector_container_type
    implicit none
    
    type gvector_container_type_ptr_type
        type(gvector_container_type), pointer :: p => NULL()
    end type gvector_container_type_ptr_type
    type gvector_type_ptr_type
        type(gvector_type), pointer :: p => NULL()
    end type gvector_type_ptr_type
    integer, intent(out) :: f90wrap_n
    integer, intent(in) :: f90wrap_this(2)
    type(gvector_container_type_ptr_type) :: this_ptr
    
    this_ptr = transfer(f90wrap_this, this_ptr)
    if (allocated(this_ptr%p%system)) then
        f90wrap_n = size(this_ptr%p%system)
    else
        f90wrap_n = 0
    end if
end subroutine f90wrap_gvector_container_type__array_len__system

subroutine f90wrap_evolver__gvector_container_type_initialise(this)
    use evolver, only: gvector_container_type
    implicit none
    
    type gvector_container_type_ptr_type
        type(gvector_container_type), pointer :: p => NULL()
    end type gvector_container_type_ptr_type
    type(gvector_container_type_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_evolver__gvector_container_type_initialise

subroutine f90wrap_evolver__gvector_container_type_finalise(this)
    use evolver, only: gvector_container_type
    implicit none
    
    type gvector_container_type_ptr_type
        type(gvector_container_type), pointer :: p => NULL()
    end type gvector_container_type_ptr_type
    type(gvector_container_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_evolver__gvector_container_type_finalise

subroutine f90wrap_evolver__set_width__binding__gvector_container_type(this, width)
    use evolver, only: gvector_container_type
    implicit none
    
    type gvector_container_type_ptr_type
        type(gvector_container_type), pointer :: p => NULL()
    end type gvector_container_type_ptr_type
    type(gvector_container_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    real(4), dimension(3), intent(in) :: width
    this_ptr = transfer(this, this_ptr)
    call this_ptr%p%set_width(width=width)
end subroutine f90wrap_evolver__set_width__binding__gvector_container_type

subroutine f90wrap_evolver__set_sigma__binding__gvector_container_type(this, sigma)
    use evolver, only: gvector_container_type
    implicit none
    
    type gvector_container_type_ptr_type
        type(gvector_container_type), pointer :: p => NULL()
    end type gvector_container_type_ptr_type
    type(gvector_container_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    real(4), dimension(3), intent(in) :: sigma
    this_ptr = transfer(this, this_ptr)
    call this_ptr%p%set_sigma(sigma=sigma)
end subroutine f90wrap_evolver__set_sigma__binding__gvector_container_type

subroutine f90wrap_evolver__set_cutoff_min__binding__gvector_container7007(this, cutoff_min)
    use evolver, only: gvector_container_type
    implicit none
    
    type gvector_container_type_ptr_type
        type(gvector_container_type), pointer :: p => NULL()
    end type gvector_container_type_ptr_type
    type(gvector_container_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    real(4), dimension(3), intent(in) :: cutoff_min
    this_ptr = transfer(this, this_ptr)
    call this_ptr%p%set_cutoff_min(cutoff_min=cutoff_min)
end subroutine f90wrap_evolver__set_cutoff_min__binding__gvector_container7007

subroutine f90wrap_evolver__set_cutoff_max__binding__gvector_container047c(this, cutoff_max)
    use evolver, only: gvector_container_type
    implicit none
    
    type gvector_container_type_ptr_type
        type(gvector_container_type), pointer :: p => NULL()
    end type gvector_container_type_ptr_type
    type(gvector_container_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    real(4), dimension(3), intent(in) :: cutoff_max
    this_ptr = transfer(this, this_ptr)
    call this_ptr%p%set_cutoff_max(cutoff_max=cutoff_max)
end subroutine f90wrap_evolver__set_cutoff_max__binding__gvector_container047c

subroutine f90wrap_evolver__set_radius_distance_tol__binding__gvector_1dda(this, radius_distance_tol)
    use evolver, only: gvector_container_type
    implicit none
    
    type gvector_container_type_ptr_type
        type(gvector_container_type), pointer :: p => NULL()
    end type gvector_container_type_ptr_type
    type(gvector_container_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    real(4), dimension(4), intent(in) :: radius_distance_tol
    this_ptr = transfer(this, this_ptr)
    call this_ptr%p%set_radius_distance_tol(radius_distance_tol=radius_distance_tol)
end subroutine f90wrap_evolver__set_radius_distance_tol__binding__gvector_1dda

subroutine f90wrap_evolver__create__binding__gvector_container_type(this, basis_list, deallocate_systems)
    use rw_geom, only: basis_type
    use evolver, only: gvector_container_type
    implicit none
    
    type gvector_container_type_ptr_type
        type(gvector_container_type), pointer :: p => NULL()
    end type gvector_container_type_ptr_type

    type basis_type_xnum_array
        type(basis_type), dimension(:), allocatable :: items
    end type basis_type_xnum_array

    type basis_type_xnum_array_ptr_type
        type(basis_type_xnum_array), pointer :: p => NULL()
    end type basis_type_xnum_array_ptr_type
    type(gvector_container_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    type(basis_type_xnum_array_ptr_type) :: basis_list_ptr
    integer, intent(in), dimension(2) :: basis_list
    logical, intent(in), optional :: deallocate_systems

    this_ptr = transfer(this, this_ptr)
    basis_list_ptr = transfer(basis_list, basis_list_ptr)
    if(present(deallocate_systems)) then
        call this_ptr%p%create(basis_list=basis_list_ptr%p%items, deallocate_systems=deallocate_systems)
    else
        call this_ptr%p%create(basis_list=basis_list_ptr%p%items)
    end if
end subroutine f90wrap_evolver__create__binding__gvector_container_type

subroutine f90wrap_evolver__update__binding__gvector_container_type(this, basis_list, from_host, deallocate_systems)
    use rw_geom, only: basis_type
    use evolver, only: gvector_container_type
    implicit none
    
    type gvector_container_type_ptr_type
        type(gvector_container_type), pointer :: p => NULL()
    end type gvector_container_type_ptr_type

    type basis_type_xnum_array
        type(basis_type), dimension(:), allocatable :: items
    end type basis_type_xnum_array

    type basis_type_xnum_array_ptr_type
        type(basis_type_xnum_array), pointer :: p => NULL()
    end type basis_type_xnum_array_ptr_type
    type(gvector_container_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    type(basis_type_xnum_array_ptr_type) :: basis_list_ptr
    integer, intent(in), dimension(2) :: basis_list
    logical, intent(in) :: from_host
    logical, intent(in) :: deallocate_systems

    this_ptr = transfer(this, this_ptr)
    basis_list_ptr = transfer(basis_list, basis_list_ptr)
    call this_ptr%p%update(basis_list=basis_list_ptr%p%items, &
         from_host=from_host, &
         deallocate_systems=deallocate_systems)
end subroutine f90wrap_evolver__update__binding__gvector_container_type

subroutine f90wrap_evolver__deallocate_systems__binding__gvector_conta8f02(this)
    use evolver, only: gvector_container_type
    implicit none
    
    type gvector_container_type_ptr_type
        type(gvector_container_type), pointer :: p => NULL()
    end type gvector_container_type_ptr_type
    type(gvector_container_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    call this_ptr%p%deallocate_systems()
end subroutine f90wrap_evolver__deallocate_systems__binding__gvector_conta8f02

subroutine f90wrap_evolver__add_basis__binding__gvector_container_type(this, basis)
    use rw_geom, only: basis_type
    use evolver, only: gvector_container_type
    implicit none
    
    type gvector_container_type_ptr_type
        type(gvector_container_type), pointer :: p => NULL()
    end type gvector_container_type_ptr_type
    type basis_type_ptr_type
        type(basis_type), pointer :: p => NULL()
    end type basis_type_ptr_type
    type(gvector_container_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    type(basis_type_ptr_type) :: basis_ptr
    integer, intent(in), dimension(2) :: basis
    this_ptr = transfer(this, this_ptr)
    basis_ptr = transfer(basis, basis_ptr)
    call this_ptr%p%add_basis(basis=basis_ptr%p)
end subroutine f90wrap_evolver__add_basis__binding__gvector_container_type

subroutine f90wrap_evolver__get__num_elements(this, ret_num_elements)
    use evolver, only: gvector_container_type
    implicit none
    type gvector_container_type_ptr_type
        type(gvector_container_type), pointer :: p => NULL()
    end type gvector_container_type_ptr_type
    integer, intent(in)   :: this(2)
    type(gvector_container_type_ptr_type) :: this_ptr
    integer, intent(out) :: ret_num_elements
    
    this_ptr = transfer(this, this_ptr)
    if(.not.allocated(this_ptr%p%element_info)) then
        ret_num_elements = 0
    else
        ret_num_elements = size(this_ptr%p%element_info,1)
    end if
end subroutine f90wrap_evolver__get__num_elements

subroutine f90wrap_evolver__set_element_energies__binding__gvector_con0537(this, elements, energies, n0)
    use evolver, only: gvector_container_type
    implicit none
    
    type gvector_container_type_ptr_type
        type(gvector_container_type), pointer :: p => NULL()
    end type gvector_container_type_ptr_type
    type(gvector_container_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    character(3), dimension(n0), intent(in) :: elements
    real(4), dimension(n0), intent(in) :: energies
    integer :: n0
    !f2py intent(hide), depend(elements) :: n0 = shape(elements,0)
    this_ptr = transfer(this, this_ptr)
    call this_ptr%p%set_element_energies(elements=elements, energies=energies)
end subroutine f90wrap_evolver__set_element_energies__binding__gvector_con0537

subroutine f90wrap_evolver__get_element_energies_staticmem__binding__g4f53(this, elements, energies, n0)
    use evolver, only: gvector_container_type
    implicit none
    
    type gvector_container_type_ptr_type
        type(gvector_container_type), pointer :: p => NULL()
    end type gvector_container_type_ptr_type
    type(gvector_container_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    character(3), intent(inout), dimension(n0) :: elements
    real(4), intent(inout), dimension(n0) :: energies
    integer :: n0
    !f2py intent(hide), depend(elements) :: n0 = shape(elements,0)
    this_ptr = transfer(this, this_ptr)
    call this_ptr%p%get_element_energies_staticmem(elements=elements, energies=energies)
end subroutine f90wrap_evolver__get_element_energies_staticmem__binding__g4f53

subroutine f90wrap_evolver__set_bond_radius__binding__gvector_containe7df9(this, elements, radius)
    use evolver, only: gvector_container_type
    implicit none
    
    type gvector_container_type_ptr_type
        type(gvector_container_type), pointer :: p => NULL()
    end type gvector_container_type_ptr_type
    type(gvector_container_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    character(3), dimension(2), intent(in) :: elements
    real(4), intent(in) :: radius
    this_ptr = transfer(this, this_ptr)
    call this_ptr%p%set_bond_radius(elements=elements, radius=radius)
end subroutine f90wrap_evolver__set_bond_radius__binding__gvector_containe7df9

subroutine f90wrap_evolver__set_bond_radii__binding__gvector_container83c5(this, elements, radii, n0, n1, n2)
    use evolver, only: gvector_container_type
    implicit none
    
    type gvector_container_type_ptr_type
        type(gvector_container_type), pointer :: p => NULL()
    end type gvector_container_type_ptr_type
    type(gvector_container_type_ptr_type) :: this_ptr
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
end subroutine f90wrap_evolver__set_bond_radii__binding__gvector_container83c5

subroutine f90wrap_evolver__get_bond_radii_staticmem__binding__gvectord2e1(this, elements, radii, n0)
    use evolver, only: gvector_container_type
    implicit none
    
    type gvector_container_type_ptr_type
        type(gvector_container_type), pointer :: p => NULL()
    end type gvector_container_type_ptr_type
    type(gvector_container_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    character(3), intent(inout), dimension(n0,2) :: elements
    real(4), intent(inout), dimension(n0) :: radii
    integer :: n0
    !f2py intent(hide), depend(elements) :: n0 = shape(elements,0)
    this_ptr = transfer(this, this_ptr)
    call this_ptr%p%get_bond_radii_staticmem(elements=elements, radii=radii)
end subroutine f90wrap_evolver__get_bond_radii_staticmem__binding__gvectord2e1

! subroutine f90wrap_evolver__set_best_energy__binding__gvector_containe4680(this)
!     use evolver, only: gvector_container_type
!     implicit none
    
!     type gvector_container_type_ptr_type
!         type(gvector_container_type), pointer :: p => NULL()
!     end type gvector_container_type_ptr_type
!     type(gvector_container_type_ptr_type) :: this_ptr
!     integer, intent(in), dimension(2) :: this
!     this_ptr = transfer(this, this_ptr)
!     call this_ptr%p%set_best_energy()
! end subroutine f90wrap_evolver__set_best_energy__binding__gvector_containe4680

subroutine f90wrap_evolver__initialise_gvectors__binding__gvector_contc1f2(this)
    use evolver, only: gvector_container_type
    implicit none
    
    type gvector_container_type_ptr_type
        type(gvector_container_type), pointer :: p => NULL()
    end type gvector_container_type_ptr_type
    type(gvector_container_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    call this_ptr%p%initialise_gvectors()
end subroutine f90wrap_evolver__initialise_gvectors__binding__gvector_contc1f2

subroutine f90wrap_evolver__evolve__binding__gvector_container_type(this, system)
    use evolver, only: gvector_type, gvector_container_type
    implicit none
    
    type gvector_type_ptr_type
        type(gvector_type), pointer :: p => NULL()
    end type gvector_type_ptr_type
    type gvector_container_type_ptr_type
        type(gvector_container_type), pointer :: p => NULL()
    end type gvector_container_type_ptr_type
    type(gvector_container_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    type(gvector_type_ptr_type) :: system_ptr
    integer, optional, intent(in), dimension(2) :: system
    this_ptr = transfer(this, this_ptr)
    if (present(system)) then
        system_ptr = transfer(system, system_ptr)
    else
        system_ptr%p => null()
    end if
    call this_ptr%p%evolve(system=system_ptr%p)
end subroutine f90wrap_evolver__evolve__binding__gvector_container_type

subroutine f90wrap_evolver__write__binding__gvector_container_type(this, file)
    use evolver, only: gvector_container_type
    implicit none
    
    type gvector_container_type_ptr_type
        type(gvector_container_type), pointer :: p => NULL()
    end type gvector_container_type_ptr_type
    type(gvector_container_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    character*(*), intent(in) :: file
    this_ptr = transfer(this, this_ptr)
    call this_ptr%p%write(file=file)
end subroutine f90wrap_evolver__write__binding__gvector_container_type

subroutine f90wrap_evolver__read__binding__gvector_container_type(this, file)
    use evolver, only: gvector_container_type
    implicit none
    
    type gvector_container_type_ptr_type
        type(gvector_container_type), pointer :: p => NULL()
    end type gvector_container_type_ptr_type
    type(gvector_container_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    character*(*), intent(in) :: file
    this_ptr = transfer(this, this_ptr)
    call this_ptr%p%read(file=file)
end subroutine f90wrap_evolver__read__binding__gvector_container_type

subroutine f90wrap_evolver__write_2body__binding__gvector_container_type(this, file)
    use evolver, only: gvector_container_type
    implicit none
    
    type gvector_container_type_ptr_type
        type(gvector_container_type), pointer :: p => NULL()
    end type gvector_container_type_ptr_type
    type(gvector_container_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    character*(*), intent(in) :: file
    this_ptr = transfer(this, this_ptr)
    call this_ptr%p%write_2body(file=file)
end subroutine f90wrap_evolver__write_2body__binding__gvector_container_type

subroutine f90wrap_evolver__write_3body__binding__gvector_container_type(this, file)
    use evolver, only: gvector_container_type
    implicit none
    
    type gvector_container_type_ptr_type
        type(gvector_container_type), pointer :: p => NULL()
    end type gvector_container_type_ptr_type
    type(gvector_container_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    character*(*), intent(in) :: file
    this_ptr = transfer(this, this_ptr)
    call this_ptr%p%write_3body(file=file)
end subroutine f90wrap_evolver__write_3body__binding__gvector_container_type

subroutine f90wrap_evolver__write_4body__binding__gvector_container_type(this, file)
    use evolver, only: gvector_container_type
    implicit none
    
    type gvector_container_type_ptr_type
        type(gvector_container_type), pointer :: p => NULL()
    end type gvector_container_type_ptr_type
    type(gvector_container_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    character*(*), intent(in) :: file
    this_ptr = transfer(this, this_ptr)
    call this_ptr%p%write_4body(file=file)
end subroutine f90wrap_evolver__write_4body__binding__gvector_container_type

subroutine f90wrap_evolver__get_pair_index__binding__gvector_container4618(this, species1, ret_idx, species2)
    use evolver, only: gvector_container_type
    implicit none
    
    type gvector_container_type_ptr_type
        type(gvector_container_type), pointer :: p => NULL()
    end type gvector_container_type_ptr_type
    type(gvector_container_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    character(3), intent(in) :: species1
    integer, intent(out) :: ret_idx
    character(3), intent(in) :: species2
    this_ptr = transfer(this, this_ptr)
    ret_idx = this_ptr%p%get_pair_index(species1=species1, species2=species2)
end subroutine f90wrap_evolver__get_pair_index__binding__gvector_container4618

subroutine f90wrap_evolver__get_bin__binding__gvector_container_type(this, value, ret_bin, dim)
    use evolver, only: gvector_container_type
    implicit none
    
    type gvector_container_type_ptr_type
        type(gvector_container_type), pointer :: p => NULL()
    end type gvector_container_type_ptr_type
    type(gvector_container_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    real(4), intent(in) :: value
    integer, intent(out) :: ret_bin
    integer, intent(in) :: dim
    this_ptr = transfer(this, this_ptr)
    ret_bin = this_ptr%p%get_bin(value=value, dim=dim)
end subroutine f90wrap_evolver__get_bin__binding__gvector_container_type

! End of module evolver defined in file ../src/lib/mod_evolver.f90

