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

subroutine f90wrap_gvector_type__array__species(this, nd, dtype, dshape, dloc)
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
    if (allocated(this_ptr%p%species)) then
        dshape(1:2) = (/len(this_ptr%p%species(1)), shape(this_ptr%p%species)/)
        dloc = loc(this_ptr%p%species)
    else
        dloc = 0
    end if
end subroutine f90wrap_gvector_type__array__species

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

subroutine f90wrap_evolver__calculate__binding__gvector_type(this, lattice, basis, nbins, width, sigma, cutoff_min, &
    cutoff_max)
    use evolver, only: gvector_type
    use rw_geom, only: bas_type
    implicit none
    
    type bas_type_ptr_type
        type(bas_type), pointer :: p => NULL()
    end type bas_type_ptr_type
    type gvector_type_ptr_type
        type(gvector_type), pointer :: p => NULL()
    end type gvector_type_ptr_type
    type(gvector_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    real(4), dimension(3,3), intent(in) :: lattice
    type(bas_type_ptr_type) :: basis_ptr
    integer, intent(in), dimension(2) :: basis
    integer, dimension(3), intent(in), optional :: nbins
    real(4), dimension(3), intent(in), optional :: width
    real(4), dimension(3), intent(in), optional :: sigma
    real(4), dimension(3), intent(in), optional :: cutoff_min
    real(4), dimension(3), intent(in), optional :: cutoff_max
    this_ptr = transfer(this, this_ptr)
    basis_ptr = transfer(basis, basis_ptr)
    call this_ptr%p%calculate(lattice=lattice, basis=basis_ptr%p, nbins=nbins, width=width, sigma=sigma, &
        cutoff_min=cutoff_min, cutoff_max=cutoff_max)
end subroutine f90wrap_evolver__calculate__binding__gvector_type

subroutine f90wrap_gvector_container_type__get__best_system(this, f90wrap_best_system)
    use evolver, only: gvector_container_type
    implicit none
    type gvector_container_type_ptr_type
        type(gvector_container_type), pointer :: p => NULL()
    end type gvector_container_type_ptr_type
    integer, intent(in)   :: this(2)
    type(gvector_container_type_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_best_system
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_best_system = this_ptr%p%best_system
end subroutine f90wrap_gvector_container_type__get__best_system

subroutine f90wrap_gvector_container_type__set__best_system(this, f90wrap_best_system)
    use evolver, only: gvector_container_type
    implicit none
    type gvector_container_type_ptr_type
        type(gvector_container_type), pointer :: p => NULL()
    end type gvector_container_type_ptr_type
    integer, intent(in)   :: this(2)
    type(gvector_container_type_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_best_system
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%best_system = f90wrap_best_system
end subroutine f90wrap_gvector_container_type__set__best_system

subroutine f90wrap_gvector_container_type__get__best_energy(this, f90wrap_best_energy)
    use evolver, only: gvector_container_type
    implicit none
    type gvector_container_type_ptr_type
        type(gvector_container_type), pointer :: p => NULL()
    end type gvector_container_type_ptr_type
    integer, intent(in)   :: this(2)
    type(gvector_container_type_ptr_type) :: this_ptr
    real(4), intent(out) :: f90wrap_best_energy
    
    this_ptr = transfer(this, this_ptr)
    f90wrap_best_energy = this_ptr%p%best_energy
end subroutine f90wrap_gvector_container_type__get__best_energy

subroutine f90wrap_gvector_container_type__set__best_energy(this, f90wrap_best_energy)
    use evolver, only: gvector_container_type
    implicit none
    type gvector_container_type_ptr_type
        type(gvector_container_type), pointer :: p => NULL()
    end type gvector_container_type_ptr_type
    integer, intent(in)   :: this(2)
    type(gvector_container_type_ptr_type) :: this_ptr
    real(4), intent(in) :: f90wrap_best_energy
    
    this_ptr = transfer(this, this_ptr)
    this_ptr%p%best_energy = f90wrap_best_energy
end subroutine f90wrap_gvector_container_type__set__best_energy

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
        endif
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

subroutine f90wrap_evolver__add_basis__binding__gvector_container_type(this, lattice, basis)
    use rw_geom, only: bas_type
    use evolver, only: gvector_container_type
    implicit none
    
    type gvector_container_type_ptr_type
        type(gvector_container_type), pointer :: p => NULL()
    end type gvector_container_type_ptr_type
    type bas_type_ptr_type
        type(bas_type), pointer :: p => NULL()
    end type bas_type_ptr_type
    type(gvector_container_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    real(4), dimension(3,3), intent(in) :: lattice
    type(bas_type_ptr_type) :: basis_ptr
    integer, intent(in), dimension(2) :: basis
    this_ptr = transfer(this, this_ptr)
    basis_ptr = transfer(basis, basis_ptr)
    call this_ptr%p%add_basis(lattice=lattice, basis=basis_ptr%p)
end subroutine f90wrap_evolver__add_basis__binding__gvector_container_type

subroutine f90wrap_evolver__set_element_info__binding__gvector_containbcb0(this, element_file, element_list, n0)
    use evolver, only: gvector_container_type
    implicit none
    
    type gvector_container_type_ptr_type
        type(gvector_container_type), pointer :: p => NULL()
    end type gvector_container_type_ptr_type
    type(gvector_container_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    character*(*), intent(in), optional :: element_file
    character(3), intent(in), optional, dimension(n0) :: element_list
    integer :: n0
    !f2py intent(hide), depend(element_list) :: n0 = shape(element_list,0)
    this_ptr = transfer(this, this_ptr)
    call this_ptr%p%set_element_info(element_file=element_file, element_list=element_list)
end subroutine f90wrap_evolver__set_element_info__binding__gvector_containbcb0

subroutine f90wrap_evolver__set_bond_info__binding__gvector_container_type(this, bond_file)
    use evolver, only: gvector_container_type
    implicit none
    
    type gvector_container_type_ptr_type
        type(gvector_container_type), pointer :: p => NULL()
    end type gvector_container_type_ptr_type
    type(gvector_container_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    character*(*), intent(in), optional :: bond_file
    this_ptr = transfer(this, this_ptr)
    call this_ptr%p%set_bond_info(bond_file=bond_file)
end subroutine f90wrap_evolver__set_bond_info__binding__gvector_container_type

subroutine f90wrap_evolver__set_best_energy__binding__gvector_containe4680(this)
    use evolver, only: gvector_container_type
    implicit none
    
    type gvector_container_type_ptr_type
        type(gvector_container_type), pointer :: p => NULL()
    end type gvector_container_type_ptr_type
    type(gvector_container_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    call this_ptr%p%set_best_energy()
end subroutine f90wrap_evolver__set_best_energy__binding__gvector_containe4680

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

subroutine f90wrap_evolver__evolve__binding__gvector_container_type(this, system, deallocate_systems_after_evolve)
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
    logical, intent(in), optional :: deallocate_systems_after_evolve
    this_ptr = transfer(this, this_ptr)
    if (present(system)) then
        system_ptr = transfer(system, system_ptr)
    else
        system_ptr%p => null()
    end if
    call this_ptr%p%evolve(system=system_ptr%p, deallocate_systems_after_evolve=deallocate_systems_after_evolve)
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

