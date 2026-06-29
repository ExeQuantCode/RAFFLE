! Module raffle__graph_builder defined in file ../fortran/lib/mod_graph_builder.f90

subroutine f90wrap_topology_type__get__num_species( &
     this, ret_num_species &
)
    use raffle__graph_builder, only: topology_type
    implicit none
    type topology_type_ptr_type
        type(topology_type), pointer :: p => NULL()
    end type topology_type_ptr_type
    integer, intent(in)   :: this(2)
    type(topology_type_ptr_type) :: this_ptr
    integer, intent(out) :: ret_num_species

    this_ptr = transfer(this, this_ptr)
    if(.not.allocated(this_ptr%p%symbols)) then
        ret_num_species = 0
    else
        ret_num_species = size(this_ptr%p%symbols,1)
    end if
end subroutine f90wrap_topology_type__get__num_species

subroutine f90wrap_topology_type__get_symbols(this, symbols, n0)
    use raffle__graph_builder, only: topology_type
    implicit none
    type topology_type_ptr_type
        type(topology_type), pointer :: p => NULL()
    end type topology_type_ptr_type
    integer, intent(in) :: this(2)
    type(topology_type_ptr_type) :: this_ptr
    character(3), intent(inout), dimension(n0) :: symbols
    integer :: n0
    !f2py intent(hide), depend(symbols) :: n0 = shape(symbols,1)
    integer :: i

    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%symbols)) then
        do i = 1, min(size(symbols), size(this_ptr%p%symbols))
            symbols(i) = this_ptr%p%symbols(i)
        end do
    end if
end subroutine f90wrap_topology_type__get_symbols

subroutine f90wrap_topology_type__set_symbols(this, symbols, n0)
    use raffle__graph_builder, only: topology_type
    implicit none
    type topology_type_ptr_type
        type(topology_type), pointer :: p => NULL()
    end type topology_type_ptr_type
    integer, intent(in) :: this(2)
    type(topology_type_ptr_type) :: this_ptr
    character(3), intent(in), dimension(n0) :: symbols
    integer :: n0
    !f2py intent(hide), depend(symbols) :: n0 = shape(symbols,1)
    integer :: i

    this_ptr = transfer(this, this_ptr)
    if (.not.allocated(this_ptr%p%symbols)) then
        allocate(this_ptr%p%symbols(n0))
    end if
    do i = 1, min(size(symbols), size(this_ptr%p%symbols))
        this_ptr%p%symbols(i) = symbols(i)
    end do
end subroutine f90wrap_topology_type__set_symbols

subroutine f90wrap_topology_type__array__species_index(this, nd, dtype, dshape, dloc)
    use raffle__graph_builder, only: topology_type
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type topology_type_ptr_type
        type(topology_type), pointer :: p => NULL()
    end type topology_type_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(topology_type_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc

    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%species_index)) then
        dshape(1:1) = shape(this_ptr%p%species_index)
        dloc = loc(this_ptr%p%species_index)
    else
        dloc = 0
    end if
end subroutine f90wrap_topology_type__array__species_index

subroutine f90wrap_topology_type__array__atomic_numbers(this, nd, dtype, dshape, dloc)
    use raffle__graph_builder, only: topology_type
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type topology_type_ptr_type
        type(topology_type), pointer :: p => NULL()
    end type topology_type_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(topology_type_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc

    nd = 1
    dtype = 11
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%atomic_numbers)) then
        dshape(1:1) = shape(this_ptr%p%atomic_numbers)
        dloc = loc(this_ptr%p%atomic_numbers)
    else
        dloc = 0
    end if
end subroutine f90wrap_topology_type__array__atomic_numbers

subroutine f90wrap_topology_type__array__covalent_radii(this, nd, dtype, dshape, dloc)
    use raffle__graph_builder, only: topology_type
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type topology_type_ptr_type
        type(topology_type), pointer :: p => NULL()
    end type topology_type_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(topology_type_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc

    nd = 1
    dtype = 11
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%covalent_radii)) then
        dshape(1:1) = shape(this_ptr%p%covalent_radii)
        dloc = loc(this_ptr%p%covalent_radii)
    else
        dloc = 0
    end if
end subroutine f90wrap_topology_type__array__covalent_radii

subroutine f90wrap_topology_type__array__pair_image_shift(this, nd, dtype, dshape, dloc)
    use raffle__graph_builder, only: topology_type
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type topology_type_ptr_type
        type(topology_type), pointer :: p => NULL()
    end type topology_type_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(topology_type_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc

    nd = 2
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%pair_image_shift)) then
        dshape(1:2) = shape(this_ptr%p%pair_image_shift)
        dloc = loc(this_ptr%p%pair_image_shift)
    else
        dloc = 0
    end if
end subroutine f90wrap_topology_type__array__pair_image_shift

subroutine f90wrap_topology_type__array__pair_target_species_index(this, nd, dtype, dshape, dloc)
    use raffle__graph_builder, only: topology_type
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type topology_type_ptr_type
        type(topology_type), pointer :: p => NULL()
    end type topology_type_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(topology_type_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc

    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%pair_target_species_index)) then
        dshape(1:1) = shape(this_ptr%p%pair_target_species_index)
        dloc = loc(this_ptr%p%pair_target_species_index)
    else
        dloc = 0
    end if
end subroutine f90wrap_topology_type__array__pair_target_species_index

subroutine f90wrap_topology_type__array__pair_index(this, nd, dtype, dshape, dloc)
    use raffle__graph_builder, only: topology_type
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type topology_type_ptr_type
        type(topology_type), pointer :: p => NULL()
    end type topology_type_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(topology_type_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc

    nd = 2
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%pair_index)) then
        dshape(1:2) = shape(this_ptr%p%pair_index)
        dloc = loc(this_ptr%p%pair_index)
    else
        dloc = 0
    end if
end subroutine f90wrap_topology_type__array__pair_index

subroutine f90wrap_topology_type__array__pair_type_index(this, nd, dtype, dshape, dloc)
    use raffle__graph_builder, only: topology_type
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type topology_type_ptr_type
        type(topology_type), pointer :: p => NULL()
    end type topology_type_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(topology_type_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc

    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%pair_type_index)) then
        dshape(1:1) = shape(this_ptr%p%pair_type_index)
        dloc = loc(this_ptr%p%pair_type_index)
    else
        dloc = 0
    end if
end subroutine f90wrap_topology_type__array__pair_type_index

subroutine f90wrap_topology_type__array__pair_cutoff_weight_3body(this, nd, dtype, dshape, dloc)
    use raffle__graph_builder, only: topology_type
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type topology_type_ptr_type
        type(topology_type), pointer :: p => NULL()
    end type topology_type_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(topology_type_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc

    nd = 1
    dtype = 11
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%pair_cutoff_weight_3body)) then
        dshape(1:1) = shape(this_ptr%p%pair_cutoff_weight_3body)
        dloc = loc(this_ptr%p%pair_cutoff_weight_3body)
    else
        dloc = 0
    end if
end subroutine f90wrap_topology_type__array__pair_cutoff_weight_3body

subroutine f90wrap_topology_type__array__pair_cutoff_weight_4body(this, nd, dtype, dshape, dloc)
    use raffle__graph_builder, only: topology_type
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type topology_type_ptr_type
        type(topology_type), pointer :: p => NULL()
    end type topology_type_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(topology_type_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc

    nd = 1
    dtype = 11
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%pair_cutoff_weight_4body)) then
        dshape(1:1) = shape(this_ptr%p%pair_cutoff_weight_4body)
        dloc = loc(this_ptr%p%pair_cutoff_weight_4body)
    else
        dloc = 0
    end if
end subroutine f90wrap_topology_type__array__pair_cutoff_weight_4body

subroutine f90wrap_topology_type__array__triplet_species_index(this, nd, dtype, dshape, dloc)
    use raffle__graph_builder, only: topology_type
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type topology_type_ptr_type
        type(topology_type), pointer :: p => NULL()
    end type topology_type_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(topology_type_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc

    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%triplet_species_index)) then
        dshape(1:1) = shape(this_ptr%p%triplet_species_index)
        dloc = loc(this_ptr%p%triplet_species_index)
    else
        dloc = 0
    end if
end subroutine f90wrap_topology_type__array__triplet_species_index

subroutine f90wrap_topology_type__array__triplet_index(this, nd, dtype, dshape, dloc)
    use raffle__graph_builder, only: topology_type
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type topology_type_ptr_type
        type(topology_type), pointer :: p => NULL()
    end type topology_type_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(topology_type_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc

    nd = 2
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%triplet_index)) then
        dshape(1:2) = shape(this_ptr%p%triplet_index)
        dloc = loc(this_ptr%p%triplet_index)
    else
        dloc = 0
    end if
end subroutine f90wrap_topology_type__array__triplet_index

subroutine f90wrap_topology_type__array__triplet_pair_ids(this, nd, dtype, dshape, dloc)
    use raffle__graph_builder, only: topology_type
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type topology_type_ptr_type
        type(topology_type), pointer :: p => NULL()
    end type topology_type_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(topology_type_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc

    nd = 2
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%triplet_pair_ids)) then
        dshape(1:2) = shape(this_ptr%p%triplet_pair_ids)
        dloc = loc(this_ptr%p%triplet_pair_ids)
    else
        dloc = 0
    end if
end subroutine f90wrap_topology_type__array__triplet_pair_ids

subroutine f90wrap_topology_type__array__triplet_center_index(this, nd, dtype, dshape, dloc)
    use raffle__graph_builder, only: topology_type
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type topology_type_ptr_type
        type(topology_type), pointer :: p => NULL()
    end type topology_type_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(topology_type_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc

    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%triplet_center_index)) then
        dshape(1:1) = shape(this_ptr%p%triplet_center_index)
        dloc = loc(this_ptr%p%triplet_center_index)
    else
        dloc = 0
    end if
end subroutine f90wrap_topology_type__array__triplet_center_index

subroutine f90wrap_topology_type__array__quadruplet_pair_ids(this, nd, dtype, dshape, dloc)
    use raffle__graph_builder, only: topology_type
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type topology_type_ptr_type
        type(topology_type), pointer :: p => NULL()
    end type topology_type_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(topology_type_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc

    nd = 2
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%quadruplet_pair_ids)) then
        dshape(1:2) = shape(this_ptr%p%quadruplet_pair_ids)
        dloc = loc(this_ptr%p%quadruplet_pair_ids)
    else
        dloc = 0
    end if
end subroutine f90wrap_topology_type__array__quadruplet_pair_ids

subroutine f90wrap_topology_type__array__quadruplet_species_index(this, nd, dtype, dshape, dloc)
    use raffle__graph_builder, only: topology_type
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type topology_type_ptr_type
        type(topology_type), pointer :: p => NULL()
    end type topology_type_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(topology_type_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc

    nd = 1
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%quadruplet_species_index)) then
        dshape(1:1) = shape(this_ptr%p%quadruplet_species_index)
        dloc = loc(this_ptr%p%quadruplet_species_index)
    else
        dloc = 0
    end if
end subroutine f90wrap_topology_type__array__quadruplet_species_index

subroutine f90wrap_topology_type__get__num_atoms(this, f90wrap_num_atoms)
    use raffle__graph_builder, only: topology_type
    implicit none
    type topology_type_ptr_type
        type(topology_type), pointer :: p => NULL()
    end type topology_type_ptr_type
    integer, intent(in)   :: this(2)
    type(topology_type_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_num_atoms

    this_ptr = transfer(this, this_ptr)
    f90wrap_num_atoms = this_ptr%p%num_atoms
end subroutine f90wrap_topology_type__get__num_atoms

subroutine f90wrap_topology_type__set__num_atoms(this, f90wrap_num_atoms)
    use raffle__graph_builder, only: topology_type
    implicit none
    type topology_type_ptr_type
        type(topology_type), pointer :: p => NULL()
    end type topology_type_ptr_type
    integer, intent(in)   :: this(2)
    type(topology_type_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_num_atoms

    this_ptr = transfer(this, this_ptr)
    this_ptr%p%num_atoms = f90wrap_num_atoms
end subroutine f90wrap_topology_type__set__num_atoms

subroutine f90wrap_topology_type__get__num_pairs(this, f90wrap_num_pairs)
    use raffle__graph_builder, only: topology_type
    implicit none
    type topology_type_ptr_type
        type(topology_type), pointer :: p => NULL()
    end type topology_type_ptr_type
    integer, intent(in)   :: this(2)
    type(topology_type_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_num_pairs

    this_ptr = transfer(this, this_ptr)
    f90wrap_num_pairs = this_ptr%p%num_pairs
end subroutine f90wrap_topology_type__get__num_pairs

subroutine f90wrap_topology_type__set__num_pairs(this, f90wrap_num_pairs)
    use raffle__graph_builder, only: topology_type
    implicit none
    type topology_type_ptr_type
        type(topology_type), pointer :: p => NULL()
    end type topology_type_ptr_type
    integer, intent(in)   :: this(2)
    type(topology_type_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_num_pairs

    this_ptr = transfer(this, this_ptr)
    this_ptr%p%num_pairs = f90wrap_num_pairs
end subroutine f90wrap_topology_type__set__num_pairs

subroutine f90wrap_topology_type__get__num_triplets(this, f90wrap_num_triplets)
    use raffle__graph_builder, only: topology_type
    implicit none
    type topology_type_ptr_type
        type(topology_type), pointer :: p => NULL()
    end type topology_type_ptr_type
    integer, intent(in)   :: this(2)
    type(topology_type_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_num_triplets

    this_ptr = transfer(this, this_ptr)
    f90wrap_num_triplets = this_ptr%p%num_triplets
end subroutine f90wrap_topology_type__get__num_triplets

subroutine f90wrap_topology_type__set__num_triplets(this, f90wrap_num_triplets)
    use raffle__graph_builder, only: topology_type
    implicit none
    type topology_type_ptr_type
        type(topology_type), pointer :: p => NULL()
    end type topology_type_ptr_type
    integer, intent(in)   :: this(2)
    type(topology_type_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_num_triplets

    this_ptr = transfer(this, this_ptr)
    this_ptr%p%num_triplets = f90wrap_num_triplets
end subroutine f90wrap_topology_type__set__num_triplets

subroutine f90wrap_topology_type__get__num_quadruplets(this, f90wrap_num_quadruplets)
    use raffle__graph_builder, only: topology_type
    implicit none
    type topology_type_ptr_type
        type(topology_type), pointer :: p => NULL()
    end type topology_type_ptr_type
    integer, intent(in)   :: this(2)
    type(topology_type_ptr_type) :: this_ptr
    integer, intent(out) :: f90wrap_num_quadruplets

    this_ptr = transfer(this, this_ptr)
    f90wrap_num_quadruplets = this_ptr%p%num_quadruplets
end subroutine f90wrap_topology_type__get__num_quadruplets

subroutine f90wrap_topology_type__set__num_quadruplets(this, f90wrap_num_quadruplets)
    use raffle__graph_builder, only: topology_type
    implicit none
    type topology_type_ptr_type
        type(topology_type), pointer :: p => NULL()
    end type topology_type_ptr_type
    integer, intent(in)   :: this(2)
    type(topology_type_ptr_type) :: this_ptr
    integer, intent(in) :: f90wrap_num_quadruplets

    this_ptr = transfer(this, this_ptr)
    this_ptr%p%num_quadruplets = f90wrap_num_quadruplets
end subroutine f90wrap_topology_type__set__num_quadruplets

subroutine f90wrap_raffle__graph_builder__topology_type_initialise(this)
    use raffle__graph_builder, only: topology_type
    implicit none

    type topology_type_ptr_type
        type(topology_type), pointer :: p => NULL()
    end type topology_type_ptr_type
    type(topology_type_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_raffle__graph_builder__topology_type_initialise

subroutine f90wrap_raffle__graph_builder__topology_type_finalise(this)
    use raffle__graph_builder, only: topology_type
    implicit none

    type topology_type_ptr_type
        type(topology_type), pointer :: p => NULL()
    end type topology_type_ptr_type
    type(topology_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_raffle__graph_builder__topology_type_finalise

subroutine f90wrap_raffle__graph_builder__allocate_arrays__binding__toda88(this, num_atoms, num_pairs, &
    num_triplets, num_quadruplets)
    use raffle__graph_builder, only: topology_type
    implicit none

    type topology_type_ptr_type
        type(topology_type), pointer :: p => NULL()
    end type topology_type_ptr_type
    type(topology_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    integer, intent(in) :: num_atoms
    integer, intent(in) :: num_pairs
    integer, intent(in) :: num_triplets
    integer, intent(in) :: num_quadruplets
    this_ptr = transfer(this, this_ptr)
    call this_ptr%p%allocate_arrays(num_atoms=num_atoms, num_pairs=num_pairs, &
        num_triplets=num_triplets, num_quadruplets=num_quadruplets)
end subroutine f90wrap_raffle__graph_builder__allocate_arrays__binding__toda88

subroutine f90wrap_raffle__graph_builder__finalize__binding__topology_type(this)
    use raffle__graph_builder, only: topology_type
    implicit none

    type topology_type_ptr_type
        type(topology_type), pointer :: p => NULL()
    end type topology_type_ptr_type
    type(topology_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    call this_ptr%p%finalize()
end subroutine f90wrap_raffle__graph_builder__finalize__binding__topology_type

subroutine f90wrap_graph_tensors_type__array__global_features(this, nd, dtype, dshape, dloc)
    use raffle__graph_builder, only: graph_tensors_type
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type graph_tensors_type_ptr_type
        type(graph_tensors_type), pointer :: p => NULL()
    end type graph_tensors_type_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(graph_tensors_type_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc

    nd = 2
    dtype = 11
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%global_features)) then
        dshape(1:2) = shape(this_ptr%p%global_features)
        dloc = loc(this_ptr%p%global_features)
    else
        dloc = 0
    end if
end subroutine f90wrap_graph_tensors_type__array__global_features

subroutine f90wrap_graph_tensors_type__array__atom_node_features(this, nd, dtype, dshape, dloc)
    use raffle__graph_builder, only: graph_tensors_type
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type graph_tensors_type_ptr_type
        type(graph_tensors_type), pointer :: p => NULL()
    end type graph_tensors_type_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(graph_tensors_type_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc

    nd = 2
    dtype = 11
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%atom_node_features)) then
        dshape(1:2) = shape(this_ptr%p%atom_node_features)
        dloc = loc(this_ptr%p%atom_node_features)
    else
        dloc = 0
    end if
end subroutine f90wrap_graph_tensors_type__array__atom_node_features

subroutine f90wrap_graph_tensors_type__array__pair_node_features(this, nd, dtype, dshape, dloc)
    use raffle__graph_builder, only: graph_tensors_type
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type graph_tensors_type_ptr_type
        type(graph_tensors_type), pointer :: p => NULL()
    end type graph_tensors_type_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(graph_tensors_type_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc

    nd = 2
    dtype = 11
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%pair_node_features)) then
        dshape(1:2) = shape(this_ptr%p%pair_node_features)
        dloc = loc(this_ptr%p%pair_node_features)
    else
        dloc = 0
    end if
end subroutine f90wrap_graph_tensors_type__array__pair_node_features

subroutine f90wrap_graph_tensors_type__array__atom_edge_index(this, nd, dtype, dshape, dloc)
    use raffle__graph_builder, only: graph_tensors_type
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type graph_tensors_type_ptr_type
        type(graph_tensors_type), pointer :: p => NULL()
    end type graph_tensors_type_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(graph_tensors_type_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc

    nd = 2
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%atom_edge_index)) then
        dshape(1:2) = shape(this_ptr%p%atom_edge_index)
        dloc = loc(this_ptr%p%atom_edge_index)
    else
        dloc = 0
    end if
end subroutine f90wrap_graph_tensors_type__array__atom_edge_index

subroutine f90wrap_graph_tensors_type__array__pair_edge_index(this, nd, dtype, dshape, dloc)
    use raffle__graph_builder, only: graph_tensors_type
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type graph_tensors_type_ptr_type
        type(graph_tensors_type), pointer :: p => NULL()
    end type graph_tensors_type_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(graph_tensors_type_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc

    nd = 2
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%pair_edge_index)) then
        dshape(1:2) = shape(this_ptr%p%pair_edge_index)
        dloc = loc(this_ptr%p%pair_edge_index)
    else
        dloc = 0
    end if
end subroutine f90wrap_graph_tensors_type__array__pair_edge_index

subroutine f90wrap_graph_tensors_type__array__atom_edge_attr(this, nd, dtype, dshape, dloc)
    use raffle__graph_builder, only: graph_tensors_type
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type graph_tensors_type_ptr_type
        type(graph_tensors_type), pointer :: p => NULL()
    end type graph_tensors_type_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(graph_tensors_type_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc

    nd = 2
    dtype = 11
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%atom_edge_attr)) then
        dshape(1:2) = shape(this_ptr%p%atom_edge_attr)
        dloc = loc(this_ptr%p%atom_edge_attr)
    else
        dloc = 0
    end if
end subroutine f90wrap_graph_tensors_type__array__atom_edge_attr

subroutine f90wrap_graph_tensors_type__array__pair_edge_attr(this, nd, dtype, dshape, dloc)
    use raffle__graph_builder, only: graph_tensors_type
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type graph_tensors_type_ptr_type
        type(graph_tensors_type), pointer :: p => NULL()
    end type graph_tensors_type_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(graph_tensors_type_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc

    nd = 2
    dtype = 11
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%pair_edge_attr)) then
        dshape(1:2) = shape(this_ptr%p%pair_edge_attr)
        dloc = loc(this_ptr%p%pair_edge_attr)
    else
        dloc = 0
    end if
end subroutine f90wrap_graph_tensors_type__array__pair_edge_attr

subroutine f90wrap_graph_tensors_type__array__atom_edge_weight(this, nd, dtype, dshape, dloc)
    use raffle__graph_builder, only: graph_tensors_type
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type graph_tensors_type_ptr_type
        type(graph_tensors_type), pointer :: p => NULL()
    end type graph_tensors_type_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(graph_tensors_type_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc

    nd = 1
    dtype = 11
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%atom_edge_weight)) then
        dshape(1:1) = shape(this_ptr%p%atom_edge_weight)
        dloc = loc(this_ptr%p%atom_edge_weight)
    else
        dloc = 0
    end if
end subroutine f90wrap_graph_tensors_type__array__atom_edge_weight

subroutine f90wrap_graph_tensors_type__array__pair_edge_weight(this, nd, dtype, dshape, dloc)
    use raffle__graph_builder, only: graph_tensors_type
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type graph_tensors_type_ptr_type
        type(graph_tensors_type), pointer :: p => NULL()
    end type graph_tensors_type_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(graph_tensors_type_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc

    nd = 1
    dtype = 11
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%pair_edge_weight)) then
        dshape(1:1) = shape(this_ptr%p%pair_edge_weight)
        dloc = loc(this_ptr%p%pair_edge_weight)
    else
        dloc = 0
    end if
end subroutine f90wrap_graph_tensors_type__array__pair_edge_weight

subroutine f90wrap_graph_tensors_type__array__hyperedge_index(this, nd, dtype, dshape, dloc)
    use raffle__graph_builder, only: graph_tensors_type
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type graph_tensors_type_ptr_type
        type(graph_tensors_type), pointer :: p => NULL()
    end type graph_tensors_type_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(graph_tensors_type_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc

    nd = 2
    dtype = 5
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%hyperedge_index)) then
        dshape(1:2) = shape(this_ptr%p%hyperedge_index)
        dloc = loc(this_ptr%p%hyperedge_index)
    else
        dloc = 0
    end if
end subroutine f90wrap_graph_tensors_type__array__hyperedge_index

subroutine f90wrap_graph_tensors_type__array__hyperedge_weight(this, nd, dtype, dshape, dloc)
    use raffle__graph_builder, only: graph_tensors_type
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type graph_tensors_type_ptr_type
        type(graph_tensors_type), pointer :: p => NULL()
    end type graph_tensors_type_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(graph_tensors_type_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc

    nd = 1
    dtype = 11
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%hyperedge_weight)) then
        dshape(1:1) = shape(this_ptr%p%hyperedge_weight)
        dloc = loc(this_ptr%p%hyperedge_weight)
    else
        dloc = 0
    end if
end subroutine f90wrap_graph_tensors_type__array__hyperedge_weight

subroutine f90wrap_graph_tensors_type__array__hyperedge_attr(this, nd, dtype, dshape, dloc)
    use raffle__graph_builder, only: graph_tensors_type
    use, intrinsic :: iso_c_binding, only : c_int
    implicit none
    type graph_tensors_type_ptr_type
        type(graph_tensors_type), pointer :: p => NULL()
    end type graph_tensors_type_ptr_type
    integer(c_int), intent(in) :: this(2)
    type(graph_tensors_type_ptr_type) :: this_ptr
    integer(c_int), intent(out) :: nd
    integer(c_int), intent(out) :: dtype
    integer(c_int), dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc

    nd = 2
    dtype = 11
    this_ptr = transfer(this, this_ptr)
    if (allocated(this_ptr%p%hyperedge_attr)) then
        dshape(1:2) = shape(this_ptr%p%hyperedge_attr)
        dloc = loc(this_ptr%p%hyperedge_attr)
    else
        dloc = 0
    end if
end subroutine f90wrap_graph_tensors_type__array__hyperedge_attr

subroutine f90wrap_raffle__graph_builder__graph_tensors_type_initialise(this)
    use raffle__graph_builder, only: graph_tensors_type
    implicit none

    type graph_tensors_type_ptr_type
        type(graph_tensors_type), pointer :: p => NULL()
    end type graph_tensors_type_ptr_type
    type(graph_tensors_type_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    allocate(this_ptr%p)
    this = transfer(this_ptr, this)
end subroutine f90wrap_raffle__graph_builder__graph_tensors_type_initialise

subroutine f90wrap_raffle__graph_builder__graph_tensors_type_finalise(this)
    use raffle__graph_builder, only: graph_tensors_type
    implicit none

    type graph_tensors_type_ptr_type
        type(graph_tensors_type), pointer :: p => NULL()
    end type graph_tensors_type_ptr_type
    type(graph_tensors_type_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    deallocate(this_ptr%p)
end subroutine f90wrap_raffle__graph_builder__graph_tensors_type_finalise

! End of module raffle__graph_builder defined in file ../fortran/lib/mod_graph_builder.f90
