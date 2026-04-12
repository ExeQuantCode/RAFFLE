! Module raffle__gnn_fingerprint defined in file ../src/lib/mod_gnn_fingerprint.f90
! Wrapper subroutines for gnn_fingerprint_type

!###############################################################################
! gnn_fingerprint_type initialise/finalise
!###############################################################################
subroutine f90wrap_gnn_fingerprint_type_initialise(this)
	use raffle__gnn_fingerprint, only: gnn_fingerprint_type
	implicit none

	type gnn_fingerprint_type_ptr_type
		type(gnn_fingerprint_type), pointer :: p => NULL()
	end type gnn_fingerprint_type_ptr_type
	type(gnn_fingerprint_type_ptr_type) :: this_ptr
	integer, intent(out), dimension(2) :: this
	allocate(this_ptr%p)
	this = transfer(this_ptr, this)
end subroutine f90wrap_gnn_fingerprint_type_initialise

subroutine f90wrap_gnn_fingerprint_type_finalise(this)
	use raffle__gnn_fingerprint, only: gnn_fingerprint_type
	implicit none

	type gnn_fingerprint_type_ptr_type
		type(gnn_fingerprint_type), pointer :: p => NULL()
	end type gnn_fingerprint_type_ptr_type
	type(gnn_fingerprint_type_ptr_type) :: this_ptr
	integer, intent(in), dimension(2) :: this
	this_ptr = transfer(this, this_ptr)
	deallocate(this_ptr%p)
end subroutine f90wrap_gnn_fingerprint_type_finalise


!###############################################################################
! Property getters
!###############################################################################
subroutine f90wrap_gnn_fingerprint_type__get__is_initialised( &
	 this, f90wrap_is_initialised)
	use raffle__gnn_fingerprint, only: gnn_fingerprint_type
	implicit none
	type gnn_fingerprint_type_ptr_type
		type(gnn_fingerprint_type), pointer :: p => NULL()
	end type gnn_fingerprint_type_ptr_type
	integer, intent(in) :: this(2)
	type(gnn_fingerprint_type_ptr_type) :: this_ptr
	logical, intent(out) :: f90wrap_is_initialised
	this_ptr = transfer(this, this_ptr)
	f90wrap_is_initialised = this_ptr%p%is_initialised
end subroutine f90wrap_gnn_fingerprint_type__get__is_initialised

subroutine f90wrap_gnn_fingerprint_type__get__is_trained( &
	 this, f90wrap_is_trained)
	use raffle__gnn_fingerprint, only: gnn_fingerprint_type
	implicit none
	type gnn_fingerprint_type_ptr_type
		type(gnn_fingerprint_type), pointer :: p => NULL()
	end type gnn_fingerprint_type_ptr_type
	integer, intent(in) :: this(2)
	type(gnn_fingerprint_type_ptr_type) :: this_ptr
	logical, intent(out) :: f90wrap_is_trained
	this_ptr = transfer(this, this_ptr)
	f90wrap_is_trained = this_ptr%p%is_trained
end subroutine f90wrap_gnn_fingerprint_type__get__is_trained

subroutine f90wrap_gnn_fingerprint_type__get__fingerprint_dim( &
	 this, f90wrap_fingerprint_dim)
	use raffle__gnn_fingerprint, only: gnn_fingerprint_type
	implicit none
	type gnn_fingerprint_type_ptr_type
		type(gnn_fingerprint_type), pointer :: p => NULL()
	end type gnn_fingerprint_type_ptr_type
	integer, intent(in) :: this(2)
	type(gnn_fingerprint_type_ptr_type) :: this_ptr
	integer, intent(out) :: f90wrap_fingerprint_dim
	this_ptr = transfer(this, this_ptr)
	f90wrap_fingerprint_dim = this_ptr%p%fingerprint_dim
end subroutine f90wrap_gnn_fingerprint_type__get__fingerprint_dim

subroutine f90wrap_gnn_fingerprint_type__get__num_species( &
	 this, f90wrap_num_species)
	use raffle__gnn_fingerprint, only: gnn_fingerprint_type
	implicit none
	type gnn_fingerprint_type_ptr_type
		type(gnn_fingerprint_type), pointer :: p => NULL()
	end type gnn_fingerprint_type_ptr_type
	integer, intent(in) :: this(2)
	type(gnn_fingerprint_type_ptr_type) :: this_ptr
	integer, intent(out) :: f90wrap_num_species
	this_ptr = transfer(this, this_ptr)
	f90wrap_num_species = this_ptr%p%num_species
end subroutine f90wrap_gnn_fingerprint_type__get__num_species

subroutine f90wrap_gnn_fingerprint_type__get__num_vertex_features( &
	 this, f90wrap_num_vertex_features)
	use raffle__gnn_fingerprint, only: gnn_fingerprint_type
	implicit none
	type gnn_fingerprint_type_ptr_type
		type(gnn_fingerprint_type), pointer :: p => NULL()
	end type gnn_fingerprint_type_ptr_type
	integer, intent(in) :: this(2)
	type(gnn_fingerprint_type_ptr_type) :: this_ptr
	integer, intent(out) :: f90wrap_num_vertex_features
	this_ptr = transfer(this, this_ptr)
	f90wrap_num_vertex_features = this_ptr%p%num_vertex_features
end subroutine f90wrap_gnn_fingerprint_type__get__num_vertex_features

subroutine f90wrap_gnn_fingerprint_type__get__num_edge_features( &
	 this, f90wrap_num_edge_features)
	use raffle__gnn_fingerprint, only: gnn_fingerprint_type
	implicit none
	type gnn_fingerprint_type_ptr_type
		type(gnn_fingerprint_type), pointer :: p => NULL()
	end type gnn_fingerprint_type_ptr_type
	integer, intent(in) :: this(2)
	type(gnn_fingerprint_type_ptr_type) :: this_ptr
	integer, intent(out) :: f90wrap_num_edge_features
	this_ptr = transfer(this, this_ptr)
	f90wrap_num_edge_features = this_ptr%p%num_edge_features
end subroutine f90wrap_gnn_fingerprint_type__get__num_edge_features

subroutine f90wrap_gnn_fingerprint_type__get__gnn_output_dim( &
	 this, f90wrap_gnn_output_dim)
	use raffle__gnn_fingerprint, only: gnn_fingerprint_type
	implicit none
	type gnn_fingerprint_type_ptr_type
		type(gnn_fingerprint_type), pointer :: p => NULL()
	end type gnn_fingerprint_type_ptr_type
	integer, intent(in) :: this(2)
	type(gnn_fingerprint_type_ptr_type) :: this_ptr
	integer, intent(out) :: f90wrap_gnn_output_dim
	this_ptr = transfer(this, this_ptr)
	f90wrap_gnn_output_dim = this_ptr%p%gnn_output_dim
end subroutine f90wrap_gnn_fingerprint_type__get__gnn_output_dim

subroutine f90wrap_gnn_fingerprint_type__get__use_mlip_layer( &
	 this, f90wrap_use_mlip_layer)
	use raffle__gnn_fingerprint, only: gnn_fingerprint_type
	implicit none
	type gnn_fingerprint_type_ptr_type
		type(gnn_fingerprint_type), pointer :: p => NULL()
	end type gnn_fingerprint_type_ptr_type
	integer, intent(in) :: this(2)
	type(gnn_fingerprint_type_ptr_type) :: this_ptr
	logical, intent(out) :: f90wrap_use_mlip_layer
	this_ptr = transfer(this, this_ptr)
	f90wrap_use_mlip_layer = this_ptr%p%use_mlip_layer
end subroutine f90wrap_gnn_fingerprint_type__get__use_mlip_layer


!###############################################################################
! Initialise the network
!###############################################################################
subroutine f90wrap_gnn_fingerprint_type__initialise( &
	 this, species_list, n_species, num_time_steps, gnn_output_dim, &
	 max_degree, hidden_sizes, n_hidden, learning_rate, bond_cutoff, &
	 use_mlip_layer, n_rbf, kernel_hidden)
	use raffle__constants, only: real32
	use raffle__gnn_fingerprint, only: gnn_fingerprint_type
	implicit none
	type gnn_fingerprint_type_ptr_type
		type(gnn_fingerprint_type), pointer :: p => NULL()
	end type gnn_fingerprint_type_ptr_type
	integer, intent(in) :: this(2)
	type(gnn_fingerprint_type_ptr_type) :: this_ptr
	integer, intent(in) :: n_species
	character(len=3), dimension(n_species), intent(in) :: species_list
	integer, intent(in) :: num_time_steps
	integer, intent(in) :: gnn_output_dim
	integer, intent(in) :: max_degree
	integer, intent(in) :: n_hidden
	integer, dimension(n_hidden), intent(in) :: hidden_sizes
	real(real32), intent(in) :: learning_rate
	real(real32), intent(in) :: bond_cutoff
	logical, intent(in) :: use_mlip_layer
	integer, intent(in) :: n_rbf
	integer, intent(in) :: kernel_hidden

	this_ptr = transfer(this, this_ptr)
	call this_ptr%p%initialise( &
		 species_list = species_list, &
		 num_time_steps = num_time_steps, &
		 gnn_output_dim = gnn_output_dim, &
		 max_degree = max_degree, &
		 hidden_layer_sizes = hidden_sizes, &
		 learning_rate = learning_rate, &
		 bond_cutoff = bond_cutoff, &
		 use_mlip_layer = use_mlip_layer, &
		 n_rbf = n_rbf, &
		 kernel_hidden = kernel_hidden)
end subroutine f90wrap_gnn_fingerprint_type__initialise


!###############################################################################
! Compute fingerprint
!###############################################################################
subroutine f90wrap_gnn_fingerprint_type__compute_fingerprint( &
	 this, basis, fingerprint_out, fp_dim)
	use raffle__constants, only: real32
	use raffle__gnn_fingerprint, only: gnn_fingerprint_type
	use raffle__geom_rw, only: basis_type
	implicit none
	type gnn_fingerprint_type_ptr_type
		type(gnn_fingerprint_type), pointer :: p => NULL()
	end type gnn_fingerprint_type_ptr_type
	type basis_type_ptr_type
		type(basis_type), pointer :: p => NULL()
	end type basis_type_ptr_type
	integer, intent(in) :: this(2)
	integer, intent(in) :: basis(2)
	integer, intent(in) :: fp_dim
	real(real32), dimension(fp_dim), intent(out) :: fingerprint_out
	type(gnn_fingerprint_type_ptr_type) :: this_ptr
	type(basis_type_ptr_type) :: basis_ptr

	this_ptr = transfer(this, this_ptr)
	basis_ptr = transfer(basis, basis_ptr)
	call this_ptr%p%compute_fingerprint(basis_ptr%p, fingerprint_out)
end subroutine f90wrap_gnn_fingerprint_type__compute_fingerprint


!###############################################################################
! Predict (forward inference)
!###############################################################################
subroutine f90wrap_gnn_fingerprint_type__predict( &
	 this, basis, fingerprint_out, fp_dim)
	use raffle__constants, only: real32
	use raffle__gnn_fingerprint, only: gnn_fingerprint_type
	use raffle__geom_rw, only: basis_type
	implicit none
	type gnn_fingerprint_type_ptr_type
		type(gnn_fingerprint_type), pointer :: p => NULL()
	end type gnn_fingerprint_type_ptr_type
	type basis_type_ptr_type
		type(basis_type), pointer :: p => NULL()
	end type basis_type_ptr_type
	integer, intent(in) :: this(2)
	integer, intent(in) :: basis(2)
	integer, intent(in) :: fp_dim
	real(real32), dimension(fp_dim), intent(out) :: fingerprint_out
	type(gnn_fingerprint_type_ptr_type) :: this_ptr
	type(basis_type_ptr_type) :: basis_ptr

	this_ptr = transfer(this, this_ptr)
	basis_ptr = transfer(basis, basis_ptr)
	call this_ptr%p%predict(basis_ptr%p, fingerprint_out)
end subroutine f90wrap_gnn_fingerprint_type__predict


!###############################################################################
! Train
!###############################################################################
subroutine f90wrap_gnn_fingerprint_type__train( &
	 this, basis_handles, n_structures, num_epochs, batch_size, verbose)
	use raffle__gnn_fingerprint, only: gnn_fingerprint_type
	use raffle__geom_rw, only: basis_type
	implicit none
	type gnn_fingerprint_type_ptr_type
		type(gnn_fingerprint_type), pointer :: p => NULL()
	end type gnn_fingerprint_type_ptr_type
	type basis_type_ptr_type
		type(basis_type), pointer :: p => NULL()
	end type basis_type_ptr_type
	integer, intent(in) :: this(2)
	integer, intent(in) :: n_structures
	integer, dimension(2, n_structures), intent(in) :: basis_handles
	integer, intent(in) :: num_epochs, batch_size, verbose
	type(gnn_fingerprint_type_ptr_type) :: this_ptr
	type(basis_type_ptr_type) :: basis_ptr
	type(basis_type), dimension(:), allocatable :: structures
	integer :: i

	this_ptr = transfer(this, this_ptr)

	allocate(structures(n_structures))
	do i = 1, n_structures
	   basis_ptr = transfer(basis_handles(:, i), basis_ptr)
	   structures(i) = basis_ptr%p
	end do

	call this_ptr%p%train(structures, num_epochs, batch_size, verbose)

	deallocate(structures)
end subroutine f90wrap_gnn_fingerprint_type__train


!###############################################################################
! Inverse design
!###############################################################################
subroutine f90wrap_gnn_fingerprint_type__inverse_design( &
	 this, target_fp, fp_dim, basis, fixed_atoms, n_atoms, &
	 num_steps, step_size, verbose)
	use raffle__constants, only: real32
	use raffle__gnn_fingerprint, only: gnn_fingerprint_type
	use raffle__geom_rw, only: basis_type
	implicit none
	type gnn_fingerprint_type_ptr_type
		type(gnn_fingerprint_type), pointer :: p => NULL()
	end type gnn_fingerprint_type_ptr_type
	type basis_type_ptr_type
		type(basis_type), pointer :: p => NULL()
	end type basis_type_ptr_type
	integer, intent(in) :: this(2)
	integer, intent(in) :: fp_dim
	real(real32), dimension(fp_dim), intent(in) :: target_fp
	integer, intent(in) :: basis(2)
	integer, intent(in) :: n_atoms
	logical, dimension(n_atoms), intent(in) :: fixed_atoms
	integer, intent(in) :: num_steps
	real(real32), intent(in) :: step_size
	integer, intent(in) :: verbose
	type(gnn_fingerprint_type_ptr_type) :: this_ptr
	type(basis_type_ptr_type) :: basis_ptr

	this_ptr = transfer(this, this_ptr)
	basis_ptr = transfer(basis, basis_ptr)
	call this_ptr%p%inverse_design( &
		 target_fingerprint = target_fp, &
		 basis = basis_ptr%p, &
		 fixed_atoms = fixed_atoms, &
		 num_steps = num_steps, &
		 step_size = step_size, &
		 verbose = verbose)
end subroutine f90wrap_gnn_fingerprint_type__inverse_design
