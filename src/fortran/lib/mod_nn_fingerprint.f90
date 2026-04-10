module raffle__nn_fingerprint
  !! Module for neural network-based descriptor fingerprint prediction.
  !!
  !! This module implements a neural network framework using the ATHENA library
  !! that:
  !!   1. Learns a RAFFLE descriptor fingerprint from atomic structures
  !!   2. Supports forward inference mode for predicting descriptors
  !!   3. Enables inverse design (generating structures from target descriptors)
  !!   4. Allows partial atomic optimisation via boolean atom masks
  use raffle__constants, only: real32, pi
  use raffle__io_utils, only: stop_program, print_warning
  use raffle__geom_rw, only: basis_type
  use raffle__distribs, only: distribs_base_type, distribs_type
  use raffle__distribs_container, only: distribs_container_type
  use athena, only: &
       network_type, &
       full_layer_type, &
       adam_optimiser_type, &
       base_optimiser_type
  implicit none


  private

  public :: nn_fingerprint_type


  !-----------------------------------------------------------------------------
  ! Maximum species for encoding
  !-----------------------------------------------------------------------------
  integer, parameter :: MAX_SPECIES = 10
  !! Maximum number of distinct species supported for one-hot encoding.

  !-----------------------------------------------------------------------------
  ! Neural network fingerprint type
  !-----------------------------------------------------------------------------
  type :: nn_fingerprint_type
     !! Neural network for learning and predicting RAFFLE descriptor
     !! fingerprints from atomic structures, and for inverse design.
     logical :: is_initialised = .false.
     !! Whether the network has been initialised.
     logical :: is_trained = .false.
     !! Whether the network has been trained.

     integer :: num_species = 0
     !! Number of distinct species in the training set.
     character(len=3), dimension(:), allocatable :: species_list
     !! List of species symbols.

     integer :: max_atoms = 0
     !! Maximum number of atoms per structure (for padding).
     integer :: input_dim = 0
     !! Dimension of the flattened input vector.
     integer :: fingerprint_dim = 0
     !! Dimension of the output fingerprint vector.

     integer, dimension(3) :: nbins = [-1, -1, -1]
     !! Number of bins for 2-body, 3-body, 4-body distributions.
     integer :: num_pairs = 0
     !! Number of element pairs for 2-body.

     real(real32), dimension(3) :: width = &
          [0.025_real32, pi/64._real32, pi/64._real32]
     !! Bin widths for 2/3/4-body distributions.
     real(real32), dimension(3) :: sigma = &
          [0.1_real32, 0.1_real32, 0.1_real32]
     !! Gaussian widths for 2/3/4-body distributions.
     real(real32), dimension(3) :: cutoff_min = &
          [0.5_real32, 0._real32, 0._real32]
     !! Minimum cutoffs for 2/3/4-body distributions.
     real(real32), dimension(3) :: cutoff_max = &
          [6._real32, pi, pi]
     !! Maximum cutoffs for 2/3/4-body distributions.
     real(real32), dimension(4) :: radius_distance_tol = &
          [1.5_real32, 2.5_real32, 3._real32, 6._real32]
     !! Radius distance tolerance factors.

     type(network_type) :: network
     !! ATHENA neural network for fingerprint prediction.

   contains
     procedure, pass(this) :: initialise
     !! Initialise the neural network architecture.
     procedure, pass(this) :: train
     !! Train the neural network on a set of structures and their fingerprints.
     procedure, pass(this) :: predict
     !! Forward inference: predict fingerprint from atomic structure.
     procedure, pass(this) :: inverse_design
     !! Inverse design: generate structure from target descriptor.
     procedure, pass(this) :: basis_to_input
     !! Convert a basis_type to a flat input vector for the network.
     procedure, pass(this) :: compute_fingerprint
     !! Compute the RAFFLE descriptor fingerprint for a given structure.
     procedure, pass(this) :: fingerprint_to_distribs
     !! Convert a flat fingerprint vector back to distribs_base_type.
     procedure, pass(this) :: set_distribution_params
     !! Set distribution function parameters from a container.
  end type nn_fingerprint_type


contains


!###############################################################################
  subroutine set_distribution_params(this, container)
    !! Set distribution parameters from an existing distribs_container_type.
    implicit none

    ! Arguments
    class(nn_fingerprint_type), intent(inout) :: this
    !! Parent. Instance of nn_fingerprint_type.
    type(distribs_container_type), intent(in) :: container
    !! Distribution container to copy parameters from.

    this%nbins = container%nbins
    this%width = container%width
    this%sigma = container%sigma
    this%cutoff_min = container%cutoff_min
    this%cutoff_max = container%cutoff_max
    this%radius_distance_tol = container%radius_distance_tol

  end subroutine set_distribution_params
!###############################################################################


!###############################################################################
  subroutine initialise(this, species_list, max_atoms, &
       hidden_layer_sizes, learning_rate)
    !! Initialise the neural network for fingerprint prediction.
    !!
    !! Sets up the network architecture with dense layers.
    !! Input dimension = max_atoms * (3 + num_species) [positions + one-hot]
    !! Output dimension = fingerprint_dim (flattened 2/3/4-body distributions)
    implicit none

    ! Arguments
    class(nn_fingerprint_type), intent(inout) :: this
    !! Parent. Instance of nn_fingerprint_type.
    character(len=3), dimension(:), intent(in) :: species_list
    !! List of element species to handle.
    integer, intent(in) :: max_atoms
    !! Maximum number of atoms per structure (for input padding).
    integer, dimension(:), intent(in), optional :: hidden_layer_sizes
    !! Sizes of hidden layers. Default: [128, 64].
    real(real32), intent(in), optional :: learning_rate
    !! Learning rate for the optimizer. Default: 0.001.

    ! Local variables
    integer :: i, num_hidden
    integer, dimension(:), allocatable :: h_sizes
    real(real32) :: lr
    integer :: features_per_atom

    ! Set species
    this%num_species = size(species_list)
    if (allocated(this%species_list)) deallocate(this%species_list)
    allocate(this%species_list(this%num_species))
    this%species_list = species_list

    ! Set max atoms
    this%max_atoms = max_atoms

    ! Features per atom: 3 coordinates + one-hot species encoding
    features_per_atom = 3 + this%num_species
    this%input_dim = max_atoms * features_per_atom

    ! Calculate number of element pairs for 2-body
    this%num_pairs = nint( &
         gamma(real(this%num_species + 2, real32)) / &
         ( gamma(real(this%num_species, real32)) * gamma(3._real32) ) &
    )

    ! Set default nbins if not yet set
    if (this%nbins(1) .lt. 0) then
       this%nbins(1) = 1 + nint( &
            (this%cutoff_max(1) - this%cutoff_min(1)) / this%width(1) )
       this%nbins(2) = 1 + nint( &
            (this%cutoff_max(2) - this%cutoff_min(2)) / this%width(2) )
       this%nbins(3) = 1 + nint( &
            (this%cutoff_max(3) - this%cutoff_min(3)) / this%width(3) )
    end if

    ! Fingerprint dimension: 2body + 3body + 4body
    this%fingerprint_dim = &
         this%nbins(1) * this%num_pairs + &
         this%nbins(2) * this%num_species + &
         this%nbins(3) * this%num_species

    ! Hidden layer sizes
    if (present(hidden_layer_sizes)) then
       num_hidden = size(hidden_layer_sizes)
       allocate(h_sizes(num_hidden))
       h_sizes = hidden_layer_sizes
    else
       num_hidden = 2
       allocate(h_sizes(2))
       h_sizes = [128, 64]
    end if

    ! Learning rate
    lr = 0.001_real32
    if (present(learning_rate)) lr = learning_rate

    ! Build network: input -> hidden layers -> output
    ! First hidden layer (takes input_dim)
    call this%network%add(full_layer_type( &
         num_inputs = this%input_dim, &
         num_outputs = h_sizes(1), &
         activation = "relu" &
    ))

    ! Intermediate hidden layers
    do i = 2, num_hidden
       call this%network%add(full_layer_type( &
            num_outputs = h_sizes(i), &
            activation = "relu" &
       ))
    end do

    ! Output layer (linear activation for regression)
    call this%network%add(full_layer_type( &
         num_outputs = this%fingerprint_dim, &
         activation = "linear" &
    ))

    ! Compile network
    call this%network%compile( &
         optimiser = adam_optimiser_type(learning_rate = lr), &
         loss_method = "mse", &
         metrics = ["loss"], &
         verbose = 0 &
    )

    this%is_initialised = .true.

    if (allocated(h_sizes)) deallocate(h_sizes)

  end subroutine initialise
!###############################################################################


!###############################################################################
  subroutine basis_to_input(this, basis, input_vector)
    !! Convert a basis_type atomic structure to a flat input vector.
    !!
    !! The input vector is structured as:
    !!   [atom_1_x, atom_1_y, atom_1_z, one_hot_species_1, ...,
    !!    atom_2_x, atom_2_y, atom_2_z, one_hot_species_2, ...,
    !!    ..., zero_padding]
    !! with total length = max_atoms * (3 + num_species)
    implicit none

    ! Arguments
    class(nn_fingerprint_type), intent(in) :: this
    !! Parent. Instance of nn_fingerprint_type.
    type(basis_type), intent(in) :: basis
    !! Atomic structure to convert.
    real(real32), dimension(:), intent(out) :: input_vector
    !! Flat input vector of length input_dim.

    ! Local variables
    integer :: is, ia, atom_idx, species_idx, offset
    integer :: features_per_atom

    features_per_atom = 3 + this%num_species
    input_vector = 0._real32

    atom_idx = 0
    do is = 1, basis%nspec
       ! Find species index in our species list
       species_idx = 0
       do ia = 1, this%num_species
          if (this%species_list(ia) == basis%spec(is)%name) then
             species_idx = ia
             exit
          end if
       end do
       if (species_idx == 0) cycle

       do ia = 1, basis%spec(is)%num
          atom_idx = atom_idx + 1
          if (atom_idx > this%max_atoms) exit

          offset = (atom_idx - 1) * features_per_atom

          ! Fractional coordinates (or Cartesian if lcart)
          input_vector(offset + 1) = basis%spec(is)%atom(ia, 1)
          input_vector(offset + 2) = basis%spec(is)%atom(ia, 2)
          input_vector(offset + 3) = basis%spec(is)%atom(ia, 3)

          ! One-hot species encoding
          input_vector(offset + 3 + species_idx) = 1._real32
       end do
    end do

  end subroutine basis_to_input
!###############################################################################


!###############################################################################
  subroutine compute_fingerprint(this, basis, fingerprint)
    !! Compute the RAFFLE descriptor fingerprint for a given structure.
    !!
    !! Uses the distribs_type%calculate() to compute 2/3/4-body distributions
    !! and flattens them into a single vector.
    implicit none

    ! Arguments
    class(nn_fingerprint_type), intent(in) :: this
    !! Parent. Instance of nn_fingerprint_type.
    type(basis_type), intent(in) :: basis
    !! Atomic structure.
    real(real32), dimension(:), intent(out) :: fingerprint
    !! Output fingerprint vector of length fingerprint_dim.

    ! Local variables
    type(distribs_type) :: distribs
    integer :: offset, i, j

    ! Calculate distributions
    call distribs%calculate(basis, &
         nbins = this%nbins, &
         width = this%width, &
         sigma = this%sigma, &
         cutoff_min = this%cutoff_min, &
         cutoff_max = this%cutoff_max, &
         radius_distance_tol = this%radius_distance_tol)

    fingerprint = 0._real32
    offset = 0

    ! Flatten 2-body
    do j = 1, size(distribs%df_2body, 2)
       do i = 1, size(distribs%df_2body, 1)
          offset = offset + 1
          fingerprint(offset) = distribs%df_2body(i, j)
       end do
    end do

    ! Flatten 3-body
    do j = 1, size(distribs%df_3body, 2)
       do i = 1, size(distribs%df_3body, 1)
          offset = offset + 1
          fingerprint(offset) = distribs%df_3body(i, j)
       end do
    end do

    ! Flatten 4-body
    do j = 1, size(distribs%df_4body, 2)
       do i = 1, size(distribs%df_4body, 1)
          offset = offset + 1
          fingerprint(offset) = distribs%df_4body(i, j)
       end do
    end do

  end subroutine compute_fingerprint
!###############################################################################


!###############################################################################
  subroutine fingerprint_to_distribs(this, fp_vector, distribs)
    !! Convert a flat fingerprint vector back to distribs_base_type.
    implicit none

    ! Arguments
    class(nn_fingerprint_type), intent(in) :: this
    !! Parent. Instance of nn_fingerprint_type.
    real(real32), dimension(:), intent(in) :: fp_vector
    !! Fingerprint vector.
    type(distribs_base_type), intent(out) :: distribs
    !! Output distribution functions.

    ! Local variables
    integer :: offset, i, j

    allocate(distribs%df_2body(this%nbins(1), this%num_pairs))
    allocate(distribs%df_3body(this%nbins(2), this%num_species))
    allocate(distribs%df_4body(this%nbins(3), this%num_species))

    offset = 0
    do j = 1, this%num_pairs
       do i = 1, this%nbins(1)
          offset = offset + 1
          distribs%df_2body(i, j) = fp_vector(offset)
       end do
    end do
    do j = 1, this%num_species
       do i = 1, this%nbins(2)
          offset = offset + 1
          distribs%df_3body(i, j) = fp_vector(offset)
       end do
    end do
    do j = 1, this%num_species
       do i = 1, this%nbins(3)
          offset = offset + 1
          distribs%df_4body(i, j) = fp_vector(offset)
       end do
    end do

  end subroutine fingerprint_to_distribs
!###############################################################################


!###############################################################################
  subroutine train(this, structures, num_epochs, batch_size, verbose)
    !! Train the neural network on a set of structures.
    !!
    !! For each structure, computes the input vector (atomic features) and
    !! the target fingerprint (RAFFLE descriptor), then trains the network.
    implicit none

    ! Arguments
    class(nn_fingerprint_type), intent(inout) :: this
    !! Parent. Instance of nn_fingerprint_type.
    type(basis_type), dimension(:), intent(in) :: structures
    !! Training structures.
    integer, intent(in), optional :: num_epochs
    !! Number of training epochs. Default: 100.
    integer, intent(in), optional :: batch_size
    !! Batch size for training. Default: 1.
    integer, intent(in), optional :: verbose
    !! Verbosity level. Default: 0.

    ! Local variables
    integer :: n, num_strucs, epochs, bs, verb
    real(real32), dimension(:,:), allocatable :: input_data, target_data
    real(real32), dimension(:), allocatable :: input_vec, fp_vec

    if (.not. this%is_initialised) then
       call stop_program("nn_fingerprint: network not initialised")
       return
    end if

    num_strucs = size(structures)
    epochs = 100
    bs = min(num_strucs, 32)
    verb = 0
    if (present(num_epochs)) epochs = num_epochs
    if (present(batch_size)) bs = batch_size
    if (present(verbose)) verb = verbose

    ! Prepare training data
    allocate(input_data(this%input_dim, num_strucs))
    allocate(target_data(this%fingerprint_dim, num_strucs))
    allocate(input_vec(this%input_dim))
    allocate(fp_vec(this%fingerprint_dim))

    do n = 1, num_strucs
       call this%basis_to_input(structures(n), input_vec)
       input_data(:, n) = input_vec

       call this%compute_fingerprint(structures(n), fp_vec)
       target_data(:, n) = fp_vec
    end do

    ! Train with ATHENA
    call this%network%train( &
         input = input_data, &
         output = target_data, &
         num_epochs = epochs, &
         batch_size = bs, &
         verbose = verb &
    )

    this%is_trained = .true.

    deallocate(input_data, target_data, input_vec, fp_vec)

  end subroutine train
!###############################################################################


!###############################################################################
  subroutine predict(this, basis, fingerprint)
    !! Forward inference: predict a RAFFLE descriptor fingerprint from
    !! an atomic structure using the trained neural network.
    implicit none

    ! Arguments
    class(nn_fingerprint_type), intent(inout) :: this
    !! Parent. Instance of nn_fingerprint_type.
    type(basis_type), intent(in) :: basis
    !! Atomic structure input.
    real(real32), dimension(:), intent(out) :: fingerprint
    !! Predicted fingerprint vector.

    ! Local variables
    real(real32), dimension(:,:), allocatable :: input_2d, output_2d
    real(real32), dimension(:), allocatable :: input_vec

    if (.not. this%is_initialised) then
       call stop_program("nn_fingerprint: network not initialised")
       return
    end if

    allocate(input_vec(this%input_dim))
    call this%basis_to_input(basis, input_vec)

    ! Reshape for ATHENA: (features, num_samples)
    allocate(input_2d(this%input_dim, 1))
    input_2d(:, 1) = input_vec

    ! Set inference mode and predict
    call this%network%set_inference_mode()
    output_2d = this%network%predict(input = input_2d)

    fingerprint(1:this%fingerprint_dim) = output_2d(:, 1)

    deallocate(input_vec, input_2d, output_2d)

  end subroutine predict
!###############################################################################


!###############################################################################
  subroutine inverse_design(this, target_fingerprint, basis, &
       fixed_atoms, num_steps, step_size, verbose)
    !! Inverse design: generate or optimise an atomic structure to match
    !! a target RAFFLE descriptor fingerprint.
    !!
    !! Starting from an initial structure (basis), optimises the positions
    !! of non-fixed atoms to minimise the MSE between the computed descriptor
    !! and the target descriptor.
    !!
    !! The optimisation uses gradient-free finite-difference approach on
    !! atomic coordinates, respecting the atom mask.
    implicit none

    ! Arguments
    class(nn_fingerprint_type), intent(inout) :: this
    !! Parent. Instance of nn_fingerprint_type.
    real(real32), dimension(:), intent(in) :: target_fingerprint
    !! Target descriptor fingerprint to match.
    type(basis_type), intent(inout) :: basis
    !! Atomic structure to optimise. Modified in place.
    logical, dimension(:), intent(in) :: fixed_atoms
    !! Boolean mask: .true. = atom is fixed, .false. = atom is optimisable.
    !! Length must equal basis%natom.
    integer, intent(in), optional :: num_steps
    !! Number of optimisation steps. Default: 200.
    real(real32), intent(in), optional :: step_size
    !! Step size for coordinate perturbation. Default: 0.01.
    integer, intent(in), optional :: verbose
    !! Verbosity level. Default: 0.

    ! Local variables
    integer :: nsteps, verb
    real(real32) :: lr
    integer :: step, is, ia, coord, atom_idx, global_idx
    real(real32) :: loss_current, loss_perturbed, grad
    real(real32) :: delta
    real(real32), dimension(:), allocatable :: current_fp
    real(real32), dimension(:), allocatable :: perturbed_fp
    type(basis_type) :: basis_perturbed
    real(real32) :: best_loss
    type(basis_type) :: best_basis

    nsteps = 200
    lr = 0.01_real32
    verb = 0
    if (present(num_steps)) nsteps = num_steps
    if (present(step_size)) lr = step_size
    if (present(verbose)) verb = verbose

    delta = 1.E-4_real32

    allocate(current_fp(this%fingerprint_dim))
    allocate(perturbed_fp(this%fingerprint_dim))

    ! Compute initial loss
    call this%compute_fingerprint(basis, current_fp)
    loss_current = sum((current_fp - target_fingerprint)**2) / &
         real(this%fingerprint_dim, real32)
    best_loss = loss_current
    best_basis = basis

    if (verb > 0) write(*,'(A,I6,A,E12.5)') &
         ' Inverse design step ', 0, ' loss = ', loss_current

    !---------------------------------------------------------------------------
    ! Gradient descent loop with finite differences
    !---------------------------------------------------------------------------
    do step = 1, nsteps
       global_idx = 0

       do is = 1, basis%nspec
          do ia = 1, basis%spec(is)%num
             global_idx = global_idx + 1

             ! Skip fixed atoms
             if (global_idx <= size(fixed_atoms)) then
                if (fixed_atoms(global_idx)) cycle
             end if

             ! Optimise each coordinate (x, y, z)
             do coord = 1, 3
                ! Create perturbed structure
                basis_perturbed = basis
                basis_perturbed%spec(is)%atom(ia, coord) = &
                     basis_perturbed%spec(is)%atom(ia, coord) + delta

                ! Compute perturbed fingerprint and loss
                call this%compute_fingerprint(basis_perturbed, perturbed_fp)
                loss_perturbed = sum( &
                     (perturbed_fp - target_fingerprint)**2 &
                ) / real(this%fingerprint_dim, real32)

                ! Finite-difference gradient
                grad = (loss_perturbed - loss_current) / delta

                ! Update coordinate
                basis%spec(is)%atom(ia, coord) = &
                     basis%spec(is)%atom(ia, coord) - lr * grad
             end do
          end do
       end do

       ! Recompute loss after full update
       call this%compute_fingerprint(basis, current_fp)
       loss_current = sum((current_fp - target_fingerprint)**2) / &
            real(this%fingerprint_dim, real32)

       ! Track best
       if (loss_current < best_loss) then
          best_loss = loss_current
          best_basis = basis
       end if

       if (verb > 0 .and. mod(step, 10) == 0) then
          write(*,'(A,I6,A,E12.5)') &
               ' Inverse design step ', step, ' loss = ', loss_current
       end if

       ! Convergence check
       if (loss_current < 1.E-8_real32) then
          if (verb > 0) write(*,'(A,I6)') &
               ' Inverse design converged at step ', step
          exit
       end if
    end do

    ! Restore best structure
    basis = best_basis

    deallocate(current_fp, perturbed_fp)

    if (verb > 0) write(*,'(A,E12.5)') &
         ' Inverse design final loss = ', best_loss

  end subroutine inverse_design
!###############################################################################


end module raffle__nn_fingerprint
