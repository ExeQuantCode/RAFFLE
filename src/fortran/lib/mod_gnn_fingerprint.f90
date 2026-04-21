module raffle__gnn_fingerprint
  !! Module for graph neural network-based descriptor fingerprint prediction.
  !!
  !! This module implements a GNN fingerprint framework using the ATHENA
  !! Duvenaud message-passing neural network. The GNN naturally respects the
  !! graph topology of atomic structures (atoms = vertices, bonds = edges),
  !! providing a permutation-invariant, variable-size representation.
  !!
  !! Capabilities:
  !!   1. Converts atomic structures to molecular graphs
  !!   2. Learns RAFFLE descriptor fingerprints via message passing
  !!   3. Supports forward inference for predicting descriptors
  !!   4. Enables inverse design (generating structures from target descriptors)
  !!   5. Allows partial atomic optimisation via boolean atom masks
  use raffle__constants, only: real32, pi
  use raffle__io_utils, only: stop_program, print_warning
  use raffle__geom_rw, only: basis_type, geom_write
  use raffle__distribs, only: distribs_base_type, distribs_type
  use raffle__distribs_container, only: distribs_container_type
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use raffle__misc, only: strip_null
  use athena, only: &
       network_type, &
       full_layer_type, &
       adam_optimiser_type, sgd_optimiser_type, &
       clip_type, &
       duvenaud_msgpass_layer_type, &
       graph_type, &
       edge_type, &
       exp_lr_decay_type, &
       random_setup
  use diffstruc, only: array_type
  use raffle__msgpass_layer, only: raffle_msgpass_layer_type
  use raffle__schnet_msgpass_layer, only: schnet_msgpass_layer_type
  use raffle__dimenet_msgpass_layer, only: dimenet_msgpass_layer_type
  use raffle__hybrid_msgpass_layer, only: hybrid_msgpass_layer_type
  implicit none


  private

  public :: gnn_fingerprint_type


  !-----------------------------------------------------------------------------
  ! Default parameters
  !-----------------------------------------------------------------------------
  real(real32), parameter :: DEFAULT_BOND_CUTOFF = 6.0_real32
  !! Default bond cutoff distance in Angstroms for graph edge construction.
  integer, parameter :: DEFAULT_NUM_TIME_STEPS = 3
  !! Default number of message-passing iterations.
  integer, parameter :: DEFAULT_GNN_OUTPUT_DIM = 32
  !! Default dimension of graph-level output from message passing.
  integer, parameter :: DEFAULT_MAX_DEGREE = 12
  !! Default maximum vertex degree (coordination number).


  !-----------------------------------------------------------------------------
  ! GNN fingerprint type
  !-----------------------------------------------------------------------------
  type :: gnn_fingerprint_type
     !! Graph neural network for learning and predicting RAFFLE descriptor
     !! fingerprints from atomic structures represented as molecular graphs.
     logical :: is_initialised = .false.
     !! Whether the network has been initialised.
     logical :: is_trained = .false.
     !! Whether the network has been trained.
     logical :: use_mlip_layer = .false.
     !! Whether to use MLIP-style message passing instead of Duvenaud.
     integer :: layer_type = 0
     !! Layer type: 0=duvenaud, 1=raffle_mlip, 2=schnet, 3=dimenet, 4=hybrid

     integer :: num_species = 0
     !! Number of distinct species in the training set.
     character(len=3), dimension(:), allocatable :: species_list
     !! List of species symbols.

     integer :: num_vertex_features = 0
     !! Number of vertex features: 3 (coords) + num_species (one-hot).
     integer :: num_edge_features = 1
     !! Number of edge features (default: 1 = bond distance).
     integer :: gnn_output_dim = DEFAULT_GNN_OUTPUT_DIM
     !! Dimension of graph-level vector from message passing readout.
     integer :: fingerprint_dim = 0
     !! Dimension of the output fingerprint vector.

     integer :: num_time_steps = DEFAULT_NUM_TIME_STEPS
     !! Number of message-passing iterations.
     integer :: max_degree = DEFAULT_MAX_DEGREE
     !! Maximum vertex degree for Duvenaud bucketing.
     real(real32) :: bond_cutoff = DEFAULT_BOND_CUTOFF
     !! Bond cutoff distance for graph construction.

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

     real(real32), dimension(:), allocatable :: target_mean
     !! Per-dimension fingerprint mean used to normalise training targets.
     real(real32), dimension(:), allocatable :: target_scale
     !! Per-dimension fingerprint scale used to normalise training targets.
     real(real32), dimension(3) :: component_weight = &
          [4.0_real32, 1.0_real32, 1.0_real32]
     !! Relative weight for 2-body, 3-body, 4-body components in loss.
     !! Default heavily biases towards 2-body.
     real(real32), dimension(:), allocatable :: target_comp_weights
     !! Per-dimension component weights applied during training.

     type(network_type) :: network
     !! ATHENA neural network (Duvenaud MPNN + dense head).

   contains
     procedure, pass(this) :: initialise
     !! Initialise the GNN architecture.
     procedure, pass(this) :: basis_to_graph
     !! Convert a basis_type to a graph_type.
     procedure, pass(this) :: compute_fingerprint
     !! Compute the RAFFLE descriptor fingerprint for a structure.
     procedure, pass(this) :: train
     !! Train on structures and their descriptor fingerprints.
     procedure, pass(this) :: predict
     !! Forward inference: predict fingerprint from structure.
     procedure, pass(this) :: inverse_design
     !! Inverse design: optimise structure to match target descriptor.
     procedure, pass(this) :: fingerprint_to_distribs
     !! Convert a flat fingerprint to distribs_base_type.
     procedure, pass(this) :: set_distribution_params
     !! Copy distribution parameters from a container.
  end type gnn_fingerprint_type


contains


!###############################################################################
  subroutine get_element_props(name, atomic_number, covalent_radius)
    !! Look up atomic number and covalent radius for an element.
    !! Returns scaled values suitable for use as vertex features.
    implicit none
    character(len=3), intent(in) :: name
    real(real32), intent(out) :: atomic_number
    real(real32), intent(out) :: covalent_radius

    character(len=3) :: trimmed
    trimmed = adjustl(name)

    select case (trim(trimmed))
    case ('H');   atomic_number = 1;  covalent_radius = 0.31
    case ('He');  atomic_number = 2;  covalent_radius = 0.28
    case ('Li');  atomic_number = 3;  covalent_radius = 1.28
    case ('Be');  atomic_number = 4;  covalent_radius = 0.96
    case ('B');   atomic_number = 5;  covalent_radius = 0.84
    case ('C');   atomic_number = 6;  covalent_radius = 0.76
    case ('N');   atomic_number = 7;  covalent_radius = 0.71
    case ('O');   atomic_number = 8;  covalent_radius = 0.66
    case ('F');   atomic_number = 9;  covalent_radius = 0.57
    case ('Ne');  atomic_number = 10; covalent_radius = 0.58
    case ('Na');  atomic_number = 11; covalent_radius = 1.66
    case ('Mg');  atomic_number = 12; covalent_radius = 1.41
    case ('Al');  atomic_number = 13; covalent_radius = 1.21
    case ('Si');  atomic_number = 14; covalent_radius = 1.11
    case ('P');   atomic_number = 15; covalent_radius = 1.07
    case ('S');   atomic_number = 16; covalent_radius = 1.05
    case ('Cl');  atomic_number = 17; covalent_radius = 1.02
    case ('Ar');  atomic_number = 18; covalent_radius = 1.06
    case ('K');   atomic_number = 19; covalent_radius = 2.03
    case ('Ca');  atomic_number = 20; covalent_radius = 1.76
    case ('Sc');  atomic_number = 21; covalent_radius = 1.70
    case ('Ti');  atomic_number = 22; covalent_radius = 1.60
    case ('V');   atomic_number = 23; covalent_radius = 1.53
    case ('Cr');  atomic_number = 24; covalent_radius = 1.39
    case ('Mn');  atomic_number = 25; covalent_radius = 1.39
    case ('Fe');  atomic_number = 26; covalent_radius = 1.32
    case ('Co');  atomic_number = 27; covalent_radius = 1.26
    case ('Ni');  atomic_number = 28; covalent_radius = 1.24
    case ('Cu');  atomic_number = 29; covalent_radius = 1.32
    case ('Zn');  atomic_number = 30; covalent_radius = 1.22
    case ('Ga');  atomic_number = 31; covalent_radius = 1.22
    case ('Ge');  atomic_number = 32; covalent_radius = 1.20
    case ('As');  atomic_number = 33; covalent_radius = 1.19
    case ('Se');  atomic_number = 34; covalent_radius = 1.20
    case ('Br');  atomic_number = 35; covalent_radius = 1.20
    case ('Mo');  atomic_number = 42; covalent_radius = 1.54
    case ('Ba');  atomic_number = 56; covalent_radius = 2.15
    case ('W');   atomic_number = 74; covalent_radius = 1.62
    case default
       atomic_number = 0
       covalent_radius = 1.0
    end select

    ! Scale: Z/100 puts it in [0, ~1] range; radius already in Angstroms
    atomic_number = atomic_number / 100._real32
  end subroutine get_element_props
!###############################################################################


!###############################################################################
  subroutine set_distribution_params(this, container)
    !! Set distribution parameters from an existing distribs_container_type.
    implicit none

    ! Arguments
    class(gnn_fingerprint_type), intent(inout) :: this
    !! Parent. Instance of gnn_fingerprint_type.
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
  subroutine initialise(this, species_list, &
       num_time_steps, gnn_output_dim, max_degree, &
       hidden_layer_sizes, learning_rate, lr_decay_rate, bond_cutoff, &
       use_mlip_layer, n_rbf, kernel_hidden, layer_type_in, seed)
    !! Initialise the GNN for fingerprint prediction.
    !!
    !! Architecture:
    !!   Duvenaud MPNN or MLIP MPNN (message-passing on atom graph)
    !!     → graph-level vector (gnn_output_dim)
    !!     → Dense hidden layers
    !!     → Dense output (fingerprint_dim, linear)
    implicit none

    ! Arguments
    class(gnn_fingerprint_type), intent(inout) :: this
    !! Parent. Instance of gnn_fingerprint_type.
    character(len=3), dimension(:), intent(in) :: species_list
    !! List of element species to handle.
    integer, intent(in), optional :: num_time_steps
    !! Number of message-passing iterations. Default: 3.
    integer, intent(in), optional :: gnn_output_dim
    !! Graph-level output dimension from readout. Default: 32.
    integer, intent(in), optional :: max_degree
    !! Maximum vertex degree. Default: 12.
    integer, dimension(:), intent(in), optional :: hidden_layer_sizes
    !! Dense head hidden layer sizes. Default: [64].
    real(real32), intent(in), optional :: learning_rate
    !! Learning rate for Adam optimizer. Default: 0.001.
    real(real32), intent(in), optional :: lr_decay_rate
    !! Learning rate decay rate. Default: 1.E-2.
    real(real32), intent(in), optional :: bond_cutoff
    !! Cutoff distance for bonds. Default: 6.0 Angstrom.
    logical, intent(in), optional :: use_mlip_layer
    !! Use MLIP-style message passing. Default: .false.
    integer, intent(in), optional :: n_rbf
    !! Number of RBF basis functions (MLIP only). Default: 20.
    integer, intent(in), optional :: kernel_hidden
    !! Kernel MLP hidden width (MLIP only). Default: 64.
    integer, intent(in), optional :: layer_type_in
    !! Layer type: 0=duvenaud, 1=raffle_mlip, 2=schnet, 3=dimenet, 4=hybrid.
    integer, intent(in), optional :: seed

    ! Local variables
    integer :: i, num_hidden
    integer :: seed_
    integer, dimension(:), allocatable :: h_sizes
    real(real32) :: lr, lr_decay_rate_
    class(clip_type), allocatable :: clip
    type(exp_lr_decay_type) :: lr_decay


    seed_ = 42
    if(present(seed)) seed_ = seed
    call random_setup(seed_, restart=.false.)

    ! Set species
    this%num_species = size(species_list)
    if (allocated(this%species_list)) deallocate(this%species_list)
    allocate(this%species_list(this%num_species))
    this%species_list = species_list

    ! Vertex features: 3 coordinates + one-hot species + atomic_number + covalent_radius
    this%num_vertex_features = 3 + this%num_species + 2
    this%num_edge_features = 1  ! bond distance

    ! Optional parameters
    this%num_time_steps = DEFAULT_NUM_TIME_STEPS
    if (present(num_time_steps)) this%num_time_steps = num_time_steps

    this%gnn_output_dim = DEFAULT_GNN_OUTPUT_DIM
    if (present(gnn_output_dim)) this%gnn_output_dim = gnn_output_dim

    this%max_degree = DEFAULT_MAX_DEGREE
    if (present(max_degree)) this%max_degree = max_degree

    this%bond_cutoff = DEFAULT_BOND_CUTOFF
    if (present(bond_cutoff)) this%bond_cutoff = bond_cutoff

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

    ! Hidden layer sizes for the dense head
    if (present(hidden_layer_sizes)) then
       num_hidden = size(hidden_layer_sizes)
       allocate(h_sizes(num_hidden))
       h_sizes = hidden_layer_sizes
    else
       num_hidden = 1
       allocate(h_sizes(1))
       h_sizes = [64]
    end if

    ! Learning rate
    lr = 0.001_real32
    lr_decay_rate_ = 1.E-2_real32
    if(present(learning_rate)) lr = learning_rate
    if(present(lr_decay_rate)) lr_decay_rate_ = lr_decay_rate

    ! MLIP layer flag
    this%use_mlip_layer = .false.
    if (present(use_mlip_layer)) this%use_mlip_layer = use_mlip_layer

    ! Layer type selection
    this%layer_type = 0  ! default: Duvenaud
    if (this%use_mlip_layer) this%layer_type = 1
    if (present(layer_type_in)) this%layer_type = layer_type_in

    !--------------------------------------------------------------------------
    ! Build network: Message-passing layer → Dense head
    !--------------------------------------------------------------------------
    block
      integer :: n_rbf_, kernel_hidden_
      n_rbf_ = 20
      kernel_hidden_ = 64
      if (present(n_rbf)) n_rbf_ = n_rbf
      if (present(kernel_hidden)) kernel_hidden_ = kernel_hidden

      select case (this%layer_type)
      case (1)
         ! MLIP-style (raffle) message-passing layer
         call this%network%add(raffle_msgpass_layer_type( &
              num_time_steps = this%num_time_steps, &
              num_vertex_features = [this%num_vertex_features], &
              num_edge_features = [this%num_edge_features], &
              num_outputs = this%gnn_output_dim, &
              n_rbf = n_rbf_, &
              rbf_cutoff = this%bond_cutoff, &
              kernel_hidden = kernel_hidden_, &
              message_activation = 'swish', &
              readout_activation = 'none', &
              kernel_initialiser = 'glorot_normal' &
         ))
      case (2)
         ! SchNet-style continuous filter
         call this%network%add(schnet_msgpass_layer_type( &
              num_time_steps = this%num_time_steps, &
              num_vertex_features = [this%num_vertex_features], &
              num_edge_features = [this%num_edge_features], &
              num_outputs = this%gnn_output_dim, &
              n_rbf = n_rbf_, &
              rbf_cutoff = this%bond_cutoff, &
              kernel_hidden = kernel_hidden_, &
              message_activation = 'swish', &
              readout_activation = 'none', &
              kernel_initialiser = 'glorot_normal' &
         ))
      case (3)
         ! DimeNet-inspired angular layer
         call this%network%add(dimenet_msgpass_layer_type( &
              num_time_steps = this%num_time_steps, &
              num_vertex_features = [this%num_vertex_features], &
              num_edge_features = [this%num_edge_features], &
              num_outputs = this%gnn_output_dim, &
              n_rbf = n_rbf_, &
              rbf_cutoff = this%bond_cutoff, &
              kernel_hidden = kernel_hidden_, &
              message_activation = 'swish', &
              readout_activation = 'none', &
              kernel_initialiser = 'glorot_normal' &
         ))
      case (4)
         ! Hybrid MLIP-inspired with gating
         call this%network%add(hybrid_msgpass_layer_type( &
              num_time_steps = this%num_time_steps, &
              num_vertex_features = [this%num_vertex_features], &
              num_edge_features = [this%num_edge_features], &
              num_outputs = this%gnn_output_dim, &
              n_rbf = n_rbf_, &
              rbf_cutoff = this%bond_cutoff, &
              kernel_hidden = kernel_hidden_, &
              message_activation = 'swish', &
              readout_activation = 'none', &
              kernel_initialiser = 'glorot_normal' &
         ))
      case default
         ! Duvenaud message-passing layer
         call this%network%add(duvenaud_msgpass_layer_type( &
              num_time_steps = this%num_time_steps, &
              num_vertex_features = [this%num_vertex_features], &
              num_edge_features = [this%num_edge_features], &
              num_outputs = this%gnn_output_dim, &
              kernel_initialiser = 'glorot_normal', &
              readout_activation = 'none', &
              min_vertex_degree = 0, &
              max_vertex_degree = this%max_degree &
         ))
      end select
    end block

    ! 2. Dense hidden layers
    do i = 1, num_hidden
       if (i == 1) then
          call this%network%add(full_layer_type( &
               num_inputs = this%gnn_output_dim, &
               num_outputs = h_sizes(i), &
               activation = 'leaky_relu', &
               kernel_initialiser = 'he_normal' &
          ))
       else
          call this%network%add(full_layer_type( &
               num_outputs = h_sizes(i), &
               activation = 'leaky_relu', &
               kernel_initialiser = 'he_normal' &
          ))
       end if
    end do

    ! 3. Output layer (linear activation for regression)
    call this%network%add(full_layer_type( &
         num_outputs = this%fingerprint_dim, &
         activation = 'none', &
         kernel_initialiser = 'glorot_normal' &
    ))

    ! 4. Compile network
    allocate(clip, source=clip_type( &
         clip_min = -1.E-1_real32, &
         clip_max = 1.E-1_real32, &
         clip_norm = 1.E-1_real32 &
    ))
    lr_decay = exp_lr_decay_type(lr_decay_rate_)
    lr_decay%iterate_per_epoch = .true.
    call this%network%compile( &
         optimiser = adam_optimiser_type( &
              learning_rate = lr, &
              clip_dict = clip, &
              lr_decay = lr_decay &
         ), &
         loss_method = 'mse', &
         metrics = ['loss'], &
         verbose = 0 &
    )

    this%is_initialised = .true.

    ! Initialise default normalisation so predict() works before training
    if (.not. allocated(this%target_mean)) then
       allocate(this%target_mean(this%fingerprint_dim))
       this%target_mean = 0._real32
    end if
    if (.not. allocated(this%target_scale)) then
       allocate(this%target_scale(this%fingerprint_dim))
       this%target_scale = 1._real32
    end if
    if (.not. allocated(this%target_comp_weights)) then
       allocate(this%target_comp_weights(this%fingerprint_dim))
       this%target_comp_weights = 1._real32
    end if

    if (allocated(h_sizes)) deallocate(h_sizes)

  end subroutine initialise
!###############################################################################


!###############################################################################
  subroutine basis_to_graph(this, basis, graph)
    !! Convert a basis_type atomic structure to a graph_type.
    !!
    !! Atoms become vertices with features [x, y, z, one_hot_species, Z/100, cov_radius].
    !! Bonds (pairs within cutoff) become edges with feature [distance].
    !! Self-loops are added and the graph is converted to sparse (CSR).
    implicit none

    ! Arguments
    class(gnn_fingerprint_type), intent(in) :: this
    !! Parent. Instance of gnn_fingerprint_type.
    type(basis_type), intent(in) :: basis
    !! Atomic structure to convert.
    type(graph_type), intent(out) :: graph
    !! Output molecular graph.

    ! Local variables
    integer :: is, ia, js, ja, atom_i, atom_j, species_idx
    integer :: num_atoms
    real(real32) :: coord_scale
    real(real32) :: dx, dy, dz, dist
    real(real32), dimension(3) :: pos_i, pos_j, shift
    real(real32), dimension(3) :: cart_i, cart_j
    real(real32) :: elem_z, elem_cov_r

    num_atoms = basis%natom
    coord_scale = max(this%bond_cutoff, 1.E-6_real32)

    ! Set up vertices (dense mode: allocates vertex(:) array)
    call graph%set_num_vertices(num_atoms, &
         num_vertex_features = this%num_vertex_features)

    ! Allocate feature array for each vertex
    do atom_i = 1, num_atoms
       allocate(graph%vertex(atom_i)%feature(this%num_vertex_features))
       graph%vertex(atom_i)%feature = 0._real32
    end do

    ! Fill vertex features: [x, y, z, one_hot_species, Z/100, cov_radius]
    ! In dense mode, vertex features are stored in vertex(i)%feature(:)
    atom_i = 0
    do is = 1, basis%nspec
       ! Find species index
       species_idx = 0
       do ia = 1, this%num_species
          if (this%species_list(ia) == basis%spec(is)%name) then
             species_idx = ia
             exit
          end if
       end do

       ! Look up element properties for this species
       call get_element_props(basis%spec(is)%name, elem_z, elem_cov_r)

       do ia = 1, basis%spec(is)%num
          atom_i = atom_i + 1

          ! Convert fractional to Cartesian if needed
          if (basis%lcart) then
             cart_i(1) = basis%spec(is)%atom(ia, 1)
             cart_i(2) = basis%spec(is)%atom(ia, 2)
             cart_i(3) = basis%spec(is)%atom(ia, 3)
          else
             cart_i(1) = basis%spec(is)%atom(ia,1) * basis%lat(1,1) &
                  + basis%spec(is)%atom(ia,2) * basis%lat(2,1) &
                  + basis%spec(is)%atom(ia,3) * basis%lat(3,1)
             cart_i(2) = basis%spec(is)%atom(ia,1) * basis%lat(1,2) &
                  + basis%spec(is)%atom(ia,2) * basis%lat(2,2) &
                  + basis%spec(is)%atom(ia,3) * basis%lat(3,2)
             cart_i(3) = basis%spec(is)%atom(ia,1) * basis%lat(1,3) &
                  + basis%spec(is)%atom(ia,2) * basis%lat(2,3) &
                  + basis%spec(is)%atom(ia,3) * basis%lat(3,3)
          end if

          ! Dense mode: write to vertex(atom_i)%feature(:)
          graph%vertex(atom_i)%feature(1) = cart_i(1) / coord_scale
          graph%vertex(atom_i)%feature(2) = cart_i(2) / coord_scale
          graph%vertex(atom_i)%feature(3) = cart_i(3) / coord_scale

          ! One-hot species encoding
          graph%vertex(atom_i)%feature(4:3+this%num_species) = 0._real32
          if (species_idx > 0) then
             graph%vertex(atom_i)%feature(3 + species_idx) = 1._real32
          end if

          ! Atomic number (scaled) and covalent radius
          graph%vertex(atom_i)%feature(3 + this%num_species + 1) = elem_z
          graph%vertex(atom_i)%feature(3 + this%num_species + 2) = elem_cov_r
       end do
    end do

    ! Build edges: all pairs within bond cutoff
    ! First pass: count edges
    atom_i = 0
    do is = 1, basis%nspec
       do ia = 1, basis%spec(is)%num
          atom_i = atom_i + 1
          if (basis%lcart) then
             pos_i = basis%spec(is)%atom(ia, 1:3)
          else
             pos_i(1) = basis%spec(is)%atom(ia,1) * basis%lat(1,1) &
                  + basis%spec(is)%atom(ia,2) * basis%lat(2,1) &
                  + basis%spec(is)%atom(ia,3) * basis%lat(3,1)
             pos_i(2) = basis%spec(is)%atom(ia,1) * basis%lat(1,2) &
                  + basis%spec(is)%atom(ia,2) * basis%lat(2,2) &
                  + basis%spec(is)%atom(ia,3) * basis%lat(3,2)
             pos_i(3) = basis%spec(is)%atom(ia,1) * basis%lat(1,3) &
                  + basis%spec(is)%atom(ia,2) * basis%lat(2,3) &
                  + basis%spec(is)%atom(ia,3) * basis%lat(3,3)
          end if

          atom_j = 0
          do js = 1, basis%nspec
             do ja = 1, basis%spec(js)%num
                atom_j = atom_j + 1
                if (atom_j <= atom_i) cycle  ! avoid double-counting

                if (basis%lcart) then
                   pos_j = basis%spec(js)%atom(ja, 1:3)
                else
                   pos_j(1) = basis%spec(js)%atom(ja,1) * basis%lat(1,1) &
                        + basis%spec(js)%atom(ja,2) * basis%lat(2,1) &
                        + basis%spec(js)%atom(ja,3) * basis%lat(3,1)
                   pos_j(2) = basis%spec(js)%atom(ja,1) * basis%lat(1,2) &
                        + basis%spec(js)%atom(ja,2) * basis%lat(2,2) &
                        + basis%spec(js)%atom(ja,3) * basis%lat(3,2)
                   pos_j(3) = basis%spec(js)%atom(ja,1) * basis%lat(1,3) &
                        + basis%spec(js)%atom(ja,2) * basis%lat(2,3) &
                        + basis%spec(js)%atom(ja,3) * basis%lat(3,3)
                end if

                ! Minimum image convention for periodic systems
                shift = pos_j - pos_i
                if (basis%pbc(1)) &
                     shift(1) = shift(1) - basis%lat(1,1) * &
                     nint(shift(1) / basis%lat(1,1))
                if (basis%pbc(2)) &
                     shift(2) = shift(2) - basis%lat(2,2) * &
                     nint(shift(2) / basis%lat(2,2))
                if (basis%pbc(3)) &
                     shift(3) = shift(3) - basis%lat(3,3) * &
                     nint(shift(3) / basis%lat(3,3))

                dist = sqrt(shift(1)**2 + shift(2)**2 + shift(3)**2)

                if (dist > 0._real32 .and. dist <= this%bond_cutoff) then
                   ! Add undirected edge (both directions)
                   call graph%add_edge( &
                        index = [atom_i, atom_j], &
                        feature = [dist] &
                   )
                end if
             end do
          end do
       end do
    end do

    ! Add self-loops and convert to sparse CSR
    call graph%add_self_loops(features = [0._real32])
    if (.not. graph%is_sparse) call graph%convert_to_sparse()

  end subroutine basis_to_graph
!###############################################################################


!###############################################################################
  subroutine compute_fingerprint(this, basis, fingerprint)
    !! Compute the RAFFLE descriptor fingerprint for a given structure.
    !!
    !! Uses distribs_type%calculate() to compute 2/3/4-body distributions
    !! and flattens them into a single vector.
    implicit none

    ! Arguments
    class(gnn_fingerprint_type), intent(in) :: this
    !! Parent. Instance of gnn_fingerprint_type.
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
    class(gnn_fingerprint_type), intent(in) :: this
    !! Parent. Instance of gnn_fingerprint_type.
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
    !! Train the GNN on a set of structures and their descriptor fingerprints.
    !!
    !! For each structure:
    !!   1. Convert to graph via basis_to_graph
    !!   2. Compute target fingerprint via compute_fingerprint
    !!   3. Train network: graph → fingerprint
    implicit none

    ! Arguments
    class(gnn_fingerprint_type), intent(inout) :: this
    !! Parent. Instance of gnn_fingerprint_type.
    type(basis_type), dimension(:), intent(in) :: structures
    !! Training structures.
    integer, intent(in), optional :: num_epochs
    !! Number of training epochs. Default: 100.
    integer, intent(in), optional :: batch_size
    !! Batch size for optimisation. Default: min(16, num_strucs).
    integer, intent(in), optional :: verbose
    !! Verbosity level. Default: 0.

    ! Local variables
    integer :: n, num_strucs, epochs, batch_size_, verb
    type(graph_type), dimension(:,:), allocatable :: graphs_in
    type(array_type), dimension(1,1) :: output_array
    real(real32), dimension(:,:), allocatable :: target_data
    real(real32), dimension(:), allocatable :: fp_vec
    real(real32), dimension(:,:), allocatable :: centred_target_data
    real(real32), dimension(:), allocatable :: comp_weights
    integer :: offset, dim2, dim3, dim4, i

    if (.not. this%is_initialised) then
       call stop_program("gnn_fingerprint: network not initialised")
       return
    end if

    num_strucs = size(structures)
    epochs = 100
    batch_size_ = min(num_strucs, 32)
    verb = 0
    if (present(num_epochs)) epochs = num_epochs
    if (present(batch_size)) batch_size_ = max(1, min(num_strucs, batch_size))
    if (present(verbose)) verb = verbose

    ! Build graph inputs: dimension(1, num_strucs)
    allocate(graphs_in(1, num_strucs))
    allocate(target_data(this%fingerprint_dim, num_strucs))
    allocate(fp_vec(this%fingerprint_dim))

    do n = 1, num_strucs
       call this%basis_to_graph(structures(n), graphs_in(1, n))
       call this%compute_fingerprint(structures(n), fp_vec)
       if (.not. all(ieee_is_finite(fp_vec))) then
          call stop_program("gnn_fingerprint: non-finite target fingerprint")
          return
       end if
       target_data(:, n) = fp_vec
    end do

    if (allocated(this%target_mean)) deallocate(this%target_mean)
    if (allocated(this%target_scale)) deallocate(this%target_scale)
    allocate(this%target_mean(this%fingerprint_dim))
    allocate(this%target_scale(this%fingerprint_dim))
    this%target_mean = sum(target_data, dim = 2) / real(num_strucs, real32)

    allocate(centred_target_data(this%fingerprint_dim, num_strucs))
    centred_target_data = &
         target_data - spread(this%target_mean, dim = 2, ncopies = num_strucs)
    this%target_scale = sqrt( &
         sum(centred_target_data**2, dim = 2) / real(max(num_strucs, 1), real32) &
    )
    where (this%target_scale < 1.E-6_real32)
       this%target_scale = 1._real32
    end where
    target_data = &
         centred_target_data / spread(this%target_scale, dim = 2, ncopies = num_strucs)
    deallocate(centred_target_data)

    ! Apply per-component weighting (2-body, 3-body, 4-body)
    dim2 = this%nbins(1) * this%num_pairs
    dim3 = this%nbins(2) * this%num_species
    dim4 = this%nbins(3) * this%num_species
    allocate(comp_weights(this%fingerprint_dim))
    comp_weights = 1._real32
    offset = 0
    ! 2-body block
    do i = 1, dim2
       offset = offset + 1
       comp_weights(offset) = sqrt(this%component_weight(1))
    end do
    ! 3-body block
    do i = 1, dim3
       offset = offset + 1
       comp_weights(offset) = sqrt(this%component_weight(2))
    end do
    ! 4-body block
    do i = 1, dim4
       offset = offset + 1
       comp_weights(offset) = sqrt(this%component_weight(3))
    end do
    ! Apply weights to targets (network learns weighted targets)
    target_data = target_data * spread(comp_weights, dim=2, ncopies=num_strucs)
    ! Store weights for de-weighting at prediction time
    if (allocated(this%target_comp_weights)) deallocate(this%target_comp_weights)
    allocate(this%target_comp_weights, source=comp_weights)
    deallocate(comp_weights)

    if (.not. all(ieee_is_finite(target_data))) then
       call stop_program("gnn_fingerprint: non-finite normalised training targets")
       return
    end if

    ! Prepare output as array_type
    call output_array(1,1)%allocate(array_shape = &
         [this%fingerprint_dim, num_strucs])
    output_array(1,1)%val = target_data

    ! Train: single call to ATHENA with all epochs.
    call this%network%train( &
         graphs_in, &
         output_array, &
         num_epochs = epochs, &
         batch_size = batch_size_, &
         shuffle_batches = .true., &
         verbose = verb &
    )

    this%is_trained = .true.

    deallocate(graphs_in, target_data, fp_vec)

  end subroutine train
!###############################################################################


!###############################################################################
  subroutine predict(this, basis, fingerprint)
    !! Forward inference: predict a RAFFLE descriptor fingerprint from
    !! an atomic structure using the trained GNN.
    implicit none

    ! Arguments
    class(gnn_fingerprint_type), intent(inout) :: this
    !! Parent. Instance of gnn_fingerprint_type.
    type(basis_type), intent(in) :: basis
    !! Atomic structure input.
    real(real32), dimension(:), intent(out) :: fingerprint
    !! Predicted fingerprint vector.

    ! Local variables
    type(graph_type), dimension(1,1) :: graphs_in
    integer :: leaf_id

    if (.not. this%is_initialised) then
       call stop_program("gnn_fingerprint: network not initialised")
       return
    end if

    ! Convert structure to graph
    call this%basis_to_graph(basis, graphs_in(1,1))

    ! Set up for inference
    call this%network%set_batch_size(1)
    call this%network%set_inference_mode()

    ! Forward pass
    call this%network%forward(graphs_in)

    ! Extract output from leaf layer
    leaf_id = this%network%leaf_vertices(1)
    fingerprint(1:this%fingerprint_dim) = &
         this%network%model(leaf_id)%layer%output(1,1)%val(:, 1)

    if (allocated(this%target_scale) .and. allocated(this%target_mean)) then
       ! Undo component weighting, then undo normalisation
       if (allocated(this%target_comp_weights)) then
          where (this%target_comp_weights > 1.E-12_real32)
             fingerprint(1:this%fingerprint_dim) = &
                  fingerprint(1:this%fingerprint_dim) / &
                  this%target_comp_weights
          end where
       end if
       fingerprint = fingerprint * this%target_scale + this%target_mean
    end if

    if (.not. all(ieee_is_finite(fingerprint))) then
       call stop_program("gnn_fingerprint: non-finite prediction")
       return
    end if

  end subroutine predict
!###############################################################################


!###############################################################################
  subroutine inverse_design(this, target_fingerprint, basis, &
       fixed_atoms, num_steps, step_size, verbose, use_predict)
    !! Inverse design: optimise atomic positions to match a target
    !! fingerprint.
    !!
    !! Uses central-difference gradients with momentum.
    !! When use_predict is true, the GNN predict() is used instead of
    !! the analytical RAFFLE descriptor. This lets the trained GNN
    !! contribute to the inverse design landscape.
    implicit none

    ! Arguments
    class(gnn_fingerprint_type), intent(inout) :: this
    !! Parent. Instance of gnn_fingerprint_type.
    real(real32), dimension(:), intent(in) :: target_fingerprint
    !! Target descriptor fingerprint to match.
    type(basis_type), intent(inout) :: basis
    !! Atomic structure to optimise. Modified in place.
    logical, dimension(:), intent(in) :: fixed_atoms
    !! Boolean mask: .true. = fixed, .false. = optimisable.
    integer, intent(in), optional :: num_steps
    !! Number of optimisation steps. Default: 500.
    real(real32), intent(in), optional :: step_size
    !! Initial step size for gradient descent. Default: 1.0.
    integer, intent(in), optional :: verbose
    !! Verbosity level. Default: 0.
    logical, intent(in), optional :: use_predict
    !! If .true., use GNN predict() instead of compute_fingerprint().

    ! Local variables
    integer :: nsteps, verb
    real(real32) :: lr
    integer :: step, is, ia, atom_i, coord, global_idx, num_movable
    real(real32) :: loss_current, grad_norm
    real(real32) :: delta
    real(real32), dimension(:), allocatable :: current_fp
    real(real32), dimension(:), allocatable :: fp_fwd, fp_bwd, residual
    type(basis_type) :: basis_perturbed
    real(real32) :: best_loss
    type(basis_type) :: best_basis
    ! Momentum storage: (3, num_movable_atoms)
    real(real32), dimension(:,:), allocatable :: velocity
    real(real32), dimension(:,:), allocatable :: grad_store
    integer, dimension(:,:), allocatable :: movable_map
    integer :: mi
    real(real32), parameter :: momentum = 0.7_real32
    real(real32), parameter :: lr_decay = 1.0_real32
    logical :: do_predict

    integer :: l
    integer :: root_id, num_elements, ifeature
    type(array_type) :: y(1,1)
    type(array_type), pointer :: loss
    type(adam_optimiser_type) :: opt
    real(real32), dimension(:), allocatable :: x_flat, x_grad
    type(graph_type) :: graphs_in(1,1)

    nsteps = 500
    lr = 1.0_real32
    verb = 0
    do_predict = .false.
    if (present(num_steps)) nsteps = num_steps
    if (present(step_size)) lr = step_size
    if (present(verbose)) verb = verbose
    if (present(use_predict)) do_predict = use_predict

    num_elements = basis%natom * this%num_vertex_features

    opt = adam_optimiser_type( &
         learning_rate=lr, &
         lr_decay = exp_lr_decay_type(lr_decay), &
         clip_dict=clip_type( &
              clip_min = -1.E-1_real32, &
              clip_max = 1.E-1_real32, &
              clip_norm = 1.E-1_real32 &
         ) )
    call opt%init(num_params=num_elements)

    ! convert target fingerprint to array_type
    call y(1,1)%allocate(array_shape = [size(target_fingerprint), 1])
    y(1,1)%val(:,1) = target_fingerprint


    ! get the input layer id
    root_id = this%network%auto_graph%vertex(this%network%root_vertices(1))%id

    call this%network%set_batch_size(1)
    call this%network%set_inference_mode()

    ! convert basis to cartesian coordinates
    if(.not. basis%lcart) then
       call basis%convert()
    end if
    do is = 1, basis%nspec

       ! strip null from the species names
       basis%spec(is)%name = strip_null(basis%spec(is)%name)
    end do

    !--------------------------------------------------------------------------
    ! Gradient descent loop with momentum
    !--------------------------------------------------------------------------
    do step = 1, nsteps

       ! recalculate the edge features at every step
       call this%basis_to_graph(basis, graphs_in(1,1))

       call this%network%forward(graphs_in, input_requires_grad = .true.)

       this%network%expected_array = y
       loss => this%network%loss_eval(1,1)
       call loss%grad_reverse()

       if(associated(this%network%model(root_id)%layer%output(1,1)%grad))then
          x_grad = reshape( &
               this%network%model(root_id)%layer%output(1,1)%grad%val, &
               [ num_elements ] &
          )
       else
          if(.not.allocated(x_grad)) allocate(x_grad(num_elements))
          x_grad = 0._real32
       end if

       ! clip gradients
       call opt%clip_dict%apply(num_elements, x_grad)

       write(*,'(A,I5,A,ES12.4)') "  step=",step, " loss=",sum(loss%val)

       ! gradient descent update on x
       x_flat = reshape( graphs_in(1,1)%vertex_features, [ num_elements ] )

       call opt%minimise(param=x_flat, gradient=x_grad)

       atom_i = 0
       do is = 1, basis%nspec
          do ia = 1, basis%spec(is)%num
             atom_i = atom_i + 1
             ifeature = (atom_i-1)*this%num_vertex_features
             if(.not. fixed_atoms(atom_i))then
                basis%spec(is)%atom(ia,1:3) = &
                     x_flat( ifeature + 1: ifeature + 3 ) * this%bond_cutoff
             end if
          end do
       end do

       ! clean up computation graph
       call loss%nullify_graph()
       deallocate(loss)
       nullify(loss)

       ! reset network gradients (do NOT call network%update)
       call this%network%reset_gradients()

    end do

    basis%lcart = .true.
    open(unit=100, file="POSCAR_optimized_structure", status="replace")
    call geom_write(100, basis)
    close(100)

  end subroutine inverse_design
!###############################################################################


end module raffle__gnn_fingerprint
