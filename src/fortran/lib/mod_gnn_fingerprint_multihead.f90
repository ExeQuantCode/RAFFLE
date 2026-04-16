module raffle__gnn_fingerprint_multihead
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
  use raffle__misc_linalg, only: cross, get_angle, get_improper_dihedral_angle, &
       inverse_3x3
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
  use diffstruc, only: array_type, ds_sum => sum, operator(+), operator(-), &
       operator(*), operator(/), operator(**)
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
  type :: multigraph_data_type
     !! Cached deterministic graph metadata used for exact position gradients.
     real(real32), dimension(:,:), allocatable :: positions
     integer, dimension(:), allocatable :: species_index
     real(real32), dimension(:), allocatable :: atomic_numbers
     real(real32), dimension(:), allocatable :: cov_radii
     integer, dimension(:,:), allocatable :: pair_atoms
     real(real32), dimension(:,:), allocatable :: pair_shift
     real(real32), dimension(:), allocatable :: pair_distance
     integer, dimension(:,:), allocatable :: angle_triplets
     integer, dimension(:,:), allocatable :: triplet_atoms
     integer, dimension(:,:), allocatable :: improper_quads
  end type multigraph_data_type

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
     integer :: num_pair_vertex_features = 0
     !! Number of 3-body graph vertex features derived from atomic pairs.
     integer :: num_triplet_vertex_features = 0
     !! Number of 4-body graph vertex features derived from atomic triplets.
     integer :: gnn_output_dim = DEFAULT_GNN_OUTPUT_DIM
     !! Dimension of graph-level vector from message passing readout.
     integer :: fingerprint_dim_2body = 0
     !! Dimension of the 2-body fingerprint block.
     integer :: fingerprint_dim_3body = 0
     !! Dimension of the 3-body fingerprint block.
     integer :: fingerprint_dim_4body = 0
     !! Dimension of the 4-body fingerprint block.
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
     !! ATHENA neural network for the 2-body branch.
     type(network_type) :: network_3body
     !! ATHENA neural network for the 3-body branch.
     type(network_type) :: network_4body
     !! ATHENA neural network for the 4-body branch.

   contains
     procedure, pass(this) :: initialise
     !! Initialise the GNN architecture.
     procedure, pass(this) :: basis_to_graph
     !! Convert a basis_type to a graph_type.
     procedure, pass(this) :: basis_to_graph_3body
     !! Convert a basis_type to the fixed 3-body graph.
     procedure, pass(this) :: basis_to_graph_4body
     !! Convert a basis_type to the fixed 4-body graph.
     procedure, pass(this) :: compute_fingerprint
     !! Compute the RAFFLE descriptor fingerprint for a structure.
     procedure, pass(this) :: compute_fingerprint_components
     !! Compute separate 2-body, 3-body, and 4-body fingerprints.
     procedure, pass(this) :: compute_gradients
     !! Compute fingerprint gradients with respect to atomic positions.
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
  subroutine extract_structure_data( &
       this, basis, positions, species_index, atomic_numbers, cov_radii)
    !! Convert a basis into flat Cartesian arrays used by all graph builders.
    implicit none

    class(gnn_fingerprint_type), intent(in) :: this
    type(basis_type), intent(in) :: basis
    real(real32), dimension(:,:), allocatable, intent(out) :: positions
    integer, dimension(:), allocatable, intent(out) :: species_index
    real(real32), dimension(:), allocatable, intent(out) :: atomic_numbers
    real(real32), dimension(:), allocatable, intent(out) :: cov_radii

    integer :: is, ia, atom_i, species_idx
    real(real32) :: elem_z, elem_cov_r
    character(len=3) :: species_name

    allocate(positions(3, basis%natom))
    allocate(species_index(basis%natom))
    allocate(atomic_numbers(basis%natom))
    allocate(cov_radii(basis%natom))

    atom_i = 0
    do is = 1, basis%nspec
       species_idx = 0
       species_name = strip_null(basis%spec(is)%name)
       do ia = 1, this%num_species
          if (trim(strip_null(this%species_list(ia))) == trim(species_name)) then
             species_idx = ia
             exit
          end if
       end do
       call get_element_props(species_name, elem_z, elem_cov_r)

       do ia = 1, basis%spec(is)%num
          atom_i = atom_i + 1
          species_index(atom_i) = species_idx
          atomic_numbers(atom_i) = elem_z
          cov_radii(atom_i) = elem_cov_r
          if (basis%lcart) then
             positions(:, atom_i) = basis%spec(is)%atom(ia, 1:3)
          else
             positions(1, atom_i) = basis%spec(is)%atom(ia,1) * basis%lat(1,1) &
                  + basis%spec(is)%atom(ia,2) * basis%lat(2,1) &
                  + basis%spec(is)%atom(ia,3) * basis%lat(3,1)
             positions(2, atom_i) = basis%spec(is)%atom(ia,1) * basis%lat(1,2) &
                  + basis%spec(is)%atom(ia,2) * basis%lat(2,2) &
                  + basis%spec(is)%atom(ia,3) * basis%lat(3,2)
             positions(3, atom_i) = basis%spec(is)%atom(ia,1) * basis%lat(1,3) &
                  + basis%spec(is)%atom(ia,2) * basis%lat(2,3) &
                  + basis%spec(is)%atom(ia,3) * basis%lat(3,3)
          end if
       end do
    end do

  end subroutine extract_structure_data
!###############################################################################


!###############################################################################
  function minimum_image_delta(basis, delta_cart) result(delta_mic)
    !! Apply a general minimum-image convention using the lattice matrix.
    implicit none

    type(basis_type), intent(in) :: basis
    real(real32), dimension(3), intent(in) :: delta_cart
    real(real32), dimension(3) :: delta_mic

    real(real32), dimension(3,3) :: lattice_matrix, lattice_inverse
    real(real32), dimension(3) :: delta_frac
    integer :: i

    delta_mic = delta_cart
    if (.not. any(basis%pbc)) return

    lattice_matrix = transpose(basis%lat)
    lattice_inverse = inverse_3x3(lattice_matrix)
    delta_frac = matmul(lattice_inverse, delta_cart)
    do i = 1, 3
       if (basis%pbc(i)) delta_frac(i) = delta_frac(i) - nint(delta_frac(i))
    end do
    delta_mic = matmul(lattice_matrix, delta_frac)

  end function minimum_image_delta
!###############################################################################


!###############################################################################
  subroutine build_pair_data(this, basis, positions, pair_atoms, pair_shift, &
       pair_distance)
    !! Build the fixed set of 2-body pairs used by all graph branches.
    implicit none

    class(gnn_fingerprint_type), intent(in) :: this
    type(basis_type), intent(in) :: basis
    real(real32), dimension(:,:), intent(in) :: positions
    integer, dimension(:,:), allocatable, intent(out) :: pair_atoms
    real(real32), dimension(:,:), allocatable, intent(out) :: pair_shift
    real(real32), dimension(:), allocatable, intent(out) :: pair_distance

    integer :: atom_i, atom_j, num_pairs_local
    real(real32), dimension(3) :: delta
    real(real32) :: dist

    num_pairs_local = 0
    do atom_i = 1, size(positions, 2) - 1
       do atom_j = atom_i + 1, size(positions, 2)
          delta = minimum_image_delta(basis, positions(:, atom_j) - &
               positions(:, atom_i))
          dist = sqrt(sum(delta**2))
          if (dist > 0._real32 .and. dist <= this%bond_cutoff) then
             num_pairs_local = num_pairs_local + 1
          end if
       end do
    end do

    allocate(pair_atoms(2, num_pairs_local))
    allocate(pair_shift(3, num_pairs_local))
    allocate(pair_distance(num_pairs_local))

    num_pairs_local = 0
    do atom_i = 1, size(positions, 2) - 1
       do atom_j = atom_i + 1, size(positions, 2)
          delta = minimum_image_delta(basis, positions(:, atom_j) - &
               positions(:, atom_i))
          dist = sqrt(sum(delta**2))
          if (dist > 0._real32 .and. dist <= this%bond_cutoff) then
             num_pairs_local = num_pairs_local + 1
             pair_atoms(:, num_pairs_local) = [atom_i, atom_j]
             pair_shift(:, num_pairs_local) = delta
             pair_distance(num_pairs_local) = dist
          end if
       end do
    end do

  end subroutine build_pair_data
!###############################################################################


!###############################################################################
  subroutine build_pair_vertex_features(this, positions, species_index, pair_atoms, &
       pair_shift, pair_distance, pair_index, features)
    !! Construct deterministic pair-vertex features for the 3-body graph.
    implicit none

    class(gnn_fingerprint_type), intent(in) :: this
    real(real32), dimension(:,:), intent(in) :: positions
    integer, dimension(:), intent(in) :: species_index
    integer, dimension(:,:), intent(in) :: pair_atoms
    real(real32), dimension(:,:), intent(in) :: pair_shift
    real(real32), dimension(:), intent(in) :: pair_distance
    integer, intent(in) :: pair_index
    real(real32), dimension(:), intent(out) :: features

    integer :: atom_i, atom_j, offset
    real(real32), dimension(3) :: midpoint, unit_vec
    real(real32) :: coord_scale

    coord_scale = max(this%bond_cutoff, 1.E-6_real32)
    atom_i = pair_atoms(1, pair_index)
    atom_j = pair_atoms(2, pair_index)
    midpoint = 0.5_real32 * (positions(:, atom_i) + positions(:, atom_j))
    unit_vec = 0._real32
    if (pair_distance(pair_index) > 1.E-6_real32) then
       unit_vec = pair_shift(:, pair_index) / pair_distance(pair_index)
    end if

    features = 0._real32
    features(1:3) = midpoint / coord_scale
    features(4:6) = unit_vec
    features(7) = pair_distance(pair_index) / coord_scale
    offset = 7
    if (species_index(atom_i) > 0) features(offset + species_index(atom_i)) = 1._real32
    offset = offset + this%num_species
    if (species_index(atom_j) > 0) features(offset + species_index(atom_j)) = 1._real32

  end subroutine build_pair_vertex_features
!###############################################################################


!###############################################################################
  subroutine build_triplet_vertex_features(this, basis, positions, species_index, &
       triplet_atoms, triplet_index, features)
    !! Construct deterministic triplet-vertex features for the 4-body graph.
    implicit none

    class(gnn_fingerprint_type), intent(in) :: this
    type(basis_type), intent(in) :: basis
    real(real32), dimension(:,:), intent(in) :: positions
    integer, dimension(:), intent(in) :: species_index
    integer, dimension(:,:), intent(in) :: triplet_atoms
    integer, intent(in) :: triplet_index
    real(real32), dimension(:), intent(out) :: features

    integer :: atom_i, atom_j, atom_k, offset
    real(real32), dimension(3) :: rij, rjk, centroid
    real(real32) :: dij, djk, angle_ijk, coord_scale

    coord_scale = max(this%bond_cutoff, 1.E-6_real32)
    atom_i = triplet_atoms(1, triplet_index)
    atom_j = triplet_atoms(2, triplet_index)
    atom_k = triplet_atoms(3, triplet_index)

    rij = minimum_image_delta(basis, positions(:, atom_i) - positions(:, atom_j))
    rjk = minimum_image_delta(basis, positions(:, atom_k) - positions(:, atom_j))
    dij = sqrt(sum(rij**2))
    djk = sqrt(sum(rjk**2))
    angle_ijk = get_angle(rij, rjk)
    centroid = (positions(:, atom_i) + positions(:, atom_j) + &
         positions(:, atom_k)) / 3._real32

    features = 0._real32
    features(1:3) = centroid / coord_scale
    features(4) = dij / coord_scale
    features(5) = djk / coord_scale
    features(6) = angle_ijk / pi
    offset = 6
    if (species_index(atom_i) > 0) features(offset + species_index(atom_i)) = 1._real32
    offset = offset + this%num_species
    if (species_index(atom_j) > 0) features(offset + species_index(atom_j)) = 1._real32
    offset = offset + this%num_species
    if (species_index(atom_k) > 0) features(offset + species_index(atom_k)) = 1._real32

  end subroutine build_triplet_vertex_features
!###############################################################################


!###############################################################################
  subroutine find_shared_pair_atom(pair_a, pair_b, shared_atom, other_a, other_b)
    !! Identify the shared atom between two unique pairs, if one exists.
    implicit none

    integer, dimension(2), intent(in) :: pair_a
    integer, dimension(2), intent(in) :: pair_b
    integer, intent(out) :: shared_atom
    integer, intent(out) :: other_a
    integer, intent(out) :: other_b

    shared_atom = 0
    other_a = 0
    other_b = 0

    if (pair_a(1) == pair_b(1)) then
       shared_atom = pair_a(1)
       other_a = pair_a(2)
       other_b = pair_b(2)
    elseif (pair_a(1) == pair_b(2)) then
       shared_atom = pair_a(1)
       other_a = pair_a(2)
       other_b = pair_b(1)
    elseif (pair_a(2) == pair_b(1)) then
       shared_atom = pair_a(2)
       other_a = pair_a(1)
       other_b = pair_b(2)
    elseif (pair_a(2) == pair_b(2)) then
       shared_atom = pair_a(2)
       other_a = pair_a(1)
       other_b = pair_b(1)
    end if

  end subroutine find_shared_pair_atom
!###############################################################################


!###############################################################################
  subroutine get_branch_layout(this, branch_id, offset, dim)
    !! Return the flat fingerprint offset and size for one branch.
    implicit none

    class(gnn_fingerprint_type), intent(in) :: this
    integer, intent(in) :: branch_id
    integer, intent(out) :: offset
    integer, intent(out) :: dim

    select case (branch_id)
    case (1)
       offset = 0
       dim = this%fingerprint_dim_2body
    case (2)
       offset = this%fingerprint_dim_2body
       dim = this%fingerprint_dim_3body
    case (3)
       offset = this%fingerprint_dim_2body + this%fingerprint_dim_3body
       dim = this%fingerprint_dim_4body
    case default
       call stop_program('gnn_fingerprint: invalid branch id')
    end select

  end subroutine get_branch_layout
   !###############################################################################


   !###############################################################################
  subroutine build_triplet_data(basis, pair_atoms, pair_distance, triplet_atoms)
    !! Build the ordered triplet-vertex list used by the 4-body graph.
    implicit none

    type(basis_type), intent(in) :: basis
    integer, dimension(:,:), intent(in) :: pair_atoms
    real(real32), dimension(:), intent(in) :: pair_distance
    integer, dimension(:,:), allocatable, intent(out) :: triplet_atoms

    integer :: atom_j, pair_i, pair_k, atom_i, atom_k, num_triplets

    num_triplets = 0
    do atom_j = 1, basis%natom
       do pair_i = 1, size(pair_distance)
          if (all(pair_atoms(:, pair_i) /= atom_j)) cycle
          atom_i = pair_atoms(1, pair_i)
          if (atom_i == atom_j) atom_i = pair_atoms(2, pair_i)
          do pair_k = 1, size(pair_distance)
             if (pair_k == pair_i) cycle
             if (all(pair_atoms(:, pair_k) /= atom_j)) cycle
             atom_k = pair_atoms(1, pair_k)
             if (atom_k == atom_j) atom_k = pair_atoms(2, pair_k)
             if (atom_i == atom_k) cycle
             num_triplets = num_triplets + 1
          end do
       end do
    end do

    allocate(triplet_atoms(3, num_triplets))

    num_triplets = 0
    do atom_j = 1, basis%natom
       do pair_i = 1, size(pair_distance)
          if (all(pair_atoms(:, pair_i) /= atom_j)) cycle
          atom_i = pair_atoms(1, pair_i)
          if (atom_i == atom_j) atom_i = pair_atoms(2, pair_i)
          do pair_k = 1, size(pair_distance)
             if (pair_k == pair_i) cycle
             if (all(pair_atoms(:, pair_k) /= atom_j)) cycle
             atom_k = pair_atoms(1, pair_k)
             if (atom_k == atom_j) atom_k = pair_atoms(2, pair_k)
             if (atom_i == atom_k) cycle
             num_triplets = num_triplets + 1
             triplet_atoms(:, num_triplets) = [atom_i, atom_j, atom_k]
          end do
       end do
    end do

  end subroutine build_triplet_data
   !###############################################################################


   !###############################################################################
  subroutine build_angle_triplet_metadata( &
       basis, positions, pair_atoms, pair_distance, angle_triplets)
    !! Store the atom triplets that define each 3-body graph edge.
    implicit none

    type(basis_type), intent(in) :: basis
    real(real32), dimension(:,:), intent(in) :: positions
    integer, dimension(:,:), intent(in) :: pair_atoms
    real(real32), dimension(:), intent(in) :: pair_distance
    integer, dimension(:,:), allocatable, intent(out) :: angle_triplets

    integer :: pair_i, pair_j, shared_atom, other_i, other_j, num_angles
    real(real32), dimension(3) :: vec_i, vec_j
    real(real32) :: angle_ijk

    num_angles = 0
    do pair_i = 1, size(pair_distance) - 1
       do pair_j = pair_i + 1, size(pair_distance)
          call find_shared_pair_atom(pair_atoms(:, pair_i), pair_atoms(:, pair_j), &
               shared_atom, other_i, other_j)
          if (shared_atom == 0) cycle
          vec_i = minimum_image_delta(basis, positions(:, other_i) - &
               positions(:, shared_atom))
          vec_j = minimum_image_delta(basis, positions(:, other_j) - &
               positions(:, shared_atom))
          if (sqrt(sum(vec_i**2)) <= 1.E-6_real32) cycle
          if (sqrt(sum(vec_j**2)) <= 1.E-6_real32) cycle
          angle_ijk = get_angle(vec_i, vec_j)
          if (.not. ieee_is_finite(angle_ijk)) cycle
          num_angles = num_angles + 1
       end do
    end do

    allocate(angle_triplets(3, num_angles))

    num_angles = 0
    do pair_i = 1, size(pair_distance) - 1
       do pair_j = pair_i + 1, size(pair_distance)
          call find_shared_pair_atom(pair_atoms(:, pair_i), pair_atoms(:, pair_j), &
               shared_atom, other_i, other_j)
          if (shared_atom == 0) cycle
          vec_i = minimum_image_delta(basis, positions(:, other_i) - &
               positions(:, shared_atom))
          vec_j = minimum_image_delta(basis, positions(:, other_j) - &
               positions(:, shared_atom))
          if (sqrt(sum(vec_i**2)) <= 1.E-6_real32) cycle
          if (sqrt(sum(vec_j**2)) <= 1.E-6_real32) cycle
          angle_ijk = get_angle(vec_i, vec_j)
          if (.not. ieee_is_finite(angle_ijk)) cycle
          num_angles = num_angles + 1
          angle_triplets(:, num_angles) = [other_i, shared_atom, other_j]
       end do
    end do

  end subroutine build_angle_triplet_metadata
   !###############################################################################


   !###############################################################################
  subroutine build_improper_quad_metadata( &
       basis, positions, triplet_atoms, improper_quads)
    !! Store the atom quadruplets that define each 4-body graph edge.
    implicit none

    type(basis_type), intent(in) :: basis
    real(real32), dimension(:,:), intent(in) :: positions
    integer, dimension(:,:), intent(in) :: triplet_atoms
    integer, dimension(:,:), allocatable, intent(out) :: improper_quads

    integer :: triplet_i, triplet_j, atom_i, atom_j, atom_k, atom_l, num_dihedrals
    real(real32) :: dihedral_angle

    num_dihedrals = 0
    do triplet_i = 1, size(triplet_atoms, 2)
       do triplet_j = 1, size(triplet_atoms, 2)
          if (triplet_i == triplet_j) cycle
          if (triplet_atoms(2, triplet_i) /= triplet_atoms(1, triplet_j)) cycle
          if (triplet_atoms(3, triplet_i) /= triplet_atoms(2, triplet_j)) cycle
          atom_i = triplet_atoms(1, triplet_i)
          atom_j = triplet_atoms(2, triplet_i)
          atom_k = triplet_atoms(3, triplet_i)
          atom_l = triplet_atoms(3, triplet_j)
          if (atom_i == atom_l) cycle
          dihedral_angle = get_improper_dihedral_angle(positions(:, atom_i), &
               positions(:, atom_j), positions(:, atom_k), positions(:, atom_l))
          if (.not. ieee_is_finite(dihedral_angle)) cycle
          num_dihedrals = num_dihedrals + 1
       end do
    end do

    allocate(improper_quads(4, num_dihedrals))

    num_dihedrals = 0
    do triplet_i = 1, size(triplet_atoms, 2)
       do triplet_j = 1, size(triplet_atoms, 2)
          if (triplet_i == triplet_j) cycle
          if (triplet_atoms(2, triplet_i) /= triplet_atoms(1, triplet_j)) cycle
          if (triplet_atoms(3, triplet_i) /= triplet_atoms(2, triplet_j)) cycle
          atom_i = triplet_atoms(1, triplet_i)
          atom_j = triplet_atoms(2, triplet_i)
          atom_k = triplet_atoms(3, triplet_i)
          atom_l = triplet_atoms(3, triplet_j)
          if (atom_i == atom_l) cycle
          dihedral_angle = get_improper_dihedral_angle(positions(:, atom_i), &
               positions(:, atom_j), positions(:, atom_k), positions(:, atom_l))
          if (.not. ieee_is_finite(dihedral_angle)) cycle
          num_dihedrals = num_dihedrals + 1
          improper_quads(:, num_dihedrals) = [atom_i, atom_j, atom_k, atom_l]
       end do
    end do

  end subroutine build_improper_quad_metadata
   !###############################################################################


   !###############################################################################
  subroutine prepare_multigraph_data(this, basis, data)
    !! Prepare deterministic graph metadata for exact chain-rule gradients.
    implicit none

    class(gnn_fingerprint_type), intent(in) :: this
    type(basis_type), intent(in) :: basis
    type(multigraph_data_type), intent(out) :: data

    call extract_structure_data(this, basis, data%positions, data%species_index, &
         data%atomic_numbers, data%cov_radii)
    call build_pair_data( &
         this, basis, data%positions, data%pair_atoms, data%pair_shift, &
         data%pair_distance)
    call build_triplet_data( &
         basis, data%pair_atoms, data%pair_distance, data%triplet_atoms)
    call build_angle_triplet_metadata(basis, data%positions, data%pair_atoms, &
         data%pair_distance, data%angle_triplets)
    call build_improper_quad_metadata(basis, data%positions, data%triplet_atoms, &
         data%improper_quads)

  end subroutine prepare_multigraph_data
   !###############################################################################


   !###############################################################################
  subroutine fill_branch_affine_arrays(this, branch_id, scale_factor, offset)
    !! Build the affine transform that maps raw network outputs to physical fingerprints.
    implicit none

    class(gnn_fingerprint_type), intent(in) :: this
    integer, intent(in) :: branch_id
    real(real32), dimension(:,:), allocatable, intent(out) :: scale_factor
    real(real32), dimension(:,:), allocatable, intent(out) :: offset

    integer :: branch_offset, branch_dim, i

    call get_branch_layout(this, branch_id, branch_offset, branch_dim)
    allocate(scale_factor(branch_dim, 1), source = 1._real32)
    allocate(offset(branch_dim, 1), source = 0._real32)

    if (allocated(this%target_comp_weights)) then
       do i = 1, branch_dim
          if (this%target_comp_weights(branch_offset + i) > 1.E-12_real32) then
             scale_factor(i, 1) = scale_factor(i, 1) / &
                  this%target_comp_weights(branch_offset + i)
          else
             scale_factor(i, 1) = 0._real32
          end if
       end do
    end if

    if (allocated(this%target_scale)) then
       scale_factor(:, 1) = scale_factor(:, 1) * &
            this%target_scale(branch_offset + 1:branch_offset + branch_dim)
    end if
    if (allocated(this%target_mean)) then
       offset(:, 1) = this%target_mean(branch_offset + 1:branch_offset + branch_dim)
    end if

  end subroutine fill_branch_affine_arrays
   !###############################################################################


   !###############################################################################
  subroutine accumulate_distance_gradient( &
       grad_scalar, delta, distance, atom_i, atom_j, grad_positions)
    !! Accumulate the position gradient contribution from one distance feature.
    implicit none

    real(real32), intent(in) :: grad_scalar
    real(real32), dimension(3), intent(in) :: delta
    real(real32), intent(in) :: distance
    integer, intent(in) :: atom_i
    integer, intent(in) :: atom_j
    real(real32), dimension(:,:), intent(inout) :: grad_positions

    real(real32), dimension(3) :: direction

    if (distance <= 1.E-6_real32) return
    direction = delta / distance
    grad_positions(atom_i, :) = grad_positions(atom_i, :) - grad_scalar * direction
    grad_positions(atom_j, :) = grad_positions(atom_j, :) + grad_scalar * direction

  end subroutine accumulate_distance_gradient
   !###############################################################################


   !###############################################################################
  subroutine apply_unit_vector_gradient(grad_unit, delta, distance, atom_i, atom_j, &
       grad_positions)
    !! Accumulate gradients for a normalised pair-direction feature.
    implicit none

    real(real32), dimension(3), intent(in) :: grad_unit
    real(real32), dimension(3), intent(in) :: delta
    real(real32), intent(in) :: distance
    integer, intent(in) :: atom_i
    integer, intent(in) :: atom_j
    real(real32), dimension(:,:), intent(inout) :: grad_positions

    real(real32), dimension(3) :: unit_vec, grad_delta
    real(real32), dimension(3,3) :: projector
    integer :: i, j

    if (distance <= 1.E-6_real32) return

    unit_vec = delta / distance
    projector = 0._real32
    do i = 1, 3
       projector(i, i) = 1._real32
       do j = 1, 3
          projector(i, j) = projector(i, j) - unit_vec(i) * unit_vec(j)
       end do
    end do
    grad_delta = matmul(projector, grad_unit) / distance

    grad_positions(atom_i, :) = grad_positions(atom_i, :) - grad_delta
    grad_positions(atom_j, :) = grad_positions(atom_j, :) + grad_delta

  end subroutine apply_unit_vector_gradient
   !###############################################################################


   !###############################################################################
  subroutine get_angle_vector_gradients(vec_a, vec_b, grad_a, grad_b)
    !! Return d angle(vec_a, vec_b) / d vec_a and d vec_b.
    implicit none

    real(real32), dimension(3), intent(in) :: vec_a
    real(real32), dimension(3), intent(in) :: vec_b
    real(real32), dimension(3), intent(out) :: grad_a
    real(real32), dimension(3), intent(out) :: grad_b

    real(real32) :: norm_a, norm_b, cos_theta, sin_theta

    grad_a = 0._real32
    grad_b = 0._real32

    norm_a = sqrt(sum(vec_a**2))
    norm_b = sqrt(sum(vec_b**2))
    if (norm_a <= 1.E-6_real32 .or. norm_b <= 1.E-6_real32) return

    cos_theta = dot_product(vec_a, vec_b) / (norm_a * norm_b)
    cos_theta = min(1._real32 - 1.E-6_real32, max(-1._real32 + 1.E-6_real32, cos_theta))
    sin_theta = sqrt(max(1._real32 - cos_theta**2, 1.E-12_real32))

    grad_a = -(vec_b / (norm_a * norm_b) - cos_theta * vec_a / (norm_a**2)) / sin_theta
    grad_b = -(vec_a / (norm_a * norm_b) - cos_theta * vec_b / (norm_b**2)) / sin_theta

  end subroutine get_angle_vector_gradients
   !###############################################################################


   !###############################################################################
  subroutine accumulate_angle_gradient(grad_scalar, vec_a, vec_b, atom_a, atom_center, &
       atom_b, grad_positions)
    !! Accumulate the position gradient for one bond-angle feature.
    implicit none

    real(real32), intent(in) :: grad_scalar
    real(real32), dimension(3), intent(in) :: vec_a
    real(real32), dimension(3), intent(in) :: vec_b
    integer, intent(in) :: atom_a
    integer, intent(in) :: atom_center
    integer, intent(in) :: atom_b
    real(real32), dimension(:,:), intent(inout) :: grad_positions

    real(real32), dimension(3) :: grad_vec_a, grad_vec_b

    call get_angle_vector_gradients(vec_a, vec_b, grad_vec_a, grad_vec_b)
    grad_positions(atom_a, :) = grad_positions(atom_a, :) + grad_scalar * grad_vec_a
    grad_positions(atom_center, :) = grad_positions(atom_center, :) - &
         grad_scalar * (grad_vec_a + grad_vec_b)
    grad_positions(atom_b, :) = grad_positions(atom_b, :) + grad_scalar * grad_vec_b

  end subroutine accumulate_angle_gradient
   !###############################################################################


   !###############################################################################
  subroutine accumulate_improper_dihedral_gradient( &
       grad_scalar, point1, point2, point3, &
       point4, atom1, atom2, atom3, atom4, grad_positions)
    !! Accumulate the position gradient for one improper-dihedral feature.
    implicit none

    real(real32), intent(in) :: grad_scalar
    real(real32), dimension(3), intent(in) :: point1
    real(real32), dimension(3), intent(in) :: point2
    real(real32), dimension(3), intent(in) :: point3
    real(real32), dimension(3), intent(in) :: point4
    integer, intent(in) :: atom1
    integer, intent(in) :: atom2
    integer, intent(in) :: atom3
    integer, intent(in) :: atom4
    real(real32), dimension(:,:), intent(inout) :: grad_positions

    real(real32), dimension(3) :: vec_a, vec_b, vec_c
    real(real32), dimension(3) :: normal_1, normal_2
    real(real32), dimension(3) :: grad_normal_1, grad_normal_2
    real(real32), dimension(3) :: grad_a, grad_b, grad_c
    real(real32), dimension(3) :: &
         grad_point_1, grad_point_2, grad_point_3, grad_point_4

    vec_a = point2 - point1
    vec_b = point3 - point1
    vec_c = point4 - point1
    normal_1 = cross(vec_a, vec_b)
    normal_2 = cross(vec_b, vec_c)

    if (sqrt(sum(normal_1**2)) <= 1.E-6_real32) return
    if (sqrt(sum(normal_2**2)) <= 1.E-6_real32) return

    call get_angle_vector_gradients(normal_1, normal_2, grad_normal_1, grad_normal_2)

    grad_a = cross(vec_b, grad_normal_1)
    grad_b = cross(grad_normal_1, vec_a) + cross(vec_c, grad_normal_2)
    grad_c = cross(grad_normal_2, vec_b)

    grad_point_1 = -(grad_a + grad_b + grad_c)
    grad_point_2 = grad_a
    grad_point_3 = grad_b
    grad_point_4 = grad_c

    grad_positions(atom1, :) = grad_positions(atom1, :) + grad_scalar * grad_point_1
    grad_positions(atom2, :) = grad_positions(atom2, :) + grad_scalar * grad_point_2
    grad_positions(atom3, :) = grad_positions(atom3, :) + grad_scalar * grad_point_3
    grad_positions(atom4, :) = grad_positions(atom4, :) + grad_scalar * grad_point_4

  end subroutine accumulate_improper_dihedral_gradient
   !###############################################################################


   !###############################################################################
  subroutine accumulate_2body_vertex_gradients(this, vertex_grad, grad_positions)
    !! Chain vertex-feature gradients on the 2-body graph back to positions.
    implicit none

    class(gnn_fingerprint_type), intent(in) :: this
    real(real32), dimension(:,:), intent(in) :: vertex_grad
    real(real32), dimension(:,:), intent(inout) :: grad_positions

    integer :: atom_i
    real(real32) :: coord_scale

    coord_scale = max(this%bond_cutoff, 1.E-6_real32)
    do atom_i = 1, min(size(vertex_grad, 2), size(grad_positions, 1))
       grad_positions(atom_i, :) = &
            grad_positions(atom_i, :) + vertex_grad(1:3, atom_i) / coord_scale
    end do

  end subroutine accumulate_2body_vertex_gradients
   !###############################################################################


   !###############################################################################
  subroutine accumulate_2body_edge_gradients(edge_grad, data, grad_positions)
    !! Chain edge-feature gradients on the 2-body graph back to positions.
    implicit none

    real(real32), dimension(:,:), intent(in) :: edge_grad
    type(multigraph_data_type), intent(in) :: data
    real(real32), dimension(:,:), intent(inout) :: grad_positions

    integer :: edge_i, num_edges

    num_edges = min(size(edge_grad, 2), size(data%pair_distance))
    do edge_i = 1, num_edges
       call accumulate_distance_gradient( &
            edge_grad(1, edge_i), data%pair_shift(:, edge_i), &
            data%pair_distance(edge_i), data%pair_atoms(1, edge_i), &
            data%pair_atoms(2, edge_i), grad_positions)
    end do

  end subroutine accumulate_2body_edge_gradients
   !###############################################################################


   !###############################################################################
  subroutine accumulate_3body_vertex_gradients(this, data, vertex_grad, grad_positions)
    !! Chain pair-vertex gradients on the 3-body graph back to positions.
    implicit none

    class(gnn_fingerprint_type), intent(in) :: this
    type(multigraph_data_type), intent(in) :: data
    real(real32), dimension(:,:), intent(in) :: vertex_grad
    real(real32), dimension(:,:), intent(inout) :: grad_positions

    integer :: pair_i, atom_i, atom_j, num_pairs
    real(real32) :: coord_scale

    coord_scale = max(this%bond_cutoff, 1.E-6_real32)
    num_pairs = min(size(vertex_grad, 2), size(data%pair_distance))
    do pair_i = 1, num_pairs
       atom_i = data%pair_atoms(1, pair_i)
       atom_j = data%pair_atoms(2, pair_i)
       grad_positions(atom_i, :) = grad_positions(atom_i, :) + 0.5_real32 * &
            vertex_grad(1:3, pair_i) / coord_scale
       grad_positions(atom_j, :) = grad_positions(atom_j, :) + 0.5_real32 * &
            vertex_grad(1:3, pair_i) / coord_scale
       call apply_unit_vector_gradient( &
            vertex_grad(4:6, pair_i), data%pair_shift(:, pair_i), &
            data%pair_distance(pair_i), atom_i, atom_j, grad_positions)
       call accumulate_distance_gradient( &
            vertex_grad(7, pair_i) / coord_scale, &
            data%pair_shift(:, pair_i), data%pair_distance(pair_i), atom_i, atom_j, &
            grad_positions)
    end do

  end subroutine accumulate_3body_vertex_gradients
!###############################################################################


!###############################################################################
  subroutine accumulate_3body_edge_gradients(basis, data, edge_grad, grad_positions)
    !! Chain angle-edge gradients on the 3-body graph back to positions.
    implicit none

    type(basis_type), intent(in) :: basis
    type(multigraph_data_type), intent(in) :: data
    real(real32), dimension(:,:), intent(in) :: edge_grad
    real(real32), dimension(:,:), intent(inout) :: grad_positions

    integer :: edge_i, atom_i, atom_j, atom_k, num_edges
    real(real32), dimension(3) :: vec_i, vec_k

    num_edges = min(size(edge_grad, 2), size(data%angle_triplets, 2))
    do edge_i = 1, num_edges
       atom_i = data%angle_triplets(1, edge_i)
       atom_j = data%angle_triplets(2, edge_i)
       atom_k = data%angle_triplets(3, edge_i)
       vec_i = minimum_image_delta(basis, data%positions(:, atom_i) - &
            data%positions(:, atom_j))
       vec_k = minimum_image_delta(basis, data%positions(:, atom_k) - &
            data%positions(:, atom_j))
       call accumulate_angle_gradient( &
            edge_grad(1, edge_i), vec_i, vec_k, atom_i, atom_j, &
            atom_k, grad_positions)
    end do

  end subroutine accumulate_3body_edge_gradients
!###############################################################################


!###############################################################################
  subroutine accumulate_4body_vertex_gradients(this, basis, data, vertex_grad, &
       grad_positions)
    !! Chain triplet-vertex gradients on the 4-body graph back to positions.
    implicit none

    class(gnn_fingerprint_type), intent(in) :: this
    type(basis_type), intent(in) :: basis
    type(multigraph_data_type), intent(in) :: data
    real(real32), dimension(:,:), intent(in) :: vertex_grad
    real(real32), dimension(:,:), intent(inout) :: grad_positions

    integer :: triplet_i, atom_i, atom_j, atom_k, num_triplets
    real(real32) :: coord_scale
    real(real32), dimension(3) :: vec_ij, vec_kj

    coord_scale = max(this%bond_cutoff, 1.E-6_real32)
    num_triplets = min(size(vertex_grad, 2), size(data%triplet_atoms, 2))
    do triplet_i = 1, num_triplets
       atom_i = data%triplet_atoms(1, triplet_i)
       atom_j = data%triplet_atoms(2, triplet_i)
       atom_k = data%triplet_atoms(3, triplet_i)
       grad_positions(atom_i, :) = grad_positions(atom_i, :) + &
            vertex_grad(1:3, triplet_i) / (3._real32 * coord_scale)
       grad_positions(atom_j, :) = grad_positions(atom_j, :) + &
            vertex_grad(1:3, triplet_i) / (3._real32 * coord_scale)
       grad_positions(atom_k, :) = grad_positions(atom_k, :) + &
            vertex_grad(1:3, triplet_i) / (3._real32 * coord_scale)
       vec_ij = minimum_image_delta(basis, data%positions(:, atom_i) - &
            data%positions(:, atom_j))
       vec_kj = minimum_image_delta(basis, data%positions(:, atom_k) - &
            data%positions(:, atom_j))
       call accumulate_distance_gradient( &
            vertex_grad(4, triplet_i) / coord_scale, vec_ij, &
            sqrt(sum(vec_ij**2)), atom_i, atom_j, grad_positions)
       call accumulate_distance_gradient( &
            vertex_grad(5, triplet_i) / coord_scale, vec_kj, &
            sqrt(sum(vec_kj**2)), atom_j, atom_k, grad_positions)
       call accumulate_angle_gradient( &
            vertex_grad(6, triplet_i) / pi, vec_ij, vec_kj, &
            atom_i, atom_j, atom_k, grad_positions)
    end do

  end subroutine accumulate_4body_vertex_gradients
!###############################################################################


!###############################################################################
  subroutine accumulate_4body_edge_gradients(data, edge_grad, grad_positions)
    !! Chain improper-dihedral edge gradients on the 4-body graph back to positions.
    implicit none

    type(multigraph_data_type), intent(in) :: data
    real(real32), dimension(:,:), intent(in) :: edge_grad
    real(real32), dimension(:,:), intent(inout) :: grad_positions

    integer :: edge_i, atom_i, atom_j, atom_k, atom_l, num_edges

    num_edges = min(size(edge_grad, 2), size(data%improper_quads, 2))
    do edge_i = 1, num_edges
       atom_i = data%improper_quads(1, edge_i)
       atom_j = data%improper_quads(2, edge_i)
       atom_k = data%improper_quads(3, edge_i)
       atom_l = data%improper_quads(4, edge_i)
       call accumulate_improper_dihedral_gradient(edge_grad(1, edge_i), &
            data%positions(:, atom_i), data%positions(:, atom_j), &
            data%positions(:, atom_k), data%positions(:, atom_l), atom_i, atom_j, &
            atom_k, atom_l, grad_positions)
    end do

  end subroutine accumulate_4body_edge_gradients
!###############################################################################


!###############################################################################
  subroutine accumulate_branch_feature_gradients( &
       this, branch_id, basis, data, vertex_grad, edge_grad, grad_positions)
    !! Push branch input-feature gradients back to Cartesian positions.
    implicit none

    class(gnn_fingerprint_type), intent(in) :: this
    integer, intent(in) :: branch_id
    type(basis_type), intent(in) :: basis
    type(multigraph_data_type), intent(in) :: data
    real(real32), dimension(:,:), intent(in) :: vertex_grad
    real(real32), dimension(:,:), intent(in) :: edge_grad
    real(real32), dimension(:,:), intent(inout) :: grad_positions

    select case (branch_id)
    case (1)
       call accumulate_2body_vertex_gradients(this, vertex_grad, grad_positions)
       if (size(edge_grad, 2) > 0) &
            call accumulate_2body_edge_gradients(edge_grad, data, grad_positions)
    case (2)
       call accumulate_3body_vertex_gradients(this, data, vertex_grad, grad_positions)
       if (size(edge_grad, 2) > 0) &
            call accumulate_3body_edge_gradients(basis, data, edge_grad, grad_positions)
    case (3)
       call accumulate_4body_vertex_gradients( &
            this, basis, data, vertex_grad, grad_positions)
       if (size(edge_grad, 2) > 0) &
            call accumulate_4body_edge_gradients(data, edge_grad, grad_positions)
    case default
       call stop_program('gnn_fingerprint: invalid branch id for feature gradients')
    end select

  end subroutine accumulate_branch_feature_gradients
!###############################################################################


!###############################################################################
  subroutine compute_branch_output_jacobian( &
       this, network, branch_id, graph, basis, data, grad_block)
    !! Compute exact position Jacobians for one branch output vector.
    implicit none

    class(gnn_fingerprint_type), intent(inout) :: this
    type(network_type), intent(inout) :: network
    integer, intent(in) :: branch_id
    type(graph_type), dimension(1,1), intent(in) :: graph
    type(basis_type), intent(in) :: basis
    type(multigraph_data_type), intent(in) :: data
    real(real32), dimension(:,:,:), intent(out) :: grad_block

    integer :: component_i, root_id, leaf_id
    real(real32), dimension(:,:), allocatable :: scale_factor, offset
    real(real32), dimension(:,:), allocatable :: target_array
    real(real32), dimension(:,:), allocatable :: grad_positions, vertex_grad, edge_grad
    type(array_type), pointer :: loss_array
    real(real32) :: branch_dim_real

    call fill_branch_affine_arrays(this, branch_id, scale_factor, offset)
    allocate(target_array(size(scale_factor, 1), 1))
    allocate(grad_positions(size(data%positions, 2), 3))
    grad_block = 0._real32
    branch_dim_real = real(max(size(scale_factor, 1), 1), real32)

    do component_i = 1, size(scale_factor, 1)
       call network%set_batch_size(1)
       call network%set_inference_mode()
       call network%forward(graph, input_requires_grad = .true.)
       root_id = network%auto_graph%vertex(network%root_vertices(1))%id
       leaf_id = network%auto_graph%vertex(network%leaf_vertices(1))%id

       target_array(:, 1) = network%model(leaf_id)%layer%output(1,1)%val(:, 1)
       target_array(component_i, 1) = target_array(component_i, 1) - &
            0.5_real32 * branch_dim_real * scale_factor(component_i, 1)
       call network%save_output(target_array)
       loss_array => network%loss_eval(1, 1)
       call loss_array%grad_reverse()

       grad_positions = 0._real32
       if (.not. associated(network%model(root_id)%layer%output(1,1)%grad)) then
          call stop_program( &
               'gnn_fingerprint: missing vertex gradients for branch jacobian')
       end if
       vertex_grad = network%model(root_id)%layer%output(1,1)%grad%val
       allocate(edge_grad(1, 0))
       if (associated(network%model(root_id)%layer%output(2,1)%grad)) then
          deallocate(edge_grad)
          edge_grad = network%model(root_id)%layer%output(2,1)%grad%val
       end if
       call accumulate_branch_feature_gradients( &
            this, branch_id, basis, data, vertex_grad, &
            edge_grad, grad_positions)
       grad_block(component_i, :, :) = grad_positions

       call loss_array%nullify_graph()
       deallocate(loss_array)
       nullify(loss_array)
       call network%reset_gradients()
       if (allocated(edge_grad)) deallocate(edge_grad)
    end do

  end subroutine compute_branch_output_jacobian
   !###############################################################################


   !###############################################################################
  subroutine compute_branch_loss_and_gradient( &
       this, network, branch_id, graph, basis, data, &
       target_component, current_component, grad_positions, branch_loss)
    !! Evaluate one branch loss in physical space and accumulate its exact gradient.
    implicit none

    class(gnn_fingerprint_type), intent(inout) :: this
    type(network_type), intent(inout) :: network
    integer, intent(in) :: branch_id
    type(graph_type), dimension(1,1), intent(in) :: graph
    type(basis_type), intent(in) :: basis
    type(multigraph_data_type), intent(in) :: data
    real(real32), dimension(:), intent(in) :: target_component
    real(real32), dimension(:), intent(out) :: current_component
    real(real32), dimension(:,:), intent(inout) :: grad_positions
    real(real32), intent(out) :: branch_loss

    integer :: component_i, root_id, leaf_id
    real(real32), dimension(:,:), allocatable :: scale_factor, offset
    logical, dimension(:,:), allocatable :: mask
    real(real32), dimension(:,:), allocatable :: vertex_grad, edge_grad
    type(array_type), pointer :: projected_output, residual_array, loss_array
    real(real32) :: branch_weight

    call fill_branch_affine_arrays(this, branch_id, scale_factor, offset)
    allocate(mask(size(scale_factor, 1), 1), source = .false.)

    branch_weight = this%component_weight(branch_id) / &
         real(max(size(scale_factor, 1), 1), real32)
    current_component = 0._real32
    branch_loss = 0._real32

    do component_i = 1, size(scale_factor, 1)
       mask = .false.
       mask(component_i, 1) = .true.

       call network%set_batch_size(1)
       call network%set_inference_mode()
       call network%forward(graph, input_requires_grad = .true.)
       root_id = network%auto_graph%vertex(network%root_vertices(1))%id
       leaf_id = network%auto_graph%vertex(network%leaf_vertices(1))%id

       current_component(component_i) = &
            network%model(leaf_id)%layer%output(1,1)%val(component_i, 1) * &
            scale_factor(component_i, 1) + offset(component_i, 1)

       projected_output => &
            ds_sum(network%model(leaf_id)%layer%output(1,1) * mask, dim = 1) * &
            scale_factor(component_i, 1) + offset(component_i, 1)
       residual_array => projected_output - target_component(component_i)
       loss_array => (residual_array * residual_array) * branch_weight
       call loss_array%grad_reverse(reset_graph = .true.)

       if (.not. associated(network%model(root_id)%layer%output(1,1)%grad)) then
          call stop_program('gnn_fingerprint: missing vertex gradients for branch loss')
       end if
       vertex_grad = network%model(root_id)%layer%output(1,1)%grad%val
       allocate(edge_grad(1, 0))
       if (associated(network%model(root_id)%layer%output(2,1)%grad)) then
          deallocate(edge_grad)
          edge_grad = network%model(root_id)%layer%output(2,1)%grad%val
       end if
       call accumulate_branch_feature_gradients( &
            this, branch_id, basis, data, vertex_grad, &
            edge_grad, grad_positions)
       branch_loss = branch_loss + loss_array%val(1, 1)

       call loss_array%nullify_graph(ignore_ownership = .false.)
       deallocate(loss_array)
       nullify(projected_output)
       nullify(residual_array)
       nullify(loss_array)
       if (allocated(edge_grad)) deallocate(edge_grad)
    end do

  end subroutine compute_branch_loss_and_gradient
   !###############################################################################


   !###############################################################################
  subroutine evaluate_predictive_loss_and_gradients(this, basis, target_2body, &
       target_3body, target_4body, current_2body, current_3body, current_4body, &
       grad_positions, loss_value)
    !! Evaluate learned component losses and exact gradients with respect to positions.
    implicit none

    class(gnn_fingerprint_type), intent(inout) :: this
    type(basis_type), intent(in) :: basis
    real(real32), dimension(:), intent(in) :: target_2body
    real(real32), dimension(:), intent(in) :: target_3body
    real(real32), dimension(:), intent(in) :: target_4body
    real(real32), dimension(:), intent(out) :: current_2body
    real(real32), dimension(:), intent(out) :: current_3body
    real(real32), dimension(:), intent(out) :: current_4body
    real(real32), dimension(:,:), intent(out) :: grad_positions
    real(real32), intent(out) :: loss_value

    real(real32), dimension(:,:,:), allocatable :: grad_2body, grad_3body, grad_4body
    integer :: atom_i, coord

    call evaluate_components( &
         this, basis, .true., current_2body, current_3body, current_4body)

    allocate(grad_2body(size(target_2body), basis%natom, 3))
    allocate(grad_3body(size(target_3body), basis%natom, 3))
    allocate(grad_4body(size(target_4body), basis%natom, 3))
    call this%compute_gradients(basis, grad_2body, grad_3body, grad_4body)

    grad_positions = 0._real32
    loss_value = this%component_weight(1) * sum((current_2body - target_2body)**2) / &
         real(max(size(target_2body), 1), real32)
    loss_value = loss_value + this%component_weight(2) * &
         sum((current_3body - target_3body)**2) / &
         real(max(size(target_3body), 1), real32)
    loss_value = loss_value + this%component_weight(3) * &
         sum((current_4body - target_4body)**2) / &
         real(max(size(target_4body), 1), real32)

    do atom_i = 1, basis%natom
       do coord = 1, 3
          grad_positions(atom_i, coord) = 2._real32 * this%component_weight(1) * &
               sum((current_2body - target_2body) * grad_2body(:, atom_i, coord)) / &
               real(max(size(target_2body), 1), real32)
          grad_positions(atom_i, coord) = grad_positions(atom_i, coord) + &
               2._real32 * this%component_weight(2) * &
               sum((current_3body - target_3body) * grad_3body(:, atom_i, coord)) / &
               real(max(size(target_3body), 1), real32)
          grad_positions(atom_i, coord) = grad_positions(atom_i, coord) + &
               2._real32 * this%component_weight(3) * &
               sum((current_4body - target_4body) * grad_4body(:, atom_i, coord)) / &
               real(max(size(target_4body), 1), real32)
       end do
    end do

    deallocate(grad_2body, grad_3body, grad_4body)

  end subroutine evaluate_predictive_loss_and_gradients
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
    integer :: num_hidden
    integer :: seed_
    integer :: n_rbf_, kernel_hidden_
    integer, dimension(:), allocatable :: h_sizes
    real(real32) :: lr, lr_decay_rate_


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
    this%num_pair_vertex_features = 7 + 2 * this%num_species
    this%num_triplet_vertex_features = 6 + 3 * this%num_species

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

    this%fingerprint_dim_2body = this%nbins(1) * this%num_pairs
    this%fingerprint_dim_3body = this%nbins(2) * this%num_species
    this%fingerprint_dim_4body = this%nbins(3) * this%num_species
    this%fingerprint_dim = &
         this%fingerprint_dim_2body + &
         this%fingerprint_dim_3body + &
         this%fingerprint_dim_4body

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
    n_rbf_ = 20
    kernel_hidden_ = 64
    if (present(n_rbf)) n_rbf_ = n_rbf
    if (present(kernel_hidden)) kernel_hidden_ = kernel_hidden

    call setup_branch_network( &
         this, this%network, this%num_vertex_features, this%num_edge_features, &
         this%fingerprint_dim_2body, h_sizes, lr, lr_decay_rate_, &
         n_rbf_, kernel_hidden_)
    call setup_branch_network( &
         this, this%network_3body, this%num_pair_vertex_features, 1, &
         this%fingerprint_dim_3body, h_sizes, lr, lr_decay_rate_, &
         n_rbf_, kernel_hidden_)
    call setup_branch_network( &
         this, this%network_4body, this%num_triplet_vertex_features, 1, &
         this%fingerprint_dim_4body, h_sizes, lr, lr_decay_rate_, &
         n_rbf_, kernel_hidden_)

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
  subroutine setup_branch_network( &
       this, network, num_vertex_features, num_edge_features, &
       output_dim, hidden_layer_sizes, learning_rate, lr_decay_rate, &
       n_rbf, kernel_hidden)
    !! Set up one body-specific message-passing branch.
    implicit none

    class(gnn_fingerprint_type), intent(in) :: this
    type(network_type), intent(inout) :: network
    integer, intent(in) :: num_vertex_features
    integer, intent(in) :: num_edge_features
    integer, intent(in) :: output_dim
    integer, dimension(:), intent(in) :: hidden_layer_sizes
    real(real32), intent(in) :: learning_rate
    real(real32), intent(in) :: lr_decay_rate
    integer, intent(in) :: n_rbf
    integer, intent(in) :: kernel_hidden

    integer :: i
    class(clip_type), allocatable :: clip
    type(exp_lr_decay_type) :: lr_decay

    call network%reset()

    select case (this%layer_type)
    case (1)
       call network%add(raffle_msgpass_layer_type( &
            num_time_steps = this%num_time_steps, &
            num_vertex_features = [num_vertex_features], &
            num_edge_features = [num_edge_features], &
            num_outputs = this%gnn_output_dim, &
            n_rbf = n_rbf, &
            rbf_cutoff = this%bond_cutoff, &
            kernel_hidden = kernel_hidden, &
            message_activation = 'swish', &
            readout_activation = 'none', &
            kernel_initialiser = 'glorot_normal' &
       ))
    case (2)
       call network%add(schnet_msgpass_layer_type( &
            num_time_steps = this%num_time_steps, &
            num_vertex_features = [num_vertex_features], &
            num_edge_features = [num_edge_features], &
            num_outputs = this%gnn_output_dim, &
            n_rbf = n_rbf, &
            rbf_cutoff = this%bond_cutoff, &
            kernel_hidden = kernel_hidden, &
            message_activation = 'swish', &
            readout_activation = 'none', &
            kernel_initialiser = 'glorot_normal' &
       ))
    case (3)
       call network%add(dimenet_msgpass_layer_type( &
            num_time_steps = this%num_time_steps, &
            num_vertex_features = [num_vertex_features], &
            num_edge_features = [num_edge_features], &
            num_outputs = this%gnn_output_dim, &
            n_rbf = n_rbf, &
            rbf_cutoff = this%bond_cutoff, &
            kernel_hidden = kernel_hidden, &
            message_activation = 'swish', &
            readout_activation = 'none', &
            kernel_initialiser = 'glorot_normal' &
       ))
    case (4)
       call network%add(hybrid_msgpass_layer_type( &
            num_time_steps = this%num_time_steps, &
            num_vertex_features = [num_vertex_features], &
            num_edge_features = [num_edge_features], &
            num_outputs = this%gnn_output_dim, &
            n_rbf = n_rbf, &
            rbf_cutoff = this%bond_cutoff, &
            kernel_hidden = kernel_hidden, &
            message_activation = 'swish', &
            readout_activation = 'none', &
            kernel_initialiser = 'glorot_normal' &
       ))
    case default
       call network%add(duvenaud_msgpass_layer_type( &
            num_time_steps = this%num_time_steps, &
            num_vertex_features = [num_vertex_features], &
            num_edge_features = [num_edge_features], &
            num_outputs = this%gnn_output_dim, &
            kernel_initialiser = 'glorot_normal', &
            readout_activation = 'none', &
            min_vertex_degree = 0, &
            max_vertex_degree = this%max_degree &
       ))
    end select

    do i = 1, size(hidden_layer_sizes)
       if (i == 1) then
          call network%add(full_layer_type( &
               num_inputs = this%gnn_output_dim, &
               num_outputs = hidden_layer_sizes(i), &
               activation = 'leaky_relu', &
               kernel_initialiser = 'he_normal' &
          ))
       else
          call network%add(full_layer_type( &
               num_outputs = hidden_layer_sizes(i), &
               activation = 'leaky_relu', &
               kernel_initialiser = 'he_normal' &
          ))
       end if
    end do

    call network%add(full_layer_type( &
         num_outputs = output_dim, &
         activation = 'none', &
         kernel_initialiser = 'glorot_normal' &
    ))

    allocate(clip, source=clip_type( &
         clip_min = -1.E-1_real32, &
         clip_max = 1.E-1_real32, &
         clip_norm = 1.E-1_real32 &
    ))
    lr_decay = exp_lr_decay_type(lr_decay_rate)
    lr_decay%iterate_per_epoch = .true.
    call network%compile( &
         optimiser = adam_optimiser_type( &
              learning_rate = learning_rate, &
              clip_dict = clip, &
              lr_decay = lr_decay &
         ), &
         loss_method = 'mse', &
         metrics = ['loss'], &
         verbose = 0 &
    )

  end subroutine setup_branch_network
!###############################################################################


!###############################################################################
  subroutine basis_to_graph(this, basis, graph)
    !! Convert a basis_type atomic structure to a graph_type.
    implicit none

    class(gnn_fingerprint_type), intent(in) :: this
    type(basis_type), intent(in) :: basis
    type(graph_type), intent(out) :: graph

    integer :: atom_i, pair_i
    real(real32), dimension(:,:), allocatable :: positions
    integer, dimension(:), allocatable :: species_index
    real(real32), dimension(:), allocatable :: atomic_numbers, cov_radii
    integer, dimension(:,:), allocatable :: pair_atoms
    real(real32), dimension(:,:), allocatable :: pair_shift
    real(real32), dimension(:), allocatable :: pair_distance
    real(real32) :: coord_scale
    call extract_structure_data( &
         this, basis, positions, species_index, atomic_numbers, cov_radii)
    call build_pair_data(this, basis, positions, pair_atoms, pair_shift, pair_distance)

    coord_scale = max(this%bond_cutoff, 1.E-6_real32)
    call graph%set_num_vertices(basis%natom, &
         num_vertex_features = this%num_vertex_features)

    do atom_i = 1, basis%natom
       allocate(graph%vertex(atom_i)%feature(this%num_vertex_features))
       graph%vertex(atom_i)%feature = 0._real32
       graph%vertex(atom_i)%feature(1:3) = positions(:, atom_i) / coord_scale
       if (species_index(atom_i) > 0) then
          graph%vertex(atom_i)%feature(3 + species_index(atom_i)) = 1._real32
       end if
       graph%vertex(atom_i)%feature(3 + this%num_species + 1) = atomic_numbers(atom_i)
       graph%vertex(atom_i)%feature(3 + this%num_species + 2) = cov_radii(atom_i)
    end do

    do pair_i = 1, size(pair_distance)
       call graph%add_edge(index = pair_atoms(:, pair_i), &
            feature = [pair_distance(pair_i)])
    end do

    call graph%add_self_loops(features = [0._real32])
    if (.not. graph%is_sparse) call graph%convert_to_sparse()

  end subroutine basis_to_graph
!###############################################################################


!###############################################################################
  subroutine basis_to_graph_3body(this, basis, graph)
    !! Convert a structure to the fixed 3-body graph (pair vertices, angle edges).
    implicit none

    class(gnn_fingerprint_type), intent(in) :: this
    type(basis_type), intent(in) :: basis
    type(graph_type), intent(out) :: graph

    integer :: pair_i, pair_j, num_vertices, shared_atom, other_i, other_j
    real(real32), dimension(:,:), allocatable :: positions
    integer, dimension(:), allocatable :: species_index
    real(real32), dimension(:), allocatable :: atomic_numbers, cov_radii
    integer, dimension(:,:), allocatable :: pair_atoms
    real(real32), dimension(:,:), allocatable :: pair_shift
    real(real32), dimension(:), allocatable :: pair_distance
    real(real32), dimension(:), allocatable :: pair_features
    real(real32), dimension(3) :: vec_i, vec_j
    real(real32) :: angle_ijk

    call extract_structure_data( &
         this, basis, positions, species_index, atomic_numbers, cov_radii)
    call build_pair_data(this, basis, positions, pair_atoms, pair_shift, pair_distance)

    num_vertices = max(1, size(pair_distance))
    call graph%set_num_vertices(num_vertices, &
         num_vertex_features = this%num_pair_vertex_features)

    allocate(pair_features(this%num_pair_vertex_features))
    do pair_i = 1, num_vertices
       allocate(graph%vertex(pair_i)%feature(this%num_pair_vertex_features))
       graph%vertex(pair_i)%feature = 0._real32
       if (pair_i <= size(pair_distance)) then
          call build_pair_vertex_features( &
               this, positions, species_index, pair_atoms, pair_shift, pair_distance, &
               pair_i, pair_features)
          graph%vertex(pair_i)%feature = pair_features
       end if
    end do
    deallocate(pair_features)

    if (size(pair_distance) >= 2) then
       do pair_i = 1, size(pair_distance) - 1
          do pair_j = pair_i + 1, size(pair_distance)
             call find_shared_pair_atom( &
                  pair_atoms(:, pair_i), pair_atoms(:, pair_j), &
                  shared_atom, other_i, other_j)
             if (shared_atom == 0) cycle
             vec_i = minimum_image_delta(basis, positions(:, other_i) - &
                  positions(:, shared_atom))
             vec_j = minimum_image_delta(basis, positions(:, other_j) - &
                  positions(:, shared_atom))
             if (sqrt(sum(vec_i**2)) <= 1.E-6_real32) cycle
             if (sqrt(sum(vec_j**2)) <= 1.E-6_real32) cycle
             angle_ijk = get_angle(vec_i, vec_j)
             if (.not. ieee_is_finite(angle_ijk)) cycle
             call graph%add_edge(index = [pair_i, pair_j], feature = [angle_ijk])
          end do
       end do
    end if

    call graph%add_self_loops(features = [0._real32])
    if (.not. graph%is_sparse) call graph%convert_to_sparse()

  end subroutine basis_to_graph_3body
!###############################################################################


!###############################################################################
  subroutine basis_to_graph_4body(this, basis, graph)
    !! Convert a structure to the fixed 4-body graph (triplet vertices, dihedral edges).
    implicit none

    class(gnn_fingerprint_type), intent(in) :: this
    type(basis_type), intent(in) :: basis
    type(graph_type), intent(out) :: graph

    integer :: atom_j, pair_i, pair_k, num_triplets, triplet_i, triplet_j
    integer :: num_vertices, atom_i, atom_k, atom_l
    real(real32), dimension(:,:), allocatable :: positions
    integer, dimension(:), allocatable :: species_index
    real(real32), dimension(:), allocatable :: atomic_numbers, cov_radii
    integer, dimension(:,:), allocatable :: pair_atoms
    real(real32), dimension(:,:), allocatable :: pair_shift
    real(real32), dimension(:), allocatable :: pair_distance
    integer, dimension(:,:), allocatable :: triplet_atoms
    real(real32), dimension(:), allocatable :: triplet_features
    real(real32) :: dihedral_angle

    call extract_structure_data( &
         this, basis, positions, species_index, atomic_numbers, cov_radii)
    call build_pair_data(this, basis, positions, pair_atoms, pair_shift, pair_distance)

    num_triplets = 0
    do atom_j = 1, basis%natom
       do pair_i = 1, size(pair_distance)
          if (all(pair_atoms(:, pair_i) /= atom_j)) cycle
          atom_i = pair_atoms(1, pair_i)
          if (atom_i == atom_j) atom_i = pair_atoms(2, pair_i)
          do pair_k = 1, size(pair_distance)
             if (pair_k == pair_i) cycle
             if (all(pair_atoms(:, pair_k) /= atom_j)) cycle
             atom_k = pair_atoms(1, pair_k)
             if (atom_k == atom_j) atom_k = pair_atoms(2, pair_k)
             if (atom_i == atom_k) cycle
             num_triplets = num_triplets + 1
          end do
       end do
    end do

    num_vertices = max(1, num_triplets)
    call graph%set_num_vertices(num_vertices, &
         num_vertex_features = this%num_triplet_vertex_features)

    allocate(triplet_atoms(3, num_vertices), source = 0)
    allocate(triplet_features(this%num_triplet_vertex_features))

    num_triplets = 0
    do atom_j = 1, basis%natom
       do pair_i = 1, size(pair_distance)
          if (all(pair_atoms(:, pair_i) /= atom_j)) cycle
          atom_i = pair_atoms(1, pair_i)
          if (atom_i == atom_j) atom_i = pair_atoms(2, pair_i)
          do pair_k = 1, size(pair_distance)
             if (pair_k == pair_i) cycle
             if (all(pair_atoms(:, pair_k) /= atom_j)) cycle
             atom_k = pair_atoms(1, pair_k)
             if (atom_k == atom_j) atom_k = pair_atoms(2, pair_k)
             if (atom_i == atom_k) cycle
             num_triplets = num_triplets + 1
             triplet_atoms(:, num_triplets) = [atom_i, atom_j, atom_k]
          end do
       end do
    end do

    do triplet_i = 1, num_vertices
       allocate(graph%vertex(triplet_i)%feature(this%num_triplet_vertex_features))
       graph%vertex(triplet_i)%feature = 0._real32
       if (triplet_i <= num_triplets) then
          call build_triplet_vertex_features( &
               this, basis, positions, species_index, triplet_atoms, triplet_i, &
               triplet_features)
          graph%vertex(triplet_i)%feature = triplet_features
       end if
    end do
    deallocate(triplet_features)

    if (num_triplets >= 2) then
       do triplet_i = 1, num_triplets
          do triplet_j = 1, num_triplets
             if (triplet_i == triplet_j) cycle
             if (triplet_atoms(2, triplet_i) /= triplet_atoms(1, triplet_j)) cycle
             if (triplet_atoms(3, triplet_i) /= triplet_atoms(2, triplet_j)) cycle
             atom_i = triplet_atoms(1, triplet_i)
             atom_j = triplet_atoms(2, triplet_i)
             atom_k = triplet_atoms(3, triplet_i)
             atom_l = triplet_atoms(3, triplet_j)
             if (atom_i == atom_l) cycle
             dihedral_angle = get_improper_dihedral_angle( &
                  positions(:, atom_i), positions(:, atom_j), &
                  positions(:, atom_k), positions(:, atom_l))
             if (.not. ieee_is_finite(dihedral_angle)) cycle
             call graph%add_edge(index = [triplet_i, triplet_j], &
                  feature = [dihedral_angle])
          end do
       end do
    end if

    call graph%add_self_loops(features = [0._real32])
    if (.not. graph%is_sparse) call graph%convert_to_sparse()

  end subroutine basis_to_graph_4body
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
  subroutine compute_fingerprint_components(this, basis, fingerprint_2body, &
       fingerprint_3body, fingerprint_4body)
    !! Compute the analytical RAFFLE fingerprint split into 2/3/4-body blocks.
    implicit none

    class(gnn_fingerprint_type), intent(in) :: this
    type(basis_type), intent(in) :: basis
    real(real32), dimension(:), intent(out) :: fingerprint_2body
    real(real32), dimension(:), intent(out) :: fingerprint_3body
    real(real32), dimension(:), intent(out) :: fingerprint_4body

    type(distribs_type) :: distribs
    integer :: i, j, offset

    call distribs%calculate(basis, &
         nbins = this%nbins, &
         width = this%width, &
         sigma = this%sigma, &
         cutoff_min = this%cutoff_min, &
         cutoff_max = this%cutoff_max, &
         radius_distance_tol = this%radius_distance_tol)

    fingerprint_2body = 0._real32
    fingerprint_3body = 0._real32
    fingerprint_4body = 0._real32

    offset = 0
    do j = 1, size(distribs%df_2body, 2)
       do i = 1, size(distribs%df_2body, 1)
          offset = offset + 1
          fingerprint_2body(offset) = distribs%df_2body(i, j)
       end do
    end do

    offset = 0
    do j = 1, size(distribs%df_3body, 2)
       do i = 1, size(distribs%df_3body, 1)
          offset = offset + 1
          fingerprint_3body(offset) = distribs%df_3body(i, j)
       end do
    end do

    offset = 0
    do j = 1, size(distribs%df_4body, 2)
       do i = 1, size(distribs%df_4body, 1)
          offset = offset + 1
          fingerprint_4body(offset) = distribs%df_4body(i, j)
       end do
    end do

  end subroutine compute_fingerprint_components
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
  subroutine predict_components_internal(this, basis, fingerprint_2body, &
       fingerprint_3body, fingerprint_4body)
    !! Evaluate the three trained branch networks and return denormalised outputs.
    implicit none

    class(gnn_fingerprint_type), intent(inout) :: this
    type(basis_type), intent(in) :: basis
    real(real32), dimension(:), intent(out) :: fingerprint_2body
    real(real32), dimension(:), intent(out) :: fingerprint_3body
    real(real32), dimension(:), intent(out) :: fingerprint_4body

    type(graph_type), dimension(1,1) :: graph_2body, graph_3body, graph_4body
    integer :: leaf_id, offset

    call this%basis_to_graph(basis, graph_2body(1,1))
    call this%basis_to_graph_3body(basis, graph_3body(1,1))
    call this%basis_to_graph_4body(basis, graph_4body(1,1))

    call this%network%set_batch_size(1)
    call this%network%set_inference_mode()
    call this%network%forward(graph_2body)
    leaf_id = this%network%leaf_vertices(1)
    fingerprint_2body = this%network%model(leaf_id)%layer%output(1,1)%val(:, 1)

    call this%network_3body%set_batch_size(1)
    call this%network_3body%set_inference_mode()
    call this%network_3body%forward(graph_3body)
    leaf_id = this%network_3body%leaf_vertices(1)
    fingerprint_3body = this%network_3body%model(leaf_id)%layer%output(1,1)%val(:, 1)

    call this%network_4body%set_batch_size(1)
    call this%network_4body%set_inference_mode()
    call this%network_4body%forward(graph_4body)
    leaf_id = this%network_4body%leaf_vertices(1)
    fingerprint_4body = this%network_4body%model(leaf_id)%layer%output(1,1)%val(:, 1)

    if (allocated(this%target_comp_weights)) then
       where(this%target_comp_weights(1:this%fingerprint_dim_2body) > &
            1.E-12_real32)
          fingerprint_2body = fingerprint_2body / &
               this%target_comp_weights(1:this%fingerprint_dim_2body)
       end where
       offset = this%fingerprint_dim_2body
       where(this%target_comp_weights(offset+1:offset+this%fingerprint_dim_3body) > &
            1.E-12_real32)
          fingerprint_3body = fingerprint_3body / &
               this%target_comp_weights(offset+1:offset+this%fingerprint_dim_3body)
       end where
       offset = offset + this%fingerprint_dim_3body
       where(this%target_comp_weights(offset+1:offset+this%fingerprint_dim_4body) > &
            1.E-12_real32)
          fingerprint_4body = fingerprint_4body / &
               this%target_comp_weights(offset+1:offset+this%fingerprint_dim_4body)
       end where
    end if

    if (allocated(this%target_scale) .and. allocated(this%target_mean)) then
       fingerprint_2body = fingerprint_2body * &
            this%target_scale(1:this%fingerprint_dim_2body) + &
            this%target_mean(1:this%fingerprint_dim_2body)
       offset = this%fingerprint_dim_2body
       fingerprint_3body = fingerprint_3body * &
            this%target_scale(offset+1:offset+this%fingerprint_dim_3body) + &
            this%target_mean(offset+1:offset+this%fingerprint_dim_3body)
       offset = offset + this%fingerprint_dim_3body
       fingerprint_4body = fingerprint_4body * &
            this%target_scale(offset+1:offset+this%fingerprint_dim_4body) + &
            this%target_mean(offset+1:offset+this%fingerprint_dim_4body)
    end if

  end subroutine predict_components_internal
!###############################################################################


!###############################################################################
  subroutine evaluate_components(this, basis, use_predict, fingerprint_2body, &
       fingerprint_3body, fingerprint_4body)
    !! Evaluate either analytical or learned component fingerprints.
    implicit none

    class(gnn_fingerprint_type), intent(inout) :: this
    type(basis_type), intent(in) :: basis
    logical, intent(in) :: use_predict
    real(real32), dimension(:), intent(out) :: fingerprint_2body
    real(real32), dimension(:), intent(out) :: fingerprint_3body
    real(real32), dimension(:), intent(out) :: fingerprint_4body

    if (use_predict) then
       call predict_components_internal( &
            this, basis, fingerprint_2body, fingerprint_3body, fingerprint_4body)
    else
       call this%compute_fingerprint_components( &
            basis, fingerprint_2body, fingerprint_3body, fingerprint_4body)
    end if

  end subroutine evaluate_components
!###############################################################################


!###############################################################################
  subroutine perturb_basis_coordinate(basis, atom_index, coord, delta)
    !! Apply a Cartesian perturbation to a flat atom index inside a basis.
    implicit none

    type(basis_type), intent(inout) :: basis
    integer, intent(in) :: atom_index
    integer, intent(in) :: coord
    real(real32), intent(in) :: delta

    integer :: is, ia, atom_i

    atom_i = 0
    do is = 1, basis%nspec
       do ia = 1, basis%spec(is)%num
          atom_i = atom_i + 1
          if (atom_i /= atom_index) cycle
          basis%spec(is)%atom(ia, coord) = basis%spec(is)%atom(ia, coord) + delta
          return
       end do
    end do

  end subroutine perturb_basis_coordinate
!###############################################################################


!###############################################################################
  subroutine compute_component_gradients_numerical(this, basis, use_predict, &
       grad_2body, grad_3body, grad_4body)
    !! Numerical component Jacobians used by the public gradient API and inverse design.
    implicit none

    class(gnn_fingerprint_type), intent(inout) :: this
    type(basis_type), intent(in) :: basis
    logical, intent(in) :: use_predict
    real(real32), dimension(:,:,:), intent(out) :: grad_2body
    real(real32), dimension(:,:,:), intent(out) :: grad_3body
    real(real32), dimension(:,:,:), intent(out) :: grad_4body

    type(basis_type) :: basis_fwd, basis_bwd
    real(real32), dimension(:), allocatable :: fp2_fwd, fp2_bwd
    real(real32), dimension(:), allocatable :: fp3_fwd, fp3_bwd
    real(real32), dimension(:), allocatable :: fp4_fwd, fp4_bwd
    real(real32) :: delta
    integer :: atom_i, coord

    delta = 1.E-3_real32
    allocate(fp2_fwd(this%fingerprint_dim_2body), fp2_bwd(this%fingerprint_dim_2body))
    allocate(fp3_fwd(this%fingerprint_dim_3body), fp3_bwd(this%fingerprint_dim_3body))
    allocate(fp4_fwd(this%fingerprint_dim_4body), fp4_bwd(this%fingerprint_dim_4body))

    grad_2body = 0._real32
    grad_3body = 0._real32
    grad_4body = 0._real32

    do atom_i = 1, basis%natom
       do coord = 1, 3
          basis_fwd = basis
          basis_bwd = basis
          if (.not. basis_fwd%lcart) call basis_fwd%convert()
          if (.not. basis_bwd%lcart) call basis_bwd%convert()
          call perturb_basis_coordinate(basis_fwd, atom_i, coord, delta)
          call perturb_basis_coordinate(basis_bwd, atom_i, coord, -delta)
          call evaluate_components( &
               this, basis_fwd, use_predict, fp2_fwd, fp3_fwd, fp4_fwd)
          call evaluate_components( &
               this, basis_bwd, use_predict, fp2_bwd, fp3_bwd, fp4_bwd)
          grad_2body(:, atom_i, coord) = (fp2_fwd - fp2_bwd) / (2._real32 * delta)
          grad_3body(:, atom_i, coord) = (fp3_fwd - fp3_bwd) / (2._real32 * delta)
          grad_4body(:, atom_i, coord) = (fp4_fwd - fp4_bwd) / (2._real32 * delta)
       end do
    end do

    deallocate(fp2_fwd, fp2_bwd, fp3_fwd, fp3_bwd, fp4_fwd, fp4_bwd)

  end subroutine compute_component_gradients_numerical
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
    type(graph_type), dimension(:,:), allocatable :: graphs_2body
    type(graph_type), dimension(:,:), allocatable :: graphs_3body
    type(graph_type), dimension(:,:), allocatable :: graphs_4body
    type(array_type), dimension(1,1) :: output_array_2body
    type(array_type), dimension(1,1) :: output_array_3body
    type(array_type), dimension(1,1) :: output_array_4body
    real(real32), dimension(:,:), allocatable :: target_data_2body
    real(real32), dimension(:,:), allocatable :: target_data_3body
    real(real32), dimension(:,:), allocatable :: target_data_4body
    real(real32), dimension(:,:), allocatable :: centred_target_data
    real(real32), dimension(:), allocatable :: fp_2body, fp_3body, fp_4body
    real(real32), dimension(:), allocatable :: comp_weights
    integer :: offset, i

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

    allocate(graphs_2body(1, num_strucs))
    allocate(graphs_3body(1, num_strucs))
    allocate(graphs_4body(1, num_strucs))
    allocate(target_data_2body(this%fingerprint_dim_2body, num_strucs))
    allocate(target_data_3body(this%fingerprint_dim_3body, num_strucs))
    allocate(target_data_4body(this%fingerprint_dim_4body, num_strucs))
    allocate(fp_2body(this%fingerprint_dim_2body))
    allocate(fp_3body(this%fingerprint_dim_3body))
    allocate(fp_4body(this%fingerprint_dim_4body))

    do n = 1, num_strucs
       call this%basis_to_graph(structures(n), graphs_2body(1, n))
       call this%basis_to_graph_3body(structures(n), graphs_3body(1, n))
       call this%basis_to_graph_4body(structures(n), graphs_4body(1, n))
       call this%compute_fingerprint_components( &
            structures(n), fp_2body, fp_3body, fp_4body)
       if (.not. all(ieee_is_finite(fp_2body)) .or. &
            .not. all(ieee_is_finite(fp_3body)) .or. &
            .not. all(ieee_is_finite(fp_4body))) then
          call stop_program("gnn_fingerprint: non-finite target fingerprint")
          return
       end if
       target_data_2body(:, n) = fp_2body
       target_data_3body(:, n) = fp_3body
       target_data_4body(:, n) = fp_4body
    end do

    if (allocated(this%target_mean)) deallocate(this%target_mean)
    if (allocated(this%target_scale)) deallocate(this%target_scale)
    allocate(this%target_mean(this%fingerprint_dim))
    allocate(this%target_scale(this%fingerprint_dim))
    this%target_mean(1:this%fingerprint_dim_2body) = &
         sum(target_data_2body, dim = 2) / real(num_strucs, real32)
    offset = this%fingerprint_dim_2body
    this%target_mean(offset+1:offset+this%fingerprint_dim_3body) = &
         sum(target_data_3body, dim = 2) / real(num_strucs, real32)
    offset = offset + this%fingerprint_dim_3body
    this%target_mean(offset+1:offset+this%fingerprint_dim_4body) = &
         sum(target_data_4body, dim = 2) / real(num_strucs, real32)

    allocate(centred_target_data(this%fingerprint_dim_2body, num_strucs))
    centred_target_data = target_data_2body - spread( &
         this%target_mean(1:this%fingerprint_dim_2body), dim = 2, ncopies = num_strucs)
    this%target_scale(1:this%fingerprint_dim_2body) = sqrt( &
         sum(centred_target_data**2, dim = 2) / real(max(num_strucs, 1), real32))
    where (this%target_scale(1:this%fingerprint_dim_2body) < 1.E-6_real32)
       this%target_scale(1:this%fingerprint_dim_2body) = 1._real32
    end where
    target_data_2body = centred_target_data / spread( &
         this%target_scale(1:this%fingerprint_dim_2body), dim = 2, ncopies = num_strucs)
    deallocate(centred_target_data)

    allocate(centred_target_data(this%fingerprint_dim_3body, num_strucs))
    offset = this%fingerprint_dim_2body
    centred_target_data = target_data_3body - spread( &
         this%target_mean(offset+1:offset+this%fingerprint_dim_3body), &
         dim = 2, ncopies = num_strucs)
    this%target_scale(offset+1:offset+this%fingerprint_dim_3body) = sqrt( &
         sum(centred_target_data**2, dim = 2) / real(max(num_strucs, 1), real32))
    where (this%target_scale(offset+1:offset+this%fingerprint_dim_3body) < 1.E-6_real32)
       this%target_scale(offset+1:offset+this%fingerprint_dim_3body) = 1._real32
    end where
    target_data_3body = centred_target_data / spread( &
         this%target_scale(offset+1:offset+this%fingerprint_dim_3body), &
         dim = 2, ncopies = num_strucs)
    deallocate(centred_target_data)

    allocate(centred_target_data(this%fingerprint_dim_4body, num_strucs))
    offset = this%fingerprint_dim_2body + this%fingerprint_dim_3body
    centred_target_data = target_data_4body - spread( &
         this%target_mean(offset+1:offset+this%fingerprint_dim_4body), &
         dim = 2, ncopies = num_strucs)
    this%target_scale(offset+1:offset+this%fingerprint_dim_4body) = sqrt( &
         sum(centred_target_data**2, dim = 2) / real(max(num_strucs, 1), real32))
    where (this%target_scale(offset+1:offset+this%fingerprint_dim_4body) < 1.E-6_real32)
       this%target_scale(offset+1:offset+this%fingerprint_dim_4body) = 1._real32
    end where
    target_data_4body = centred_target_data / spread( &
         this%target_scale(offset+1:offset+this%fingerprint_dim_4body), &
         dim = 2, ncopies = num_strucs)
    deallocate(centred_target_data)

    allocate(comp_weights(this%fingerprint_dim))
    comp_weights = 1._real32
    offset = 0
    do i = 1, this%fingerprint_dim_2body
       offset = offset + 1
       comp_weights(offset) = sqrt(this%component_weight(1))
    end do
    do i = 1, this%fingerprint_dim_3body
       offset = offset + 1
       comp_weights(offset) = sqrt(this%component_weight(2))
    end do
    do i = 1, this%fingerprint_dim_4body
       offset = offset + 1
       comp_weights(offset) = sqrt(this%component_weight(3))
    end do
    target_data_2body = target_data_2body * spread( &
         comp_weights(1:this%fingerprint_dim_2body), dim = 2, ncopies = num_strucs)
    offset = this%fingerprint_dim_2body
    target_data_3body = target_data_3body * spread( &
         comp_weights(offset+1:offset+this%fingerprint_dim_3body), &
         dim = 2, ncopies = num_strucs)
    offset = offset + this%fingerprint_dim_3body
    target_data_4body = target_data_4body * spread( &
         comp_weights(offset+1:offset+this%fingerprint_dim_4body), &
         dim = 2, ncopies = num_strucs)
    if (allocated(this%target_comp_weights)) deallocate(this%target_comp_weights)
    allocate(this%target_comp_weights, source=comp_weights)
    deallocate(comp_weights)

    if (.not. all(ieee_is_finite(target_data_2body)) .or. &
         .not. all(ieee_is_finite(target_data_3body)) .or. &
         .not. all(ieee_is_finite(target_data_4body))) then
       call stop_program("gnn_fingerprint: non-finite normalised training targets")
       return
    end if

    call output_array_2body(1,1)%allocate(array_shape = &
         [this%fingerprint_dim_2body, num_strucs])
    output_array_2body(1,1)%val = target_data_2body
    call output_array_3body(1,1)%allocate(array_shape = &
         [this%fingerprint_dim_3body, num_strucs])
    output_array_3body(1,1)%val = target_data_3body
    call output_array_4body(1,1)%allocate(array_shape = &
         [this%fingerprint_dim_4body, num_strucs])
    output_array_4body(1,1)%val = target_data_4body

    call this%network%train( &
         graphs_2body, &
         output_array_2body, &
         num_epochs = epochs, &
         batch_size = batch_size_, &
         shuffle_batches = .true., &
         verbose = verb &
    )
    call this%network_3body%train( &
         graphs_3body, &
         output_array_3body, &
         num_epochs = epochs, &
         batch_size = batch_size_, &
         shuffle_batches = .true., &
         verbose = verb &
    )
    call this%network_4body%train( &
         graphs_4body, &
         output_array_4body, &
         num_epochs = epochs, &
         batch_size = batch_size_, &
         shuffle_batches = .true., &
         verbose = verb &
    )

    this%is_trained = .true.

    deallocate(graphs_2body, graphs_3body, graphs_4body)
    deallocate(target_data_2body, target_data_3body, target_data_4body)
    deallocate(fp_2body, fp_3body, fp_4body)

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

    real(real32), dimension(:), allocatable :: fingerprint_2body
    real(real32), dimension(:), allocatable :: fingerprint_3body
    real(real32), dimension(:), allocatable :: fingerprint_4body

    if (.not. this%is_initialised) then
       call stop_program("gnn_fingerprint: network not initialised")
       return
    end if

    allocate(fingerprint_2body(this%fingerprint_dim_2body))
    allocate(fingerprint_3body(this%fingerprint_dim_3body))
    allocate(fingerprint_4body(this%fingerprint_dim_4body))

    call predict_components_internal( &
         this, basis, fingerprint_2body, fingerprint_3body, fingerprint_4body)

    fingerprint(1:this%fingerprint_dim_2body) = fingerprint_2body
    fingerprint(this%fingerprint_dim_2body+1: &
         this%fingerprint_dim_2body+this%fingerprint_dim_3body) = fingerprint_3body
    fingerprint(this%fingerprint_dim_2body+this%fingerprint_dim_3body+1: &
         this%fingerprint_dim) = fingerprint_4body

    if (.not. all(ieee_is_finite(fingerprint))) then
       call stop_program("gnn_fingerprint: non-finite prediction")
       return
    end if

    deallocate(fingerprint_2body, fingerprint_3body, fingerprint_4body)

  end subroutine predict
!###############################################################################


!###############################################################################
  subroutine compute_gradients(this, basis, grad_2body, grad_3body, grad_4body)
    !! Compute learned fingerprint gradients with respect to Cartesian positions.
    implicit none

    class(gnn_fingerprint_type), intent(inout) :: this
    type(basis_type), intent(in) :: basis
    real(real32), dimension(:,:,:), intent(out) :: grad_2body
    real(real32), dimension(:,:,:), intent(out) :: grad_3body
    real(real32), dimension(:,:,:), intent(out) :: grad_4body

    type(basis_type) :: working_basis
    type(multigraph_data_type) :: data
    type(graph_type), dimension(1,1) :: graph_2body, graph_3body, graph_4body

    if (.not. this%is_initialised) then
       call stop_program("gnn_fingerprint: network not initialised")
       return
    end if

    working_basis = basis
    if (.not. working_basis%lcart) call working_basis%convert()

    call prepare_multigraph_data(this, working_basis, data)
    call this%basis_to_graph(working_basis, graph_2body(1,1))
    call this%basis_to_graph_3body(working_basis, graph_3body(1,1))
    call this%basis_to_graph_4body(working_basis, graph_4body(1,1))

    call compute_branch_output_jacobian(this, this%network, 1, graph_2body, &
         working_basis, data, grad_2body)
    call compute_branch_output_jacobian(this, this%network_3body, 2, graph_3body, &
         working_basis, data, grad_3body)
    call compute_branch_output_jacobian(this, this%network_4body, 3, graph_4body, &
         working_basis, data, grad_4body)

  end subroutine compute_gradients
!###############################################################################


!###############################################################################
  subroutine inverse_design(this, target_fingerprint, basis, &
       fixed_atoms, num_steps, step_size, verbose, use_predict)
    !! Inverse design: optimise atomic positions to match a target
    !! fingerprint.
    !!
    !! Uses exact predictive gradients for the learned model and a numerical
    !! fallback only when the analytical RAFFLE descriptor path is requested.
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

    integer :: nsteps, verb, step, atom_i, coord, num_movable, offset
    integer :: is, ia, atom_j
    real(real32) :: lr, best_loss, loss_current
    real(real32) :: overlap_scale, overlap_penalty, dist, min_dist, penalty_grad
    logical :: do_predict
    type(sgd_optimiser_type) :: opt
    type(basis_type) :: working_basis, best_basis
    real(real32), dimension(:), allocatable :: &
         current_2body, current_3body, current_4body
    real(real32), dimension(:), allocatable :: &
         target_2body, target_3body, target_4body
    real(real32), dimension(:), allocatable :: x_flat, grad_flat
    real(real32), dimension(:,:,:), allocatable :: grad_2body, grad_3body, grad_4body
    real(real32), dimension(:,:), allocatable :: grad_positions
    real(real32), dimension(:,:), allocatable :: positions
    integer, dimension(:), allocatable :: species_index
    real(real32), dimension(:), allocatable :: atomic_numbers, cov_radii
    real(real32), dimension(3) :: delta_vec

    if (.not. this%is_initialised) then
       call stop_program("gnn_fingerprint: network not initialised")
       return
    end if
    if (size(fixed_atoms) /= basis%natom) then
       call stop_program("gnn_fingerprint: fixed_atoms mask has wrong size")
       return
    end if

    nsteps = 500
    lr = 1.0_real32
    verb = 0
    do_predict = .false.
    if (present(num_steps)) nsteps = num_steps
    if (present(step_size)) lr = step_size
    if (present(verbose)) verb = verbose
    if (present(use_predict)) do_predict = use_predict

    working_basis = basis
    if (.not. working_basis%lcart) call working_basis%convert()
    best_basis = working_basis
    best_loss = huge(1._real32)
    overlap_scale = 10._real32
    do is = 1, working_basis%nspec
       do ia = 1, working_basis%spec(is)%num
          atom_i = atom_i + 1
          write(*,*) "start", working_basis%spec(is)%atom(ia, 1:3)
       end do
    end do

    allocate(target_2body(this%fingerprint_dim_2body))
    allocate(target_3body(this%fingerprint_dim_3body))
    allocate(target_4body(this%fingerprint_dim_4body))
    target_2body = target_fingerprint(1:this%fingerprint_dim_2body)
    offset = this%fingerprint_dim_2body
    target_3body = target_fingerprint(offset+1:offset+this%fingerprint_dim_3body)
    offset = offset + this%fingerprint_dim_3body
    target_4body = target_fingerprint(offset+1:offset+this%fingerprint_dim_4body)

    allocate(current_2body(this%fingerprint_dim_2body))
    allocate(current_3body(this%fingerprint_dim_3body))
    allocate(current_4body(this%fingerprint_dim_4body))
    allocate(grad_2body(this%fingerprint_dim_2body, basis%natom, 3))
    allocate(grad_3body(this%fingerprint_dim_3body, basis%natom, 3))
    allocate(grad_4body(this%fingerprint_dim_4body, basis%natom, 3))
    allocate(grad_positions(basis%natom, 3))

    num_movable = count(.not. fixed_atoms)
    allocate(x_flat(3 * max(num_movable, 1)), source = 0._real32)
    allocate(grad_flat(3 * max(num_movable, 1)), source = 0._real32)

    opt = sgd_optimiser_type( &
         learning_rate = lr, &
         lr_decay = exp_lr_decay_type(1.E-2_real32), &
         clip_dict = clip_type( &
              clip_min = -1.E-1_real32, &
              clip_max = 1.E-1_real32, &
              clip_norm = 1.E-1_real32 &
         ) )
    call opt%init(num_params = size(x_flat))

    do step = 1, nsteps
       if (do_predict) then
          call evaluate_predictive_loss_and_gradients( &
               this, working_basis, target_2body, &
               target_3body, target_4body, &
               current_2body, current_3body, current_4body, &
               grad_positions, loss_current)
       else
          call evaluate_components(this, working_basis, .false., current_2body, &
               current_3body, current_4body)
          call compute_component_gradients_numerical(this, working_basis, .false., &
               grad_2body, grad_3body, grad_4body)

          grad_positions = 0._real32
          loss_current = this%component_weight(1) * &
               sum((current_2body - target_2body)**2) / &
               real(max(this%fingerprint_dim_2body, 1), real32)
          loss_current = loss_current + this%component_weight(2) * &
               sum((current_3body - target_3body)**2) / &
               real(max(this%fingerprint_dim_3body, 1), real32)
          loss_current = loss_current + this%component_weight(3) * &
               sum((current_4body - target_4body)**2) / &
               real(max(this%fingerprint_dim_4body, 1), real32)

          do atom_i = 1, basis%natom
             do coord = 1, 3
                grad_positions(atom_i, coord) = &
                     2._real32 * this%component_weight(1) * &
                     sum((current_2body - target_2body) * &
                          grad_2body(:, atom_i, coord)) / &
                     real(max(this%fingerprint_dim_2body, 1), real32)
                grad_positions(atom_i, coord) = grad_positions(atom_i, coord) + &
                     2._real32 * this%component_weight(2) * &
                     sum((current_3body - target_3body) * &
                          grad_3body(:, atom_i, coord)) / &
                     real(max(this%fingerprint_dim_3body, 1), real32)
                grad_positions(atom_i, coord) = grad_positions(atom_i, coord) + &
                     2._real32 * this%component_weight(3) * &
                     sum((current_4body - target_4body) * &
                          grad_4body(:, atom_i, coord)) / &
                     real(max(this%fingerprint_dim_4body, 1), real32)
             end do
          end do
       end if
       call extract_structure_data( &
            this, working_basis, positions, species_index, atomic_numbers, cov_radii)
       overlap_penalty = 0._real32
       do atom_i = 1, basis%natom - 1
          do atom_j = atom_i + 1, basis%natom
             delta_vec = minimum_image_delta( &
                  working_basis, positions(:, atom_j) - positions(:, atom_i))
             dist = sqrt(sum(delta_vec**2))
             min_dist = 0.75_real32 * (cov_radii(atom_i) + cov_radii(atom_j))
             if (dist >= min_dist .or. min_dist <= 1.E-6_real32) cycle
             overlap_penalty = (min_dist - dist) / min_dist
             loss_current = loss_current + overlap_scale * overlap_penalty**2
             if (dist > 1.E-6_real32) then
                penalty_grad = &
                     2._real32 * overlap_scale * overlap_penalty / (min_dist * dist)
                grad_positions(atom_i, :) = &
                     grad_positions(atom_i, :) + penalty_grad * delta_vec
                grad_positions(atom_j, :) = &
                     grad_positions(atom_j, :) - penalty_grad * delta_vec
             end if
          end do
       end do
       deallocate(positions, species_index, atomic_numbers, cov_radii)

       if (loss_current < best_loss) then
          best_loss = loss_current
          best_basis = working_basis
       end if

       if (verb > 0) then
          write(*,'(A,I5,A,ES12.4)') "  step=", step, " loss=", loss_current
       end if

       if (num_movable == 0) cycle
       offset = 0
       atom_i = 0
       do is = 1, working_basis%nspec
          do ia = 1, working_basis%spec(is)%num
             atom_i = atom_i + 1
             if (fixed_atoms(atom_i)) cycle
             x_flat(offset+1:offset+3) = working_basis%spec(is)%atom(ia, 1:3)
             grad_flat(offset+1:offset+3) = grad_positions(atom_i, 1:3)
             offset = offset + 3
          end do
       end do

       call opt%clip_dict%apply(offset, grad_flat(1:offset))
       call opt%minimise(param = x_flat(1:offset), gradient = grad_flat(1:offset))

       offset = 0
       atom_i = 0
       do is = 1, working_basis%nspec
          do ia = 1, working_basis%spec(is)%num
             atom_i = atom_i + 1
             if (fixed_atoms(atom_i)) cycle
             working_basis%spec(is)%atom(ia, 1:3) = x_flat(offset+1:offset+3)
             offset = offset + 3
          end do
       end do
    end do

    basis = best_basis
    basis%lcart = .true.
    open(unit=100, file="POSCAR_optimized_structure", status="replace")
    call geom_write(100, basis)
    close(100)

  end subroutine inverse_design
!###############################################################################


end module raffle__gnn_fingerprint_multihead
