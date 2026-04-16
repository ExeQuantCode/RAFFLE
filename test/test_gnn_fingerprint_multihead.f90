program test_gnn_fingerprint
  use raffle__constants, only: real32
  use raffle__io_utils
  use raffle__geom_rw, only: basis_type, geom_write
  use raffle__gnn_fingerprint_multihead, only: gnn_fingerprint_type
  use athena, only: graph_type
  implicit none

  type(gnn_fingerprint_type) :: gnn, gnn_mlip
  type(basis_type) :: basis1, basis2, basis2_before
  type(basis_type), dimension(2) :: training_set
  type(graph_type) :: graph1, graph3, graph4
  real(real32), dimension(:), allocatable :: fingerprint1, fingerprint2
  real(real32), dimension(:), allocatable :: fingerprint1_2body
  real(real32), dimension(:), allocatable :: fingerprint1_3body
  real(real32), dimension(:), allocatable :: fingerprint1_4body
  real(real32), dimension(:), allocatable :: predicted_fp
  real(real32), dimension(:), allocatable :: predicted_fp_mlip
  real(real32), dimension(:), allocatable :: target_fp
  real(real32), dimension(:,:,:), allocatable :: grad_2body
  real(real32), dimension(:,:,:), allocatable :: grad_3body
  real(real32), dimension(:,:,:), allocatable :: grad_4body
  logical, dimension(:), allocatable :: fixed_mask
  real(real32) :: loss
  integer :: i, atom_i
  logical :: success

  success = .true.
  test_error_handling = .true.

  !-----------------------------------------------------------------------------
  ! Set up test structures (diamond cubic carbon, 8 atoms)
  !-----------------------------------------------------------------------------
  call setup_carbon_diamond(basis1)
  call setup_carbon_diamond_perturbed(basis2)

  !-----------------------------------------------------------------------------
  ! Test 1: Initialisation
  !-----------------------------------------------------------------------------
  write(*,*) "Test 1: Initialise GNN fingerprint module"
  call gnn%initialise( &
       species_list = [character(len=3) :: 'C  '], &
       num_time_steps = 2, &
       gnn_output_dim = 16, &
       max_degree = 8, &
       hidden_layer_sizes = [32], &
       learning_rate = 0.001_real32 &
  )
  call assert(gnn%is_initialised, &
       'GNN fingerprint failed to initialise', success)
  call assert(gnn%num_vertex_features == 6, &
       'GNN num_vertex_features should be 6 (3 coords + 1 species + Z + radius)', &
       success)
  call assert(gnn%num_edge_features == 1, &
       'GNN num_edge_features should be 1 (distance)', success)
  call assert(gnn%num_pair_vertex_features > gnn%num_vertex_features, &
       '3-body pair vertex features should exceed 2-body vertex features', success)
  call assert(gnn%num_triplet_vertex_features >= gnn%num_pair_vertex_features, &
       '4-body triplet vertex features should be at least as large as &
       &3-body pair features', success)
  call assert(gnn%fingerprint_dim > 0, &
       'GNN fingerprint_dim should be positive', success)
  call assert(gnn%num_species == 1, &
       'GNN num_species should be 1', success)
  write(*,*) "  num_vertex_features = ", gnn%num_vertex_features
  write(*,*) "  gnn_output_dim = ", gnn%gnn_output_dim
  write(*,*) "  fingerprint_dim = ", gnn%fingerprint_dim
  write(*,*) "  num_pairs = ", gnn%num_pairs

  !-----------------------------------------------------------------------------
  ! Test 2: Basis to graph conversion
  !-----------------------------------------------------------------------------
  write(*,*) "Test 2: basis_to_graph conversion"
  call gnn%basis_to_graph(basis1, graph1)
  call gnn%basis_to_graph_3body(basis1, graph3)
  call gnn%basis_to_graph_4body(basis1, graph4)
  call assert(graph1%num_vertices == 8, &
       'Graph should have 8 vertices (atoms)', success)
  call assert(graph1%num_vertex_features == 6, &
       'Graph vertex features should be 6', success)
  call assert(graph1%num_edges > 0, &
       'Graph should have edges (bonds within cutoff)', success)
  call assert(graph1%is_sparse, &
       'Graph should be sparse (CSR format)', success)
  call assert(graph3%num_vertices > 0, &
       '3-body graph should have pair vertices', success)
  call assert(graph3%num_edges > 0, &
       '3-body graph should have angle edges', success)
  call assert(graph4%num_vertices > 0, &
       '4-body graph should have triplet vertices', success)
  write(*,*) "  num_vertices = ", graph1%num_vertices
  write(*,*) "  num_edges = ", graph1%num_edges
  write(*,*) "  3-body vertices / edges = ", graph3%num_vertices, graph3%num_edges
  write(*,*) "  4-body vertices / edges = ", graph4%num_vertices, graph4%num_edges

  !-----------------------------------------------------------------------------
  ! Test 3: Compute fingerprint (RAFFLE descriptor)
  !-----------------------------------------------------------------------------
  write(*,*) "Test 3: Compute RAFFLE descriptor fingerprint"
  allocate(fingerprint1(gnn%fingerprint_dim))
  allocate(fingerprint2(gnn%fingerprint_dim))
  allocate(fingerprint1_2body(gnn%fingerprint_dim_2body))
  allocate(fingerprint1_3body(gnn%fingerprint_dim_3body))
  allocate(fingerprint1_4body(gnn%fingerprint_dim_4body))
  call gnn%compute_fingerprint(basis1, fingerprint1)
  call gnn%compute_fingerprint(basis2, fingerprint2)
  call gnn%compute_fingerprint_components( &
       basis1, fingerprint1_2body, fingerprint1_3body, fingerprint1_4body)
  call assert(any(abs(fingerprint1) > 0._real32), &
       'Fingerprint should be non-zero for valid structure', success)
  call assert(any(abs(fingerprint1_2body) > 0._real32), &
       '2-body fingerprint block should be non-zero', success)
  call assert(size(fingerprint1_2body) + size(fingerprint1_3body) + &
       size(fingerprint1_4body) == gnn%fingerprint_dim, &
       'Component fingerprints should reconstruct total fingerprint dimension', success)
  loss = sum((fingerprint1 - fingerprint2)**2)
  call assert(loss > 0._real32, &
       'Different structures should have different fingerprints', success)
  write(*,*) "  Fingerprint L2 distance between structures: ", loss

  !-----------------------------------------------------------------------------
  ! Test 4: Training with GNN
  !-----------------------------------------------------------------------------
  write(*,*) "Test 4: Train GNN on structures"
  training_set(1) = basis1
  training_set(2) = basis2
  call gnn%train(training_set, num_epochs = 5, batch_size = 1, verbose = 0)
  call assert(gnn%is_trained, 'GNN should be marked as trained', success)

  !-----------------------------------------------------------------------------
  ! Test 4b: Training with MLIP-style message passing
  !-----------------------------------------------------------------------------
  write(*,*) "Test 4b: Train GNN with MLIP-style message passing"
  call gnn_mlip%initialise( &
       species_list = [character(len=3) :: 'C  '], &
       num_time_steps = 2, &
       gnn_output_dim = 16, &
       max_degree = 8, &
       hidden_layer_sizes = [32], &
       learning_rate = 0.0001_real32, &
       use_mlip_layer = .true., &
       n_rbf = 12, &
       kernel_hidden = 32 &
  )
  call gnn_mlip%train(training_set, num_epochs = 100, batch_size = 1, verbose = 0)
  call assert(gnn_mlip%is_trained, 'MLIP GNN should be marked as trained', success)

  !-----------------------------------------------------------------------------
  ! Test 5: Forward inference (predict)
  !-----------------------------------------------------------------------------
  write(*,*) "Test 5: Forward inference via GNN"
  allocate(predicted_fp(gnn%fingerprint_dim))
  call gnn%predict(basis1, predicted_fp)
  call assert(size(predicted_fp) == gnn%fingerprint_dim, &
       'Predicted fingerprint should have correct dimension', success)
  call assert(all(predicted_fp == predicted_fp), &
       'Predicted fingerprint should not contain NaNs', success)
  write(*,*) "  Predicted fingerprint max value: ", maxval(abs(predicted_fp))

  write(*,*) "Test 5b: Forward inference via MLIP GNN"
  allocate(predicted_fp_mlip(gnn_mlip%fingerprint_dim))
  call gnn_mlip%predict(basis1, predicted_fp_mlip)
  call assert(size(predicted_fp_mlip) == gnn_mlip%fingerprint_dim, &
       'MLIP predicted fingerprint should have correct dimension', success)
  call assert(all(predicted_fp_mlip == predicted_fp_mlip), &
       'MLIP predicted fingerprint should not contain NaNs', success)

  !-----------------------------------------------------------------------------
  ! Test 5c: Gradient evaluation
  !-----------------------------------------------------------------------------
  write(*,*) "Test 5c: Gradient evaluation"
  allocate(grad_2body(gnn%fingerprint_dim_2body, basis1%natom, 3))
  allocate(grad_3body(gnn%fingerprint_dim_3body, basis1%natom, 3))
  allocate(grad_4body(gnn%fingerprint_dim_4body, basis1%natom, 3))
  call gnn%compute_gradients(basis1, grad_2body, grad_3body, grad_4body)
  call assert(all(grad_2body == grad_2body), &
       '2-body gradients should not contain NaNs', success)
  call assert(all(grad_3body == grad_3body), &
       '3-body gradients should not contain NaNs', success)
  call assert(all(grad_4body == grad_4body), &
       '4-body gradients should not contain NaNs', success)

  !-----------------------------------------------------------------------------
  ! Test 6: Inverse design with atom masking
  !-----------------------------------------------------------------------------
  write(*,*) "Test 6: Inverse design with atom mask"
  allocate(target_fp(gnn%fingerprint_dim))
  target_fp = fingerprint1

  ! Fix first 4 atoms, allow last 4 to move
  allocate(fixed_mask(basis2%natom))
  fixed_mask = .false.
  do i = 1, min(1, basis2%natom)
     fixed_mask(i) = .true.
  end do

  basis2_before = basis2
  if (.not. basis2_before%lcart) call basis2_before%convert()

  call gnn%inverse_design( &
       target_fingerprint = target_fp, &
       basis = basis2, &
       fixed_atoms = fixed_mask, &
       num_steps = 100, &
       step_size = 0.1_real32, &
       verbose = 1, &
       use_predict = .true. &
  )
  do atom_i = 1, min(4, basis2%natom)
     call assert(all(abs(basis2%spec(1)%atom(atom_i,1:3) - &
          basis2_before%spec(1)%atom(atom_i,1:3)) < 1.E-5_real32), &
     'Fixed atoms should remain unchanged during inverse design', success)
  end do
  write(*,*) "  Inverse design completed"

  open(unit=9, file='POSCAR_original', status='replace')
  call geom_write(9, basis2_before)
  close(9)

  open(unit=10, file='POSCAR_inverse', status='replace')
  call geom_write(10, basis2)
  close(10)

  !-----------------------------------------------------------------------------
  ! Test 7: Fingerprint to distribs conversion
  !-----------------------------------------------------------------------------
  write(*,*) "Test 7: Fingerprint to distribs conversion"
  block
    use raffle__distribs, only: distribs_base_type
    type(distribs_base_type) :: reconverted
    call gnn%fingerprint_to_distribs(fingerprint1, reconverted)
    call assert(allocated(reconverted%df_2body), &
         'Reconverted 2-body should be allocated', success)
    call assert(allocated(reconverted%df_3body), &
         'Reconverted 3-body should be allocated', success)
    call assert(allocated(reconverted%df_4body), &
         'Reconverted 4-body should be allocated', success)
    call assert( &
         size(reconverted%df_2body, 1) == gnn%nbins(1), &
         'Reconverted 2-body should have correct nbins', success)
  end block

  !-----------------------------------------------------------------------------
  ! Summary
  !-----------------------------------------------------------------------------
  deallocate(fingerprint1, fingerprint2, fingerprint1_2body, fingerprint1_3body)
  deallocate(fingerprint1_4body, predicted_fp, predicted_fp_mlip, target_fp)
  deallocate(fixed_mask, grad_2body, grad_3body, grad_4body)

  if (success) then
     write(*,*) "All gnn_fingerprint tests PASSED"
  else
     write(*,*) "Some gnn_fingerprint tests FAILED"
     stop 1
  end if


contains


  subroutine setup_carbon_diamond(basis)
    !! Set up a simple 8-atom carbon diamond cubic cell.
    type(basis_type), intent(out) :: basis
    real(real32) :: a

    a = 3.567_real32

    basis%sysname = "C_diamond"
    basis%nspec = 1
    basis%natom = 8
    basis%energy = -72.0_real32
    basis%lcart = .false.
    basis%pbc = [.true., .true., .true.]
    basis%lat(1,:) = [a, 0._real32, 0._real32]
    basis%lat(2,:) = [0._real32, a, 0._real32]
    basis%lat(3,:) = [0._real32, 0._real32, a]

    allocate(basis%spec(1))
    basis%spec(1)%name = 'C  '
    basis%spec(1)%num = 8
    allocate(basis%spec(1)%atom(8, 4))
    basis%spec(1)%atom(1,:) = [0.000_real32, 0.000_real32, 0.000_real32, 0._real32]
    basis%spec(1)%atom(2,:) = [0.500_real32, 0.500_real32, 0.000_real32, 0._real32]
    basis%spec(1)%atom(3,:) = [0.500_real32, 0.000_real32, 0.500_real32, 0._real32]
    basis%spec(1)%atom(4,:) = [0.000_real32, 0.500_real32, 0.500_real32, 0._real32]
    basis%spec(1)%atom(5,:) = [0.250_real32, 0.250_real32, 0.250_real32, 0._real32]
    basis%spec(1)%atom(6,:) = [0.750_real32, 0.750_real32, 0.250_real32, 0._real32]
    basis%spec(1)%atom(7,:) = [0.750_real32, 0.250_real32, 0.750_real32, 0._real32]
    basis%spec(1)%atom(8,:) = [0.250_real32, 0.750_real32, 0.750_real32, 0._real32]
  end subroutine setup_carbon_diamond


  subroutine setup_carbon_diamond_perturbed(basis)
    !! Set up a slightly perturbed 8-atom carbon diamond cell.
    type(basis_type), intent(out) :: basis

    call setup_carbon_diamond(basis)
    basis%sysname = "C_diamond_pert"
    basis%energy = -71.5_real32

    basis%spec(1)%atom(5,1) = basis%spec(1)%atom(5,1) + 0.2_real32
    basis%spec(1)%atom(6,2) = basis%spec(1)%atom(6,2) - 0.1_real32
    basis%spec(1)%atom(7,3) = basis%spec(1)%atom(7,3) + 0.15_real32
    basis%spec(1)%atom(8,1) = basis%spec(1)%atom(8,1) - 0.1_real32
  end subroutine setup_carbon_diamond_perturbed


  subroutine assert(condition, message, success)
    implicit none
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    logical, intent(inout) :: success
    if (.not. condition) then
       write(0,*) "Test failed: ", message
       success = .false.
    end if
  end subroutine assert


end program test_gnn_fingerprint
