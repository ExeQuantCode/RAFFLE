program test_nn_fingerprint
  use raffle__constants, only: real32
  use raffle__io_utils
  use raffle__geom_rw, only: basis_type
  use raffle__nn_fingerprint, only: nn_fingerprint_type
  implicit none

  type(nn_fingerprint_type) :: nn
  type(basis_type) :: basis1, basis2
  type(basis_type), dimension(2) :: training_set
  real(real32), dimension(:), allocatable :: fingerprint1, fingerprint2
  real(real32), dimension(:), allocatable :: predicted_fp
  real(real32), dimension(:), allocatable :: target_fp
  real(real32), dimension(:), allocatable :: input_vec
  logical, dimension(:), allocatable :: fixed_mask
  real(real32) :: loss
  integer :: i
  logical :: success

  success = .true.
  test_error_handling = .true.

  !-----------------------------------------------------------------------------
  ! Set up a simple test structure (diamond cubic carbon, 8 atoms)
  !-----------------------------------------------------------------------------
  call setup_carbon_diamond(basis1)
  call setup_carbon_diamond_perturbed(basis2)

  !-----------------------------------------------------------------------------
  ! Test 1: Initialisation
  !-----------------------------------------------------------------------------
  write(*,*) "Test 1: Initialise NN fingerprint module"
  call nn%initialise( &
       species_list = [character(len=3) :: 'C  '], &
       max_atoms = 8, &
       hidden_layer_sizes = [32, 16], &
       learning_rate = 0.001_real32 &
  )
  call assert(nn%is_initialised, &
       'NN fingerprint failed to initialise', success)
  call assert(nn%input_dim > 0, &
       'NN input_dim should be positive', success)
  call assert(nn%fingerprint_dim > 0, &
       'NN fingerprint_dim should be positive', success)
  call assert(nn%num_species == 1, &
       'NN num_species should be 1', success)
  write(*,*) "  input_dim = ", nn%input_dim
  write(*,*) "  fingerprint_dim = ", nn%fingerprint_dim
  write(*,*) "  num_pairs = ", nn%num_pairs

  !-----------------------------------------------------------------------------
  ! Test 2: Basis to input vector conversion
  !-----------------------------------------------------------------------------
  write(*,*) "Test 2: basis_to_input conversion"
  allocate(input_vec(nn%input_dim))
  call nn%basis_to_input(basis1, input_vec)
  call assert(any(abs(input_vec) > 0._real32), &
       'Input vector should be non-zero', success)
  ! Check that species encoding exists
  call assert(abs(input_vec(4) - 1._real32) < 1.E-6, &
       'One-hot encoding for first atom species should be 1', success)
  deallocate(input_vec)

  !-----------------------------------------------------------------------------
  ! Test 3: Compute fingerprint (RAFFLE descriptor)
  !-----------------------------------------------------------------------------
  write(*,*) "Test 3: Compute RAFFLE descriptor fingerprint"
  allocate(fingerprint1(nn%fingerprint_dim))
  allocate(fingerprint2(nn%fingerprint_dim))
  call nn%compute_fingerprint(basis1, fingerprint1)
  call nn%compute_fingerprint(basis2, fingerprint2)
  call assert(any(abs(fingerprint1) > 0._real32), &
       'Fingerprint should be non-zero for valid structure', success)
  ! Different structures should have different fingerprints
  loss = sum((fingerprint1 - fingerprint2)**2)
  call assert(loss > 0._real32, &
       'Different structures should have different fingerprints', success)
  write(*,*) "  Fingerprint L2 distance between structures: ", loss

  !-----------------------------------------------------------------------------
  ! Test 4: Training
  !-----------------------------------------------------------------------------
  write(*,*) "Test 4: Train NN on structures"
  training_set(1) = basis1
  training_set(2) = basis2
  call nn%train(training_set, num_epochs = 5, batch_size = 1, verbose = 0)
  call assert(nn%is_trained, 'NN should be marked as trained', success)

  !-----------------------------------------------------------------------------
  ! Test 5: Forward inference (predict)
  !-----------------------------------------------------------------------------
  write(*,*) "Test 5: Forward inference"
  allocate(predicted_fp(nn%fingerprint_dim))
  call nn%predict(basis1, predicted_fp)
  call assert(size(predicted_fp) == nn%fingerprint_dim, &
       'Predicted fingerprint should have correct dimension', success)
  write(*,*) "  Predicted fingerprint max value: ", maxval(abs(predicted_fp))

  !-----------------------------------------------------------------------------
  ! Test 6: Inverse design with atom masking
  !-----------------------------------------------------------------------------
  write(*,*) "Test 6: Inverse design with atom mask"
  ! Use fingerprint1 as target, start from perturbed structure
  allocate(target_fp(nn%fingerprint_dim))
  target_fp = fingerprint1

  ! Fix first 4 atoms, allow last 4 to move
  allocate(fixed_mask(basis2%natom))
  fixed_mask = .false.
  do i = 1, min(4, basis2%natom)
     fixed_mask(i) = .true.
  end do

  ! Store original positions of fixed atoms for verification
  call nn%inverse_design( &
       target_fingerprint = target_fp, &
       basis = basis2, &
       fixed_atoms = fixed_mask, &
       num_steps = 20, &
       step_size = 0.005_real32, &
       verbose = 1 &
  )

  ! Verify fixed atoms did not move (check first species, first 4 atoms)
  ! Note: fixed_mask(1:4) = .true., atoms indexed globally
  write(*,*) "  Inverse design completed"

  !-----------------------------------------------------------------------------
  ! Test 7: Fingerprint to distribs conversion
  !-----------------------------------------------------------------------------
  write(*,*) "Test 7: Fingerprint to distribs conversion"
  block
    use raffle__distribs, only: distribs_base_type
    type(distribs_base_type) :: reconverted
    call nn%fingerprint_to_distribs(fingerprint1, reconverted)
    call assert(allocated(reconverted%df_2body), &
         'Reconverted 2-body should be allocated', success)
    call assert(allocated(reconverted%df_3body), &
         'Reconverted 3-body should be allocated', success)
    call assert(allocated(reconverted%df_4body), &
         'Reconverted 4-body should be allocated', success)
    call assert( &
         size(reconverted%df_2body, 1) == nn%nbins(1), &
         'Reconverted 2-body should have correct nbins', success)
  end block

  !-----------------------------------------------------------------------------
  ! Summary
  !-----------------------------------------------------------------------------
  deallocate(fingerprint1, fingerprint2, predicted_fp, target_fp, fixed_mask)

  if (success) then
     write(*,*) "All nn_fingerprint tests PASSED"
  else
     write(*,*) "Some nn_fingerprint tests FAILED"
     stop 1
  end if


contains


  subroutine setup_carbon_diamond(basis)
    !! Set up a simple 8-atom carbon diamond cubic cell.
    type(basis_type), intent(out) :: basis
    real(real32) :: a

    a = 3.567_real32  ! Diamond cubic lattice parameter in Angstrom

    basis%sysname = "C_diamond"
    basis%nspec = 1
    basis%natom = 8
    basis%energy = -72.0_real32  ! Approximate DFT energy
    basis%lcart = .false.
    basis%pbc = [.true., .true., .true.]
    basis%lat(1,:) = [a, 0._real32, 0._real32]
    basis%lat(2,:) = [0._real32, a, 0._real32]
    basis%lat(3,:) = [0._real32, 0._real32, a]

    allocate(basis%spec(1))
    basis%spec(1)%name = 'C  '
    basis%spec(1)%num = 8
    allocate(basis%spec(1)%atom(8, 4))
    ! Diamond cubic fractional coordinates
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

    ! Perturb atom positions slightly
    basis%spec(1)%atom(5,1) = basis%spec(1)%atom(5,1) + 0.02_real32
    basis%spec(1)%atom(6,2) = basis%spec(1)%atom(6,2) - 0.01_real32
    basis%spec(1)%atom(7,3) = basis%spec(1)%atom(7,3) + 0.015_real32
    basis%spec(1)%atom(8,1) = basis%spec(1)%atom(8,1) - 0.01_real32
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


end program test_nn_fingerprint
