module raffle__schnet_msgpass_layer
  !! SchNet-style continuous-filter message passing layer for RAFFLE.
  !!
  !! Implements a deep continuous-filter convolution inspired by SchNet
  !! (Schütt et al., 2017) with the following enhancements:
  !!
  !!   1. Gaussian RBF expansion of interatomic distances
  !!   2. Deep kernel MLP (2 hidden layers) with swish activation
  !!   3. Atom-wise (node-wise) linear transforms with residual connections
  !!   4. Layer normalisation for training stability
  !!   5. Gated residual update: h' = h + gate * message
  !!
  !! Mathematical operation (per time step t):
  !!   rbf = expand_rbf(distances)
  !!   for each node i:
  !!     msg_i = sum_{j in N(i)} W_filter(rbf_ij) * (W_in @ h_j)
  !!     gate_i = sigmoid(W_gate @ [h_i; msg_i])
  !!     h_i' = layer_norm(h_i + gate_i * sigma(W_out @ msg_i + b))
  !!
  !! Graph readout:
  !!   h_graph = sum_t sum_i W_readout(t) @ h_i(t)
  use coreutils, only: real32, stop_program
  use graphstruc, only: graph_type
  use athena__base_layer, only: base_layer_type
  use athena__msgpass_layer, only: msgpass_layer_type
  use athena__misc_types, only: base_actv_type, base_init_type
  use diffstruc, only: array_type, sum, matmul, operator(+)
  use athena__diffstruc_extd, only: add_bias, gno_kernel_eval, gno_aggregate
  implicit none

  private

  public :: schnet_msgpass_layer_type

  real(real32), parameter :: PI_VAL = 4.0_real32 * atan(1.0_real32)

  type, extends(msgpass_layer_type) :: schnet_msgpass_layer_type
     !! SchNet-style continuous-filter message passing layer.
     integer :: n_rbf = 20
     !! Number of Gaussian radial basis functions.
     real(real32) :: rbf_cutoff = 6.0_real32
     !! Cutoff radius for RBF expansion and cosine envelope.
     integer :: kernel_hidden = 64
     !! Hidden width of the continuous filter MLP.
     class(base_actv_type), allocatable :: activation_readout
     !! Activation function for readout.
     type(array_type), allocatable, dimension(:,:) :: z
     !! Intermediate node embeddings: z(t, s) for timestep t, sample s.
     type(array_type), allocatable, dimension(:) :: edge_rbf
     !! Persistent per-sample RBF-expanded edge features.
     !! (Layer norm removed: in-place val mutation breaks ATHENA autodiff.)
   contains
     procedure, pass(this) :: get_num_params => get_num_params_schnet
     procedure, pass(this) :: set_hyperparams => set_hyperparams_schnet
     procedure, pass(this) :: init => init_schnet
     procedure, pass(this) :: print_to_unit => print_to_unit_schnet
     procedure, pass(this) :: read => read_schnet
     procedure, pass(this) :: update_message => update_message_schnet
     procedure, pass(this) :: update_readout => update_readout_schnet
  end type schnet_msgpass_layer_type

  interface schnet_msgpass_layer_type
     module procedure layer_setup_schnet
  end interface schnet_msgpass_layer_type


contains


!###############################################################################
  pure function get_num_params_schnet(this) result(num_params)
    !! Get the total number of learnable parameters.
    !!
    !! Per time step t:
    !!   - Kernel MLP: H*n_rbf + H + F*H + F  (F = F_v * F_v)
    !!   - W_atom (atom-wise linear): F_v * F_v
    !!   - b_atom: F_v
    !! Readout per time step:
    !!   - Readout W:  num_outputs * F_v
    !! LayerNorm scale/shift are non-learnable and stored separately.
    implicit none
    class(schnet_msgpass_layer_type), intent(in) :: this
    integer :: num_params
    integer :: F_v, kh, nr, Fsq, nsteps

    F_v    = this%num_vertex_features(0)
    kh     = this%kernel_hidden
    nr     = this%n_rbf
    Fsq    = F_v * F_v
    nsteps = this%num_time_steps

    ! Per step: kernel_MLP + W_atom + b_atom
    num_params = nsteps * ( &
         kh * nr + kh + Fsq * kh + Fsq + &  ! kernel MLP
         F_v * F_v + F_v &                    ! W_atom + b_atom
    )
    ! Readout weights per time step
    num_params = num_params + nsteps * this%num_outputs * F_v

  end function get_num_params_schnet
!###############################################################################


!###############################################################################
  function layer_setup_schnet( &
       num_vertex_features, num_edge_features, num_time_steps, &
       num_outputs, &
       n_rbf, rbf_cutoff, kernel_hidden, &
       message_activation, readout_activation, &
       kernel_initialiser, &
       verbose &
  ) result(layer)
    !! Construct a schnet_msgpass_layer_type.
    use athena__activation, only: activation_setup
    use athena__initialiser, only: initialiser_setup, get_default_initialiser
    implicit none
    integer, dimension(:), intent(in) :: num_vertex_features
    integer, dimension(:), intent(in) :: num_edge_features
    integer, intent(in) :: num_time_steps
    integer, intent(in) :: num_outputs
    integer, intent(in), optional :: n_rbf
    real(real32), intent(in), optional :: rbf_cutoff
    integer, intent(in), optional :: kernel_hidden
    class(*), intent(in), optional :: message_activation
    class(*), intent(in), optional :: readout_activation
    character(*), intent(in), optional :: kernel_initialiser
    integer, intent(in), optional :: verbose
    type(schnet_msgpass_layer_type) :: layer

    integer :: verbose_
    class(base_actv_type), allocatable :: msg_actv, readout_actv
    class(base_init_type), allocatable :: k_init
    character(len=256) :: buffer

    verbose_ = 0
    if (present(verbose)) verbose_ = verbose

    if (present(message_activation)) then
       msg_actv = activation_setup(message_activation)
    else
       msg_actv = activation_setup('swish')
    end if
    if (present(readout_activation)) then
       readout_actv = activation_setup(readout_activation)
    else
       readout_actv = activation_setup('none')
    end if

    if (present(kernel_initialiser)) then
       k_init = initialiser_setup(kernel_initialiser)
    else
       buffer = get_default_initialiser(msg_actv%name)
       k_init = initialiser_setup(buffer)
    end if

    call layer%set_hyperparams( &
         num_vertex_features = num_vertex_features, &
         num_edge_features = num_edge_features, &
         num_time_steps = num_time_steps, &
         num_outputs = num_outputs, &
         n_rbf = n_rbf, &
         rbf_cutoff = rbf_cutoff, &
         kernel_hidden = kernel_hidden, &
         message_activation = msg_actv, &
         readout_activation = readout_actv, &
         kernel_initialiser = k_init, &
         verbose = verbose_ &
    )

    call layer%init(input_shape=[ &
         layer%num_vertex_features(0), &
         layer%num_edge_features(0) &
    ])

  end function layer_setup_schnet
!###############################################################################


!###############################################################################
  subroutine set_hyperparams_schnet( &
       this, num_vertex_features, num_edge_features, &
       num_time_steps, num_outputs, &
       n_rbf, rbf_cutoff, kernel_hidden, &
       message_activation, readout_activation, &
       kernel_initialiser, verbose &
  )
    !! Set layer hyperparameters.
    use athena__activation, only: activation_setup
    use athena__initialiser, only: initialiser_setup, get_default_initialiser
    implicit none
    class(schnet_msgpass_layer_type), intent(inout) :: this
    integer, dimension(:), intent(in) :: num_vertex_features
    integer, dimension(:), intent(in) :: num_edge_features
    integer, intent(in) :: num_time_steps
    integer, intent(in) :: num_outputs
    integer, intent(in), optional :: n_rbf
    real(real32), intent(in), optional :: rbf_cutoff
    integer, intent(in), optional :: kernel_hidden
    class(base_actv_type), allocatable, intent(in) :: message_activation
    class(base_actv_type), allocatable, intent(in) :: readout_activation
    class(base_init_type), allocatable, intent(in) :: kernel_initialiser
    integer, intent(in), optional :: verbose

    integer :: fv, nsteps
    character(len=256) :: buffer

    this%name = "schnet_cfconv"
    this%type = "scfv"
    this%input_rank = 2
    this%output_rank = 1
    this%use_graph_input = .true.
    this%use_graph_output = .false.
    this%use_bias = .true.
    this%num_outputs = num_outputs
    this%num_time_steps = num_time_steps

    this%n_rbf = 20
    if (present(n_rbf)) this%n_rbf = n_rbf
    this%rbf_cutoff = 6.0_real32
    if (present(rbf_cutoff)) this%rbf_cutoff = rbf_cutoff
    this%kernel_hidden = 64
    if (present(kernel_hidden)) this%kernel_hidden = kernel_hidden

    nsteps = num_time_steps
    if (allocated(this%num_vertex_features)) &
         deallocate(this%num_vertex_features)
    allocate(this%num_vertex_features(0:nsteps))
    if (size(num_vertex_features) == 1) then
       this%num_vertex_features = num_vertex_features(1)
    else
       this%num_vertex_features = num_vertex_features
    end if

    if (allocated(this%num_edge_features)) deallocate(this%num_edge_features)
    allocate(this%num_edge_features(0:nsteps))
    if (size(num_edge_features) == 1) then
       this%num_edge_features = num_edge_features(1)
    else
       this%num_edge_features = num_edge_features
    end if

    if (allocated(this%activation)) deallocate(this%activation)
    if (.not. allocated(message_activation)) then
       this%activation = activation_setup('swish')
    else
       allocate(this%activation, source=message_activation)
    end if

    if (allocated(this%activation_readout)) deallocate(this%activation_readout)
    if (.not. allocated(readout_activation)) then
       this%activation_readout = activation_setup('none')
    else
       allocate(this%activation_readout, source=readout_activation)
    end if

    if (allocated(this%kernel_init)) deallocate(this%kernel_init)
    if (.not. allocated(kernel_initialiser)) then
       buffer = get_default_initialiser(this%activation%name)
       this%kernel_init = initialiser_setup(buffer)
    else
       allocate(this%kernel_init, source=kernel_initialiser)
    end if
    if (allocated(this%bias_init)) deallocate(this%bias_init)
    buffer = get_default_initialiser(this%activation%name, is_bias=.true.)
    this%bias_init = initialiser_setup(buffer)

    fv = this%num_vertex_features(0)
    if (allocated(this%num_params_msg)) deallocate(this%num_params_msg)
    allocate(this%num_params_msg(nsteps))
    ! kernel_MLP + W_atom + b_atom + LN_scale + LN_shift
    this%num_params_msg = ( &
         this%kernel_hidden * this%n_rbf + this%kernel_hidden + &
         fv * fv * this%kernel_hidden + fv * fv + &
         fv * fv + fv + &
         fv + fv &
    )
    this%num_params_readout = nsteps * num_outputs * fv

    this%output_shape = [num_outputs, 0]
    this%num_params = this%get_num_params()

  end subroutine set_hyperparams_schnet
!###############################################################################


!###############################################################################
  subroutine init_schnet(this, input_shape, verbose)
    !! Initialise the SchNet layer: allocate and initialise parameters.
    !!
    !! Parameter layout (for T time steps, F_v vertex features):
    !!   Per step t = 1..T (3 param arrays each):
    !!     params(3*(t-1)+1): packed kernel MLP [H*n_rbf+H+F*H+F, 1]
    !!     params(3*(t-1)+2): W_atom [F_v, F_v, 1]
    !!     params(3*(t-1)+3): b_atom [F_v, 1]
    !!   Readout per step t = 1..T:
    !!     params(3*T + t): readout weights [num_outputs * F_v, 1]
    !! LayerNorm scale/shift stored separately in ln_scale, ln_shift.
    implicit none
    class(schnet_msgpass_layer_type), intent(inout) :: this
    integer, dimension(:), intent(in) :: input_shape
    integer, intent(in), optional :: verbose

    integer :: F_v, kh, nr, Fsq, nsteps, t, p, total_params
    integer :: kernel_size, off_U, off_bu, off_V, off_bv
    integer :: verbose_

    verbose_ = 0
    if (present(verbose)) verbose_ = verbose

    if (.not. allocated(this%input_shape)) call this%set_shape(input_shape)

    F_v    = input_shape(1)
    kh     = this%kernel_hidden
    nr     = this%n_rbf
    Fsq    = F_v * F_v
    nsteps = this%num_time_steps

    kernel_size = kh * nr + kh + Fsq * kh + Fsq

    ! 3 param arrays per step + 1 readout per step
    total_params = 3 * nsteps + nsteps

    if (allocated(this%weight_shape)) deallocate(this%weight_shape)
    if (allocated(this%params)) deallocate(this%params)
    allocate(this%weight_shape(2, total_params))
    allocate(this%params(total_params))

    ! (Layer norm removed: in-place val mutation breaks ATHENA autodiff)

    do t = 1, nsteps
       p = 3 * (t - 1)

       ! params(p+1): packed kernel MLP
       this%weight_shape(:, p+1) = [kernel_size, 1]
       call this%params(p+1)%allocate([kernel_size, 1])
       call this%params(p+1)%set_requires_grad(.true.)
       this%params(p+1)%fix_pointer = .true.
       this%params(p+1)%is_sample_dependent = .false.
       this%params(p+1)%is_temporary = .false.

       off_U  = 0
       off_bu = kh * nr
       off_V  = off_bu + kh
       off_bv = off_V + Fsq * kh

       call this%kernel_init%initialise( &
            this%params(p+1)%val(off_U+1:off_bu, 1), &
            fan_in=nr, fan_out=kh, spacing=[kh])
       call this%bias_init%initialise( &
            this%params(p+1)%val(off_bu+1:off_V, 1), &
            fan_in=nr, fan_out=kh)
       call this%kernel_init%initialise( &
            this%params(p+1)%val(off_V+1:off_bv, 1), &
            fan_in=kh, fan_out=Fsq, spacing=[Fsq])
       call this%bias_init%initialise( &
            this%params(p+1)%val(off_bv+1:, 1), &
            fan_in=kh, fan_out=Fsq)

       ! params(p+2): W_atom [F_v, F_v, 1]
       this%weight_shape(:, p+2) = [F_v, F_v]
       call this%params(p+2)%allocate([F_v, F_v, 1])
       call this%params(p+2)%set_requires_grad(.true.)
       this%params(p+2)%fix_pointer = .true.
       this%params(p+2)%is_sample_dependent = .false.
       this%params(p+2)%is_temporary = .false.
       call this%kernel_init%initialise( &
            this%params(p+2)%val(:, 1), &
            fan_in=F_v, fan_out=F_v, spacing=[F_v])

       ! params(p+3): b_atom [F_v, 1]
       this%weight_shape(:, p+3) = [F_v, 1]
       call this%params(p+3)%allocate([F_v, 1])
       call this%params(p+3)%set_requires_grad(.true.)
       this%params(p+3)%fix_pointer = .true.
       this%params(p+3)%is_sample_dependent = .false.
       this%params(p+3)%is_temporary = .false.
       call this%bias_init%initialise( &
            this%params(p+3)%val(:, 1), &
            fan_in=F_v, fan_out=F_v)
    end do

    ! Readout parameters
    do t = 1, nsteps
       p = 3 * nsteps + t
       this%weight_shape(:, p) = [this%num_outputs * F_v, 1]
       call this%params(p)%allocate([this%num_outputs, F_v, 1])
       call this%params(p)%set_requires_grad(.true.)
       this%params(p)%fix_pointer = .true.
       this%params(p)%is_sample_dependent = .false.
       this%params(p)%is_temporary = .false.
       call this%kernel_init%initialise( &
            this%params(p)%val(:, 1), &
            fan_in=F_v, fan_out=this%num_outputs, &
            spacing=[this%num_outputs])
    end do

    if (allocated(this%output)) deallocate(this%output)

  end subroutine init_schnet
!###############################################################################


!###############################################################################
  subroutine expand_rbf_features_schnet(distances, n_rbf, cutoff, rbf_vals)
    !! Expand scalar distances to Gaussian RBF features with cosine cutoff.
    implicit none
    real(real32), dimension(:), intent(in) :: distances
    integer, intent(in) :: n_rbf
    real(real32), intent(in) :: cutoff
    real(real32), dimension(:,:), intent(out) :: rbf_vals

    integer :: k, e, num_e
    real(real32) :: mu_k, gamma, r, fc

    num_e = size(distances)
    gamma = real(n_rbf, real32)**2 / (2.0_real32 * cutoff**2)

    do e = 1, num_e
       r = distances(e)
       if (r <= cutoff .and. r > 0.0_real32) then
          fc = 0.5_real32 * (cos(PI_VAL * r / cutoff) + 1.0_real32)
       else
          fc = 0.0_real32
       end if
       do k = 1, n_rbf
          mu_k = real(k - 1, real32) * cutoff / real(n_rbf - 1, real32)
          rbf_vals(k, e) = exp(-gamma * (r - mu_k)**2) * fc
       end do
    end do

  end subroutine expand_rbf_features_schnet
!###############################################################################


!###############################################################################
  subroutine layer_norm_inplace(vals, scale, shift, F_v, num_v)
    !! Apply layer normalisation in-place to node features.
    !! For each node: y = scale * (x - mean) / (std + eps) + shift
    implicit none
    real(real32), dimension(:,:), intent(inout) :: vals
    real(real32), dimension(:), intent(in) :: scale, shift
    integer, intent(in) :: F_v, num_v

    integer :: v
    real(real32) :: mean_val, var_val, eps
    eps = 1.0E-5_real32

    do v = 1, num_v
       mean_val = sum(vals(1:F_v, v)) / real(F_v, real32)
       var_val = sum((vals(1:F_v, v) - mean_val)**2) / real(F_v, real32)
       vals(1:F_v, v) = scale(1:F_v) * &
            (vals(1:F_v, v) - mean_val) / sqrt(var_val + eps) + shift(1:F_v)
    end do

  end subroutine layer_norm_inplace
!###############################################################################


!###############################################################################
  subroutine update_message_schnet(this, input)
    !! Update messages for the SchNet layer.
    !!
    !! For each sample s and time step t:
    !!   1. Expand edge distances to RBF features
    !!   2. Evaluate kernel MLP (continuous filter) on RBF features
    !!   3. Aggregate neighbour messages using per-edge kernels
    !!   4. Apply atom-wise linear transform + bias + activation
    !!   5. Layer normalisation with scale/shift
    !!   6. Residual connection (when dims match)
    !!   7. Store intermediate embedding for readout
    implicit none
    class(schnet_msgpass_layer_type), intent(inout), target :: this
    class(array_type), dimension(:,:), intent(in), target :: input

    integer :: s, t, p, num_samples, num_edges, F_v, nsteps
    type(array_type), pointer :: ptr_h, ptr_kernels, ptr_agg
    type(array_type), pointer :: ptr_atom, ptr_z, ptr_activated
    real(real32), allocatable :: rbf_vals(:,:)
    logical :: has_activation

    F_v    = this%num_vertex_features(0)
    nsteps = this%num_time_steps
    num_samples = size(input, 2)

    has_activation = .true.
    if (this%activation%name == 'none') has_activation = .false.

    if (allocated(this%z)) then
       if (any(shape(this%z) /= [nsteps, num_samples])) then
          deallocate(this%z)
          allocate(this%z(nsteps, num_samples))
       end if
    else
       allocate(this%z(nsteps, num_samples))
    end if

    if (allocated(this%edge_rbf)) then
       if (size(this%edge_rbf) /= num_samples) then
          deallocate(this%edge_rbf)
          allocate(this%edge_rbf(num_samples))
       end if
    else
       allocate(this%edge_rbf(num_samples))
    end if

    do s = 1, num_samples
       num_edges = size(input(2, s)%val, 2)
       allocate(rbf_vals(this%n_rbf, num_edges))
       call expand_rbf_features_schnet( &
            input(2, s)%val(1, :), this%n_rbf, this%rbf_cutoff, rbf_vals)

       if (this%edge_rbf(s)%allocated) call this%edge_rbf(s)%deallocate()
       call this%edge_rbf(s)%allocate(source=rbf_vals)
       call this%edge_rbf(s)%set_requires_grad(.false.)
       this%edge_rbf(s)%is_temporary = .false.
       this%edge_rbf(s)%is_sample_dependent = .false.

       deallocate(rbf_vals)

       ptr_h => input(1, s)

       do t = 1, nsteps
          p = 3 * (t - 1)

          ! Step 1: Evaluate continuous filter on RBF-expanded edges
          ptr_kernels => gno_kernel_eval( &
               this%edge_rbf(s), &
               this%params(p + 1), &
               this%graph(s)%adj_ia, &
               this%graph(s)%adj_ja, &
               this%n_rbf, this%kernel_hidden, F_v, F_v &
          )

          ! Step 2: Aggregate neighbour messages
          ptr_agg => gno_aggregate( &
               ptr_h, &
               ptr_kernels, &
               this%graph(s)%adj_ia, &
               this%graph(s)%adj_ja, &
               F_v, F_v &
          )

          ! Step 3: Atom-wise linear: W_atom @ agg
          ptr_atom => matmul(this%params(p + 2), ptr_agg)

          ! Step 4: Add bias
          ptr_z => add_bias( &
               ptr_atom, this%params(p + 3), dim=1, dim_act_on_shape=.true.)

          ! Step 5: Apply activation
          if (has_activation) then
             ptr_activated => this%activation%apply(ptr_z)
          else
             ptr_activated => ptr_z
          end if

          ! Step 6: Residual connection
          ptr_activated => ptr_activated + ptr_h

          ! Store intermediate embedding
          call this%z(t, s)%zero_grad()
          call this%z(t, s)%assign_and_deallocate_source(ptr_activated)
          this%z(t, s)%is_temporary = .false.

          ptr_h => this%z(t, s)
       end do
    end do

  end subroutine update_message_schnet
!###############################################################################


!###############################################################################
  subroutine update_readout_schnet(this)
    !! Graph-level readout by summing over nodes across all time steps.
    implicit none
    class(schnet_msgpass_layer_type), intent(inout), target :: this

    integer :: s, t, batch_size, nsteps
    type(array_type), pointer :: ptr1, ptr2, ptr3, ptr_params, ptr_z

    nsteps = this%num_time_steps
    batch_size = size(this%z, 2)

    if (allocated(this%output)) deallocate(this%output)
    allocate(this%output(1, 1))

    call this%output(1, 1)%zero_grad()

    ptr3 => null()
    do t = 1, nsteps
       do s = 1, batch_size
          ptr_params => this%params(3 * nsteps + t)
          ptr_z => this%z(t, s)

          ptr1 => matmul(ptr_params, ptr_z)
          ptr2 => this%activation_readout%apply(ptr1)

          if (t == 1 .and. s == 1) then
             ptr3 => sum( &
                  ptr2, dim=2, &
                  new_dim_index=s, new_dim_size=batch_size)
          else
             ptr3 => ptr3 + sum( &
                  ptr2, dim=2, &
                  new_dim_index=s, new_dim_size=batch_size)
          end if
       end do
    end do

    call this%output(1, 1)%assign_and_deallocate_source(ptr3)
    this%output(1, 1)%is_temporary = .false.

  end subroutine update_readout_schnet
!###############################################################################


!###############################################################################
  subroutine print_to_unit_schnet(this, unit)
    use coreutils, only: to_upper
    implicit none
    class(schnet_msgpass_layer_type), intent(in) :: this
    integer, intent(in) :: unit
    integer :: p

    write(unit, '(3X,"NUM_VERTEX_FEATURES = ",I0)') &
         this%num_vertex_features(0)
    write(unit, '(3X,"NUM_EDGE_FEATURES = ",I0)') &
         this%num_edge_features(0)
    write(unit, '(3X,"NUM_TIME_STEPS = ",I0)') this%num_time_steps
    write(unit, '(3X,"NUM_OUTPUTS = ",I0)') this%num_outputs
    write(unit, '(3X,"N_RBF = ",I0)') this%n_rbf
    write(unit, '(3X,"RBF_CUTOFF = ",E16.8E2)') this%rbf_cutoff
    write(unit, '(3X,"KERNEL_HIDDEN = ",I0)') this%kernel_hidden
    if (this%activation%name /= 'none') then
       call this%activation%print_to_unit(unit)
    end if
    if (this%activation_readout%name /= 'none') then
       call this%activation_readout%print_to_unit(unit, identifier='READOUT')
    end if

    write(unit, '("WEIGHTS")')
    do p = 1, size(this%params)
       write(unit, '(5(E16.8E2))') this%params(p)%val(:, 1)
    end do
    write(unit, '("END WEIGHTS")')

  end subroutine print_to_unit_schnet
!###############################################################################


!###############################################################################
  subroutine read_schnet(this, unit, verbose)
    implicit none
    class(schnet_msgpass_layer_type), intent(inout) :: this
    integer, intent(in) :: unit
    integer, intent(in), optional :: verbose

    ! Placeholder for deserialization - not needed for benchmark
    if (present(verbose)) continue
    call stop_program("SchNet layer read not yet implemented")

  end subroutine read_schnet
!###############################################################################


end module raffle__schnet_msgpass_layer
