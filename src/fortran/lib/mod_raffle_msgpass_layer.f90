module raffle__msgpass_layer
  !! Module implementing an MLIP-style message passing layer for RAFFLE.
  !!
  !! This module implements a continuous-filter message passing neural network
  !! inspired by SchNet / MLIP architectures. It extends ATHENA's
  !! msgpass_layer_type to provide:
  !!
  !!   1. Gaussian RBF expansion of interatomic distances as edge features
  !!   2. Learned continuous filter (kernel MLP) on expanded edge features
  !!   3. Message aggregation using per-edge kernels
  !!   4. Residual node updates (when input/output dimensions match)
  !!   5. Multi-step message passing with graph-level readout
  !!
  !! Mathematical operation (per time step t):
  !! \[ h_i^{(t+1)} = h_i^{(t)} + \sigma\!\left(
  !!    \mathbf{W}^{(t)} h_i^{(t)} + b^{(t)}
  !!    + \sum_{j \in \mathcal{N}(i)}
  !!        \kappa_\theta^{(t)}(\phi_{\mathrm{RBF}}(r_{ij})) \, h_j^{(t)}
  !! \right) \]
  !!
  !! where:
  !!   - \(\phi_{\mathrm{RBF}}(r) = \exp(-\gamma_k (r - \mu_k)^2) \cdot f_c(r)\)
  !!     is a Gaussian radial basis function with cosine cutoff envelope
  !!   - \(\kappa_\theta\) is a learnable kernel MLP mapping RBF features
  !!     to per-edge transformation matrices
  !!   - \(\mathbf{W}^{(t)}\) is a linear (bypass) transform
  !!
  !! Graph readout:
  !! \[ h_{\mathrm{graph}} = \sum_{t=1}^{T} \sigma_r\!\left(
  !!    \mathbf{W}_r^{(t)} \sum_{i} h_i^{(t)} \right) \]
  !!
  !! The layer uses ATHENA's diffstruc autodiff framework through the
  !! gno_kernel_eval and gno_aggregate differentiable operations.
  use coreutils, only: real32, stop_program
  use graphstruc, only: graph_type
  use athena__base_layer, only: base_layer_type
  use athena__msgpass_layer, only: msgpass_layer_type
  use athena__misc_types, only: base_actv_type, base_init_type
  use diffstruc, only: array_type, sum, matmul, operator(+)
  use athena__diffstruc_extd, only: add_bias, gno_kernel_eval, gno_aggregate
  implicit none


  private

  public :: raffle_msgpass_layer_type


  real(real32), parameter :: PI_VAL = 4.0_real32 * atan(1.0_real32)


  type, extends(msgpass_layer_type) :: raffle_msgpass_layer_type
     !! MLIP-style message passing layer with RBF edge features
     !! and learned continuous filters.
     integer :: n_rbf = 20
     !! Number of Gaussian radial basis functions.
     real(real32) :: rbf_cutoff = 6.0_real32
     !! Cutoff radius for RBF expansion and cosine envelope.
     integer :: kernel_hidden = 64
     !! Hidden width of the kernel MLP.
     class(base_actv_type), allocatable :: activation_readout
     !! Activation function for readout.
     type(array_type), allocatable, dimension(:,:) :: z
     !! Intermediate node embeddings: z(t, s) for timestep t, sample s.
     type(array_type), allocatable, dimension(:) :: edge_rbf
     !! Persistent per-sample RBF-expanded edge features used during backprop.
   contains
     procedure, pass(this) :: get_num_params => get_num_params_raffle
     procedure, pass(this) :: set_hyperparams => set_hyperparams_raffle
     procedure, pass(this) :: init => init_raffle
     procedure, pass(this) :: print_to_unit => print_to_unit_raffle
     procedure, pass(this) :: read => read_raffle
     procedure, pass(this) :: update_message => update_message_raffle
     procedure, pass(this) :: update_readout => update_readout_raffle
  end type raffle_msgpass_layer_type

  interface raffle_msgpass_layer_type
     module procedure layer_setup_raffle
  end interface raffle_msgpass_layer_type


contains


!###############################################################################
  pure function get_num_params_raffle(this) result(num_params)
    !! Get the total number of learnable parameters.
    !!
    !! Per time step t (T total):
    !!   - Kernel MLP: H*n_rbf + H + F*H + F  (where F = F_v * F_v)
    !!   - Bypass W:   F_v * F_v
    !!   - Bias b:     F_v
    !! Readout per time step:
    !!   - Readout W:  num_outputs * F_v
    implicit none

    ! Arguments
    class(raffle_msgpass_layer_type), intent(in) :: this
    !! Layer instance.
    integer :: num_params
    !! Total number of learnable parameters.

    ! Local variables
    integer :: F_v, kh, nr, Fsq, nsteps
    !! Vertex feature count, hidden width, RBF count, kernel output, time steps.

    F_v    = this%num_vertex_features(0)
    kh     = this%kernel_hidden
    nr     = this%n_rbf
    Fsq    = F_v * F_v
    nsteps = this%num_time_steps

    ! Per time step: kernel_MLP + bypass_W + bias_b
    num_params = nsteps * (kh * nr + kh + Fsq * kh + Fsq + F_v * F_v + F_v)
    ! Readout weights per time step
    num_params = num_params + nsteps * this%num_outputs * F_v

  end function get_num_params_raffle
!###############################################################################


!###############################################################################
  function layer_setup_raffle( &
       num_vertex_features, num_edge_features, num_time_steps, &
       num_outputs, &
       n_rbf, rbf_cutoff, kernel_hidden, &
       message_activation, readout_activation, &
       kernel_initialiser, &
       verbose &
  ) result(layer)
    !! Construct a raffle_msgpass_layer_type.
    use athena__activation, only: activation_setup
    use athena__initialiser, only: initialiser_setup, get_default_initialiser
    implicit none

    ! Arguments
    integer, dimension(:), intent(in) :: num_vertex_features
    !! Number of vertex features (scalar array).
    integer, dimension(:), intent(in) :: num_edge_features
    !! Number of edge features (scalar array; typically [1] for distance).
    integer, intent(in) :: num_time_steps
    !! Number of message-passing iterations.
    integer, intent(in) :: num_outputs
    !! Dimension of graph-level output (fingerprint).
    integer, intent(in), optional :: n_rbf
    !! Number of RBF basis functions. Default: 20.
    real(real32), intent(in), optional :: rbf_cutoff
    !! Cutoff radius for RBF. Default: 6.0.
    integer, intent(in), optional :: kernel_hidden
    !! Hidden width of kernel MLP. Default: 64.
    class(*), intent(in), optional :: message_activation
    !! Message activation function specification (string or actv_type).
    class(*), intent(in), optional :: readout_activation
    !! Readout activation function specification.
    character(*), intent(in), optional :: kernel_initialiser
    !! Kernel initialiser name. Default: inferred from activation.
    integer, intent(in), optional :: verbose
    !! Verbosity level.
    type(raffle_msgpass_layer_type) :: layer
    !! Constructed layer.

    ! Local variables
    integer :: verbose_
    class(base_actv_type), allocatable :: msg_actv, readout_actv
    class(base_init_type), allocatable :: k_init
    character(len=256) :: buffer

    verbose_ = 0
    if (present(verbose)) verbose_ = verbose

    ! Activation functions
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

    ! Kernel initialiser
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

  end function layer_setup_raffle
!###############################################################################


!###############################################################################
  subroutine set_hyperparams_raffle( &
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

    ! Arguments
    class(raffle_msgpass_layer_type), intent(inout) :: this
    !! Layer instance.
    integer, dimension(:), intent(in) :: num_vertex_features
    !! Number of vertex features.
    integer, dimension(:), intent(in) :: num_edge_features
    !! Number of edge features.
    integer, intent(in) :: num_time_steps
    !! Number of message-passing time steps.
    integer, intent(in) :: num_outputs
    !! Dimension of graph-level output.
    integer, intent(in), optional :: n_rbf
    !! Number of RBF basis functions.
    real(real32), intent(in), optional :: rbf_cutoff
    !! Cutoff radius for RBF.
    integer, intent(in), optional :: kernel_hidden
    !! Hidden width of kernel MLP.
    class(base_actv_type), allocatable, intent(in) :: message_activation
    !! Message activation function.
    class(base_actv_type), allocatable, intent(in) :: readout_activation
    !! Readout activation function.
    class(base_init_type), allocatable, intent(in) :: kernel_initialiser
    !! Kernel initialiser.
    integer, intent(in), optional :: verbose
    !! Verbosity level.

    ! Local variables
    integer :: fv, nsteps
    character(len=256) :: buffer

    this%name = "raffle_mlip"
    this%type = "rmlip"
    this%input_rank = 2
    this%output_rank = 1
    this%use_graph_input = .true.
    this%use_graph_output = .false.
    this%use_bias = .true.
    this%num_outputs = num_outputs
    this%num_time_steps = num_time_steps

    ! MLIP-specific hyperparameters
    this%n_rbf = 20
    if (present(n_rbf)) this%n_rbf = n_rbf
    this%rbf_cutoff = 6.0_real32
    if (present(rbf_cutoff)) this%rbf_cutoff = rbf_cutoff
    this%kernel_hidden = 64
    if (present(kernel_hidden)) this%kernel_hidden = kernel_hidden

    ! Vertex and edge features
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

    ! Activation functions
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

    ! Kernel initialiser
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

    ! Compute num_params_msg and num_params_readout
    fv = this%num_vertex_features(0)
    if (allocated(this%num_params_msg)) deallocate(this%num_params_msg)
    allocate(this%num_params_msg(nsteps))
    this%num_params_msg = ( &
         this%kernel_hidden * this%n_rbf + this%kernel_hidden + &
         fv * fv * this%kernel_hidden + fv * fv + &
         fv * fv + fv &
    )
    this%num_params_readout = nsteps * num_outputs * fv

    this%output_shape = [num_outputs, 0]
    this%num_params = this%get_num_params()

    if (present(verbose)) then
       if (abs(verbose) > 0) then
          write(*,'("RAFFLE_MLIP message activation: ",A)') &
               trim(this%activation%name)
          write(*,'("RAFFLE_MLIP readout activation: ",A)') &
               trim(this%activation_readout%name)
          write(*,'("RAFFLE_MLIP n_rbf = ",I0)') this%n_rbf
          write(*,'("RAFFLE_MLIP rbf_cutoff = ",F8.3)') this%rbf_cutoff
          write(*,'("RAFFLE_MLIP kernel_hidden = ",I0)') this%kernel_hidden
       end if
    end if

  end subroutine set_hyperparams_raffle
!###############################################################################


!###############################################################################
  subroutine init_raffle(this, input_shape, verbose)
    !! Initialise the RAFFLE MLIP layer: allocate and initialise parameters.
    !!
    !! Parameter layout (for T time steps, F_v vertex features):
    !!   Message params per step t = 1..T (3 params each):
    !!     params(3*(t-1)+1): packed kernel MLP
    !!       [H*n_rbf + H + F*H + F, 1] where F = F_v * F_v
    !!     params(3*(t-1)+2): W bypass [F_v, F_v, 1]
    !!     params(3*(t-1)+3): b bias   [F_v, 1]
    !!   Readout params per step t = 1..T:
    !!     params(3*T + t): readout weights [num_outputs * F_v, 1]
    implicit none

    ! Arguments
    class(raffle_msgpass_layer_type), intent(inout) :: this
    !! Layer instance.
    integer, dimension(:), intent(in) :: input_shape
    !! Input shape: [num_vertex_features, 0].
    integer, intent(in), optional :: verbose
    !! Verbosity level.

    ! Local variables
    integer :: F_v, kh, nr, Fsq, nsteps, t, p, total_params
    integer :: kernel_size, off_U, off_bu, off_V, off_bv
    integer :: verbose_

    verbose_ = 0
    if (present(verbose)) verbose_ = verbose

    ! Set shapes
    if (.not. allocated(this%input_shape)) call this%set_shape(input_shape)

    F_v    = input_shape(1)
    kh     = this%kernel_hidden
    nr     = this%n_rbf
    Fsq    = F_v * F_v
    nsteps = this%num_time_steps

    kernel_size = kh * nr + kh + Fsq * kh + Fsq

    ! Total parameter arrays: 3 per step (kernel, W, b) + 1 readout per step
    total_params = 4 * nsteps

    if (allocated(this%weight_shape)) deallocate(this%weight_shape)
    if (allocated(this%params)) deallocate(this%params)
    allocate(this%weight_shape(2, total_params))
    allocate(this%params(total_params))

    ! Allocate and initialise per-step message parameters
    do t = 1, nsteps
       p = 3 * (t - 1)

       ! params(p+1): packed kernel MLP [kernel_size, 1]
       this%weight_shape(:, p+1) = [kernel_size, 1]
       call this%params(p+1)%allocate([kernel_size, 1])
       call this%params(p+1)%set_requires_grad(.true.)
       this%params(p+1)%fix_pointer = .true.
       this%params(p+1)%is_sample_dependent = .false.
       this%params(p+1)%is_temporary = .false.

       ! Initialise kernel MLP sub-blocks
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

       ! params(p+2): W bypass [F_v, F_v, 1]
       this%weight_shape(:, p+2) = [F_v, F_v]
       call this%params(p+2)%allocate([F_v, F_v, 1])
       call this%params(p+2)%set_requires_grad(.true.)
       this%params(p+2)%fix_pointer = .true.
       this%params(p+2)%is_sample_dependent = .false.
       this%params(p+2)%is_temporary = .false.
       call this%kernel_init%initialise( &
            this%params(p+2)%val(:, 1), &
            fan_in=F_v+1, fan_out=F_v, spacing=[F_v])

       ! params(p+3): bias [F_v, 1]
       this%weight_shape(:, p+3) = [F_v, 1]
       call this%params(p+3)%allocate([F_v, 1])
       call this%params(p+3)%set_requires_grad(.true.)
       this%params(p+3)%fix_pointer = .true.
       this%params(p+3)%is_sample_dependent = .false.
       this%params(p+3)%is_temporary = .false.
       call this%bias_init%initialise( &
            this%params(p+3)%val(:, 1), &
            fan_in=F_v+1, fan_out=F_v)
    end do

    ! Allocate and initialise readout parameters
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

    ! Clear any existing output
    if (allocated(this%output)) deallocate(this%output)

  end subroutine init_raffle
!###############################################################################


!###############################################################################
  subroutine expand_rbf_features(distances, n_rbf, cutoff, rbf_vals)
    !! Expand scalar distances to Gaussian RBF features with cosine cutoff.
    !!
    !! phi_k(r) = exp(-gamma * (r - mu_k)^2) * f_c(r)
    !! where mu_k = k * cutoff / (n_rbf - 1), k = 0, ..., n_rbf-1
    !!       gamma = n_rbf^2 / (2 * cutoff^2)
    !!       f_c(r) = 0.5 * (cos(pi * r / cutoff) + 1)  for r <= cutoff
    implicit none

    ! Arguments
    real(real32), dimension(:), intent(in) :: distances
    !! Interatomic distances [num_edges].
    integer, intent(in) :: n_rbf
    !! Number of RBF basis functions.
    real(real32), intent(in) :: cutoff
    !! Cutoff radius.
    real(real32), dimension(:,:), intent(out) :: rbf_vals
    !! Output RBF features [n_rbf, num_edges].

    ! Local variables
    integer :: k, e, num_e
    real(real32) :: mu_k, gamma, r, fc

    num_e = size(distances)
    gamma = real(n_rbf, real32)**2 / (2.0_real32 * cutoff**2)

    do e = 1, num_e
       r = distances(e)
       ! Cosine cutoff envelope
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

  end subroutine expand_rbf_features
!###############################################################################


!###############################################################################
  subroutine update_message_raffle(this, input)
    !! Update messages for the RAFFLE MLIP layer.
    !!
    !! For each sample s and time step t:
    !!   1. Expand edge distances to RBF features
    !!   2. Evaluate kernel MLP on RBF features
    !!   3. Aggregate neighbour messages using per-edge kernels
    !!   4. Add linear bypass and bias
    !!   5. Apply activation
    !!   6. Apply residual connection
    !!   7. Store intermediate embedding for readout
    !!
    !! input(1, s) = node features  [F_v, num_vertices]
    !! input(2, s) = edge features  [1, num_edges] (raw distances)
    implicit none

    ! Arguments
    class(raffle_msgpass_layer_type), intent(inout), target :: this
    !! Layer instance.
    class(array_type), dimension(:,:), intent(in), target :: input
    !! Input node-feature and edge-feature tensors.

    ! Local variables
    integer :: s, t, p, num_samples, num_edges, F_v, nsteps
    type(array_type), pointer :: ptr_h, ptr_kernels, ptr_agg
    type(array_type), pointer :: ptr_bypass, ptr_z, ptr_activated
    real(real32), allocatable :: rbf_vals(:,:)
    logical :: has_activation

    F_v    = this%num_vertex_features(0)
    nsteps = this%num_time_steps
    num_samples = size(input, 2)

    has_activation = .true.
    if (this%activation%name == 'none') has_activation = .false.

    ! Allocate intermediate storage
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
       ! Get edge distances and expand to RBF
       num_edges = size(input(2, s)%val, 2)
       allocate(rbf_vals(this%n_rbf, num_edges))
       call expand_rbf_features( &
            input(2, s)%val(1, :), this%n_rbf, this%rbf_cutoff, rbf_vals)

       ! Persist the expanded edge features for the lifetime of the autodiff graph.
       if (this%edge_rbf(s)%allocated) call this%edge_rbf(s)%deallocate()
       call this%edge_rbf(s)%allocate(source=rbf_vals)
       call this%edge_rbf(s)%set_requires_grad(.false.)
       this%edge_rbf(s)%is_temporary = .false.
       this%edge_rbf(s)%is_sample_dependent = .false.

       deallocate(rbf_vals)

       ! Initial node features
       ptr_h => input(1, s)

       do t = 1, nsteps
          p = 3 * (t - 1)

          ! Step 1: Evaluate kernel MLP on RBF-expanded edge features
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

          ! Step 3: Linear bypass W @ h
          ptr_bypass => matmul(this%params(p + 2), ptr_h)

          ! Step 4: Combine aggregation + bypass
          ptr_z => ptr_agg + ptr_bypass

          ! Step 5: Add bias
          ptr_z => add_bias( &
               ptr_z, this%params(p + 3), dim=1, dim_act_on_shape=.true.)

          ! Step 6: Apply activation
          if (has_activation) then
             ptr_activated => this%activation%apply(ptr_z)
          else
             ptr_activated => ptr_z
          end if

          ! Step 7: Propagate the updated embedding.
          ptr_h => ptr_activated

          ! Store intermediate embedding
          call this%z(t, s)%zero_grad()
          call this%z(t, s)%assign_and_deallocate_source(ptr_h)
          this%z(t, s)%is_temporary = .false.

          ! Point to stored embedding for next step
          ptr_h => this%z(t, s)
       end do
    end do

  end subroutine update_message_raffle
!###############################################################################


!###############################################################################
  subroutine update_readout_raffle(this)
    !! Graph-level readout by summing over nodes across all time steps.
    !!
    !! For each time step t and sample s:
    !!   readout += activation_readout(W_readout(t) @ z(t,s)) summed over nodes
    implicit none

    ! Arguments
    class(raffle_msgpass_layer_type), intent(inout), target :: this
    !! Layer instance.

    ! Local variables
    integer :: s, t, batch_size, nsteps
    type(array_type), pointer :: ptr1, ptr2, ptr3, ptr_params, ptr_z

    nsteps = this%num_time_steps
    batch_size = size(this%z, 2)

    ! Allocate output
    if (allocated(this%output)) deallocate(this%output)
    allocate(this%output(1, 1))

    call this%output(1, 1)%zero_grad()

    ptr3 => null()
    do t = 1, nsteps
       do s = 1, batch_size
          ptr_params => this%params(3 * nsteps + t)
          ptr_z => this%z(t, s)

          ! Apply readout weight matrix
          ptr1 => matmul(ptr_params, ptr_z)

          ! Apply readout activation
          ptr2 => this%activation_readout%apply(ptr1)

          ! Sum over nodes (dim=2), building batch dimension
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

  end subroutine update_readout_raffle
!###############################################################################


!###############################################################################
  subroutine print_to_unit_raffle(this, unit)
    !! Print layer settings and parameters to a unit for serialisation.
    use coreutils, only: to_upper
    implicit none

    ! Arguments
    class(raffle_msgpass_layer_type), intent(in) :: this
    !! Layer instance.
    integer, intent(in) :: unit
    !! Output unit number.

    ! Local variables
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

  end subroutine print_to_unit_raffle
!###############################################################################


!###############################################################################
  subroutine read_raffle(this, unit, verbose)
    !! Read layer from a file unit.
    use athena__tools_infile, only: assign_val, move
    use coreutils, only: to_lower, to_upper, icount
    use athena__activation, only: read_activation, activation_setup
    use athena__initialiser, only: initialiser_setup
    implicit none

    ! Arguments
    class(raffle_msgpass_layer_type), intent(inout) :: this
    !! Layer instance.
    integer, intent(in) :: unit
    !! Input unit number.
    integer, intent(in), optional :: verbose
    !! Verbosity level.

    ! Local variables
    integer :: stat, verbose_
    integer :: j, k, c, itmp1, iline
    integer :: num_vertex_features, num_edge_features
    integer :: num_time_steps, num_outputs
    integer :: n_rbf, kernel_hidden
    real(real32) :: rbf_cutoff
    class(base_actv_type), allocatable :: msg_activation, readout_activation
    class(base_init_type), allocatable :: k_init
    character(256) :: buffer, tag, err_msg
    real(real32), allocatable, dimension(:) :: data_list
    integer :: param_line, final_line, num_vals, p

    verbose_ = 0
    if (present(verbose)) verbose_ = verbose

    ! Defaults
    n_rbf = 20
    rbf_cutoff = 6.0_real32
    kernel_hidden = 64

    iline = 0
    param_line = 0
    final_line = 0

    tag_loop: do
       read(unit, '(A)', iostat=stat) buffer
       if (stat /= 0) then
          write(err_msg, &
               '("file encountered error (EoF?) before END ",A)') &
               to_upper(this%name)
          call stop_program(err_msg)
          return
       end if
       if (trim(adjustl(buffer)) == "") cycle tag_loop

       if (trim(adjustl(buffer)) == &
            "END "//to_upper(trim(this%name))) then
          final_line = iline
          backspace(unit)
          exit tag_loop
       end if
       iline = iline + 1

       tag = trim(adjustl(buffer))
       if (scan(buffer, "=") /= 0) tag = trim(tag(:scan(tag, "=") - 1))

       select case(trim(tag))
       case("NUM_VERTEX_FEATURES")
          call assign_val(buffer, num_vertex_features, itmp1)
       case("NUM_EDGE_FEATURES")
          call assign_val(buffer, num_edge_features, itmp1)
       case("NUM_TIME_STEPS")
          call assign_val(buffer, num_time_steps, itmp1)
       case("NUM_OUTPUTS")
          call assign_val(buffer, num_outputs, itmp1)
       case("N_RBF")
          call assign_val(buffer, n_rbf, itmp1)
       case("RBF_CUTOFF")
          call assign_val(buffer, rbf_cutoff, itmp1)
       case("KERNEL_HIDDEN")
          call assign_val(buffer, kernel_hidden, itmp1)
       case("ACTIVATION")
          iline = iline - 1
          backspace(unit)
          msg_activation = read_activation(unit, iline)
       case("READOUT")
          iline = iline - 1
          backspace(unit)
          readout_activation = read_activation(unit, iline)
       case("WEIGHTS")
          param_line = iline
       case default
          if (scan(to_lower(trim(adjustl(buffer))), &
               'abcdfghijklmnopqrstuvwxyz') == 0) then
             cycle tag_loop
          elseif (tag(:3) == 'END') then
             cycle tag_loop
          end if
          write(err_msg, '("Unrecognised line in input file: ",A)') &
               trim(adjustl(buffer))
          call stop_program(err_msg)
          return
       end select
    end do tag_loop

    ! Set up the layer with parsed hyperparameters
    if (.not. allocated(msg_activation)) &
         msg_activation = activation_setup('swish')
    if (.not. allocated(readout_activation)) &
         readout_activation = activation_setup('none')

    k_init = initialiser_setup('zeros')

    call this%set_hyperparams( &
         num_vertex_features=[num_vertex_features], &
         num_edge_features=[num_edge_features], &
         num_time_steps=num_time_steps, &
         num_outputs=num_outputs, &
         n_rbf=n_rbf, &
         rbf_cutoff=rbf_cutoff, &
         kernel_hidden=kernel_hidden, &
         message_activation=msg_activation, &
         readout_activation=readout_activation, &
         kernel_initialiser=k_init, &
         verbose=verbose_)
    call this%init(input_shape=[num_vertex_features, 0])

    ! Read weights
    if (param_line /= 0) then
       call move(unit, param_line - iline, iostat=stat)
       do p = 1, size(this%params)
          num_vals = size(this%params(p)%val(:, 1))
          allocate(data_list(num_vals), source=0._real32)
          c = 1
          do while (c <= num_vals)
             read(unit, '(A)', iostat=stat) buffer
             if (stat /= 0) exit
             k = icount(buffer)
             read(buffer, *, iostat=stat) (data_list(j), j=c, c+k-1)
             c = c + k
          end do
          this%params(p)%val(:, 1) = data_list
          deallocate(data_list)
       end do

       read(unit, '(A)') buffer
       if (trim(adjustl(buffer)) /= "END WEIGHTS") then
          call stop_program("END WEIGHTS not where expected")
          return
       end if
    end if

    call move(unit, final_line - iline, iostat=stat)
    read(unit, '(A)') buffer
    if (trim(adjustl(buffer)) /= &
         "END "//to_upper(trim(this%name))) then
       write(err_msg, '("END ",A," not where expected")') &
            to_upper(this%name)
       call stop_program(err_msg)
       return
    end if

  end subroutine read_raffle
!###############################################################################


end module raffle__msgpass_layer
