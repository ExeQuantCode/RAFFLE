module machine_learning
  use raffle__constants, only: real32
  use raffle__misc_linalg, only: modu
  use raffle__geom_rw, only: basis_type
#ifdef ENABLE_ATHENA
  use athena
  use athena, only: graph_type, edge_type
#endif
  implicit none


  private

#ifdef ENABLE_ATHENA

  public :: network_setup
  public :: network_train, network_train_graph
  public :: network_predict, network_predict_graph
  public :: get_graph_from_basis

  
  type(network_type) :: network
#endif


contains

#ifdef ENABLE_ATHENA
  subroutine network_setup(num_inputs, num_outputs)
    implicit none
    integer, intent(in) :: num_inputs, num_outputs

    type(metric_dict_type), dimension(2) :: metric_dict

    ! call network%add(conv1d_layer_type( &
    !      input_shape = [num_inputs, 1], &
    !      num_filters = 16, &
    !      kernel_size = 3, &
    !      activation_function = 'tanh' &
    !      ))
    ! ! call network%add(full_layer_type( &
    ! !      num_inputs = num_inputs, &
    ! !      num_outputs = 500, &
    ! !      activation_function = 'tanh' &
    ! !      ))
    ! call network%add(full_layer_type( &
    !      num_outputs = 200, &
    !      activation_function = 'tanh' &
    !      ))
    ! call network%add(full_layer_type( &
    !      num_outputs = num_outputs, &
    !      activation_function = 'tanh' &
    !      ))
    ! call network%compile( &
    !      optimiser = sgd_optimiser_type(learning_rate=1.E-3_real32, momentum=0.99_real32,&
    !                                     clip_dict=clip_type(clip_min=0._real32, &
    !                                                         clip_max=1.E-4_real32), &
    !                                     regulariser=l1l2_regulariser_type(l1=1.E-3_real32, &
    !                                                                       l2=1.E-3_real32) &
    !                                     ), &
    !      loss_method = 'hubber', accuracy_method = 'rmse', verbose=1, &
    !      metrics = ['accuracy'] &
    !      )
    ! call network%set_batch_size(12)
    
     call network%add(conv_mpnn_layer_type( &
          num_time_steps=4, &
          num_vertex_features=2, num_edge_features=1, &
          num_outputs=10, &
          max_vertex_degree = 8, &
          batch_size=1 ))
     call network%add(full_layer_type( &
          num_inputs  = 10, &
          num_outputs = 1, &
          batch_size  = 1, &
          activation_function='leaky_relu', &
          kernel_initialiser='he_normal', &
          bias_initialiser='ones' &
          ))

    metric_dict%active = .false.
    metric_dict(1)%key = "loss"
    metric_dict(2)%key = "accuracy"
    metric_dict%threshold = 1.E-1_real32
    call network%compile(optimiser=sgd_optimiser_type(), &
         loss_method="mse", metrics=metric_dict, &
         batch_size = 1, verbose = 0)
    call network%set_batch_size(1)

    network%metrics%threshold = 1.E-3_real32

  end subroutine network_setup

  subroutine network_train(x, y, num_epochs)
    implicit none
    real(real32), dimension(:,:), intent(in) :: x
    real(real32), dimension(:), intent(in) :: y
    integer, intent(in) :: num_epochs
    real(real32), dimension(:,:), allocatable :: y_renorm


    allocate(y_renorm(1,size(y,1)))
    y_renorm(1,:) = -1._real32 * y

    call renormalise_norm(y_renorm(1,:))
    write(*,*) y_renorm(1,:)
    call network%train(x, y_renorm, num_epochs=num_epochs, plateau_threshold=1.E-3_real32)
    !write(*,*) y_renorm


  end subroutine network_train

  subroutine network_train_graph(graphs, labels, num_epochs)
    implicit none
    type(graph_type), dimension(:), intent(in) :: graphs
    real(real32), dimension(:), intent(in) :: labels
    integer, intent(in) :: num_epochs

    integer :: n, s
    integer, dimension(:), allocatable :: sample_list

    do n = 1, num_epochs
      sample_list = [(s, s = 1, size(labels))]
      call shuffle(sample_list) 
      do s = 1, size(sample_list)
         
         !write(*,*) n, s
         select type(layer => network%model(2)%layer)
         type is (conv_mpnn_layer_type)
            call layer%set_graph(graphs(sample_list(s):sample_list(s)))
         end select
         call network%forward(reshape([1._real32], [1,1]))
         call network%backward(reshape([labels(sample_list(s))], [1,1]))
         !call network%model(size(network%model,1))%layer%get_output(output_tmp)
         !write(*,*) "predicted",output_tmp(1,1), labels(sample_list(s))
 
         call network%update()
  
      end do
   end do


  end subroutine network_train_graph

  function network_predict(x) result(y)
    implicit none
    real(real32), dimension(:,:), intent(in) :: x
    real(real32), dimension(1,size(x,2)) :: y

    y = network%predict(x)

  end function network_predict

  function network_predict_graph(graphs) result(y)
    implicit none
    type(graph_type), dimension(:), intent(in) :: graphs
    real(real32), dimension(size(graphs)) :: y

    integer :: s
    real(real32), dimension(:,:), allocatable :: output_tmp

    do s = 1, size(graphs)
       select type(layer => network%model(2)%layer)
       type is (conv_mpnn_layer_type)
          call layer%set_graph(graphs(s:s))
       end select
       call network%forward(reshape([1._real32], [1,1]))
       call network%model(size(network%model,1))%layer%get_output(output_tmp)
       y(s) = output_tmp(1,1)
   end do

  end function network_predict_graph


!###############################################################################
  function get_graph_from_basis(basis) result(graph)
    !! Get a graph representation of a basis.
    implicit none

    ! Arguments
    type(basis_type), intent(in) :: basis
    !! The basis to be converted to a graph.
    type(graph_type) :: graph
    !! The graph representation of the basis.

    ! Local variables
    integer :: is, ia, js, ja, i, j, k, iatom, jatom
    !! Loop indices.
    integer :: amax, bmax, cmax
    !! Maximum number of lattice translations.
    type(edge_type) :: edge
    !! An edge in the graph.
    real(real32) :: rtmp1
    !! Temporary real.
    real(real32) :: cutoff_min, cutoff_max
    !! Cutoff radii.
    real(real32), dimension(3) :: diff, vtmp1
    !! Difference vector and temporary vector.

    
    graph%num_vertices = basis%natom
    graph%num_vertex_features = 2
    graph%num_edge_features = 1

    allocate(graph%vertex(graph%num_vertices))

    iatom = 0
    do is = 1, basis%nspec
       do ia = 1, basis%spec(is)%num
          iatom = iatom + 1
          allocate(graph%vertex(iatom)%feature(graph%num_vertex_features))
          graph%vertex(iatom)%feature = [ basis%spec(is)%charge / 100._real32, &
               basis%spec(is)%mass / 52._real32 ]
       end do
    end do

    cutoff_min = 0.5_real32
    cutoff_max = 6.0_real32
    amax = ceiling(cutoff_max/modu(basis%lat(1,:)))
    bmax = ceiling(cutoff_max/modu(basis%lat(2,:)))
    cmax = ceiling(cutoff_max/modu(basis%lat(3,:)))

    iatom = 0
    allocate(graph%edge(0))
    spec_loop1: do is=1,basis%nspec
       atom_loop1: do ia=1,basis%spec(is)%num
          iatom = iatom + 1
          jatom = 0
          spec_loop2: do js=is,basis%nspec
             atom_loop2: do ja=1,basis%spec(js)%num
                jatom = jatom + 1
                if(is.eq.js.and.ja.lt.ia) cycle atom_loop2
                diff = basis%spec(is)%atom(ia,:3) -  basis%spec(js)%atom(ja,:3)
                diff = diff - ceiling(diff - 0.5_real32)
                do i=-amax,amax+1,1
                   vtmp1(1) = diff(1) + real(i, real32)
                   do j=-bmax,bmax+1,1
                      vtmp1(2) = diff(2) + real(j, real32)
                      do k=-cmax,cmax+1,1
                         vtmp1(3) = diff(3) + real(k, real32)
                         rtmp1 = modu(matmul(vtmp1,basis%lat))
                         if( rtmp1 .gt. cutoff_min .and. &
                             rtmp1 .lt. cutoff_max )then
                            edge%index = [iatom,jatom]
                            edge%feature = [rtmp1]
                            graph%edge = [ graph%edge, edge ]
                         end if
                      end do
                   end do
                end do
             end do atom_loop2
          end do spec_loop2
       end do atom_loop1
    end do spec_loop1
    graph%num_edges = size(graph%edge)
    call graph%generate_adjacency()
    call graph%calculate_degree()


  end function get_graph_from_basis
!###############################################################################

#endif
end module machine_learning