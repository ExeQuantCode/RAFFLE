module machine_learning
  use constants, only: real12
  use athena
  implicit none


  private

  public :: network_setup
  public :: network_train, network_train_graph
  public :: network_predict, network_predict_graph

  type(network_type) :: network


contains

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
    !      optimiser = sgd_optimiser_type(learning_rate=1.E-3_real12, momentum=0.99_real12,&
    !                                     clip_dict=clip_type(clip_min=0._real12, &
    !                                                         clip_max=1.E-4_real12), &
    !                                     regulariser=l1l2_regulariser_type(l1=1.E-3_real12, &
    !                                                                       l2=1.E-3_real12) &
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
    metric_dict%threshold = 1.E-1_real12
    call network%compile(optimiser=sgd_optimiser_type(), &
         loss_method="mse", metrics=metric_dict, &
         batch_size = 1, verbose = 0)
    call network%set_batch_size(1)

    network%metrics%threshold = 1.E-3_real12

  end subroutine network_setup

  subroutine network_train(x, y, num_epochs)
    implicit none
    real(real12), dimension(:,:), intent(in) :: x
    real(real12), dimension(:), intent(in) :: y
    integer, intent(in) :: num_epochs
    real(real12), dimension(:,:), allocatable :: y_renorm


    allocate(y_renorm(1,size(y,1)))
    y_renorm(1,:) = -1._real12 * y

    call renormalise_norm(y_renorm(1,:))
    write(*,*) y_renorm(1,:)
    call network%train(x, y_renorm, num_epochs=num_epochs, plateau_threshold=1.E-3_real12)
    !write(*,*) y_renorm


  end subroutine network_train

  subroutine network_train_graph(graphs, labels, num_epochs)
    implicit none
    type(graph_type), dimension(:), intent(in) :: graphs
    real(real12), dimension(:), intent(in) :: labels
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
         call network%forward(reshape([1._real12], [1,1]))
         call network%backward(reshape([labels(sample_list(s))], [1,1]))
         !call network%model(size(network%model,1))%layer%get_output(output_tmp)
         !write(*,*) "predicted",output_tmp(1,1), labels(sample_list(s))
 
         call network%update()
  
      end do
   end do


  end subroutine network_train_graph

  function network_predict(x) result(y)
    implicit none
    real(real12), dimension(:,:), intent(in) :: x
    real(real12), dimension(1,size(x,2)) :: y

    y = network%predict(x)

  end function network_predict

  function network_predict_graph(graphs) result(y)
    implicit none
    type(graph_type), dimension(:), intent(in) :: graphs
    real(real12), dimension(size(graphs)) :: y

    integer :: s
    real(real12), dimension(:,:), allocatable :: output_tmp

    do s = 1, size(graphs)
       select type(layer => network%model(2)%layer)
       type is (conv_mpnn_layer_type)
          call layer%set_graph(graphs(s:s))
       end select
       call network%forward(reshape([1._real12], [1,1]))
       call network%model(size(network%model,1))%layer%get_output(output_tmp)
       y(s) = output_tmp(1,1)
   end do

  end function network_predict_graph

end module machine_learning