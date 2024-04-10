module machine_learning
  use constants, only: real12
  use athena
  implicit none


  private

  public :: network_setup
  public :: network_train, network_predict

  type(network_type) :: network


contains

  subroutine network_setup(num_inputs, num_outputs)
    implicit none
    integer, intent(in) :: num_inputs, num_outputs

    call network%add(full_layer_type( &
         num_inputs = num_inputs, &
         num_outputs = 10, &
         activation_function = 'relu' &
         ))
    call network%add(full_layer_type( &
         num_outputs = num_outputs, &
         activation_function = 'sigmoid' &
         ))
    call network%compile( &
         optimiser = sgd_optimiser_type(learning_rate=1._real12), &
         loss_method = 'mse', &
         metrics = ['accuracy'] &
         )
    call network%set_batch_size(20)

  end subroutine network_setup

  subroutine network_train(x, y, num_epochs)
    implicit none
    real(real12), dimension(:,:), intent(in) :: x
    real(real12), dimension(:,:), intent(in) :: y
    integer, intent(in) :: num_epochs

    call network%train(x, y, num_epochs=num_epochs)

  end subroutine network_train

  function network_predict(x) result(y)
    implicit none
    real(real12), dimension(:,:), intent(in) :: x
    real(real12), dimension(size(x,1),1) :: y

    y = network%predict(x)

  end function network_predict

end module machine_learning