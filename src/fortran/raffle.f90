module raffle
  use raffle__constants, only: real12
  use generator, only: raffle_generator_type
  use raffle__distribs_container, only: distribs_container_type
  implicit none


  private
  public :: real12
  public :: distribs_container_type
  public :: raffle_generator_type


end module raffle