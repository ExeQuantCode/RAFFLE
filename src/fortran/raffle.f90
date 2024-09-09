module raffle
  use constants, only: real12
  use generator, only: raffle_generator_type
  use evolver, only: gvector_container_type
  implicit none


  private
  public :: real12
  public :: gvector_container_type
  public :: raffle_generator_type


end module raffle