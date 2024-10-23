module raffle
  use raffle__constants, only: real32
  use raffle__io_utils, only: raffle__version__
  use raffle__generator, only: raffle_generator_type
  use raffle__distribs_container, only: distribs_container_type
  implicit none


  private
  public :: real32
  public :: distribs_container_type
  public :: raffle_generator_type


end module raffle