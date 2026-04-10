module raffle
  use raffle__constants, only: real32
  use raffle__io_utils, only: raffle__version__
  use raffle__generator, only: raffle_generator_type
  use raffle__distribs_container, only: distribs_container_type
#ifdef ENABLE_ATHENA
  use raffle__nn_fingerprint, only: nn_fingerprint_type
  use raffle__gnn_fingerprint, only: gnn_fingerprint_type
#endif
  use raffle__cache, only: &
       store_probability_density, retrieve_probability_density
  implicit none


  private
  public :: real32
  public :: distribs_container_type
  public :: raffle_generator_type
#ifdef ENABLE_ATHENA
  public :: nn_fingerprint_type
  public :: gnn_fingerprint_type
#endif


end module raffle
