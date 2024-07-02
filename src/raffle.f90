module raffle
  use constants, only: real12
  use gen, only: generation
  use rw_geom, only: bas_type
  use evolver, only: gvector_container_type
  implicit none



  private
  public :: raffle_generator_type


  type :: stoichiometry_type
    character(len=3) :: element
    integer :: num
  end type stoichiometry_type


  type :: raffle_generator_type
    real(real12), dimension(3,3) :: lattice_host
    type(bas_type) :: basis_host
    type(gvector_container_type) :: distributions
    real(real12), dimension(3) :: method_probab
   contains
    procedure :: generate
    !procedure :: get_structures
    !procedure :: evaluate
  end type raffle_generator_type

  interface raffle_generator_type
    module function init_raffle_generator( &
         lattice_host, basis_host, &
         width, sigma, cutoff_min, cutoff_max) result(generator)
      real(real12), dimension(3,3), intent(in) :: lattice_host
      type(bas_type), intent(in) :: basis_host
      real(real12), dimension(3), intent(in), optional :: width
      real(real12), dimension(3), intent(in), optional :: sigma
      real(real12), dimension(3), intent(in), optional :: cutoff_min
      real(real12), dimension(3), intent(in), optional :: cutoff_max
      type(raffle_generator_type) :: generator
    end function init_raffle_generator
  end interface raffle_generator_type

  interface
    module subroutine generate( this, &
         num_structures, stoichiometry, method_probab )
      class(raffle_generator_type), intent(inout) :: this
      integer, intent(in) :: num_structures
      type(stoichiometry_type), dimension(:), intent(in) :: stoichiometry
      real(real12), dimension(:), intent(in) :: method_probab
    end subroutine generate
    end interface

contains
  module function init_raffle_generator( &
       lattice_host, basis_host, width, sigma, cutoff_min, cutoff_max ) &
       result(generator)
    !! Initialise an instance of the raffle generator.
    !! Set up run-independent parameters.
    implicit none
    ! Arguments
    real(real12), dimension(3,3), intent(in) :: lattice_host
    !! Lattice vectors of the host structure.
    type(bas_type), intent(in) :: basis_host
    !! Basis of the host structure.
    real(real12), dimension(3), intent(in), optional :: width
    !! Width of the gaussians used in the 2-, 3-, and 4-body 
    !! distribution functions.
    real(real12), dimension(3), intent(in), optional :: sigma
    !! Width of the gaussians used in the 2-, 3-, and 4-body
    !! distribution functions.
    real(real12), dimension(3), intent(in), optional :: cutoff_min
    !! Minimum cutoff for the 2-, 3-, and 4-body distribution functions.
    real(real12), dimension(3), intent(in), optional :: cutoff_max
    !! Maximum cutoff for the 2-, 3-, and 4-body distribution functions.

    type(raffle_generator_type) :: generator


    generator%lattice_host = lattice_host
    generator%basis_host = basis_host

    if( present(width) ) &
         call generator%distributions%set_width(width)
    if( present(sigma) ) &
         call generator%distributions%set_sigma(sigma)
    if( present(cutoff_min) ) &
         call generator%distributions%set_cutoff_min(cutoff_min)
    if( present(cutoff_max) ) &
         call generator%distributions%set_cutoff_max(cutoff_max)


  end function init_raffle_generator


  module subroutine generate( this, &
       num_structures, stoichiometry, method_probab )
    !! Generate random structures.
    implicit none
    ! Arguments
    class(raffle_generator_type), intent(inout) :: this
    !! Instance of the raffle generator.
    integer, intent(in) :: num_structures
    !! Number of structures to generate.
    type(stoichiometry_type), dimension(:), intent(in) :: stoichiometry
    !! Stoichiometry of the structures to generate.
    real(real12), dimension(:), intent(in) :: method_probab
    !! Probability of each placement method.

    call generation( this%distributions, num_structures, &
        stoichiometry(:)%element, stoichiometry(:)%num, &
        method_probab )

  end subroutine generate

end module raffle