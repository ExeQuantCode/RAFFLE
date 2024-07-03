module generator
  use constants, only: real12
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
    integer, dimension(3) :: bins
    real(real12), dimension(3,3) :: lattice_host
    type(bas_type) :: basis_host
    type(gvector_container_type) :: distributions
    real(real12), dimension(3) :: method_probab
   contains
    procedure, pass(this) :: generate
    procedure, pass(this), private :: generate_structure
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
      real(real12), dimension(:), intent(in), optional :: method_probab
    end subroutine generate

    module function generate_structure( &
         this, &
         basis_initial, &
         placement_list, method_probab ) result(basis)
      class(raffle_generator_type), intent(in) :: this
      type(bas_type), intent(in) :: basis_initial
      integer, dimension(:,:), intent(in) :: placement_list
      real(real12), dimension(3) :: method_probab
      type(bas_type) :: basis
    end function generate_structure
  end interface



end module generator