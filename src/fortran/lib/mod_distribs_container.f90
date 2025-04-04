module raffle__distribs_container
  !! Module for handling the distribution function container.
  !!
  !! This module defines the distribution function container and associated
  !! procedures.
  !! The container holds the distribution functions for a set of atomic
  !! structures, alongside parameters for initialising the distributions.
  !! The container also holds the generalised distribution functions, built
  !! from the distribution functions of the individual systems.
  !! The generalised distribution functions are used to evaluate the viability
  !! of a new structure.
  use raffle__constants, only: real32, pi
  use raffle__io_utils, only: stop_program, print_warning, suppress_warnings
  use raffle__misc, only: set, icount, strip_null, sort_str
  use raffle__misc_maths, only: triangular_number, set_difference
  use raffle__geom_rw, only: basis_type, get_element_properties
  use raffle__element_utils, only: &
       element_type, element_bond_type, &
       element_database, element_bond_database
  use raffle__distribs, only: distribs_base_type, distribs_type, get_distrib
  use raffle__distribs_host, only: distribs_host_type
  implicit none


  private

  public :: distribs_container_type


  type :: distribs_container_type
     !! Container for distribution functions.
     !!
     !! This type contains the distribution functions for a set of atomic
     !! structures, alongside parameters for initialising the distributions.
     integer :: num_evaluated = 0
     !! Number of evaluated systems.
     integer :: num_evaluated_allocated = 0
     !! Number of evaluated systems still allocated.
     real(real32) :: kBT = 0.2_real32
     !! Boltzmann constant times temperature.
     logical :: weight_by_hull = .false.
     !! Boolean whether to weight the distribution functions by the energy
     !! above the hull. If false, the formation energy from the element
     !! reference energies is used.
     real(real32) :: &
          viability_3body_default = 0.1_real32, &
          viability_4body_default = 0.1_real32
     !! Default viability for the 3- and 4-body distribution functions.
     logical, dimension(:), allocatable :: &
          in_dataset_2body, in_dataset_3body, in_dataset_4body
     !! Whether the 2-, 3-, and 4-body distribution functions are in
     !! the dataset.
     real(real32), dimension(:), allocatable :: &
          best_energy_pair, &
          best_energy_per_species
     !! Best energy for the 2-body and species distribution functions.
     integer, dimension(3) :: nbins = -1
     !! Number of bins for the 2-, 3-, and 4-body distribution functions.
     real(real32), dimension(3) :: &
          sigma = [ 0.1_real32, 0.1_real32, 0.1_real32 ]
     !! Sigma of the gaussians used in the 2-, 3-, and 4-body
     !! distribution functions.
     real(real32), dimension(3) :: &
          width = [ 0.025_real32, pi/64._real32, pi/64._real32 ]
     !! Width of the bins used in the 2-, 3-, and 4-body distribution functions.
     real(real32), dimension(3) :: &
          cutoff_min = [ 0.5_real32, 0._real32, 0._real32 ]
     !! Minimum cutoff for the 2-, 3-, and 4-body distribution functions.
     real(real32), dimension(3) :: &
          cutoff_max = [ 6._real32, pi, pi ]
     !! Maximum cutoff for the 2-, 3-, and 4-body distribution functions.
     real(real32), dimension(4) :: &
          radius_distance_tol = [ 1.5_real32, 2.5_real32, 3._real32, 6._real32 ]
     !! Tolerance for the distance between atoms for 3- and 4-body.
     !! index 1 = lower bound for 3-body
     !! index 2 = upper bound for 3-body
     !! index 3 = lower bound for 4-body
     !! index 4 = upper bound for 4-body
     real(real32), dimension(:), allocatable :: &
          norm_2body, norm_3body, norm_4body
     !! Normalisation factors for the 2-, 3-, and 4-body distribution functions.
     type(distribs_base_type) :: gdf
     !! Generalised distribution functions for all systems.
     !! Generated from combining the energy-weighted distribution functions
     !! of all systems
     integer, dimension(:), allocatable :: element_map
     !! Mapping of host elements to distribution function elements.
     type(distribs_host_type) :: host_system
     !! Host structure for the distribution functions.
     type(distribs_type), dimension(:), allocatable :: system
     !! Distribution functions for each system.
     type(element_type), dimension(:), allocatable :: element_info
     !! Information about the elements in the container.
     type(element_bond_type), dimension(:), allocatable :: bond_info
     !! Information about the 2-body bonds in the container.
   contains
     procedure, pass(this) :: set_width
     !! Set the width of the bins used in the 2-, 3-, and 4-body.
     procedure, pass(this) :: set_sigma
     !! Set the sigma of the gaussians used in the 2-, 3-, and 4-body.
     procedure, pass(this) :: set_cutoff_min
     !! Set the minimum cutoff for the 2-, 3-, and 4-body.
     procedure, pass(this) :: set_cutoff_max
     !! Set the maximum cutoff for the 2-, 3-, and 4-body.
     procedure, pass(this) :: set_radius_distance_tol
     !! Set the tolerance for the distance between atoms for 3- and 4-body.

     procedure, pass(this) :: create
     !! Create the distribution functions for all systems, and the learned one.
     procedure, pass(this) :: update
     !! Update the distribution functions for all systems, and the learned one.

     procedure, pass(this) :: deallocate_systems
     !! Deallocate the systems in the container.

     procedure, pass(this) :: add, add_basis
     !! Add a system to the container.

     procedure, pass(this), private :: set_element_info
     !! Set the list of elements for the container.
     procedure, pass(this), private :: update_element_info
     !! Update the element information in the container.
     procedure, pass(this) :: set_element_energy
     !! Set the energy of an element in the container.
     procedure, pass(this) :: set_element_energies
     !! Set the energies of elements in the container.
     procedure, pass(this) :: get_element_energies
     !! Return the energies of elements in the container.
     procedure, pass(this) :: get_element_energies_staticmem
     !! Return the energies of elements in the container.
     !! Used in Python interface.


     procedure, pass(this) :: set_element_map
     !! Set the mapping of elements to distribution function elements.
     procedure, pass(this), private :: set_bond_info
     !! Set the 2-body bond information for the container.
     procedure, pass(this), private :: update_bond_info
     !! Update the bond information in the container.
     procedure, pass(this) :: set_bond_radius
     !! Set the radius of a bond in the container.
     procedure, pass(this) :: set_bond_radii
     !! Set the radii of multiple bonds in the container.
     procedure, pass(this) :: get_bond_radii
     !! Return the radii of all bonds in the container.
     procedure, pass(this) :: get_bond_radii_staticmem
     !! Return the radii of all bonds in the container.
     !! Used in Python interface.


     procedure, pass(this) :: set_best_energy
     !! Set the best energy and system in the container.
     procedure, pass(this) :: initialise_gdfs
     !! Initialise the distribution functions in the container.
     procedure, pass(this) :: set_gdfs_to_default
     !! Set the generalised distribution function to the default value.
     procedure, pass(this) :: evolve
     !! Evolve the learned distribution function.


     procedure, pass(this) :: write_gdfs
     !! Write the generalised distribution functions to a file.
     procedure, pass(this) :: read_gdfs
     !! Read the generalised distribution functions from a file.
     procedure, pass(this) :: write_dfs
     !! Write all distribution functions to a file.
     procedure, pass(this) :: read_dfs
     !! Read all distribution functions from a file.
     procedure, pass(this) :: write_2body
     !! Write the learned 2-body distribution function to a file.
     procedure, pass(this) :: write_3body
     !! Write the learned 3-body distribution function to a file.
     procedure, pass(this) :: write_4body
     !! Write the learned 4-body distribution function to a file.
     procedure, pass(this) :: get_pair_index
     !! Return the index for bond_info given two elements.
     procedure, pass(this) :: get_element_index
     !! Return the index for element_info given one element.
     procedure, pass(this) :: get_bin
     !! Return the bin index for a given distance.
  end type distribs_container_type

  interface distribs_container_type
     !! Interface for the distribution functions container.
     module function init_distribs_container( &
          nbins, width, sigma, cutoff_min, cutoff_max &
     ) result(distribs_container)
       !! Initialise the distribution functions container.
       integer, dimension(3), intent(in), optional :: nbins
       !! Optional. Number of bins for the 2-, 3-, and 4-body distribution
       !! functions.
       real(real32), dimension(3), intent(in), optional :: width, sigma
       !! Optional. Width and sigma of the gaussians used in the 2-, 3-, and
       !! 4-body.
       real(real32), dimension(3), intent(in), optional :: &
            cutoff_min, cutoff_max
       !! Optional. Minimum and maximum cutoff for the 2-, 3-, and 4-body.
       type(distribs_container_type) :: distribs_container
       !! Instance of the distribution functions container.
     end function init_distribs_container
  end interface distribs_container_type


contains

!###############################################################################
  module function init_distribs_container( &
       nbins, width, sigma, &
       cutoff_min, cutoff_max &
  ) result(distribs_container)
    !! Initialise the distribution functions container.
    implicit none

    ! Arguments
    integer, dimension(3), intent(in), optional :: nbins
    !! Optional. Number of bins for the 2-, 3-, and 4-body distribution
    !! functions.
    real(real32), dimension(3), intent(in), optional :: width, sigma
    !! Optional. Width and sigma of the gaussians used in the 2-, 3-, and
    !! 4-body.
    real(real32), dimension(3), intent(in), optional :: cutoff_min, cutoff_max
    !! Optional. Minimum and maximum cutoff for the 2-, 3-, and 4-body.
    type(distribs_container_type) :: distribs_container
    !! Instance of the distribution functions container.

    ! Local variables
    character(256) :: stop_msg
    !! Error message.


    if(present(nbins))then
       if(all(nbins .gt. 0)) distribs_container%nbins = nbins
    end if

    if(present(width))then
       if(all(width.ge.0._real32)) distribs_container%width = width
    end if

    if(present(sigma))then
       if(all(sigma.ge.0._real32)) distribs_container%sigma = sigma
    end if

    if(present(cutoff_min))then
       if(any(cutoff_min.ge.0._real32)) &
            distribs_container%cutoff_min = cutoff_min
    end if
    if(present(cutoff_max))then
       if(all(cutoff_max.ge.0._real32)) &
            distribs_container%cutoff_max = cutoff_max
    end if
    if( &
         any(distribs_container%cutoff_max .le. distribs_container%cutoff_min) &
    )then
       write(stop_msg,*) &
            "cutoff_max <= cutoff_min" // &
            achar(13) // achar(10) // &
            "cutoff min: ", distribs_container%cutoff_min, &
            achar(13) // achar(10) // &
            "cutoff max: ", distribs_container%cutoff_max
       call stop_program( stop_msg )
       return
    end if

  end function init_distribs_container
!###############################################################################


!###############################################################################
  subroutine set_width(this, width)
    !! Set the width of the gaussians used in the 2-, 3-, and 4-body
    !! distribution functions.
    implicit none

    ! Arguments
    class(distribs_container_type), intent(inout) :: this
    !! Parent. Instance of distribution functions container.
    real(real32), dimension(3), intent(in) :: width
    !! Width of the gaussians used in the 2-, 3-, and 4-body
    !! distribution functions.

    this%width = width

  end subroutine set_width
!###############################################################################


!###############################################################################
  subroutine set_sigma(this, sigma)
    !! Set the sigma of the gaussians used in the 2-, 3-, and 4-body
    !! distribution functions.
    implicit none

    ! Arguments
    class(distribs_container_type), intent(inout) :: this
    !! Parent. Instance of distribution functions container.
    real(real32), dimension(3), intent(in) :: sigma
    !! Sigma of the gaussians used in the 2-, 3-, and 4-body distribution
    !! functions.

    this%sigma = sigma

  end subroutine set_sigma
!###############################################################################


!###############################################################################
  subroutine set_cutoff_min(this, cutoff_min)
    !! Set the minimum cutoff for the 2-, 3-, and 4-body distribution functions.
    implicit none

    ! Arguments
    class(distribs_container_type), intent(inout) :: this
    !! Parent. Instance of distribution functions container.
    real(real32), dimension(3), intent(in) :: cutoff_min
    !! Minimum cutoff for the 2-, 3-, and 4-body distribution functions.

    this%cutoff_min = cutoff_min

  end subroutine set_cutoff_min
!###############################################################################


!###############################################################################
  subroutine set_cutoff_max(this, cutoff_max)
    !! Set the maximum cutoff for the 2-, 3-, and 4-body distribution functions.
    implicit none

    ! Arguments
    class(distribs_container_type), intent(inout) :: this
    !! Parent. Instance of distribution functions container.
    real(real32), dimension(3), intent(in) :: cutoff_max
    !! Maximum cutoff for the 2-, 3-, and 4-body distribution functions.

    this%cutoff_max = cutoff_max

  end subroutine set_cutoff_max
!###############################################################################


!###############################################################################
  subroutine set_radius_distance_tol(this, radius_distance_tol)
    !! Set the tolerance for the distance between atoms for 3- and 4-body.
    implicit none

    ! Arguments
    class(distribs_container_type), intent(inout) :: this
    !! Parent. Instance of distribution functions container.
    real(real32), dimension(4), intent(in) :: radius_distance_tol
    !! Tolerance for the distance between atoms for 3- and 4-body.

    this%radius_distance_tol = radius_distance_tol

  end subroutine set_radius_distance_tol
!###############################################################################


!###############################################################################
  subroutine create( &
       this, basis_list, energy_above_hull_list, deallocate_systems, &
       verbose &
  )
    !! create the distribution functions from the input file
    implicit none
    ! Arguments
    class(distribs_container_type), intent(inout) :: this
    !! Parent. Instance of distribution functions container.
    type(basis_type), dimension(:), intent(in) :: basis_list
    !! List of basis structures.
    real(real32), dimension(:), intent(in), optional :: energy_above_hull_list
    !! List of energies above the hull for the structures.
    logical, intent(in), optional :: deallocate_systems
    !! Boolean whether to deallocate the systems after the
    !! distribution functions are created.
    integer, intent(in), optional :: verbose
    !! Verbosity level.

    ! Local variables
    logical :: deallocate_systems_
    !! Boolean whether to deallocate the systems after the distribution
    character(256) :: stop_msg
    !! Error message.
    integer :: verbose_
    !! Verbosity level.
    logical :: suppress_warnings_store
    !! Boolean to store the suppress_warnings value.


    ! Set the verbosity level
    verbose_ = 0
    if(present(verbose)) verbose_ = verbose
    if(verbose_ .eq. 0)then
       suppress_warnings_store = suppress_warnings
       suppress_warnings = .true.
    end if

    ! Check if element_database is allocated
    if(.not.allocated(element_database))then
       write(stop_msg,*) "element_database not allocated" // &
            achar(13) // achar(10) // &
            "Run the set_element_energies() procedure of " // &
            "distribs_container_type before calling create()"
       call stop_program( stop_msg )
       return
    end if

    ! Check if energy_above_hull_list and basis_list are the same size
    if(present(energy_above_hull_list))then
       if(size(energy_above_hull_list).eq.0)then
          this%weight_by_hull = .false.
       elseif(size(energy_above_hull_list) .ne. size(basis_list))then
          this%weight_by_hull = .true.
          write(stop_msg,*) "energy_above_hull_list and basis_list " // &
               "not the same size" // &
               achar(13) // achar(10) // &
               "energy_above_hull_list: ", size(energy_above_hull_list), &
               achar(13) // achar(10) // &
               "basis_list: ", size(basis_list)
          call stop_program( stop_msg )
          return
       end if
    end if
    if(this%weight_by_hull.and..not.present(energy_above_hull_list))then
       write(stop_msg,*) "energy_above_hull_list not present" // &
            achar(13) // achar(10) // &
            "energy_above_hull_list must be present when using hull weighting"
       call stop_program( stop_msg )
       return
    end if


    ! Check if deallocate_systems is present
    deallocate_systems_ = .true.
    if(present(deallocate_systems)) deallocate_systems_ = deallocate_systems

    this%num_evaluated = 0
    this%num_evaluated_allocated = 0
    if(allocated(this%gdf%df_2body)) deallocate(this%gdf%df_2body)
    if(allocated(this%gdf%df_3body)) deallocate(this%gdf%df_3body)
    if(allocated(this%gdf%df_4body)) deallocate(this%gdf%df_4body)
    if(allocated(this%norm_2body)) deallocate(this%norm_2body)
    if(allocated(this%norm_3body)) deallocate(this%norm_3body)
    if(allocated(this%norm_4body)) deallocate(this%norm_4body)
    if(allocated(this%system)) deallocate(this%system)
    if(allocated(this%in_dataset_2body)) deallocate(this%in_dataset_2body)
    if(allocated(this%in_dataset_3body)) deallocate(this%in_dataset_3body)
    if(allocated(this%in_dataset_4body)) deallocate(this%in_dataset_4body)
    if(allocated(this%best_energy_pair)) deallocate(this%best_energy_pair)
    if(allocated(this%best_energy_per_species)) &
         deallocate(this%best_energy_per_species)
    allocate(this%system(0))
    call this%add(basis_list)
    if(present(energy_above_hull_list).and.this%weight_by_hull)then
       this%system(:)%energy_above_hull = energy_above_hull_list(:)
    end if
    call this%set_bond_info()
    call this%evolve()
    if(deallocate_systems_) call this%deallocate_systems()
    if(this%host_system%defined) &
         call this%host_system%set_element_map(this%element_info)

    if(verbose_ .eq. 0)then
       suppress_warnings = suppress_warnings_store
    end if

  end subroutine create
!###############################################################################


!###############################################################################
  subroutine update( &
       this, basis_list, energy_above_hull_list, from_host, &
       deallocate_systems, &
       verbose &
  )
    !! update the distribution functions from the input file
    implicit none
    ! Arguments
    class(distribs_container_type), intent(inout) :: this
    !! Parent. Instance of distribution functions container.
    type(basis_type), dimension(:), intent(in) :: basis_list
    !! List of basis structures.
    real(real32), dimension(:), intent(in), optional :: energy_above_hull_list
    !! List of energies above the hull for the structures.
    logical, intent(in), optional :: from_host
    !! Optional. Boolean whether structures are derived from the host.
    logical, intent(in), optional :: deallocate_systems
    !! Boolean whether to deallocate the systems after the
    !! distribution functions are created.
    integer, intent(in), optional :: verbose
    !! Verbosity level.

    ! Local variables
    integer :: i
    !! Loop index.
    logical :: deallocate_systems_
    !! Boolean whether to deallocate the systems after the distribution
    logical :: from_host_
    !! Boolean whether structures are derived from the host.
    character(256) :: stop_msg
    !! Error message.
    integer :: verbose_
    !! Verbosity level.
    logical :: suppress_warnings_store
    !! Boolean to store the suppress_warnings value.


    ! Set the verbosity level
    verbose_ = 0
    if(present(verbose)) verbose_ = verbose
    if(verbose_ .eq. 0)then
       suppress_warnings_store = suppress_warnings
       suppress_warnings = .true.
    end if

    ! Check if energy_above_hull_list and basis_list are the same size
    if(present(energy_above_hull_list))then
       if(size(energy_above_hull_list).eq.0 .and. .not. this%weight_by_hull)then
       elseif(size(energy_above_hull_list) .ne. size(basis_list) .and. &
            this%weight_by_hull &
       )then
          write(stop_msg,*) "energy_above_hull_list and basis_list " // &
               "not the same size whilst using hull weighting" // &
               achar(13) // achar(10) // &
               "energy_above_hull_list: ", size(energy_above_hull_list), &
               achar(13) // achar(10) // &
               "basis_list: ", size(basis_list)
          call stop_program( stop_msg )
          return
       end if
    end if
    if(this%weight_by_hull.and..not.present(energy_above_hull_list))then
       write(stop_msg,*) "energy_above_hull_list not present" // &
            achar(13) // achar(10) // &
            "energy_above_hull_list must be present when using hull weighting"
       call stop_program( stop_msg )
       return
    end if

    ! Check if from_host is present
    if(present(from_host))then
       from_host_ = from_host
       if(this%weight_by_hull.and.from_host_)then
          write(stop_msg,*) "Hull weighting and from_host are incompatible" // &
               achar(13) // achar(10) // &
               "Set from_host = .false. to use hull weighting"
          call stop_program( stop_msg )
          return
       end if
    else
       if(this%weight_by_hull)then
          from_host_ = .false.
       else
          from_host_ = .true.
       end if
    end if

    ! Check if deallocate_systems is present
    deallocate_systems_ = .true.
    if(present(deallocate_systems)) deallocate_systems_ = deallocate_systems

    ! Add the new basis structures
    call this%add(basis_list)
    call this%update_bond_info()
    if(present(energy_above_hull_list).and.this%weight_by_hull)then
       this%system(this%num_evaluated_allocated + 1:)%energy_above_hull = &
            energy_above_hull_list(:)
    end if

    ! If the structures are derived from the host, subtract the interface energy
    if(from_host_)then
       if(.not.this%host_system%defined)then
          write(stop_msg,*) "host not set" // &
               achar(13) // achar(10) // &
               "Run the set_host() procedure of parent of" // &
               "distribs_container_type before calling create()"
          call stop_program( stop_msg )
          return
       else
          if(.not.allocated(this%host_system%df_2body))then
             call this%host_system%calculate( &
                  this%host_system%basis, &
                  width = this%width, &
                  sigma = this%sigma, &
                  cutoff_min = this%cutoff_min, &
                  cutoff_max = this%cutoff_max, &
                  radius_distance_tol = this%radius_distance_tol &
             )
          end if
          call this%host_system%calculate_interface_energy(this%element_info)
          do i = this%num_evaluated_allocated + 1, size(this%system), 1
             this%system(i)%from_host = .true.
             this%system(i)%energy = this%system(i)%energy - &
                  this%host_system%interface_energy
             this%system(i)%num_atoms = this%system(i)%num_atoms - &
                  this%host_system%num_atoms
          end do
       end if
    end if

    call this%evolve()
    if(deallocate_systems_) call this%deallocate_systems()
    if(this%host_system%defined) &
         call this%host_system%set_element_map(this%element_info)

    if(verbose_ .eq. 0)then
       suppress_warnings = suppress_warnings_store
    end if

  end subroutine update
!###############################################################################


!###############################################################################
  subroutine deallocate_systems(this)
    !! Deallocate the systems in the container.
    implicit none

    ! Arguments
    class(distribs_container_type), intent(inout) :: this
    !! Parent. Instance of distribution functions container.

    deallocate(this%system)
    this%num_evaluated_allocated = 0

  end subroutine deallocate_systems
!###############################################################################


!###############################################################################
  subroutine write_gdfs(this, file)
    !! Write the generalised distribution functions to a file.
    implicit none

    ! Arguments
    class(distribs_container_type), intent(in) :: this
    !! Parent. Instance of distribution functions container.
    character(*), intent(in) :: file
    !! Filename to write the generalised distribution functions to.

    ! Local variables
    integer :: unit
    !! File unit.
    integer :: i
    !! Loop index.
    character(256) :: fmt
    !! Format string.
    character(256) :: stop_msg
    !! Error message.

    if(.not.allocated(this%gdf%df_2body))then
       write(stop_msg,*) &
            "Generalised distribution functions are not allocated." // &
            achar(13) // achar(10) // &
            "create() or read() must be called before writing the " // &
            "generalised distribution functions."
       call stop_program( stop_msg )
       return
    end if
    open(newunit=unit, file=file)
    write(unit, '("# nbins",3(1X,I0))') this%nbins
    write(unit, '("# width",3(1X,ES0.4))') this%width
    write(unit, '("# sigma",3(1X,ES0.4))') this%sigma
    write(unit, '("# cutoff_min",3(1X,ES0.4))') this%cutoff_min
    write(unit, '("# cutoff_max",3(1X,ES0.4))') this%cutoff_max
    write(unit, '("# radius_distance_tol",4(1X,ES0.4))') &
         this%radius_distance_tol
    write(fmt, '("(""# "",A,",I0,"(1X,A))")') size(this%element_info)
    write(unit, fmt) "elements", this%element_info(:)%name
    write(fmt, '("(""# "",A,",I0,"(1X,ES0.4))")') size(this%element_info)
    write(unit, fmt) "energies", this%element_info(:)%energy
    write(unit, fmt) "best_energy_per_element", this%best_energy_per_species
    write(unit, fmt) "3-body_norm", this%norm_3body
    write(unit, fmt) "4-body_norm", this%norm_4body
    write(fmt, '("(""# "",A,",I0,"(1X,L1))")') size(this%element_info)
    write(unit, fmt) "in_dataset_3body", this%in_dataset_3body
    write(unit, fmt) "in_dataset_4body", this%in_dataset_4body
    write(fmt, '("(""# "",A,",I0,"(1X,A))")') size(this%bond_info)
    write(unit, fmt) "element_pairs", &
         ( &
              trim(this%bond_info(i)%element(1)) // "-" // &
              trim(this%bond_info(i)%element(2)), &
              i = 1, size(this%bond_info) &
         )
    write(fmt, '("(""# "",A,",I0,"(1X,ES0.4))")') size(this%bond_info)
    write(unit, fmt) "radii", this%bond_info(:)%radius_covalent
    write(unit, fmt) "best_energy_per_pair", this%best_energy_pair
    write(unit, fmt) "2-body_norm", this%norm_2body
    write(fmt, '("(""# "",A,",I0,"(1X,L1))")') size(this%bond_info)
    write(unit, fmt) "in_dataset_2body", this%in_dataset_2body
    write(unit, *)
    write(unit, '("# 2-body")')
    write(fmt,'("(""# bond-length "",",I0,"(1X,A))")') size(this%bond_info)
    write(unit, fmt) &
         ( &
              trim(this%bond_info(i)%element(1)) // "-" // &
              trim(this%bond_info(i)%element(2)), &
              i = 1, size(this%bond_info) &
         )
    do i = 1, this%nbins(1)
       write(unit, *) &
            this%cutoff_min(1) + this%width(1) * ( i - 1 ), &
            this%gdf%df_2body(i,:)
    end do
    write(unit, *)
    write(unit, '("# 3-body")')
    write(fmt,'("(""# bond-angle "",",I0,"(1X,A))")') size(this%bond_info)
    write(unit, fmt) this%element_info(:)%name
    do i = 1, this%nbins(2)
       write(unit, *) &
            this%cutoff_min(2) + this%width(2) * ( i - 1 ), &
            this%gdf%df_3body(i,:)
    end do
    write(unit, *)
    write(unit, '("# 4-body")')
    write(fmt,'("(""# dihedral-angle "",",I0,"(1X,A))")') size(this%bond_info)
    write(unit, fmt) this%element_info(:)%name
    do i = 1, this%nbins(3)
       write(unit, *) &
            this%cutoff_min(2) + this%width(2) * ( i - 1 ), &
            this%gdf%df_4body(i,:)
    end do
    close(unit)

  end subroutine write_gdfs
!###############################################################################


!###############################################################################
  subroutine read_gdfs(this, file)
    !! Read the generalised distribution functions from a file.
    implicit none

    ! Arguments
    class(distribs_container_type), intent(inout) :: this
    !! Parent. Instance of distribution functions container.
    character(*), intent(in) :: file
    !! Filename to read the generalised distribution functions from.

    ! Local variables
    integer :: unit
    !! File unit.

    integer :: i
    !! Loop index.
    integer :: nspec
    !! Number of species.
    logical :: exist
    !! Boolean whether the file exists.
    character(256) :: buffer, buffer1, buffer2
    !! Buffer for reading lines.

    ! check if file exists
    inquire(file=file, exist=exist)
    if(.not.exist)then
       call stop_program( "File does not exist" )
       return
    end if

    ! read the file
    open(newunit=unit, file=file)
    read(unit, *) buffer1, buffer2, this%nbins
    read(unit, *) buffer1, buffer2, this%width
    read(unit, *) buffer1, buffer2, this%sigma
    read(unit, *) buffer1, buffer2, this%cutoff_min
    read(unit, *) buffer1, buffer2, this%cutoff_max
    read(unit, *) buffer1, buffer2, this%radius_distance_tol
    read(unit, '(A)') buffer
    nspec = icount(buffer(index(buffer,"elements")+8:))
    allocate(this%element_info(nspec))
    read(buffer, *) buffer1, buffer2, this%element_info(:)%name
    read(unit, *) buffer1, buffer2, this%element_info(:)%energy
    do i = 1, nspec
       call this%set_element_energy( &
            this%element_info(i)%name, &
            this%element_info(i)%energy &
       )
       call this%element_info(i)%set(this%element_info(i)%name)
    end do
    call this%update_bond_info()
    allocate(this%best_energy_per_species(nspec))
    allocate(this%norm_3body(nspec))
    allocate(this%norm_4body(nspec))
    allocate(this%in_dataset_3body(nspec))
    allocate(this%in_dataset_4body(nspec))
    read(unit, *) buffer1, buffer2, this%best_energy_per_species
    read(unit, *) buffer1, buffer2, this%norm_3body
    read(unit, *) buffer1, buffer2, this%norm_4body
    read(unit, *) buffer1, buffer2, this%in_dataset_3body
    read(unit, *) buffer1, buffer2, this%in_dataset_4body
    read(unit, *)
    allocate(this%best_energy_pair(size(this%bond_info)))
    allocate(this%norm_2body(size(this%bond_info)))
    allocate(this%in_dataset_2body(size(this%bond_info)))
    read(unit, *) buffer1, buffer2, this%bond_info(:)%radius_covalent
    read(unit, *) buffer1, buffer2, this%best_energy_pair
    read(unit, *) buffer1, buffer2, this%norm_2body
    read(unit, *) buffer1, buffer2, this%in_dataset_2body
    read(unit, *)
    read(unit, *)
    read(unit, *)
    allocate(this%gdf%df_2body(this%nbins(1),size(this%bond_info)))
    do i = 1, this%nbins(1)
       read(unit, *) buffer, this%gdf%df_2body(i,:)
    end do
    read(unit, *)
    read(unit, *)
    read(unit, *)
    allocate(this%gdf%df_3body(this%nbins(2),nspec))
    do i = 1, this%nbins(2)
       read(unit, *) buffer, this%gdf%df_3body(i,:)
    end do
    read(unit, *)
    read(unit, *)
    read(unit, *)
    allocate(this%gdf%df_4body(this%nbins(3),nspec))
    do i = 1, this%nbins(3)
       read(unit, *) buffer, this%gdf%df_4body(i,:)
    end do
    close(unit)

  end subroutine read_gdfs
!###############################################################################


!###############################################################################
  subroutine write_dfs(this, file)
    !! Write all distribution functions for each system to a file.
    implicit none

    ! Arguments
    class(distribs_container_type), intent(in) :: this
    !! Parent. Instance of distribution functions container.
    character(*), intent(in) :: file
    !! Filename to write the distribution functions to.

    ! Local variables
    integer :: unit
    !! File unit.
    integer :: i, j
    !! Loop indices.
    character(256) :: stop_msg
    !! Error message.

    if(.not.allocated(this%system))then
       write(stop_msg,*) "No systems to write" // &
            achar(13) // achar(10) // &
            "Systems either not created or deallocated after evolve" // &
            achar(13) // achar(10) // &
            "To stop automatic deallocation, " // &
            "use the following flag in create()" // &
            achar(13) // achar(10) // &
            "   deallocate_systems = .false."
       call stop_program( stop_msg )
       return
    end if
    open(newunit=unit, file=file)
    write(unit, *) "nbins", this%nbins
    write(unit, *) "width", this%width
    write(unit, *) "sigma", this%sigma
    write(unit, *) "cutoff_min", this%cutoff_min
    write(unit, *) "cutoff_max", this%cutoff_max
    write(unit, *)
    do i = 1, size(this%system,1)
       write(unit, *) this%system(i)%energy
       write(unit, *) this%system(i)%element_symbols
       write(unit, *) this%system(i)%stoichiometry
       do j = 1, this%nbins(1)
          write(unit, *) this%system(i)%df_2body(j,:)
       end do
       do j = 1, this%nbins(2)
          write(unit, *) this%system(i)%df_3body(j,:)
       end do
       do j = 1, this%nbins(3)
          write(unit, *) this%system(i)%df_4body(j,:)
       end do
       write(unit, *)
    end do
    close(unit)

  end subroutine write_dfs
!###############################################################################


!###############################################################################
  subroutine read_dfs(this, file)
    !! Read all distribution functions for each system from a file.
    implicit none

    ! Arguments
    class(distribs_container_type), intent(inout) :: this
    !! Parent. Instance of distribution functions container.
    character(*), intent(in) :: file
    !! Filename to read the distribution functions from.

    ! Local variables
    integer :: unit
    !! File unit.
    integer :: j
    !! Loop indices.
    integer :: iostat
    !! I/O status.
    integer :: num_species, num_pairs
    !! Number of species and pairs.
    character(256) :: buffer
    !! Buffer for reading lines.
    type(distribs_type) :: system
    !! System to read distribution functions into.


    open(newunit=unit, file=file)
    read(unit, *) buffer, this%nbins
    read(unit, *) buffer, this%width
    read(unit, *) buffer, this%sigma
    read(unit, *) buffer, this%cutoff_min
    read(unit, *) buffer, this%cutoff_max
    do
       read(unit, '(A)', iostat=iostat) buffer
       if(iostat.ne.0) exit
       if(trim(buffer).eq.''.or.trim(buffer).eq.'#') cycle
       read(buffer, *) system%energy
       read(unit, '(A)') buffer
       num_species = icount(buffer)
       allocate(system%element_symbols(num_species))
       allocate(system%stoichiometry(num_species))
       read(buffer, *) system%element_symbols
       read(unit, *) system%stoichiometry
       system%num_atoms = sum(system%stoichiometry)
       num_pairs = nint( &
            gamma(real(num_species + 2, real32)) / &
            ( gamma(real(num_species, real32)) * gamma( 3._real32 ) ) &
       )
       allocate(system%df_2body(this%nbins(1),num_pairs))
       do j = 1, this%nbins(1)
          read(unit, *) system%df_2body(j,:)
       end do
       allocate(system%df_3body(this%nbins(2),num_species))
       do j = 1, this%nbins(2)
          read(unit, *) system%df_3body(j,:)
       end do
       allocate(system%df_4body(this%nbins(3),num_species))
       do j = 1, this%nbins(3)
          read(unit, *) system%df_4body(j,:)
       end do

       this%system = [ this%system, system ]
       deallocate(&
            system%element_symbols, system%stoichiometry, &
            system%df_2body, system%df_3body, system%df_4body &
       )
    end do
    close(unit)

  end subroutine read_dfs
!###############################################################################


!###############################################################################
  subroutine write_2body(this, file)
    !! Write the learned 2-body distribution functions to a file.
    implicit none

    ! Arguments
    class(distribs_container_type), intent(in) :: this
    !! Parent. Instance of distribution functions container.
    character(*), intent(in) :: file
    !! Filename to write the 2-body distribution functions to.

    ! Local variables
    integer :: unit
    !! File unit.
    integer :: i, j, is, js
    !! Loop indices.
    integer :: num_pairs
    !! Number of pairs.
    integer, allocatable, dimension(:,:) :: idx
    !! Pair indices.


    num_pairs = nint( &
         gamma(real(size(this%element_info) + 2, real32)) / &
         ( gamma(real(size(this%element_info), real32)) * gamma( 3._real32 ) ) &
    )
    allocate(idx(2,num_pairs))
    i = 0
    do is = 1, size(this%element_info)
       do js = is, size(this%element_info), 1
          i = i + 1
          idx(:,i) = [is, js]
       end do
    end do

    open(newunit=unit, file=file)
    do i = 1,  size(this%gdf%df_2body, dim=2)
       write(unit,'("# ",A,2X,A)') &
            this%element_info(idx(1,i))%name, &
            this%element_info(idx(2,i))%name
       do j = 1, size(this%gdf%df_2body, dim=1)
          write(unit,*) &
               this%cutoff_min(1) + this%width(1) * ( j - 1 ), &
               this%gdf%df_2body(j,i)
       end do
       write(unit,*)
    end do
    close(unit)

  end subroutine write_2body
!###############################################################################


!###############################################################################
  subroutine write_3body(this, file)
    !! Write the learned 3-body distribution functions to a file.
    implicit none

    ! Arguments
    class(distribs_container_type), intent(in) :: this
    !! Parent. Instance of distribution functions container.
    character(*), intent(in) :: file
    !! Filename to write the 3-body distribution functions to.

    ! Local variables
    integer :: unit
    !! File unit.
    integer :: i, j
    !! Loop indices.


    open(newunit=unit, file=file)
    do i = 1,  size(this%gdf%df_3body, dim=2)
       write(unit,'("# ",A)') this%element_info(i)%name
       do j = 1, size(this%gdf%df_3body, dim=1)
          write(unit,*) &
               this%cutoff_min(2) + this%width(2) * ( j - 1 ), &
               this%gdf%df_3body(j,i)
       end do
       write(unit,*)
    end do
    close(unit)

  end subroutine write_3body
!###############################################################################


!###############################################################################
  subroutine write_4body(this, file)
    !! Write the learned 4-body distribution functions to a file.
    implicit none

    ! Arguments
    class(distribs_container_type), intent(in) :: this
    !! Parent. Instance of distribution functions container.
    character(*), intent(in) :: file
    !! Filename to write the 4-body distribution functions to.

    ! Local variables
    integer :: unit
    !! File unit.
    integer :: i, j
    !! Loop indices.


    open(newunit=unit, file=file)
    do i = 1,  size(this%gdf%df_4body, dim=2)
       write(unit,'("# ",A)') this%element_info(i)%name
       do j = 1, size(this%gdf%df_4body, dim=1)
          write(unit,*) &
               this%cutoff_min(3) + this%width(3) * ( j - 1 ), &
               this%gdf%df_4body(j,i)
       end do
       write(unit,*)
    end do
    close(unit)

  end subroutine write_4body
!###############################################################################


!###############################################################################
  subroutine add(this, system)
    !! Add a system to the container.
    implicit none

    ! Arguments
    class(distribs_container_type), intent(inout) :: this
    !! Parent. Instance of distribution functions container.
    class(*), dimension(..), intent(in) :: system
    !! System to add to the container.

    ! Local variables
    integer :: i
    !! Loop index.
    integer :: num_structures_previous
    !! Number of structures in the container before adding the system.
    character(256) :: stop_msg
    !! Error message.

    select rank(rank_ptr => system)
    rank(0)
       select type(type_ptr => rank_ptr)
       type is (distribs_type)
          this%system = [ this%system, type_ptr ]
       class is (basis_type)
#if defined(GFORTRAN)
          call this%add_basis(type_ptr)
#else
          block
            type(basis_type), dimension(1) :: basis

            basis = type_ptr
            call this%add_basis(basis(1))
          end block
#endif
       class default
          write(stop_msg,*) "Invalid type for system" // &
               achar(13) // achar(10) // &
               "Expected type distribs_type or basis_type"
          call stop_program( stop_msg )
          return
       end select
    rank(1)
       num_structures_previous = size(this%system)
       if(.not.allocated(this%system))then
          allocate(this%system(0))
       end if
       select type(type_ptr => rank_ptr)
       type is (distribs_type)
          this%system = [ this%system, type_ptr ]
       class is (basis_type)
          do i = 1, size(type_ptr)
             call this%add_basis(type_ptr(i))
          end do
       class default
          write(stop_msg,*) "Invalid type for system" // &
               achar(13) // achar(10) // &
               "Expected type distribs_type or basis_type"
          call stop_program( stop_msg )
          return
       end select
    rank default
       write(stop_msg,*) "Invalid rank for system" // &
            achar(13) // achar(10) // &
            "Expected rank 0 or 1, got ", rank(rank_ptr)
       call stop_program( stop_msg )
       return
    end select
    call this%update_element_info()
    call this%update_bond_info()

  end subroutine add
!###############################################################################


!###############################################################################
  subroutine add_basis(this, basis)
    !! Add a basis to the container.
    implicit none

    ! Arguments
    class(distribs_container_type), intent(inout) :: this
    !! Parent. Instance of distribution functions container.
    type(basis_type), intent(in) :: basis
    !! Basis to add to the container.

    ! Local variables
    type(distribs_type) :: system
    !! System to add to the container.

    call system%calculate( &
         basis, &
         width = this%width, &
         sigma = this%sigma, &
         cutoff_min = this%cutoff_min, &
         cutoff_max = this%cutoff_max, &
         radius_distance_tol = this%radius_distance_tol &
    )

    if(.not.allocated(this%system))then
       this%system = [ system ]
    else
       this%system = [ this%system, system ]
    end if
  end subroutine add_basis
!###############################################################################


!###############################################################################
  subroutine set_element_map(this, element_list)
    !! Set the element map for the container.
    implicit none

    ! Arguments
    class(distribs_container_type), intent(inout) :: this
    !! Parent of the procedure. Instance of distribution functions container.
    character(3), dimension(:), intent(in) :: element_list
    !! Element information.

    ! Local variables
    integer :: is
    !! Loop index.
    integer :: idx
    !! Index of the element in the element_info array.
    character(256) :: stop_msg
    !! Error message.

    if(allocated(this%element_map)) deallocate(this%element_map)
    allocate(this%element_map(size(element_list)))
    do is = 1, size(element_list)
       idx = findloc( &
            [ this%element_info(:)%name ], element_list(is), dim=1 &
       )
       if(idx.eq.0)then
          write(stop_msg,*) "Element not found in element_info array" // &
               achar(13) // achar(10) // &
               "Element: ", element_list(is)
          call stop_program( stop_msg )
          return
       end if
       this%element_map(is) = idx
    end do

  end subroutine set_element_map
!###############################################################################


!###############################################################################
  subroutine set_element_info(this)
    !! Set the list of elements for the container.
    implicit none

    ! Arguments
    class(distribs_container_type), intent(inout) :: this
    !! Parent of the procedure. Instance of distribution functions container.

    ! Local variables
    integer :: i
    !! Loop index.
    real(real32) :: radius
    !! Element radii.
    character(len=3), dimension(:), allocatable :: element_list
    !! List of elements in the container.


    !---------------------------------------------------------------------------
    ! get list of species in dataset
    !---------------------------------------------------------------------------
    allocate(element_list(0))
    element_list = [ this%system(1)%element_symbols ]
    do i = 2, size(this%system),1
       element_list = [ element_list, this%system(i)%element_symbols ]
    end do
    call set(element_list)
    if(allocated(this%element_info)) deallocate(this%element_info)
    if(.not.allocated(element_database)) allocate(element_database(0))
    allocate(this%element_info(size(element_list)))
    do i = 1, size(element_list)
       if( findloc( &
            [ element_database(:)%name ], element_list(i), dim=1 ) .lt. 1 )then
          call get_element_properties(element_list(i), radius=radius)
          element_database = [ &
               element_database(:), &
               element_type(name = element_list(i), radius = radius) &
          ]
       end if
       call this%element_info(i)%set(element_list(i))
    end do

  end subroutine set_element_info
!###############################################################################


!###############################################################################
  subroutine update_element_info(this)
    !! Update the element information in the container.
    implicit none

    ! Arguments
    class(distribs_container_type), intent(inout) :: this
    !! Parent of the procedure. Instance of distribution functions container.

    ! Local variables
    integer :: i
    !! Loop index.
    real(real32) :: radius
    !! Element radii.
    character(len=3), dimension(:), allocatable :: element_list
    !! List of elements in the container.


    !---------------------------------------------------------------------------
    ! check if element_info is allocated, if not, set it and return
    !---------------------------------------------------------------------------
    if(.not.allocated(this%element_info))then
       call this%set_element_info()
       return
    elseif(size(this%element_info).eq.0)then
       call this%set_element_info()
       return
    end if


    !---------------------------------------------------------------------------
    ! get list of species in dataset
    !---------------------------------------------------------------------------
    allocate(element_list(0))
    element_list = [ this%system(1)%element_symbols ]
    do i = 2, size(this%system),1
       element_list = [ element_list, this%system(i)%element_symbols ]
    end do
    call set(element_list)


    !---------------------------------------------------------------------------
    ! check if all elements are in the element_info array
    !---------------------------------------------------------------------------
    if(.not.allocated(element_database)) allocate(element_database(0))
    do i = 1, size(element_list)
       if(findloc( &
            [ element_database(:)%name ], &
            element_list(i), dim=1 ) .lt. 1 )then
          call get_element_properties(element_list(i), radius=radius)
          element_database = [ &
               element_database(:), &
               element_type(name = element_list(i), radius = radius) &
          ]
       end if
       if( findloc( &
            [ this%element_info(:)%name ], &
            element_list(i), dim=1 ) .lt. 1 )then
          this%element_info = [ &
               this%element_info(:), &
               element_type(element_list(i)) &
          ]
          call this%element_info(size(this%element_info))%set(element_list(i))
       end if
    end do

  end subroutine update_element_info
!###############################################################################


!###############################################################################
  subroutine set_element_energy(this, element, energy)
    !! Set the energy of an element in the container.
    implicit none

    ! Arguments
    class(distribs_container_type), intent(inout) :: this
    !! Parent of the procedure. Instance of distribution functions container.
    character(len=3), intent(in) :: element
    !! Element name.
    real(real32), intent(in) :: energy
    !! Energy of the element.

    ! Local variables
    integer :: idx, idx_db
    !! Index of the element in the element_info array.
    real(real32) :: radius
    !! Element radius.
    character(len=3) :: element_
    !! Element name without null characters.


    !---------------------------------------------------------------------------
    ! remove python formatting
    !---------------------------------------------------------------------------
    element_ = strip_null(element)


    !---------------------------------------------------------------------------
    ! if element_info is allocated, update energy of associated index
    !---------------------------------------------------------------------------
    if(allocated(this%element_info))then
       idx = findloc( [ this%element_info(:)%name ], element_, dim=1 )
       if(idx.ge.1) this%element_info(idx)%energy = energy
    end if


    !---------------------------------------------------------------------------
    ! if element_database is allocated, update energy of associated index
    !---------------------------------------------------------------------------
    if(.not.allocated(element_database)) allocate(element_database(0))
    idx_db = findloc( [ element_database(:)%name ], element_, dim=1 )
    if(idx_db.lt.1)then
       call get_element_properties( element_, radius = radius )
       element_database = [ &
            element_database(:), &
            element_type(name = element_, energy = energy, radius = radius) &
       ]
    else
       element_database(idx_db)%energy = energy
    end if

  end subroutine set_element_energy
!###############################################################################


!###############################################################################
  subroutine set_element_energies(this, elements, energies)
    !! Set the energies of elements in the container.
    implicit none

    ! Arguments
    class(distribs_container_type), intent(inout) :: this
    !! Parent of the procedure. Instance of distribution functions container.
    character(len=3), dimension(:), intent(in) :: elements
    !! Element names.
    real(real32), dimension(:), intent(in) :: energies
    !! Energies of the elements.

    ! Local variables
    integer :: i

    do i = 1, size(elements)
       call this%set_element_energy(elements(i), energies(i))
    end do

  end subroutine set_element_energies
!###############################################################################


!###############################################################################
  subroutine get_element_energies(this, elements, energies)
    !! Return the energies of elements in the container.
    implicit none

    ! Arguments
    class(distribs_container_type), intent(in) :: this
    !! Parent of the procedure. Instance of distribution functions container.
    character(len=3), dimension(:), allocatable, intent(out) :: elements
    !! Element names.
    real(real32), dimension(:), allocatable, intent(out) :: energies
    !! Energies of the elements.

    ! Local variables
    integer :: i
    !! Loop index.


    allocate(elements(size(this%element_info)))
    allocate(energies(size(this%element_info)))
    do i = 1, size(this%element_info)
       elements(i) = this%element_info(i)%name
       energies(i) = this%element_info(i)%energy
    end do

  end subroutine get_element_energies
!###############################################################################


!###############################################################################
  subroutine get_element_energies_staticmem(this, elements, energies)
    !! Return the energies of elements in the container.
    !!
    !! This subroutine is used when the memory for the output arrays is
    !! allocated outside of the subroutine. Used in Python interface.
    implicit none

    ! Arguments
    class(distribs_container_type), intent(in) :: this
    !! Parent of the procedure. Instance of distribution functions container.
    character(len=3), dimension(size(this%element_info,1)), intent(out) :: &
         elements
    !! Element names.
    real(real32), dimension(size(this%element_info,1)), intent(out) :: energies
    !! Energies of the elements.

    ! Local variables
    integer :: i
    !! Loop index.


    do i = 1, size(this%element_info,1)
       elements(i) = this%element_info(i)%name
       energies(i) = this%element_info(i)%energy
    end do

  end subroutine get_element_energies_staticmem
!###############################################################################


!###############################################################################
  subroutine set_bond_info(this)
    !! Set the 2-body bond information for the container.
    implicit none

    ! Arguments
    class(distribs_container_type), intent(inout) :: this
    !! Parent of the procedure. Instance of distribution functions container.

    ! Local variables
    integer :: i, j
    !! Loop index.
    integer :: num_elements, num_pairs
    !! Number of elements and pairs.
    logical :: success
    !! Success flag.


    !---------------------------------------------------------------------------
    ! allocate the bond information array
    !---------------------------------------------------------------------------
    num_elements = size(this%element_info)
    num_pairs = nint(gamma(real(num_elements + 2, real32)) / &
         ( gamma(real(num_elements, real32)) * gamma( 3._real32 ) ) )
    if(allocated(this%bond_info)) deallocate(this%bond_info)
    allocate(this%bond_info(num_pairs))


    !---------------------------------------------------------------------------
    ! loop over all pairs of elements to set the bond information
    !---------------------------------------------------------------------------
    num_pairs = 0
    pair_loop1: do i = 1, num_elements
       pair_loop2: do j = i, num_elements
          num_pairs = num_pairs + 1
          call this%bond_info(num_pairs)%set( &
               this%element_info(i)%name, &
               this%element_info(j)%name, &
               success &
          )
          if(success) cycle pair_loop2
          call set_bond_radius_to_default( [ &
               this%element_info(i)%name, &
               this%element_info(j)%name &
          ] )
          call this%bond_info(num_pairs)%set( &
               this%element_info(i)%name, &
               this%element_info(j)%name, &
               success &
          )
       end do pair_loop2
    end do pair_loop1

  end subroutine set_bond_info
!###############################################################################


!###############################################################################
  subroutine set_bond_radius_to_default(elements)
    !! Set the bond radius to the default value.
    !!
    !! The default value is the average of the covalent radii of the elements.
    implicit none

    ! Arguments
    character(len=3), dimension(2), intent(in) :: elements
    !! Element symbols.

    ! Local variables
    integer :: idx1, idx2
    !! Index of the elements in the element database.
    real(real32) :: radius, radius1, radius2
    !! Average of covalent radii.
    character(256) :: warn_msg


    write(warn_msg,'("No bond data for element pair ",A," and ",A)') &
         elements(1), elements(2)
    warn_msg = trim(warn_msg) // &
         achar(13) // achar(10) // &
         "Setting bond to average of covalent radii"
    call print_warning(warn_msg)
    if(.not.allocated(element_database)) allocate(element_database(0))
    idx1 = findloc([ element_database(:)%name ], &
         elements(1), dim=1)
    if(idx1.lt.1)then
       call get_element_properties(elements(1), radius=radius1)
       element_database = [ element_database, &
            element_type(name=elements(1), radius=radius1) ]
       idx1 = size(element_database)
    end if
    idx2 = findloc([ element_database(:)%name ], &
         elements(2), dim=1)
    if(idx2.lt.1)then
       call get_element_properties(elements(2), radius=radius2)
       element_database = [ element_database, &
            element_type(name=elements(2), radius=radius2) ]
       idx2 = size(element_database)
    end if
    radius = ( element_database(idx1)%radius + &
         element_database(idx2)%radius ) / 2._real32
    if(.not.allocated(element_bond_database)) &
         allocate(element_bond_database(0))
    element_bond_database = [ element_bond_database, &
         element_bond_type(elements=[ &
              elements(1), &
              elements(2) &
         ], radius=radius) &
    ]
    call sort_str( &
         element_bond_database(size(element_bond_database))%element &
    )

  end subroutine set_bond_radius_to_default
!###############################################################################


!###############################################################################
  subroutine update_bond_info(this)
    !! Update the element information in the container.
    implicit none

    ! Arguments
    class(distribs_container_type), intent(inout) :: this
    !! Parent of the procedure. Instance of distribution functions container.

    ! Local variables
    integer :: i, j, k, is, js
    !! Loop index.
    real(real32) :: radius, radius1, radius2
    !! Average of covalent radii.
    character(len=3), dimension(:), allocatable :: element_list
    !! List of elements in the container.
    character(len=3), dimension(:,:), allocatable :: pair_list
    !! List of element pairs in the container.


    !---------------------------------------------------------------------------
    ! check if bond_info is allocated, if not, set it and return
    !---------------------------------------------------------------------------
    if(.not.allocated(this%bond_info))then
       call this%set_bond_info()
       return
    elseif(size(this%bond_info).eq.0)then
       call this%set_bond_info()
       return
    end if


    !---------------------------------------------------------------------------
    ! get list of element pairs in dataset
    !---------------------------------------------------------------------------
    allocate(element_list(0))
    element_list = [ this%system(1)%element_symbols ]
    do i = 2, size(this%system),1
       element_list = [ element_list, this%system(i)%element_symbols ]
    end do
    call set(element_list)
    allocate(pair_list(triangular_number(size(element_list)),2))
    k = 0
    do i = 1, size(element_list)
       do j = i, size(element_list)
          k = k + 1
          pair_list(k,:) = [ element_list(i), element_list(j) ]
          call sort_str(pair_list(k,:))
       end do
    end do

    !---------------------------------------------------------------------------
    ! check if all element pairs are in the database
    !---------------------------------------------------------------------------
    if(.not.allocated(element_bond_database)) allocate(element_bond_database(0))
    pair_loop1: do i = 1, size(pair_list,1)
       do j = 1, size(element_bond_database)
          if( ( &
               element_bond_database(j)%element(1) .eq. pair_list(i,1) .and. &
               element_bond_database(j)%element(2) .eq. pair_list(i,2) &
          ) .or. ( &
               element_bond_database(j)%element(1) .eq. pair_list(i,2) .and. &
               element_bond_database(j)%element(2) .eq. pair_list(i,1) &
          ) ) cycle pair_loop1
       end do
       ! if not found, add to the database
       is = findloc([ this%element_info(:)%name ], pair_list(i,1), dim=1)
       js = findloc([ this%element_info(:)%name ], pair_list(i,2), dim=1)
       radius1 = this%element_info(is)%radius
       if(radius1.le.1.E-6) &
            call get_element_properties(pair_list(i,1), radius = radius1)
       radius2 = this%element_info(js)%radius
       if(radius2.le.1.E-6) &
            call get_element_properties(pair_list(i,2), radius = radius2)
       radius = ( radius1 + radius2 ) / 2._real32
       element_bond_database = [ element_bond_database, &
            element_bond_type(elements=[pair_list(i,:)], radius=radius) ]
       call sort_str(element_bond_database(size(element_bond_database))%element)
    end do pair_loop1


    ! --------------------------------------------------------------------------
    ! check if all element pairs are in the bond_info array
    ! --------------------------------------------------------------------------
    pair_loop2: do i = 1, size(pair_list,1)
       info_loop: do j = 1, size(this%bond_info)
          if( ( &
               this%bond_info(j)%element(1) .eq. pair_list(i,1) .and. &
               this%bond_info(j)%element(2) .eq. pair_list(i,2) &
          ) .or. ( &
               this%bond_info(j)%element(1) .eq. pair_list(i,2) .and. &
               this%bond_info(j)%element(2) .eq. pair_list(i,1) &
          ) ) cycle pair_loop2
       end do info_loop
       this%bond_info = [ &
            this%bond_info(:), &
            element_bond_type( [ pair_list(i,1:2) ] ) &
       ]
       call this%bond_info(size(this%bond_info))%set( &
            pair_list(i,1), &
            pair_list(i,2) )
    end do pair_loop2

  end subroutine update_bond_info
!###############################################################################


!###############################################################################
  subroutine set_bond_radius(this, elements, radius)
    !! Set the bond radius for a pair of elements in the container.
    implicit none

    ! Arguments
    class(distribs_container_type), intent(inout) :: this
    !! Parent of the procedure. Instance of distribution functions container.
    character(len=3), dimension(2), intent(in) :: elements
    !! Element name.
    real(real32), intent(in) :: radius
    !! Bond radius.

    ! Local variables
    integer :: idx, i
    !! Index of the bond in the bond_info array.
    character(len=3) :: element_1, element_2
    !! Element names.


    !---------------------------------------------------------------------------
    ! remove python formatting
    !---------------------------------------------------------------------------
    element_1 = strip_null(elements(1))
    element_2 = strip_null(elements(2))


    !---------------------------------------------------------------------------
    ! if bond_info is allocated, update radius of associated index
    !---------------------------------------------------------------------------
    if(allocated(this%bond_info))then
       idx = this%get_pair_index(element_1, element_2)
       if(idx.ge.1) this%bond_info(idx)%radius_covalent = radius
    end if


    !---------------------------------------------------------------------------
    ! if element_bond_database is allocated, update radius of associated index
    !---------------------------------------------------------------------------
    if(.not.allocated(element_bond_database))then
       allocate(element_bond_database(1))
       element_bond_database(1)%element = [ element_1, element_2 ]
       call sort_str(element_bond_database(1)%element)
       element_bond_database(1)%radius_covalent = radius
       return
    end if
    do i = 1, size(element_bond_database)
       if( ( &
            element_bond_database(i)%element(1) .eq. element_1 .and. &
            element_bond_database(i)%element(2) .eq. element_2 &
       ) .or. ( &
            element_bond_database(i)%element(1) .eq. element_2 .and. &
            element_bond_database(i)%element(2) .eq. element_1 &
       ) )then
          element_bond_database(i)%radius_covalent = radius
          return
       end if
    end do
    ! if allocated and not found, add to the database
    element_bond_database = [ element_bond_database, &
         element_bond_type([ element_1, element_2 ], radius) ]
    call sort_str(element_bond_database(size(element_bond_database))%element)

  end subroutine set_bond_radius
!###############################################################################


!###############################################################################
  subroutine set_bond_radii(this, elements, radii)
    !! Set the bond radii for a pair of elements in the container.
    implicit none

    ! Arguments
    class(distribs_container_type), intent(inout) :: this
    !! Parent of the procedure. Instance of distribution functions container.
    character(len=3), dimension(:,:), intent(in) :: elements
    !! Element names.
    real(real32), dimension(:), intent(in) :: radii
    !! Bond radii.

    ! Local variables
    integer :: i
    !! Loop index.


    do i = 1, size(elements,1)
       call this%set_bond_radius(elements(i,:), radii(i))
    end do

  end subroutine set_bond_radii
!###############################################################################


!###############################################################################
  subroutine get_bond_radii(this, elements, radii)
    !! Return the bond radii of elements in the container.
    implicit none

    ! Arguments
    class(distribs_container_type), intent(in) :: this
    !! Parent of the procedure. Instance of distribution functions container.
    character(len=3), dimension(:,:), allocatable, intent(out) :: elements
    !! Element pair names.
    real(real32), dimension(:), allocatable, intent(out) :: radii
    !! Radii of the bond pairs.

    ! Local variables
    integer :: i
    !! Loop index.


    allocate(elements(size(this%bond_info),2))
    allocate(radii(size(this%bond_info)))
    do i = 1, size(this%bond_info)
       elements(i,:) = this%bond_info(i)%element
       radii(i) = this%bond_info(i)%radius_covalent
    end do

  end subroutine get_bond_radii
!###############################################################################


!###############################################################################
  subroutine get_bond_radii_staticmem(this, elements, radii)
    !! Return the energies of elements in the container.
    !!
    !! This subroutine is used when the memory for the output arrays is
    !! allocated outside of the subroutine. Used in Python interface.
    implicit none

    ! Arguments
    class(distribs_container_type), intent(in) :: this
    !! Parent of the procedure. Instance of distribution functions container.
    character(len=3), dimension(size(this%bond_info,1),2), intent(out) :: &
         elements
    !! Element pair names.
    real(real32), dimension(size(this%bond_info,1)), intent(out) :: radii
    !! Radii of the bond pairs.

    ! Local variables
    integer :: i
    !! Loop index.


    do i = 1, size(this%bond_info)
       elements(i,:) = this%bond_info(i)%element
       radii(i) = this%bond_info(i)%radius_covalent
    end do

  end subroutine get_bond_radii_staticmem
!###############################################################################


!###############################################################################
  subroutine set_best_energy(this)
    !! Set the best energy in the container.
    implicit none

    ! Arguments
    class(distribs_container_type), intent(inout) :: this
    !! Parent of the procedure. Instance of distribution functions container.

    ! Local variables
    integer :: i, j, is, js, idx1, idx2
    !! Loop index.
    real(real32) :: energy, energy_per_species, energy_pair
    !! Energy of the system.
    integer, dimension(:,:), allocatable :: idx_list
    !! Index list for pairs of elements.
    character(len=256) :: warn_msg
    !! Warning message.

    if(.not.allocated(this%best_energy_pair))then
       allocate( &
            this%best_energy_pair(size(this%bond_info,1)), &
            source = 0._real32 &
       )
    elseif(size(this%best_energy_pair).ne.size(this%bond_info))then
       deallocate(this%best_energy_pair)
       allocate( &
            this%best_energy_pair(size(this%bond_info,1)), &
            source = 0._real32 &
       )
    end if

    if(.not.allocated(this%best_energy_per_species))then
       allocate( &
            this%best_energy_per_species(size(this%element_info,1)), &
            source = 0._real32 &
       )
    elseif(size(this%best_energy_per_species).ne.size(this%element_info))then
       deallocate(this%best_energy_per_species)
       allocate( &
            this%best_energy_per_species(size(this%element_info,1)), &
            source = 0._real32 &
       )
    end if

    do i = 1, size(this%system)
       j = 0
       allocate( idx_list( &
            size(this%system(i)%element_symbols), &
            size(this%system(i)%element_symbols) &
       ) )
       do is = 1, size(this%system(i)%element_symbols)
          do js = is, size(this%system(i)%element_symbols), 1
             j = j + 1
             idx_list(is,js) = j
             idx_list(js,is) = j
          end do
       end do

       energy = this%system(i)%energy
       do is = 1, size(this%system(i)%element_symbols)
          idx1 = findloc( &
               [ this%element_info(:)%name ], &
               this%system(i)%element_symbols(is), &
               dim = 1 &
          )
          if(idx1.lt.1)then
             call stop_program( "Species not found in element_info" )
             return
          end if
          energy = energy - &
               this%system(i)%stoichiometry(is) * this%element_info(idx1)%energy
       end do
       energy = energy / this%system(i)%num_atoms

       do is = 1, size(this%system(i)%element_symbols)
          if(this%system(i)%num_per_species(is).eq.0)then
             write(warn_msg, &
                  '("No neighbours found for species ",A," (",I0,") &
                  &in system ",I0)' &
             ) trim(this%system(i)%element_symbols(is)), is, i
             call print_warning(warn_msg)
             cycle
          end if
          idx1 = findloc( &
               [ this%element_info(:)%name ], &
               this%system(i)%element_symbols(is), &
               dim = 1 &
          )
          energy_per_species = &
               energy * this%system(i)%weight_per_species(is) / &
               real( sum( this%system(i)%num_per_species(:) ), real32 )

          if( energy_per_species .lt. this%best_energy_per_species(idx1) )then
             this%best_energy_per_species(idx1) = energy_per_species
          end if
          do js = 1, size(this%system(i)%element_symbols)
             idx2 = findloc( &
                  [ this%element_info(:)%name ], &
                  this%system(i)%element_symbols(js), &
                  dim = 1 &
             )
             j = nint( &
                  ( &
                       size(this%element_info) - &
                       min( idx1, idx2 ) / 2._real32 &
                  ) * ( min( idx1, idx2 ) - 1._real32 ) + max( idx1, idx2 ) &
             )

             energy_pair = &
                  energy * this%system(i)%weight_pair(idx_list(is,js)) / &
                  real( sum( this%system(i)%num_per_species(:) ), real32 )

             if( energy_pair .lt. this%best_energy_pair(j) )then
                this%best_energy_pair(j) = energy_pair
             end if

          end do
       end do
       deallocate(idx_list)
       if( this%best_energy_pair(j) .lt. -1.E1 )then
          write(warn_msg, &
               '("Best energy pair is less than -10 eV, &
               &this is likely to be unphysical. Check the energy values.")' &
          )
          call print_warning(warn_msg)
       end if

    end do

  end subroutine set_best_energy
!###############################################################################


!###############################################################################
  pure function get_pair_index(this, species1, species2) result(idx)
    !! Get the index of a pair of elements in the container.
    implicit none

    ! Arguments
    class(distribs_container_type), intent(in) :: this
    !! Parent of the procedure. Instance of distribution functions container.
    character(len=3), intent(in) :: species1, species2
    !! Element names.
    integer :: idx
    !! Index of the pair in the bond_info array.

    ! Local variables
    integer :: is, js
    !! Index of the elements in the element_info array.

    is = findloc([ this%element_info(:)%name ], species1, dim=1)
    js = findloc([ this%element_info(:)%name ], species2, dim=1)
    !! This comes from subtraction of nth triangular numbers
    !! nth triangular number: N_n = n(n+1)/2 = ( n^2 + n ) / 2
    !! idx = N_n - N_{n-is+1} + ( js - is + 1)
    !! idx = ( n - is/2 ) * ( is - 1 ) + js
    !idx = nint( ( size(this%element_info) - min( is, js ) / 2._real32 ) * &
    !      ( is - 1._real32 ) + js )
    idx = nint( ( size(this%element_info) - min( is, js ) / 2._real32 ) * &
         ( min( is, js ) - 1._real32 ) + max( is, js ) )

  end function get_pair_index
!###############################################################################


!###############################################################################
  pure function get_element_index(this, species) result(idx)
    !! Get the index of an element in the container.
    implicit none

    ! Arguments
    class(distribs_container_type), intent(in) :: this
    !! Parent of the procedure. Instance of distribution functions container.
    character(len=3), intent(in) :: species
    !! Element name.
    integer :: idx
    !! Index of the element in the element_info array.

    ! Local variables

    idx = findloc([ this%element_info(:)%name ], species, dim=1)

  end function get_element_index
!###############################################################################


!###############################################################################
  pure function get_bin(this, value, dim) result(bin)
    !! Get the bin index for a value in a dimension.
    implicit none

    ! Arguments
    class(distribs_container_type), intent(in) :: this
    !! Parent of the procedure. Instance of distribution functions container.
    real(real32), intent(in) :: value
    !! Value to get the bin index for.
    integer, intent(in) :: dim
    !! Dimension to get the bin index for.
    integer :: bin
    !! Bin index for the value.

    if(value .lt. this%cutoff_min(dim) - this%width(dim) .or. &
         value .gt. this%cutoff_max(dim) + this%width(dim))then
       bin = 0
    else
       bin = nint( &
            ( this%nbins(dim) - 1 ) * &
            ( value - this%cutoff_min(dim) ) / &
            ( this%cutoff_max(dim) - this%cutoff_min(dim) ) &
       ) + 1
    end if

  end function get_bin
!###############################################################################


!###############################################################################
  subroutine initialise_gdfs(this)
    !! Initialise the distribution functions for the container.
    implicit none

    ! Arguments
    class(distribs_container_type), intent(inout) :: this
    !! Parent of the procedure. Instance of distribution functions container.

    ! Local variables
    integer :: num_pairs
    !! Number of pairs.


    num_pairs = nint( gamma(real(size(this%element_info) + 2, real32)) / &
         ( gamma(real(size(this%element_info), real32)) * gamma( 3._real32 ) ) )
    allocate(this%gdf%df_2body(this%nbins(1),num_pairs), &
         source = 0._real32 )
    allocate(this%gdf%df_3body(this%nbins(2),size(this%element_info)), &
         source = 0._real32 )
    allocate(this%gdf%df_4body(this%nbins(3),size(this%element_info)), &
         source = 0._real32 )
    allocate(this%in_dataset_2body(num_pairs), source = .false. )
    allocate(this%in_dataset_3body(size(this%element_info)), source = .false. )
    allocate(this%in_dataset_4body(size(this%element_info)), source = .false. )

  end subroutine initialise_gdfs
!###############################################################################


!###############################################################################
  subroutine set_gdfs_to_default(this, body, index)
    !! Initialise the gdfs for index of body distribution function.
    implicit none

    ! Arguments
    class(distribs_container_type), intent(inout) :: this
    !! Parent of the procedure. Instance of distribution functions container.
    integer, intent(in) :: body
    !! Body distribution function to initialise.
    integer, intent(in) :: index
    !! Index of the pair in the bond_info array.

    ! Local variables
    real(real32) :: eta, weight, height
    !! Parameters for the distribution functions.
    real(real32), dimension(1) :: bonds


    if( body .eq. 2 )then
       weight = exp( -4._real32 )
       height = 1._real32 / this%nbins(1)
       eta = 1._real32 / ( 2._real32 * ( this%sigma(1) )**2._real32 )
       if(size(this%bond_info).eq.0)then
          call set_bond_radius_to_default( [ &
               this%bond_info(index)%element(1), &
               this%bond_info(index)%element(2) &
          ] )
       end if
       bonds = [ 2._real32 * this%bond_info(index)%radius_covalent ]
       if(abs(bonds(1)).lt.1.E-6)then
          call stop_program( "Bond radius is zero" )
       end if
       this%gdf%df_2body(:,index) = weight * height * get_distrib( &
            bonds , &
            this%nbins(1), eta, this%width(1), &
            this%cutoff_min(1), &
            scale_list = [ 1._real32 ] &
       )
    elseif( body .eq. 3 )then
       this%gdf%df_3body(:,index) = 1._real32/this%nbins(2)
    elseif( body .eq. 4 )then
       this%gdf%df_4body(:,index) = 1._real32/this%nbins(3)
    end if

  end subroutine set_gdfs_to_default
!###############################################################################


!###############################################################################
  subroutine evolve(this, system)
    !! Evolve the generalised distribution functions for the container.
    implicit none

    ! Arguments
    class(distribs_container_type), intent(inout) :: this
    !! Parent of the procedure. Instance of distribution functions container.
    type(distribs_type), dimension(..), intent(in), optional :: system
    !! Optional. System to add to the container.

    ! Local variables
    integer :: i, j, is, js
    !! Loop index.
    integer :: idx1, idx2
    !! Index of the element in the element_info array.
    integer :: num_evaluated
    !! Number of systems evaluated this iteration.
    real(real32) :: weight, energy
    !! Energy and weight variables for a system.
    real(real32), dimension(:), allocatable :: &
         best_energy_pair_old, &
         best_energy_per_species_old
    !! Old best energies.
    integer, dimension(:,:), allocatable :: idx_list
    !! Index list for the element pairs in a system.
    real(real32), dimension(:,:), allocatable :: tmp_df
    !! Temporary array for the distribution functions.
    logical, dimension(:), allocatable :: tmp_in_dataset

    character(256) :: err_msg
    integer, dimension(:), allocatable :: host_idx_list


    weight = 1._real32

    !---------------------------------------------------------------------------
    ! if present, add the system to the container
    !---------------------------------------------------------------------------
    if(present(system)) call this%add(system)
    do i = 1, 3
       if(this%nbins(i).eq.-1)then
          this%nbins(i) = 1 + &
               nint( &
                    ( this%cutoff_max(i) - this%cutoff_min(i) ) / &
                    this%width(i) &
               )
       end if
    end do


    !---------------------------------------------------------------------------
    ! initialise the generalised distribution functions and get
    ! best energies from lowest formation energy system
    !---------------------------------------------------------------------------
    if(.not.allocated(this%gdf%df_2body))then
       call this%set_best_energy()
       call this%initialise_gdfs()
    else
       best_energy_pair_old = this%best_energy_pair
       best_energy_per_species_old = this%best_energy_per_species
       call this%set_best_energy()
       do i = 1, size(this%gdf%df_2body,2)
          this%gdf%df_2body(:,i) = this%gdf%df_2body(:,i) * &
               exp( this%best_energy_pair(i) / this%kBT ) / &
               exp( best_energy_pair_old(i) / this%kBT )
       end do
       do i = 1, size(this%gdf%df_3body,2)
          this%gdf%df_3body(:,i) = &
               this%gdf%df_3body(:,i) * exp( &
                    this%best_energy_per_species(i) / this%kBT &
               ) / exp( best_energy_per_species_old(i) / this%kBT )
       end do
       do i = 1, size(this%gdf%df_4body,2)
          this%gdf%df_4body(:,i) = &
               this%gdf%df_4body(:,i) * exp( &
                    this%best_energy_per_species(i) / this%kBT &
               ) / exp( best_energy_per_species_old(i) / this%kBT )
       end do
       if(size(this%gdf%df_2body,2).ne.size(this%bond_info))then
          allocate(tmp_df(this%nbins(1),size(this%bond_info)), &
               source = 0._real32 )
          tmp_df(:,1:size(this%gdf%df_2body,2)) = this%gdf%df_2body
          deallocate(this%gdf%df_2body)
          call move_alloc( tmp_df, this%gdf%df_2body )
          allocate(tmp_in_dataset(size(this%bond_info)), source = .false. )
          tmp_in_dataset(1:size(this%in_dataset_2body)) = this%in_dataset_2body
          deallocate(this%in_dataset_2body)
          call move_alloc( tmp_in_dataset, this%in_dataset_2body )
       end if
       if(size(this%gdf%df_3body,2).ne.size(this%element_info))then
          allocate(tmp_df(this%nbins(2),size(this%element_info)), &
               source = 0._real32 )
          tmp_df(:,1:size(this%gdf%df_3body,2)) = this%gdf%df_3body
          deallocate(this%gdf%df_3body)
          call move_alloc( tmp_df, this%gdf%df_3body )
          allocate(tmp_in_dataset(size(this%element_info)), source = .false. )
          tmp_in_dataset(1:size(this%in_dataset_3body)) = this%in_dataset_3body
          deallocate(this%in_dataset_3body)
          call move_alloc( tmp_in_dataset, this%in_dataset_3body )
       end if
       if(size(this%gdf%df_4body,2).ne.size(this%element_info))then
          allocate(tmp_df(this%nbins(3),size(this%element_info)), &
               source = 0._real32 )
          tmp_df(:,1:size(this%gdf%df_4body,2)) = this%gdf%df_4body
          deallocate(this%gdf%df_4body)
          call move_alloc( tmp_df, this%gdf%df_4body )
          allocate(tmp_in_dataset(size(this%element_info)), source = .false. )
          tmp_in_dataset(1:size(this%in_dataset_4body)) = this%in_dataset_4body
          deallocate(this%in_dataset_4body)
          call move_alloc( tmp_in_dataset, this%in_dataset_4body )
       end if
       do j = 1, size(this%gdf%df_2body,2)
          if(.not.this%in_dataset_2body(j))then
             this%gdf%df_2body(:,j) = 0._real32
          else
             this%gdf%df_2body(:,j) = &
                  this%gdf%df_2body(:,j) * this%norm_2body(j)
          end if
       end do
       do is = 1, size(this%element_info)
          if(.not.this%in_dataset_3body(is))then
             this%gdf%df_3body(:,is) = 0._real32
          else
             this%gdf%df_3body(:,is) = &
                  this%gdf%df_3body(:,is) * this%norm_3body(is)
          end if
          if(.not.this%in_dataset_4body(is))then
             this%gdf%df_4body(:,is) = 0._real32
          else
             this%gdf%df_4body(:,is) = &
                  this%gdf%df_4body(:,is) * this%norm_4body(is)
          end if
       end do
       deallocate(this%norm_2body)
       deallocate(this%norm_3body)
       deallocate(this%norm_4body)
    end if

    if( &
         any(this%system(this%num_evaluated_allocated+1:)%from_host) .and. &
         this%host_system%defined &
    )then
       ! set host_idx_list
       allocate(host_idx_list(size(this%element_info)))
       host_idx_list = 0
       do is = 1, size(this%host_system%element_symbols)
          idx1 = findloc( &
               [ this%element_info(:)%name ], &
               this%host_system%element_symbols(is), &
               dim = 1 &
          )
          if(idx1.lt.1)then
             call stop_program( "Host species not found in species list" )
             return
          end if
          host_idx_list(idx1) = is
       end do
    end if

    !---------------------------------------------------------------------------
    ! loop over all systems to calculate the generalised distribution functions
    !---------------------------------------------------------------------------
    num_evaluated = 0
    do i = this%num_evaluated_allocated + 1, size(this%system), 1
       num_evaluated = num_evaluated + 1
       if(this%weight_by_hull)then
          weight = exp( this%system(i)%energy_above_hull / this%kBT )
          if(weight.lt.1.E-6) cycle
       end if
       !------------------------------------------------------------------------
       ! get the list of 2-body species pairs the system
       !------------------------------------------------------------------------
       j = 0
       allocate( idx_list( &
            size(this%system(i)%element_symbols), &
            size(this%system(i)%element_symbols) &
       ) )
       do is = 1, size(this%system(i)%element_symbols)
          do js = is, size(this%system(i)%element_symbols), 1
             j = j + 1
             idx_list(is,js) = j
             idx_list(js,is) = j
          end do
       end do


       !------------------------------------------------------------------------
       ! calculate the weight for the system
       !------------------------------------------------------------------------
       energy = this%system(i)%energy
       do is = 1, size(this%system(i)%element_symbols)
          idx1 = findloc( &
               [ this%element_info(:)%name ], &
               this%system(i)%element_symbols(is), &
               dim = 1 &
          )
          if(idx1.lt.1)then
             call stop_program( "Species not found in species list" )
             return
          end if
          energy = energy - this%system(i)%stoichiometry(is) * &
               this%element_info(idx1)%energy
       end do
       energy = energy / this%system(i)%num_atoms
       j = 0
       !------------------------------------------------------------------------
       ! loop over all species in the system to add the distributions
       !------------------------------------------------------------------------
       do is = 1, size(this%system(i)%element_symbols)

          if( this%system(i)%num_per_species(is).eq.0 )cycle

          idx1 = findloc( &
               [ this%element_info(:)%name ], &
               this%system(i)%element_symbols(is), &
               dim = 1 &
          )

          if(.not.this%weight_by_hull)then
             weight = exp( &
                  ( &
                       this%best_energy_per_species(is) - &
                       energy * ( &
                            this%system(i)%weight_per_species(is) / &
                            real( &
                                 sum( this%system(i)%num_per_species(:) ), &
                                 real32 &
                            ) &
                       ) &
                  ) / this%kBT &
             )
             if(weight.lt.1.E-6) cycle
          end if

          this%gdf%df_3body(:,idx1) = this%gdf%df_3body(:,idx1) + &
               set_difference( &
                    weight * this%system(i)%df_3body(:,is), &
                    this%gdf%df_3body(:,idx1), &
                    set_min_zero = .true. &
               )

          this%gdf%df_4body(:,idx1) = this%gdf%df_4body(:,idx1) + &
               set_difference( &
                    weight * this%system(i)%df_4body(:,is), &
                    this%gdf%df_4body(:,idx1), &
                    set_min_zero = .true. &
               )

          do js = is, size(this%system(i)%element_symbols), 1
             idx2 = findloc( &
                  [ this%element_info(:)%name ], &
                  this%system(i)%element_symbols(js), &
                  dim = 1 &
             )
             j = nint( ( &
                  size(this%element_info) - min( idx1, idx2 ) / 2._real32 &
             ) * ( min( idx1, idx2 ) - 1._real32 ) + max( idx1, idx2 ) )

             if(.not.this%weight_by_hull)then
                weight = exp( &
                     ( &
                          this%best_energy_pair(j) - &
                          energy * ( &
                               this%system(i)%weight_pair(idx_list(is,js)) / &
                               real( &
                                    sum( this%system(i)%num_per_species(:) ), &
                                    real32 &
                               ) &
                          ) &
                     ) / this%kBT &
                )
                if(weight.lt.1.E-6) cycle
             end if

             this%gdf%df_2body(:,j) = this%gdf%df_2body(:,j) + &
                  set_difference( &
                       weight * this%system(i)%df_2body(:,idx_list(is,js)), &
                       this%gdf%df_2body(:,j), &
                       set_min_zero = .true. &
                  )

          end do
       end do
       deallocate(idx_list)
    end do

    !---------------------------------------------------------------------------
    ! if not in the dataset, set distribution functions to default
    !---------------------------------------------------------------------------
    do j = 1, size(this%gdf%df_2body,2)
       if(all(abs(this%gdf%df_2body(:,j)).lt.1.E-6))then
          call this%set_gdfs_to_default(2, j)
       else
          this%in_dataset_2body(j) = .true.
       end if
    end do
    do is = 1, size(this%element_info)
       if(all(abs(this%gdf%df_3body(:,is)).lt.1.E-6))then
          call this%set_gdfs_to_default(3, is)
       else
          this%in_dataset_3body(is) = .true.
       end if
       if(all(abs(this%gdf%df_4body(:,is)).lt.1.E-6))then
          call this%set_gdfs_to_default(4, is)
       else
          this%in_dataset_4body(is) = .true.
       end if
    end do

    allocate(this%norm_2body(size(this%gdf%df_2body,2)))
    do j = 1, size(this%gdf%df_2body,2)
       this%norm_2body(j) = maxval(this%gdf%df_2body(:,j))
       if(abs(this%norm_2body(j)).lt.1.E-6)then
          call stop_program( "Zero norm for 2-body distribution function" )
          return
       end if
       this%gdf%df_2body(:,j) = &
            this%gdf%df_2body(:,j) / this%norm_2body(j)
       if(any(isnan(this%gdf%df_2body(:,j))))then
          write(err_msg, &
               '("NaN in 2-body distribution function for ",A,"-",A,&
               &" with norm of ",F0.3)' &
          ) &
               this%bond_info(j)%element(1), &
               this%bond_info(j)%element(2), &
               this%norm_2body(j)
          call stop_program( err_msg )
          return
       end if
    end do
    allocate(this%norm_3body(size(this%element_info)))
    allocate(this%norm_4body(size(this%element_info)))
    do is = 1, size(this%element_info)
       this%norm_3body(is) = maxval(this%gdf%df_3body(:,is))
       if(abs(this%norm_3body(is)).lt.1.E-6)then
          call stop_program( "Zero norm for 3-body distribution function" )
          return
       end if
       this%norm_4body(is) = maxval(this%gdf%df_4body(:,is))
       if(abs(this%norm_4body(is)).lt.1.E-6)then
          call stop_program( "Zero norm for 4-body distribution function" )
          return
       end if
       this%gdf%df_3body(:,is) = &
            this%gdf%df_3body(:,is) / this%norm_3body(is)
       this%gdf%df_4body(:,is) = &
            this%gdf%df_4body(:,is) / this%norm_4body(is)

       if(any(isnan(this%gdf%df_3body(:,is))))then
          write(err_msg,'("NaN in 3-body distribution function for ",A,&
               &" with norm of ",F0.3)' &
          ) &
               this%element_info(is)%name, this%norm_3body(is)
          call stop_program( err_msg )
          return
       elseif(any(isnan(this%gdf%df_4body(:,is))))then
          write(err_msg,'("NaN in 4-body distribution function for ",A,&
               &" with norm of ",F0.3)' &
          ) &
               this%element_info(is)%name, this%norm_4body(is)
          call stop_program( err_msg )
          return
       end if
    end do

    this%num_evaluated_allocated = size(this%system)
    this%num_evaluated = this%num_evaluated + num_evaluated

    this%viability_3body_default = sum( this%gdf%df_3body ) / &
         real( size( this%gdf%df_3body ), real32 )
    this%viability_4body_default = sum( this%gdf%df_4body ) / &
         real( size( this%gdf%df_4body ), real32 )

  end subroutine evolve
!###############################################################################

end module raffle__distribs_container
