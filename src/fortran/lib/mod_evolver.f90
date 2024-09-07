module evolver
  !! Module for handling distribution functions.
  !!
  !! This module contains the types and subroutines for handling distribution
  !! distributions. The distribution functions are used as fingerprints for
  !! atomic structures to identify similarities and differences between
  !! structures.
  use constants, only: real12, pi
  use misc_raffle, only: set, icount, strip_null, sort_str
  use misc_maths, only: triangular_number, set_difference
  use misc_linalg, only: get_angle, get_vol, get_improper_dihedral_angle, &
       cross, modu
  use rw_geom, only: basis_type, get_element_properties
  use extended_geom, only: extended_basis_type
  use elements, only: &
       element_type, element_bond_type, &
       element_database, element_bond_database
  implicit none

  
  private

  public :: gvector_container_type, gvector_base_type, gvector_type
  
  
  type :: gvector_base_type
     !! Base type for distribution functions.
     real(real12), dimension(:,:), allocatable :: df_2body
     !! 2-body distribution function.
     real(real12), dimension(:,:), allocatable :: df_3body
     !! 3-body distribution function.
     real(real12), dimension(:,:), allocatable :: df_4body
     !! 4-body distribution function.
  end type gvector_base_type

  type, extends(gvector_base_type) :: gvector_type
     !! Type for distribution functions.
     !!
     !! This type contains the distribution functions for a single atomic
     !! structure. It also contains other structure properties, including:
     !! - energy
     !! - stoichiometry
     !! - elements
     !! - number of atoms
     integer :: num_atoms = 0
     !! Number of atoms in the structure.
     real(real12) :: energy = 0.0_real12
     !! Energy of the structure.
     integer, dimension(:), allocatable :: stoichiometry
     !! Stoichiometry of the structure.
     character(len=3), dimension(:), allocatable :: element_symbols
     !! Elements contained within the structure.
   contains
     procedure, pass(this) :: calculate
  end type gvector_type


  !! should gvector be for the entire prediction, or one for each system?
  !! if one for each system, do not contain nbins, width, as this would ...
  !! ... result in lots of duplicated data
  type :: gvector_container_type
     !! Container for distribution functions.
     !!
     !! This type contains the distribution functions for a set of atomic
     !! structures, alongside parameters for initialising the distributions.
     integer :: num_evaluated = 0
     !! Number of evaluated systems.
     integer :: num_evaluated_allocated = 0
     !! Number of evaluated systems still allocated.
     integer :: best_system = 0
     !! Index of the best system.
     real(real12) :: best_energy = 0.0_real12
     !! Energy of the best system.
     real(real12) :: kbt = 0.2_real12
     !! Boltzmann constant times temperature.
     real(real12) :: &
          viability_3body_default = 0.1_real12, &
          viability_4body_default = 0.1_real12
     !! Default viability for the 3- and 4-body distribution functions.
     integer, dimension(3) :: nbins = -1
     !! Number of bins for the 2-, 3-, and 4-body distribution functions.
     real(real12), dimension(3) :: &
          sigma = [ 0.1_real12, 0.1_real12, 0.1_real12 ]
     !! Sigma of the gaussians used in the 2-, 3-, and 4-body 
     !! distribution functions.
     real(real12), dimension(3) :: &
          width = [ 0.025_real12, pi/64._real12, pi/64._real12 ]
     !! Width of the bins used in the 2-, 3-, and 4-body distribution functions.
     real(real12), dimension(3) :: &
          cutoff_min = [ 0.5_real12, 0._real12, 0._real12 ]
     !! Minimum cutoff for the 2-, 3-, and 4-body distribution functions.
     real(real12), dimension(3) :: &
          cutoff_max = [ 6._real12, pi, pi ]
     !! Maximum cutoff for the 2-, 3-, and 4-body distribution functions.
     real(real12), dimension(4) :: &
          radius_distance_tol = [ 1.5_real12, 2.5_real12, 3._real12, 6._real12 ]
     !! Tolerance for the distance between atoms for 3- and 4-body.
     !! index 1 = lower bound for 3-body
     !! index 2 = upper bound for 3-body
     !! index 3 = lower bound for 4-body
     !! index 4 = upper bound for 4-body
     real(real12), dimension(:), allocatable :: &
          norm_2body, norm_3body, norm_4body
     !! Normalisation factors for the 2-, 3-, and 4-body distribution functions.
     type(gvector_base_type) :: total !! name it best instead?
     !! Total distribution functions for all systems.
     !! Generated from combining the energy-weighted distribution functions
     !! of all systems
     type(gvector_type), dimension(:), allocatable :: system
     !! Distribution functions for each system.
     type(element_type), dimension(:), allocatable :: element_info
     !! Information about the elements in the container.
     type(element_bond_type), dimension(:), allocatable :: bond_info
     !! Information about the 2-body bonds in the container.

     !! @note
     !! Defaults for distribution function parametsr are randomly chosen for now.
     !! @endnote
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
     procedure, pass(this) :: initialise_gvectors
     !! Initialise the distribution functions in the container.
     procedure, pass(this) :: evolve
     !! Evolve the learned distribution function.
     procedure, pass(this) :: write
     !! Write all distribution functions to a file.
     procedure, pass(this) :: read
     !! Read all distribution functions from a file.
     procedure, pass(this) :: write_2body
     !! Write the learned 2-body distribution function to a file.
     procedure, pass(this) :: write_3body
     !! Write the learned 3-body distribution function to a file.
     procedure, pass(this) :: write_4body
     !! Write the learned 4-body distribution function to a file.
     procedure, pass(this) :: get_pair_index
     !! Return the index for bond_info given two elements.
     procedure, pass(this) :: get_bin
     !! Return the bin index for a given distance.
  end type gvector_container_type

  interface gvector_container_type
    !! Interface for the distribution functions container.
    module function init_gvector_container( &
         nbins, width, sigma, cutoff_min, cutoff_max &
         ) result(gvector_container)
         !! Initialise the distribution functions container.
         integer, dimension(3), intent(in), optional :: nbins
         !! Optional. Number of bins for the 2-, 3-, and 4-body distribution
         !! functions.
         real(real12), dimension(3), intent(in), optional :: width, sigma
         !! Optional. Width and sigma of the gaussians used in the 2-, 3-, and 
         !! 4-body.
         real(real12), dimension(3), intent(in), optional :: &
              cutoff_min, cutoff_max
         !! Optional. Minimum and maximum cutoff for the 2-, 3-, and 4-body.
         type(gvector_container_type) :: gvector_container
         !! Instance of the distribution functions container.
    end function init_gvector_container
  end interface gvector_container_type


  contains
  
!###############################################################################
  module function init_gvector_container( &
       nbins, width, sigma, &
       cutoff_min, cutoff_max ) &
       result(gvector_container)
    !! Initialise the distribution functions container.
    implicit none

    ! Arguments
    integer, dimension(3), intent(in), optional :: nbins
    !! Optional. Number of bins for the 2-, 3-, and 4-body distribution 
    !! functions.
    real(real12), dimension(3), intent(in), optional :: width, sigma
    !! Optional. Width and sigma of the gaussians used in the 2-, 3-, and 
    !! 4-body.
    real(real12), dimension(3), intent(in), optional :: cutoff_min, cutoff_max
    !! Optional. Minimum and maximum cutoff for the 2-, 3-, and 4-body.
    type(gvector_container_type) :: gvector_container
    !! Instance of the distribution functions container.


    if(present(nbins))then
       if(all(nbins .gt. 0)) gvector_container%nbins = nbins
    end if

    if(present(width))then
       if(all(width.ge.0._real12)) gvector_container%width = width
    end if

    if(present(sigma))then
       if(all(sigma.ge.0._real12)) gvector_container%sigma = sigma
    end if

    if(present(cutoff_min))then
       if(any(cutoff_min.ge.0._real12)) &
            gvector_container%cutoff_min = cutoff_min
    end if
    if(present(cutoff_max))then
       if(all(cutoff_max.ge.0._real12)) &
            gvector_container%cutoff_max = cutoff_max
    end if
    if(any(gvector_container%cutoff_max .le. gvector_container%cutoff_min))then
       write(0,*) "ERROR: cutoff_max <= cutoff_min"
       write(0,*) "cutoff min: ", gvector_container%cutoff_min
       write(0,*) "cutoff max: ", gvector_container%cutoff_max
       stop 1
    end if

  end function init_gvector_container
!###############################################################################


!###############################################################################
  subroutine set_width(this, width)
    !! Set the width of the gaussians used in the 2-, 3-, and 4-body
    !! distribution functions.
    implicit none

    ! Arguments
    class(gvector_container_type), intent(inout) :: this
    !! Parent. Instance of distribution functions container.
    real(real12), dimension(3), intent(in) :: width
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
    class(gvector_container_type), intent(inout) :: this
    !! Parent. Instance of distribution functions container.
    real(real12), dimension(3), intent(in) :: sigma
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
    class(gvector_container_type), intent(inout) :: this
    !! Parent. Instance of distribution functions container.
    real(real12), dimension(3), intent(in) :: cutoff_min
    !! Minimum cutoff for the 2-, 3-, and 4-body distribution functions.

     this%cutoff_min = cutoff_min

  end subroutine set_cutoff_min
!###############################################################################


!###############################################################################
  subroutine set_cutoff_max(this, cutoff_max)
    !! Set the maximum cutoff for the 2-, 3-, and 4-body distribution functions.
    implicit none
   
    ! Arguments
    class(gvector_container_type), intent(inout) :: this
    !! Parent. Instance of distribution functions container.
    real(real12), dimension(3), intent(in) :: cutoff_max
    !! Maximum cutoff for the 2-, 3-, and 4-body distribution functions.
   
    this%cutoff_max = cutoff_max
   
  end subroutine set_cutoff_max
!###############################################################################


!###############################################################################
  subroutine set_radius_distance_tol(this, radius_distance_tol)
    !! Set the tolerance for the distance between atoms for 3- and 4-body.
    implicit none

    ! Arguments
    class(gvector_container_type), intent(inout) :: this
    !! Parent. Instance of distribution functions container.
    real(real12), dimension(4), intent(in) :: radius_distance_tol
    !! Tolerance for the distance between atoms for 3- and 4-body.

    this%radius_distance_tol = radius_distance_tol

  end subroutine set_radius_distance_tol
!###############################################################################


!###############################################################################
  subroutine create(this, basis_list, deallocate_systems)
    !! create the distribution functions from the input file
    implicit none
    ! Arguments
    class(gvector_container_type), intent(inout) :: this
    !! Parent. Instance of distribution functions container.
    type(basis_type), dimension(:), intent(in) :: basis_list
    !! List of basis structures.
    logical, intent(in), optional :: deallocate_systems
    !! Optional. Boolean whether to deallocate the systems after the 
    !! distribution functions are created.

    ! Local variables
    logical :: deallocate_systems_
    
    
    if(.not.allocated(element_database))then
       write(0,*) "ERROR: element_database not allocated"
       write(0,*) "Run the set_element_energies() procedure of &
            &gvector_container_type before calling create()"
       stop 1
    end if

    deallocate_systems_ = .true.
    if(present(deallocate_systems)) deallocate_systems_ = deallocate_systems

    this%num_evaluated = 0
    this%num_evaluated_allocated = 0
    if(allocated(this%total%df_2body)) deallocate(this%total%df_2body)
    if(allocated(this%total%df_3body)) deallocate(this%total%df_3body)
    if(allocated(this%total%df_4body)) deallocate(this%total%df_4body)
    if(allocated(this%norm_2body)) deallocate(this%norm_2body)
    if(allocated(this%norm_3body)) deallocate(this%norm_3body)
    if(allocated(this%norm_4body)) deallocate(this%norm_4body)
    if(allocated(this%system)) deallocate(this%system)
    allocate(this%system(0))
    call this%add(basis_list)
    call this%set_bond_info()
    call this%evolve()
    if(deallocate_systems_) call this%deallocate_systems()
    
  end subroutine create
!###############################################################################


!###############################################################################
  subroutine update(this, basis_list, deallocate_systems)
    !! update the distribution functions from the input file
    implicit none
    ! Arguments
    class(gvector_container_type), intent(inout) :: this
    !! Parent. Instance of distribution functions container.
    type(basis_type), dimension(:), intent(in) :: basis_list
    !! List of basis structures.
    logical, intent(in), optional :: deallocate_systems
    !! Optional. Boolean whether to deallocate the systems after the
    !! distribution functions are created.

    ! Local variables
    logical :: deallocate_systems_


    deallocate_systems_ = .true.
    if(present(deallocate_systems)) deallocate_systems_ = deallocate_systems

    call this%add(basis_list)
    call this%update_bond_info()
    call this%evolve()
    if(deallocate_systems_) call this%deallocate_systems()
    
  end subroutine update
!###############################################################################

!###############################################################################
  subroutine deallocate_systems(this)
    !! Deallocate the systems in the container.
    implicit none

    ! Arguments
    class(gvector_container_type), intent(inout) :: this
    !! Parent. Instance of distribution functions container.

    deallocate(this%system)
    this%best_system = 0
    this%num_evaluated_allocated = 0

  end subroutine deallocate_systems
!###############################################################################


!###############################################################################
  subroutine write(this, file)
    !! Write all distribution functions for each system to a file.
    implicit none

    ! Arguments
    class(gvector_container_type), intent(in) :: this
    !! Parent. Instance of distribution functions container.
    character(*), intent(in) :: file
    !! Filename to write the distribution functions to.

    ! Local variables
    integer :: unit
    !! File unit.
    integer :: i, j
    !! Loop indices.


    if(.not.allocated(this%system))then
       write(0,*) "ERROR: No systems to write"
       write(0,*) "Systems either not created or deallocated after evolve"
       write(0,*) "To stop automatic deallocation, &
            &use the following flag in create()"
       write(0,*)
       write(0,*) "   deallocate_systems = .false."
       write(0,*)
       stop 1
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

  end subroutine write
!###############################################################################


!###############################################################################
  subroutine read(this, file)
    !! Read all distribution functions for each system from a file.
    implicit none

    ! Arguments
    class(gvector_container_type), intent(inout) :: this
    !! Parent. Instance of distribution functions container.
    character(*), intent(in) :: file
    !! Filename to read the distribution functions from.

    ! Local variables
    integer :: unit
    !! File unit.
    integer :: i, j
    !! Loop indices.
    integer :: iostat
    !! I/O status.
    integer :: num_species, num_pairs
    !! Number of species and pairs.
    character(256) :: buffer
    !! Buffer for reading lines.
    type(gvector_type) :: system
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
       num_pairs = nint( gamma(real(num_species + 2, real12)) / &
                   ( gamma(real(num_species, real12)) * gamma( 3._real12 ) ) )
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
       deallocate(system%element_symbols, system%stoichiometry, &
                  system%df_2body, system%df_3body, system%df_4body)
    end do
    close(unit)

  end subroutine read
!###############################################################################


!###############################################################################
  subroutine write_2body(this, file)
    !! Write the learned 2-body distribution functions to a file.
    implicit none

    ! Arguments
    class(gvector_container_type), intent(in) :: this
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


    num_pairs = nint( gamma(real(size(this%element_info) + 2, real12)) / &
                ( gamma(real(size(this%element_info), real12)) * &
                  gamma( 3._real12 ) ) )
    allocate(idx(2,num_pairs))
    i = 0 
    do is = 1, size(this%element_info)
       do js = is, size(this%element_info), 1
          i = i + 1
          idx(:,i) = [is, js]
       end do
    end do

    open(newunit=unit, file=file)
    do i = 1,  size(this%total%df_2body, dim=2)
       write(unit,'("# ",A,2X,A)') &
            this%element_info(idx(1,i))%name, &
            this%element_info(idx(2,i))%name
       do j = 1, size(this%total%df_2body, dim=1)
          write(unit,*) this%cutoff_min(1) + this%width(1) * ( j - 1 ), &
                        this%total%df_2body(j,i)
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
    class(gvector_container_type), intent(in) :: this
    !! Parent. Instance of distribution functions container.
    character(*), intent(in) :: file
    !! Filename to write the 3-body distribution functions to.

    ! Local variables
    integer :: unit
    !! File unit.
    integer :: i, j
    !! Loop indices.


    open(newunit=unit, file=file)
    do i = 1,  size(this%total%df_3body, dim=2)
       write(unit,'("# ",A)') this%element_info(i)%name
       do j = 1, size(this%total%df_3body, dim=1)
          write(unit,*) this%cutoff_min(2) + this%width(2) * ( j - 1 ), &
                        this%total%df_3body(j,i)
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
    class(gvector_container_type), intent(in) :: this
    !! Parent. Instance of distribution functions container.
    character(*), intent(in) :: file
    !! Filename to write the 4-body distribution functions to.

    ! Local variables
    integer :: unit
    !! File unit.
    integer :: i, j
    !! Loop indices.


    open(newunit=unit, file=file)
    do i = 1,  size(this%total%df_4body, dim=2)
       write(unit,'("# ",A)') this%element_info(i)%name
       do j = 1, size(this%total%df_4body, dim=1)
          write(unit,*) this%cutoff_min(3) + this%width(3) * ( j - 1 ), &
                        this%total%df_4body(j,i)
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
    class(gvector_container_type), intent(inout) :: this
    !! Parent. Instance of distribution functions container.
    class(*), dimension(..), intent(in) :: system
    !! System to add to the container.

    ! Local variables
    integer :: i
    !! Loop index.
    integer :: num_structures_previous
    !! Number of structures in the container before adding the system.
    character(128) :: buffer
    !! Buffer for writing messages.


    select rank(system)
    rank(0)
       select type(system)
       type is (gvector_type)
          this%system = [ this%system, system ]
       type is (basis_type)
          call this%add_basis(system)
       class default
          write(0,*) "ERROR: Invalid type for system"
          write(0,*) "Expected type gvector_type or basis_type"
          stop 1
       end select
    rank(1)
       num_structures_previous = size(this%system)
       select type(system)
       type is (gvector_type)
          this%system = [ this%system, system ]
       type is (basis_type)
          do i = 1, size(system)
             call this%add_basis(system(i))
          end do
       class default
          write(0,*) "ERROR: Invalid type for system"
          write(0,*) "Expected type gvector_type or basis_type"
          stop 1
       end select
    rank default
       write(0,*) "ERROR: Invalid rank for system"
       write(buffer,*) rank(system)
       write(0,*) "Expected rank 0 or 1, got ", trim(buffer)
       stop 1
    end select
    call this%update_element_info()
    call this%update_bond_info()
    call this%set_best_energy()

  end subroutine add
!###############################################################################


!###############################################################################
  subroutine add_basis(this, basis)
    !! Add a basis to the container.
    implicit none

    ! Arguments
    class(gvector_container_type), intent(inout) :: this
    !! Parent. Instance of distribution functions container.
    type(basis_type), intent(in) :: basis
    !! Basis to add to the container.

    ! Local variables
    type(gvector_type) :: system
    !! System to add to the container.

    call system%calculate(basis, width = this%width, &
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
  subroutine set_element_info(this)
    !! Set the list of elements for the container.
    implicit none

    ! Arguments
    class(gvector_container_type), intent(inout) :: this
    !! Parent of the procedure. Instance of distribution functions container.

    ! Local variables
    integer :: i
    !! Loop index.
    real(real12) :: radius
    !! Element radii.
    character(len=3), dimension(:), allocatable :: element_list
    !! List of elements in the container.


    !---------------------------------------------------------------------------
    ! get list of species in dataset
    !---------------------------------------------------------------------------
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
    class(gvector_container_type), intent(inout) :: this
    !! Parent of the procedure. Instance of distribution functions container.

    ! Local variables
    integer :: i
    !! Loop index.
    real(real12) :: radius
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
    class(gvector_container_type), intent(inout) :: this
    !! Parent of the procedure. Instance of distribution functions container.
    character(len=3), intent(in) :: element
    !! Element name.
    real(real12), intent(in) :: energy
    !! Energy of the element.

    ! Local variables
    integer :: idx, idx_db
    !! Index of the element in the element_info array.
    real(real12) :: radius
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
    class(gvector_container_type), intent(inout) :: this
    !! Parent of the procedure. Instance of distribution functions container.
    character(len=3), dimension(:), intent(in) :: elements
    !! Element names.
    real(real12), dimension(:), intent(in) :: energies
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
   class(gvector_container_type), intent(in) :: this
   !! Parent of the procedure. Instance of distribution functions container.
   character(len=3), dimension(:), allocatable, intent(out) :: elements
   !! Element names.
   real(real12), dimension(:), allocatable, intent(out) :: energies
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
    class(gvector_container_type), intent(in) :: this
    !! Parent of the procedure. Instance of distribution functions container.
    character(len=3), dimension(size(this%element_info,1)), intent(out) :: &
         elements
    !! Element names.
    real(real12), dimension(size(this%element_info,1)), intent(out) :: energies
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
    class(gvector_container_type), intent(inout) :: this
    !! Parent of the procedure. Instance of distribution functions container.

    ! Local variables
    integer :: i, j, k, idx1, idx2
    !! Loop index.
    integer :: num_elements, num_pairs
    !! Number of elements and pairs.
    real(real12) :: radius
    !! Average of covalent radii.
    logical :: success
    !! Success flag.


    !---------------------------------------------------------------------------
    ! allocate the bond information array
    !---------------------------------------------------------------------------
    num_elements = size(this%element_info)
    num_pairs = nint(gamma(real(num_elements + 2, real12)) / &
         ( gamma(real(num_elements, real12)) * gamma( 3._real12 ) ) )
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
               this%element_info(j)%name ] &
          )
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
    real(real12) :: radius
    !! Average of covalent radii.


    write(0,*) 'WARNING: No bond data for element pair ', &
               elements(1), ' and ', &
               elements(2)
    write(0,*) 'WARNING: Setting bond to average of covalent radii'
    if(.not.allocated(element_database))then
       write(0,*) "ERROR: Element database not initialised"
       stop 1
    end if
    idx1 = findloc([ element_database(:)%name ], &
         elements(1), dim=1)
    idx2 = findloc([ element_database(:)%name ], &
         elements(2), dim=1)
    if(idx1.lt.1.or.idx2.lt.1)then
       write(0,*) "ERROR: Element not found in database"
       if(idx1.lt.1) write(0,*) "Element: ", elements(1)
       if(idx2.lt.1) write(0,*) "Element: ", elements(2)
       write(0,*) "Indices: ", idx1, idx2
       stop 1
    end if
    radius = ( element_database(idx1)%radius + &
         element_database(idx2)%radius ) / 2._real12
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
   class(gvector_container_type), intent(inout) :: this
   !! Parent of the procedure. Instance of distribution functions container.

   ! Local variables
   integer :: i, j, k, is, js
   !! Loop index.
   real(real12) :: radius, radius1, radius2
   !! Average of covalent radii.
   character(len=3), dimension(:), allocatable :: element_list
   !! List of elements in the container.
   character(len=3), dimension(:,:), allocatable :: pair_list
   !! List of element pairs in the container.


   !----------------------------------------------------------------------------
   ! check if bond_info is allocated, if not, set it and return
   !----------------------------------------------------------------------------
   if(.not.allocated(this%bond_info))then
      call this%set_bond_info()
      return
   elseif(size(this%bond_info).eq.0)then
      call this%set_bond_info()
      return
   end if


   !----------------------------------------------------------------------------
   ! get list of element pairs in dataset
   !----------------------------------------------------------------------------
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

   !----------------------------------------------------------------------------
   ! check if all element pairs are in the database
   !----------------------------------------------------------------------------
   if(.not.allocated(element_bond_database)) allocate(element_bond_database(0))
   pair_loop1: do i = 1, size(pair_list,1)
      do j = 1, size(element_bond_database)
         if( ( element_bond_database(j)%element(1) .eq. pair_list(i,1) .and. &
               element_bond_database(j)%element(2) .eq. pair_list(i,2) ) .or. &
             ( element_bond_database(j)%element(1) .eq. pair_list(i,2) .and. &
               element_bond_database(j)%element(2) .eq. pair_list(i,1) ) &
              ) cycle pair_loop1
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
      radius = ( radius1 + radius2 ) / 2._real12
      element_bond_database = [ element_bond_database, &
           element_bond_type(elements=[pair_list(i,:)], radius=radius) ]
      call sort_str(element_bond_database(size(element_bond_database))%element)
   end do pair_loop1


   ! ---------------------------------------------------------------------------
   ! check if all element pairs are in the bond_info array
   ! ---------------------------------------------------------------------------
   pair_loop2: do i = 1, size(pair_list,1)
       info_loop: do j = 1, size(this%bond_info)
          if( ( this%bond_info(j)%element(1) .eq. pair_list(i,1) .and. &
                this%bond_info(j)%element(2) .eq. pair_list(i,2) ) .or. &
              ( this%bond_info(j)%element(1) .eq. pair_list(i,2) .and. &
                this%bond_info(j)%element(2) .eq. pair_list(i,1) ) &
               )then
             cycle pair_loop2
          end if
       end do info_loop
       this%bond_info = [ this%bond_info(:), &
                          element_bond_type([ pair_list(i,1:2) ]) ]
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
    class(gvector_container_type), intent(inout) :: this
    !! Parent of the procedure. Instance of distribution functions container.
    character(len=3), dimension(2), intent(in) :: elements
    !! Element name.
    real(real12), intent(in) :: radius
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
       if( ( element_bond_database(i)%element(1) .eq. element_1 .and. &
             element_bond_database(i)%element(2) .eq. element_2 ) .or. &
           ( element_bond_database(i)%element(1) .eq. element_2 .and. &
             element_bond_database(i)%element(2) .eq. element_1 ) &
            )then
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
    class(gvector_container_type), intent(inout) :: this
    !! Parent of the procedure. Instance of distribution functions container.
    character(len=3), dimension(:,:), intent(in) :: elements
    !! Element names.
    real(real12), dimension(:), intent(in) :: radii
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
   class(gvector_container_type), intent(in) :: this
   !! Parent of the procedure. Instance of distribution functions container.
   character(len=3), dimension(:,:), allocatable, intent(out) :: elements
   !! Element pair names.
   real(real12), dimension(:), allocatable, intent(out) :: radii
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
    class(gvector_container_type), intent(in) :: this
    !! Parent of the procedure. Instance of distribution functions container.
    character(len=3), dimension(size(this%bond_info,1),2), intent(out) :: &
         elements
    !! Element pair names.
    real(real12), dimension(size(this%bond_info,1)), intent(out) :: radii
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
    class(gvector_container_type), intent(inout) :: this
    !! Parent of the procedure. Instance of distribution functions container.

    ! Local variables
    integer :: i, is, idx
    !! Loop index.
    real(real12) :: energy
    !! Energy of the system.

    do i = 1, size(this%system)
       energy = this%system(i)%energy
       do is = 1, size(this%system(i)%element_symbols)
          idx = findloc( [ this%element_info(:)%name ], &
                           this%system(i)%element_symbols(is), dim=1 )
          if(idx.lt.1)then
             write(0,*) "ERROR: Species not found in element_info"
             stop 1
          end if
          energy = energy - this%system(i)%stoichiometry(is) * &
                            this%element_info(idx)%energy
       end do
       energy = energy / this%system(i)%num_atoms
       if( energy .lt. this%best_energy ) then
          this%best_energy = energy
          this%best_system = i
       end if
    end do

  end subroutine set_best_energy
!###############################################################################


!###############################################################################
  pure function get_pair_index(this, species1, species2) result(idx)
    !! Get the index of a pair of elements in the container.
    implicit none

    ! Arguments
    class(gvector_container_type), intent(in) :: this
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
    !idx = nint( ( size(this%element_info) - min( is, js ) / 2._real12 ) * &
    !      ( is - 1._real12 ) + js )
    idx = nint( ( size(this%element_info) - min( is, js ) / 2._real12 ) * &
          ( min( is, js ) - 1._real12 ) + max( is, js ) )

  end function get_pair_index
!###############################################################################


!###############################################################################
  pure function get_bin(this, value, dim) result(bin)
    !! Get the bin index for a value in a dimension.
    implicit none

    ! Arguments
    class(gvector_container_type), intent(in) :: this
    !! Parent of the procedure. Instance of distribution functions container.
    real(real12), intent(in) :: value
    !! Value to get the bin index for.
    integer, intent(in) :: dim
    !! Dimension to get the bin index for.
    integer :: bin
    !! Bin index for the value.

    if(value .lt. this%cutoff_min(dim) - this%width(dim) .or. &
         value .gt. this%cutoff_max(dim) + this%width(dim))then
       bin = 0
    else
       bin = nint( ( this%nbins(dim) - 1 ) * &
                   ( value - this%cutoff_min(dim) ) / &
                   ( this%cutoff_max(dim) - this%cutoff_min(dim) ) ) + 1
    end if

  end function get_bin
!###############################################################################


!###############################################################################
  subroutine initialise_gvectors(this)
    !! Initialise the g-vectors for the container.
    implicit none

    ! Arguments
    class(gvector_container_type), intent(inout) :: this
    !! Parent of the procedure. Instance of distribution functions container.

    ! Local variables
    integer :: i
    !! Loop index.
    integer :: num_pairs
    !! Number of pairs.
    real(real12) :: eta, weight, height
    !! Parameters for the g-vectors.
    !real(real12), dimension(42) :: bonds_cubic

    num_pairs = nint( gamma(real(size(this%element_info) + 2, real12)) / &
         ( gamma(real(size(this%element_info), real12)) * gamma( 3._real12 ) ) )
    allocate(this%total%df_2body(this%nbins(1),num_pairs), &
         source = 0._real12 )
    allocate(this%total%df_3body(this%nbins(2),size(this%element_info)), &
         source = 0._real12 )
    allocate(this%total%df_4body(this%nbins(3),size(this%element_info)), &
         source = 0._real12 )
   !  allocate(this%total%df_3body(this%nbins(2),size(this%element_info)), &
   !       source = 1._real12/this%nbins(2))
   !  allocate(this%total%df_4body(this%nbins(3),size(this%element_info)), &
   !       source = 1._real12/this%nbins(3))

    !  this%total%df_2body(:,:) = 1._real12 / this%nbins(1)
    !! make it extra broad
    !bonds_cubic(:6) = 1._real12
    !bonds_cubic(7:14) = sqrt(2._real12)
    !bonds_cubic(15:22) = sqrt(3._real12)
    !bonds_cubic(23:30) = 2._real12
    !bonds_cubic(31:42) = sqrt(5._real12)
    !weight = exp( this%best_energy )
    !height = 1._real12 / this%nbins(1)
    !eta = 1._real12 / ( 2._real12 * ( this%sigma(1) )**2._real12 )
    !do i = 1, num_pairs
    !   this%total%df_2body(:,i) = weight * height * get_gvector( &
    !                      bonds_cubic * this%bond_info(i)%radius_covalent , &
    !                      this%nbins(1), eta, this%width(1), &
    !                      this%cutoff_min(1), &
    !                      ( this%cutoff_max(1) - this%cutoff_min(1) ) &
    !   )
    !end do

  end subroutine initialise_gvectors
!###############################################################################


!###############################################################################
  subroutine evolve(this, system)
    !! Evolve the g-vectors for the container.
    implicit none

    ! Arguments
    class(gvector_container_type), intent(inout) :: this
    !! Parent of the procedure. Instance of distribution functions container.
    type(gvector_type), dimension(..), intent(in), optional :: system
    !! Optional. System to add to the container.

    ! Local variables
    integer :: i, j, is, js
    !! Loop index.
    integer :: idx1, idx2
    !! Index of the element in the element_info array.
    integer :: num_evaluated
    !! Number of systems evaluated this iteration.
    real(real12) :: weight, energy, best_energy_old
    !! Energy and weight variables for a system.
    real(real12), dimension(:), allocatable :: height
    !! Height of the g-vectors.
    integer, dimension(:,:), allocatable :: idx_list
    !! Index list for the element pairs in a system.


    !---------------------------------------------------------------------------
    ! if present, add the system to the container
    !---------------------------------------------------------------------------
    if(present(system)) call this%add(system)
    do i = 1, 3
       if(this%nbins(i).eq.-1)then
          this%nbins(i) = 1 + &
               ( this%cutoff_max(i) - this%cutoff_min(i) ) / this%width(i)
       end if
    end do


    !---------------------------------------------------------------------------
    ! get the energy from the lowest formation energy system
    !---------------------------------------------------------------------------
    best_energy_old = this%best_energy
    call this%set_best_energy()


    !---------------------------------------------------------------------------
    ! initialise the total gvectors
    !---------------------------------------------------------------------------
    if(.not.allocated(this%total%df_2body))then
       call this%initialise_gvectors()
    else
      this%total%df_2body = this%total%df_2body * &
                              exp( this%best_energy / this%kbt ) / &
                              exp( best_energy_old / this%kbt )
      this%total%df_3body = this%total%df_3body * &
                              exp( this%best_energy / this%kbt ) / &
                              exp( best_energy_old / this%kbt )
      this%total%df_4body = this%total%df_4body * &
                              exp( this%best_energy / this%kbt ) / &
                              exp( best_energy_old / this%kbt )
      do j = 1, size(this%total%df_2body,2)
         this%total%df_2body(:,j) = &
              this%total%df_2body(:,j) * this%norm_2body(j)
      end do
      do is = 1, size(this%element_info)
         this%total%df_3body(:,is) = &
              this%total%df_3body(:,is) * this%norm_3body(is)
         this%total%df_4body(:,is) = &
              this%total%df_4body(:,is) * this%norm_4body(is)
      end do
      deallocate(this%norm_2body)
      deallocate(this%norm_3body)
      deallocate(this%norm_4body)
    end if


    !---------------------------------------------------------------------------
    ! loop over all systems to calculate the total gvectors
    !---------------------------------------------------------------------------
    num_evaluated = 0
    do i = this%num_evaluated_allocated + 1, size(this%system), 1
       num_evaluated = num_evaluated + 1
       !------------------------------------------------------------------------
       ! get the list of 2-body species pairs the system
       !------------------------------------------------------------------------
       j = 0
       allocate(idx_list(size(this%system(i)%element_symbols),&
                         size(this%system(i)%element_symbols)))
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
          idx1 = findloc( [ this%element_info(:)%name ], &
                          this%system(i)%element_symbols(is), dim=1)
          if(idx1.lt.1)then
             write(0,*) "ERROR: Species not found in species list"
             stop 1
          end if
          energy = energy - this%system(i)%stoichiometry(is) * &
               this%element_info(idx1)%energy
       end do
       energy = energy / this%system(i)%num_atoms
       weight = exp( ( this%best_energy - energy ) / this%kbt )
       j = 0

       !------------------------------------------------------------------------
       ! loop over all species in the system to add the gvectors
       !------------------------------------------------------------------------
       do is = 1, size(this%system(i)%element_symbols)

          idx1 = findloc( [ this%element_info(:)%name ], &
                          this%system(i)%element_symbols(is), dim=1)

          ! height = 1._real12 / ( 1._real12 + this%total%df_3body(:,idx1) )
          this%total%df_3body(:,idx1) = this%total%df_3body(:,idx1) + &
               ! weight * &
               set_difference( weight * this%system(i)%df_3body(:,is), &
                               this%total%df_3body(:,idx1), & !/ max( &
                              !       1._real12, &
                              !       maxval(this%total%df_3body(:,idx1)) &
                              !  ), &
                               set_min_zero = .true. &
               )
               ! height * this%system(i)%df_4body(:,is)
          
          ! height = 1._real12 / ( 1._real12 + this%total%df_4body(:,idx1) )
          this%total%df_4body(:,idx1) = this%total%df_4body(:,idx1) + &
               ! weight * &
               set_difference( weight * this%system(i)%df_4body(:,is), &
                               this%total%df_4body(:,idx1), &! !/ max( &
                              !       1._real12, &
                              !       maxval(this%total%df_4body(:,idx1)) &
                              !  ), &
                               set_min_zero = .true. &
               )
               ! height * this%system(i)%df_4body(:,is)
          
          do js = is, size(this%system(i)%element_symbols), 1
             idx2 = findloc( [ this%element_info(:)%name ], &
                             this%system(i)%element_symbols(js), dim=1)
             j = nint( ( size(this%element_info) - &
                         min( idx1, idx2 ) / 2._real12 ) * &
                         ( min( idx1, idx2 ) - 1._real12 ) + max( idx1, idx2 ) )
             ! height = 1._real12 / &
             !      ( 1._real12 + this%total%df_2body(:,j) ) ** 2._real12
             this%total%df_2body(:,j) = this%total%df_2body(:,j) + &
                  ! weight * &
                  set_difference( weight * this%system(i)%df_2body(:,idx_list(is,js)), &
                                  this%total%df_2body(:,j), & !/ max( &
                                 !       1._real12, &
                                 !       maxval(this%total%df_2body(:,j)) &
                                 !  ), &
                                  set_min_zero = .true. &
                  )
                  ! height * this%system(i)%df_2body(:,idx_list(is,js))

          end do
       end do
       deallocate(idx_list)
   end do
   
   allocate(this%norm_2body(size(this%total%df_2body,2)))
   do j = 1, size(this%total%df_2body,2)
      this%norm_2body(j) = maxval(this%total%df_2body(:,j))
      ! this%norm_2body(j) = sum(this%total%df_2body(:,j))/size(this%total%df_2body,1)
      if(abs(this%norm_2body(j)).lt.1.E-6)then
         write(0,*) "ERROR: Zero norm for 2-body g-vector"
         stop 1
      end if
      this%total%df_2body(:,j) = &
           this%total%df_2body(:,j) / this%norm_2body(j)
   end do
   allocate(this%norm_3body(size(this%element_info)))
   allocate(this%norm_4body(size(this%element_info)))
   do is = 1, size(this%element_info)
      this%norm_3body(is) = maxval(this%total%df_3body(:,is))
      ! this%norm_3body(is) = sum(this%total%df_3body(:,is))/size(this%total%df_3body,1)
      if(abs(this%norm_3body(is)).lt.1.E-6)then
         write(0,*) "ERROR: Zero norm for 3-body g-vector"
         stop 1
      end if
      this%norm_4body(is) = maxval(this%total%df_4body(:,is))
      ! this%norm_4body(is) = sum(this%total%df_4body(:,is))/size(this%total%df_4body,1)
      if(abs(this%norm_4body(is)).lt.1.E-6)then
         write(0,*) "ERROR: Zero norm for 4-body g-vector"
         stop 1
      end if
      this%total%df_3body(:,is) = &
           this%total%df_3body(:,is) / this%norm_3body(is)
      this%total%df_4body(:,is) = &
           this%total%df_4body(:,is) / this%norm_4body(is)
   end do

   this%num_evaluated_allocated = size(this%system)
   this%num_evaluated = this%num_evaluated + num_evaluated

   this%viability_3body_default = sum(this%total%df_3body)/size(this%total%df_3body)
   this%viability_4body_default = sum(this%total%df_4body)/size(this%total%df_4body)

  end subroutine evolve
!###############################################################################


!###############################################################################
  subroutine calculate(this, basis, &
       nbins, width, sigma, cutoff_min, cutoff_max, radius_distance_tol)
    !! Calculate the distribution functions for the container.
    !!
    !! This procedure calculates the 2-, 3-, and 4-body distribution function 
    !! for a given atomic structure (i.e. basis).
    implicit none

    ! Arguments
    class(gvector_type), intent(inout) :: this
    !! Parent of the procedure. Instance of distribution functions container.
    type(basis_type), intent(in) :: basis
    !! Atomic structure.
    integer, dimension(3), intent(in), optional :: nbins
    !! Optional. Number of bins for the distribution functions.
    real(real12), dimension(3), intent(in), optional :: width, sigma
    !! Optional. Width and sigma for the distribution functions.
    real(real12), dimension(3), intent(in), optional :: cutoff_min, cutoff_max
    !! Optional. Cutoff minimum and maximum for the distribution functions.
    real(real12), dimension(4), intent(in), optional :: radius_distance_tol
    !! Tolerance for the distance between atoms for 3- and 4-body.

    ! Local variables
    integer, dimension(3) :: nbins_
    !! Number of bins for the distribution functions.
    real(real12), dimension(3) :: sigma_
    !! Sigma for the distribution functions.
    real(real12), dimension(3) :: width_
    !! Width of the bins for the distribution functions.
    real(real12), dimension(3) :: cutoff_min_
    !! Cutoff minimum for the distribution functions.
    real(real12), dimension(3) :: cutoff_max_
    !! Cutoff maximum for the distribution functions.
    type(element_bond_type), dimension(:), allocatable :: bond_info
    !! Bond information for radii.
    real(real12), dimension(4) :: radius_distance_tol_
    !! Tolerance for the distance between atoms for 3- and 4-body.


    !! @note
    !! Defaults for distribution function parametsr are randomly chosen for now.
    !! @endnote

    integer :: bin
    !! Bin index and maximum number of steps.
    integer :: i, j, b, itmp1, idx
    !! Loop index.
    integer :: is, js, ia, ja, ka, la
    !! Loop index.
    integer :: num_pairs!, num_angles
    !! Number of pairs and angles.
    real(real12) :: bondlength
    !! Temporary real variables.
    logical :: success
    !! Boolean for success.
    type(extended_basis_type) :: basis_extd
    !! Extended basis of the system.
    type(extended_basis_type) :: neighbour_basis
    !! Basis for storing neighbour data.
    real(real12), dimension(3) :: eta, limit
    !! Parameters for the distribution functions.
    real(real12), dimension(3) :: vtmp1, vtmp2, vtmp3, diff
    !! Temporary real arrays.
    real(real12), allocatable, dimension(:) :: angle_list, bondlength_list, &
         distance
    !! Temporary real arrays.
    integer, allocatable, dimension(:,:) :: pair_index
    !! Index of element pairs.


    !---------------------------------------------------------------------------
    ! initialise optional variables
    !---------------------------------------------------------------------------
    if(present(cutoff_min))then
       cutoff_min_ = cutoff_min
    else
       cutoff_min_ = [0.5_real12, 0._real12, 0._real12]
    end if
    if(present(cutoff_max))then
       cutoff_max_ = cutoff_max
    else
       cutoff_max_ = [6._real12, pi, pi]
    end if
    if(present(width))then
       width_ = width
    else
       width_ = [0.25_real12, pi/64._real12, pi/64._real12]
    end if
    if(present(sigma))then
       sigma_ = sigma
    else
       sigma_ = [0.1_real12, 0.1_real12, 0.1_real12]
    end if
    if(present(nbins))then
       nbins_ = nbins
       width_ = ( cutoff_max_ - cutoff_min_ )/real( nbins_ - 1, real12 )
    else
       nbins_ = 1 + nint( (cutoff_max_ - cutoff_min_)/width_ )
    end if
    limit = cutoff_max_ - cutoff_min_
    if(present(radius_distance_tol))then
       radius_distance_tol_ = radius_distance_tol
    else
       radius_distance_tol_ = [1.5_real12, 2.5_real12, 3._real12, 6._real12]
    end if
       


    !---------------------------------------------------------------------------
    ! get the number of pairs of species
    ! (this uses a combination calculator with repetition)
    !---------------------------------------------------------------------------
    num_pairs = gamma(real(basis%nspec + 2, real12)) / &
                ( gamma(real(basis%nspec, real12)) * gamma( 3._real12 ) )
    allocate(this%element_symbols(basis%nspec))
    do is = 1, basis%nspec
       this%element_symbols(is) = strip_null(basis%spec(is)%name)
    end do
    i = 0
    allocate(bond_info(num_pairs))
    allocate(pair_index(basis%nspec,basis%nspec))
    do is = 1, basis%nspec
       do js = is, basis%nspec, 1
          i = i + 1
          pair_index(js,is) = i
          pair_index(is,js) = i
          call bond_info(i)%set( this%element_symbols(is), &
                                 this%element_symbols(js), success &
          )
          if(success) cycle
          call set_bond_radius_to_default( [ &
               this%element_symbols(is), &
               this%element_symbols(js) ] &
          )
          call bond_info(i)%set( this%element_symbols(is), &
                                 this%element_symbols(js), success &
          )
       end do
    end do


    !---------------------------------------------------------------------------
    ! get the stoichiometry, energy, and number of atoms
    !---------------------------------------------------------------------------
    this%stoichiometry = basis%spec(:)%num
    this%energy = basis%energy
    this%num_atoms = basis%natom


    !---------------------------------------------------------------------------
    ! calculate the gaussian width and allocate the distribution functions
    !---------------------------------------------------------------------------
    eta = 1._real12 / ( 2._real12 * sigma_**2._real12 )
    allocate(this%df_2body(nbins_(1), num_pairs), source = 0._real12)
    allocate(this%df_3body(nbins_(2), basis%nspec), source = 0._real12)
    allocate(this%df_4body(nbins_(3), basis%nspec), source = 0._real12)


    !---------------------------------------------------------------------------
    ! create the extended basis and neighbour basis
    !---------------------------------------------------------------------------
    call basis_extd%copy(basis)
    call basis_extd%create_images( max_bondlength = cutoff_max_(1) )
    allocate(bondlength_list(basis_extd%natom+basis_extd%num_images))

    allocate(neighbour_basis%spec(1))
    allocate(neighbour_basis%image_spec(1))
    allocate(neighbour_basis%spec(1)%atom( &
         sum(basis_extd%spec(:)%num)+sum(basis_extd%image_spec(:)%num), 3 &
    ) )
    allocate(neighbour_basis%image_spec(1)%atom( &
         sum(basis_extd%spec(:)%num)+sum(basis_extd%image_spec(:)%num), 3 &
    ) )
    neighbour_basis%nspec = basis%nspec
    neighbour_basis%natom = 0
    neighbour_basis%num_images = 0
    neighbour_basis%lat = basis%lat


    !---------------------------------------------------------------------------
    ! calculate the distribution functions
    !---------------------------------------------------------------------------
    do is = 1, basis%nspec
       do ia = 1, basis%spec(is)%num
          allocate(distance(basis_extd%natom+basis_extd%num_images)) !!! ALLOCATE THIS ONCE AND JUST WRITE OVER ?
          neighbour_basis%spec(1)%num = 0
          neighbour_basis%image_spec(1)%num = 0
          do js = 1, basis%nspec
             itmp1 = 0

             !------------------------------------------------------------------
             ! loop over all atoms inside the unit cell
             !------------------------------------------------------------------
             atom_loop: do ja = 1, basis%spec(js)%num

                associate( vector =>  matmul( [ &
                          basis_extd%spec(js)%atom(ja,1:3) - &
                          basis_extd%spec(is)%atom(ia,1:3) &
                     ], basis%lat ) &
                )
                   bondlength = modu( vector )
                   
                   if( bondlength .lt. cutoff_min_(1) .or. &
                       bondlength .gt. cutoff_max_(1) ) cycle atom_loop
                  
                   ! add 2-body bond to store if within tolerances for 3-body
                   ! distance
                   if( &
                        bondlength .ge. &
                             bond_info(pair_index(is, js))%radius_covalent * &
                             radius_distance_tol(1) .and. &
                        bondlength .le. &
                             bond_info(pair_index(is, js))%radius_covalent * &
                             radius_distance_tol(2) &
                   ) then
                      neighbour_basis%spec(1)%num = &
                           neighbour_basis%spec(1)%num + 1
                      neighbour_basis%spec(1)%atom( &
                           neighbour_basis%spec(1)%num,1:3) = vector
                   end if

                   ! add 2-body bond to store if within tolerances for 4-body
                   ! distance
                   if( bondlength .ge. ( & 
                        bond_info(pair_index(is, js))%radius_covalent * &
                        radius_distance_tol(3) ) .and. &
                       bondlength .le. ( &
                        bond_info(pair_index(is, js))%radius_covalent * &
                        radius_distance_tol(4) ) &
                   ) then
                      neighbour_basis%image_spec(1)%num = &
                           neighbour_basis%image_spec(1)%num + 1
                      neighbour_basis%image_spec(1)%atom( &
                           neighbour_basis%image_spec(1)%num,1:3) = vector
                   end if

                   !if(js.lt.js.or.(is.eq.js.and.ja.le.ia)) cycle
                   itmp1 = itmp1 + 1
                   bondlength_list(itmp1) = bondlength
                   distance(itmp1) = 1._real12
                
                end associate
             end do atom_loop


             !------------------------------------------------------------------
             ! loop over all image atoms outside of the unit cell
             !------------------------------------------------------------------
             image_loop: do ja = 1, basis_extd%image_spec(js)%num
                associate( vector =>  matmul( [ &
                          basis_extd%image_spec(js)%atom(ja,1:3) - &
                          basis_extd%spec(is)%atom(ia,1:3) &
                     ], basis_extd%lat ) &
                )

                   bondlength = modu( vector )
                   
                   if( bondlength .lt. cutoff_min_(1) .or. &
                       bondlength .gt. cutoff_max_(1) ) cycle image_loop
                  
                   ! add 2-body bond to store if within tolerances for 3-body
                   ! distance
                   if( &
                        bondlength .ge. &
                             bond_info(pair_index(is, js))%radius_covalent * &
                             radius_distance_tol(1) .and. &
                        bondlength .le. &
                             bond_info(pair_index(is, js))%radius_covalent * &
                             radius_distance_tol(2) &
                   ) then
                      neighbour_basis%spec(1)%num = &
                           neighbour_basis%spec(1)%num + 1
                      neighbour_basis%spec(1)%atom( &
                           neighbour_basis%spec(1)%num,1:3 &
                      ) = vector
                   end if

                   ! add 2-body bond to store if within tolerances for 4-body
                   ! distance
                   if( bondlength .ge. ( & 
                        bond_info(pair_index(is, js))%radius_covalent * &
                        radius_distance_tol(3) ) .and. &
                       bondlength .le. ( &
                        bond_info(pair_index(is, js))%radius_covalent * &
                        radius_distance_tol(4) ) &
                   ) then
                      neighbour_basis%image_spec(1)%num = &
                           neighbour_basis%image_spec(1)%num + 1
                      neighbour_basis%image_spec(1)%atom( &
                           neighbour_basis%image_spec(1)%num,1:3 &
                      ) = vector
                   end if

                   itmp1 = itmp1 + 1
                   bondlength_list(itmp1) = bondlength
                   distance(itmp1) = 1._real12
                
                end associate
             end do image_loop

             !------------------------------------------------------------------
             ! calculate the 2-body distribution function contributions from
             ! atom (is,ia) for species pair (is,js)
             !------------------------------------------------------------------
             if(itmp1.gt.0)then
                this%df_2body(:,pair_index(is, js)) = &
                     this%df_2body(:,pair_index(is, js)) + &
                     get_gvector( &
                          bondlength_list(:itmp1), &
                          nbins_(1), eta(1), width_(1), &
                          cutoff_min_(1), &
                          limit(1), scale = distance(:itmp1) &
                     )
             end if

          end do
          deallocate(distance)


          !---------------------------------------------------------------------
          ! calculate the 3-body distribution function for atom (is,ia)
          !---------------------------------------------------------------------
          if(neighbour_basis%spec(1)%num.le.1) cycle
          associate( &
               num_angles => &
               triangular_number( neighbour_basis%spec(1)%num - 1 ) &
          )
             allocate( angle_list(num_angles), distance(num_angles) )
          end associate
          do concurrent ( ja = 1:neighbour_basis%spec(1)%num:1 )
             do concurrent ( ka = ja + 1:neighbour_basis%spec(1)%num:1 )
                idx = nint( &
                     (ja - 1) * (neighbour_basis%spec(1)%num - ja / 2.0) + &
                     (ka - ja) &
                )
                angle_list(idx) = get_angle( &
                     [ neighbour_basis%spec(1)%atom(ja,:3) ], &
                     [ neighbour_basis%spec(1)%atom(ka,:3) ] &
                )
                distance(idx) = &
                     ( &
                          modu(neighbour_basis%spec(1)%atom(ja,:3)) ** 2 * &
                          modu(neighbour_basis%spec(1)%atom(ka,:3)) ** 2 &
                     )
             end do
          end do
          this%df_3body(:,is) = this%df_3body(:,is) + &
               get_gvector( angle_list, &
                            nbins_(2), eta(2), width_(2), &
                            cutoff_min_(2), limit(2), &
                            scale = distance &
               )
          deallocate( angle_list, distance )


          !---------------------------------------------------------------------
          ! calculate the 4-body distribution function for atom (is,ia)
          !---------------------------------------------------------------------
          if(neighbour_basis%image_spec(1)%num.eq.0) cycle
          associate( &
               num_angles => &
               triangular_number( neighbour_basis%spec(1)%num - 1 ) * &
               neighbour_basis%image_spec(1)%num &
          )
             allocate( angle_list(num_angles), distance(num_angles) )
          end associate
          idx = 0
          do concurrent ( &
                 ja = 1:neighbour_basis%spec(1)%num:1, &
                 la = 1:neighbour_basis%image_spec(1)%num:1 &
          )
             do concurrent ( ka = ja + 1:neighbour_basis%spec(1)%num:1 )
                idx = nint( &
                     (ja - 1) * (neighbour_basis%spec(1)%num - ja / 2.0) + &
                     (ka - ja - 1) &
                ) * neighbour_basis%image_spec(1)%num + la
                angle_list(idx) = &
                     get_improper_dihedral_angle( &
                          [ neighbour_basis%spec(1)%atom(ja,:3) ], &
                          [ neighbour_basis%spec(1)%atom(ka,:3) ], &
                          [ neighbour_basis%image_spec(1)%atom(la,:3) ] &
                     )
                distance(idx) = &
                    modu(neighbour_basis%spec(1)%atom(ja,:3)) ** 2 * &
                    modu(neighbour_basis%spec(1)%atom(ka,:3)) ** 2 * &
                    modu(neighbour_basis%image_spec(1)%atom(la,:3)) ** 2
             end do
          end do
          this%df_4body(:,is) = this%df_4body(:,is) + &
               get_gvector( angle_list, &
                            nbins_(3), eta(3), width_(3), &
                            cutoff_min_(3), limit(3), &
                            scale = distance &
               )
          deallocate( angle_list, distance )

       end do
    end do

    !---------------------------------------------------------------------------
    ! apply the cutoff function to the 2-body distribution function
    !---------------------------------------------------------------------------
    do b = 1, nbins_(1)
       this%df_2body(b,:) = this%df_2body(b,:) / ( cutoff_min_(1) + &
            width_(1) * real(b-1, real12) ) ** 2
    end do


    !---------------------------------------------------------------------------
    ! renormalise the distribution functions so that area under the curve is 1
    !---------------------------------------------------------------------------
    do i = 1, num_pairs
       if(all(abs(this%df_2body(:,i)).lt.1.E-6))then
          this%df_2body(:,i) = 1._real12 / nbins_(1)
       end if
       this%df_2body(:,i) = this%df_2body(:,i) / sum(this%df_2body(:,i))
    end do
    do is = 1, basis%nspec
       if(all(abs(this%df_3body(:,is)).lt.1.E-6))then
          this%df_3body(:,is) = 1._real12 / nbins_(2)
       end if
       this%df_3body(:,is) = this%df_3body(:,is) / sum(this%df_3body(:,is))
       if(all(abs(this%df_4body(:,is)).lt.1.E-6))then
          this%df_4body(:,is) = 1._real12 / nbins_(3)
       end if
       this%df_4body(:,is) = this%df_4body(:,is) / sum(this%df_4body(:,is))
    end do

  end subroutine calculate
!###############################################################################


!###############################################################################
  function get_gvector(value_list, nbins, eta, width, cutoff_min, limit, &
       scale ) result(gvector)
    !! Calculate the angular distribution function for a list of values.
    implicit none

    ! Arguments
    integer, intent(in) :: nbins
    !! Number of bins for the distribution functions.
    real(real12), intent(in) :: eta, width, cutoff_min, limit
    !! Parameters for the distribution functions.
    real(real12), dimension(:), intent(in) :: value_list
    !! List of angles.
    real(real12), dimension(:), intent(in) :: scale
    !! List of scaling for each angle (distance**3 or distance**4)
    real(real12), dimension(nbins) :: gvector
    !! Distribution function for the list of values.

    ! Local variables
    integer :: i, j, b, bin
    !! Loop index.
    integer :: max_num_steps
    !! Maximum number of steps.
    integer, dimension(3,2) :: loop_limits
    !! Loop limits for the 3-body distribution function.


    max_num_steps = ceiling( sqrt(16._real12/eta) / width )
    gvector = 0._real12

    !---------------------------------------------------------------------------
    ! calculate the gvector for a list of values
    !---------------------------------------------------------------------------
    do i = 1, size(value_list)

       !------------------------------------------------------------------------
       ! get the bin closest to the value
       !------------------------------------------------------------------------
       bin = nint( ( value_list(i) - cutoff_min ) / width ) + 1


       !------------------------------------------------------------------------
       ! calculate the gaussian for this bond
       !------------------------------------------------------------------------
       ! gvector_tmp = 0._real12
       loop_limits(:,1) = &
            [ min(nbins, bin), min(nbins, bin + max_num_steps), 1 ]
       loop_limits(:,2) = &
            [ max(1, bin - 1), max(1, bin - max_num_steps), -1 ]


       !------------------------------------------------------------------------
       ! do forward and backward loops to add gaussian from its centre
       !------------------------------------------------------------------------
       do concurrent ( j = 1:2 )
          do concurrent ( &
                 b = loop_limits(1,j):loop_limits(2,j):loop_limits(3,j) )
             gvector(b) = gvector(b) + &
                  exp( -eta * ( value_list(i) - &
                                   ( width * real(b-1, real12) + &
                                     cutoff_min ) ) ** 2._real12 &
                  ) / scale(i)
          end do
       end do
    end do
    gvector = gvector * sqrt( eta / pi ) / real(size(value_list),real12)

  end function get_gvector
!###############################################################################

end module evolver