module evolver
  !! Module for handling distribution functions.
  !!
  !! This module contains the types and subroutines for handling distribution
  !! distributions. The distribution functions are used as fingerprints for
  !! atomic structures to identify similarities and differences between
  !! structures.
  use constants, only: real12, pi
  use misc_raffle, only: set, icount, strip_null, sort_str
  use misc_maths, only: lnsum, triangular_number
  use misc_linalg, only: get_angle, get_vol, get_improper_dihedral_angle, &
       cross, modu
  use rw_geom, only: basis_type, get_element_properties
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
    
    
    deallocate_systems_ = .true.
    if(present(deallocate_systems)) deallocate_systems_ = deallocate_systems

    this%num_evaluated = 0
    this%num_evaluated_allocated = 0
    if(allocated(this%total%df_2body)) deallocate(this%total%df_2body)
    if(allocated(this%total%df_3body)) deallocate(this%total%df_3body)
    if(allocated(this%total%df_4body)) deallocate(this%total%df_4body)
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
   if(.not.allocated(this%bond_info).or.size(this%bond_info).eq.0)then
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
      this%total%df_2body = this%total%df_2body * exp( this%best_energy ) / &
                              exp( best_energy_old )
      this%total%df_3body = this%total%df_3body * exp( this%best_energy ) / &
                              exp( best_energy_old )
      this%total%df_4body = this%total%df_4body * exp( this%best_energy ) / &
                              exp( best_energy_old )
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
       weight = exp( this%best_energy - energy )
       j = 0

       !------------------------------------------------------------------------
       ! loop over all species in the system to add the gvectors
       !------------------------------------------------------------------------
       do is = 1, size(this%system(i)%element_symbols)

          idx1 = findloc( [ this%element_info(:)%name ], &
                          this%system(i)%element_symbols(is), dim=1)

          height = 1._real12 / ( 1._real12 + this%total%df_3body(:,idx1) )
          this%total%df_3body(:,idx1) = this%total%df_3body(:,idx1) + &
               height * weight * this%system(i)%df_3body(:,is)
          
          height = 1._real12 / ( 1._real12 + this%total%df_4body(:,idx1) )
          this%total%df_4body(:,idx1) = this%total%df_4body(:,idx1) + &
               height * weight * this%system(i)%df_4body(:,is)
          
          do js = is, size(this%system(i)%element_symbols), 1
             idx2 = findloc( [ this%element_info(:)%name ], &
                             this%system(i)%element_symbols(js), dim=1)
             j = nint( ( size(this%element_info) - &
                         min( idx1, idx2 ) / 2._real12 ) * &
                         ( min( idx1, idx2 ) - 1._real12 ) + max( idx1, idx2 ) )
             height = 1._real12 / &
                  ( 1._real12 + this%total%df_2body(:,j) ) ** 2._real12
             this%total%df_2body(:,j) = this%total%df_2body(:,j) + &
                  height * weight * this%system(i)%df_2body(:,idx_list(is,js))
 
          end do
       end do
       deallocate(idx_list)
   end do
   
   allocate(this%norm_2body(size(this%total%df_2body,2)))
   do j = 1, size(this%total%df_2body,2)
      this%norm_2body(j) = maxval(this%total%df_2body(:,j))
      this%total%df_2body(:,j) = &
           this%total%df_2body(:,j) / this%norm_2body(j)
   end do
   allocate(this%norm_3body(size(this%element_info)))
   allocate(this%norm_4body(size(this%element_info)))
   do is = 1, size(this%element_info)
      this%norm_3body(is) = maxval(this%total%df_3body(:,is))
      this%norm_4body(is) = maxval(this%total%df_4body(:,is))
      this%total%df_3body(:,is) = &
           this%total%df_3body(:,is) / this%norm_3body(is)
      this%total%df_4body(:,is) = &
           this%total%df_4body(:,is) / this%norm_4body(is)
   end do

   this%num_evaluated_allocated = size(this%system)
   this%num_evaluated = this%num_evaluated + num_evaluated

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

    integer :: bin, max_num_steps
    !! Bin index and maximum number of steps.
    integer :: i, j, k, b, itmp1
    !! Loop index.
    integer :: is, js, ia, ja
    !! Loop index.
    integer :: num_pairs, num_angles
    !! Number of pairs and angles.
    integer :: amax, bmax, cmax
    !! Maximum number of lattice vectors to consider.
    real(real12) :: rtmp1, rtmp2
    !! Temporary real variables.
    real(real12) :: fc, weight, scale
    !! Cutoff, weight and scale for the distribution functions.
    logical :: success
    !! Boolean for success.
    real(real12), dimension(3) :: eta, limit
    !! Parameters for the distribution functions.
    real(real12), dimension(3) :: vtmp1, vtmp2, vtmp3, diff
    !! Temporary real arrays.
    real(real12), allocatable, dimension(:) :: gvector_tmp, angle, distance
    !! Temporary real arrays.
    integer, dimension(:), allocatable :: idx_list!, count_list
    !! Index list for the element pairs in a system.

    integer, dimension(3,2) :: loop_limits
    !! Loop limits for the 3-body distribution function.
    integer, allocatable, dimension(:,:) :: idx
    !! Index list for the element pairs in a system.
    integer, allocatable, dimension(:,:) :: pair_index
    !! Index of element pairs.

    type :: bond_type
       !! Derived type for a bond.
       integer, dimension(2) :: species, atom
       !! Species and atom indices.
       logical :: skip = .false.
       !! Boolean whether to skip the bond.
       real(real12), dimension(3) :: vector
       !! Vector of the bond.
    end type bond_type
    type(bond_type), dimension(:), allocatable :: bond_list
    !! List of bonds in the system.

    !type :: plane_type
    !   real(real12), dimension(3) :: vector
    !   integer :: count
    !end type plane_type
    !type(plane_type), dimension(:), allocatable :: plane_list


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
    ! get the maximum number of lattice vectors to consider
    ! NOTE: this is not perfect
    !       won't work for extremely acute/obtuse angle cells
    !       (due to diagonal path being shorter than individual lattice vectors)
    !---------------------------------------------------------------------------
    amax = ceiling(cutoff_max_(1)/modu(basis%lat(1,:)))
    bmax = ceiling(cutoff_max_(1)/modu(basis%lat(2,:)))
    cmax = ceiling(cutoff_max_(1)/modu(basis%lat(3,:)))


    !---------------------------------------------------------------------------
    ! build the bond list
    !---------------------------------------------------------------------------
    allocate(bond_list(0)) !if doesn't work, allocate a dummy bond first
    spec_loop1: do is=1,basis%nspec
       atom_loop1: do ia=1,basis%spec(is)%num
          spec_loop2: do js=is,basis%nspec
             atom_loop2: do ja=1,basis%spec(js)%num
                if(is.eq.js.and.ja.lt.ia) cycle atom_loop2
                diff = basis%spec(is)%atom(ia,:3) -  basis%spec(js)%atom(ja,:3)
                diff = diff - ceiling(diff - 0.5_real12)
                do i=-amax,amax+1,1
                   vtmp1(1) = diff(1) + real(i, real12)
                   do j=-bmax,bmax+1,1
                      vtmp1(2) = diff(2) + real(j, real12)
                      do k=-cmax,cmax+1,1
                         vtmp1(3) = diff(3) + real(k, real12)
                         rtmp1 = modu(matmul(vtmp1,basis%lat))
                         if( rtmp1 .gt. cutoff_min_(1) - &
                                        width_(1)/2._real12 .and. &
                             rtmp1 .lt. cutoff_max_(1) + &
                                        width_(1)/2._real12 )then
                            bond_list = [ bond_list, bond_type( &
                                 species=[is,js], &
                                 atom=[ia,ja], skip=.false., &
                                 vector=matmul(vtmp1,basis%lat)) ]
                         end if
                      end do
                   end do
                end do
             end do atom_loop2
          end do spec_loop2
       end do atom_loop1
    end do spec_loop1


    !---------------------------------------------------------------------------
    ! calculate the gaussian width
    !---------------------------------------------------------------------------
    eta = 1._real12 / ( 2._real12 * sigma_**2._real12 )
    max_num_steps = ceiling( sqrt(16._real12/eta(1)) / width_(1) )


    !---------------------------------------------------------------------------
    ! build the 2-body gvectors (radial distribution functions)
    !---------------------------------------------------------------------------
    allocate(this%df_2body(nbins_(1),num_pairs), source = 0._real12)
    allocate(gvector_tmp(nbins_(1)),             source = 0._real12)
    do i = 1, size(bond_list)
       if(bond_list(i)%skip) cycle
       if(abs(modu(bond_list(i)%vector)).lt.1.E-3) cycle
       is = bond_list(i)%species(1)
       js = bond_list(i)%species(2)
       rtmp1 = modu(bond_list(i)%vector)
       !------------------------------------------------------------------------
       ! get number of equivalent bonds
       !------------------------------------------------------------------------
       scale = 1._real12
       do j = i+1, size(bond_list)
          !! don't need to look at reverse of species, as the ordering will ...
          !! always be enforced by the loop above, i.e. is <= js
          if( is .ne. bond_list(j)%species(1) .or. &
              js .ne. bond_list(j)%species(2) ) cycle
          if(abs(modu(bond_list(j)%vector)-rtmp1).lt.1.E-3)then !!MAKE THIS LINKED TO WIDTH?
             bond_list(j)%skip = .true.
             scale = scale + 1._real12
          end if          
       end do
       !------------------------------------------------------------------------
       bin = nint( ( rtmp1 - cutoff_min_(1) ) / width_(1) ) + 1
       if(bin.gt.nbins_(1).or.bin.lt.1) cycle

       fc = 0.5_real12 * cos( pi * ( rtmp1 - cutoff_min_(1) ) / limit(1) ) + &
            0.5_real12

       !------------------------------------------------------------------------
       ! calculate the gaussian for this bond
       !------------------------------------------------------------------------
       gvector_tmp = 0._real12
       loop_limits(:,1) = [ bin, min(nbins_(1), (bin + max_num_steps) ), 1 ]
       loop_limits(:,2) = [ bin - 1, max(1, bin - max_num_steps), -1 ]
       
       ! do forward and backward loops to add gaussian for larger distances
       do concurrent ( j = 1:2 )
          do concurrent ( b = &
                            loop_limits(1,j):loop_limits(2,j):loop_limits(3,j) )
             gvector_tmp(b) = gvector_tmp(b) + &
                  exp( -eta(1) * ( rtmp1 - &
                                   ( width_(1) * real(b-1, real12) + &
                                     cutoff_min_(1) ) ) ** 2._real12 )
          end do
       end do
       itmp1 = count( [ ( ( bond_list(j)%species(1) .eq. is .and. &
                            bond_list(j)%species(2) .eq. js ) .or. &
                          ( bond_list(j)%species(2) .eq. is .and. &
                            bond_list(j)%species(1) .eq. js ), &
                              j = 1, size(bond_list), 1 ) ] )
       this%df_2body(:,pair_index(is, js)) = &
            this%df_2body(:,pair_index(is, js)) + &
            gvector_tmp * scale * sqrt( eta(1) / pi ) / real(itmp1,real12) ! / width_(1)
    end do
    do b = 1, nbins_(1)
       this%df_2body(b,:) = this%df_2body(b,:) / ( cutoff_min_(1) + &
            width_(1) * real(b-1, real12) ) ** 2
    end do
    ! renormalise so that the area under the curve is 1
    do k = 1, num_pairs
       if(all(abs(this%df_2body(:,k)).lt.1.E-6))then
          this%df_2body(:,k) = 1._real12 / nbins_(1)
       end if
       this%df_2body(:,k) = this%df_2body(:,k) / sum(this%df_2body(:,k))
    end do


    !---------------------------------------------------------------------------
    ! build the 3-body gvectors (angular distribution functions)
    !---------------------------------------------------------------------------
    deallocate(gvector_tmp)
    allocate(this%df_3body(nbins_(2), basis%nspec), source = 0._real12)
    allocate(gvector_tmp(nbins_(2)),                source = 0._real12)
    do is = 1, basis%nspec
       num_angles = 0
       ! number of comibnations without repetitions:
       !    = n! / (n - r)! r!
       !    as r = 1, this simplifies to n
       do i = 1, size(bond_list), 1
          if( is .eq. bond_list(i)%species(1) )then
             ia = bond_list(i)%atom(1)
          elseif( is .eq. bond_list(i)%species(2) )then
             ia = bond_list(i)%atom(2)
          else
             cycle
          end if
          num_angles = num_angles + count( &
             [ ( ( bond_list(j)%species(1) .eq. is .and. &
                   bond_list(j)%atom(1) .eq. ia ) .or. &
                 ( bond_list(j)%species(2) .eq. is .and. &
                   bond_list(j)%atom(2) .eq. ia ), &
                     j = i + 1, size(bond_list), 1 ) ] )
       end do
       allocate(angle(num_angles))
       allocate(distance(num_angles))
       num_angles = 0


       !------------------------------------------------------------------------
       ! loop over all bonds to find the first bond
       !------------------------------------------------------------------------
       do i = 1, size(bond_list)
          if( is .eq. bond_list(i)%species(1) )then
             vtmp1 = -bond_list(i)%vector
             ia = bond_list(i)%atom(1)
             js = bond_list(i)%species(2)
          elseif( is .eq. bond_list(i)%species(2) )then
             vtmp1 = bond_list(i)%vector
             ia = bond_list(i)%atom(2)
             js = bond_list(i)%species(1)
            else
             cycle
          end if
          if( &
               modu(vtmp1) .lt. &
                    bond_info(pair_index(is, js))%radius_covalent * &
                    radius_distance_tol_(1) .or. &
               modu(vtmp1) .gt. &
                    bond_info(pair_index(is, js))%radius_covalent * &
                    radius_distance_tol_(2) &
          ) cycle
        
          !---------------------------------------------------------------------
          ! loop over all bonds to find the second bond
          !---------------------------------------------------------------------
          do j = i + 1, size(bond_list)
             if( is .eq. bond_list(j)%species(1) .and. &
                 ia .eq. bond_list(j)%atom(1) )then
                vtmp2 = -bond_list(j)%vector
                js = bond_list(j)%species(2)
             elseif( is .eq. bond_list(j)%species(2) .and. &
                     ia .eq. bond_list(j)%atom(2) )then
                vtmp2 = bond_list(j)%vector
                js = bond_list(j)%species(1)
               else
                cycle
             end if
             if( &
                  modu(vtmp2) .lt. &
                       bond_info(pair_index(is, js))%radius_covalent * &
                       radius_distance_tol_(1) .or. &
                  modu(vtmp2) .gt. &
                       bond_info(pair_index(is, js))%radius_covalent * &
                       radius_distance_tol_(2) &
             ) cycle
 
             if( abs( &
                  dot_product(vtmp1, vtmp1) - &
                  dot_product(vtmp1, vtmp2) &
             ) .lt. 1.E-3 ) cycle

             !------------------------------------------------------------------
             ! calculate the angle between the two bonds
             !------------------------------------------------------------------
             num_angles = num_angles + 1
             angle(num_angles) = get_angle( vtmp1, vtmp2 )
             distance(num_angles) = &
                  ( modu(vtmp1) ** 2 * modu(vtmp2) ** 2 )

          end do
       end do
       if(num_angles.gt.size(angle))then
          write(0,*) "ERROR: Number of 3-body angles exceeds allocated array"
          write(0,'("Expected ",I0," got ",I0)') size(angle), num_angles
          stop 1
       elseif(num_angles.gt.0)then
          this%df_3body(:,is) = this%df_3body(:,is) + &
               get_gvector( angle(:num_angles), nbins_(2), eta(2), width_(2), &
                                  cutoff_min_(2), &
                                  limit(2), scale = distance(:num_angles) &
               )
       end if
       deallocate(angle)
       deallocate(distance)
    end do
    ! renormalise so that the area under the curve is 1
    do is = 1, basis%nspec
       if(all(abs(this%df_3body(:,is)).lt.1.E-6))then
          this%df_3body(:,is) = 1._real12 / nbins_(2)
       end if
       this%df_3body(:,is) = this%df_3body(:,is) / sum(this%df_3body(:,is))
    end do

  
    !---------------------------------------------------------------------------
    ! build the 4-body gvectors (angular distribution functions)
    !---------------------------------------------------------------------------
    deallocate(gvector_tmp)
    allocate(this%df_4body(nbins_(3),basis%nspec), source = 0._real12)
    allocate(gvector_tmp(nbins_(3)),               source = 0._real12)
    do is = 1, basis%nspec
       num_angles = 0
       ! number of comibnations without repetitions:
       ! = n! / (n - r)! r!
       do i = 1, size(bond_list), 1
          if( is .eq. bond_list(i)%species(1) )then
             ia = bond_list(i)%atom(1)
          elseif( is .eq. bond_list(i)%species(2) )then
             ia = bond_list(i)%atom(2)
          else
             cycle
          end if
          itmp1 = count( &
             [ ( ( bond_list(j)%species(1) .eq. is .and. &
                   bond_list(j)%atom(1) .eq. ia ) .or. &
                 ( bond_list(j)%species(2) .eq. is .and. &
                   bond_list(j)%atom(2) .eq. ia ), &
                     j = i+1, size(bond_list), 1 ) ] )
          num_angles = num_angles + &
               nint(exp(lnsum(itmp1) - lnsum(itmp1 - 2) - lnsum(2)))
       end do
       allocate(angle(num_angles))
       allocate(distance(num_angles))
       !allocate(count_list(num_angles))
       num_angles = 0

       !------------------------------------------------------------------------
       ! loop over all bonds to find the first bond
       !------------------------------------------------------------------------
       do i = 1, size(bond_list)
          if(abs(modu(bond_list(i)%vector)).lt.1.E-3) cycle
 
          if( is .eq. bond_list(i)%species(1) )then
             vtmp1 = -bond_list(i)%vector
             ia = bond_list(i)%atom(1)
             js = bond_list(i)%species(2)
          elseif( is .eq. bond_list(i)%species(2) )then
             vtmp1 = bond_list(i)%vector
             ia = bond_list(i)%atom(2)
             js = bond_list(i)%species(1)
          else
             cycle
          end if
          if( &
               modu(vtmp1) .lt. &
                    bond_info(pair_index(is, js))%radius_covalent * &
                    radius_distance_tol(1) .or. &
               modu(vtmp1) .gt. &
                    bond_info(pair_index(is, js))%radius_covalent * &
                    radius_distance_tol(2) &
          ) cycle

          ! ! make list of indices where species is in bond_list(i)%species and atom is in bond_list(i)%atom
          allocate(idx_list(count( [ ( ( bond_list(j)%species(1) .eq. is .and. &
                                         bond_list(j)%atom(1) .eq. ia ) .or. &
                                       ( bond_list(j)%species(2) .eq. is .and. &
                                         bond_list(j)%atom(2) .eq. ia ), &
                                           j = i + 1, size(bond_list), 1 ) ] )))
          k = 0
          index_list_loop: do j = i + 1, size(bond_list)
             if(abs(modu(bond_list(j)%vector)).lt.1.E-3) cycle
             if( ( is .eq. bond_list(j)%species(1) .and. &
                   ia .eq. bond_list(j)%atom(1)  ) )then
                k = k + 1
                idx_list(k) = -j
             elseif( is .eq. bond_list(j)%species(2) .and. &
                     ia .eq. bond_list(j)%atom(2) ) then
                k = k + 1
                idx_list(k) = j
             end if
             ! ! get a list of unique plane normal vectors
             ! vtmp2 = cross( vtmp1, bond_list(j)%vector )
             ! vtmp2 = vtmp2 / modu(vtmp2)
             ! plane_check_loop: do k = 1, size(plane_list)
             !    if(abs(dot_product(vtmp2, plane_list(k)%vector)) .lt. 1.E-3)then
             !       plane_list(k)%count = plane_list(k)%count + 1
             !       cycle index_list_loop
             !    end if
             ! end do plane_check_loop
             ! plane_list = [ plane_list, plane_type(vector = vtmp2, count = 1) ]
          end do index_list_loop

          !---------------------------------------------------------------------
          ! loop over all bonds to find the second bond
          !---------------------------------------------------------------------
          do j = 1, size(idx_list)!size(plane_list)!size(idx_list)
 
             if( idx_list(j) .lt. 0 )then
                vtmp2 = -bond_list(-idx_list(j))%vector
                js = bond_list(-idx_list(j))%species(2)
             else
                vtmp2 = bond_list(idx_list(j))%vector
                js = bond_list(idx_list(j))%species(1)
             end if
             if( &
                  modu(vtmp2) .lt. &
                       bond_info(pair_index(is, js))%radius_covalent * &
                       radius_distance_tol(1) .or. &
                  modu(vtmp2) .gt. &
                       bond_info(pair_index(is, js))%radius_covalent * &
                       radius_distance_tol(2) &
             ) cycle
             if( abs( &
                  dot_product(vtmp1, vtmp1) - &
                  dot_product(vtmp1, vtmp2) &
             ) .lt. 1.E-3 ) cycle
             !------------------------------------------------------------------
             ! loop over all bonds to find the third bond
             !------------------------------------------------------------------
             !count_list(num_angles+1:num_angles+(size(idx_list)-j)) = plane_list(j)%count
             do k = abs(j) + 1, size(idx_list)

                if( idx_list(k) .lt. 0 )then
                   vtmp3 = -bond_list(-idx_list(k))%vector
                   js = bond_list(-idx_list(k))%species(2)
                else
                   vtmp3 = bond_list(idx_list(k))%vector
                   js = bond_list(idx_list(k))%species(1)
                end if
                if( &
                     modu(vtmp3) .lt. &
                          bond_info(pair_index(is, js))%radius_covalent * &
                          radius_distance_tol(3) .or. &
                     modu(vtmp3) .gt. &
                          bond_info(pair_index(is, js))%radius_covalent * &
                          radius_distance_tol(4) &
                ) cycle
 
                !---------------------------------------------------------------
                ! calculate the angle between the two bonds
                !---------------------------------------------------------------
                ! rtmp1 = get_angle( vtmp2, &
                !      bond_list(idx_list(k))%vector )
                ! if(rtmp1 .gt. pi/2._real12) rtmp1 = pi - rtmp1
                ! if(any(abs(rtmp1 - angle(:num_angles)) .lt. 1.E-3)) cycle
                num_angles = num_angles + 1
                ! angle(num_angles) = rtmp1
                angle(num_angles) = &
                     get_improper_dihedral_angle( &
                          vtmp1, &
                          vtmp2, &
                          vtmp3 &
                     )

                distance(num_angles) = &
                    modu(vtmp1) ** 2 * &
                    modu(vtmp2) ** 2 * &
                    modu(vtmp3) ** 2
                
                ! count_list(num_angles) = plane_list(j)%count

             end do
             ! num_angles = num_angles + size(idx_list) - j
          end do
          deallocate(idx_list)
          ! deallocate(plane_list)
       end do
       if(num_angles.gt.size(angle))then
          write(0,*) "ERROR: Number of 4-body angles exceeds allocated array"
          write(0,'("Expected ",I0," got ",I0)') size(angle), num_angles
          stop 1
       elseif(num_angles.gt.0)then
          this%df_4body(:,is) = this%df_4body(:,is) + &
               get_gvector( angle(:num_angles), nbins_(3), eta(3), width_(3), &
                                  cutoff_min_(3), &
                                  limit(3), &
                                  scale = distance(:num_angles) &
               )!, count_list(:num_angles) )
       end if
       deallocate(angle)
       deallocate(distance)
       ! deallocate(count_list)
    end do
    ! renormalise so that the area under the curve is 1
    do is = 1, basis%nspec
      if(all(abs(this%df_4body(:,is)).lt.1.E-6))then
         this%df_4body(:,is) = 1._real12 / nbins_(3)
      end if
      this%df_4body(:,is) = this%df_4body(:,is) / sum(this%df_4body(:,is))
   end do

  end subroutine calculate
!###############################################################################


!###############################################################################
  function get_gvector(vector, nbins, eta, width, cutoff_min, limit, &
       scale ) result(gvector)
    !! Calculate the angular distribution function for a list of vectors.
    implicit none

    ! Arguments
    integer, intent(in) :: nbins
    !! Number of bins for the distribution functions.
    real(real12), intent(in) :: eta, width, cutoff_min, limit
    !! Parameters for the distribution functions.
    real(real12), dimension(:), intent(in) :: vector
    !! List of angles.
    real(real12), dimension(:), intent(in) :: scale
    !! List of scaling for each angle (distance**3 or distance**4)
    real(real12), dimension(nbins) :: gvector
    !! Distribution function for the list of vectors.
    ! integer, dimension(:), intent(in), optional :: count_list

    ! Local variables
    integer :: i, j, b, bin
    !! Loop index.
    integer :: max_num_steps
    !! Maximum number of steps.
    !integer, dimension(:), allocatable :: scale_list
    !real(real12), dimension(nbins) :: gvector_tmp
    !real(real12), dimension(:), allocatable :: vector_copy
    integer, dimension(3,2) :: loop_limits
    !! Loop limits for the 3-body distribution function.


    ! max_num_steps = ceiling( abs( vector_copy(i) - sqrt(16._real12/eta) ) / width )
    max_num_steps = ceiling( sqrt(16._real12/eta) / width )
    gvector = 0._real12

    !---------------------------------------------------------------------------
    ! calculate the gvector for a list of vectors
    !---------------------------------------------------------------------------
    ! vector_copy = vector
    ! itmp1 = 0
    ! !  order the vector list, and remove duplicates within a tolerance
    ! call set(vector_copy, 1.E-3, scale_list)
    do i = 1, size(vector)
        ! if( vector_copy(i) .lt. -1.E-3 ) cycle
        ! !-----------------------------------------------------------------------
        ! ! remove duplicates
        ! !-----------------------------------------------------------------------
        ! scale = 1._real12
        ! do j = i + 1, size(vector_copy)
        !    if(abs(vector_copy(i)-vector_copy(j)) .lt. 1.E-3 )then
        !       !vector_copy(j) = -1.E3_real12
        !       !scale = scale + 1._real12
        !       !itmp1 = itmp1 + 1
        !    end if
        ! end do


       !------------------------------------------------------------------------
       ! get the bin closest to the value
       !------------------------------------------------------------------------
       bin = nint( ( vector(i) - cutoff_min ) / width ) + 1
       ! if(bin.gt.nbins.or.bin.lt.1) cycle


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
                  exp( -eta * ( vector(i) - &
                                   ( width * real(b-1, real12) + &
                                     cutoff_min ) ) ** 2._real12 &
                  ) / scale(i)
          end do
       end do
       ! if(present(count_list)) gvector_tmp = gvector_tmp * count_list(i)
       ! gvector(:) = gvector(:) + gvector_tmp !* scale !real(scale_list(i), real12)
    end do
    gvector = gvector * sqrt( eta / pi ) / real(size(vector),real12) ! / width

  end function get_gvector
!###############################################################################

end module evolver