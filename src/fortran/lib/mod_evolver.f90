module evolver
  use constants, only: real12, pi
  use misc_raffle, only: set, icount, strip_null
  use misc_maths, only: lnsum, triangular_number
  use misc_linalg, only: get_angle, get_vol, cross, modu
  use rw_geom, only: bas_type
  use elements, only: &
       element_type, element_bond_type, &
       load_elements, load_element_bonds, &
       element_database, element_bond_database
  implicit none

  
  private

  public :: gvector_container_type, gvector_base_type, gvector_type
  
  
  type :: gvector_base_type
     real(real12), dimension(:,:), allocatable :: df_2body
     real(real12), dimension(:,:), allocatable :: df_3body
     real(real12), dimension(:,:), allocatable :: df_4body
  end type gvector_base_type

  type, extends(gvector_base_type) :: gvector_type
     integer :: num_atoms = 0
     real(real12) :: energy = 0.0_real12 !! should be formation energy
     integer, dimension(:), allocatable :: stoichiometry
     character(len=3), dimension(:), allocatable :: species
   contains
     procedure, pass(this) :: calculate
  end type gvector_type


  !! should gvector be for the entire prediction, or one for each system?
  !! if one for each system, do not contain nbins, width, as this would ...
  !! ... result in lots of duplicated data
  type :: gvector_container_type
     integer :: best_system = 0
     real(real12) :: best_energy = 0.0_real12
     integer, dimension(3) :: nbins = -1
     real(real12), dimension(3) :: sigma = [ 0.1_real12, 0.05_real12, 0.05_real12 ]
     real(real12), dimension(3) :: width = [ 0.025_real12, pi/24._real12, pi/32._real12 ]
     real(real12), dimension(3) :: cutoff_min = [ 0.5_real12, 0._real12, 0._real12 ]
     real(real12), dimension(3) :: cutoff_max = [ 6._real12, pi, pi/2._real12 ]
     type(gvector_base_type) :: total !! name it best instead?
     type(gvector_type), dimension(:), allocatable :: system
     type(element_type), dimension(:), allocatable :: element_info
     type(element_bond_type), dimension(:), allocatable :: bond_info
   contains
     procedure, pass(this) :: set_width
     procedure, pass(this) :: set_sigma
     procedure, pass(this) :: set_cutoff_min
     procedure, pass(this) :: set_cutoff_max

     procedure, pass(this) :: create
     procedure, pass(this) :: update
     
     procedure, pass(this) :: add, add_basis

     procedure, pass(this), private :: set_element_info
     procedure, pass(this), private :: update_element_info
     procedure, pass(this) :: set_element_energy
     procedure, pass(this) :: set_element_energies
     procedure, pass(this) :: get_element_energies
     procedure, pass(this) :: get_element_energies_staticmem

     procedure, pass(this) :: set_bond_info
     
     procedure, pass(this) :: set_best_energy
     procedure, pass(this) :: initialise_gvectors
     procedure, pass(this) :: evolve
     procedure, pass(this) :: write
     procedure, pass(this) :: read
     procedure, pass(this) :: write_2body
     procedure, pass(this) :: write_3body
     procedure, pass(this) :: write_4body
     procedure, pass(this) :: get_pair_index
     procedure, pass(this) :: get_bin
  !   procedure :: read
  end type gvector_container_type

  interface gvector_container_type
    module function init_gvector_container( &
         nbins, width, sigma, cutoff_min, cutoff_max &
         ) result(gvector_container)
         integer, dimension(3), intent(in), optional :: nbins
         real(real12), dimension(3), intent(in), optional :: width, sigma
         real(real12), dimension(3), intent(in), optional :: cutoff_min, cutoff_max
         type(gvector_container_type) :: gvector_container
    end function init_gvector_container
  end interface gvector_container_type


  contains
!!!#############################################################################
!!! initialise gvector container type
!!!#############################################################################
  module function init_gvector_container(nbins, width, sigma, cutoff_min, cutoff_max) &
       result(gvector_container)
    implicit none
    integer, dimension(3), intent(in), optional :: nbins
    real(real12), dimension(3), intent(in), optional :: width, sigma
    real(real12), dimension(3), intent(in), optional :: cutoff_min, cutoff_max
    type(gvector_container_type) :: gvector_container

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
       if(all(cutoff_min.ge.0._real12)) gvector_container%cutoff_min = cutoff_min
    end if
    if(present(cutoff_max))then
       if(all(cutoff_max.ge.0._real12)) gvector_container%cutoff_max = cutoff_max
    end if
    if(any(gvector_container%cutoff_max .le. gvector_container%cutoff_min))then
       write(0,*) "ERROR: cutoff_max <= cutoff_min"
       write(0,*) "cutoff min: ", gvector_container%cutoff_min
       write(0,*) "cutoff max: ", gvector_container%cutoff_max
       stop 1
    end if

  end function init_gvector_container
!!!#############################################################################

  subroutine set_width(this, width)
    !! Set the width of the gaussians used in the 2-, 3-, and 4-body
    !! distribution functions.
    implicit none

    ! Arguments
    class(gvector_container_type), intent(inout) :: this
    !! Self, parent of the procedure
    real(real12), dimension(3), intent(in) :: width
    !! Width of the gaussians used in the 2-, 3-, and 4-body
    !! distribution functions.

    this%width = width

  end subroutine set_width


  subroutine set_sigma(this, sigma)
    !! Set the sigma of the gaussians used in the 2-, 3-, and 4-body
    !! distribution functions.
    implicit none

    ! Arguments
    class(gvector_container_type), intent(inout) :: this
    !! Self, parent of the procedure.
    real(real12), dimension(3), intent(in) :: sigma
    !! Sigma of the gaussians used in the 2-, 3-, and 4-body
    !! distribution functions.

    this%sigma = sigma

  end subroutine set_sigma


  subroutine set_cutoff_min(this, cutoff_min)
    !! Set the minimum cutoff for the 2-, 3-, and 4-body distribution functions.
    implicit none

    ! Arguments
    class(gvector_container_type), intent(inout) :: this
    !! Self, parent of the procedure.
    real(real12), dimension(3), intent(in) :: cutoff_min
    !! Minimum cutoff for the 2-, 3-, and 4-body distribution functions.

     this%cutoff_min = cutoff_min

  end subroutine set_cutoff_min


  subroutine set_cutoff_max(this, cutoff_max)
    !! Set the maximum cutoff for the 2-, 3-, and 4-body distribution functions.
    implicit none
   
    ! Arguments
    class(gvector_container_type), intent(inout) :: this
    !! Self, parent of the procedure.
    real(real12), dimension(3), intent(in) :: cutoff_max
    !! Maximum cutoff for the 2-, 3-, and 4-body distribution functions.
   
    this%cutoff_max = cutoff_max
   
  end subroutine set_cutoff_max


  subroutine create(this, basis_list)
    !! create the distribution functions from the input file
    implicit none
    ! Arguments
    class(gvector_container_type), intent(inout) :: this
    !! Self, parent of the procedure.
    type(bas_type), dimension(:), intent(in) :: basis_list
    !! List of basis structures.

    if(allocated(this%total%df_2body)) deallocate(this%total%df_2body)
    if(allocated(this%total%df_3body)) deallocate(this%total%df_3body)
    if(allocated(this%total%df_4body)) deallocate(this%total%df_4body)
    if(allocated(this%system)) deallocate(this%system)
    allocate(this%system(0))
    call this%add(basis_list)
    call this%set_bond_info()
    call this%evolve()
    
  end subroutine create


  subroutine update(this, basis_list)
    !! update the distribution functions from the input file
    implicit none
    ! Arguments
    class(gvector_container_type), intent(inout) :: this
    !! Self, parent of the procedure.
    type(bas_type), dimension(:), intent(in) :: basis_list
    !! List of basis structures.

    
    call this%add(basis_list)
    call this%set_bond_info()
    call this%evolve()
    
  end subroutine update


!!!#############################################################################
!!! write all systems
!!!#############################################################################
  subroutine write(this, file)
    implicit none
    class(gvector_container_type), intent(in) :: this
    character(*), intent(in) :: file

    integer :: unit, i, j
    integer, allocatable, dimension(:,:) :: idx

    open(newunit=unit, file=file)
    write(unit, *) "nbins", this%nbins
    write(unit, *) "width", this%width
    write(unit, *) "sigma", this%sigma
    write(unit, *) "cutoff_min", this%cutoff_min
    write(unit, *) "cutoff_max", this%cutoff_max
    write(unit, *)
    do i = 1, size(this%system,1)
       write(unit, *) this%system(i)%energy
       write(unit, *) this%system(i)%species
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
!!!#############################################################################


!!!#############################################################################
!!! read all systems
!!!#############################################################################
  subroutine read(this, file)
    implicit none
    class(gvector_container_type), intent(inout) :: this
    character(*), intent(in) :: file

    integer :: unit, i, j, ierror
    integer :: num_species, num_pairs
    character(256) :: buffer
    type(gvector_type) :: system

   
    open(newunit=unit, file=file)
    read(unit, *) buffer, this%nbins
    read(unit, *) buffer, this%width
    read(unit, *) buffer, this%sigma
    read(unit, *) buffer, this%cutoff_min
    read(unit, *) buffer, this%cutoff_max
    do
       read(unit, '(A)', iostat=ierror) buffer
       if(ierror.ne.0) exit
       if(trim(buffer).eq.''.or.trim(buffer).eq.'#') cycle
       read(buffer, *) system%energy
       read(unit, '(A)') buffer
       num_species = icount(buffer)
       allocate(system%species(num_species))
       allocate(system%stoichiometry(num_species))
       read(buffer, *) system%species
       read(unit, *) system%stoichiometry
       system%num_atoms = sum(system%stoichiometry)
       num_pairs = gamma(real(num_species + 2, real12)) / &
                   ( gamma(real(num_species, real12)) * gamma( 3._real12 ) )
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
       deallocate(system%species, system%stoichiometry, &
                  system%df_2body, system%df_3body, system%df_4body)
    end do
    close(unit)

  end subroutine read
!!!#############################################################################


!!!#############################################################################
!!! write the 2body gvectors to a file
!!!#############################################################################
  subroutine write_2body(this, file)
    implicit none
    class(gvector_container_type), intent(in) :: this
    character(*), intent(in) :: file

    integer :: unit, i, j, is, js
    integer :: num_pairs
    integer, allocatable, dimension(:,:) :: idx


    num_pairs = gamma(real(size(this%element_info) + 2, real12)) / &
                ( gamma(real(size(this%element_info), real12)) * &
                  gamma( 3._real12 ) )
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
!!!#############################################################################


!!!#############################################################################
!!! write the 3body gvectors to a file
!!!#############################################################################
    subroutine write_3body(this, file)
    implicit none
    class(gvector_container_type), intent(in) :: this
    character(*), intent(in) :: file

    integer :: unit, i, j


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
!!!#############################################################################


!!!#############################################################################
!!! write the 3body gvectors to a file
!!!#############################################################################
    subroutine write_4body(this, file)
    implicit none
    class(gvector_container_type), intent(in) :: this
    character(*), intent(in) :: file

    integer :: unit, i, j


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
!!!#############################################################################


!!!#############################################################################
!!! add system (basis or gvector) to the container
!!!#############################################################################
  subroutine add(this, system)
    implicit none
    class(gvector_container_type), intent(inout) :: this
    class(*), dimension(..), intent(in) :: system

    integer :: i, num_structures_previous
    character(128) :: buffer


    select rank(system)
    rank(0)
       select type(system)
       type is (gvector_type)
          this%system = [ this%system, system ]
       type is (bas_type)
          call this%add_basis(system)
       class default
          write(0,*) "ERROR: Invalid type for system"
          write(0,*) "Expected type gvector_type or bas_type"
          stop 1
       end select
    rank(1)
       num_structures_previous = size(this%system)
       select type(system)
       type is (gvector_type)
          this%system = [ this%system, system ]
       type is (bas_type)
          do i = 1, size(system)
             call this%add_basis(system(i))
          end do
       class default
          write(0,*) "ERROR: Invalid type for system"
          write(0,*) "Expected type gvector_type or bas_type"
          stop 1
       end select
    rank default
       write(0,*) "ERROR: Invalid rank for system"
       write(buffer,*) rank(system)
       write(0,*) "Expected rank 0 or 1, got ", trim(buffer)
       stop 1
    end select
    call this%update_element_info()
    call this%set_best_energy()

  end subroutine add
!!!#############################################################################


!!!#############################################################################
!!! generate gvectors from basis and add to the container
!!!#############################################################################
  subroutine add_basis(this, basis)
    implicit none
    class(gvector_container_type), intent(inout) :: this
    type(bas_type), intent(in) :: basis

    type(gvector_type) :: system

    call system%calculate(basis, width = this%width, &
                     sigma = this%sigma, &
                     cutoff_min = this%cutoff_min, &
                     cutoff_max = this%cutoff_max)

    if(.not.allocated(this%system))then
       this%system = [ system ]
    else
       this%system = [ this%system, system ]
    end if
  end subroutine add_basis
!!!#############################################################################


!!!#############################################################################
!!! set the species list for the container
!!!#############################################################################
  subroutine set_element_info(this)
    implicit none
    class(gvector_container_type), intent(inout) :: this

    integer :: i, unit
    character(len=3), dimension(:), allocatable :: element_list


    !!--------------------------------------------------------------------------
    !! get list of species in dataset
    !!--------------------------------------------------------------------------
    element_list = [ this%system(1)%species ]
    do i = 2, size(this%system),1
       element_list = [ element_list, this%system(i)%species ]
    end do
    call set(element_list)
    if(allocated(this%element_info)) deallocate(this%element_info)
    allocate(this%element_info(size(element_list)))
    do i = 1, size(element_list)
       call this%element_info(i)%set(element_list(i))
    end do
    
  end subroutine set_element_info
!!!#############################################################################


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
    element_list = [ this%system(1)%species ]
    do i = 2, size(this%system),1
       element_list = [ element_list, this%system(i)%species ]
    end do
    call set(element_list)


    !---------------------------------------------------------------------------
    ! check if all elements are in the element_info array
    !---------------------------------------------------------------------------
    do i = 1, size(element_list)
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
    character(len=3) :: element_


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
       element_database = [ &
            element_database(:), &
            element_type(name = element_, energy = energy) &
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


! !###############################################################################
  subroutine get_element_energies_staticmem(this, elements, energies)
    !! Return the energies of elements in the container.
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
! !###############################################################################


!!!#############################################################################
!!! set the element bond list for the container
!!!#############################################################################
  subroutine set_bond_info(this, bond_file)
    implicit none
    class(gvector_container_type), intent(inout) :: this
    character(*), intent(in), optional :: bond_file

    integer :: i, j, k
    integer :: num_elements, num_pairs


    !!--------------------------------------------------------------------------
    !! load the element bonds database
    !!--------------------------------------------------------------------------
    if(present(bond_file))then
       call load_element_bonds(bond_file)
    elseif(.not.allocated(element_bond_database))then
       call load_element_bonds()
    end if


    !!--------------------------------------------------------------------------
    !! allocate the bond information array
    !!--------------------------------------------------------------------------
    num_elements = size(this%element_info)
    num_pairs = nint(gamma(real(num_elements + 2, real12)) / &
         ( gamma(real(num_elements, real12)) * gamma( 3._real12 ) ) )
    allocate(this%bond_info(num_pairs))


    !!--------------------------------------------------------------------------
    !! loop over all pairs of elements to set the bond information
    !!--------------------------------------------------------------------------
    num_pairs = 0
    pair_loop1: do i = 1, num_elements 
       pair_loop2: do j = i, num_elements
          num_pairs = num_pairs + 1
          do k = 1, size(element_bond_database)
             if( this%element_info(i)%name .eq. &
                      element_bond_database(k)%element(1) .and. &
                 this%element_info(j)%name .eq. &
                      element_bond_database(k)%element(2) ) then
                this%bond_info(num_pairs) = element_bond_database(k)
                cycle pair_loop2
             elseif( this%element_info(i)%name .eq. &
                      element_bond_database(k)%element(2) .and. &
                 this%element_info(j)%name .eq. &
                      element_bond_database(k)%element(1) ) then
                this%bond_info(num_pairs) = element_bond_database(k)
                this%bond_info(num_pairs)%coordination = &
                     this%bond_info(num_pairs)%coordination(2:1:-1)
                this%bond_info(num_pairs)%element = &
                     this%bond_info(num_pairs)%element(2:1:-1)
                cycle pair_loop2
             end if
          end do
          !! check if all pairs were found
          write(0,*) 'Error reading element bond data'
          write(0,*) 'Could not find bond data for ', &
                     this%element_info(i)%name, ' and ', &
                     this%element_info(j)%name
          stop 1
       end do pair_loop2
    end do pair_loop1

  end subroutine set_bond_info
!!!#############################################################################


!!!#############################################################################
!!! set the best energy and system
!!!#############################################################################
  subroutine set_best_energy(this)
    implicit none
    class(gvector_container_type), intent(inout) :: this

    integer :: i, is, idx
    real(real12) :: energy

    do i = 1, size(this%system)
       energy = this%system(i)%energy
       do is = 1, size(this%system(i)%species)
          idx = findloc( [ this%element_info(:)%name ], &
                           this%system(i)%species(is), dim=1 )
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
!!!#############################################################################


!!!#############################################################################
!!! get index corresponding to element pair
!!!#############################################################################
  pure function get_pair_index(this, species1, species2) result(idx)
    implicit none
    class(gvector_container_type), intent(in) :: this
    character(len=3), intent(in) :: species1, species2
    integer :: idx

    integer :: is, js

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
!!!#############################################################################


!!!#############################################################################
!!! get bin number associated with input value
!!!#############################################################################
  pure function get_bin(this, value, dim) result(bin)
    implicit none
    class(gvector_container_type), intent(in) :: this
    real(real12), intent(in) :: value
    integer, intent(in) :: dim
    integer :: bin

    if(value .lt. this%cutoff_min(dim) - this%width(dim) .or. &
         value .gt. this%cutoff_max(dim) + this%width(dim))then
       bin = 0
    else
       bin = nint( ( this%nbins(dim) - 1 ) * &
                   ( value - this%cutoff_min(dim) ) / &
                   ( this%cutoff_max(dim) - this%cutoff_min(dim) ) ) + 1
    end if

  end function get_bin
!!!#############################################################################


!!!#############################################################################
!!! initialise the gvectors
!!!#############################################################################
  subroutine initialise_gvectors(this)
    implicit none
    class(gvector_container_type), intent(inout) :: this

    integer :: i, num_pairs
    real(real12) :: eta, weight, height
    !real(real12), dimension(42) :: bonds_cubic

    num_pairs = nint( gamma(real(size(this%element_info) + 2, real12)) / &
         ( gamma(real(size(this%element_info), real12)) * gamma( 3._real12 ) ) )
    allocate(this%total%df_2body(this%nbins(1),num_pairs), source = 0._real12)
    allocate(this%total%df_3body(this%nbins(2),size(this%element_info)), &
         source = 1._real12/this%nbins(2))
    allocate(this%total%df_4body(this%nbins(3),size(this%element_info)), &
         source = 1._real12/this%nbins(3))

    this%total%df_2body(:,:) = 1._real12 / this%nbins(1)
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
!!!#############################################################################


!!!#############################################################################
!!! generate the evolved gvectors
!!!#############################################################################
!!! EASIER TO STORE THE LIST OF LENGTHS, AND ANGLES, OR THE INDIVIDUAL SYSTEM GVECTORS?
  subroutine evolve(this, system, deallocate_systems_after_evolve)
    implicit none
    class(gvector_container_type), intent(inout) :: this
    type(gvector_type), dimension(..), intent(in), optional :: system
    logical, intent(in), optional :: deallocate_systems_after_evolve

    integer :: idx1, idx2
    integer :: i, j, is, js, num_structures_previous
    real(real12) :: weight, energy, best_energy_old
    logical :: deallocate_systems_after_evolve_ = .true.
    real(real12), dimension(:), allocatable :: height
    integer, dimension(:,:), allocatable :: idx_list
    

    !!--------------------------------------------------------------------------
    !! if present, set the deallocate flag
    !!--------------------------------------------------------------------------
    if(present(deallocate_systems_after_evolve)) &
       deallocate_systems_after_evolve_ = deallocate_systems_after_evolve


    !!--------------------------------------------------------------------------
    !! if present, add the system to the container
    !!--------------------------------------------------------------------------
    if(present(system)) call this%add(system)
    do i = 1, 3
       if(this%nbins(i).eq.-1)then
          this%nbins(i) = 1 + &
               ( this%cutoff_max(i) - this%cutoff_min(i) ) / this%width(i)
       end if
    end do


    !!--------------------------------------------------------------------------
    !! get the energy from the lowest formation energy system
    !!--------------------------------------------------------------------------
    best_energy_old = this%best_energy
    call this%set_best_energy()


    !!--------------------------------------------------------------------------
    !! initialise the total gvectors
    !!--------------------------------------------------------------------------
    if(.not.allocated(this%total%df_2body))then
       call this%initialise_gvectors()
    else
      this%total%df_3body = this%total%df_3body * exp( this%best_energy ) / &
                              exp( best_energy_old )
    end if


    !!--------------------------------------------------------------------------
    !! loop over all systems to calculate the total gvectors
    !!--------------------------------------------------------------------------
    do i = 1, size(this%system)
       !!-----------------------------------------------------------------------
       !! get the list of 2-body species pairs the system
       !!-----------------------------------------------------------------------
       j = 0
       allocate(idx_list(size(this%system(i)%species),&
                         size(this%system(i)%species)))
       do is = 1, size(this%system(i)%species)
          do js = is, size(this%system(i)%species), 1
             j = j + 1
             idx_list(is,js) = j
             idx_list(js,is) = j
          end do
       end do
       
       
       !!-----------------------------------------------------------------------
       !! calculate the weight for the system
       !!-----------------------------------------------------------------------
       energy = this%system(i)%energy
       do is = 1, size(this%system(i)%species)
          idx1 = findloc( [ this%element_info(:)%name ], &
                          this%system(i)%species(is), dim=1)
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

       !!-----------------------------------------------------------------------
       !! loop over all species in the system to add the gvectors
       !!-----------------------------------------------------------------------
       do is = 1, size(this%system(i)%species)

          idx1 = findloc( [ this%element_info(:)%name ], &
                          this%system(i)%species(is), dim=1)

          height = 1._real12 / ( 1._real12 + this%total%df_3body(:,idx1) )
          this%total%df_3body(:,idx1) = this%total%df_3body(:,idx1) + &
               height * weight * this%system(i)%df_3body(:,is)
          
          height = 1._real12 / ( 1._real12 + this%total%df_4body(:,idx1) )
          this%total%df_4body(:,idx1) = this%total%df_4body(:,idx1) + &
               height * weight * this%system(i)%df_4body(:,is)
          
          do js = is, size(this%system(i)%species), 1
             idx2 = findloc( [ this%element_info(:)%name ], &
                             this%system(i)%species(js), dim=1)
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


    !!--------------------------------------------------------------------------
    !! deallocate the individual gvectors
    !!--------------------------------------------------------------------------
    if(deallocate_systems_after_evolve_) deallocate(this%system)

  end subroutine evolve
!!!#############################################################################


!!!#############################################################################
!!! calculate the gvectors for a system from its basis
!!!#############################################################################
  subroutine calculate(this, basis, &
       nbins, width, sigma, cutoff_min, cutoff_max)
    implicit none
    class(gvector_type), intent(inout) :: this
    type(bas_type), intent(in) :: basis

    integer, dimension(3), intent(in), optional :: nbins
    real(real12), dimension(3), intent(in), optional :: width, sigma
    real(real12), dimension(3), intent(in), optional :: cutoff_min, cutoff_max

    integer, dimension(3) :: nbins_
    real(real12), dimension(3) :: sigma_ = [0.1_real12, 0.05_real12, 0.05_real12] !!! RANDOMLY CHOSEN DEFAULTS FOR NOW
    real(real12), dimension(3) :: width_ = [0.25_real12, pi/24._real12, pi/32._real12] !!! RANDOMLY CHOSEN DEFAULTS FOR NOW
    real(real12), dimension(3) :: cutoff_min_ = [0._real12, 0._real12, 0._real12]
    real(real12), dimension(3) :: cutoff_max_ = [6._real12, pi, pi/2._real12]

    integer :: bin, max_num_steps
    integer :: i, j, k, b, itmp1
    integer :: is, js, ia, ja
    integer :: num_pairs, num_angles
    integer :: amax, bmax, cmax
    real(real12) :: rtmp1, rtmp2, fc, weight, scale
    real(real12), dimension(3) :: eta, limit
    real(real12), dimension(3) :: vtmp1, vtmp2, vtmp3, diff
    real(real12), allocatable, dimension(:) :: gvector_tmp, angle
    integer, dimension(:), allocatable :: idx_list!, count_list

    integer, dimension(3,2) :: loop_limits
    integer, allocatable, dimension(:,:) :: idx

    type :: bond_type
       integer, dimension(2) :: species, atom
       logical :: skip = .false.
       real(real12), dimension(3) :: vector
    end type bond_type
    type(bond_type), dimension(:), allocatable :: bond_list

    !type :: plane_type
    !   real(real12), dimension(3) :: vector
    !   integer :: count
    !end type plane_type
    !type(plane_type), dimension(:), allocatable :: plane_list


    !!--------------------------------------------------------------------------
    !! initialise optional variables
    !!--------------------------------------------------------------------------
    if(present(cutoff_min)) cutoff_min_ = cutoff_min
    if(present(cutoff_max)) cutoff_max_ = cutoff_max
    if(present(width)) width_ = width
    if(present(sigma)) sigma_ = sigma
    if(present(nbins))then
       nbins_ = nbins
       width_ = ( cutoff_max_ - cutoff_min_ )/real( nbins_ - 1, real12 )
    else
       nbins_ = 1 + (cutoff_max_ - cutoff_min_)/width_
    end if
    limit = cutoff_max_ - cutoff_min_


    !!--------------------------------------------------------------------------
    !! get the number of pairs of species
    !!--------------------------------------------------------------------------
    !! combination calculator with repetition
    i = 0
    num_pairs = gamma(real(basis%nspec + 2, real12)) / &
                ( gamma(real(basis%nspec, real12)) * gamma( 3._real12 ) )
    allocate(idx(2,num_pairs))
    allocate(this%species(basis%nspec))
    do is = 1, basis%nspec
       this%species(is) = strip_null(basis%spec(is)%name)
    end do
    do is = 1, basis%nspec
       do js = is, basis%nspec, 1
          i = i + 1
          idx(:,i) = [is, js]
       end do
    end do


    !!--------------------------------------------------------------------------
    !! get the stoichiometry, energy, and number of atoms
    !!--------------------------------------------------------------------------
    this%stoichiometry = basis%spec(:)%num
    this%energy = basis%energy
    this%num_atoms = basis%natom


    !!--------------------------------------------------------------------------
    !! get the maximum number of lattice vectors to consider
    !!--------------------------------------------------------------------------
    !! this is not perfect
    !! won't work for extremely acute/obtuse angle cells
    !! (due to diagonal path being shorter than individual lattice vectors)
    amax = ceiling(cutoff_max_(1)/modu(basis%lat(1,:)))
    bmax = ceiling(cutoff_max_(1)/modu(basis%lat(2,:)))
    cmax = ceiling(cutoff_max_(1)/modu(basis%lat(3,:)))


    !!--------------------------------------------------------------------------
    !! build the bond list
    !!--------------------------------------------------------------------------
    !write(*,*) amax, bmax, cmax
    !! estimate number of bonds
    !write(*,*) "estimated number of bonds: ", triangular_number(basis%natom) * &
    !     ceiling( (pi * 4._real12/3._real12) * &
    !     (cutoff_max(1)**3 - cutoff_min(1)**3)/ get_vol(basis%lat) )

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

    !write(*,*) "number of bonds: ", size(bond_list)
    !write(*,*)
    !return


    !!--------------------------------------------------------------------------
    !! calculate the gaussian width
    !!--------------------------------------------------------------------------
    eta = 1._real12 / ( 2._real12 * sigma_**2._real12 )
    max_num_steps = ceiling( sqrt(16._real12/eta(1)) / width_(1) )

    !write(*,*) "starting 2-body gvectors"
    !!--------------------------------------------------------------------------
    !! build the 2-body gvectors (radial distribution functions)
    !!--------------------------------------------------------------------------
    allocate(this%df_2body(nbins_(1),num_pairs), source = 0._real12)
    allocate(gvector_tmp(nbins_(1)),             source = 0._real12)
    do i = 1, size(bond_list)
       if(bond_list(i)%skip) cycle
       if(abs(modu(bond_list(i)%vector)).lt.1.E-3) cycle
       is = bond_list(i)%species(1)
       js = bond_list(i)%species(2)
       rtmp1 = modu(bond_list(i)%vector)
       !!-----------------------------------------------------------------------
       !! get number of equivalent bonds
       !!-----------------------------------------------------------------------
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
       !!-----------------------------------------------------------------------
       bin = nint( ( rtmp1 - cutoff_min_(1) ) / width_(1) ) + 1
       if(bin.gt.nbins_(1).or.bin.lt.1) cycle

       fc = 0.5_real12 * cos( pi * ( rtmp1 - cutoff_min_(1) ) / limit(1) ) + &
            0.5_real12

       !!-----------------------------------------------------------------------
       !! calculate the gaussian for this bond
       !!-----------------------------------------------------------------------
       gvector_tmp = 0._real12
       loop_limits(:,1) = [ bin, min(nbins_(1), (bin + max_num_steps) ), 1 ]
       loop_limits(:,2) = [ bin - 1, max(1, bin - max_num_steps), -1 ]
       
       !! do forward and backward loops to add gaussian for larger distances
       do concurrent ( j = 1:2 )
          do concurrent ( b = &
                            loop_limits(1,j):loop_limits(2,j):loop_limits(3,j) )
             gvector_tmp(b) = gvector_tmp(b) + &
                  exp( -eta(1) * ( rtmp1 - &
                                   ( width_(1) * real(b-1, real12) + &
                                     cutoff_min_(1) ) ) ** 2._real12 )
          end do
       end do
       get_pair_index_loop: do j = 1, num_pairs
          if( is .eq. idx(1,j) .and. js .eq. idx(2,j) )then
             k = j
             exit get_pair_index_loop
          end if
       end do get_pair_index_loop
       itmp1 = count( [ ( ( bond_list(j)%species(1) .eq. is .and. &
                            bond_list(j)%species(2) .eq. js ) .or. &
                          ( bond_list(j)%species(2) .eq. is .and. &
                            bond_list(j)%species(1) .eq. js ), &
                              j = 1, size(bond_list), 1 ) ] )
       this%df_2body(:,k) = this%df_2body(:,k) + &
            gvector_tmp * scale * sqrt( eta(1) / pi ) / real(itmp1,real12) ! / width_(1)
    end do


    !write(*,*) "starting 3-body gvectors"
    !!--------------------------------------------------------------------------
    !! build the 3-body gvectors (angular distribution functions)
    !!--------------------------------------------------------------------------
    deallocate(gvector_tmp)
    allocate(this%df_3body(nbins_(2), basis%nspec), source = 0._real12)
    allocate(gvector_tmp(nbins_(2)),                source = 0._real12)
    do is = 1, basis%nspec
       num_angles = 0
       !! number of comibnations without repetitions:
       !! = n! / (n - r)! r!
       !! as r = 1, this simplifies to n
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
       num_angles = 0

       !!-----------------------------------------------------------------------
       !! loop over all bonds to find the first bond
       !!-----------------------------------------------------------------------
       do i = 1, size(bond_list)
          if( is .eq. bond_list(i)%species(1) )then
             vtmp1 = bond_list(i)%vector
             ia = bond_list(i)%atom(1)
          elseif( is .eq. bond_list(i)%species(2) )then
             vtmp1 = -bond_list(i)%vector
             ia = bond_list(i)%atom(2)
          else
             !write(0,*) "ERROR: Species not found in bond list1", is, i
             cycle
          end if
        
          !!--------------------------------------------------------------------
          !! loop over all bonds to find the second bond
          !!--------------------------------------------------------------------
          do j = i + 1, size(bond_list)
             if( is .eq. bond_list(j)%species(1) .and. &
                 ia .eq. bond_list(j)%atom(1) )then
                vtmp2 = bond_list(j)%vector
             elseif( is .eq. bond_list(j)%species(2) .and. &
                     ia .eq. bond_list(j)%atom(2) )then
                vtmp2 = -bond_list(j)%vector
             else
                !write(0,*) "ERROR: Species not found in bond list2", is, i
                cycle
             end if
 
             !!-----------------------------------------------------------------
             !! calculate the angle between the two bonds
             !!-----------------------------------------------------------------
             num_angles = num_angles + 1
             angle(num_angles) = get_angle( vtmp1, vtmp2 )

          end do
       end do
       if(num_angles.ne.size(angle))then
          write(0,*) "ERROR: Number of 3-body angles exceeds allocated array"
          write(0,'("Expected ",I0," got ",I0)') size(angle), num_angles
          stop 1
       end if
       this%df_3body(:,is) = this%df_3body(:,is) + &
            get_gvector( angle, nbins_(2), eta(2), width_(2), &
                               cutoff_min(2), &
                               limit(2) )
       deallocate(angle)
    end do

  
    !write(*,*) "starting 4-body gvectors"
    !!--------------------------------------------------------------------------
    !! build the 4-body gvectors (angular distribution functions)
    !!--------------------------------------------------------------------------
    deallocate(gvector_tmp)
    allocate(this%df_4body(nbins_(3),basis%nspec), source = 0._real12)
    allocate(gvector_tmp(nbins_(3)),               source = 0._real12)
    do is = 1, basis%nspec
       num_angles = 0
       !! number of comibnations without repetitions:
       !! = n! / (n - r)! r!
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
       !allocate(count_list(num_angles))
       num_angles = 0

       !!-----------------------------------------------------------------------
       !! loop over all bonds to find the first bond
       !!-----------------------------------------------------------------------
       do i = 1, size(bond_list)
          if(abs(modu(bond_list(i)%vector)).lt.1.E-3) cycle
 
          if( is .eq. bond_list(i)%species(1) )then
             vtmp1 = bond_list(i)%vector
             ia = bond_list(i)%atom(1)
          elseif( is .eq. bond_list(i)%species(2) )then
             vtmp1 = -bond_list(i)%vector
             ia = bond_list(i)%atom(2)
          else
             cycle
          end if
          !!! make list of indices where species is in bond_list(i)%species and atom is in bond_list(i)%atom
          !allocate(idx_list(0))
          !write(*,*) "START"
          !allocate(plane_list(0))
          allocate(idx_list(count( [ ( ( bond_list(j)%species(1) .eq. is .and. &
                                         bond_list(j)%atom(1) .eq. ia ) .or. &
                                       ( bond_list(j)%species(2) .eq. is .and. &
                                         bond_list(j)%atom(2) .eq. ia ), &
                                           j = i + 1, size(bond_list), 1 ) ] )))
          k = 0
          index_list_loop: do j = i + 1, size(bond_list)
             if(abs(modu(bond_list(j)%vector)).lt.1.E-3) cycle
             if( ( is .eq. bond_list(j)%species(1) .and. &
                   ia .eq. bond_list(j)%atom(1)  ) .or. &
                 ( is .eq. bond_list(j)%species(2) .and. &
                   ia .eq. bond_list(j)%atom(2) ) ) then
                k = k + 1
                idx_list(k) = j
             end if
             !! get a list of unique plane normal vectors
             !vtmp2 = cross( vtmp1, bond_list(j)%vector )
             !vtmp2 = vtmp2 / modu(vtmp2)
             !plane_check_loop: do k = 1, size(plane_list)
             !   if(abs(dot_product(vtmp2, plane_list(k)%vector)) .lt. 1.E-3)then
             !      plane_list(k)%count = plane_list(k)%count + 1
             !      cycle index_list_loop
             !   end if
             !end do plane_check_loop
             !plane_list = [ plane_list, plane_type(vector = vtmp2, count = 1) ]
          end do index_list_loop
          !write(*,*) "END", idx_list

          !!--------------------------------------------------------------------
          !! loop over all bonds to find the second bond
          !!--------------------------------------------------------------------
          do j = 1, size(idx_list)!size(plane_list)!size(idx_list)
 
             vtmp2 = cross(vtmp1, bond_list(idx_list(j))%vector)
             !!-----------------------------------------------------------------
             !! loop over all bonds to find the third bond
             !!-----------------------------------------------------------------
             !count_list(num_angles+1:num_angles+(size(idx_list)-j)) = plane_list(j)%count
             do k = j + 1, size(idx_list)
 
                !!--------------------------------------------------------------
                !! calculate the angle between the two bonds
                !!--------------------------------------------------------------
                !rtmp1 = get_angle( vtmp2, &
                !     bond_list(idx_list(k))%vector )
                !if(rtmp1 .gt. pi/2._real12) rtmp1 = pi - rtmp1
                !if(any(abs(rtmp1 - angle(:num_angles)) .lt. 1.E-3)) cycle
                num_angles = num_angles + 1
                !angle(num_angles) = rtmp1
                angle(num_angles) = &
                          !get_angle( plane_list(j)%vector, &
                          !           bond_list(idx_list(k))%vector )
                           get_angle( vtmp2, &
                                      bond_list(idx_list(k))%vector )

                !! angle never greater than 90, as this corresponds to ...
                !! ... perpendicular to plane
                !if(angle(num_angles) .gt. pi/2._real12) &
                !     angle(num_angles) = pi - angle(num_angles)
                angle(num_angles) = abs( anint( angle(num_angles)/pi )*pi - angle(num_angles) )
                
                !count_list(num_angles) = plane_list(j)%count

             end do
             !num_angles = num_angles + size(idx_list) - j
          end do
          deallocate(idx_list)
          !deallocate(plane_list)
       end do
       !write(*,'("Expected ",I0," got ",I0)') size(angle), num_angles
       if(num_angles.ne.size(angle))then
          write(0,*) "ERROR: Number of 4-body angles exceeds allocated array"
          write(0,'("Expected ",I0," got ",I0)') size(angle), num_angles
          stop 1
       end if
       this%df_4body(:,is) = this%df_4body(:,is) + &
            get_gvector( angle(:num_angles), nbins_(3), eta(3), width_(3), &
                               cutoff_min(3), &
                               limit(3))!, count_list(:num_angles) )
       deallocate(angle)
       !deallocate(count_list)
    end do
    !write(*,*) "finished 4-body gvectors"

  end subroutine calculate
!!!#############################################################################


!!!#############################################################################
!!! get the gvector for a distribution of values
!!!#############################################################################
  function get_gvector(vector, nbins, eta, width, cutoff_min, limit) &!, count_list) &
       result(gvector)
    implicit none
    integer, intent(in) :: nbins
    real(real12), intent(in) :: eta, width, cutoff_min, limit
    real(real12), dimension(:), intent(in) :: vector
    real(real12), dimension(nbins) :: gvector
    !integer, dimension(:), intent(in), optional :: count_list

    integer :: i, j, b, bin, max_num_steps
    !integer, dimension(:), allocatable :: scale_list
    !real(real12), dimension(nbins) :: gvector_tmp
    !real(real12), dimension(:), allocatable :: vector_copy
    integer, dimension(3,2) :: loop_limits


    ! max_num_steps = ceiling( abs( vector_copy(i) - sqrt(16._real12/eta) ) / width )
    max_num_steps = ceiling( sqrt(16._real12/eta) / width )
    gvector = 0._real12

    !!--------------------------------------------------------------------------
    !! calculate the gvector for a list of vectors
    !!--------------------------------------------------------------------------
    !vector_copy = vector
    !itmp1 = 0
    !!  order the vector list, and remove duplicates within a tolerance
    !call set(vector_copy, 1.E-3, scale_list)
    do i = 1, size(vector)
        !if( vector_copy(i) .lt. -1.E-3 ) cycle
        !!-----------------------------------------------------------------------
        !! remove duplicates
        !!-----------------------------------------------------------------------
        !scale = 1._real12
        !do j = i + 1, size(vector_copy)
        !   if(abs(vector_copy(i)-vector_copy(j)) .lt. 1.E-3 )then
        !      !vector_copy(j) = -1.E3_real12
        !      !scale = scale + 1._real12
        !      !itmp1 = itmp1 + 1
        !   end if
        !end do


       !!-----------------------------------------------------------------------
       !! get the bin closest to the value
       !!-----------------------------------------------------------------------
       bin = nint( ( vector(i) - cutoff_min ) / width ) + 1
       !if(bin.gt.nbins.or.bin.lt.1) cycle


       !!-----------------------------------------------------------------------
       !! calculate the gaussian for this bond
       !!-----------------------------------------------------------------------
       !gvector_tmp = 0._real12
       loop_limits(:,1) = [ min(nbins, bin), min(nbins, bin + max_num_steps), 1 ]
       loop_limits(:,2) = [ min(1, bin - 1), max(1, bin - max_num_steps), -1 ]


       !!-----------------------------------------------------------------------
       !! do forward and backward loops to add gaussian from its centre
       !!-----------------------------------------------------------------------
       do concurrent ( j = 1:2 )
          do concurrent ( b = loop_limits(1,j):loop_limits(2,j):loop_limits(3,j) )
             gvector(b) = gvector(b) + &
                  exp( -eta * ( vector(i) - &
                                   ( width * real(b-1, real12) + &
                                     cutoff_min ) ) ** 2._real12 )
          end do
       end do
       !if(present(count_list)) gvector_tmp = gvector_tmp * count_list(i)
       !gvector(:) = gvector(:) + gvector_tmp !* scale !real(scale_list(i), real12)
    end do
    gvector = gvector * sqrt( eta / pi ) / real(size(vector),real12) ! / width
    !write(*,*) "LOOK", size(vector_copy), itmp1

  end function get_gvector
!!!#############################################################################

end module evolver