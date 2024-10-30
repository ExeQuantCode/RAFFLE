module raffle__distribs_host
  !! Module for handling the host distribution functions.
  !!
  !! This module contains the type and procedures for handling the host
  !! distribution function. Procedures are also provided to calculate the
  !! interface energy of the host and to set the mapping of host elements to
  !! the element database.
  use raffle__constants, only: real32
  use raffle__io_utils, only: stop_program
  use raffle__geom_rw, only: basis_type
  use raffle__element_utils, only: element_type
  use raffle__distribs, only: distribs_type
  implicit none

  
  private

  public :: distribs_host_type


  type, extends(distribs_type) :: distribs_host_type
     !! Type for host information.
     !!
     !! This type contains the information regarding the host structure that
     !! will be used in the grandparent generator type.
     logical :: defined = .false.
     !! Boolean whether the host structure has been set.
     real(real32) :: interface_energy = 0.0_real32
     !! Energy associated with the formation of the interface in the host.
     type(basis_type) :: basis
     !! Host structure.
     integer, dimension(:,:), allocatable :: pair_index
     !! Index for the 2-body distribution function.
     integer, dimension(:), allocatable :: element_map
     !! Mapping of host elements to distribution function elements.
   contains
     procedure, pass(this) :: calculate_interface_energy
     !! Calculate the interface formation energy of the host.
     procedure, pass(this) :: set => set_host
     !! Set the host structure for the distribution functions.
     procedure, pass(this) :: set_element_map => set_host_element_map
     !! Set the mapping of host elements to distribution function elements.
  end type distribs_host_type


contains

!###############################################################################
  subroutine set_host(this, host)
    !! Set the host structure for the distribution functions.
    !!
    !! distribution function not needed for host
    implicit none

    ! Arguments
    class(distribs_host_type), intent(inout) :: this
    !! Parent. Instance of distribution functions container.
    type(basis_type), intent(in) :: host
    !! Host structure for the distribution functions.

    ! Local variables
    integer :: i, is, js
    !! Loop indices.

    call this%basis%copy(host)
    this%defined = .true.
    if(allocated(this%pair_index)) deallocate(this%pair_index)
    allocate(this%pair_index(this%basis%nspec, this%basis%nspec))
    i = 0
    do is = 1, this%basis%nspec
       do js = is, this%basis%nspec, 1
          i = i + 1
          this%pair_index(js,is) = i
          this%pair_index(is,js) = i
       end do
    end do
    if(allocated(this%df_2body)) deallocate(this%df_2body)
    if(allocated(this%df_3body)) deallocate(this%df_3body)
    if(allocated(this%df_4body)) deallocate(this%df_4body)

    if(allocated(this%stoichiometry)) deallocate(this%stoichiometry)
    if(allocated(this%element_symbols)) deallocate(this%element_symbols)
    if(allocated(this%num_pairs)) deallocate(this%num_pairs)
    if(allocated(this%num_per_species)) deallocate(this%num_per_species)
    if(allocated(this%weight_pair)) deallocate(this%weight_pair)
    if(allocated(this%weight_per_species)) deallocate(this%weight_per_species)

  end subroutine set_host
!###############################################################################


!###############################################################################
  subroutine calculate_interface_energy(this, element_info)
    !! Calculate the interface formation energy of the host.
    implicit none

    ! Arguments
    class(distribs_host_type), intent(inout) :: this
    !! Parent. Instance of host type.
    type(element_type), dimension(:), intent(in) :: element_info
    !! List of elements and properties.

    ! Local variables
    integer :: is, idx1
    !! Loop indices.

    this%interface_energy = this%energy
    do is = 1, size(this%element_symbols)
       idx1 = findloc( &
            [ element_info(:)%name ], &
            this%element_symbols(is), dim=1 &
       )
       if(idx1.lt.1)then
          call stop_program( "Species not found in species list" )
          return
       end if
       this%interface_energy = this%interface_energy - &
            this%stoichiometry(is) * element_info(idx1)%energy
    end do

  end subroutine calculate_interface_energy
!###############################################################################


!###############################################################################
  subroutine set_host_element_map(this, element_info)
    !! Set the host element map for the container.
    implicit none

    ! Arguments
    class(distribs_host_type), intent(inout) :: this
    !! Parent of the procedure. Instance of distribution functions container.
    type(element_type), dimension(:), intent(in) :: element_info
    !! Element information.

    ! Local variables
    integer :: is
    !! Index of the elements in the element_info array.

    if(.not.this%defined)then
       call stop_program( "Host not defined" )
       return
    end if
    if(allocated(this%element_map)) deallocate(this%element_map)
    allocate(this%element_map(this%basis%nspec))
    do is = 1, this%basis%nspec
       this%element_map(is) = findloc(&
            [ element_info(:)%name ], &
            this%basis%spec(is)%name, dim=1 &
       )
    end do

  end subroutine set_host_element_map
!###############################################################################

end module raffle__distribs_host