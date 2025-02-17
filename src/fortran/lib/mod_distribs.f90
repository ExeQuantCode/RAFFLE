module raffle__distribs
  !! Module for handling distribution functions.
  !!
  !! This module contains the types and subroutines for generating distribution
  !! fucntions for individual materials.
  !! The distribution functions are used as fingerprints for atomic structures
  !! to identify similarities and differences between structures.
  use raffle__constants, only: real32, pi
  use raffle__io_utils, only: stop_program, print_warning
  use raffle__misc, only: strip_null, sort_str
  use raffle__misc_maths, only: triangular_number
  use raffle__misc_linalg, only: get_angle, get_improper_dihedral_angle, modu
  use raffle__geom_rw, only: basis_type, get_element_properties
  use raffle__geom_extd, only: extended_basis_type
  use raffle__element_utils, only: &
       element_type, element_bond_type, &
       element_database, element_bond_database
  implicit none

  
  private

  public :: distribs_base_type, distribs_type, get_distrib
  
  
  type :: distribs_base_type
     !! Base type for distribution functions.
     real(real32), dimension(:,:), allocatable :: df_2body
     !! 2-body distribution function.
     real(real32), dimension(:,:), allocatable :: df_3body
     !! 3-body distribution function.
     real(real32), dimension(:,:), allocatable :: df_4body
     !! 4-body distribution function.
  end type distribs_base_type

  type, extends(distribs_base_type) :: distribs_type
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
     real(real32) :: energy = 0.0_real32
     !! Energy of the structure.
     real(real32) :: energy_above_hull = 0.0_real32
     !! Energy above the hull of the structure.
     logical :: from_host = .false.
     !! Boolean whether the structure is derived from the host.
     integer, dimension(:), allocatable :: stoichiometry
     !! Stoichiometry of the structure.
     character(len=3), dimension(:), allocatable :: element_symbols
     !! Elements contained within the structure.
     integer, dimension(:), allocatable :: num_pairs, num_per_species
     !! Number of pairs and number of pairs per species.
     real(real32), dimension(:), allocatable :: weight_pair, weight_per_species
     !! Weights for the 2-body and species distribution functions.
   contains
     procedure, pass(this) :: calculate
  end type distribs_type


contains

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
  subroutine calculate(this, basis, &
       nbins, width, sigma, cutoff_min, cutoff_max, radius_distance_tol)
    !! Calculate the distribution functions for the container.
    !!
    !! This procedure calculates the 2-, 3-, and 4-body distribution function 
    !! for a given atomic structure (i.e. basis).
    implicit none

    ! Arguments
    class(distribs_type), intent(inout) :: this
    !! Parent of the procedure. Instance of distribution functions container.
    type(basis_type), intent(in) :: basis
    !! Atomic structure.
    integer, dimension(3), intent(in), optional :: nbins
    !! Optional. Number of bins for the distribution functions.
    real(real32), dimension(3), intent(in), optional :: width, sigma
    !! Optional. Width and sigma for the distribution functions.
    real(real32), dimension(3), intent(in), optional :: cutoff_min, cutoff_max
    !! Optional. Cutoff minimum and maximum for the distribution functions.
    real(real32), dimension(4), intent(in), optional :: radius_distance_tol
    !! Tolerance for the distance between atoms for 3- and 4-body.

    ! Local variables
    integer, dimension(3) :: nbins_
    !! Number of bins for the distribution functions.
    real(real32), dimension(3) :: sigma_
    !! Sigma for the distribution functions.
    real(real32), dimension(3) :: width_
    !! Width of the bins for the distribution functions.
    real(real32), dimension(3) :: cutoff_min_
    !! Cutoff minimum for the distribution functions.
    real(real32), dimension(3) :: cutoff_max_
    !! Cutoff maximum for the distribution functions.
    type(element_bond_type), dimension(:), allocatable :: bond_info
    !! Bond information for radii.
    real(real32), dimension(4) :: radius_distance_tol_
    !! Tolerance for the distance between atoms for 3- and 4-body.


    integer :: i, b, itmp1, idx
    !! Loop index.
    integer :: is, js, ia, ja, ka, la
    !! Loop index.
    integer :: num_pairs
    !! Number of pairs and angles.
    real(real32) :: bondlength
    !! Temporary real variables.
    logical :: success
    !! Boolean for success.
    type(extended_basis_type) :: basis_extd
    !! Extended basis of the system.
    type(extended_basis_type) :: neighbour_basis
    !! Basis for storing neighbour data.
    real(real32), dimension(3) :: eta
    !! Parameters for the distribution functions.
    real(real32), allocatable, dimension(:) :: angle_list, bondlength_list, &
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
       cutoff_min_ = [0.5_real32, 0._real32, 0._real32]
    end if
    if(present(cutoff_max))then
       cutoff_max_ = cutoff_max
    else
       cutoff_max_ = [6._real32, pi, pi]
    end if
    if(present(width))then
       width_ = width
    else
       width_ = [0.25_real32, pi/64._real32, pi/64._real32]
    end if
    if(present(sigma))then
       sigma_ = sigma
    else
       sigma_ = [0.1_real32, 0.1_real32, 0.1_real32]
    end if
    if(present(nbins))then
       nbins_ = nbins
       width_ = ( cutoff_max_ - cutoff_min_ )/real( nbins_ - 1, real32 )
    else
       nbins_ = 1 + nint( (cutoff_max_ - cutoff_min_)/width_ )
    end if
    if(present(radius_distance_tol))then
       radius_distance_tol_ = radius_distance_tol
    else
       radius_distance_tol_ = [1.5_real32, 2.5_real32, 3._real32, 6._real32]
    end if
       


    !---------------------------------------------------------------------------
    ! get the number of pairs of species
    ! (this uses a combination calculator with repetition)
    !---------------------------------------------------------------------------
    num_pairs = nint( &
         gamma(real(basis%nspec + 2, real32)) / &
         ( gamma(real(basis%nspec, real32)) * gamma( 3._real32 ) ) &
    )
    allocate(this%element_symbols(basis%nspec))
    do is = 1, basis%nspec
       this%element_symbols(is) = strip_null(basis%spec(is)%name)
    end do
    i = 0
    allocate(bond_info(num_pairs))
    allocate(pair_index(basis%nspec,basis%nspec))
    do is = 1, basis%nspec, 1
       do js = is, basis%nspec, 1
          i = i + 1
          pair_index(js,is) = i
          pair_index(is,js) = i
          call bond_info(i)%set( &
               this%element_symbols(is), &
               this%element_symbols(js), success &
          )
          if(success) cycle
          call set_bond_radius_to_default( [ &
               this%element_symbols(is), &
               this%element_symbols(js) &
          ] )
          call bond_info(i)%set( &
               this%element_symbols(is), &
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
    eta = 1._real32 / ( 2._real32 * sigma_**2._real32 )
    allocate(this%num_pairs(num_pairs), source = 0)
    allocate(this%num_per_species(basis%nspec), source = 0)
    allocate(this%weight_pair(num_pairs), source = 0._real32)
    allocate(this%weight_per_species(basis%nspec), source = 0._real32)
    allocate(this%df_2body(nbins_(1), num_pairs), source = 0._real32)
    allocate(this%df_3body(nbins_(2), basis%nspec), source = 0._real32)
    allocate(this%df_4body(nbins_(3), basis%nspec), source = 0._real32)


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
          allocate(distance(basis_extd%natom+basis_extd%num_images))
          neighbour_basis%spec(1)%num = 0
          neighbour_basis%image_spec(1)%num = 0
          do js = 1, basis%nspec
             itmp1 = 0

             !------------------------------------------------------------------
             ! loop over all atoms inside the unit cell
             !------------------------------------------------------------------
             atom_loop: do ja = 1, basis_extd%spec(js)%num

                associate( vector =>  matmul( &
                     [ &
                          basis_extd%spec(js)%atom(ja,1:3) - &
                          basis_extd%spec(is)%atom(ia,1:3) &
                     ], basis_extd%lat ) &
                )
                   bondlength = modu( vector )
                   
                   if( &
                        bondlength .lt. cutoff_min_(1) .or. &
                        bondlength .gt. cutoff_max_(1) &
                   ) cycle atom_loop
                  
                   ! add 2-body bond to store if within tolerances for 3-body
                   ! distance
                   if( &
                        bondlength .ge. ( &
                             bond_info(pair_index(is, js))%radius_covalent * &
                             radius_distance_tol_(1) &
                        ) .and. &
                        bondlength .le. ( &
                             bond_info(pair_index(is, js))%radius_covalent * &
                             radius_distance_tol_(2) &
                        ) &
                   ) then
                      neighbour_basis%spec(1)%num = &
                           neighbour_basis%spec(1)%num + 1
                      neighbour_basis%spec(1)%atom( &
                           neighbour_basis%spec(1)%num,1:3) = vector
                   end if

                   ! add 2-body bond to store if within tolerances for 4-body
                   ! distance
                   if( &
                        bondlength .ge. ( & 
                             bond_info(pair_index(is, js))%radius_covalent * &
                             radius_distance_tol_(3) &
                        ) .and. &
                        bondlength .le. ( &
                             bond_info(pair_index(is, js))%radius_covalent * &
                             radius_distance_tol_(4) &
                        ) &
                   ) then
                      neighbour_basis%image_spec(1)%num = &
                           neighbour_basis%image_spec(1)%num + 1
                      neighbour_basis%image_spec(1)%atom( &
                           neighbour_basis%image_spec(1)%num,1:3) = vector
                   end if

                   !if(js.lt.js.or.(is.eq.js.and.ja.le.ia)) cycle
                   itmp1 = itmp1 + 1
                   bondlength_list(itmp1) = bondlength
                   distance(itmp1) = 1._real32
                
                end associate
             end do atom_loop


             !------------------------------------------------------------------
             ! loop over all image atoms outside of the unit cell
             !------------------------------------------------------------------
             image_loop: do ja = 1, basis_extd%image_spec(js)%num
                associate( vector =>  matmul( &
                     [ &
                          basis_extd%image_spec(js)%atom(ja,1:3) - &
                          basis_extd%spec(is)%atom(ia,1:3) &
                     ], basis_extd%lat ) &
                )

                   bondlength = modu( vector )
                   
                   if( &
                        bondlength .lt. cutoff_min_(1) .or. &
                        bondlength .gt. cutoff_max_(1) &
                   ) cycle image_loop
                  
                   ! add 2-body bond to store if within tolerances for 3-body
                   ! distance
                   if( &
                        bondlength .ge. ( &
                             bond_info(pair_index(is, js))%radius_covalent * &
                             radius_distance_tol_(1) &
                        ) .and. &
                        bondlength .le. ( &
                             bond_info(pair_index(is, js))%radius_covalent * &
                             radius_distance_tol_(2) &
                        ) &
                   ) then
                      neighbour_basis%spec(1)%num = &
                           neighbour_basis%spec(1)%num + 1
                      neighbour_basis%spec(1)%atom( &
                           neighbour_basis%spec(1)%num,1:3 &
                      ) = vector
                   end if

                   ! add 2-body bond to store if within tolerances for 4-body
                   ! distance
                   if( &
                        bondlength .ge. ( & 
                             bond_info(pair_index(is, js))%radius_covalent * &
                             radius_distance_tol_(3) &
                        ) .and. &
                        bondlength .le. ( &
                             bond_info(pair_index(is, js))%radius_covalent * &
                             radius_distance_tol_(4) &
                        ) &
                   ) then
                      neighbour_basis%image_spec(1)%num = &
                           neighbour_basis%image_spec(1)%num + 1
                      neighbour_basis%image_spec(1)%atom( &
                           neighbour_basis%image_spec(1)%num,1:3 &
                      ) = vector
                   end if

                   itmp1 = itmp1 + 1
                   bondlength_list(itmp1) = bondlength
                   distance(itmp1) = 1._real32
                
                end associate
             end do image_loop

             !------------------------------------------------------------------
             ! calculate the 2-body distribution function contributions from
             ! atom (is,ia) for species pair (is,js)
             !------------------------------------------------------------------
             if(itmp1.gt.0)then
                this%df_2body(:,pair_index(is, js)) = &
                     this%df_2body(:,pair_index(is, js)) + &
                     get_distrib( &
                          bondlength_list(:itmp1), &
                          nbins_(1), eta(1), width_(1), &
                          cutoff_min_(1), &
                          scale_list = distance(:itmp1) &
                     )
                this%weight_pair(pair_index(is, js)) = &
                     this%weight_pair(pair_index(is, js)) + &
                     4._real32 * sum( &
                          ( &
                               bond_info(pair_index(is, js))%radius_covalent / &
                               bondlength_list(:itmp1) ) ** 2 &
                     )
                this%num_pairs(pair_index(is, js)) = &
                     this%num_pairs(pair_index(is, js)) + itmp1
                this%weight_per_species(is) = &
                     this%weight_per_species(is) + &
                     4._real32 * sum( &
                          ( &
                               bond_info(pair_index(is, js))%radius_covalent / &
                               bondlength_list(:itmp1) ) ** 2 &
                     )
                this%num_per_species(is) = this%num_per_species(is) + itmp1
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
          ! a NaN in the angle refers to one where two of the vectors are
          ! parallel, so the angle is undefined
          do i = 1, size(angle_list)
             if(isnan(angle_list(i)))then
                angle_list(i) = -huge(1._real32)
                distance(i) =  1._real32
             end if
          end do
          this%df_3body(:,is) = this%df_3body(:,is) + &
               get_distrib( &
                    angle_list, &
                    nbins_(2), eta(2), width_(2), &
                    cutoff_min_(2), &
                    scale_list = distance &
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
          ! a NaN in the angle refers to one where two of the vectors are
          ! parallel, so the angle is undefined
          do i = 1, size(angle_list)
             if(isnan(angle_list(i)))then
                angle_list(i) = -huge(1._real32)
                distance(i) =  1._real32
             end if
          end do
          this%df_4body(:,is) = this%df_4body(:,is) + &
               get_distrib( &
                    angle_list, &
                    nbins_(3), eta(3), width_(3), &
                    cutoff_min_(3), &
                    scale_list = distance &
               )
          deallocate( angle_list, distance )

       end do
    end do

    !---------------------------------------------------------------------------
    ! apply the cutoff function to the 2-body distribution function
    !---------------------------------------------------------------------------
    do b = 1, nbins_(1)
       this%df_2body(b,:) = this%df_2body(b,:) / ( cutoff_min_(1) + &
            width_(1) * real(b-1, real32) ) ** 2
    end do


    !---------------------------------------------------------------------------
    ! renormalise the distribution functions so that area under the curve is 1
    !---------------------------------------------------------------------------
    do i = 1, num_pairs
       if(any(abs(this%df_2body(:,i)).gt.1.E-6))then
          this%df_2body(:,i) = this%df_2body(:,i) / sum(this%df_2body(:,i))
       end if
    end do
    do is = 1, basis%nspec
       if(any(abs(this%df_3body(:,is)).gt.1.E-6))then
          this%df_3body(:,is) = this%df_3body(:,is) / sum(this%df_3body(:,is))
       end if
       if(any(abs(this%df_4body(:,is)).gt.1.E-6))then
          this%df_4body(:,is) = this%df_4body(:,is) / sum(this%df_4body(:,is))
       end if
    end do


    !---------------------------------------------------------------------------
    ! check for NaN in the distribution functions
    !---------------------------------------------------------------------------
    if(any(isnan(this%df_2body)))then
       call stop_program('NaN in 2-body distribution function')
    end if
    if(any(isnan(this%df_3body)))then
       call stop_program('NaN in 3-body distribution function')
    end if
    if(any(isnan(this%df_4body)))then
       call stop_program('NaN in 4-body distribution function')
    end if

  end subroutine calculate
!###############################################################################


!###############################################################################
  function get_distrib(value_list, nbins, eta, width, cutoff_min, &
       scale_list ) result(output)
    !! Calculate the angular distribution function for a list of values.
    implicit none

    ! Arguments
    integer, intent(in) :: nbins
    !! Number of bins for the distribution functions.
    real(real32), intent(in) :: eta, width, cutoff_min
    !! Parameters for the distribution functions.
    real(real32), dimension(:), intent(in) :: value_list
    !! List of angles.
    real(real32), dimension(:), intent(in) :: scale_list
    !! List of scaling for each angle (distance**3 or distance**4)
    real(real32), dimension(nbins) :: output
    !! Distribution function for the list of values.

    ! Local variables
    integer :: i, j, b, bin
    !! Loop index.
    integer :: max_num_steps
    !! Maximum number of steps.
    integer, dimension(3,2) :: loop_limits
    !! Loop limits for the 3-body distribution function.


    max_num_steps = ceiling( sqrt(16._real32/eta) / width )
    output = 0._real32

    !---------------------------------------------------------------------------
    ! calculate the distribution function for a list of values
    !---------------------------------------------------------------------------
    do i = 1, size(value_list), 1

       !------------------------------------------------------------------------
       ! get the bin closest to the value
       !------------------------------------------------------------------------
       bin = nint( ( value_list(i) - cutoff_min ) / width ) + 1
       if( &
            bin .lt. 1 - max_num_steps .or. &
            bin .gt. nbins + max_num_steps &
       ) cycle


       !------------------------------------------------------------------------
       ! calculate the gaussian for this bond
       !------------------------------------------------------------------------
       loop_limits(:,1) = &
            [ min(nbins, bin), min(nbins, bin + max_num_steps), 1 ]
       loop_limits(:,2) = &
            [ max(0, bin - 1), max(1, bin - max_num_steps), -1 ]


       !------------------------------------------------------------------------
       ! do forward and backward loops to add gaussian from its centre
       !------------------------------------------------------------------------
       do concurrent ( j = 1:2 )
          do concurrent ( &
               b = loop_limits(1,j):loop_limits(2,j):loop_limits(3,j) &
          )
             output(b) = output(b) + &
                  exp( &
                       -eta * ( &
                            value_list(i) - &
                            ( width * real(b-1, real32) + cutoff_min ) &
                       ) ** 2._real32 &
                  ) / scale_list(i)
          end do
       end do
    end do
    output = output * sqrt( eta / pi ) / real(size(value_list,1),real32)

  end function get_distrib
!###############################################################################

end module raffle__distribs