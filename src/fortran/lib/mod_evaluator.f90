module evaluator
  !! Module to build viability map of a structure
  !!
  !! This module handles the viability map for a structure, which is a map of
  !! the system with each point in the map representing the suitability of
  !! that point for a new atom. The map is built by checking the bond lengths,
  !! bond angles and dihedral angles between the test point and all atoms.
  use constants, only: real12, pi
  use misc_linalg, only: modu, get_distance, get_angle, get_dihedral_angle, &
       get_improper_dihedral_angle
  use rw_geom, only: basis_type
  use extended_geom, only: extended_basis_type
  use edit_geom, only: get_min_dist_between_point_and_atom
  use evolver, only: gvector_container_type
  implicit none


  private
  public :: evaluate_point


contains

!###############################################################################
  function evaluate_point( gvector_container, &
       position, species, basis, atom_ignore_list, &
       radius_list ) &
       result(output)
    !! Build a map of basis and returns the value of the map at a given point
    implicit none

    ! Arguments
    integer, intent(in) :: species
    !! Index of the query element.
    type(gvector_container_type), intent(in) :: gvector_container
    !! Distribution function (gvector) container.
    type(extended_basis_type), intent(in) :: basis
    !! Basis of the system.
    real(real12), dimension(3), intent(in) :: position
    !! Position of the test point.
    integer, dimension(:,:), intent(in) :: atom_ignore_list
    !! List of atoms to ignore (i.e. indices of atoms not yet placed).
    real(real12), dimension(:), intent(in) :: radius_list
    !! List of radii for each pair of elements.
    real(real12) :: output
    !! Suitability of the test point.
    
    ! Local variables
    integer :: i, is, js, ia, ja
    !! Loop counters.
    integer :: num_2body, num_3body, num_4body
    !! Number of 2-, 3- and 4-body interactions.
    real(real12) :: contribution
    !! Contribution to the viability map
    real(real12) :: viability_2body
    !! Viability of the test point for 2-body interactions.
    real(real12) :: viability_3body, viability_4body
    !! Viability of the test point for 3- and 4-body interactions.
    real(real12) :: bondlength
    integer, dimension(:,:), allocatable :: pair_index
    !! Index of element pairs.
    type(extended_basis_type) :: neighbour_basis
    !! Basis of the neighbouring atoms.
 
    
    ! Initialisation
    output = 0._real12
    viability_2body = 0._real12


    !---------------------------------------------------------------------------
    ! get list of element pair indices
    !---------------------------------------------------------------------------
    allocate(pair_index(basis%nspec, basis%nspec), source = 0)
    do is = 1, basis%nspec
       do js = 1, basis%nspec
          pair_index(is, js) = gvector_container%get_pair_index( &
               basis%spec(is)%name, basis%spec(js)%name )
       end do
    end do


    !---------------------------------------------------------------------------
    ! loop over all atoms in the system
    !---------------------------------------------------------------------------
    allocate(neighbour_basis%spec(basis%nspec))
    allocate(neighbour_basis%image_spec(basis%nspec))
    neighbour_basis%nspec = basis%nspec
    neighbour_basis%lat = basis%lat
    num_2body = 0
    species_loop: do is = 1, basis%nspec
       ! 2-body map
       ! check bondlength between test point and all other atoms
       !------------------------------------------------------------------------
       allocate(neighbour_basis%spec(is)%atom( &
            basis%spec(is)%num+basis%image_spec(is)%num, &
            size(basis%spec(is)%atom,2) &
       ) )
       allocate(neighbour_basis%image_spec(is)%atom( &
            basis%spec(is)%num+basis%image_spec(is)%num, &
            size(basis%spec(is)%atom,2) &
       ) )
       neighbour_basis%spec(is)%num = 0
       neighbour_basis%image_spec(is)%num = 0
       atom_loop: do ia = 1, basis%spec(is)%num
          ! Check if the atom is in the ignore list
          ! If it is, skip the atom.
          do i = 1, size(atom_ignore_list,dim=1), 1
             if(all(atom_ignore_list(i,:).eq.[is,ia])) cycle atom_loop
          end do
          associate( position_store => [ basis%spec(is)%atom(ia,1:3) ] )
             bondlength = modu( matmul(position - position_store, basis%lat) )
             if( bondlength .gt. gvector_container%cutoff_max(1) ) &
                 cycle atom_loop
             if( bondlength .lt. ( &
                  radius_list(pair_index(species,is)) * &
                  gvector_container%radius_distance_tol(1) ) &
             )then
                ! If the bond length is less than the minimum allowed bond,
                ! return 0 (i.e. the point is not viable).
                return
             elseif( bondlength .le. ( &
                    radius_list(pair_index(species,is)) * &
                    gvector_container%radius_distance_tol(2) ) &
             )then
                ! If the bond length is within the tolerance of the covalent
                ! radius of the pair, add the atom to the list of
                ! atoms to consider for 3-body interactions.
                neighbour_basis%spec(is)%num = neighbour_basis%spec(is)%num + 1
                neighbour_basis%spec(is)%atom( &
                     neighbour_basis%spec(is)%num,:3 &
                ) = matmul(position_store, basis%lat)
             end if

             if( bondlength .ge. ( &
                         radius_list(pair_index(species,is)) * &
                         gvector_container%radius_distance_tol(3) &
                    ) .and. &
                    bondlength .le. min( &
                         gvector_container%cutoff_max(1), &
                         ( &
                              radius_list(pair_index(species,is)) * &
                              gvector_container%radius_distance_tol(4) &
                         ) &
                    ) &
             )then
                ! If the bond length is within the min and max allowed size,
                ! add the atom to the list of atoms to consider for 4-body.
                neighbour_basis%image_spec(is)%num = &
                     neighbour_basis%image_spec(is)%num + 1
                neighbour_basis%image_spec(is)%atom( &
                     neighbour_basis%image_spec(is)%num,:3 &
                ) = matmul(position_store, basis%lat)
             end if
        
             ! Add the contribution of the bond length to the viability map.
             viability_2body = viability_2body + &
                  gvector_container%total%df_2body( &
                       gvector_container%get_bin(bondlength, dim = 1), &
                       pair_index(species,is) &
                  )
             num_2body = num_2body + 1
          end associate
       end do atom_loop

       ! Repeat the process for the image atoms.
       ! i.e. atoms that are not in the unit cell but are within the cutoff
       ! distance.
       image_loop: do ia = 1, basis%image_spec(is)%num, 1
          associate( position_store => [ basis%image_spec(is)%atom(ia,1:3) ] )
             bondlength = modu( matmul(position - position_store, basis%lat) )
             if( bondlength .gt. gvector_container%cutoff_max(1) ) &
                  cycle image_loop
             if( bondlength .lt. ( &
                  radius_list(pair_index(species,is)) * &
                  gvector_container%radius_distance_tol(1) ) &
             )then
                return
             elseif( bondlength .le. ( &
                    radius_list(pair_index(species,is)) * &
                    gvector_container%radius_distance_tol(2) ) &
             )then
                neighbour_basis%spec(is)%num = neighbour_basis%spec(is)%num + 1
                neighbour_basis%spec(is)%atom( &
                     neighbour_basis%spec(is)%num,:3 &
                ) = matmul(position_store, basis%lat)
             end if

             if( bondlength .ge. ( &
                         radius_list(pair_index(species,is)) * &
                         gvector_container%radius_distance_tol(3) &
                    ) .and. &
                    bondlength .le. min( &
                         gvector_container%cutoff_max(1), &
                         ( &
                              radius_list(pair_index(species,is)) * &
                              gvector_container%radius_distance_tol(4) &
                         ) &
                    ) &
             )then
                neighbour_basis%image_spec(is)%num = &
                     neighbour_basis%image_spec(is)%num + 1
                neighbour_basis%image_spec(is)%atom( &
                     neighbour_basis%image_spec(is)%num,:3 &
                ) = matmul(position_store, basis%lat)
             end if
        
             viability_2body = viability_2body + &
                  gvector_container%total%df_2body( &
                       gvector_container%get_bin(bondlength, dim = 1), &
                       pair_index(species,is) &
                  )
             num_2body = num_2body + 1
          end associate
       end do image_loop
    end do species_loop
    neighbour_basis%natom = sum(neighbour_basis%spec(:)%num)
    ! Normalise the viability map
    if(num_2body.eq.0)then
       ! This does not matter as, if there are no 2-body bonds, the point is
       ! not meant to be included in the viability set.
       ! The evaluator cannot comment on the viability of the point.
       viability_2body = 0.5_real12
    else
       viability_2body = viability_2body / real( num_2body, real12 )
    end if



    ! store 3-body viable atoms in neighbour_basis%spec
    ! store 4-body viable atoms in neighbour_basis%image_spec
    ! for 3-bdoy, just cycle over neighbour_basis%spec
    ! for 4-body, cycle over neighbour_basis%spec for first two atoms,
    ! and neighbour_basis%image_spec for the third atom
    num_3body = 0
    num_4body = 0
    viability_3body = 1._real12
    viability_4body = 1._real12
    do is = 1, neighbour_basis%nspec
       do ia = 1, neighbour_basis%spec(is)%num
          ! 3-body map
          ! check bondangle between test point and all other atoms
          !---------------------------------------------------------------------
          associate( &
               position_1 => matmul(position, basis%lat), &
               position_2 => [neighbour_basis%spec(is)%atom(ia,1:3)] &
          )
             num_3body = num_3body + 1
             viability_3body = viability_3body * &
                   evaluate_3body_contributions( gvector_container, &
                      position_1, &
                      position_2, &
                      neighbour_basis, species, [is, ia], num_3body &
                   )
             do js = is, neighbour_basis%nspec
            !  do js = 1, neighbour_basis%nspec
                do ja = 1, neighbour_basis%spec(js)%num
                   if(js.eq.is .and. ja.le.ia) cycle
                   if(all(neighbour_basis%image_spec(:)%num.eq.0))cycle
                   num_4body = num_4body + 1
                   ! 4-body map
                   ! check improperdihedral angle between test point and all
                   ! other atoms
                   !------------------------------------------------------------
                   viability_4body = viability_4body * &
                        evaluate_4body_contributions( gvector_container, &
                           position_1, &
                           position_2, &
                           [neighbour_basis%spec(js)%atom(ja,1:3)], &
                           neighbour_basis, species &
                        )
                end do
             end do
          end associate
       end do
    end do
    ! Normalise the viability map
    if(num_3body.eq.0)then
       viability_3body = gvector_container%viability_3body_default
    else
       viability_3body = viability_3body ** ( &
            1._real12 / real(num_3body,real12) &
       )
    end if
    if(num_4body.eq.0)then
       viability_4body = gvector_container%viability_4body_default
    else
       viability_4body = viability_4body ** ( &
            1._real12 / real(num_4body,real12) &
       )
    end if
    
    ! Combine the 2-, 3- and 4-body maps
    output = viability_2body * viability_3body * viability_4body
    
  end function evaluate_point
!###############################################################################


!###############################################################################
  function evaluate_3body_contributions( gvector_container, &
       position_1, position_2, basis, species, current_idx, num_3body &
  ) result(output)
    !! Return the contribution to the viability map from 3-body interactions
    implicit none

    ! Arguments
    type(gvector_container_type), intent(in) :: gvector_container
    !! Distribution function (gvector) container.
    real(real12), dimension(3), intent(in) :: position_1, position_2
    !! Positions of the atoms.
    type(extended_basis_type), intent(in) :: basis
    !! Basis of the system.
    integer, intent(in) :: species
    !! Index of the query element.
    integer, dimension(2), intent(in) :: current_idx
    !! Index of the 1st-atom query element.
    integer, intent(inout) :: num_3body
    !! Number of 3-body interactions.
    real(real12) :: output
    !! Contribution to the viability map.

    ! Local variables
    integer :: js, ja
    !! Loop indices.
    integer :: bin
    !! Bin for the distribution function.
    integer :: num_3body_local
    !! Number of 3-body interactions local to the current atom pair.


    output = 1._real12
    num_3body_local = 0
    species_loop: do js = current_idx(1), basis%nspec, 1
       atom_loop: do ja = 1, basis%spec(js)%num
          if(all([js,ja].eq.current_idx))cycle
          num_3body_local = num_3body_local + 1
          associate( position_store => [ basis%spec(js)%atom(ja,1:3) ] )
             bin = gvector_container%get_bin( &
                  get_angle( position_2, &
                             position_1, &
                             position_store ), &
                  dim = 2 &
             )
             !!! THIS IS A TEMPORARY CHECK
             !!! IF IT IS NEVER ENCOUNTERED, WE CAN REMOVE IT
             !!! AND WHEN REMOVING, WE CAN MAKE ALL PROCEDURES HERE PURE
             if(bin.eq.0)then
                write(0,*) "Error: bin = 0, IF NOT TRIGGERED, WE CAN REMOVE THIS IF"
                stop 1
             end if
             output = output * gvector_container%total%df_3body(bin,species)
          end associate
       end do atom_loop
    end do species_loop
    if(num_3body_local.eq.0)then
       output = 1._real12
       num_3body = num_3body - 1
    else
       output = output ** ( &
            1._real12 / real( num_3body_local, real12 ) &
       )
    end if

  end function evaluate_3body_contributions
!###############################################################################


!###############################################################################
  function evaluate_4body_contributions( gvector_container, &
       position_1, position_2, position_3, basis, species ) result(output)
    !! Return the contribution to the viability map from 4-body interactions
    implicit none

    ! Arguments
    type(gvector_container_type), intent(in) :: gvector_container
    !! Distribution function (gvector) container.
    real(real12), dimension(3), intent(in) :: position_1, position_2, position_3
    !! Positions of the atoms.
    type(extended_basis_type), intent(in) :: basis
    !! Basis of the system.
    integer, intent(in) :: species
    !! Index of the query element.
    real(real12) :: output
    !! Contribution to the viability map.

    ! Local variables
    integer :: ks, ka
    !! Loop indices.
    integer :: bin
    !! Bin for the distribution function.
    real(real12) :: contribution
    !! Contribution to the viability map.
    integer :: num_4body_local
    !! Number of 4-body interactions local to the current atom triplet.


    output = 1._real12
    num_4body_local = sum(basis%image_spec(:)%num)
    if(num_4body_local.eq.0) return
    species_loop: do ks = 1, basis%nspec, 1
       atom_loop: do ka = 1, basis%image_spec(ks)%num
          associate( position_store => [ basis%image_spec(ks)%atom(ka,1:3) ] )
             bin = gvector_container%get_bin( &
                  get_improper_dihedral_angle( &
                                 position_1, &
                                 position_2, &
                                 position_3, &
                                 position_store &
                  ), &
                  dim = 3 &
             )
             !!! THIS IS A TEMPORARY CHECK
             !!! IF IT IS NEVER ENCOUNTERED, WE CAN REMOVE IT
             !!! AND WHEN REMOVING, WE CAN MAKE ALL PROCEDURES HERE PURE
             if(bin.eq.0)then
                write(0,*) "Error: bin = 0, IF NOT TRIGGERED, WE CAN REMOVE THIS IF"
                stop 1
             end if
             output = output * &
                  gvector_container%total%df_4body(bin,species) ** ( &
                       1._real12 / real( num_4body_local, real12 ) &
                  )
          end associate
       end do atom_loop
    end do species_loop

  end function evaluate_4body_contributions
!###############################################################################

end module evaluator