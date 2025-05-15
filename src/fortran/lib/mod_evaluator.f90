module raffle__evaluator
  !! Module to build viability map of a structure
  !!
  !! This module handles the viability map for a structure, which is a map of
  !! the system with each point in the map representing the suitability of
  !! that point for a new atom. The map is built by checking the bond lengths,
  !! bond angles and dihedral angles between the test point and all atoms.
  use raffle__constants, only: real32, tau
  use raffle__misc_linalg, only: modu, get_angle, get_improper_dihedral_angle
  use raffle__geom_extd, only: extended_basis_type
  use raffle__distribs_container, only: distribs_container_type
  implicit none


  private
  public :: evaluate_point


contains

!###############################################################################
  pure function evaluate_point( distribs_container, &
       position, species, basis, atom_ignore_list, &
       radius_list &
  ) result(output)
    !! Return the viability of a point in a basis for a specified species
    !!
    !! This function evaluates the viability of a point in a basis for a
    !! specified species. The viability is determined by the bond lengths,
    !! bond angles and dihedral angles between the test point and all atoms.
    implicit none

    ! Arguments
    integer, intent(in) :: species
    !! Index of the query element.
    type(distribs_container_type), intent(in) :: distribs_container
    !! Distribution function (gvector) container.
    type(extended_basis_type), intent(in) :: basis
    !! Basis of the system.
    real(real32), dimension(3), intent(in) :: position
    !! Position of the test point.
    integer, dimension(:,:), intent(in) :: atom_ignore_list
    !! List of atoms to ignore (i.e. indices of atoms not yet placed).
    real(real32), dimension(:), intent(in) :: radius_list
    !! List of radii for each pair of elements.
    real(real32) :: output
    !! Suitability of the test point.

    ! Local variables
    integer :: i, is, js, ia, ja
    !! Loop counters.
    integer :: num_2body, num_3body, num_4body
    !! Number of 2-, 3- and 4-body interactions.
    real(real32) :: viability_2body
    !! Viability of the test point for 2-body interactions.
    real(real32) :: viability_3body, viability_4body
    !! Viability of the test point for 3- and 4-body interactions.
    real(real32) :: bondlength, rtmp1
    !! Temporary variables.
    real(real32), dimension(2) :: cos_scales
    !! Cosine scales for the 3- and 4-body interactions.
    real(real32), dimension(4) :: tolerances
    !! Tolerance for the distance between atoms for 3- and 4-body.
    integer, dimension(:,:), allocatable :: pair_index
    !! Index of element pairs.
    type(extended_basis_type) :: neighbour_basis
    !! Basis of the neighbouring atoms.


    ! Initialisation
    output = 0._real32
    viability_2body = 0._real32


    !---------------------------------------------------------------------------
    ! get list of element pair indices
    ! (i.e. the index for bond_info for each element pair)
    !---------------------------------------------------------------------------
    allocate(pair_index(basis%nspec, basis%nspec), source = 0)
    do is = 1, basis%nspec
       do js = 1, basis%nspec
          pair_index(is, js) = distribs_container%get_pair_index( &
               basis%spec(is)%name, basis%spec(js)%name &
          )
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
       allocate(neighbour_basis%spec(is)%atom( &
            basis%spec(is)%num+basis%image_spec(is)%num, 4 &
       ) )
       allocate(neighbour_basis%image_spec(is)%atom( &
            basis%spec(is)%num+basis%image_spec(is)%num, 4 &
       ) )
       neighbour_basis%spec(is)%num = 0
       neighbour_basis%image_spec(is)%num = 0
       tolerances = distribs_container%radius_distance_tol * &
            radius_list(pair_index(species,is))
       tolerances(2) = min( distribs_container%cutoff_max(1), tolerances(2) )
       tolerances(4) = min( distribs_container%cutoff_max(1), tolerances(4) )
       cos_scales(1) = tau / ( tolerances(2) - tolerances(1) )
       cos_scales(2) = tau / ( tolerances(4) - tolerances(3) )
       !------------------------------------------------------------------------
       ! 2-body map
       ! check bondlength between test point and all other atoms
       !------------------------------------------------------------------------
       atom_loop: do ia = 1, basis%spec(is)%num
          ! Check if the atom is in the ignore list
          ! If it is, skip the atom.
          do i = 1, size(atom_ignore_list,dim=2), 1
             if(all(atom_ignore_list(:,i).eq.[is,ia])) cycle atom_loop
          end do
          associate( position_store => [ basis%spec(is)%atom(ia,1:3) ] )
             bondlength = modu( matmul(position - position_store, basis%lat) )
             if( bondlength .gt. distribs_container%cutoff_max(1) ) &
                  cycle atom_loop
             if( bondlength .lt. tolerances(1) )then
                ! If the bond length is less than the minimum allowed bond,
                ! return 0 (i.e. the point is not viable).
                return
             elseif( bondlength .le. tolerances(2) )then
                ! If the bond length is within the tolerance of the covalent
                ! radius of the pair, add the atom to the list of
                ! atoms to consider for 3-body interactions.
                neighbour_basis%spec(is)%num = neighbour_basis%spec(is)%num + 1
                neighbour_basis%spec(is)%atom( &
                     neighbour_basis%spec(is)%num,:3 &
                ) = matmul(position_store, basis%lat)
                neighbour_basis%spec(is)%atom( &
                     neighbour_basis%spec(is)%num,4 &
                ) = 0.5_real32 * ( 1._real32 - &
                     cos( cos_scales(1) * ( bondlength - tolerances(1) ) ) )
             end if

             if( bondlength .ge. tolerances(3) .and. &
                  bondlength .le. tolerances(4) &
             )then
                ! If the bond length is within the min and max allowed size,
                ! add the atom to the list of atoms to consider for 4-body.
                neighbour_basis%image_spec(is)%num = &
                     neighbour_basis%image_spec(is)%num + 1
                neighbour_basis%image_spec(is)%atom( &
                     neighbour_basis%image_spec(is)%num,:3 &
                ) = matmul(position_store, basis%lat)
                neighbour_basis%image_spec(is)%atom( &
                     neighbour_basis%image_spec(is)%num,4 &
                ) = 0.5_real32 * ( 1._real32 - &
                     cos( cos_scales(2) * ( bondlength - tolerances(3) ) ) )
             end if

             !------------------------------------------------------------------
             ! Add the contribution of the bond length to the viability
             !------------------------------------------------------------------
             viability_2body = viability_2body + &
                  evaluate_2body_contributions( &
                       distribs_container, bondlength, pair_index(species,is) &
                  )
             num_2body = num_2body + 1
          end associate
       end do atom_loop

       !------------------------------------------------------------------------
       ! Repeat the process for the image atoms.
       ! i.e. atoms that are not in the unit cell but are within the cutoff
       ! distance.
       !------------------------------------------------------------------------
       image_loop: do ia = 1, basis%image_spec(is)%num, 1
          associate( position_store => [ basis%image_spec(is)%atom(ia,1:3) ] )
             bondlength = modu( matmul(position - position_store, basis%lat) )
             if( bondlength .gt. distribs_container%cutoff_max(1) ) &
                  cycle image_loop
             if( bondlength .lt. tolerances(1) )then
                return
             elseif( bondlength .le. tolerances(2) )then
                neighbour_basis%spec(is)%num = neighbour_basis%spec(is)%num + 1
                neighbour_basis%spec(is)%atom( &
                     neighbour_basis%spec(is)%num,:3 &
                ) = matmul(position_store, basis%lat)
                neighbour_basis%spec(is)%atom( &
                     neighbour_basis%spec(is)%num,4 &
                ) = 0.5_real32 * ( 1._real32 - &
                     cos( cos_scales(1) * ( bondlength - tolerances(1) ) ) )
             end if

             if( bondlength .ge. tolerances(3) .and. &
                  bondlength .le. tolerances(4) &
             )then
                neighbour_basis%image_spec(is)%num = &
                     neighbour_basis%image_spec(is)%num + 1
                neighbour_basis%image_spec(is)%atom( &
                     neighbour_basis%image_spec(is)%num,:3 &
                ) = matmul(position_store, basis%lat)
                neighbour_basis%image_spec(is)%atom( &
                     neighbour_basis%image_spec(is)%num,4 &
                ) =  0.5_real32 * ( 1._real32 - &
                     cos( cos_scales(2) * ( bondlength - tolerances(3) ) ) )
             end if

             !------------------------------------------------------------------
             ! Add the contribution of the bond length to the viability
             !------------------------------------------------------------------
             viability_2body = viability_2body + &
                  evaluate_2body_contributions( &
                       distribs_container, bondlength, pair_index(species,is) &
                  )
             num_2body = num_2body + 1
          end associate
       end do image_loop
    end do species_loop
    neighbour_basis%natom = sum(neighbour_basis%spec(:)%num)


    !---------------------------------------------------------------------------
    ! Normalise the bond length viability
    !---------------------------------------------------------------------------
    if(num_2body.eq.0)then
       ! This does not matter as, if there are no 2-body bonds, the point is
       ! not meant to be included in the viability set.
       ! The evaluator cannot comment on the viability of the point.
       viability_2body = 0.5_real32
    else
       viability_2body = viability_2body / real( num_2body, real32 )
    end if


    ! store 3-body viable atoms in neighbour_basis%spec
    ! store 4-body viable atoms in neighbour_basis%image_spec
    ! for 3-bdoy, just cycle over neighbour_basis%spec
    ! for 4-body, cycle over neighbour_basis%spec for first two atoms,
    ! and neighbour_basis%image_spec for the third atom
    num_3body = 0
    num_4body = 0
    viability_3body = 1._real32
    viability_4body = 1._real32
    do is = 1, neighbour_basis%nspec
       do ia = 1, neighbour_basis%spec(is)%num
          !---------------------------------------------------------------------
          ! 3-body map
          ! check bondangle between test point and all other atoms
          !---------------------------------------------------------------------
          associate( &
               position_1 => matmul(position, basis%lat), &
               position_2 => [neighbour_basis%spec(is)%atom(ia,1:4)] &
          )
             if(sum(basis%spec(is:)%num).eq.ia) cycle
             rtmp1 = evaluate_3body_contributions( distribs_container, &
                  position_1, &
                  position_2, &
                  neighbour_basis, species, [is, ia] &
             )
             if(rtmp1.lt.-999._real32) cycle
             num_3body = num_3body + 1
             viability_3body = viability_3body * rtmp1
             do js = is, neighbour_basis%nspec
                do ja = 1, neighbour_basis%spec(js)%num
                   if(js.eq.is .and. ja.le.ia) cycle
                   if(all(neighbour_basis%image_spec(:)%num.eq.0))cycle
                   !------------------------------------------------------------
                   ! 4-body map
                   ! check improperdihedral angle between test point and all
                   ! other atoms
                   !------------------------------------------------------------
                   rtmp1 = evaluate_4body_contributions( distribs_container, &
                        position_1, &
                        position_2, &
                        [neighbour_basis%spec(js)%atom(ja,1:4)], &
                        neighbour_basis, species &
                   )
                   if(rtmp1.lt.-999._real32) cycle
                   num_4body = num_4body + 1
                   viability_4body = viability_4body * rtmp1
                end do
             end do
          end associate
       end do
    end do


    !---------------------------------------------------------------------------
    ! Normalise the angular viabilities
    !---------------------------------------------------------------------------
    if(num_3body.gt.0)then
       viability_3body = viability_3body ** ( &
            1._real32 / real(num_3body,real32) &
       )
    end if
    if(num_4body.gt.0)then
       viability_4body = viability_4body ** ( &
            1._real32 / real(num_4body,real32) &
       )
    end if


    !---------------------------------------------------------------------------
    ! Combine the 2-, 3- and 4-body viabilities to get the overall viability
    !---------------------------------------------------------------------------
    output = viability_2body * viability_3body * viability_4body

  end function evaluate_point
!###############################################################################


!###############################################################################
  pure function evaluate_2body_contributions( distribs_container, &
       bondlength, pair_index &
  ) result(output)
    !! Return the contribution to the viability from 2-body interactions
    implicit none

    ! Arguments
    type(distribs_container_type), intent(in) :: distribs_container
    !! Distribution function (gvector) container.
    real(real32), intent(in) :: bondlength
    !! Bond length.
    integer, intent(in) :: pair_index
    !! Index of the element pair.
    real(real32) :: output
    !! Contribution to the viability.


    output = distribs_container%gdf%df_2body( &
         distribs_container%get_bin(bondlength, dim = 1), &
         pair_index &
    )

  end function evaluate_2body_contributions
!###############################################################################


!###############################################################################
  pure function evaluate_3body_contributions( distribs_container, &
       position_1, position_2, basis, species, current_idx &
  ) result(output)
    !! Return the contribution to the viability from 3-body interactions
    implicit none

    ! Arguments
    type(distribs_container_type), intent(in) :: distribs_container
    !! Distribution function (gvector) container.
    real(real32), dimension(3), intent(in) :: position_1
    !! Positions of the atom.
    real(real32), dimension(4), intent(in) :: position_2
    !! Positions of the second atom and its cutoff weighting.
    type(extended_basis_type), intent(in) :: basis
    !! Basis of the system.
    integer, intent(in) :: species
    !! Index of the query element.
    integer, dimension(2), intent(in) :: current_idx
    !! Index of the 1st-atom query element.
    real(real32) :: output
    !! Contribution to the viability.

    ! Local variables
    integer :: js, ja
    !! Loop indices.
    integer :: bin
    !! Bin for the distribution function.
    integer :: num_3body_local
    !! Number of 3-body interactions local to the current atom pair.
    real(real32) :: power
    !! Power for the contribution to the viability.


    num_3body_local = sum(basis%spec(current_idx(1):)%num) - current_idx(2)
    if(num_3body_local.eq.0)then
       output = -1000._real32
       return
    end if
    output = 1._real32
    power = 1._real32 / real( num_3body_local, real32 )
    species_loop: do js = current_idx(1), basis%nspec, 1
       atom_loop: do ja = 1, basis%spec(js)%num
          if(js.eq.current_idx(1) .and. ja.le.current_idx(2))cycle
          associate( position_store => [ basis%spec(js)%atom(ja,1:4) ] )
             bin = distribs_container%get_bin( &
                  get_angle( &
                       position_2(1:3), &
                       position_1, &
                       position_store(1:3) &
                  ), dim = 2 &
             )
             output = output * &
                  ( 1._real32 + ( position_2(4) * position_store(4) ) * ( &
                       distribs_container%gdf%df_3body( &
                            bin, &
                            distribs_container%element_map(species) &
                       ) - distribs_container%viability_3body_default &
                  ) ) ** power
          end associate
       end do atom_loop
    end do species_loop

  end function evaluate_3body_contributions
!###############################################################################


!###############################################################################
  pure function evaluate_4body_contributions( distribs_container, &
       position_1, position_2, position_3, basis, species ) result(output)
    !! Return the contribution to the viability from 4-body interactions
    implicit none

    ! Arguments
    type(distribs_container_type), intent(in) :: distribs_container
    !! Distribution function (gvector) container.
    real(real32), dimension(3), intent(in) :: position_1
    !! Positions of the atoms.
    real(real32), dimension(4), intent(in) :: position_2, position_3
    type(extended_basis_type), intent(in) :: basis
    !! Basis of the system.
    integer, intent(in) :: species
    !! Index of the query element.
    real(real32) :: output
    !! Contribution to the viability.

    ! Local variables
    integer :: ks, ka
    !! Loop indices.
    integer :: bin
    !! Bin for the distribution function.
    integer :: num_4body_local
    !! Number of 4-body interactions local to the current atom triplet.
    real(real32) :: power
    !! Power for the contribution to the viability.


    num_4body_local = sum(basis%image_spec(:)%num)
    if(num_4body_local.eq.0)then
       output = -1000._real32
       return
    end if
    output = 1._real32
    power = 1._real32 / real( num_4body_local, real32 )
    species_loop: do ks = 1, basis%nspec, 1
       atom_loop: do ka = 1, basis%image_spec(ks)%num
          associate( position_store => [ basis%image_spec(ks)%atom(ka,1:4) ] )
             bin = distribs_container%get_bin( &
                  get_improper_dihedral_angle( &
                       position_1, &
                       position_2(1:3), &
                       position_3(1:3), &
                       position_store(1:3) &
                  ), dim = 3 &
             )
             output = output * ( &
                  1._real32 + ( &
                       position_2(4) * position_3(4) * position_store(4) &
                  ) * ( &
                       distribs_container%gdf%df_4body( &
                            bin, &
                            distribs_container%element_map(species) &
                       ) - distribs_container%viability_4body_default &
                  ) ) ** power
          end associate
       end do atom_loop
    end do species_loop

  end function evaluate_4body_contributions
!###############################################################################

end module raffle__evaluator
