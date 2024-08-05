module evaluator
  !! Module to build viability map of a structure
  !!
  !! This module handles the viability map for a structure, which is a map of
  !! the system with each point in the map representing the suitability of
  !! that point for a new atom. The map is built by checking the bond lengths,
  !! bond angles and dihedral angles between the test point and all atoms.
  use constants, only: real12
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
  pure function evaluate_point(gvector_container, &
       position, basis, atom_ignore_list, &
       radius_list, uptol, lowtol) &
       result(output)
    !! Build a map of basis and returns the value of the map at a given point
    implicit none

    ! Arguments
    type(gvector_container_type), intent(in) :: gvector_container
    !! Distribution function (gvector) container.
    real(real12), intent(in) :: uptol, lowtol
    !! Upper and lower tolerance for bond lengths and angles.
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
    integer :: i, is, js, ia, ls
    !! Loop counters.
    integer :: bin
    !! Bin for the distribution function.
    real(real12) :: contribution
    !! Contribution to the viability map
    real(real12) :: viability_2body
    !! Viability of the test point for 2-body interactions.
    !! 2-body viabilities are summed.
    real(real12) :: viability_angles
    !! Viability of the test point for 3- and 4-body interactions.
    !! Angle viabilities are multiplied.
    real(real12), dimension(3) :: position_store
    !! Storage for atom positions.
    integer, dimension(:,:), allocatable :: pair_index
    !! Index of element pairs.
 
    
    ! Initialisation
    output = 0._real12
    viability_2body = 0._real12
    viability_angles = 1._real12
    

    !---------------------------------------------------------------------------
    ! get list of element pair indices
    !---------------------------------------------------------------------------
    ls = atom_ignore_list(1,1)
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
    species_loop: do is=1, basis%nspec
       ! 2-body map
       ! check bondlength between test point and all other atoms
       !-------------------------------------------------------------------------
       atom_loop: do ia = 1, basis%spec(is)%num
          do i = 1, size(atom_ignore_list,dim=1), 1
             if(all(atom_ignore_list(i,:).eq.[is,ia])) cycle atom_loop
          end do
          position_store = basis%spec(is)%atom(ia,:)
          contribution = get_2body_contribution( gvector_container, &
               position, position_store, basis%lat, &
               radius_list(pair_index(ls,is))*lowtol, &
               gvector_container%cutoff_max(1), &
               pair_index(ls,is) &
          )
          if(contribution .lt. -100._real12)then
             return
          elseif(contribution.lt.-50._real12)then
             cycle atom_loop
          end if
          viability_2body = viability_2body + contribution

          ! 3-body map
          ! check bondangle between test point and all other atoms
          !----------------------------------------------------------------------
          viability_angles = viability_angles * &
               evaluate_3body_contributions( gvector_container, &
                    position, position_store, basis, atom_ignore_list, &
                    radius_list, uptol, lowtol, pair_index, ls, [is, ia] &
               )
          if(viability_angles .lt. 1.E-6)then
             return
          end if
       end do atom_loop

       image_loop: do ia = 1, basis%image_spec(is)%num, 1
          position_store = basis%image_spec(is)%atom(ia,:)
          contribution = get_2body_contribution( gvector_container, &
               position, position_store, basis%lat, &
               radius_list(pair_index(ls,is))*lowtol, &
               gvector_container%cutoff_max(1), &
               pair_index(ls,is) &
          )
          if(contribution .lt. -100._real12)then
             return
          elseif(contribution.lt.-50._real12)then
             cycle image_loop
          end if
          viability_2body = viability_2body + contribution

          ! 3-body map
          ! check bondangle between test point and all other atoms
          !----------------------------------------------------------------------
          viability_angles = viability_angles * &
               evaluate_3body_contributions( gvector_container, &
                    position, position_store, basis, atom_ignore_list, &
                    radius_list, uptol, lowtol, pair_index, ls, &
                    [is, basis%spec(is)%num + ia] &
               )
          if(viability_angles .lt. 1.E-6)then
             return
          end if
       end do image_loop
    end do species_loop
    !!! CHECK FOR NaN VALUES FOR VIABILITY_ANGLES

    if(abs(viability_2body).lt.1.E-6) viability_2body = 1._real12
    output = viability_2body * viability_angles

    deallocate(pair_index)
    
  end function evaluate_point
!###############################################################################


!###############################################################################
  pure function get_2body_contribution( gvector_container, &
       position_1, position_2, lattice, lower_limit, upper_limit, pair_index ) &
       result(output)
    !! Return the contribution to the viability map for a 2-body interaction
    implicit none

    ! Arguments
    type(gvector_container_type), intent(in) :: gvector_container
    !! Distribution function (gvector) container.
    real(real12), dimension(3), intent(in) :: position_1, position_2
    !! Positions of the atoms.
    real(real12), dimension(3,3), intent(in) :: lattice
    !! Lattice vectors.
    real(real12), intent(in) :: lower_limit, upper_limit
    !! Lower and upper limits for the bond length.
    integer, intent(in) :: pair_index
    !! Index of the pair of elements.
    real(real12) :: output
    !! Contribution from the 2-body interaction.

    ! Local variables
    integer :: bin
    !! Bin for the distribution function.
    real(real12) :: bondlength
    !! Bond length between the test point and the second atom.


    output = 0._real12
    
    bondlength = modu( matmul(position_1 - position_2, lattice) )

    !! check if the bondlength is within the tolerance for bonds ...
    !! ... between its own element and the element of the current atom
    if(bondlength .lt. lower_limit)then
       output = -120._real12
       return
    else if(bondlength.gt.gvector_container%cutoff_max(1))then
       output = -60._real12
       return
    end if

    bin = gvector_container%get_bin(bondlength, dim = 1)
    if(bin.eq.0)then
       output = 0._real12
       return
    end if
    output = gvector_container%total%df_2body(bin, pair_index)

  end function get_2body_contribution
!###############################################################################


!###############################################################################
  pure function evaluate_3body_contributions( gvector_container, &
       position_1, position_2, basis, &
       atom_ignore_list, radius_list, uptol, lowtol, &
       pair_index, ls, atom_index ) result(output)
    !! Return the contribution to the viability map from 3-body interactions
    implicit none

    ! Arguments
    type(gvector_container_type), intent(in) :: gvector_container
    !! Distribution function (gvector) container.
    real(real12), dimension(3), intent(in) :: position_1, position_2
    !! Positions of the atoms.
    type(extended_basis_type), intent(in) :: basis
    !! Basis of the system.
    integer, dimension(:,:), intent(in) :: atom_ignore_list
    !! List of atoms to ignore (i.e. indices of atoms not yet placed).
    real(real12), dimension(:), intent(in) :: radius_list
    !! List of radii for each pair of elements.
    real(real12), intent(in) :: uptol, lowtol
    !! Upper and lower tolerance for bond lengths and angles.
    integer, dimension(:,:), intent(in) :: pair_index
    !! Index of the pair of elements.
    integer, intent(in) :: ls
    !! Index of the query element.
    integer, dimension(2), intent(in) :: atom_index
    !! Index of the 1st atom.
    real(real12) :: output
    !! Contribution to the viability map.

    ! Local variables
    integer :: i, js, ja
    !! Bin for the distribution function.
    real(real12) :: contribution
    !! Contribution to the viability map.
    real(real12) :: repeat_power
    !! Repeat power for 3-body interactions.
    real(real12), dimension(3) :: position_store
    !! Storage for atom positions.
    real(real12) :: viability_4body
    !! Viability of the test point for 4-body interactions.


    repeat_power = 2._real12
    output = 1._real12
    viability_4body = 1._real12
    species_loop: do js = atom_index(1), basis%nspec, 1
      atom_loop: do ja = 1, basis%spec(js)%num
         if(js.eq.atom_index(1) .and. ja.le.atom_index(2)) cycle atom_loop
         do i = 1, size(atom_ignore_list,dim=1)
            if(all(atom_ignore_list(i,:).eq.[js,ja])) cycle atom_loop
         end do
         position_store = basis%image_spec(js)%atom(ja,:)
         contribution = get_3body_contribution( gvector_container, &
              position_1, position_2, position_store, basis%lat, &
              radius_list(pair_index(ls,js))*lowtol, &
              radius_list(pair_index(ls,js))*uptol, &
              ls &
         )
         if (contribution .lt. -100._real12) then
            output = 0._real12
            return
         elseif (contribution .lt. -50._real12) then
            cycle atom_loop
         end if
         output = output * contribution
            
         ! ! 4-body map
         ! ! check improperdihedral angle between test point and all other atoms
         ! !----------------------------------------------------------------------
         ! contribution = evaluate_4body_contributions( gvector_container, &
         !      position_1, position_2, position_store, &
         !      basis, atom_ignore_list, radius_list, uptol, lowtol, &
         !      pair_index, ls, [js, ja] )
         ! if (contribution .lt. 1.E-6) then
         !    output = 0._real12
         !    return
         ! end if
         ! viability_4body = viability_4body * contribution
      end do atom_loop

      image_loop: do ja = 1, basis%image_spec(js)%num, 1
         if( js.eq.atom_index(1) .and. &
              basis%spec(js)%num + ja.le.atom_index(2)) cycle
         position_store = basis%image_spec(js)%atom(ja,:)
         contribution = get_3body_contribution( gvector_container, &
              position_1, position_2, position_store, basis%lat, &
              radius_list(pair_index(ls,js))*lowtol, &
              radius_list(pair_index(ls,js))*uptol, &
              ls &
         )
         if (contribution .lt. -100._real12) then
            output = 0._real12
            return
         elseif (contribution .lt. -50._real12) then
            cycle image_loop
         end if
         output = output * contribution

         ! ! 4-body map
         ! ! check improperdihedral angle between test point and all other atoms
         ! !----------------------------------------------------------------------
         ! contribution = evaluate_4body_contributions( gvector_container, &
         !      position_1, position_2, position_store, &
         !      basis, atom_ignore_list, radius_list, uptol, lowtol, &
         !      pair_index, ls, [js, basis%spec(js)%num + ja] )
         ! if (contribution .lt. 1.E-6) then
         !    output = 0._real12
         !    return
         ! end if
         ! viability_4body = viability_4body * contribution
      end do image_loop
    end do species_loop
    output = output ** (1._real12/repeat_power)
    output = output * viability_4body

  end function evaluate_3body_contributions
!###############################################################################


!###############################################################################
  pure function get_3body_contribution( gvector_container, &
       position_1, position_2, position_3, &
       lattice, lower_limit, upper_limit, ls ) &
       result(output)
    !! Return the contribution to the viability map for a 3-body interaction
    implicit none

    ! Arguments
    type(gvector_container_type), intent(in) :: gvector_container
    !! Distribution function (gvector) container.
    real(real12), dimension(3), intent(in) :: position_1, position_2, position_3
    !! Positions of the atoms.
    real(real12), dimension(3,3), intent(in) :: lattice
    !! Lattice vectors.
    real(real12), intent(in) :: lower_limit, upper_limit
    !! Lower and upper limits for the bond angle.
    integer, intent(in) :: ls
    !! Index of the query element.
    real(real12) :: output
    !! Contribution from the 3-body interaction.

    ! Local variables
    integer :: bin
    !! Bin for the distribution function.
    real(real12) :: bondlength
    !! Bond length between the test point and the third atom.


    output = 1._real12
    bondlength = modu( matmul(position_1 - position_3, lattice) )
    if(bondlength.lt.lower_limit)then
       output = -120._real12
       return
    else if(bondlength.lt.upper_limit) then
       bin = gvector_container%get_bin( &
            get_angle( position_2, &
                       position_1, &
                       position_3 ), &
            dim = 2 )
       if(bin.eq.0) return
       output = gvector_container%total%df_3body(bin,ls)
    else
       output = -60._real12
    end if

  end function get_3body_contribution
!###############################################################################


!###############################################################################
  pure function evaluate_4body_contributions( gvector_container, &
       position_1, position_2, position_3, basis, &
       atom_ignore_list, radius_list, uptol, lowtol, &
       pair_index, ls, atom_index ) result(output)
    !! Return the contribution to the viability map from 4-body interactions
    implicit none

    ! Arguments
    type(gvector_container_type), intent(in) :: gvector_container
    !! Distribution function (gvector) container.
    real(real12), dimension(3), intent(in) :: position_1, position_2, position_3
    !! Positions of the atoms.
    type(extended_basis_type), intent(in) :: basis
    !! Basis of the system.
    integer, dimension(:,:), intent(in) :: atom_ignore_list
    !! List of atoms to ignore (i.e. indices of atoms not yet placed).
    real(real12), dimension(:), intent(in) :: radius_list
    !! List of radii for each pair of elements.
    real(real12), intent(in) :: uptol, lowtol
    !! Upper and lower tolerance for bond lengths and angles.
    integer, dimension(:,:), intent(in) :: pair_index
    !! Index of the pair of elements.
    integer, dimension(2), intent(in) :: atom_index
    !! Index of the 1st atom.
    integer, intent(in) :: ls
    !! Index of the query element.
    real(real12) :: output
    !! Contribution to the viability map.

    ! Local variables
    integer :: i, ks, ka
    !! Loop indices.
    integer :: bin
    !! Bin for the distribution function.
    real(real12) :: contribution
    !! Contribution to the viability map.
    real(real12) :: repeat_power
    !! Repeat power for 4-body interactions.
    real(real12), dimension(3) :: position_store
    !! Storage for atom positions.


    repeat_power = 4._real12
    output = 1._real12
    species_loop: do ks = atom_index(1), basis%nspec, 1
       atom_loop: do ka = 1, basis%spec(ks)%num
          if(ks.eq.atom_index(1) .and. ka.le.atom_index(2)) cycle atom_loop
          do i = 1, size(atom_ignore_list,dim=1)
             if(all(atom_ignore_list(i,:).eq.[ks,ka])) cycle atom_loop
          end do
          position_store = basis%spec(ks)%atom(ka,:)
          contribution = get_4body_contribution( gvector_container, &
               position_1, position_2, position_3, position_store, &
               basis%lat, radius_list(pair_index(ls,ks))*lowtol, &
               radius_list(pair_index(ls,ks))*uptol, ls )
          if(contribution .lt. -100._real12) then
             output = 0._real12
             return
          elseif(contribution .lt. -50._real12) then
             cycle atom_loop
          end if
          output = output * contribution
       end do atom_loop

       image_loop: do ka = 1, basis%image_spec(ks)%num, 1
          if(ks.eq.atom_index(1) .and. &
             basis%spec(ks)%num + ka.le.atom_index(2)) cycle
          position_store = basis%image_spec(ks)%atom(ka,:)
          contribution = get_4body_contribution( gvector_container, &
               position_1, position_2, position_3, position_store, &
               basis%lat, radius_list(pair_index(ls,ks))*lowtol, &
               radius_list(pair_index(ls,ks))*uptol, ls )
          if(contribution .lt. -100._real12) then
             output = 0._real12
             return
          elseif(contribution .lt. -50._real12) then
             cycle image_loop
          end if
          output = output * contribution
       end do image_loop
    end do species_loop
    output = output ** (1._real12/repeat_power)

  end function evaluate_4body_contributions
!###############################################################################


!###############################################################################
  pure function get_4body_contribution( gvector_container, &
       position_1, position_2, position_3, position_4, &
       lattice, lower_limit, upper_limit, ls ) &
       result(output)
    !! Return the contribution to the viability map for a 4-body interaction
    implicit none

    ! Arguments
    type(gvector_container_type), intent(in) :: gvector_container
    !! Distribution function (gvector) container.
    real(real12), dimension(3), intent(in) :: position_1, position_2, position_3, position_4
    !! Positions of the atoms.
    real(real12), dimension(3,3), intent(in) :: lattice
    !! Lattice vectors.
    real(real12), intent(in) :: lower_limit, upper_limit
    !! Lower and upper limits for the dihedral angle.
    integer, intent(in) :: ls
    !! Index of the query element.
    real(real12) :: output
    !! Contribution from the 4-body interaction.

    ! Local variables
    integer :: bin
    !! Bin for the distribution function.
    real(real12) :: bondlength
    !! Bond length between the test point and the fourth atom.


    output = 1._real12
    bondlength = modu( matmul(position_1 - position_4, lattice) )
    if(bondlength.lt.lower_limit) then
      output = -120._real12
      return
    else if(bondlength.lt.upper_limit) then
      bin = gvector_container%get_bin( &
               get_improper_dihedral_angle( &
                          position_1, &
                          position_2, &
                          position_3, &
                          position_4), &
               dim = 3 )
      !!! HANDLE IF IT IS GREATER THAN 90 DEGREES!!!
      !!! Should map back into 0-90 degrees !!!
      if(bin.eq.0) return
      output = gvector_container%total%df_4body(bin,ls)
    else
      output = -60._real12
      return
    end if

   end function get_4body_contribution
!###############################################################################

end module evaluator