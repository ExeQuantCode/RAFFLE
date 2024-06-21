module buildmap
  use constants, only: real12
  use misc_linalg, only: get_distance, get_angle, get_dihedral_angle
  use rw_geom, only: bas_type
  use edit_geom, only: &
       get_min_dist_between_two_atoms, &
       get_min_dist_between_point_and_atom
  use evolver, only: gvector_container_type
  implicit none

contains

!!!#############################################################################
!!! builds a map of the system and returns the value of the map at a given point
!!!#############################################################################
!!! output = suitability of tested point
  pure function buildmap_POINT(gvector_container, &
       position, lattice, basis, atom_ignore_list, &
       radius_list, uptol, lowtol) &
       result(output)
    implicit none
    type(gvector_container_type), intent(in) :: gvector_container
    real(real12), intent(in) :: uptol, lowtol
    type(bas_type), intent(in) :: basis
    real(real12), dimension(3), intent(in) :: position
    integer, dimension(:,:), intent(in) :: atom_ignore_list
    real(real12), dimension(3,3), intent(in) :: lattice
    real(real12), dimension(:), intent(in) :: radius_list
    real(real12) :: output
  
    integer :: i
    integer :: is, ia, js, ja, ks, ka, ls
    integer :: bin
    real(real12) :: contribution, repeat_power, bondlength
    real(real12) :: viability_2body !! 2-body is addition
    real(real12) :: viability_3body !! 3-body is multiplication
    real(real12) :: viability_4body !! 4-body is multiplication
    real(real12), dimension(3) :: &
         position_storage1, position_storage2, position_storage3
    integer, dimension(:,:), allocatable :: pair_index
  

    output = 0._real12
    repeat_power = 1._real12
    viability_2body = 0._real12
    viability_3body = 1._real12
    viability_4body = 1._real12
    

    !! get list of element pair indices
    ls = atom_ignore_list(1,1)
    allocate(pair_index(basis%nspec, basis%nspec), source = 0)
    do is = 1, basis%nspec
       do js = 1, basis%nspec
          pair_index(is, js) = gvector_container%get_pair_index( &
               basis%spec(is)%name, basis%spec(js)%name )
       end do
    end do


    species_loop1: do is=1, basis%nspec
      !! loops over all atoms currently in the system
      !! 2-body map
      !! checks bondlength between the current atom and all other atoms
      atom_loop1: do ia = 1, basis%spec(is)%num
         do i = 1, size(atom_ignore_list,dim=1), 1
            if(all(atom_ignore_list(i,:).eq.[is,ia])) cycle atom_loop1
         end do
         !!! NEED TO HAVE A LOOP FOR REPEATING CELLS
         position_storage1 = basis%spec(is)%atom(ia,:)
         contribution=0
         !!! ONLY NEEDS TO CHECK FOR THE SMALLEST BONDLENGTH BETWEEN A ...
         !!! ... PERIODIC ATOM AND THE CURRENT ATOM
         !!! if not looping over periodic images explicitly, doesn't need ...
         !!! ... to check if the bondlength is larger than the upper tolerance
         !!! But cycling still needed as it fits that criteria if it is ...
         !!! ... above upper tolerance, and doesn't fall within 3- and 4-body ...
         !!! ... check requirements
    
         bondlength = get_min_dist_between_point_and_atom(lattice, basis, &
              position, [is, ia])

         !! check if the bondlength is within the tolerance for bonds ...
         !! ... between its own element and the element of the current atom
         if(bondlength .lt. radius_list(pair_index(ls,is))*lowtol)then
            return
         else if(bondlength .gt. radius_list(pair_index(ls,is))*uptol)then
            cycle atom_loop1
         end if

         bin = gvector_container%get_bin(bondlength, dim = 1)
         if(bin.eq.0) cycle
         contribution = gvector_container%total%df_2body(bin, pair_index(ls,is))
   
         viability_2body = viability_2body + contribution

         !! loops over all atoms currently in the system
         !! 3-body map
         !! checks bondangle between the current atom and all other atoms
         !! i.e. nested loop here
         !!! NEEDS TO BE SPECIES AND ATOM
         !!! SHOULD BE ITS OWN PROCEDURE
         species_loop2: do js = is, basis%nspec, 1
           atom_loop2: do ja = 1, basis%spec(js)%num
              if(js.eq.is .and. ja.lt.ia) cycle
              do i = 2, size(atom_ignore_list,dim=1)
                 if(all(atom_ignore_list(i,:).eq.[js,ja])) cycle atom_loop1
              end do
              position_storage2 = basis%spec(js)%atom(ja,:)
              if(get_distance(position,position_storage2).lt.&
                   radius_list(pair_index(ls,js))*lowtol) return
              if(get_distance(position,position_storage2).lt.&
                   radius_list(pair_index(ls,js))*uptol) then
                 bin = gvector_container%get_bin( &
                      get_angle( position_storage1, &
                                 position, &
                                 position_storage2 ), &
                      dim = 2 )
                 if(bin.eq.0) cycle
                 contribution = gvector_container%total%df_3body(bin,is)
                 viability_3body = ( viability_3body * &
                       contribution ** (1._real12/repeat_power))
              end if
              !!! CHECKS WHETHER THE BONDLENGTH BETWEEN THE CURRENT ATOM AND THE ...
              !!! ... THIRD ATOM IS WITHIN THE TOLERANCE
              !!! I have removed the second check as this, again, is just checking ...
              !!! ... the effect of a periodic image
              if((get_distance(position_storage1,position_storage2).ge.&
                   radius_list(pair_index(ls,js))*uptol)) cycle
                 
              !! loops over all atoms currently in the system
              !! 4-body map
              !! checks dihedral angle between the current atom and all other atoms
              !! i.e. nested loop here
              !!! NEEDS TO BE SPECIES AND ATOM
              !!! SHOULD BE ITS OWN PROCEDURE
              species_loop3: do ks = 1, basis%nspec, 1
                 atom_loop3: do ka = 1, basis%spec(ks)%num
                    do i = 2, size(atom_ignore_list,dim=1)
                       if(all(atom_ignore_list(i,:).eq.[ks,ka])) cycle atom_loop1
                    end do
                    position_storage3 = basis%spec(js)%atom(ja,:)
                    if(get_distance(position,position_storage3).lt.&
                         radius_list(pair_index(ls,ks))*lowtol) return
                    if(get_distance(position_storage1,position_storage3).lt.&
                         radius_list(pair_index(is,ks))*uptol) then
                       bin = gvector_container%get_bin( &
                                get_dihedral_angle( &
                                           position, &
                                           position_storage1, &
                                           position_storage2, &
                                           position_storage3), &
                                dim = 3 )
                       if(bin.eq.0) cycle
                       contribution = gvector_container%total%df_4body(bin,is)
                       if(abs(contribution).lt.1.E-6) return
                       !Here have taken a large root of value return, to ...
                       !...account for sumamtion of atoms in 3D.
                       !Will need to think further on this.
                       viability_4body = ( viability_4body * &
                            contribution ** (1._real12/repeat_power))
                    end if
                 end do atom_loop3
              end do species_loop3
           end do atom_loop2
         end do species_loop2
      end do atom_loop1
   end do species_loop1

   if(abs(viability_2body).lt.1.E-6) viability_2body = 1._real12
   output = viability_2body * viability_4body * viability_3body
    
  end function buildmap_POINT
!!!#############################################################################

end module buildmap