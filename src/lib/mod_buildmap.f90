module buildmap
  use constants, only: real12
  use misc_linalg, only: get_distance, get_angle, get_dihedral_angle
  use rw_geom, only: bas_type
  use edit_geom, only: get_min_dist_between_two_atoms
  use contributions, only: get_2body_contribution, get_3body_contribution, &
       get_4body_contribution
  implicit none

contains

!!!#############################################################################
!!! builds a map of the system and returns the value of the map at a given point
!!!#############################################################################
!!! output = suitability of tested point
  function buildmap_POINT(position, lattice, basis, atom_ignore_list, &
       radius_arr, uptol, lowtol) &
       result(output)
    implicit none
    real(real12), intent(in) :: uptol, lowtol
    type(bas_type), intent(in) :: basis
    real(real12), dimension(3), intent(inout) :: position
    integer, dimension(:,:), intent(in) :: atom_ignore_list
    real(real12), dimension(3,3), intent(in) :: lattice
    real(real12), dimension(:,:,:), intent(in) :: radius_arr
    real(real12) :: output
  
    integer :: i
    integer :: is, ia, js, ja, ks, ka
    real(real12) :: contribution, repeat_power, bondlength
    real(real12) :: viability_2body = 0._real12 !! 2-body is addition
    real(real12) :: viability_3body = 1._real12 !! 3-body is multiplication
    real(real12) :: viability_4body = 1._real12 !! 4-body is multiplication
    real(real12), dimension(3,3) :: position_storage
  

    output = 0._real12
    repeat_power = 1._real12
    
    species_loop1: do is=1, basis%nspec
      !! loops over all atoms currently in the system
      !! 2-body map
      !! checks bondlength between the current atom and all other atoms
      atom_loop1: do ia = 1, basis%spec(is)%num
         do i = 2, size(atom_ignore_list,dim=1)
            if(all(atom_ignore_list(i,:).eq.[is,ia])) cycle atom_loop1
         end do
         !!! NEED TO HAVE A LOOP FOR REPEATING CELLS
         position_storage(1,:) = basis%spec(is)%atom(ia,:)
         contribution=0
         !!! ONLY NEEDS TO CHECK FOR THE SMALLEST BONDLENGTH BETWEEN A ...
         !!! ... PERIODIC ATOM AND THE CURRENT ATOM
         !!! if not looping over periodic images explicitly, doesn't need ...
         !!! ... to check if the bondlength is larger than the upper tolerance
         !!! But cycling still needed as it fits that criteria if it is ...
         !!! ... above upper tolerance, and doesn't fall within 3- and 4-body ...
         !!! ... check requirements
    

         bondlength = get_min_dist_between_two_atoms(lattice, basis, &
              atom_ignore_list(1,:), [is,ia])
         !! check if the bondlength is within the tolerance for bonds ...
         !! ... between its own element and the element of the current atom
         if(bondlength .lt. radius_arr(1,atom_ignore_list(1,1),is)*lowtol)then
            viability_2body = 0._real12
            viability_3body = 0._real12
            viability_4body = 0._real12
            exit species_loop1
         else if(bondlength .gt. radius_arr(1,atom_ignore_list(1,1),is)*uptol)then
    !!!! This has been left out to ease on computation. Could be reimplemented but increases cost dramatically to consider ALL atoms
            !call get_contribution (trim(adjustl(&
            !     &atomlist(structures,atom_number_previous+1)%name)),&
            !     &trim(adjustl(alistrep(structures,atom_number_previous+1)%name)),get_distance(position,&
            !     &position_storage(1,:)),contribution)
            cycle
         end if
         contribution = get_2body_contribution( &
               trim(adjustl(basis%spec(atom_ignore_list(1,1))%name)),&
               trim(adjustl(basis%spec(is)%name)),&
               bondlength)
   
         viability_2body = viability_2body + contribution

         !! loops over all atoms currently in the system
         !! 3-body map
         !! checks bondangle between the current atom and all other atoms
         !! i.e. nested loop here
         !!! NEEDS TO BE SPECIES AND ATOM
         !!! SHOULD BE ITS OWN PROCEDURE
         species_loop2: do js = 1, basis%nspec
            atom_loop2: do ja = 1, basis%spec(js)%num
               !!! NAH, DON'T NEED TO LOOP ITSELF ANYMORE AS CHECKS NEED TO ...
               !!! ... BE MADE WITH PERIODIC IMAGES
               !if(all([is,l].eq.[js,ja])) cycle
               do i = 2, size(atom_ignore_list,dim=1)
                  if(all(atom_ignore_list(i,:).eq.[js,ja])) cycle atom_loop1
               end do
               position_storage(2,:) = basis%spec(js)%atom(ja,:)

               if(get_distance(position,position_storage(2,:)).lt.&
                    radius_arr(1,atom_ignore_list(1,1),js)*lowtol) then
                  viability_2body = 0._real12
                  viability_3body = 0._real12
                  viability_4body = 0._real12
                  exit species_loop1
               end if
               !!! ARE WE NOT DOUBLE COUNTING HERE!?!
               !!! by looking at the angle between p1, p2, and p3
               if(get_distance(position_storage(1,:),position_storage(2,:)).lt.&
                    radius_arr(1,is,js)*uptol) then
   
                  contribution = get_3body_contribution( &
                       trim(adjustl(basis%spec(atom_ignore_list(1,1))%name)),&
                       get_angle( &
                       position,position_storage(1,:),position_storage(2,:)))

                  viability_3body = ( viability_3body * &
                       contribution ** (1._real12/(repeat_power)))
               else if(get_distance(position,position_storage(2,:)).lt.&
                    radius_arr(1,atom_ignore_list(1,1),js)*uptol) then 
   
                  contribution = get_3body_contribution( &
                       trim(adjustl(basis%spec(atom_ignore_list(1,1))%name)),&
                       get_angle(&
                       position_storage(1,:),position,position_storage(2,:)))
      
                  viability_3body = ( viability_3body * &
                        contribution ** (1._real12/repeat_power))
               end if
   
               !!! CHECKS WHETHER THE BONDLENGTH BETWEEN THE CURRENT ATOM AND THE ...
               !!! ... THIRD ATOM IS WITHIN THE TOLERANCE
               !!! I have removed the second check as this, again, is just checking ...
               !!! ... the effect of a periodic image
               if((get_distance(position_storage(1,:),position_storage(2,:)).ge.&
                    radius_arr(1,is,js)*uptol)) cycle
                  
               !! loops over all atoms currently in the system
               !! 4-body map
               !! checks dihedral angle between the current atom and all other atoms
               !! i.e. nested loop here
               !!! NEEDS TO BE SPECIES AND ATOM
               !!! SHOULD BE ITS OWN PROCEDURE
               species_loop3: do ks = 1, basis%nspec
                  atom_loop3: do ka = 1, basis%spec(ks)%num
                     !!! NAH, DON'T NEED TO LOOP ITSELF ANYMORE AS CHECKS NEED TO ...
                     !!! ... BE MADE WITH PERIODIC IMAGES
                     !if(all([is,l].eq.[ks,ka]).or.all([js,ja].eq.[ks,ka])) cycle
                     do i = 2, size(atom_ignore_list,dim=1)
                        if(all(atom_ignore_list(i,:).eq.[ks,ka])) cycle atom_loop1
                     end do
                     position_storage(3,:) = basis%spec(js)%atom(ja,:)
                     if(get_distance(position,position_storage(3,:)).lt.&
                          radius_arr(1,atom_ignore_list(1,1),ks)*lowtol) then
                        viability_2body = 0._real12
                        viability_3body = 0._real12
                        viability_4body = 0._real12
                        exit species_loop1
                     end if
                     if(get_distance(position_storage(1,:),position_storage(3,:)).lt.&
                          radius_arr(1,is,ks)*uptol) then
                        contribution = get_4body_contribution( &
                             trim(adjustl(basis%spec(atom_ignore_list(1,1))%name)),&
                             get_dihedral_angle(&
                             position,position_storage(1,:),&
                             position_storage(2,:),position_storage(3,:)))
                        if(contribution.eq.0) then
                           viability_2body = 0._real12
                           viability_3body = 0._real12
                           viability_4body = 0._real12
                           exit species_loop1
                        end if
      
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

   output = viability_2body * viability_4body * viability_3body
    
  end function buildmap_POINT
!!!#############################################################################

end module buildmap