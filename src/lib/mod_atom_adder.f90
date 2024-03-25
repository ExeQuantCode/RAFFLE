module add_atom
  use constants, only: real12
  use misc_linalg, only: modu
  use rw_geom, only: bas_type
  use edit_geom, only: get_min_dist
  use buildmap, only: buildmap_POINT
  use evolver, only: gvector_container_type
  implicit none

contains


!!!#############################################################################
!!! add atom to unit cell using the scan method
!!!#############################################################################
  subroutine add_atom_scan (bin_size, gvector_container, &
       lattice, basis, atom_ignore_list, &
       radius_arr, placed)
    implicit none
    type(gvector_container_type), intent(in) :: gvector_container
    type(bas_type), intent(inout) :: basis
    logical, intent(out) :: placed
    integer, dimension(3), intent(in) :: bin_size
    integer, dimension(:,:), intent(in) :: atom_ignore_list
    real(real12), dimension(3,3) :: lattice
    real(real12), dimension(:,:,:) :: radius_arr

    integer :: el_correct_read,i, j, k,n,l
    integer, dimension(3) :: best_gridpoint
    real(real12), dimension(3) :: tmpvector
    real(real12), dimension(:,:,:), allocatable :: suitability_grid


    placed = .false. 
    allocate(suitability_grid(0:bin_size(1)-1,0:bin_size(2)-1,0:bin_size(3)-1))
    ! run buildmap_point for a set of points in the unit cell
    ! set up a grid and run with it
    ! easiest way to do that would be to find all atoms, find their radii of ...
    ! ... influence, and make a set of points for all that lie outside of them
    ! have it species-dependent   
    do concurrent( i = 0:bin_size(1)-1, &
                   j = 0:bin_size(2)-1, &
                   k = 0:bin_size(3)-1 )
       suitability_grid(i,j,k) = &
            buildmap_POINT( gvector_container, &
            [i, j, k] / real(bin_size,real12), lattice, basis, &
            atom_ignore_list, radius_arr, &
            1.1_real12, 0.95_real12)
    end do
    if(abs(maxval(suitability_grid)).lt.1.E-6) return

    placed = .true.
    best_gridpoint = maxloc(suitability_grid) + lbound(suitability_grid) - 1

    basis%spec(atom_ignore_list(1,1))%atom(atom_ignore_list(1,2),:) = &
         best_gridpoint / real(bin_size,real12)
   
   end subroutine add_atom_scan
!!!#############################################################################

   
!!!#############################################################################
!!! add atom to unit cell considering the void space
!!!#############################################################################
  subroutine add_atom_void (bin_size, lattice, basis, atom_ignore_list, placed)
    implicit none
    type(bas_type), intent(inout) :: basis
    integer, dimension(3), intent(in) :: bin_size
    integer, dimension(:,:), intent(in) :: atom_ignore_list
    real(real12), dimension(3,3), intent(in) :: lattice
    logical, intent(out) :: placed
    
    integer :: i, j, k, l
    real(real12), dimension(3) :: best_location
    real(real12) :: best_location_bond, smallest_bond, comparison
    real(real12), dimension(3) :: tmpvector


    best_location_bond = -huge(1._real12)
   
    do i = 0, bin_size(1) - 1, 1
       do j = 0, bin_size(2) - 1, 1
          do k = 0, bin_size(3) - 1, 1
             tmpvector = [i, j, k] / real(bin_size,real12)
             smallest_bond = modu(get_min_dist(&
                  lattice, basis, tmpvector, .false., &
                  ignore_list = atom_ignore_list))
             if( smallest_bond .gt. best_location_bond ) then
                best_location_bond = smallest_bond
                best_location = tmpvector
             end if
          end do
       end do
    end do

    basis%spec(atom_ignore_list(1,1))%atom(atom_ignore_list(1,2),:) = &
         best_location
    placed = .true.

  end subroutine add_atom_void
!!!#############################################################################


!!!#############################################################################
!!! add atom to unit cell using a pseudo-random walk method
!!!#############################################################################
  subroutine add_atom_pseudo (bin_size, gvector_container, &
       lattice, basis, atom_ignore_list, &
       radius_arr, placed)
    implicit none
    type(gvector_container_type), intent(in) :: gvector_container
    type(bas_type), intent(inout) :: basis
    logical, intent(out) :: placed
    integer, dimension(3), intent(in) :: bin_size
    integer, dimension(:,:), intent(in) :: atom_ignore_list
    real(real12), dimension(3,3), intent(in) :: lattice
    real(real12), dimension(:,:,:), intent(in) :: radius_arr

    integer :: i, j, k, l
    real(real12) :: rtmp1, crude_norm
    real(real12) :: calculated_value, calculated_test
    real(real12), dimension(3) :: tmpvector, testvector
   

    !! test a random point in the unit cell
    i = 0
    placed = .false.
    crude_norm = 0._real12
100 random_loop : do 
       i = i + 1
       write(6,'(A)',ADVANCE='NO') achar(13)
       write(6,'(I5.0, A)', ADVANCE='NO') i
      
       do j = 1, 3
          call random_number(rtmp1)
          tmpvector(j) = rtmp1
       end do

       calculated_value = buildmap_POINT( gvector_container, &
            tmpvector, lattice, basis, &
            atom_ignore_list, radius_arr, &
            1.1_real12, 0.95_real12)
     
       call random_number(rtmp1)
       if (rtmp1.lt.calculated_value) then
          placed = .TRUE.
          exit random_loop
       end if
 
       !! NOTE: HARDCODED LIMIT ON NUMBER OF TRIES. SET IN INPUT FILE
       if(i.ge.10000) return
       if(calculated_value.eq.0) cycle
    end do random_loop


    l = 0
    !! if we have found a point, we can now walk around it
    walk_loop : do
       do j=1, 3
          call random_number(rtmp1)
          if(k.gt.10) then 
             testvector(j) = tmpvector(j) + &
                  ( rtmp1 * 2._real12 - 1._real12 ) * &
                  0.01_real12
          else
             testvector(j) = tmpvector(j) + &
             ( rtmp1 * 2._real12 - 1._real12 ) * &
             0.1_real12
          end if 
       end do
       testvector = testvector - floor(testvector)

       calculated_test = buildmap_POINT( gvector_container, &
            testvector, lattice, basis, &
            atom_ignore_list, radius_arr, &
            1.1_real12, 0.95_real12)
     
       if(calculated_test.lt.calculated_value) then 
          l = l + 1
          if(l.ge.10) then
             call random_number(rtmp1)
             if(crude_norm.lt.calculated_value) crude_norm = calculated_value

             if (rtmp1.lt.calculated_value) then
                placed = .TRUE.
                exit walk_loop
             end if
             if(k.gt.10) then 
                calculated_value = calculated_value / crude_norm
                if (rtmp1.lt.calculated_value) then
                   placed = .TRUE.
                   exit walk_loop
                end if
             end if
   
             !! if we have tried 10 times, and still no luck, then we need to ...
             !! ... reduce the tolerance
             k = k + 1
          end if   
          cycle walk_loop
       end if

       l=0
       tmpvector = testvector
       calculated_value = calculated_test
       call random_number(rtmp1)
 
       if(k.gt.10) calculated_test = calculated_test / crude_norm
       if (rtmp1.lt.calculated_test) then 
          placed=.TRUE.
          tmpvector = testvector
          exit walk_loop
       end if
 
    end do walk_loop
 
    basis%spec(atom_ignore_list(1,1))%atom(atom_ignore_list(1,2),:) = tmpvector
   
  end subroutine add_atom_pseudo
!!!#############################################################################

end module add_atom