module add_atom
  use constants, only: real12
  use atomtype
  implicit none

contains


!!!#############################################################################
!!! add atom to unit cell using the scan method (v2????)
!!!#############################################################################
  subroutine add_atom_scan_2 (bin_size,formula,atomlist,alistrep,atom_number_previous,structures,elrad,&
    &leng, results_matrix,eltot,elnames,placed,num_VOID)
   character(1024) :: name
   character(3), dimension(:), allocatable :: elnames
   type(unitcell), dimension(:), allocatable :: formula
   integer :: el_correct_read,i, j, k,n,l, num_VOID,atom_number_previous, new_atom_number, structures, leng,eltot
   integer, dimension(3) :: bin_size
   real(real12), dimension(3) :: best_location
   real(real12) :: best_location_value, distribution, sigma1, sigma2, bondcutoff&
      &,agausssamp
   real(real12), dimension(3) :: tmpvector
   type (atom), dimension(:,:), allocatable :: atomlist, alistrep, predicted_positions
   real(real12), dimension(:,:,:), allocatable :: elrad, product_matrix
   logical :: new_position_needed, placed
   real(real12), dimension(:,:,:,:,:), allocatable :: results_matrix, append_matrix
   real(real12), dimension(:,:), allocatable :: sorting_matrix
   
   open(55,file="scan_history",Access='append')
   
   allocate(sorting_matrix((bin_size(1)+1)*(bin_size(2)+1)*(bin_size(3)+1),4))
   results_matrix=0
   call buildmap_WIP (bin_size, formula, atomlist, alistrep, atom_number_previous, structures, &
      &elrad, leng, results_matrix,eltot,elnames,placed,num_VOID,append_matrix)
   if(placed.eqv..FALSE.) return
   n=0
   l=0
   
   
   
   do el_correct_read=1, eltot
    do i=0, bin_size(1)
       do j=0, bin_size(2)
          do k=0, bin_size(3)
             if(el_correct_read.ne.atomlist(structures,atom_number_previous+1)%element_index) then 
                cycle
             else
                n=n+1
             end if
             sorting_matrix(n,:)=results_matrix(i+1,j+1,k+1,:,el_correct_read)
             !write(*,*) results_matrix(i+1,j+1,k+1,:,el_correct_read)
             shuttle:do L=0,n
                if(n-L.gt.1) then 
                   if(sorting_matrix(n-L,4).gt.sorting_matrix(n-L-1,4)) then
                      sorting_matrix(n-L,:)=sorting_matrix(n-1-L,:)
                      sorting_matrix(n-1-L,:)=results_matrix(i+1,j+1,k+1,:,el_correct_read)
                   else
                      sorting_matrix(n-L,:)=results_matrix(i+1,j+1,k+1,:,el_correct_read)
                      exit shuttle
                   end if
                end if
             end do shuttle
          end do
       end do
    end do
   end do
   
   
   
   
   !write(*,*) "-------------------------------------------------------------"
   do i=1, (bin_size(1)+1)*(bin_size(2)+1)*(bin_size(3)+1)
    write(55,*) structures, atom_number_previous+1, sorting_matrix(i,:)
   end do
   !write(*,*) "-------------------------------------------------------------"
   
   do i=1, 11
    do j=1, 11
       do k=1, 11
          if((append_matrix&
               &(i,j,k,4,atomlist(structures,atom_number_previous+1)%element_index)).eq.&
               &maxval(append_matrix&
               &(:,:,:,4,atomlist(structures,atom_number_previous+1)%element_index))) then 
   
   
             do L=1, 3
                atomlist(structures,atom_number_previous+1)%position(L)=&
                     &append_matrix&
                     &(i,j,k,L,atomlist(structures,atom_number_previous+1)%element_index)
             end do
   
          end if
       end do
    end do
   end do
   write(name,'(A,I0)') "atoms",atom_number_previous+1
   open(1062,file=trim(adjustl(name)))
   write(1062,*) "x,y,z,data,label,color"
   
   do i=1, atom_number_previous+1
    if(trim(adjustl(atomlist(structures,i)%name)).eq."Si") then
       write(1062,*) atomlist(structures,i)%position(1), ",",&
            &atomlist(structures,i)%position(2), ",",&
            &atomlist(structures,i)%position(3), ",",&
            &"1,",&
            &trim(adjustl(atomlist(structures,i)%name)), ",",&
            &"blue"
    else 
       write(1062,*) atomlist(structures,i)%position(1), ",",&
            &atomlist(structures,i)%position(2), ",",&
            &atomlist(structures,i)%position(3), ",",&
            &"1,",&
            &trim(adjustl(atomlist(structures,i)%name)), ",",&
            &"green"
    end if
   end do
   close(1062)
   do i=1, (bin_size(1)+1)*(bin_size(2)+1)*(bin_size(3)+1)
      !write(*,*) sorting_matrix(i,:), "MAT"
   end do
   call atomrepeater(structures,atomlist(structures,atom_number_previous+1)%position,&
      &atomlist(structures,atom_number_previous+1)%name,alistrep,&
      &formula,atom_number_previous+1,leng)
   
   close(55)
   
   
   end subroutine add_atom_scan_2
!!!#############################################################################
   
   
!!!#############################################################################
!!! add atom to unit cell using the scan method
!!!#############################################################################
   !This routine needs to consider the effects of PLACING the atom that it might interact with itself. Hard to do
   subroutine add_atom_scan (bin_size,formula,atom_number_previous, sigma1,&
    &structures,sigma2,elnames,eltot,bondcutoff,atomlist,alistrep,tmpvector,elrad,best_location)
   character(3), dimension(:), allocatable :: elnames
   type(unitcell), dimension(:), allocatable :: formula
   integer :: scan_counter, bin_size,p,q, i, j, k, best_location_index, eltot, atom_number_previous, new_atom_number, structures
   real(real12), dimension(3) :: best_location
   real(real12) :: best_location_value, distribution, sigma1, sigma2, bond_distribution, angle_distribution, bondcutoff&
      &,agausssamp
   real(real12), dimension(3) :: tmpvector
   type (atom), dimension(:,:), allocatable :: atomlist, alistrep, predicted_positions
   real(real12), dimension(:,:,:), allocatable :: elrad
   logical :: new_position_needed
   allocate(predicted_positions(1,26))
   new_atom_number=atom_number_previous+1
   call atomprojector(tmpvector,predicted_positions,formula,new_atom_number,structures)
   write(*,*) "Scanning unit cell for a good location"
   best_location_value=0
   best_location_index=0
   p=0
   do i=0, bin_size-1
    do j=0, bin_size-1
       firstloop :do k=0, bin_size-1 
          tmpvector(1)=i*abs(formula(structures)%cell(1,1)/bin_size)
          tmpvector(2)=i*abs(formula(structures)%cell(2,1)/bin_size)
          tmpvector(3)=i*abs(formula(structures)%cell(3,1)/bin_size)
          tmpvector(1)=tmpvector(1)+j*abs(formula(structures)%cell(1,2)/bin_size)
          tmpvector(2)=tmpvector(2)+j*abs(formula(structures)%cell(2,2)/bin_size)
          tmpvector(3)=tmpvector(3)+j*abs(formula(structures)%cell(3,2)/bin_size)
          tmpvector(1)=tmpvector(1)+k*abs(formula(structures)%cell(1,3)/bin_size)
          tmpvector(2)=tmpvector(2)+k*abs(formula(structures)%cell(2,3)/bin_size)
          tmpvector(3)=tmpvector(3)+k*abs(formula(structures)%cell(3,3)/bin_size)
   
          do q=1, atom_number_previous !!!This should really be a *27, although is probably redundant
             if (get_bondlength(tmpvector,alistrep(structures,q)%position).lt.1.0) cycle firstloop
   !!! Experimental section for ruling out areas of unit cell 
             !if(tmpvector(3).lt.0.375*formula(structures)%cell(3,3)) cycle firstloop
             !if(tmpvector(3).gt.0.625*formula(structures)%cell(3,3)) cycle firstloop
          end do
          call atomprojector(tmpvector,predicted_positions,formula,new_atom_number,structures)
   
          !First, we consider the bondlength match to the existing framework. 
          distribution=0
          call generate_bondlength_distribution(sigma1,sigma2,new_atom_number,distribution,&
               &atomlist,alistrep,tmpvector,elrad,eltot,structures,elnames,new_position_needed,&
               &predicted_positions)
          bond_distribution=distribution
          if(new_position_needed.eqv..TRUE.) cycle
          ! THIS WILL NOT TRIGGER IN TWO ATOM SYSTEMS; it should do though, as there are repeated atoms"
   
          !call generate_bondangle_distribution(tmpvector,atomlist,alistrep,structures,bondcutoff,new_atom_number&
          !&,new_position_needed,agausssamp,distribution,eltot,sigma1,sigma2,elnames,predicted_positions)
          !angle_distribution=distribution
          if(new_position_needed.eqv..TRUE.) cycle
          p=p+1
          distribution=bond_distribution!angle_distribution*bond_distribution
   
          if (distribution.gt.best_location_value) then 
             best_location_value=distribution
             best_location=tmpvector
          end if
       end do firstloop
    end do
   end do
   do i=1, (new_atom_number-1)*27
      !     print*,"!", get_bondlength(best_location,alistrep(structures,i)%position)
   end do
   end subroutine add_atom_scan
!!!#############################################################################
   
   
!!!#############################################################################
!!! add atom to unit cell considering the void space
!!!#############################################################################
   !This routine needs to consider the effects of PLACING the atom that it might interact with itself. Hard to do
   subroutine add_atom_void (bin_size,formula,atom_number_previous, sigma1,&
    &structures,sigma2,elnames,eltot,bondcutoff,atomlist,alistrep,tmpvector,elrad,leng)
   character(3), dimension(:), allocatable :: elnames
   type(unitcell), dimension(:), allocatable :: formula
   integer :: scan_counter,leng,p,q, i, j, k, best_location_index, eltot, atom_number_previous&
      &, new_atom_number, structures
   integer, dimension(3) :: bin_size
   real(real12), dimension(3) :: best_location
   real(real12) :: best_location_value, smallest_bond, sigma1, sigma2, bond_distribution, angle_distribution, bondcutoff&
      &,agausssamp, comparison
   real(real12), dimension(3) :: tmpvector
   type (atom), dimension(:,:), allocatable :: atomlist, alistrep, predicted_positions
   real(real12), dimension(:,:,:), allocatable :: elrad
   logical :: new_position_needed
   
   
   
   write(*,*) atom_number_previous
   allocate(predicted_positions(1,26))
   new_atom_number=atom_number_previous+1
   call atomprojector(tmpvector,predicted_positions,formula,new_atom_number,structures)
   best_location_value=0
   best_location_index=0
   p=0
   open(1001,file="void_heatmap.txt")
   
   
   do i=0, bin_size(1)
    do j=0, bin_size(2)
       firstloop :do k=0, bin_size(3)
   
          tmpvector(1)=i*dble(formula(structures)%cell(1,1)/bin_size(1))
          tmpvector(2)=i*dble(formula(structures)%cell(2,1)/bin_size(1))
          tmpvector(3)=i*dble(formula(structures)%cell(3,1)/bin_size(1))
          tmpvector(1)=tmpvector(1)+j*dble(formula(structures)%cell(1,2)/bin_size(2))
          tmpvector(2)=tmpvector(2)+j*dble(formula(structures)%cell(2,2)/bin_size(2))
          tmpvector(3)=tmpvector(3)+j*dble(formula(structures)%cell(3,2)/bin_size(2))
          tmpvector(1)=tmpvector(1)+k*dble(formula(structures)%cell(1,3)/bin_size(3))
          tmpvector(2)=tmpvector(2)+k*dble(formula(structures)%cell(2,3)/bin_size(3))
          tmpvector(3)=tmpvector(3)+k*dble(formula(structures)%cell(3,3)/bin_size(3)) 
          if((k.eq.0).and.(j.eq.0).and.(i.eq.0)) smallest_bond=get_bondlength(tmpvector,alistrep(structures,1)%position)
          p=0
          do q=1, atom_number_previous*27
             !if (get_bondlength(tmpvector,alistrep(structures,q)%position).lt.1.0) cycle firstloop
   !!! Experimental section for ruling out areas of unit cell
             !if(tmpvector(3).lt.0.375*formula(structures)%cell(3,3)) cycle firstloop
             !if(tmpvector(3).gt.0.625*formula(structures)%cell(3,3)) cycle firstloop
             if(p.eq.0) smallest_bond=get_bondlength(tmpvector,alistrep(structures,1)%position)
             p=p+1
             if (get_bondlength(tmpvector,alistrep(structures,q)%position).lt.smallest_bond) then 
                smallest_bond=get_bondlength(tmpvector,alistrep(structures,q)%position)
             end if
          end do
          comparison=100.0*dble(dble(k)/dble(bin_size(3)))
          if(comparison.ge.c_cut) then
             !write(*,*) comparison, c_cut, "HERE"
             smallest_bond=-1
   
          end if
   
          if(smallest_bond.eq.-1) then 
   
             write(1001,*) tmpvector(1), tmpvector(2), tmpvector(3), 0
          else 
             write(1001,*) tmpvector(1), tmpvector(2), tmpvector(3), smallest_bond
          end if
   
          call atomprojector(tmpvector,predicted_positions,formula,new_atom_number,structures)
   
          !  do q=atom_number_previous+2, 27*atom_number_previous
          !     if (get_bondlength(tmpvector,alistrep(structures,q)%position).lt.1.0) cycle firstloop
   !!! Experimental section for ruling out areas of unit cell
          !            if(tmpvector(3).lt.0.375*formula(structures)%cell(3,3)) cycle firstloop
          !!            if(tmpvector(3).gt.0.625*formula(structures)%cell(3,3)) cycle firstloop
          !    if (get_bondlength(tmpvector,alistrep(structures,q)%position).lt.smallest_bond) then
          !       smallest_bond=get_bondlength(tmpvector,alistrep(structures,q)%position)
          !    end if
          ! end do
   
   
   
          if (smallest_bond.gt.best_location_value) then
             best_location_value=smallest_bond
             best_location=tmpvector
          end if
       end do firstloop
    end do
   end do
   close(1001)
   atomlist(structures,atom_number_previous+1)%position=best_location
   !write(*,*) best_location, "$^%$"
   open(56,file="void_history",access='append')
   write(56,*) structures, atom_number_previous+1,atomlist(structures,atom_number_previous+1)%position
   close(56)
   call atomrepeater(structures,atomlist(structures,atom_number_previous+1)%position,&
      &atomlist(structures,atom_number_previous+1)%name,alistrep,formula,atom_number_previous+1,leng)
   !do i=1, (new_atom_number-1)*27
   !   print*,"!", get_bondlength(best_location,alistrep(structures,i)%position)
   !end do
   end subroutine add_atom_void
!!!#############################################################################


!!!#############################################################################
!!! add atom to unit cell using a pseudo-random method
!!!#############################################################################
   subroutine add_atom_random (formula,atom_number_previous, sigma1,&
    &structures,sigma2,elnames,eltot,bondcutoff,atomlist,alistrep,tmpvector,elrad)
   character(3), dimension(:), allocatable :: elnames
   type(unitcell), dimension(:), allocatable :: formula
   integer :: scan_counter, bin_size,p,q, i, j, k, best_location_index, eltot, atom_number_previous, new_atom_number, structures
   real(real12), dimension(3) :: best_location
   real(real12) :: best_location_value, distribution, sigma1, sigma2, bond_distribution, angle_distribution, bondcutoff&
      &,agausssamp,r,search_region
   real(real12), dimension(3) :: tmpvector
   type (atom), dimension(:,:), allocatable :: atomlist, alistrep, predicted_positions
   real(real12), dimension(:,:,:), allocatable :: elrad
   logical :: new_position_needed
   allocate(predicted_positions(1,26))
   new_atom_number=atom_number_previous+1
   p=0
   scan_counter=1
   search_region=0.25
   
   write(*,*) "Pseuodorandom placement"
   firstloop : do while (scan_counter.ne.0)
    do j=1, 3 
       call random_number(r)
       tmpvector(j)=r
       !        if(j.eq.3) tmpvector(j)=search_region*r
    end do
    tmpvector(:)=matmul(formula(structures)%cell,tmpvector(:))
    !     tmpvector(:)=tmpvector(:)+0.375*formula(structures)%cell(3,:)
    !     write(*,*) tmpvector
   
   
    do q=1, atom_number_previous
       if (get_bondlength(tmpvector,alistrep(structures,q)%position).lt.1.0) cycle firstloop
   !!! Experimental section for ruling out areas of unit cell
       if(tmpvector(3).lt.0.375*formula(structures)%cell(3,3)) cycle firstloop
       if(tmpvector(3).gt.0.625*formula(structures)%cell(3,3)) cycle firstloop
    end do
    call atomprojector(tmpvector,predicted_positions,formula,new_atom_number,structures)
   
    !First, we consider the bondlength match to the existing framework.
    distribution=0
    call generate_bondlength_distribution(sigma1,sigma2,new_atom_number,distribution,&
         &atomlist,alistrep,tmpvector,elrad,eltot,structures,elnames,new_position_needed,&
         &predicted_positions)
    bond_distribution=distribution
    if(new_position_needed.eqv..TRUE.) cycle firstloop
    ! THIS WILL NOT TRIGGER IN TWO ATOM SYSTEMS; it should do though, as there are repeated atoms"
   
    !   call generate_bondangle_distribution(tmpvector,atomlist,alistrep,structures,bondcutoff,new_atom_number&
    !        &,new_position_needed,agausssamp,distribution,eltot,sigma1,sigma2,elnames,predicted_positions)
    !   angle_distribution=distribution
    if(new_position_needed.eqv..TRUE.) cycle firstloop
    p=p+1
    if(p.eq.100) then 
       write(*,*) "100 attemps"
       p=p-100
    end if
    distribution=bond_distribution!angle_distribution*bond_distribution
    call random_number(r)
    if (r.gt.distribution) cycle firstloop 
    scan_counter=0
   
   end do firstloop
   
   end subroutine add_atom_random
!!!#############################################################################
   
   
!!!#############################################################################
!!! add atom to unit cell using a pseudo-random walk method
!!!#############################################################################
   subroutine add_atom_pseudo (bin_size,formula,atomlist,alistrep,atom_number_previous,structures,elrad,&
    &leng, results_matrix,eltot,elnames,placed,num_VOID)
   character(1024) :: name
   character(3), dimension(:), allocatable :: elnames
   type(unitcell), dimension(:), allocatable :: formula
   integer :: el_correct_read,i, j, k,n,l, num_VOID,atom_number_previous, new_atom_number, structures, leng,eltot
   integer, dimension(3) :: bin_size
   real(real12), dimension(3) :: best_location
   real(real12) :: best_location_value, distribution, sigma1, sigma2, bondcutoff&
      &,agausssamp, calculated_value, uptol, lowtol, r
   real(real12), dimension(3) :: tmpvector
   type (atom), dimension(:,:), allocatable :: atomlist, alistrep, predicted_positions
   real(real12), dimension(:,:,:), allocatable :: elrad, product_matrix
   logical :: new_position_needed, placed
   real(real12), dimension(:,:,:,:,:), allocatable :: results_matrix, append_matrix
   real(real12), dimension(:,:), allocatable :: sorting_matrix
   
   uptol=1.1
   lowtol=0.95
   
   results_matrix=0
   write(*,*) "here"
   infloop : do 
    do j=1, 3
       call random_number(r)
       tmpvector(j)=r 
    end do
    tmpvector(:)=matmul(formula(structures)%cell,tmpvector(:))
   
   
    call buildmap_POINT (tmpvector,formula,atomlist,alistrep&
         &,atom_number_previous,structures,elrad,leng,eltot, elnames,placed,num_VOID&
         &,uptol,lowtol,calculated_value)
    n=0
    l=0 
    call random_number(r)
   
    write(*,*) calculated_value, r, tmpvector
   
    if (r.lt.calculated_value) exit infloop
   end do infloop
   
   
   do L=1, 3
    atomlist(structures,atom_number_previous+1)%position(L)=tmpvector(L)
   
   end do
   
   write(name,'(A,I0)') "atoms",atom_number_previous+1
   open(1062,file=trim(adjustl(name)))
   write(1062,*) "x,y,z,data,label,color"
   
   do i=1, atom_number_previous+1
    if(trim(adjustl(atomlist(structures,i)%name)).eq."Si") then
       write(1062,*) atomlist(structures,i)%position(1), ",",&
            &atomlist(structures,i)%position(2), ",",&
            &atomlist(structures,i)%position(3), ",",&
            &"1,",&
            &trim(adjustl(atomlist(structures,i)%name)), ",",&
            &"blue"
    else 
       write(1062,*) atomlist(structures,i)%position(1), ",",&
            &atomlist(structures,i)%position(2), ",",&
            &atomlist(structures,i)%position(3), ",",&
            &"1,",&
            &trim(adjustl(atomlist(structures,i)%name)), ",",&
            &"green"
    end if
   end do
   close(1062)
   do i=1, (bin_size(1)+1)*(bin_size(2)+1)*(bin_size(3)+1)
      !write(*,*) sorting_matrix(i,:), "MAT"
   end do
   call atomrepeater(structures,atomlist(structures,atom_number_previous+1)%position,&
      &atomlist(structures,atom_number_previous+1)%name,alistrep,&
      &formula,atom_number_previous+1,leng)
   
   close(55)
   
   
   end subroutine add_atom_pseudo
!!!#############################################################################
   
   

end module add_atom