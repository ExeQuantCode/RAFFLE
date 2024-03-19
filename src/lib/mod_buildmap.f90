module buildmap
  use constants, only: real12
  use geom
  use vasp_file_handler, only: structurecounter, unitcell
  use contributions
  use atomtype
  implicit none

contains


subroutine buildmap_POINT(tmpvector,formula,atomlist,alistrep&
  &,atom_number_previous,structures,elrad,atom_total,eltot, elnames,placed,num_VOID&
  &,uptol,lowtol,calculated_value)

type(unitcell), dimension(:), allocatable :: formula
type (atom), dimension(:,:), allocatable :: atomlist,alistrep,predicted_positions
integer :: el,i,j,k,l,m,n, normalisation,atom_number_previous,structures,&
    &m_count, n_count, atom_total, eltot, cleanup, cleanup2, repeats, num_VOID, norm
real(real12) :: value_return, summation, max_1, max_2, max_3, comparison, uptol, lowtol, normaliser&
    &,repeat_power, i_comp, j_comp, k_comp, totbin, calculated_value
integer, dimension(3) :: sum, bin_size
real(real12), dimension(3) :: tmpvector, location_vector
real(real12), dimension(:,:,:), allocatable :: elrad
real(real12), dimension(:), allocatable :: f_body_matrix, th_body_matrix, t_body_matrix, product_matrix
real(real12), dimension(:,:), allocatable :: position_storage
integer, dimension(:,:), allocatable :: index_storage
character(3), dimension(:,:), allocatable :: name_storage, remaining_elements, remaining_2
real(real12), dimension(:,:,:,:,:), allocatable ::  results_matrix
logical :: file_exists, placed
character(3), dimension(:), allocatable :: elnames
character(5) :: APP
character(1024) :: name

allocate(f_body_matrix(maxval(atomlist(structures,:)%element_index)))
allocate(t_body_matrix(maxval(atomlist(structures,:)%element_index)))
allocate(th_body_matrix(maxval(atomlist(structures,:)%element_index)))
allocate(product_matrix(maxval(atomlist(structures,:)%element_index)))
allocate(predicted_positions(1,26))

f_body_matrix=1
t_body_matrix=0
th_body_matrix=1
product_matrix=1
calculated_value=0

allocate(position_storage(3,3))
allocate(index_storage(3,1))
allocate(name_storage(3,1))

!write(*,*) tmpvector, atomlist(structures,atom_number_previous+1)%name
repeat_power=1

elloop:do el=1, maxval(atomlist(structures,:)%element_index)
  if(atomlist(structures,atom_number_previous+1)%element_index.ne.el) cycle elloop            

  call atomprojector(tmpvector,predicted_positions,formula,atom_number_previous+1,structures)
  m_count=0
  n_count=0
  do L=1, atom_number_previous*27+26
     norm=1
     if(L.gt.atom_number_previous*27) then 
        position_storage(1,:)=predicted_positions(1,L-(atom_number_previous*27))%position(:)
        index_storage(1,1)=atomlist(structures,atom_number_previous+1)%element_index
        name_storage(1,1)=atomlist(structures,atom_number_previous+1)%name 

     else 

        position_storage(1,:)=alistrep(structures,L)%position(:)
        index_storage(1,1)=alistrep(structures,L)%element_index 
        name_storage(1,1)=alistrep(structures,L)%name
     end if


     value_return=0

     if(get_bondlength(tmpvector,&
          &position_storage(1,:)).lt.&
          &elrad(1,atomlist(structures,atom_number_previous+1)%element_index,index_storage(1,1))*lowtol) then 
        !write(*,*) tmpvector, get_bondlength(tmpvector,position_storage(1,:))
        t_body_matrix(el)=0
        f_body_matrix(el)=0
        th_body_matrix(el)=0

        !if((i.eq.2).and.(j.eq.2).and.(k.eq.5)) then
        !   write(*,*) "cycle: bondlength"
        !end if

        cycle elloop
     else if(get_bondlength(tmpvector,&
          &position_storage(1,:)).gt.&
          &elrad(1,atomlist(structures,atom_number_previous+1)%element_index,index_storage(1,1))*uptol) then
!!!! This has been left out to ease on computation. Could be reimplemented but increases cost dramatically to consider ALL atoms
        !call evaluate_contribution (trim(adjustl(&
        !     &atomlist(structures,atom_number_previous+1)%name)),&
        !     &trim(adjustl(alistrep(structures,atom_number_previous+1)%name)),get_bondlength(tmpvector,&
        !     &position_storage(1,:)),value_return)
        !if(get_bondlength(tmpvector,position_storage(1,:)).lt.4) then 
        !   write(*,*) tmpvector, get_bondlength(tmpvector,position_storage(1,:)), "!"
        !end if
        !write(*,*) "CYCLE0.5", elrad(1,atomlist(structures,atom_number_previous+1)%element_index,index_storage(1,1))*uptol&
        !     &,get_bondlength(tmpvector,&
        !     &position_storage(1,:))
        cycle
     else

        call evaluate_contribution (trim(adjustl(&
             &atomlist(structures,atom_number_previous+1)%name)),&
             &trim(adjustl(name_storage(1,1))),get_bondlength(tmpvector,&
             &position_storage(1,:)),value_return)

        t_body_matrix(el)=((t_body_matrix(el)*norm)+value_return)
        !ADD +1 to norm here to curtail contributions to bondlength. I do not think you need to, as it can dampen certain resonance points 
        norm=norm+1
        t_body_matrix(el)=t_body_matrix(el)/norm
        !write(*,*) tmpvector, t_body_matrix(el), el
        !if((i.eq.2).and.(j.eq.2).and.(k.eq.5)) then
        !   write(*,*) value_return, "2b"
        !end if


     end if
     if(get_bondlength(&
          &tmpvector,&
          &position_storage(1,:)).lt.&
          &elrad(1,atomlist(structures,atom_number_previous+1)%element_index,index_storage(1,1))*uptol) then
        do n=1, atom_number_previous*27+26
           if(n.eq.L) cycle
           if(n.gt.atom_number_previous*27) then
              position_storage(2,:)=predicted_positions(1,n-atom_number_previous*27)%position(:)
              index_storage(2,1)=atomlist(structures,atom_number_previous+1)%element_index
              name_storage(2,1)=atomlist(structures,atom_number_previous+1)%name
           else
              position_storage(2,:)=alistrep(structures,n)%position(:)
              index_storage(2,1)=alistrep(structures,n)%element_index
              name_storage(2,1)=alistrep(structures,n)%name
           end if
           if(get_bondlength(&
                &tmpvector,&
                &position_storage(2,:)).lt.&
                &elrad(1,atomlist(structures,atom_number_previous+1)%element_index,index_storage(2,1))*lowtol) then
              !if((i.eq.2).and.(j.eq.2).and.(k.eq.5)) then
              !   write(*,*) "cycle: bondlength ang"
              !end if

              t_body_matrix(el)=0
              f_body_matrix(el)=0
              th_body_matrix(el)=0
              !write(*,*) "CYCLE 2"
              cycle elloop
           end if
           if(get_bondlength(&
                &position_storage(1,:),&
                &position_storage(2,:)).lt.&
                &elrad(1,index_storage(1,1),index_storage(2,1))*uptol) then

              call evaluate_angle_contribution(trim(adjustl(&
                   &atomlist(structures,atom_number_previous+1)%name))&
                   &,get_bondangle(&
                   &tmpvector,&
                   &position_storage(1,:),&
                   &position_storage(2,:)),&
                   &value_return)

              th_body_matrix(el)=(th_body_matrix(el)*value_return**&
                   &(1.0/(repeat_power)))!*n_count)*value_return
              n_count=n_count+1
              th_body_matrix(el)=th_body_matrix(el)!/n_count
              !if((i.eq.2).and.(j.eq.2).and.(k.eq.5)) then
              !   write(*,*) value_return, "ang"
              !end if

           else if(get_bondlength(&
                &tmpvector,&
                &position_storage(2,:)).lt.&
                &elrad(1,atomlist(structures,atom_number_previous+1)%element_index,index_storage(2,1))*uptol) then 

              call evaluate_angle_contribution(trim(adjustl(&
                   &atomlist(structures,atom_number_previous+1)%name))&
                   &,get_bondangle(&
                   &position_storage(1,:),&
                   &tmpvector,&
                   &position_storage(2,:)),&
                   &value_return)

              th_body_matrix(el)=(th_body_matrix(el)*(value_return**&
                   (1.0/(repeat_power))))!n_count)+value_return
              n_count=n_count+1
              th_body_matrix(el)=th_body_matrix(el)!/n_count
              !if((i.eq.2).and.(j.eq.2).and.(k.eq.5)) then
              !   write(*,*) value_return, "ang"
              !end if

           end if
           !if((i.eq.2).and.(j.eq.2).and.(k.eq.5)) then
           !   write(*,*) get_bondlength(&
           !        &position_storage(1,:),&
           !        &position_storage(2,:)), "!"
           !
           !end if

           if((get_bondlength(&
                &position_storage(1,:),& 
                &position_storage(2,:)).lt.&
                &elrad(1,index_storage(1,1),index_storage(2,1))*uptol)&
                &.OR.(get_bondlength(&
                &tmpvector,&
                &position_storage(2,:)).lt.&
                &elrad(1,index_storage(1,1),index_storage(2,1))*uptol)) then

              do m=1, atom_number_previous*27+26 
                 if(m.eq.L) cycle 
                 if(m.eq.n) cycle
                 if(m.gt.atom_number_previous*27) then 
                    position_storage(3,:)=predicted_positions(1,m-atom_number_previous*27)%position(:)
                    index_storage(3,1)=atomlist(structures,atom_number_previous+1)%element_index
                    name_storage(3,1)=atomlist(structures,atom_number_previous+1)%name
                 else 
                    position_storage(3,:)=alistrep(structures,m)%position(:)
                    index_storage(3,1)=alistrep(structures,m)%element_index
                    name_storage(3,1)=alistrep(structures,m)%name
                 end if
                 if(get_bondlength(&
                      &tmpvector,&
                      &position_storage(3,:)).lt.&
                      elrad(1,atomlist(structures,atom_number_previous+1)%element_index,index_storage(3,1))*lowtol) then

                    f_body_matrix(el)=0
                    th_body_matrix(el)=0
                    t_body_matrix(el)=0

                    !if((i.eq.2).and.(j.eq.2).and.(k.eq.5)) then
                    !write(*,*) tmpvector, position_storage(3,:), get_bondlength(&
                    !     &tmpvector,&
                    !    &position_storage(3,:))
                    !end if
                    !write(*,*) "4 cycle"
                    cycle elloop
                 end if
                 if(get_bondlength(&
                      &position_storage(1,:),&
                      &position_storage(3,:)).lt.&
                      elrad(1,index_storage(1,1),index_storage(3,1))*uptol) then
                    call evaluate_4body_contribution(trim(adjustl(&
                         &atomlist(structures,atom_number_previous+1)%name))&
                         &,get_dihedral_angle(&
                         &tmpvector,&
                         &position_storage(1,:),&
                         &position_storage(2,:),&
                         &position_storage(3,:)&
                         &),&
                         &value_return)
                    !if((i.eq.2).and.(j.eq.2).and.(k.eq.5)) then
                    !write(*,*) f_body_matrix(el), value_return, "4b", index_storage(:,1)
                    !end if
                    !write(*,*) "!!!!!!!!!!!!!!!!!!!"
                    !write(*,*) value_return
                    !write(*,*) tmpvector
                    !write(*,*) position_storage(1,:)
                    !write(*,*) position_storage(2,:)
                    !write(*,*) position_storage(3,:)
                    !write(*,*) "!!!!!!!!!!!!!!!!!!!"     

                    if(value_return.eq.0) then

                       f_body_matrix(el)=0
                       th_body_matrix(el)=0
                       t_body_matrix(el)=0

                       !write(*,*) "CYCLE 4"
                       cycle elloop
                       !else
                       !print*,f_body_matrix(i+1,j+1,k+1,el), value_return, i,j,k
                    end if

                    !Here have taken a large root of value return, to account for sumamtion of atoms in 3D. 
                    !Will need to think further on this.

                    f_body_matrix(el)=(f_body_matrix(el)*((value_return)**&
                         &(1.0/(repeat_power))))
                    m_count=m_count+1
                    !f_body_matrix(i+1,j+1,k+1,el)=f_body_matrix(i+1,j+1,k+1,el)!/m_count
                    !write(*,*) f_body_matrix(i+1,j+1,k+1,el), value_return**(1.0/4.0),value_return,i,j,k
                 end if
              end do
           end if
        end do
     end if
  end do
  !write(*,*) f_body_matrix(i+1,j+1,k+1,el)
  if(m_count.eq.0) then 
     f_body_matrix(el)=1
     !write(*,*) "4_body set to zero on", tmpvector 
  end if
  if(n_count.eq.0) then 
     !write(*,*) "3_body set to one on", tmpvector
     th_body_matrix(el)=1
  end if
  product_matrix(el)=t_body_matrix(el)*&
       &f_body_matrix(el)*th_body_matrix(el)
  calculated_value=product_matrix(el)
  !write(*,*) product_matrix(i+1,j+1,k+1,el), f_body_matrix(i+1,j+1,k+1,el), i, j, k
  summation=product_matrix(el)
  !write(*,*) tmpvector, product_matrix(el),f_body_matrix(el), "!!!" 


end do elloop


end subroutine buildmap_POINT


subroutine buildmap_WIP (bin_size,formula,atomlist,alistrep&
  &,atom_number_previous,structures,elrad,atom_total,&
  &results_matrix,eltot, elnames,placed,num_VOID, append_matrix,c_cut,c_min)
implicit none
integer, intent(in) :: c_cut, c_min
type(unitcell), dimension(:), allocatable :: formula
type (atom), dimension(:,:), allocatable :: atomlist,alistrep,predicted_positions
integer :: el,i,j,k,l,m,n, normalisation,atom_number_previous,structures,&
    &m_count, n_count, atom_total, eltot, cleanup, cleanup2, repeats, num_VOID, norm
real(real12) :: value_return, summation, max_1, max_2, max_3, comparison, uptol, lowtol, normaliser&
    &,repeat_power, i_comp, j_comp, k_comp, totbin, calculated_value
integer, dimension(3) :: sum, bin_size, best
real(real12), dimension(3) :: tmpvector, location_vector, best_vector
real(real12), dimension(:,:,:), allocatable :: elrad
real(real12), dimension(:,:,:,:), allocatable :: update_region,f_body_matrix, th_body_matrix, t_body_matrix, product_matrix
real(real12), dimension(:,:), allocatable :: position_storage
integer, dimension(:,:), allocatable :: index_storage 
character(3), dimension(:,:), allocatable :: name_storage, remaining_elements, remaining_2 
real(real12), dimension(:,:,:,:,:), allocatable ::  results_matrix, append_matrix
logical :: file_exists, placed
character(3), dimension(:), allocatable :: elnames
character(5) :: APP
character(1024) :: name

uptol=1.1
lowtol=0.95



repeat_power=1


APP="APPEND"
write(*,*) "WELCOME TO THE NEW BUILDMAP FUNCTION; EXPERIMENTAL"
sum=bin_size
allocate(update_region(sum(1)+1,sum(2)+1,sum(3)+1, maxval(atomlist(structures,:)%element_index)))
allocate(predicted_positions(1,26))
allocate(f_body_matrix(sum(1)+1,sum(2)+1,sum(3)+1, maxval(atomlist(structures,:)%element_index)))
!!1 is value, 2 is contributions to value
allocate(t_body_matrix(sum(1)+1,sum(2)+1,sum(3)+1, maxval(atomlist(structures,:)%element_index)))
allocate(th_body_matrix(sum(1)+1,sum(2)+1,sum(3)+1, maxval(atomlist(structures,:)%element_index)))
allocate(position_storage(3,3))
allocate(index_storage(3,1))
allocate(name_storage(3,1))
allocate(product_matrix(sum(1)+1,sum(2)+1,sum(3)+1, maxval(atomlist(structures,:)%element_index)))
summation=0
allocate(remaining_elements(1, maxval(atomlist(structures,:)%element_index)))
allocate(remaining_2(1, maxval(atomlist(structures,:)%element_index)))
!write(*,*)  maxval(atomlist(structures,:)%element_index)

repeats=1
f_body_matrix=1
t_body_matrix=0
th_body_matrix=1
update_region=1
product_matrix=0
norm=0


INQUIRE(FILE="buildmap_testfile.txt", EXIST=file_exists)
if(file_exists) then 
  update_region=0
  write(*,*) "Loading in an existing buildmap. PLEASE ADD SPECIES SUPPORT"
  open(6969,file="buildmap_testfile.txt")
  n=1
  remaining_elements(1,1)=atomlist(structures,atom_number_previous+1)%name
  do cleanup=atom_number_previous+1, atom_total

     if(trim(adjustl(atomlist(structures,cleanup)%name)).ne.&
          trim(adjustl(remaining_elements(1,n)))) then
        n=n+1
        remaining_2=remaining_elements
        deallocate(remaining_elements) 
        allocate(remaining_elements(1,n))
        do i=1, n-1
           remaining_elements(1,i)=remaining_2(1,i)
        end do
        deallocate(remaining_2)
        allocate(remaining_2(1,n))
        remaining_elements(1,n)=atomlist(structures,cleanup)%name

     end if
  end do
  eloop: do m=1,  maxval(atomlist(structures,:)%element_index)
     do cleanup=1, n
        if(elnames(m).eq.remaining_elements(1,cleanup)) exit
        if(cleanup.eq.n) then 
           write(*,*) "CLEANUP CYCLE"
           cycle eloop
        end if
     end do


     do i=0, sum(1) 
        do j=0, sum(2)
           do k=0, sum(3) 
              if(atomlist(structures,atom_number_previous+1)%element_index.lt.m) cycle
              read(6969,*) location_vector,&
                   &product_matrix(i+1,j+1,k+1,atomlist(structures,atom_number_previous+1)%element_index),&
                   &f_body_matrix(i+1,j+1,k+1,atomlist(structures,atom_number_previous+1)%element_index),&
                   &t_body_matrix(i+1,j+1,k+1,atomlist(structures,atom_number_previous+1)%element_index),&
                   &th_body_Matrix(i+1,j+1,k+1,atomlist(structures,atom_number_previous+1)%element_index)
              if(atomlist(structures,atom_number_previous+1)%element_index.ne.m) then 
                 !write(*,*) "INDEX CYCLING"
                 cycle
              end if

              !Chops off the top of the unit cell, should be migrated to the Infile and extended to 3D
              !if(location_vector(3).gt.11) then 
              !   update_region(i+1,j+1,k+1,m)=0
              !end if


              logic_loop : do L=(atom_number_previous-num_VOID)*27+1,(atom_number_previous)*27
                 if(get_bondlength(location_vector,alistrep(structures,L)%position)&
                      &.gt.(elrad(3,m,&
                      &alistrep(structures,L)%element_index)*2.0)) then
                    update_region(i+1,j+1,k+1,m)=0
                    !write(*,*) i, j, k
                 end if
              end do logic_loop

              logic_loop2 : do L=(atom_number_previous-num_VOID)*27+1,(atom_number_previous)*27
                 if(get_bondlength(location_vector,alistrep(structures,L)%position)&
                      &.lt.(elrad(3,m,&
                      &alistrep(structures,L)%element_index)*2.0)) then 
                    !write(*,*) i, j, k , m,1
                    update_region(i+1,j+1,k+1,m)=1
                    exit logic_loop2

                 end if

              end do logic_loop2
           end do
        end do
     end do
  end do eloop
  rewind(6969)
else


  open(6969,file="buildmap_testfile.txt")
end if
product_matrix=0

do i=0, sum(1)
  do j=0, sum(2)
     do k=0, sum(3)
        !write(*,*) i, j, k , update_region(i+1,j+1,k+1,1)
     end do
  end do
end do

write(*,*) "Arriving here", repeats

100 if(repeats.ne.1) then
  write(*,*) "WIPING"
  !deallocate(update_region)
  !deallocate(f_body_matrix)
  !deallocate(t_body_matrix)
  !deallocate(th_body_matrix)
  !deallocate(product_matrix)
  !deallocate(results_matrix)
  deallocate(append_matrix)
  close(6969)
  call execute_command_line("rm buildmap_testfile.txt",WAIT=.TRUE.)
  open(6969,file="buildmap_testfile.txt")
  !sum=nint(sum*1.5)
  !allocate(results_matrix(sum+1,sum+1,sum+1,4, maxval(atomlist(structures,:)%element_index)))
  !allocate(update_region(sum+1,sum+1,sum+1, maxval(atomlist(structures,:)%element_index)))
  !allocate(f_body_matrix(sum+1,sum+1,sum+1, maxval(atomlist(structures,:)%element_index)))
  !!1 is value, 2 is contributions to value
  !allocate(t_body_matrix(sum+1,sum+1,sum+1, maxval(atomlist(structures,:)%element_index)))
  !allocate(th_body_matrix(sum+1,sum+1,sum+1, maxval(atomlist(structures,:)%element_index)))
  !allocate(product_matrix(sum+1,sum+1,sum+1, maxval(atomlist(structures,:)%element_index)))
  f_body_matrix=1
  t_body_matrix=0
  th_body_matrix=1
  update_region=1
  product_matrix=0

end if



!write(*,*) sum 
do i=0,sum(1); do j=0,sum(2); do k=0,sum(3)
  write(6,'(A)',ADVANCE='NO') achar(13)
  i_comp=i*(sum(2)+1.0)*(sum(3)+1.0)
  j_comp=j*(sum(3)+1)
  k_comp=k
  totbin=(sum(1)+1)*(sum(2)+1)*(sum(3)+1)
  write(6,'(I3.0, A)', ADVANCE='NO') NINT(100*(i_comp+j_comp+k_comp)/totbin), "%" 
  comparison=100.0*dble(k)/dble(sum(3))
  !write(*,*) i, j, k
  if(comparison.gt.c_cut) then 
     f_body_matrix(i+1,j+1,k+1,:)=0
     th_body_matrix(i+1,j+1,k+1,:)=0
     update_region(i+1,j+1,k+1,:)=0
     t_body_matrix(i+1,j+1,k+1,:)=0
     !print *, "CUTCYCLE"
     cycle
  end if
  if(comparison.lt.c_min) then
     f_body_matrix(i+1,j+1,k+1,:)=0
     th_body_matrix(i+1,j+1,k+1,:)=0
     update_region(i+1,j+1,k+1,:)=0
     t_body_matrix(i+1,j+1,k+1,:)=0
     !print *, "CUTCYCLE"
     cycle
  end if
  elloop: do el=1, maxval(atomlist(structures,:)%element_index)
     if((f_body_matrix(i+1,j+1,k+1,el).eq.0)) then 
        !write(*,*) "cycling: f body doesn't need updating", i,j,k
        !cycle elloop
        continue
     end if
     if((th_body_matrix(i+1,j+1,k+1,el).eq.0)) then 
        !write(*,*) "CYCLE B", i,j,k
        !cycle elloop
        continue
     end if
     if(update_region(i+1,j+1,k+1,el).ne.1) then
        !write(*,*) "CYCLE C",update_region(i+1,j+1,k+1,el),update_region(i+1,j+1,k+1,1),el, i, j, k
        cycle elloop
     end if
     if(atomlist(structures,atom_number_previous+1)%element_index.ne.el) cycle elloop            

     tmpvector(:)=i*dble(formula(structures)%cell(1,:)/(sum(1)))
     tmpvector(1)=tmpvector(1)+j*dble(formula(structures)%cell(2,1)/(sum(2)))
     tmpvector(2)=tmpvector(2)+j*dble(formula(structures)%cell(2,2)/(sum(2)))
     tmpvector(3)=tmpvector(3)+j*dble(formula(structures)%cell(2,3)/(sum(2)))
     tmpvector(1)=tmpvector(1)+k*dble(formula(structures)%cell(3,1)/(sum(3)))
     tmpvector(2)=tmpvector(2)+k*dble(formula(structures)%cell(3,2)/(sum(3)))
     tmpvector(3)=tmpvector(3)+k*dble(formula(structures)%cell(3,3)/(sum(3)))

     !This chops. Do not enable chop unless you want to chop
     !if(tmpvector(3).gt.11) then 
     !   f_body_matrix(i+1,j+1,k+1,:)=0
     !   th_body_matrix(i+1,j+1,k+1,:)=0
     !   update_region(i+1,j+1,k+1,:)=0
     !   t_body_matrix(i+1,j+1,k+1,:)=0
     !   cycle elloop
     !end if


     call atomprojector(tmpvector,predicted_positions,formula,atom_number_previous+1,structures)
     m_count=0
     n_count=0
     do L=1, atom_number_previous*27+26


        norm=1
        !cutoff here should be made sensibly for how close exclusion range is to each atom
        !write(*,*) "-----------------------------------"
        if(L.gt.atom_number_previous*27) then 
           position_storage(1,:)=predicted_positions(1,L-(atom_number_previous*27))%position(:)
           index_storage(1,1)=atomlist(structures,atom_number_previous+1)%element_index
           name_storage(1,1)=atomlist(structures,atom_number_previous+1)%name 

        else 

           position_storage(1,:)=alistrep(structures,L)%position(:)
           index_storage(1,1)=alistrep(structures,L)%element_index 
           name_storage(1,1)=alistrep(structures,L)%name
        end if


        value_return=0
        if(get_bondlength(tmpvector,&
             &position_storage(1,:)).lt.&
             &elrad(1,atomlist(structures,atom_number_previous+1)%element_index,index_storage(1,1))*lowtol) then 
           !write(*,*) tmpvector, get_bondlength(tmpvector,position_storage(1,:))
           t_body_matrix(i+1,j+1,k+1,el)=0
           f_body_matrix(i+1,j+1,k+1,el)=0
           th_body_matrix(i+1,j+1,k+1,el)=0

           !if((i.eq.2).and.(j.eq.2).and.(k.eq.5)) then
           !   write(*,*) "cycle: bondlength"
           !end if

           !write(*,*) "CYCLE 1", i, j, k
           cycle elloop
        else if(get_bondlength(tmpvector,&
             &position_storage(1,:)).gt.&
             &elrad(1,atomlist(structures,atom_number_previous+1)%element_index,index_storage(1,1))*uptol) then
!!!! This has been left out to ease on computation. Could be reimplemented but increases cost dramatically to consider ALL atoms
           !call evaluate_contribution (trim(adjustl(&
           !     &atomlist(structures,atom_number_previous+1)%name)),&
           !     &trim(adjustl(alistrep(structures,atom_number_previous+1)%name)),get_bondlength(tmpvector,&
           !     &position_storage(1,:)),value_return)
           !if(get_bondlength(tmpvector,position_storage(1,:)).lt.4) then 
           !   write(*,*) tmpvector, get_bondlength(tmpvector,position_storage(1,:)), "!"
           !end if
           !write(*,*) "CYCLE0.5", elrad(1,atomlist(structures,atom_number_previous+1)%element_index,index_storage(1,1))*uptol&
           !     &,get_bondlength(tmpvector,&
           !     &position_storage(1,:))
           cycle
        else

           call evaluate_contribution (trim(adjustl(&
                &atomlist(structures,atom_number_previous+1)%name)),&
                &trim(adjustl(name_storage(1,1))),get_bondlength(tmpvector,&
                &position_storage(1,:)),value_return)

           t_body_matrix(i+1,j+1,k+1,el)=((t_body_matrix(i+1,j+1,k+1,el)*norm)+value_return)
           !ADD +1 to norm here to curtail contributions to bondlength. I do not think you need to, as it can dampen certain resonance points 
           norm=norm+1
           t_body_matrix(i+1,j+1,k+1,el)=t_body_matrix(i+1,j+1,k+1,el)/norm
           !write(*,*) tmpvector, t_body_matrix(i+1,j+1,k+1,el),i, j, k
           !if((i.eq.2).and.(j.eq.2).and.(k.eq.5)) then
           !   write(*,*) value_return, "2b"
           !end if


        end if
        if(get_bondlength(&
             &tmpvector,&
             &position_storage(1,:)).lt.&
             &elrad(1,atomlist(structures,atom_number_previous+1)%element_index,index_storage(1,1))*uptol) then
           do n=1, atom_number_previous*27+26
              if(n.eq.L) cycle
              if(n.gt.atom_number_previous*27) then
                 position_storage(2,:)=predicted_positions(1,n-atom_number_previous*27)%position(:)
                 index_storage(2,1)=atomlist(structures,atom_number_previous+1)%element_index
                 name_storage(2,1)=atomlist(structures,atom_number_previous+1)%name
              else
                 position_storage(2,:)=alistrep(structures,n)%position(:)
                 index_storage(2,1)=alistrep(structures,n)%element_index
                 name_storage(2,1)=alistrep(structures,n)%name
              end if
              if(get_bondlength(&
                   &tmpvector,&
                   &position_storage(2,:)).lt.&
                   &elrad(1,atomlist(structures,atom_number_previous+1)%element_index,index_storage(2,1))*lowtol) then
                 !if((i.eq.2).and.(j.eq.2).and.(k.eq.5)) then
                 !   write(*,*) "cycle: bondlength ang"
                 !end if

                 t_body_matrix(i+1,j+1,k+1,el)=0
                 f_body_matrix(i+1,j+1,k+1,el)=0
                 th_body_matrix(i+1,j+1,k+1,el)=0
                 !write(*,*) "CYCLE 2"
                 cycle elloop
              end if
              if(get_bondlength(&
                   &position_storage(1,:),&
                   &position_storage(2,:)).lt.&
                   &elrad(1,index_storage(1,1),index_storage(2,1))*uptol) then

                 call evaluate_angle_contribution(trim(adjustl(&
                      &atomlist(structures,atom_number_previous+1)%name))&
                      &,get_bondangle(&
                      &tmpvector,&
                      &position_storage(1,:),&
                      &position_storage(2,:)),&
                      &value_return)

                 th_body_matrix(i+1,j+1,k+1,el)=(th_body_matrix(i+1,j+1,k+1,el)*value_return**&
                      &(1.0/(repeat_power)))!*n_count)*value_return
                 n_count=n_count+1
                 th_body_matrix(i+1,j+1,k+1,el)=th_body_matrix(i+1,j+1,k+1,el)!/n_count
                 !if((i.eq.2).and.(j.eq.2).and.(k.eq.5)) then
                 !   write(*,*) value_return, "ang"
                 !end if

              else if(get_bondlength(&
                   &tmpvector,&
                   &position_storage(2,:)).lt.&
                   &elrad(1,atomlist(structures,atom_number_previous+1)%element_index,index_storage(2,1))*uptol) then 

                 call evaluate_angle_contribution(trim(adjustl(&
                      &atomlist(structures,atom_number_previous+1)%name))&
                      &,get_bondangle(&
                      &position_storage(1,:),&
                      &tmpvector,&
                      &position_storage(2,:)),&
                      &value_return)

                 th_body_matrix(i+1,j+1,k+1,el)=(th_body_matrix(i+1,j+1,k+1,el)*(value_return**&
                      (1.0/(repeat_power))))!n_count)+value_return
                 n_count=n_count+1
                 th_body_matrix(i+1,j+1,k+1,el)=th_body_matrix(i+1,j+1,k+1,el)!/n_count
                 !if((i.eq.2).and.(j.eq.2).and.(k.eq.5)) then
                 !   write(*,*) value_return, "ang"
                 !end if

              end if
              !if((i.eq.2).and.(j.eq.2).and.(k.eq.5)) then
              !   write(*,*) get_bondlength(&
              !        &position_storage(1,:),&
              !        &position_storage(2,:)), "!"
              !
              !end if

              if((get_bondlength(&
                   &position_storage(1,:),& 
                   &position_storage(2,:)).lt.&
                   &elrad(1,index_storage(1,1),index_storage(2,1))*uptol)&
                   &.OR.(get_bondlength(&
                   &tmpvector,&
                   &position_storage(2,:)).lt.&
                   &elrad(1,index_storage(1,1),index_storage(2,1))*uptol)) then

                 do m=1, atom_number_previous*27+26 
                    !write(*,*) L, m, n
                    if(m.eq.L) cycle 
                    if(m.eq.n) cycle
                    if(m.gt.atom_number_previous*27) then 
                       position_storage(3,:)=predicted_positions(1,m-atom_number_previous*27)%position(:)
                       index_storage(3,1)=atomlist(structures,atom_number_previous+1)%element_index
                       name_storage(3,1)=atomlist(structures,atom_number_previous+1)%name
                    else 
                       position_storage(3,:)=alistrep(structures,m)%position(:)
                       index_storage(3,1)=alistrep(structures,m)%element_index
                       name_storage(3,1)=alistrep(structures,m)%name
                    end if
                    if(get_bondlength(&
                         &tmpvector,&
                         &position_storage(3,:)).lt.&
                         elrad(1,atomlist(structures,atom_number_previous+1)%element_index,index_storage(3,1))*lowtol) then

                       f_body_matrix(i+1,j+1,k+1,el)=0
                       th_body_matrix(i+1,j+1,k+1,el)=0
                       t_body_matrix(i+1,j+1,k+1,el)=0

                       !if((i.eq.2).and.(j.eq.2).and.(k.eq.5)) then
                       !   write(*,*) tmpvector, position_storage(3,:), get_bondlength(&
                       !        &tmpvector,&
                       !        &position_storage(3,:))
                       !end if

                       cycle elloop
                    end if
                    if(get_bondlength(&
                         &position_storage(1,:),&
                         &position_storage(3,:)).lt.&
                         elrad(1,index_storage(1,1),index_storage(3,1))*uptol) then
                       call evaluate_4body_contribution(trim(adjustl(&
                            &atomlist(structures,atom_number_previous+1)%name))&
                            &,get_dihedral_angle(&
                            &tmpvector,&
                            &position_storage(1,:),&
                            &position_storage(2,:),&
                            &position_storage(3,:)&
                            &),&
                            &value_return)
                       !if((i.eq.2).and.(j.eq.2).and.(k.eq.5)) then
                       !   write(*,*) value_return, "4b"
                       !end if


                       if(m_count.eq.0) then 
                          normaliser=1
                       end if
                       if(value_return.eq.0) then

                          f_body_matrix(i+1,j+1,k+1,el)=0
                          th_body_matrix(i+1,j+1,k+1,el)=0
                          t_body_matrix(i+1,j+1,k+1,el)=0

                          !write(*,*) "CYCLE 4"
                          cycle elloop
                          !else
                          !print*,f_body_matrix(i+1,j+1,k+1,el), value_return, i,j,k
                       end if

                       !Here have taken a large root of value return, to account for sumamtion of atoms in 3D. 
                       !Will need to think further on this.

                       f_body_matrix(i+1,j+1,k+1,el)=(f_body_matrix(i+1,j+1,k+1,el)*((value_return)**&
                            &(1.0/(repeat_power))))
                       m_count=m_count+1
                       !f_body_matrix(i+1,j+1,k+1,el)=f_body_matrix(i+1,j+1,k+1,el)!/m_count
                       !write(*,*) f_body_matrix(i+1,j+1,k+1,el), value_return**(1.0/4.0),value_return,i,j,k
                    end if
                 end do
              end if
           end do
        end if
     end do
     !write(*,*) f_body_matrix(i+1,j+1,k+1,el)
     if(m_count.eq.0) then 
        f_body_matrix(i+1,j+1,k+1,el)=1
        !write(*,*) "4_body set to zero on", tmpvector 
     end if
     if(n_count.eq.0) then 
        !write(*,*) "3_body set to one on", tmpvector
        th_body_matrix(i+1,j+1,k+1,el)=1
     end if
     product_matrix(i+1,j+1,k+1,el)=t_body_matrix(i+1,j+1,k+1,el)*&
          &f_body_matrix(i+1,j+1,k+1,el)*th_body_matrix(i+1,j+1,k+1,el)
     !write(*,*) product_matrix(i+1,j+1,k+1,el), f_body_matrix(i+1,j+1,k+1,el), i, j, k
     if(product_matrix(i+1,j+1,k+1,el).gt.summation) then 
        summation=product_matrix(i+1,j+1,k+1,el)
        best(1)=i
        best(2)=j
        best(3)=k
     end if

  end do elloop
end do; end do; end do
!! 3 is arbitrary maxdensity of grid
l=0
do concurrent(i=0:sum(1), j=0:sum(2), k=0:sum(3)) 
  do el=1, maxval(atomlist(structures,:)%element_index)
     l=l+1
     product_matrix(i+1,j+1,k+1,el)=t_body_matrix(i+1,j+1,k+1,el)*&
          &f_body_matrix(i+1,j+1,k+1,el)*th_body_matrix(i+1,j+1,k+1,el)
     !write(*,*) th_body_matrix(i+1,j+1,k+1,el), t_body_matrix(i+1,j+1,k+1,el), f_body_matrix(i+1,j+1,k+1,el),i,j,k

     if(product_matrix(i+1,j+1,k+1,el).gt.summation) then 
        !write(*,*) "DID IT"
        summation=product_matrix(i+1,j+1,k+1,el)
        best(1)=i
        best(2)=j
        best(3)=k
     end if
  end do
end do

!write(*,*) product_matrix(3,3,6,1), f_body_matrix(3,3,6,1), th_body_matrix(3,3,6,1)

!write(*,*) summation

best_vector(1)=best(1)*dble(formula(structures)%cell(1,1)/(sum(1)))
best_vector(2)=best(1)*dble(formula(structures)%cell(1,2)/(sum(1)))
best_vector(3)=best(1)*dble(formula(structures)%cell(1,3)/(sum(1)))
best_vector(1)=best_vector(1)+best(2)*dble(formula(structures)%cell(2,1)/(sum(2)))
best_vector(2)=best_vector(2)+best(2)*dble(formula(structures)%cell(2,2)/(sum(2)))
best_vector(3)=best_vector(3)+best(2)*dble(formula(structures)%cell(2,3)/(sum(2)))
best_vector(1)=best_vector(1)+best(3)*dble(formula(structures)%cell(3,1)/(sum(3)))
best_vector(2)=best_vector(2)+best(3)*dble(formula(structures)%cell(3,2)/(sum(3)))
best_vector(3)=best_vector(3)+best(3)*dble(formula(structures)%cell(3,3)/(sum(3)))

!write(*,*) best_vector
!write(*,*) summation
results_matrix=0
m=0
n=1
l=0

!!! Append Marix should most likely not be related to the input paramater bin size, as for high resolutions this would result in too many loops. Have left this for testing. 

!!!UPdate. Buildmaps do not take into account append matrixes, so later results are wrong thanks to wriggle. Make sure to fix ASAP

allocate(append_matrix(11,11,11,4,maxval(atomlist(structures,:)%element_index)))

do i= -5, 5
l=l+1
m=0
do j= -5, 5 
m=m+1
n=0
do k= -5, 5
  n=n+1


  tmpvector(1)=best_vector(1)+dble(i)*formula(structures)%cell(1,1)/(100.0*dble(sum(1)))
  tmpvector(2)=best_vector(2)+dble(i)*formula(structures)%cell(1,2)/(100.0*dble(sum(1)))
  tmpvector(3)=best_vector(3)+dble(i)*formula(structures)%cell(1,3)/(100.0*dble(sum(1)))
  tmpvector(1)=tmpvector(1)+dble(j)*formula(structures)%cell(2,1)/(100.0*dble(sum(2)))
  tmpvector(2)=tmpvector(2)+dble(j)*formula(structures)%cell(2,2)/(100.0*dble(sum(2)))
  tmpvector(3)=tmpvector(3)+dble(j)*formula(structures)%cell(2,3)/(100.0*dble(sum(2)))
  tmpvector(1)=tmpvector(1)+dble(k)*formula(structures)%cell(3,1)/(100.0*dble(sum(3)))
  tmpvector(2)=tmpvector(2)+dble(k)*formula(structures)%cell(3,2)/(100.0*dble(sum(3)))
  tmpvector(3)=tmpvector(3)+dble(k)*formula(structures)%cell(3,3)/(100.0*dble(sum(3)))


  call buildmap_POINT(tmpvector,formula,atomlist,alistrep&
       &,atom_number_previous,structures,elrad,atom_total,eltot, elnames,placed,num_VOID&
       &,uptol,lowtol,calculated_value)

  append_matrix(L,m,n,1,atomlist(structures,atom_number_previous+1)%element_index)=&
       &tmpvector(1)
  append_matrix(L,m,n,2,atomlist(structures,atom_number_previous+1)%element_index)=&
       &tmpvector(2)
  append_matrix(L,m,n,3,atomlist(structures,atom_number_previous+1)%element_index)=&
       &tmpvector(3)
  append_matrix(L,m,n,4,atomlist(structures,atom_number_previous+1)%element_index)=&
       &calculated_value


end do
end do
end do

do i=1, 11
do j=1, 11
do k=1, 11
  if((append_matrix&
       &(i,j,k,4,atomlist(structures,atom_number_previous+1)%element_index)).eq.&
       &maxval(append_matrix&
       &(:,:,:,4,atomlist(structures,atom_number_previous+1)%element_index))) then


     do L=1, 3
        best(L)=&
             &append_matrix&
             &(i,j,k,L,atomlist(structures,atom_number_previous+1)%element_index)
     end do
     calculated_value=append_matrix(L,m,n,4,atomlist(structures,atom_number_previous+1)%element_index)
  end if
end do
end do
end do





if(summation.lt.10D-300) then 
write(*,*) "SUMMATION EFFECTIVELY 0, USE WITH CAUTION"
summation=0
end if


write(*,*) summation, "SUMMATION"
if(summation.eq.0) then
if(repeats.lt.10) then
repeats=repeats+1
repeat_power=repeat_power+5
write(*,*) "CYCLING"

GO TO 100
else 
max_1=0
max_2=0
max_3=0
do i=0, sum(1)
  do j=0, sum(2) 
     do k=0, sum(3) 
        if(t_body_matrix(i+1,j+1,k+1,atomlist(structures,atom_number_previous+1)%element_index).gt.max_1) then 
           max_1=t_body_matrix(i+1,j+1,k+1,atomlist(structures,atom_number_previous+1)%element_index)
        end if
        if(th_body_matrix(i+1,j+1,k+1,atomlist(structures,atom_number_previous+1)%element_index).gt.max_2) then
           max_2=th_body_matrix(i+1,j+1,k+1,atomlist(structures,atom_number_previous+1)%element_index)
        end if
        if(f_body_matrix(i+1,j+1,k+1,atomlist(structures,atom_number_previous+1)%element_index).gt.max_3) then
           max_3=f_body_matrix(i+1,j+1,k+1,atomlist(structures,atom_number_previous+1)%element_index)
        end if
     end do
  end do
end do
write(*,*) "BREAKING", max_1, max_2, max_3, max_1*max_2*max_3
!Set eq.0 to activate
if(max_2.eq.-1) then 
  write(*,*) "resetting th"
  do el=1,  maxval(atomlist(structures,:)%element_index)
     do i=0, sum(1)
        do j=0, sum(2)
           do k=0, sum(3)
              th_body_matrix(i+1,j+1,k+1,atomlist(structures,atom_number_previous+1)%element_index)=1.0
              product_matrix(i+1,j+1,k+1,el)=t_body_matrix(i+1,j+1,k+1,el)*&
                   &f_body_matrix(i+1,j+1,k+1,el)*th_body_matrix(i+1,j+1,k+1,el)
              !           write(*,*) product_matrix(i+1,j+1,k+1,el)
              if(summation.lt.product_matrix(i+1,j+1,k+1,el)) summation=product_matrix(i+1,j+1,k+1,el)
           end do
        end do
     end do
  end do
end if
!Set eq.0 to activate
if(max_3.eq.-1) then 
  write(*,*) "resetting f"
  do el=1,  maxval(atomlist(structures,:)%element_index)
     do i=0, sum(1)
        do j=0, sum(2)
           do k=0, sum(3)
              f_body_matrix(i+1,j+1,k+1,atomlist(structures,atom_number_previous+1)%element_index)=1.0
              product_matrix(i+1,j+1,k+1,el)=t_body_matrix(i+1,j+1,k+1,el)*&
                   &f_body_matrix(i+1,j+1,k+1,el)*th_body_matrix(i+1,j+1,k+1,el)
              !          write(*,*) product_matrix(i+1,j+1,k+1,el)
              if(summation.lt.product_matrix(i+1,j+1,k+1,el)) summation=product_matrix(i+1,j+1,k+1,el)
           end do
        end do
     end do
  end do

end if
end if
placed=.FALSE.
write(*,*) "WAS ATOM PLACED?", placed
call execute_command_line("rm buildmap_testfile.txt",WAIT=.TRUE.)   
close(6969)
return
else
placed=.TRUE.
end if


open(1050,file="heatmap_plotter")

do el=1,  maxval(atomlist(structures,:)%element_index)
write(name,'(A,I0,A,I0,A,I0)') "pos/POSCAR_", structurecounter("pos"),"heatmap_plotter",el,"_",atom_number_previous+1
open(1060+el,file=trim(adjustl(name)))
write(1060+el,*) "x,y,z,data,label,color"
do i=0, sum(1)
do j=0, sum(2)
  do k=0,sum(3)
     tmpvector(1)=i*dble(formula(structures)%cell(1,1)/sum(1))
     tmpvector(2)=i*dble(formula(structures)%cell(1,2)/sum(1))
     tmpvector(3)=i*dble(formula(structures)%cell(1,3)/sum(1))
     tmpvector(1)=tmpvector(1)+j*dble(formula(structures)%cell(2,1)/sum(2))
     tmpvector(2)=tmpvector(2)+j*dble(formula(structures)%cell(2,2)/sum(2))
     tmpvector(3)=tmpvector(3)+j*dble(formula(structures)%cell(2,3)/sum(2))
     tmpvector(1)=tmpvector(1)+k*dble(formula(structures)%cell(3,1)/sum(3))
     tmpvector(2)=tmpvector(2)+k*dble(formula(structures)%cell(3,2)/sum(3))
     tmpvector(3)=tmpvector(3)+k*dble(formula(structures)%cell(3,3)/sum(3))


     product_matrix(i+1,j+1,k+1,el)=f_body_matrix(i+1,j+1,k+1,el)*th_body_matrix(i+1,j+1,k+1,el)&
          &*t_body_matrix(i+1,j+1,k+1,el)


     if((tmpvector(1)-best_vector(1)).lt.0.01) then 
        if((tmpvector(2)-best_vector(2)).lt.0.01) then
           if((tmpvector(3)-best_vector(3)).lt.0.01) then
              do L=1, 3
                 tmpvector(L)=best(L)
              end do
              product_matrix(i+1,j+1,k+1,el)=calculated_value
           end if
        end if
     end if

     !        product_matrix(i+1,j+1,k2*(product_matrix(i+1,j+1,k+1,el)**2/summation**2)
     ! if(isnan(2*(product_matrix(i+1,j+1,k+1,el)**2/summation**2))) then 
     !    !write(*,*) product_matrix(i+1,j+1,k+1,el), product_matrix(i+1,j+1,k+1,el)**2, summation, summation**2

     !    if(product_matrix(i+1,j+1,k+1,el)**2.eq.0) then 
     !       product_matrix(i+1,j+1,k+1,el)=(product_matrix(i+1,j+1,k+1,el))
     !    else
     !       product_matrix(i+1,j+1,k+1,el)=(product_matrix(i+1,j+1,k+1,el))**2
     !    end if
     !    if(summation**2.eq.0) then 
     !       summation=summation
     !    else 
     !       summation=summation**2
     !    end if
     ! end if
     !product_matrix(i+1,j+1,k+1,el)=(product_matrix(i+1,j+1,k+1,el)/summation)

     write(6969,*) tmpvector(1),tmpvector(2),tmpvector(3), product_matrix(i+1,j+1,k+1,el)&
          &,f_body_matrix(i+1,j+1,k+1,el), t_body_matrix(i+1,j+1,k+1,el), &
          &th_body_matrix(i+1,j+1,k+1,el)
     if(atomlist(structures,atom_number_previous+1)%element_index.ne.el) cycle

     write(1050,*) tmpvector(1),tmpvector(2),tmpvector(3), product_matrix(i+1,j+1,k+1,el)&
          &,f_body_matrix(i+1,j+1,k+1,el), th_body_matrix(i+1,j+1,k+1,el), &
          &t_body_matrix(i+1,j+1,k+1,el)
     if(el.eq.1) then 
        write(1060+el,*) tmpvector(1),", ", tmpvector(2),", ", tmpvector(3),", ",&
             &product_matrix(i+1,j+1,k+1,el),", ",elnames(el)," blue"
     else
        write(1060+el,*) tmpvector(1),", ", tmpvector(2),", ", tmpvector(3),", ",&
             &product_matrix(i+1,j+1,k+1,el),", ",elnames(el),", green"
     end if

     !write(*,*) product_matrix(i+1,j+1,k+1,el)
     results_matrix(i+1,j+1,k+1,1,el)=tmpvector(1) 
     results_matrix(i+1,j+1,k+1,2,el)=tmpvector(2)
     results_matrix(i+1,j+1,k+1,3,el)=tmpvector(3)
     results_matrix(i+1,j+1,k+1,4,el)=product_matrix(i+1,j+1,k+1,el)
  end do
end do
end do
close(1060+el)
end do
do i=1, atom_number_previous+1
write(1000,*) atomlist(structures,i)%position
!write(6969,*) atomlist(structures,i)%position
end do
close(1000)
close(6969)






write(*,*) "BUILDMAP COMPLETED"
close(1061) 





end subroutine buildmap_WIP


end module buildmap