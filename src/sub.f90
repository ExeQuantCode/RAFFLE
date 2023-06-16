module help
use atomtype
implicit none


TYPE densitymatrix 
   double precision, dimension(3) :: position 
   double precision :: density 
   integer :: checked
end type densitymatrix

TYPE unitcell 
   double precision, dimension(3,3) :: cell 
end type unitcell

TYPE nearmatrix
   double precision, dimension(3) :: position 
end type nearmatrix

TYPE sim_check
   character(3), dimension(:), allocatable :: names
   integer, dimension(:), allocatable :: stoichio   
end type sim_check

contains


subroutine Duplicate_check ()
  character(1024) :: command, name, name2, buffer 
  integer :: i,j,k,l,x,q,w,e,r, loop_structures, loop_elements, assigned_index, error_calculator,nbin, stat 
  TYPE(sim_check), dimension(:), allocatable :: database
  integer, dimension(:), allocatable :: total_atoms
  double precision, dimension(2) :: Maximum_value 
  double precision, dimension(:), allocatable :: read_in
  double precision, dimension(:,:), allocatable :: Gaussian, tmp, comp_a, comp_b
  double precision :: Gaus_in, bondcut, sigma, similarity_index
  call execute_command_line("./similarity.sh")
  
!!!Inputs 
bondcut=10  
sigma=0.1
!!!
  error_calculator=0
  open(11,file="similarity_combinations.txt")
  read(11,*) loop_structures
  
  print*, loop_structures
  allocate(database(loop_structures))
  allocate(total_atoms(loop_structures))
  do i=1, loop_structures-error_calculator
     
     read(11,*, END=99) assigned_index
     
     print*, i, assigned_index, loop_structures-error_calculator
    
     error_calculator=assigned_index-i

     read(11,*) loop_elements 
     allocate(database(i)%names(loop_elements))
     allocate(database(i)%stoichio(loop_elements))
     do j=1, loop_elements
        read(11,*) database(i)%stoichio(j)
        read(11,*) database(i)%names(j)
     end do
     read(11,*) total_atoms(i)
  end do
  99 close(11)
  assigned_index=loop_structures
  error_calculator=assigned_index-i+1
  open(12,file="test")
  print*, loop_structures-error_calculator
  do i=1, loop_structures-error_calculator
     do j=1, size(database(i)%stoichio,1)
        !print*, i,j
        !print*, database(i)%names(j), database(i)%stoichio(j)
     end do
  end do
  close(12)
  

  


  nbin=1000
  allocate(Gaussian(nbin,2))
  allocate(comp_a(nbin,2))
  allocate(comp_b(nbin,2))
 allocate(Tmp(nbin,2))
  
  do i=1, nbin 
     Gaussian(i,1)=dble(i*bondcut/nbin)
     Gaussian(i,2)=0
  end do
  print*, "HERE"
  do i=1, loop_structures-error_calculator
     do j=1, size(database(i)%stoichio,1)
        do k=1, database(i)%stoichio(j)
           do x=1, size(database(i)%stoichio,1)
              write(name,'(A,I0.3,A,A,A,I0.3,A,A,A,A,A)')"bon/BON_", i,"/",trim(adjustl(database(i)%names(j))),&
                   &"_",k,trim(adjustl(database(i)%names(j))),"_",trim(adjustl(database(i)%names(x))),"_bond_distribution"
              open(12,file=trim(adjustl(name)))
              
              do while(1.eq.1) 
                 read(12,*,IOSTAT=stat) Gaus_in
                 !print*, Gaus_in
                 IF(IS_IOSTAT_END(stat)) exit 
                 do l=1, nbin
                    Gaussian(l,2)=Gaussian(l,2)+Gaus_in*(1/(2*3.14159)**2)*exp(-0.5*((Gaussian(l,1)-Gaus_in)/sigma)**2)&
                         &/(Gaussian(L,1))**2
                 end do
              end do
              close(12)
              maximum_value=maxval(Gaussian,1)

              write(name,'(A,I0.3,A,A,A,I0.3,A,A,A,A,A)')"bon/BON_", i,"/",trim(adjustl(database(i)%names(j))),&
                   &"_",k,trim(adjustl(database(i)%names(j))),"_",trim(adjustl(database(i)%names(x))),"_gaussian"
              open(12,file=trim(adjustl(name)))
              print*, trim(adjustl(name))

              do L=1, nbin
                 Gaussian(L,2)=Gaussian(L,2)/maximum_value(2)
                 write(12,*) Gaussian(L,1), Gaussian(L,2)
              end do
              close(12)
           end do
        end do
     end do
  end do


  do i=1, loop_structures-error_calculator
     do j=1, size(database(i)%stoichio,1)
        do k=1, database(i)%stoichio(j)
           do x=1, size(database(i)%stoichio,1)

!!!! Find a bon profile of ith structure, kth atom of Jths species in relation to xth species.  
              write(name,'(A,I0.3,A,A,A,I0.3,A,A,A,A,A)')"bon/BON_", i,"/",&
                   &trim(adjustl(database(i)%names(j))),&
                   &"_",k,trim(adjustl(database(i)%names(j))),"_",&
                   &trim(adjustl(database(i)%names(x))),"_gaussian"
              open(12,file=trim(adjustl(name)))
              !print*, name 
              do L=1, nbin
                 read(12,*,IOSTAT=stat) comp_a(L,1), comp_a(L,2)
              end do
              if(i+1.gt.loop_structures-error_calculator) then 
                 close(12) 
                 cycle 
              end if
              do q=i+1, loop_structures-error_calculator
                 do w=1, size(database(q)%stoichio,1)
                    do e=1, database(q)%stoichio(w)
                       do r=1, size(database(q)%stoichio,1)

!!! As above, find something to compare it to. Jth and Xth species types should match with Wth and Rth, 
!!! but could be swapped around.  
                           if((database(i)%names(j).ne.database(q)%names(w)).OR.&
                                &(database(i)%names(j).ne.database(q)%names(r))) cycle
                           if((database(i)%names(x).ne.database(q)%names(w)).OR.&
                                (database(i)%names(x).ne.database(q)%names(r))) cycle

                           write(name,'(A,I0.3,A,A,A,I0.3,A,A,A,A,A)')"bon/BON_", q,"/",&
                                &trim(adjustl(database(q)%names(w))),&
                                &"_",e,trim(adjustl(database(q)%names(w))),"_",trim(adjustl(database(q)%names(r))),&
                                &"_gaussian"
                           open(12,file=trim(adjustl(name)))

                           !print*, name
                          do L=1, nbin
                             read(12,*,IOSTAT=stat) comp_b(L,1), comp_b(L,2)
                          end do
                          similarity_index=0
                          do L=1, nbin 
                             
                             if(L.ne.nbin) then 
                                similarity_index=similarity_index+abs(comp_a(L,2)**2-comp_b(L,2)**2)*(comp_a(L+1,1)-comp_a(L,1))
                             else
                                similarity_index=similarity_index+abs(comp_a(L,2)**2-comp_b(L,2)**2)*(comp_a(L,1)-comp_a(L-1,1))
                             end if
                          end do
                          print*, similarity_index, i,j,k,x,q,w,e,r
                          if(similarity_index.lt.0.08) then 
                              open(10,file="TEST") 
                              
                             write(10,*) i,j,k,q,w,e
                             
                           end if
                          
                          

                       end do
                    end do
                 end do
              end do
           end do
        end do
     end do
  end do






end subroutine Duplicate_check



! subroutine buildmap (bin_size,formula,atomlist,alistrep&
!      &,atom_number_previous,structures,elrad,atom_total,results_matrix) 
! type(unitcell), dimension(:), allocatable :: formula
! type (atom), dimension(:,:), allocatable :: atomlist,alistrep,predicted_positions
! integer :: i,j,k,l,m,n, normalisation,atom_number_previous,structures, bin_size, sum, m_count, n_count, atom_total
! double precision :: value_return, summation
! double precision, dimension(3) :: tmpvector, location_vector 
! ! double precision, dimension(:,:,:), allocatable :: elrad
! ! double precision, dimension(:,:,:), allocatable :: update_region,f_body_matrix, th_body_matrix, t_body_matrix, product_matrix
! ! double precision, dimension(:,:), allocatable :: position_storage
! ! double precision, dimension(:,:,:,:), allocatable ::  results_matrix
! ! logical :: file_exists
! ! character(5) :: APP

! ! APP="APPEND"
! ! print*, "WELCOME TO THE NEW BUILDMAP FUNCTION; EXPERIMENTAL"
! ! sum=bin_size
! ! allocate(update_region(sum+1,sum+1,sum+1))
! ! allocate(predicted_positions(1,26))
! ! allocate(f_body_matrix(sum+1,sum+1,sum+1))
! ! !!1 is value, 2 is contributions to value
! ! allocate(t_body_matrix(sum+1,sum+1,sum+1))
! ! allocate(th_body_matrix(sum+1,sum+1,sum+1))
! ! allocate(position_storage(3,3))
! ! allocate(product_matrix(sum+1,sum+1,sum+1))
! summation=0

! f_body_matrix=1
! t_body_matrix=0
! th_body_matrix=1
! update_region=1
! product_matrix
! INQUIRE(FILE="buildmap_testfile.txt", EXIST=file_exists)
! if(file_exists) then 
!    update_region=0
!    print*, "Loading in an existing buildmap. PLEASE ADD SPECIES SUPPORT"
!    open(6969,file="buildmap_testfile.txt")
   
!    do i=0, sum 
!       do j=0, sum
!          do k=0, sum 
!             read(6969,*) location_vector,&
!                  &product_matrix(i+1,j+1,k+1),&
!                  &f_body_matrix(i+1,j+1,k+1),&
!                  &t_body_matrix(i+1,j+1,k+1),&
!                  &th_body_Matrix(i+1,j+1,k+1)
!             logic_loop : do L=(atom_number_previous-1)*27+1,(atom_number_previous)*27
!                if(bondlength(location_vector,alistrep(structures,L)%position)&
!                     &.lt.3.2) then 
!                   update_region(i+1,j+1,k+1)=1 
!                   exit logic_loop
!                end if
!             end do logic_loop
!          end do
!       end do
!    end do
!    rewind(6969)
! else
   
   
!    open(6969,file="buildmap_testfile.txt")
! end if

! product_matrix=0

! print*, "Arriving here"

! do i=0, sum
   
!    do j=0, sum
!       firstloop :do k=0,sum  
!          print*, i,j,k
!          if(f_body_matrix(i+1,j+1,k+1).eq.0) cycle firstloop
!          if(th_body_matrix(i+1,j+1,k+1).eq.0) cycle firstloop
!          if(update_region(i+1,j+1,k+1).lt.0.5) then 
!             cycle firstloop
!          end if
!          tmpvector(1)=i*dble(formula(structures)%cell(1,1)/sum)
!          tmpvector(2)=i*dble(formula(structures)%cell(1,2)/sum)
!          tmpvector(3)=i*dble(formula(structures)%cell(1,3)/sum)
!          tmpvector(1)=tmpvector(1)+j*dble(formula(structures)%cell(2,1)/sum)
!          tmpvector(2)=tmpvector(2)+j*dble(formula(structures)%cell(2,2)/sum)
!          tmpvector(3)=tmpvector(3)+j*dble(formula(structures)%cell(2,3)/sum)
!          tmpvector(1)=tmpvector(1)+k*dble(formula(structures)%cell(3,1)/sum)
!          tmpvector(2)=tmpvector(2)+k*dble(formula(structures)%cell(3,2)/sum)
!          tmpvector(3)=tmpvector(3)+k*dble(formula(structures)%cell(3,3)/sum)
!          call atomprojector(tmpvector,predicted_positions,formula,atom_number_previous+1,structures)

!          m_count=0
!          n_count=0
!          do L=1, atom_number_previous*27+26
!             !cutoff here should be made sensibly for how close exclusion range is to each atom
!             if(L.gt.atom_number_previous*27) then 
!                position_storage(1,:)=predicted_positions(1,L-atom_number_previous*27)%position(:)
!             else 
!                position_storage(1,:)=alistrep(structures,L)%position(:)
!             end if
!             value_return=0
!             if(bondlength(tmpvector,&
!                  &position_storage(1,:)).lt.1.2) then 
!                t_body_matrix(i+1,j+1,k+1)=0
!                f_body_matrix(i+1,j+1,k+1)=0
!                th_body_matrix(i+1,j+1,k+1)=0
!                cycle firstloop
!             else if(bondlength(tmpvector,&
!                  &position_storage(1,:)).gt.5) then
!                !call evaluate_contribution (trim(adjustl(&
!                !     &atomlist(structures,atom_number_previous+1)%name)),&
!                !     &trim(adjustl(alistrep(structures,atom_number_previous+1)%name)),bondlength(tmpvector,&
!                !     &position_storage(1,:)),value_return)
!                cycle
!             else
!                call evaluate_contribution (trim(adjustl(&
!                     &atomlist(structures,atom_number_previous+1)%name)),&
!                     &trim(adjustl(alistrep(structures,atom_number_previous+1)%name)),bondlength(tmpvector,&
!                     &position_storage(1,:)),value_return)
!             end if

!             t_body_matrix(i+1,j+1,k+1)=t_body_matrix(i+1,j+1,k+1)+value_return
!             if(bondlength(&
!                  &tmpvector,&
!                  &position_storage(1,:)).lt.2) then
              
!                do n=1, atom_number_previous*27+26
!                   if(n.gt.atom_number_previous*27) then
!                      position_storage(2,:)=predicted_positions(1,n-atom_number_previous*27)%position(:)
!                   else
!                      position_storage(2,:)=alistrep(structures,n)%position(:)
!                   end if
!                   if(n.eq.L) cycle
!                   if(bondlength(&
!                        &tmpvector,&
!                        &position_storage(2,:)).lt.1.2) then
!                      t_body_matrix(i+1,j+1,k+1)=0
!                      f_body_matrix(i+1,j+1,k+1)=0
!                      th_body_matrix(i+1,j+1,k+1)=0
!                      cycle firstloop
!                   end if
                  
!                   if(bondlength(&
!                        &position_storage(1,:),&
!                        &position_storage(2,:)).lt.2) then
!                      n_count=n_count+1
!                      call evaluate_angle_contribution(trim(adjustl(&
!                           &atomlist(structures,atom_number_previous+1)%name))&
!                           &,bondangle(&
!                           &tmpvector,&
!                           &position_storage(1,:),&
!                           &position_storage(2,:)),&
!                           &value_return)
                     
!                      th_body_matrix(i+1,j+1,k+1)=th_body_matrix(i+1,j+1,k+1)*value_return
!                   else if(bondlength(&
!                        &tmpvector,&
!                        &position_storage(2,:)).lt.2) then 
!                      n_count=n_count+1
!                      call evaluate_angle_contribution(trim(adjustl(&
!                           &atomlist(structures,atom_number_previous+1)%name))&
!                           &,bondangle(&
!                           &position_storage(1,:),&
!                           &tmpvector,&
!                           &position_storage(2,:)),&
!                           &value_return)

!                      th_body_matrix(i+1,j+1,k+1)=th_body_matrix(i+1,j+1,k+1)*value_return

                        
!                   end if
                  
!                   if(bondlength(&
!                        &position_storage(1,:),& 
!                        &position_storage(2,:)).lt.2) then
                     
!                      do m=1, atom_number_previous*27+26 
!                         if(m.eq.L) cycle 
!                         if(m.eq.n) cycle
!                         if(m.gt.atom_number_previous*27) then 
!                            position_storage(3,:)=predicted_positions(1,m-atom_number_previous*27)%position(:)
!                         else 
!                            position_storage(3,:)=alistrep(structures,m)%position(:)
!                         end if
!                         if(bondlength(&
!                              &tmpvector,&
!                              &position_storage(3,:)).lt.1.2) then
!                            f_body_matrix(i+1,j+1,k+1)=0
!                            th_body_matrix(i+1,j+1,k+1)=0
!                            t_body_matrix(i+1,j+1,k+1)=0
!                            cycle firstloop
!                         end if
!                         if(bondlength(&
!                              &position_storage(1,:),&
!                              &position_storage(3,:)).lt.2) then
!                            m_count=m_count+1
!                            call evaluate_4body_contribution(trim(adjustl(&
!                                 &atomlist(structures,atom_number_previous+1)%name))&
!                                 &,fourbody(&
!                                 &position_storage(1,:),&
!                                 &position_storage(2,:),&
!                                 &position_storage(3,:),&
!                                 &tmpvector&
!                                 ),&
!                                 &value_return)
!                            f_body_matrix(i+1,j+1,k+1)=f_body_matrix(i+1,j+1,k+1)*value_return
!                            if(value_return.eq.0) then 
!                               cycle firstloop
!                               !else 
!                               print*,f_body_matrix(i+1,j+1,k+1), value_return, i,j,k
!                            end if
!                         end if
!                      end do
!                   end if
!                end do
!             end if
!          end do
!          if(m_count.eq.0) then 
!             f_body_matrix(i+1,j+1,k+1)=1
!             !print*, "4_body set to zero on", i, j, k 
!          end if
!          if(n_count.eq.0) then 
            
!             th_body_matrix(i+1,j+1,k+1)=1
!          end if
!          product_matrix(i+1,j+1,k+1)=t_body_matrix(i+1,j+1,k+1)*f_body_matrix(i+1,j+1,k+1)*th_body_matrix(i+1,j+1,k+1)
!          if(product_matrix(i+1,j+1,k+1).gt.summation) summation=product_matrix(i+1,j+1,k+1)
!       end do firstloop
!    end do
! end do

! open(1000,file="heatmap_plotter.txt",ACCESS=APP)
! results_matrix=0
! do i=0, sum
!    do j=0, sum
!       do k=0,sum
!          tmpvector(1)=i*dble(formula(structures)%cell(1,1)/sum)
!          tmpvector(2)=i*dble(formula(structures)%cell(1,2)/sum)
!          tmpvector(3)=i*dble(formula(structures)%cell(1,3)/sum)
!          tmpvector(1)=tmpvector(1)+j*dble(formula(structures)%cell(2,1)/sum)
!          tmpvector(2)=tmpvector(2)+j*dble(formula(structures)%cell(2,2)/sum)
!          tmpvector(3)=tmpvector(3)+j*dble(formula(structures)%cell(2,3)/sum)
!          tmpvector(1)=tmpvector(1)+k*dble(formula(structures)%cell(3,1)/sum)
!          tmpvector(2)=tmpvector(2)+k*dble(formula(structures)%cell(3,2)/sum)
!          tmpvector(3)=tmpvector(3)+k*dble(formula(structures)%cell(3,3)/sum)
         
!          product_matrix(i+1,j+1,k+1)=2*(product_matrix(i+1,j+1,k+1)**2/summation**2)
!          write(6969,*) tmpvector(1),tmpvector(2),tmpvector(3), product_matrix(i+1,j+1,k+1)&
!               &,f_body_matrix(i+1,j+1,k+1), t_body_matrix(i+1,j+1,k+1), &
!               &th_body_matrix(i+1,j+1,k+1)
!          write(1000,*) tmpvector(1),tmpvector(2),tmpvector(3), product_matrix(i+1,j+1,k+1)&
!               &,f_body_matrix(i+1,j+1,k+1), t_body_matrix(i+1,j+1,k+1), &
!               &th_body_matrix(i+1,j+1,k+1)

!          results_matrix(i+1,j+1,k+1,1)=tmpvector(1) 
!          results_matrix(i+1,j+1,k+1,2)=tmpvector(2)
!          results_matrix(i+1,j+1,k+1,3)=tmpvector(3)
!          results_matrix(i+1,j+1,k+1,4)=product_matrix(i+1,j+1,k+1)
!       end do
!    end do
! end do
! do i=1, atom_number_previous+1
!    write(1000,*) atomlist(structures,i)%position

!    write(6969,*) atomlist(structures,i)%position
! end do
! close(1000)
! close(6969)

! print*, "BUILDMAP COMPLETED"
! end subroutine buildmap

subroutine buildmap_POINT(tmpvector,formula,atomlist,alistrep&
     &,atom_number_previous,structures,elrad,atom_total,eltot, elnames,placed,num_VOID&
     &,uptol,lowtol,calculated_value)

type(unitcell), dimension(:), allocatable :: formula
type (atom), dimension(:,:), allocatable :: atomlist,alistrep,predicted_positions
integer :: el,i,j,k,l,m,n, normalisation,atom_number_previous,structures,&
     &m_count, n_count, atom_total, eltot, cleanup, cleanup2, repeats, num_VOID, c_cut, norm&
     &,c_min
double precision :: value_return, summation, max_1, max_2, max_3, comparison, uptol, lowtol, normaliser&
     &,repeat_power, i_comp, j_comp, k_comp, totbin, calculated_value
integer, dimension(3) :: sum, bin_size
double precision, dimension(3) :: tmpvector, location_vector
double precision, dimension(:,:,:), allocatable :: elrad
double precision, dimension(:), allocatable :: f_body_matrix, th_body_matrix, t_body_matrix, product_matrix
double precision, dimension(:,:), allocatable :: position_storage
integer, dimension(:,:), allocatable :: index_storage
character(3), dimension(:,:), allocatable :: name_storage, remaining_elements, remaining_2
double precision, dimension(:,:,:,:,:), allocatable ::  results_matrix
logical :: file_exists, placed
character(3), dimension(:), allocatable :: elnames
integer, dimension(:), allocatable :: tmpdig
character(1024), dimension(:), allocatable :: tmpels
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

!print*, tmpvector, atomlist(structures,atom_number_previous+1)%name
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
      
      if(bondlength(tmpvector,&
           &position_storage(1,:)).lt.&
           &elrad(1,atomlist(structures,atom_number_previous+1)%element_index,index_storage(1,1))*lowtol) then 
         !print*, tmpvector, bondlength(tmpvector,position_storage(1,:))
         t_body_matrix(el)=0
         f_body_matrix(el)=0
         th_body_matrix(el)=0

         !if((i.eq.2).and.(j.eq.2).and.(k.eq.5)) then
         !   print*, "cycle: bondlength"
         !end if

         cycle elloop
      else if(bondlength(tmpvector,&
           &position_storage(1,:)).gt.&
           &elrad(1,atomlist(structures,atom_number_previous+1)%element_index,index_storage(1,1))*uptol) then
!!!! This has been left out to ease on computation. Could be reimplemented but increases cost dramatically to consider ALL atoms
         !call evaluate_contribution (trim(adjustl(&
         !     &atomlist(structures,atom_number_previous+1)%name)),&
         !     &trim(adjustl(alistrep(structures,atom_number_previous+1)%name)),bondlength(tmpvector,&
         !     &position_storage(1,:)),value_return)
         !if(bondlength(tmpvector,position_storage(1,:)).lt.4) then 
         !   print*, tmpvector, bondlength(tmpvector,position_storage(1,:)), "!"
         !end if
         !print*, "CYCLE0.5", elrad(1,atomlist(structures,atom_number_previous+1)%element_index,index_storage(1,1))*uptol&
         !     &,bondlength(tmpvector,&
         !     &position_storage(1,:))
         cycle
      else

         call evaluate_contribution (trim(adjustl(&
              &atomlist(structures,atom_number_previous+1)%name)),&
              &trim(adjustl(name_storage(1,1))),bondlength(tmpvector,&
              &position_storage(1,:)),value_return)

         t_body_matrix(el)=((t_body_matrix(el)*norm)+value_return)
         !ADD +1 to norm here to curtail contributions to bondlength. I do not think you need to, as it can dampen certain resonance points 
         norm=norm+1
         t_body_matrix(el)=t_body_matrix(el)/norm
         !print*, tmpvector, t_body_matrix(el), el
         !if((i.eq.2).and.(j.eq.2).and.(k.eq.5)) then
         !   print*, value_return, "2b"
         !end if


      end if
      if(bondlength(&
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
            if(bondlength(&
                 &tmpvector,&
                 &position_storage(2,:)).lt.&
                 &elrad(1,atomlist(structures,atom_number_previous+1)%element_index,index_storage(2,1))*lowtol) then
               !if((i.eq.2).and.(j.eq.2).and.(k.eq.5)) then
               !   print*, "cycle: bondlength ang"
               !end if

               t_body_matrix(el)=0
               f_body_matrix(el)=0
               th_body_matrix(el)=0
               !print*, "CYCLE 2"
               cycle elloop
            end if
            if(bondlength(&
                 &position_storage(1,:),&
                 &position_storage(2,:)).lt.&
                 &elrad(1,index_storage(1,1),index_storage(2,1))*uptol) then

               call evaluate_angle_contribution(trim(adjustl(&
                    &atomlist(structures,atom_number_previous+1)%name))&
                    &,bondangle(&
                    &tmpvector,&
                    &position_storage(1,:),&
                    &position_storage(2,:)),&
                    &value_return)

               th_body_matrix(el)=(th_body_matrix(el)*value_return**&
                    &(1.0/(repeat_power)))!*n_count)*value_return
               n_count=n_count+1
               th_body_matrix(el)=th_body_matrix(el)!/n_count
               !if((i.eq.2).and.(j.eq.2).and.(k.eq.5)) then
               !   print*, value_return, "ang"
               !end if

            else if(bondlength(&
                 &tmpvector,&
                 &position_storage(2,:)).lt.&
                 &elrad(1,atomlist(structures,atom_number_previous+1)%element_index,index_storage(2,1))*uptol) then 

               call evaluate_angle_contribution(trim(adjustl(&
                    &atomlist(structures,atom_number_previous+1)%name))&
                    &,bondangle(&
                    &position_storage(1,:),&
                    &tmpvector,&
                    &position_storage(2,:)),&
                    &value_return)

               th_body_matrix(el)=(th_body_matrix(el)*(value_return**&
                    (1.0/(repeat_power))))!n_count)+value_return
               n_count=n_count+1
               th_body_matrix(el)=th_body_matrix(el)!/n_count
               !if((i.eq.2).and.(j.eq.2).and.(k.eq.5)) then
               !   print*, value_return, "ang"
               !end if

            end if
            !if((i.eq.2).and.(j.eq.2).and.(k.eq.5)) then
            !   print*, bondlength(&
            !        &position_storage(1,:),&
            !        &position_storage(2,:)), "!"
            !
            !end if

            if((bondlength(&
                 &position_storage(1,:),& 
                 &position_storage(2,:)).lt.&
                 &elrad(1,index_storage(1,1),index_storage(2,1))*uptol)&
                 &.OR.(bondlength(&
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
                  if(bondlength(&
                       &tmpvector,&
                       &position_storage(3,:)).lt.&
                       elrad(1,atomlist(structures,atom_number_previous+1)%element_index,index_storage(3,1))*lowtol) then

                     f_body_matrix(el)=0
                     th_body_matrix(el)=0
                     t_body_matrix(el)=0

                     !if((i.eq.2).and.(j.eq.2).and.(k.eq.5)) then
                        !print*, tmpvector, position_storage(3,:), bondlength(&
                        !     &tmpvector,&
                         !    &position_storage(3,:))
                     !end if
                     !print*, "4 cycle"
                     cycle elloop
                  end if
                  if(bondlength(&
                       &position_storage(1,:),&
                       &position_storage(3,:)).lt.&
                       elrad(1,index_storage(1,1),index_storage(3,1))*uptol) then
                     call evaluate_4body_contribution(trim(adjustl(&
                          &atomlist(structures,atom_number_previous+1)%name))&
                          &,fourbody(&
                          &tmpvector,&
                          &position_storage(1,:),&
                          &position_storage(2,:),&
                          &position_storage(3,:)&
                          &),&
                          &value_return)
                     !if((i.eq.2).and.(j.eq.2).and.(k.eq.5)) then
                        !print*, f_body_matrix(el), value_return, "4b", index_storage(:,1)
                     !end if
                     !print*, "!!!!!!!!!!!!!!!!!!!"
                     !print*, value_return
                     !print*, tmpvector
                     !print*, position_storage(1,:)
                     !print*, position_storage(2,:)
                     !print*, position_storage(3,:)
                     !print*, "!!!!!!!!!!!!!!!!!!!"     
                     
                     if(value_return.eq.0) then

                        f_body_matrix(el)=0
                        th_body_matrix(el)=0
                        t_body_matrix(el)=0

                        !print*, "CYCLE 4"
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
                     !print*, f_body_matrix(i+1,j+1,k+1,el), value_return**(1.0/4.0),value_return,i,j,k
                  end if
               end do
            end if
         end do
      end if
   end do
   !print*, f_body_matrix(i+1,j+1,k+1,el)
   if(m_count.eq.0) then 
      f_body_matrix(el)=1
      !print*, "4_body set to zero on", tmpvector 
   end if
   if(n_count.eq.0) then 
      !print*, "3_body set to one on", tmpvector
      th_body_matrix(el)=1
   end if
   product_matrix(el)=t_body_matrix(el)*&
        &f_body_matrix(el)*th_body_matrix(el)
   calculated_value=product_matrix(el)
   !print*, product_matrix(i+1,j+1,k+1,el), f_body_matrix(i+1,j+1,k+1,el), i, j, k
   summation=product_matrix(el)
   !print*, tmpvector, product_matrix(el),f_body_matrix(el), "!!!" 
   

end do elloop


end subroutine buildmap_POINT


subroutine buildmap_WIP (bin_size,formula,atomlist,alistrep&
     &,atom_number_previous,structures,elrad,atom_total,&
     &results_matrix,eltot, elnames,placed,num_VOID, append_matrix) 

type(unitcell), dimension(:), allocatable :: formula
type (atom), dimension(:,:), allocatable :: atomlist,alistrep,predicted_positions
integer :: el,i,j,k,l,m,n, normalisation,atom_number_previous,structures,&
     &m_count, n_count, atom_total, eltot, cleanup, cleanup2, repeats, num_VOID, c_cut, norm&
     &,c_min
double precision :: value_return, summation, max_1, max_2, max_3, comparison, uptol, lowtol, normaliser&
     &,repeat_power, i_comp, j_comp, k_comp, totbin, calculated_value
integer, dimension(3) :: sum, bin_size, best
double precision, dimension(3) :: tmpvector, location_vector, best_vector
double precision, dimension(:,:,:), allocatable :: elrad
double precision, dimension(:,:,:,:), allocatable :: update_region,f_body_matrix, th_body_matrix, t_body_matrix, product_matrix
double precision, dimension(:,:), allocatable :: position_storage
integer, dimension(:,:), allocatable :: index_storage 
character(3), dimension(:,:), allocatable :: name_storage, remaining_elements, remaining_2 
double precision, dimension(:,:,:,:,:), allocatable ::  results_matrix, append_matrix
logical :: file_exists, placed
character(3), dimension(:), allocatable :: elnames
integer, dimension(:), allocatable :: tmpdig
character(1024), dimension(:), allocatable :: tmpels
character(5) :: APP
character(1024) :: name

uptol=1.1
lowtol=0.95



repeat_power=1


APP="APPEND"
print*, "WELCOME TO THE NEW BUILDMAP FUNCTION; EXPERIMENTAL"
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
!print*,  maxval(atomlist(structures,:)%element_index)

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
   print*, "Loading in an existing buildmap. PLEASE ADD SPECIES SUPPORT"
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
            print*, "CLEANUP CYCLE"
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
                  !print*, "INDEX CYCLING"
                  cycle
               end if

               !Chops off the top of the unit cell, should be migrated to the Infile and extended to 3D
               !if(location_vector(3).gt.11) then 
               !   update_region(i+1,j+1,k+1,m)=0
               !end if


               logic_loop : do L=(atom_number_previous-num_VOID)*27+1,(atom_number_previous)*27
                  if(bondlength(location_vector,alistrep(structures,L)%position)&
                       &.gt.(elrad(3,m,&
                       &alistrep(structures,L)%element_index)*2.0)) then
                     update_region(i+1,j+1,k+1,m)=0
                     !print*, i, j, k
                  end if
               end do logic_loop

               logic_loop2 : do L=(atom_number_previous-num_VOID)*27+1,(atom_number_previous)*27
                  if(bondlength(location_vector,alistrep(structures,L)%position)&
                       &.lt.(elrad(3,m,&
                       &alistrep(structures,L)%element_index)*2.0)) then 
                     !print*, i, j, k , m,1
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
         !print*, i, j, k , update_region(i+1,j+1,k+1,1)
      end do
   end do
end do

print*, "Arriving here", repeats

100 if(repeats.ne.1) then
   print*, "WIPING"
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

call invar(15,tmpdig,tmpels)
c_cut=tmpdig(1)
deallocate(tmpdig)
call invar(17,tmpdig,tmpels)
c_min=tmpdig(1)
deallocate(tmpdig)


!print*, sum 
do i=0, sum(1)
   ! write(6,'(TL10, I3.0, A)', ADVANCE='NO') NINT(100*dble(i*(sum)**2+j*(sum)+k)/dble((sum+1)**3)), "%" !dble((i*sum**2 +j*sum +k)/(sum**3))
   do j=0, sum(2)
      firstloop :do k=0,sum(3)  
         write(6,'(A)',ADVANCE='NO') achar(13)
         i_comp=i*(sum(2)+1.0)*(sum(3)+1.0)
         j_comp=j*(sum(3)+1)
         k_comp=k
         totbin=(sum(1)+1)*(sum(2)+1)*(sum(3)+1)
         write(6,'(I3.0, A)', ADVANCE='NO') NINT(100*(i_comp+j_comp+k_comp)/totbin), "%" 
         comparison=100.0*dble(k)/dble(sum(3))
         !print*, i, j, k
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
         elloop:do el=1, maxval(atomlist(structures,:)%element_index)
            if((f_body_matrix(i+1,j+1,k+1,el).eq.0)) then 
               !print*, "cycling: f body doesn't need updating", i,j,k
               !cycle elloop
               continue
            end if
            if((th_body_matrix(i+1,j+1,k+1,el).eq.0)) then 
               !print*, "CYCLE B", i,j,k
               !cycle elloop
               continue
            end if
            if(update_region(i+1,j+1,k+1,el).ne.1) then
               !print*, "CYCLE C",update_region(i+1,j+1,k+1,el),update_region(i+1,j+1,k+1,1),el, i, j, k
               cycle elloop
            end if
            if(atomlist(structures,atom_number_previous+1)%element_index.ne.el) cycle elloop            

            tmpvector(1)=i*dble(formula(structures)%cell(1,1)/(sum(1)))
            tmpvector(2)=i*dble(formula(structures)%cell(1,2)/(sum(1)))
            tmpvector(3)=i*dble(formula(structures)%cell(1,3)/(sum(1)))
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
               !print*, "-----------------------------------"
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
               if(bondlength(tmpvector,&
                    &position_storage(1,:)).lt.&
                    &elrad(1,atomlist(structures,atom_number_previous+1)%element_index,index_storage(1,1))*lowtol) then 
                  !print*, tmpvector, bondlength(tmpvector,position_storage(1,:))
                  t_body_matrix(i+1,j+1,k+1,el)=0
                  f_body_matrix(i+1,j+1,k+1,el)=0
                  th_body_matrix(i+1,j+1,k+1,el)=0
                  
                  !if((i.eq.2).and.(j.eq.2).and.(k.eq.5)) then
                  !   print*, "cycle: bondlength"
                  !end if
                        
                  !print*, "CYCLE 1", i, j, k
                  cycle elloop
               else if(bondlength(tmpvector,&
                    &position_storage(1,:)).gt.&
                    &elrad(1,atomlist(structures,atom_number_previous+1)%element_index,index_storage(1,1))*uptol) then
!!!! This has been left out to ease on computation. Could be reimplemented but increases cost dramatically to consider ALL atoms
                  !call evaluate_contribution (trim(adjustl(&
                  !     &atomlist(structures,atom_number_previous+1)%name)),&
                  !     &trim(adjustl(alistrep(structures,atom_number_previous+1)%name)),bondlength(tmpvector,&
                  !     &position_storage(1,:)),value_return)
                  !if(bondlength(tmpvector,position_storage(1,:)).lt.4) then 
                  !   print*, tmpvector, bondlength(tmpvector,position_storage(1,:)), "!"
                  !end if
                  !print*, "CYCLE0.5", elrad(1,atomlist(structures,atom_number_previous+1)%element_index,index_storage(1,1))*uptol&
                  !     &,bondlength(tmpvector,&
                  !     &position_storage(1,:))
                  cycle
               else

                  call evaluate_contribution (trim(adjustl(&
                       &atomlist(structures,atom_number_previous+1)%name)),&
                       &trim(adjustl(name_storage(1,1))),bondlength(tmpvector,&
                       &position_storage(1,:)),value_return)

                  t_body_matrix(i+1,j+1,k+1,el)=((t_body_matrix(i+1,j+1,k+1,el)*norm)+value_return)
                  !ADD +1 to norm here to curtail contributions to bondlength. I do not think you need to, as it can dampen certain resonance points 
                  norm=norm+1
                  t_body_matrix(i+1,j+1,k+1,el)=t_body_matrix(i+1,j+1,k+1,el)/norm
                  !print*, tmpvector, t_body_matrix(i+1,j+1,k+1,el),i, j, k
                  !if((i.eq.2).and.(j.eq.2).and.(k.eq.5)) then
                  !   print*, value_return, "2b"
                  !end if


               end if
               if(bondlength(&
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
                     if(bondlength(&
                          &tmpvector,&
                          &position_storage(2,:)).lt.&
                          &elrad(1,atomlist(structures,atom_number_previous+1)%element_index,index_storage(2,1))*lowtol) then
                        !if((i.eq.2).and.(j.eq.2).and.(k.eq.5)) then
                        !   print*, "cycle: bondlength ang"
                        !end if

                        t_body_matrix(i+1,j+1,k+1,el)=0
                        f_body_matrix(i+1,j+1,k+1,el)=0
                        th_body_matrix(i+1,j+1,k+1,el)=0
                        !print*, "CYCLE 2"
                        cycle elloop
                     end if
                     if(bondlength(&
                          &position_storage(1,:),&
                          &position_storage(2,:)).lt.&
                          &elrad(1,index_storage(1,1),index_storage(2,1))*uptol) then
                        
                        call evaluate_angle_contribution(trim(adjustl(&
                             &atomlist(structures,atom_number_previous+1)%name))&
                             &,bondangle(&
                             &tmpvector,&
                             &position_storage(1,:),&
                             &position_storage(2,:)),&
                             &value_return)

                        th_body_matrix(i+1,j+1,k+1,el)=(th_body_matrix(i+1,j+1,k+1,el)*value_return**&
                             &(1.0/(repeat_power)))!*n_count)*value_return
                        n_count=n_count+1
                        th_body_matrix(i+1,j+1,k+1,el)=th_body_matrix(i+1,j+1,k+1,el)!/n_count
                        !if((i.eq.2).and.(j.eq.2).and.(k.eq.5)) then
                        !   print*, value_return, "ang"
                        !end if

                     else if(bondlength(&
                          &tmpvector,&
                          &position_storage(2,:)).lt.&
                          &elrad(1,atomlist(structures,atom_number_previous+1)%element_index,index_storage(2,1))*uptol) then 

                        call evaluate_angle_contribution(trim(adjustl(&
                             &atomlist(structures,atom_number_previous+1)%name))&
                             &,bondangle(&
                             &position_storage(1,:),&
                             &tmpvector,&
                             &position_storage(2,:)),&
                             &value_return)

                        th_body_matrix(i+1,j+1,k+1,el)=(th_body_matrix(i+1,j+1,k+1,el)*(value_return**&
                             (1.0/(repeat_power))))!n_count)+value_return
                        n_count=n_count+1
                        th_body_matrix(i+1,j+1,k+1,el)=th_body_matrix(i+1,j+1,k+1,el)!/n_count
                        !if((i.eq.2).and.(j.eq.2).and.(k.eq.5)) then
                        !   print*, value_return, "ang"
                        !end if
                           
                     end if
                     !if((i.eq.2).and.(j.eq.2).and.(k.eq.5)) then
                     !   print*, bondlength(&
                     !        &position_storage(1,:),&
                     !        &position_storage(2,:)), "!"
                     !
                     !end if
                     
                     if((bondlength(&
                          &position_storage(1,:),& 
                          &position_storage(2,:)).lt.&
                          &elrad(1,index_storage(1,1),index_storage(2,1))*uptol)&
                          &.OR.(bondlength(&
                          &tmpvector,&
                          &position_storage(2,:)).lt.&
                          &elrad(1,index_storage(1,1),index_storage(2,1))*uptol)) then

                        do m=1, atom_number_previous*27+26 
                           !print*, L, m, n
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
                           if(bondlength(&
                                &tmpvector,&
                                &position_storage(3,:)).lt.&
                                elrad(1,atomlist(structures,atom_number_previous+1)%element_index,index_storage(3,1))*lowtol) then

                              f_body_matrix(i+1,j+1,k+1,el)=0
                              th_body_matrix(i+1,j+1,k+1,el)=0
                              t_body_matrix(i+1,j+1,k+1,el)=0

                              !if((i.eq.2).and.(j.eq.2).and.(k.eq.5)) then
                              !   print*, tmpvector, position_storage(3,:), bondlength(&
                              !        &tmpvector,&
                              !        &position_storage(3,:))
                              !end if

                              cycle elloop
                           end if
                           if(bondlength(&
                                &position_storage(1,:),&
                                &position_storage(3,:)).lt.&
                                elrad(1,index_storage(1,1),index_storage(3,1))*uptol) then
                              call evaluate_4body_contribution(trim(adjustl(&
                                   &atomlist(structures,atom_number_previous+1)%name))&
                                   &,fourbody(&
                                   &tmpvector,&
                                   &position_storage(1,:),&
                                   &position_storage(2,:),&
                                   &position_storage(3,:)&
                                   &),&
                                   &value_return)
                              !if((i.eq.2).and.(j.eq.2).and.(k.eq.5)) then
                              !   print*, value_return, "4b"
                              !end if


                              if(m_count.eq.0) then 
                                 normaliser=1
                              end if
                              if(value_return.eq.0) then

                                 f_body_matrix(i+1,j+1,k+1,el)=0
                                 th_body_matrix(i+1,j+1,k+1,el)=0
                                 t_body_matrix(i+1,j+1,k+1,el)=0

                                 !print*, "CYCLE 4"
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
                              !print*, f_body_matrix(i+1,j+1,k+1,el), value_return**(1.0/4.0),value_return,i,j,k
                           end if
                        end do
                     end if
                  end do
               end if
            end do
            !print*, f_body_matrix(i+1,j+1,k+1,el)
            if(m_count.eq.0) then 
               f_body_matrix(i+1,j+1,k+1,el)=1
               !print*, "4_body set to zero on", tmpvector 
            end if
            if(n_count.eq.0) then 
               !print*, "3_body set to one on", tmpvector
               th_body_matrix(i+1,j+1,k+1,el)=1
            end if
            product_matrix(i+1,j+1,k+1,el)=t_body_matrix(i+1,j+1,k+1,el)*&
                 &f_body_matrix(i+1,j+1,k+1,el)*th_body_matrix(i+1,j+1,k+1,el)
            !print*, product_matrix(i+1,j+1,k+1,el), f_body_matrix(i+1,j+1,k+1,el), i, j, k
            if(product_matrix(i+1,j+1,k+1,el).gt.summation) then 
               summation=product_matrix(i+1,j+1,k+1,el)
               best(1)=i
               best(2)=j
               best(3)=k
            end if

         end do elloop
      end do firstloop
   end do
end do
!! 3 is arbitrary maxdensity of grid
l=0
do i=0, sum(1) 
   do j=0, sum(2) 
      do k=0, sum(3) 
         do el=1, maxval(atomlist(structures,:)%element_index)
            l=l+1
            product_matrix(i+1,j+1,k+1,el)=t_body_matrix(i+1,j+1,k+1,el)*&
                 &f_body_matrix(i+1,j+1,k+1,el)*th_body_matrix(i+1,j+1,k+1,el)
            !print*, th_body_matrix(i+1,j+1,k+1,el), t_body_matrix(i+1,j+1,k+1,el), f_body_matrix(i+1,j+1,k+1,el),i,j,k
            
            if(product_matrix(i+1,j+1,k+1,el).gt.summation) then 
               !print*, "DID IT"
               summation=product_matrix(i+1,j+1,k+1,el)
               best(1)=i
               best(2)=j
               best(3)=k
            end if
         end do
      end do
   end do
end do

!print*, product_matrix(3,3,6,1), f_body_matrix(3,3,6,1), th_body_matrix(3,3,6,1)

!print*, summation

best_vector(1)=best(1)*dble(formula(structures)%cell(1,1)/(sum(1)))
best_vector(2)=best(1)*dble(formula(structures)%cell(1,2)/(sum(1)))
best_vector(3)=best(1)*dble(formula(structures)%cell(1,3)/(sum(1)))
best_vector(1)=best_vector(1)+best(2)*dble(formula(structures)%cell(2,1)/(sum(2)))
best_vector(2)=best_vector(2)+best(2)*dble(formula(structures)%cell(2,2)/(sum(2)))
best_vector(3)=best_vector(3)+best(2)*dble(formula(structures)%cell(2,3)/(sum(2)))
best_vector(1)=best_vector(1)+best(3)*dble(formula(structures)%cell(3,1)/(sum(3)))
best_vector(2)=best_vector(2)+best(3)*dble(formula(structures)%cell(3,2)/(sum(3)))
best_vector(3)=best_vector(3)+best(3)*dble(formula(structures)%cell(3,3)/(sum(3)))

!print*, best_vector
!print*, summation
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
   print*, "SUMMATION EFFECTIVELY 0, USE WITH CAUTION"
   summation=0
end if


print*, summation, "SUMMATION"
if(summation.eq.0) then
   if(repeats.lt.10) then
      repeats=repeats+1
      repeat_power=repeat_power+5
      print*, "CYCLING"
      
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
      print*, "BREAKING", max_1, max_2, max_3, max_1*max_2*max_3
      !Set eq.0 to activate
      if(max_2.eq.-1) then 
         print*, "resetting th"
         do el=1,  maxval(atomlist(structures,:)%element_index)
            do i=0, sum(1)
               do j=0, sum(2)
                  do k=0, sum(3)
                     th_body_matrix(i+1,j+1,k+1,atomlist(structures,atom_number_previous+1)%element_index)=1.0
                     product_matrix(i+1,j+1,k+1,el)=t_body_matrix(i+1,j+1,k+1,el)*&
                          &f_body_matrix(i+1,j+1,k+1,el)*th_body_matrix(i+1,j+1,k+1,el)
          !           print*, product_matrix(i+1,j+1,k+1,el)
                     if(summation.lt.product_matrix(i+1,j+1,k+1,el)) summation=product_matrix(i+1,j+1,k+1,el)
                  end do
               end do
            end do
         end do
      end if
      !Set eq.0 to activate
      if(max_3.eq.-1) then 
         print*, "resetting f"
         do el=1,  maxval(atomlist(structures,:)%element_index)
            do i=0, sum(1)
               do j=0, sum(2)
                  do k=0, sum(3)
                     f_body_matrix(i+1,j+1,k+1,atomlist(structures,atom_number_previous+1)%element_index)=1.0
                     product_matrix(i+1,j+1,k+1,el)=t_body_matrix(i+1,j+1,k+1,el)*&
                          &f_body_matrix(i+1,j+1,k+1,el)*th_body_matrix(i+1,j+1,k+1,el)
           !          print*, product_matrix(i+1,j+1,k+1,el)
                     if(summation.lt.product_matrix(i+1,j+1,k+1,el)) summation=product_matrix(i+1,j+1,k+1,el)
                  end do
               end do
            end do
         end do
         
      end if
   end if
   placed=.FALSE.
   print*, "WAS ATOM PLACED?", placed
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
            !    !print*, product_matrix(i+1,j+1,k+1,el), product_matrix(i+1,j+1,k+1,el)**2, summation, summation**2
               
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

            !print*, product_matrix(i+1,j+1,k+1,el)
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






print*, "BUILDMAP COMPLETED"
close(1061) 





end subroutine buildmap_WIP




subroutine evaluate_contribution (element_a,element_b,&
&bondlength_target,return_slot)
character(*) :: element_a, element_b
character(1024) :: name 
integer :: atom_number_a, atom_number_b, stat, i
double precision, dimension(3) :: value_read
double precision, dimension(3) :: bondlength_read  
double precision :: a,b,c,d, bondlength_target, return_slot
write(name,'(A,A,A,A,A)') "Devolved/",trim(adjustl(element_a))&
     &,"_",trim(adjustl(element_b)),"_evolved_bondlength_gauss"
open(199, file=trim(adjustl(name)))
bondlength_read=0 
value_read=0
return_slot=0
bondlength_read(2)=bondlength_target


      master :do while(1.eq.1)

      bondlength_read(1)=bondlength_read(3)
      value_read(1)=value_read(3) 
      
      read(199,*,IOSTAT=stat) bondlength_read(3), value_read(3) 
      
      IF(IS_IOSTAT_END(stat)) then 
         return_slot=0
         close(199)
         exit master
      end IF
      IF((bondlength_read(2).lt.bondlength_read(3)).and.&
           &(bondlength_read(2).gt.bondlength_read(1))) then 
         !!Assuming the line between points can be modelled as linear
         a=bondlength_read(1)
         b=bondlength_read(2) 
         c=bondlength_read(3) 
         
         d=(b-a)/(c-a)
         
         
         
         

         if(value_read(3).lt.value_read(1)) d=1-d

         value_read(2)=value_read(1)+dble((value_read(3)-value_read(1))*d)
         !!Arbitrary could be added for when bondlength contribution starts to decay,
         !! could be changed to be bond specific (should be)!
         
         
         return_slot=value_read(2)
         if (bondlength_read(2).lt.1.6) then 
            return_slot=return_slot 
         else if (bondlength_read(2).lt.3.0) then  
            
            return_slot=return_slot
         else 
            return_slot=return_slot/26
         end if
         




         close(199)

         exit master
      end IF
      
   
   IF(IS_IOSTAT_END(stat)) exit
end do master
end subroutine evaluate_contribution

subroutine evaluate_angle_contribution (element_a,&
&bondangle_target,return_slot)
character(*) :: element_a
character(1024) :: name
integer :: atom_number_a, atom_number_b, stat, i
double precision, dimension(3) :: value_read
double precision, dimension(3) :: angles_read
double precision :: a,b,c,d, bondangle_target, return_slot


write(name,'(A,A,A,A,A)') "Devolved/",trim(adjustl(element_a))&
     &,"_evolved_angles_gauss"
open(199, file=trim(adjustl(name)))
angles_read=0
value_read=0
angles_read(2)=bondangle_target


      master :do while(1.eq.1)

      angles_read(1)=angles_read(3)
      value_read(1)=value_read(3)

      read(199,*,IOSTAT=stat) angles_read(3), value_read(3)

      IF(IS_IOSTAT_END(stat)) then
         return_slot=0
         close(199)
         exit
      end IF
      IF((angles_read(2).lt.angles_read(3)).and.&
           &(angles_read(2).gt.angles_read(1))) then
         !!Assuming the line between points can be modelled as linear
         a=angles_read(1)
         b=angles_read(2)
         c=angles_read(3)
         

         value_read(2)=(value_read(3)-value_read(1))/(angles_read(3)-angles_read(1))*(angles_read(2)-angles_read(1))+value_read(1)
         
         
         



         !!Arbitrary could be added for when bondlength contribution starts to decay,
         !! could be changed to be bond specific (should be)!


         return_slot=value_read(2)
         close(199)
         if(return_slot.gt.0) then
            !print*, value_read(1), value_read(3), d,b, return_slot
         end if

         exit master
      end IF


   IF(IS_IOSTAT_END(stat)) exit
end do master
end subroutine evaluate_angle_contribution

subroutine evaluate_4body_contribution (element_a,&
&bondangle_target,return_slot)
  character(*) :: element_a
  character(1024) :: name
  integer :: atom_number_a, atom_number_b, stat, i
  double precision, dimension(3) :: value_read
  double precision, dimension(3) :: angles_read
  double precision :: a,b,c,d, bondangle_target, return_slot


  write(name,'(A,A,A,A,A)') "Devolved/",trim(adjustl(element_a))&
       &,"_evolved_4body_gauss"
  open(199, file=trim(adjustl(name)))
  angles_read=0
  value_read=0
  angles_read(2)=bondangle_target
  if(bondangle_target.eq.1000) then 
     return_slot=1
     close(199)
  else
     master :do while(1.eq.1)

        angles_read(1)=angles_read(3)
        value_read(1)=value_read(3)

        read(199,*,IOSTAT=stat) angles_read(3), value_read(3)
        IF(IS_IOSTAT_END(stat)) then
           return_slot=0
           close(199)
           exit
        end IF
        IF((angles_read(2).lt.angles_read(3)).and.&
             &(angles_read(2).gt.angles_read(1))) then
           !!Assuming the line between points can be modelled as linear
           a=angles_read(1)
           b=angles_read(2)
           c=angles_read(3)

           d=(b-a)/(c-a)
           !print*, a, b, c

           if(value_read(3).lt.value_read(1)) d=1-d

           !value_read(2)=value_read(1)+dble((value_read(3)-value_read(1))*d)
           value_read(2)=(value_read(3)-value_read(1))/(angles_read(3)-angles_read(1))*(angles_read(2)-angles_read(1))+value_read(1)

           
           !if(value_read(2).lt.0.01) value_read(2)=0 
           
           !!Arbitrary could be added for when bondlength contribution starts to decay,
           !! could be changed to be bond specific (should be)!


           return_slot=value_read(2)
           !print*, "Angle read in, between plane and AD vector: value determined for that angle"
           if(return_slot.gt.0) then 
              !            print*, value_read(1), value_read(3), d,b, return_slot
           end if
           close(199)
           exit master
        end IF
        IF(IS_IOSTAT_END(stat)) exit
     end do master
  end if
end subroutine evaluate_4body_contribution


subroutine add_atom_scan_2 (bin_size,formula,atomlist,alistrep,atom_number_previous,structures,elrad,&
     &leng, results_matrix,eltot,elnames,placed,num_VOID)
  character(1024) :: name
  character(3), dimension(:), allocatable :: elnames
  type(unitcell), dimension(:), allocatable :: formula
  integer :: el_correct_read,i, j, k,n,l, num_VOID,atom_number_previous, new_atom_number, structures, leng,eltot
  integer, dimension(3) :: bin_size
  double precision, dimension(3) :: best_location
  double precision :: best_location_value, distribution, sigma1, sigma2, bondcutoff&
       &,agausssamp
  double precision, dimension(3) :: tmpvector
  type (atom), dimension(:,:), allocatable :: atomlist, alistrep, predicted_positions
  double precision, dimension(:,:,:), allocatable :: elrad, product_matrix
  logical :: new_position_needed, placed
  double precision, dimension(:,:,:,:,:), allocatable :: results_matrix, append_matrix
  double precision, dimension(:,:), allocatable :: sorting_matrix
  
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
              !print*, results_matrix(i+1,j+1,k+1,:,el_correct_read)
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




  !print*, "-------------------------------------------------------------"
  do i=1, (bin_size(1)+1)*(bin_size(2)+1)*(bin_size(3)+1)
     write(55,*) structures, atom_number_previous+1, sorting_matrix(i,:)
  end do
  !print*, "-------------------------------------------------------------"
  
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
     !print*, sorting_matrix(i,:), "MAT"
  end do
  call atomrepeater(structures,atomlist(structures,atom_number_previous+1)%position,&
       &atomlist(structures,atom_number_previous+1)%name,alistrep,&
       &formula,atom_number_previous+1,leng)
  
  close(55)
  
  
end subroutine add_atom_scan_2




!This routine needs to consider the effects of PLACING the atom that it might interact with itself. Hard to do
subroutine add_atom_scan (bin_size,formula,atom_number_previous, sigma1,&
     &structures,sigma2,elnames,eltot,bondcutoff,atomlist,alistrep,tmpvector,elrad,best_location)
  character(3), dimension(:), allocatable :: elnames
  type(unitcell), dimension(:), allocatable :: formula
  integer :: scan_counter, bin_size,p,q, i, j, k, best_location_index, eltot, atom_number_previous, new_atom_number, structures
  double precision, dimension(3) :: best_location
  double precision :: best_location_value, distribution, sigma1, sigma2, bond_distribution, angle_distribution, bondcutoff&
       &,agausssamp
  double precision, dimension(3) :: tmpvector
  type (atom), dimension(:,:), allocatable :: atomlist, alistrep, predicted_positions
  double precision, dimension(:,:,:), allocatable :: elrad
  logical :: new_position_needed
  allocate(predicted_positions(1,26))
  new_atom_number=atom_number_previous+1
  call atomprojector(tmpvector,predicted_positions,formula,new_atom_number,structures)
  print*, "Scanning unit cell for a good location"
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
              if (bondlength(tmpvector,alistrep(structures,q)%position).lt.1.0) cycle firstloop
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
     !     print*,"!", bondlength(best_location,alistrep(structures,i)%position)
  end do
end subroutine add_atom_scan


!This routine needs to consider the effects of PLACING the atom that it might interact with itself. Hard to do
subroutine add_atom_void (bin_size,formula,atom_number_previous, sigma1,&
     &structures,sigma2,elnames,eltot,bondcutoff,atomlist,alistrep,tmpvector,elrad,leng)
  character(3), dimension(:), allocatable :: elnames
  type(unitcell), dimension(:), allocatable :: formula
  integer :: scan_counter,leng,p,q, i, j, k, best_location_index, eltot, atom_number_previous&
       &, new_atom_number, structures,c_cut
  integer, dimension(3) :: bin_size
  double precision, dimension(3) :: best_location
  double precision :: best_location_value, smallest_bond, sigma1, sigma2, bond_distribution, angle_distribution, bondcutoff&
       &,agausssamp, comparison
  double precision, dimension(3) :: tmpvector
  type (atom), dimension(:,:), allocatable :: atomlist, alistrep, predicted_positions
  double precision, dimension(:,:,:), allocatable :: elrad
  logical :: new_position_needed
  integer, dimension(:), allocatable :: tmpdig
  character(1024), dimension(:), allocatable :: tmpels



  print*, atom_number_previous
  allocate(predicted_positions(1,26))
  new_atom_number=atom_number_previous+1
  call atomprojector(tmpvector,predicted_positions,formula,new_atom_number,structures)
  best_location_value=0
  best_location_index=0
  p=0
  call invar(15,tmpdig,tmpels)
  c_cut=tmpdig(1)
  deallocate(tmpdig)
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
              if((k.eq.0).and.(j.eq.0).and.(i.eq.0)) smallest_bond=bondlength(tmpvector,alistrep(structures,1)%position)
              p=0
              do q=1, atom_number_previous*27
                 !if (bondlength(tmpvector,alistrep(structures,q)%position).lt.1.0) cycle firstloop
!!! Experimental section for ruling out areas of unit cell
                 !if(tmpvector(3).lt.0.375*formula(structures)%cell(3,3)) cycle firstloop
                 !if(tmpvector(3).gt.0.625*formula(structures)%cell(3,3)) cycle firstloop
                 if(p.eq.0) smallest_bond=bondlength(tmpvector,alistrep(structures,1)%position)
                 p=p+1
                 if (bondlength(tmpvector,alistrep(structures,q)%position).lt.smallest_bond) then 
                    smallest_bond=bondlength(tmpvector,alistrep(structures,q)%position)
                 end if
              end do
              comparison=100.0*dble(dble(k)/dble(bin_size(3)))
              if(comparison.ge.c_cut) then
                 !print*, comparison, c_cut, "HERE"
                 smallest_bond=-1
                                
              end if

              if(smallest_bond.eq.-1) then 
                 
                 write(1001,*) tmpvector(1), tmpvector(2), tmpvector(3), 0
              else 
                 write(1001,*) tmpvector(1), tmpvector(2), tmpvector(3), smallest_bond
              end if
              
           call atomprojector(tmpvector,predicted_positions,formula,new_atom_number,structures)

           !  do q=atom_number_previous+2, 27*atom_number_previous
           !     if (bondlength(tmpvector,alistrep(structures,q)%position).lt.1.0) cycle firstloop
!!! Experimental section for ruling out areas of unit cell
           !            if(tmpvector(3).lt.0.375*formula(structures)%cell(3,3)) cycle firstloop
           !!            if(tmpvector(3).gt.0.625*formula(structures)%cell(3,3)) cycle firstloop
           !    if (bondlength(tmpvector,alistrep(structures,q)%position).lt.smallest_bond) then
           !       smallest_bond=bondlength(tmpvector,alistrep(structures,q)%position)
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
  !print*, best_location, "$^%$"
  open(56,file="void_history",access='append')
  write(56,*) structures, atom_number_previous+1,atomlist(structures,atom_number_previous+1)%position
  close(56)
  call atomrepeater(structures,atomlist(structures,atom_number_previous+1)%position,&
       &atomlist(structures,atom_number_previous+1)%name,alistrep,formula,atom_number_previous+1,leng)
  !do i=1, (new_atom_number-1)*27
  !   print*,"!", bondlength(best_location,alistrep(structures,i)%position)
  !end do
end subroutine add_atom_void







subroutine add_atom_random (formula,atom_number_previous, sigma1,&
     &structures,sigma2,elnames,eltot,bondcutoff,atomlist,alistrep,tmpvector,elrad)
  character(3), dimension(:), allocatable :: elnames
  type(unitcell), dimension(:), allocatable :: formula
  integer :: scan_counter, bin_size,p,q, i, j, k, best_location_index, eltot, atom_number_previous, new_atom_number, structures
  double precision, dimension(3) :: best_location
  double precision :: best_location_value, distribution, sigma1, sigma2, bond_distribution, angle_distribution, bondcutoff&
       &,agausssamp,r,search_region
  double precision, dimension(3) :: tmpvector
  type (atom), dimension(:,:), allocatable :: atomlist, alistrep, predicted_positions
  double precision, dimension(:,:,:), allocatable :: elrad
  logical :: new_position_needed
  allocate(predicted_positions(1,26))
  new_atom_number=atom_number_previous+1
  p=0
  scan_counter=1
  search_region=0.25

  print*, "Pseuodorandom placement"
  firstloop : do while (scan_counter.ne.0)
     do j=1, 3 
        call random_number(r)
        tmpvector(j)=r
        !        if(j.eq.3) tmpvector(j)=search_region*r
     end do
     tmpvector(:)=matmul(formula(structures)%cell,tmpvector(:))
     !     tmpvector(:)=tmpvector(:)+0.375*formula(structures)%cell(3,:)
     !     print*, tmpvector


     do q=1, atom_number_previous
        if (bondlength(tmpvector,alistrep(structures,q)%position).lt.1.0) cycle firstloop
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
        print*, "100 attemps"
        p=p-100
     end if
     distribution=bond_distribution!angle_distribution*bond_distribution
     call random_number(r)
     if (r.gt.distribution) cycle firstloop 
     scan_counter=0

  end do firstloop

end subroutine add_atom_random




subroutine add_atom_pseudo (bin_size,formula,atomlist,alistrep,atom_number_previous,structures,elrad,&
     &leng, results_matrix,eltot,elnames,placed,num_VOID)
  character(1024) :: name
  character(3), dimension(:), allocatable :: elnames
  type(unitcell), dimension(:), allocatable :: formula
  integer :: el_correct_read,i, j, k,n,l, num_VOID,atom_number_previous, new_atom_number, structures, leng,eltot
  integer, dimension(3) :: bin_size
  double precision, dimension(3) :: best_location
  double precision :: best_location_value, distribution, sigma1, sigma2, bondcutoff&
       &,agausssamp, calculated_value, uptol, lowtol, r
  double precision, dimension(3) :: tmpvector
  type (atom), dimension(:,:), allocatable :: atomlist, alistrep, predicted_positions
  double precision, dimension(:,:,:), allocatable :: elrad, product_matrix
  logical :: new_position_needed, placed
  double precision, dimension(:,:,:,:,:), allocatable :: results_matrix, append_matrix
  double precision, dimension(:,:), allocatable :: sorting_matrix
  
  uptol=1.1
  lowtol=0.95
  
  results_matrix=0
  print*, "here"
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
     
     print*, calculated_value, r, tmpvector

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
     !print*, sorting_matrix(i,:), "MAT"
  end do
  call atomrepeater(structures,atomlist(structures,atom_number_previous+1)%position,&
       &atomlist(structures,atom_number_previous+1)%name,alistrep,&
       &formula,atom_number_previous+1,leng)
  
  close(55)
  
  
end subroutine add_atom_pseudo





subroutine generate_bondangle_distribution (tmpvector,atomlist,alistrep,structures,bondcutoff, &
     &atom_number_new, new_position_needed&
     &,agausssamp,distribution,eltot, sigma1, sigma2,elnames, predicted_positions)
  character(3), dimension(:), allocatable :: elnames
  double precision, dimension(3) :: tmpvector
  integer :: y, atom_number_new, structures, eltot,i,j,k,L, nlines 
  logical :: new_position_needed, file_exists
  type (atom), dimension(:,:), allocatable :: atomlist, alistrep, predicted_positions
  type (nearmatrix), dimension(:), allocatable :: nearneighbourmatrix
  double precision :: agausssamp, bondcutoff, bondangle_store, distribution, read_angle, sigma1, sigma2, step_significance
  double precision :: tmpvalue, total_significance
  character(1024) :: determined_anglefile
  double precision, allocatable, dimension(:) :: initial_angle_peaks

  distribution=0
  !if(allocated(nearneighbourmatrix)) deallocate(nearneighbourmatrix) 
  allocate(nearneighbourmatrix(((atom_number_new-1)*27)+26))
  k=0

  do y=1, ((atom_number_new-1)*27)+26
     if(y.le.((atom_number_new-1)*27)) then 
        if(bondlength(tmpvector,alistrep(structures,y)%position).lt.bondcutoff) then
           k=k+1
           nearneighbourmatrix(k)%position=alistrep(structures,y)%position
        end if
     else
        if(bondlength(tmpvector,predicted_positions(1,y-(atom_number_new-1)*27)%position).lt.bondcutoff) then
           k=k+1
           nearneighbourmatrix(k)%position=predicted_positions(1,y-(atom_number_new-1)*27)%position
        end if
     end if
  end do


!!! CHANGE alistrep etc to nearneighbourmatrix%position. SHould work.
  total_significance=0
  if(new_position_needed.neqv..TRUE.) then 
     do L=1, eltot
        !print*, elnames(L), structures, atom_number_previous, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        !print*, atomlist(structures,atom_number_previous+1)%name, "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"
        if(trim(adjustl(elnames(L))).ne.trim(adjustl(atomlist(structures,atom_number_new)%name))) cycle

        do j=1, k 
           do y=1, k 
              if (j.eq.y) cycle 
              bondangle_store=bondangle(nearneighbourmatrix(j)%position,tmpvector,&
                   &nearneighbourmatrix(k)%position)
              if (bondlength(tmpvector,nearneighbourmatrix(k)%position).lt.1.0) then 
                 new_position_needed=.TRUE.
                 cycle
              end if
              if (bondlength(tmpvector,nearneighbourmatrix(k)%position).lt.1.0) then
                 new_position_needed=.TRUE.
                 cycle
              end if

              if(new_position_needed.eqv..TRUE.) cycle 

              CALL angledistribution(elnames(L),initial_angle_peaks)
              do i=1, size(initial_angle_peaks)
                 step_significance=1
                 total_significance=total_significance+step_significance
                 distribution=distribution+(1.0/2.0)*((erf((tmpvalue-initial_angle_peaks(i)+&
                      &(0.5*agausssamp))/(sigma2*sqrt(2.0)))-&
                      &erf((tmpvalue-initial_angle_peaks(i)-(0.5*agausssamp))/(sigma2*sqrt(2.0)))))
              end do

              write(determined_anglefile,'(A,A)') trim(adjustL(elnames(L))),"_evolved_angles"
              INQUIRE(FILE=determined_anglefile, EXIST=file_exists)
              if (file_exists.eqv..TRUE.) then
                 nlines=0
                 OPEN (11, file=trim(adjustl(determined_anglefile)))
                 DO
                    READ (11,*, END=10) read_angle, step_significance
                    total_significance=total_significance+step_significance
                    distribution=distribution+step_significance*(1.0/2.0)*(&
                         &abs(erf((bondangle_store-read_angle+0.5*sigma1)/((sigma1)*2**dble(1.0/2.0)))-&
                         &erf((bondangle_store-read_angle-0.5*sigma1)/((sigma1)*2**dble(1.0/2.0)))))
                    nlines=nlines+1
                 END DO
10               CLOSE (11)
              end if
           end do
        end do
     end do
     distribution=abs(distribution/total_significance) 
  end if







end subroutine generate_bondangle_distribution




subroutine generate_bondlength_distribution (sigma1,sigma2,atom_number_proposed,distribution,atomlist,alistrep,&
     &tmpvector,elrad,eltot,structures,elnames,new_position_needed, predicted_positions)
  double precision, dimension(:,:,:), allocatable :: elrad
  integer :: nlines, atom_number_proposed, atom_number_previous, atom_total_repeated,j,k,l,eltot,structures
  type (atom), dimension(:,:), allocatable :: atomlist, alistrep, predicted_positions
  double precision, dimension(3) :: tmpvector
  double precision :: distribution, total_significance, step_significance, bondlength_store, sigma1,sigma2, read_bondlength 
  character(1024) :: determined_bondfile 
  logical :: file_exists, new_position_needed
  character(3), dimension(:), allocatable :: elnames

  atom_number_previous=atom_number_proposed-1
  atom_total_repeated=atom_number_previous*27
  total_significance=0
  new_position_needed=.FALSE.
  do while(new_position_needed.eqv..FALSE.)
     do j=1, atom_total_repeated
        do k=1, eltot 
           do L=1, eltot
              !! Two eltot loops in here to make sure the atom being added and the atom being checked are the correct species 
              if(atomlist(structures,atom_number_proposed)%name.ne.elnames(k)) cycle
              if(alistrep(structures,j)%name.ne.elnames(L)) cycle
              step_significance=1
              bondlength_store=bondlength(tmpvector,alistrep(structures,j)%position)
              if (bondlength_store.lt.0.7*elrad(1,L,k)) then 
                 new_position_needed=.TRUE.
              end if
              total_significance=total_significance+step_significance
              distribution=distribution+(1.0/2.0)*(&
                   &abs(erf((bondlength_store-elrad(1,k,L)+0.5*sigma2)/((sigma1)*2**dble(1.0/2.0)))-&
                   &erf((bondlength_store-elrad(1,k,L)-0.5*sigma2)/((sigma1)*2**dble(1.0/2.0)))))
              file_exists=.FALSE.

              write(determined_bondfile,'(A,A,A,A)') &
                   &trim(adjustL(elnames(k))),"_",trim(adjustL(elnames(L))),"_evolved_bondlength"
              INQUIRE(FILE=determined_bondfile, EXIST=file_exists)           
              if (file_exists .eqv..FALSE.) then  
                 write(determined_bondfile,'(A,A,A,A)') &
                      &trim(adjustL(elnames(L))),"_",trim(adjustL(elnames(K))),"_evolved_bondlength"
                 INQUIRE(FILE=determined_bondfile, EXIST=file_exists)
              end if

              if (file_exists .eqv..TRUE.) then 
                 nlines=0
                 OPEN (11, file=determined_bondfile)
                 DO
                    READ (11,*, END=10) read_bondlength, step_significance
                    !           print*, distribution, nlines, total_significance 
                    total_significance=total_significance+step_significance
                    distribution=distribution+step_significance*(1.0/2.0)*(&
                         &abs(erf((bondlength_store-read_bondlength+0.5*sigma2)/((sigma1)*2**dble(1.0/2.0)))-&
                         &erf((bondlength_store-read_bondlength-0.5*sigma2)/((sigma1)*2**dble(1.0/2.0)))))
                    nlines=nlines+1
                 END DO
10               CLOSE (11)
              end if
           end do
        end do
     end do
     !distribution=abs(distribution/total_significance)
     exit
  end do

  do while(new_position_needed.eqv..FALSE.)
     do j=1, 26
        do k=1, eltot
           !! Two eltot loops in here to make sure the atom being added and the atom being checked are the correct species
           if(atomlist(structures,atom_number_proposed)%name.ne.elnames(k)) cycle
           step_significance=1
           bondlength_store=bondlength(tmpvector,predicted_positions(1,j)%position)
           if (bondlength_store.lt.1.0) new_position_needed=.TRUE.
           total_significance=total_significance+step_significance
           distribution=distribution+(1.0/2.0)*(&
                &abs(erf((bondlength_store-elrad(1,k,k)+0.5*sigma2)/((sigma1)*2**dble(1.0/2.0)))-&
                &erf((bondlength_store-elrad(1,k,k)-0.5*sigma2)/((sigma1)*2**dble(1.0/2.0)))))
           file_exists=.FALSE.

           write(determined_bondfile,'(A,A,A,A)') &
                &trim(adjustL(elnames(k))),"_",trim(adjustL(elnames(k))),"_evolved_bondlength"
           INQUIRE(FILE=determined_bondfile, EXIST=file_exists)
           if (file_exists .eqv..TRUE.) then
              OPEN (11, file=determined_bondfile)
              DO
                 READ (11,*, END=11) read_bondlength, step_significance
                 !           print*, distribution, nlines, total_significance
                 total_significance=total_significance+step_significance
                 distribution=distribution+step_significance*(1.0/2.0)*(&
                      &abs(erf((bondlength_store-read_bondlength+0.5*sigma2)/((sigma1)*2**dble(1.0/2.0)))-&
                      &erf((bondlength_store-read_bondlength-0.5*sigma2)/((sigma1)*2**dble(1.0/2.0)))))
                 nlines=nlines+1
              END DO
11            CLOSE (11)
           end if
        end do
     end do
     distribution=abs(distribution/total_significance)

     exit
  end do

end subroutine generate_bondlength_distribution



subroutine touchpos() 

  character(1024) :: file, command 
  logical :: dir_e

  inquire(file="pos", exist=dir_e)
  if(dir_e) then
  else
     write(command,*) "mkdir pos"
     CALL execute_command_line(command)
  end if
end subroutine touchpos


function ReadLine(aunit, InLine, trimmed) result(OK)
integer, intent(IN) :: aunit
character(LEN=:), allocatable, optional :: InLine
logical, intent(in), optional :: trimmed
integer, parameter :: line_buf_len= 1024*4
character(LEN=line_buf_len) :: InS
logical :: OK, set
integer status, size

OK = .false.
set = .true.
do
    read (aunit,'(a)',advance='NO',iostat=status, size=size) InS
    OK = .not. IS_IOSTAT_END(status)
    if (.not. OK) return
    if (present(InLine)) then
        if (set) then
            InLine = InS(1:size)
            set=.false.
        else
            InLine = InLine // InS(1:size)
        end if
    end if
    if (IS_IOSTAT_EOR(status)) exit
end do
if (present(trimmed) .and. present(InLine)) then
    if (trimmed) InLine = trim(adjustl(InLine))
end if

end function ReadLine



subroutine poscar_read(pathway,unit_cell,structure_elements,structure_stochiometry, structure_factor,&
     &coordinate_type,atomic_positions, header)
  use atomtype 
  implicit none 
  character(1024), intent(out) :: coordinate_type, header 
  character(1024) :: read_in, tmp
  character(1024), intent(in) :: pathway
  character(3), dimension(:), allocatable, intent(out) :: structure_elements
  character(3), dimension(:), allocatable :: temp_names
  integer, dimension(:), allocatable, intent(out) :: structure_stochiometry 
  integer :: loop, status, size,eltot, i, j,k,l,m
  integer, parameter :: line_buf_len= 1024*4
  character(LEN=line_buf_len) :: InS
  type(unitcell), dimension(1), intent(out) :: unit_cell 
  double precision, dimension(:,:), intent(out), allocatable :: atomic_positions
  double precision, dimension(:,:), allocatable :: temp_positions
  double precision, intent(out) :: structure_factor 
  logical :: file_e  
  
  eltot=0
  allocate(structure_elements(1))
  inquire(file=trim(adjustl(pathway)), exist=file_e)
  print*, file_e
  if(file_e) then 
     open(102,file=pathway)
     
     read(102,*) header
     read(102,*) structure_factor
     do loop=1, 3
        read(102,*) unit_cell(1)%cell(loop,1), unit_cell(1)%cell(loop,2), unit_cell(1)%cell(loop,3)
     end do
     read(102,'(A)') tmp
     print*, tmp
     write(read_in,'(X,A)') trim(adjustl(tmp))
     print*, read_in
     do i=1,len(read_in)-1
     if(i.eq.1) then
        if((scan(read_in(i:i+1)," ").eq.0).or.&
             &((scan(read_in(i:i+1)," ").eq.2))) then
           eltot=eltot+1
           if(scan(read_in(i:i+1)," ").eq.0) then
              allocate(temp_names(eltot))
              
              do j=1, eltot-1
                 temp_names(j)=structure_elements(j)
              end do
              temp_names(eltot)=read_in(i:i+1)
              deallocate(structure_elements)
              allocate(structure_elements(eltot))
              structure_elements=temp_names
              deallocate(temp_names)

           end if
           if((scan(read_in(i:i+1)," ").eq.2).and.(scan(read_in(i-1:i)," ").eq.1)) then
              allocate(temp_names(eltot))

              do j=1, eltot-1
                 temp_names(j)=structure_elements(j)
              end do
              temp_names(eltot)=read_in(i:i)
              deallocate(structure_elements)
              allocate(structure_elements(eltot))
              structure_elements=temp_names
              deallocate(temp_names)
           end if
        else
        end if
     else
        if((scan(read_in(i:i+1)," ").eq.0).or.&
             &((scan(read_in(i:i+1)," ").eq.2).and.(scan(read_in(i-1:i)," ")&
             &.eq.1))) then
           eltot=eltot+1
           !print*, tmp(i:i+1)
           if(scan(read_in(i:i+1)," ").eq.0) then
              allocate(temp_names(eltot))

              do j=1, eltot-1
                 temp_names(j)=structure_elements(j)
              end do
              temp_names(eltot)=read_in(i:i+1)
              deallocate(structure_elements)
              allocate(structure_elements(eltot))
              structure_elements=temp_names
              deallocate(temp_names)


           end if
           if((scan(read_in(i:i+1)," ").eq.2).and.(scan(read_in(i-1:i)," ").eq.1)) then
              allocate(temp_names(eltot))
              
              do j=1, eltot-1
                 temp_names(j)=structure_elements(j)
              end do
              temp_names(eltot)=read_in(i:i)
              deallocate(structure_elements)
              allocate(structure_elements(eltot))
              structure_elements=temp_names
              deallocate(temp_names)

           end if
        else
        end if
     end if
  end do
  allocate(structure_stochiometry(eltot))
  
  read(102,*) structure_stochiometry 
  read(102,*) coordinate_type
  allocate(temp_positions(sum(structure_stochiometry),3))
  allocate(atomic_positions(sum(structure_stochiometry),3))
  if(coordinate_type.eq."Cartesian") then 
     do i=1, sum(structure_stochiometry)
        read(102,*) atomic_positions(i,1),atomic_positions(i,2),atomic_positions(i,3)
     end do
  else
     print*, coordinate_type
     do i=1, sum(structure_stochiometry)
                
        read(102,*) temp_positions(i,1),temp_positions(i,2),temp_positions(i,3)
        
     end do
     do i=1, sum(structure_stochiometry)
        atomic_positions(i,:)=matmul(temp_positions(i,:),unit_cell(1)%cell) 
!temp_positions(i,1)*unit_cell(1,1)+temp_positions(i,2)*unit_cell() 
        !print*, atomic_positions(i,:)
        !print*, temp_positions(i,:) 
        !print*, unit_cell(i)%cell 
        !call sleep(5) 
     end do
  end if
else
   print*, "ERROR: NOT FOUND"
end if

end subroutine poscar_read

subroutine regenerate_distribution_files (prev_structures)
  use atomtype
  implicit none
  character(1024) :: command,name,tmp, location_string
  integer :: l,i,structures,prev_structures,stage, eltot, atomnumber, tmpint,m,j,k, structure_loop
  type (atom), dimension(:,:), allocatable :: array, repeatedarray
  integer, dimension(:), allocatable :: stochio
  logical :: dir_e, file_e, empty
  character(1024), dimension(:), allocatable :: tmpels
  integer, dimension(:), allocatable :: tmpdig


  
  character(1024) :: coordinate_type, header
  character(LEN=:), allocatable :: read_in
  character(1024) :: pathway, path_test
  character(3), dimension(:), allocatable :: structure_elements
  integer, dimension(:), allocatable :: structure_stochiometry
  integer :: loop, status, size, stat
  integer, parameter :: line_buf_len= 1024*4
  character(LEN=line_buf_len) :: InS
  type(unitcell), dimension(:), allocatable :: unit_cell
  double precision, dimension(:,:), allocatable :: atomic_positions
  double precision :: structure_factor
  double precision, dimension(:,:), allocatable :: bondcutoff
  logical ::  OK, set
  
allocate(unit_cell(1))
  write(pathway,*) "rm -r bon bad 4body bon2 bad2 4body2 bon3 bad3 4body3"
  call execute_command_line(pathway)
  strucloop:do structure_loop=1, prev_structures
  !DEBUG
  !strucloop:do structure_loop=303, 304 
    write(pathway,'(A,I0.3,A)') "pos/POSCAR_",structure_loop,"/KPOINTS" 
     !print*, pathway
     inquire(file=trim(adjustl(pathway)),exist=dir_e)
     if(dir_e.eqv..FALSE.) cycle strucloop
     write(pathway,'(A,I0.3,A)') "pos/POSCAR_",structure_loop,"/CONTCAR"
     inquire(file=trim(adjustl(pathway)),exist=dir_e)
     write(command,'(A,A,A,I0.3,A,I0.3,A,I0.3,A,I0.3,A)') &
          &"if [ -s ",trim(adjustl(pathway)), &
          &" ]; then rm -f  pos/POSCAR_",structure_loop,&
          &"/empty.txt; touch pos/POSCAR_",structure_loop&
          &,"/full.txt; else rm -f touch pos/POSCAR_",structure_loop&
          &,"/full.txt; touch pos/POSCAR_",structure_loop,&
          &"/empty.txt; fi"
     print*, command
     call execute_command_line(command)
     write(path_test,'(A,I0.3,A)') "pos/POSCAR_",structure_loop,"/full.txt"

     inquire(file=trim(adjustl(path_test)),exist=dir_e)

     if(dir_e.eqv..FALSE.) then
        write(pathway,'(A,I0.3,A)') "pos/POSCAR_",structure_loop,"/POSCAR"
        inquire(file=trim(adjustl(pathway)),exist=dir_e)
        if(dir_e.eqv..FALSE.) cycle strucloop
     end if


     print*, pathway, structure_loop 
     

     call poscar_read(pathway,unit_cell(1),structure_elements,structure_stochiometry, structure_factor,coordinate_type&
          &,atomic_positions, header)
     !print*, pathway
     allocate(array(1,sum(structure_stochiometry)))
     allocate(repeatedarray(1,27*sum(structure_stochiometry)))
     k=0
     !print*, size(array)
     do i=1, size(structure_stochiometry)
        do j=1, structure_stochiometry(i) 
           k=k+1
           array(1,k)%name=structure_elements(i)
           array(1,k)%position=atomic_positions(k,:)
        end do
        

     end do
     

     allocate(bondcutoff(size(structure_stochiometry),size(structure_stochiometry)))
     do i=1, size(structure_stochiometry)
        do j=1, size(structure_stochiometry)
           if(j.gt.i) cycle
           write(command,*) "echo $(sed -n 's|",trim(adjustl(structure_elements(i)))&
                &," ",trim(adjustl(structure_elements(j)))," ||p'<chem.in)|",&
                &"sed -n 's| .*||p'>tmp.out"
           call execute_command_line(command) 
           open(72,file="tmp.out")
           read(72,*,IOSTAT=stat) bondcutoff(i,j)
           IF(IS_IOSTAT_END(stat)) then 
              close(72)
              write(command,*) "echo $(sed -n 's|",trim(adjustl(structure_elements(j)))&
                &," ",trim(adjustl(structure_elements(i)))," ||p'<chem.in)|",&
                &"sed -n 's| .*||p'>tmp.out"
              call execute_command_line(command)
              open(72,file="tmp.out")
              read(72,*) bondcutoff(i,j)
              bondcutoff(i,j)=bondcutoff(i,j)*1.2
              bondcutoff(j,i)=bondcutoff(i,j)
              close(72)
           else

              close(72)
              bondcutoff(i,j)=bondcutoff(i,j)*1.2
              bondcutoff(j,i)=bondcutoff(i,j)
           end if
        end do
     end do
     
     


     print*, "MANUALLY CHANGE"
     do i=1, sum(structure_stochiometry)
        call atomrepeater(1,array(1,i)%position,array(1,i)%name,&
             &repeatedarray,unit_cell,i,sum(structure_stochiometry)) 
     end do

     
     


     print*, "YOU SHOULD ADD IN HERE A FUNCTION WHICH CERATES DYNAMIC BONDCUTOFF VALUES FOR DIFFERENT ELEMTNS."
     print*, structure_stochiometry, sum(structure_stochiometry)
     do i=1, sum(structure_stochiometry) 
        call generatebondfiles(1,1,structure_loop-1,array,repeatedarray,size(structure_stochiometry),&             
             &structure_stochiometry,i,structure_elements)
        call generate4files(1,1,structure_loop-1,array,repeatedarray,size(structure_stochiometry),&
             &structure_stochiometry,i,bondcutoff,structure_elements)
        
        call generateanglefiles(1,1,structure_loop-1,array,repeatedarray,size(structure_stochiometry),&
             &structure_stochiometry,i,bondcutoff,structure_elements)
      end do
     deallocate(array)
     deallocate(repeatedarray)
     deallocate(bondcutoff)
  end do strucloop
  
  strucloop2:do structure_loop=1, prev_structures
  !  strucloop2:do structure_loop=303, 304
   
     write(pathway,'(A,I0.3,A)') "pos/POSCAR_",structure_loop,"/RELAX/CONTCAR" 
    
     write(command,'(A,A,A,A,A)') " if [ -s ", trim(adjustl(pathway))," ]; then :; else rm ",trim(adjustl(pathway))," ; fi"
     call execute_command_line(command)
     inquire(file=trim(adjustl(pathway)),exist=dir_e)
    
     if(dir_e.eqv..FALSE.) then 
        write(pathway,'(A,I0.3,A)') "pos/POSCAR_",structure_loop,"/RELAX/POSCAR"
        
        write(command,'(A,A,A,A,A)') " if [ -s ", trim(adjustl(pathway))," ]; then :; else rm ",trim(adjustl(pathway))," ; fi"
        call execute_command_line(command)
        inquire(file=trim(adjustl(pathway)),exist=dir_e)
     
        if(dir_e.eqv..FALSE.) cycle strucloop2
     else
        continue 
     end if
     




     call poscar_read(pathway,unit_cell(1),structure_elements,structure_stochiometry, structure_factor,coordinate_type&
          &,atomic_positions, header)
     print*, pathway
     allocate(array(1,sum(structure_stochiometry)))
     allocate(repeatedarray(1,27*sum(structure_stochiometry)))
     k=0
     print*, size(array)
     do i=1, size(structure_stochiometry)
        do j=1, structure_stochiometry(i) 
           k=k+1
           array(1,k)%name=structure_elements(i)
           array(1,k)%position=atomic_positions(k,:)
        end do
     end do
     
     allocate(bondcutoff(size(structure_stochiometry),size(structure_stochiometry)))
     do i=1, size(structure_stochiometry)
        do j=1, size(structure_stochiometry)
           if(j.gt.i) cycle
           write(command,*) "echo $(sed -n 's|",trim(adjustl(structure_elements(i)))&
                &," ",trim(adjustl(structure_elements(j)))," ||p'<chem.in)|",&
                &"sed -n 's| .*||p'>tmp.out"
           call execute_command_line(command)
           open(72,file="tmp.out")
           read(72,*,IOSTAT=stat) bondcutoff(i,j)
           IF(IS_IOSTAT_END(stat)) then
              close(72)
              write(command,*) "echo $(sed -n 's|",trim(adjustl(structure_elements(j)))&
                &," ",trim(adjustl(structure_elements(i)))," ||p'<chem.in)|",&
                &"sed -n 's| .*||p'>tmp.out"
              call execute_command_line(command)
              open(72,file="tmp.out")
              read(72,*) bondcutoff(i,j)
              bondcutoff(j,i)=bondcutoff(i,j)
              bondcutoff(i,j)=bondcutoff(i,j)*1.1
              close(72)
           else

              close(72)
              bondcutoff(j,i)=bondcutoff(i,j)
              bondcutoff(i,j)=bondcutoff(i,j)*1.1
           end if
        end do
     end do

     do i=1, sum(structure_stochiometry)
        call atomrepeater(1,array(1,i)%position,array(1,i)%name,&
             &repeatedarray,unit_cell,i,sum(structure_stochiometry)) 
     end do

     print*, "YOU SHOULD ADD IN HERE A FUNCTION WHICH CERATES DYNAMIC BONDCUTOFF VALUES FOR DIFFERENT ELEMTNS."
     do i=1, sum(structure_stochiometry) 
        print*, prev_structures,structure_loop,prev_structures+structure_loop
        call generatebondfiles(2,1,structure_loop-1,array,repeatedarray,size(structure_stochiometry),&
             &structure_stochiometry,i,structure_elements)
        print*, "11"
        call generate4files(2,1,structure_loop-1,array,repeatedarray,size(structure_stochiometry),&
             &structure_stochiometry,i,bondcutoff,structure_elements)
        print*, "22"
        call generateanglefiles(2,1,structure_loop-1,array,repeatedarray,size(structure_stochiometry),&
             &structure_stochiometry,i,bondcutoff,structure_elements)
        print*, "33"
     end do
     deallocate(array)
     deallocate(repeatedarray)
     deallocate(bondcutoff)
  end do strucloop2
  strucloop3:do structure_loop=1, prev_structures
  ! strucloop3:do structure_loop=303, 304

     write(pathway,'(A,I0.3,A)') "pos/POSCAR_",structure_loop,"/RELAX/RELAX2/CONTCAR" 
     print*, pathway
      
     write(command,'(A,A,A,A,A)') " if [ -s ", trim(adjustl(pathway))," ]; then :; else rm ",trim(adjustl(pathway))," ; fi"
     call execute_command_line (command)

     inquire(file=trim(adjustl(pathway)),exist=dir_e)
     
     
     if(dir_e.eqv..FALSE.) then 
        write(pathway,'(A,I0.3,A)') "pos/POSCAR_",structure_loop,"/RELAX/RELAX2/POSCAR"
        
        write(command,'(A,A,A,A,A)') " if [ -s ", trim(adjustl(pathway))," ]; then :; else rm ",trim(adjustl(pathway))," ; fi"
        call execute_command_line(command)
        inquire(file=trim(adjustl(pathway)),exist=dir_e)
     
        if(dir_e.eqv..FALSE.) cycle strucloop3
     else
        continue 
     end if

     call poscar_read(pathway,unit_cell(1),structure_elements,structure_stochiometry, structure_factor,coordinate_type&
          &,atomic_positions, header)
     print*, pathway
     allocate(array(1,sum(structure_stochiometry)))
     allocate(repeatedarray(1,27*sum(structure_stochiometry)))
     k=0
     print*, size(array)
     do i=1, size(structure_stochiometry)
        do j=1, structure_stochiometry(i) 
           k=k+1
           array(1,k)%name=structure_elements(i)
           array(1,k)%position=atomic_positions(k,:)
        end do
     end do
     
     allocate(bondcutoff(size(structure_stochiometry),size(structure_stochiometry)))
     do i=1, size(structure_stochiometry)
        do j=1, size(structure_stochiometry)
           if(j.gt.i) cycle
           write(command,*) "echo $(sed -n 's|",trim(adjustl(structure_elements(i)))&
                &," ",trim(adjustl(structure_elements(j)))," ||p'<chem.in)|",&
                &"sed -n 's| .*||p'>tmp.out"
           call execute_command_line(command)
           open(72,file="tmp.out")
           read(72,*,IOSTAT=stat) bondcutoff(i,j)
           IF(IS_IOSTAT_END(stat)) then
              close(72)
              write(command,*) "echo $(sed -n 's|",trim(adjustl(structure_elements(j)))&
                &," ",trim(adjustl(structure_elements(i)))," ||p'<chem.in)|",&
                &"sed -n 's| .*||p'>tmp.out"
              call execute_command_line(command)
              open(72,file="tmp.out")
              read(72,*) bondcutoff(i,j)
              bondcutoff(i,j)=bondcutoff(i,j)*1.1
                            
              bondcutoff(j,i)=bondcutoff(i,j)
              
              close(72)
           else

              close(72)
              bondcutoff(i,j)=bondcutoff(i,j)*1.1
              bondcutoff(j,i)=bondcutoff(i,j)
             end if
        end do
     end do
     print*, bondcutoff
     stop

     do i=1, sum(structure_stochiometry)
        call atomrepeater(1,array(1,i)%position,array(1,i)%name,&
             &repeatedarray,unit_cell,i,sum(structure_stochiometry)) 
     end do

     print*, "YOU SHOULD ADD IN HERE A FUNCTION WHICH CERATES DYNAMIC BONDCUTOFF VALUES FOR DIFFERENT ELEMTNS."
     do i=1, sum(structure_stochiometry) 
        print*, prev_structures,structure_loop,prev_structures+structure_loop
        call generatebondfiles(3,1,structure_loop-1,array,repeatedarray,size(structure_stochiometry),&
             &structure_stochiometry,i,structure_elements)
        print*, "11"
        call generate4files(3,1,structure_loop-1,array,repeatedarray,size(structure_stochiometry),&
             &structure_stochiometry,i,bondcutoff,structure_elements)
        print*, "22"
        call generateanglefiles(3,1,structure_loop-1,array,repeatedarray,size(structure_stochiometry),&
             &structure_stochiometry,i,bondcutoff,structure_elements)
        print*, "33"
     end do
     deallocate(array)
     deallocate(repeatedarray)
     deallocate(bondcutoff)
  end do strucloop3
  


end subroutine regenerate_distribution_files

subroutine bond_evolution(mode)
  character(1024) :: name, read_element_pairing, read_element
  integer :: prev_structures, mode,i, nbin, stat, exitst, exitst2, exitst3
  character(1024), dimension(:), allocatable :: tmpels
  integer, dimension(:), allocatable :: tmpdig
  logical :: dir_e
  double precision, dimension(:,:), allocatable :: gaussian
  double precision, dimension(2) :: read_in, norma_vector
  double precision :: sigma, bondcut, dist_height

  if(mode.eq.1) then 
     write(name,'(A)') "rm tmp_energies.txt"
     call execute_command_line(name,wait=.TRUE.) 
     write(name,'(A)') "prevstructures.txt" 
     open(13,file=name,status="old") 
     read(13, *) prev_structures
     print*, prev_structures, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
     do i=1, prev_structures
        write(name,'(A,I0.3,A,I0.3,A)') "echo ",i,"'.1 '$(echo $(echo $(grep 'e  e' pos/POSCAR_",i,"/OUTCAR) &
             | sed -n 's|^.*= *||p') | sed -n 's| .*||p')>>tmp_energies.txt"
        call execute_command_line(name,wait=.TRUE.,exitstat=exitst)
        do  
           if(exitst.eq.0) exit
        end do
        write(name,'(A,I0.3,A,I0.3,A)') "echo ",i,"'.2 '$(echo $(echo $(grep 'e  e' pos/POSCAR_",i,"/RELAX/OUTCAR) &
             | sed -n 's|^.*= *||p') | sed -n 's| .*||p')>>tmp_energies.txt"
        call execute_command_line(name,wait=.TRUE.,exitstat=exitst2)
        !do 
        !   if(exitst2.eq.0) exit
        !end do!
        !write(name,'(A,I0.3,A)') "pos/POSCAR_",i,"/RELAX"
        
        write(name,'(A,I0.3,A,I0.3,A)') "echo ",i,"'.3 '$(echo $(echo $(grep 'e  e' pos/POSCAR_",i,"/RELAX/RELAX2/&
             &OUTCAR) | sed -n 's|^.*= *||p') | sed -n 's| .*||p')>>tmp_energies.txt"
        call execute_command_line(name,wait=.TRUE.,exitstat=exitst3)
        
        do 
           if(exitst3.eq.0) exit
        end do!
        
        
     end do
     close(13)
  end if
  if((mode.eq.1).or.(mode.eq.2)) then
     print*, "#################################################"
     call execute_command_line("./bond_evolution.sh tmp_energies.txt prevstructures.txt",wait=.TRUE.)
  end if
  !  The arguements in the next function are bond cap, bins and element 1 and 2. These should be changed from carbon et in this file when finished 
  !  The sigma here should be related to the sampling of the placement map; too small compared to map, and expanding regions to study in more 
  !  detail will be difficult
  bondcut=5
  nbin=1000
  sigma=0.1
  !dist_height=1/(sigma*(sqrt(2*3.141592654)))
  allocate(gaussian(2,nbin))
  call execute_command_line("mkdir Devolved; rm Devolved/*evolved_*_gauss; ./bond_database_builder.sh")
  write(name,*) "bond_element_tempfile"
  open(103,file=trim(adjustl(name)))
  write(name,'(A,A,A)') "Devolved/ALL_evolved_bondlength_gauss"
  open(105,file=trim(adjustl(name)))

  do while(1.eq.1) 
     read(103,*,IOSTAT=stat) read_element_pairing
     IF(IS_IOSTAT_END(stat)) exit 
     do i=1, nbin
        gaussian(1,i)=dble(i*bondcut/nbin)
        gaussian(2,i)=0
     end do


     write(name,'(A,A)') trim(adjustl(read_element_pairing)), "_evolved_bondlength"
     open(101,file=trim(adjustl(name)))
     print*, trim(adjustl(name)), "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
     write(name,'(A,A,A)') "Devolved/", trim(adjustl(read_element_pairing)), "_evolved_bondlength_gauss"
     open(102,file=trim(adjustl(name)))

     do while(1.eq.1)
        read(101,*,IOSTAT=stat) read_in
        !print*, read_in
        IF(IS_IOSTAT_END(stat)) exit

        do i=1, nbin
           dist_height=1/((1+gaussian(2,i)))**2
           gaussian(2,i)=gaussian(2,i)+dist_height*read_in(2)*(1/(sigma*(2*3.14159)**0.5))*&
                exp(-0.5*((gaussian(1,i)-read_in(1))/sigma)**2)
        end do
     end do
     write(105,*) trim(adjustl(read_element_pairing))
     

     dist_height=maxval(gaussian(2,:))
     
     do i=1, nbin 
        write(102,*) gaussian(1,i), gaussian(2,i)/(dist_height)
        write(105,*) gaussian(1,i), gaussian(2,i)/(dist_height)
        
     end do
     close(102) 
     close(101)
    end do
  close(103)

  sigma=0.05
  write(name,*) "angle_element_tempfile"
  open(103,file=trim(adjustl(name)))

  do while(1.eq.1)
     read(103,*,IOSTAT=stat) read_element
     IF(IS_IOSTAT_END(stat)) exit
     do i=1, nbin
        gaussian(1,i)=dble(i*3.141592653/(nbin))
        gaussian(2,i)=0
     end do
     
     
     write(name,'(A,A)') trim(adjustl(read_element)), "_evolved_angles"
     open(101,file=trim(adjustl(name)))
     print*, trim(adjustl(name)), "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
     write(name,'(A,A,A)') "Devolved/", trim(adjustl(read_element)), "_evolved_angles_gauss"
     open(102,file=trim(adjustl(name)))

     do while(1.eq.1)
        read(101,*,IOSTAT=stat) read_in
      
        if(isnan(read_in(1))) cycle 
        IF(IS_IOSTAT_END(stat)) exit
        
        do i=1, nbin
           dist_height=1/(1+gaussian(2,i))
           gaussian(2,i)=gaussian(2,i)+dist_height*read_in(2)*(1/(sigma*(2*3.14159)**0.5))*&
                exp(-0.5*((gaussian(1,i)-read_in(1))/sigma)**2)
        end do
     end do
     norma_vector=maxval(gaussian,2)
     do i=1, nbin
        write(102,*) gaussian(1,i), gaussian(2,i)/norma_vector(2)
        
     end do
     close(102)
     close(101)
  end do

  rewind(103)
  sigma=0.05
  do while(1.eq.1)
     read(103,*,IOSTAT=stat) read_element
     IF(IS_IOSTAT_END(stat)) exit
     do i=1, nbin
        gaussian(1,i)=dble(i*3.141592653/(nbin))
        gaussian(2,i)=0
     end do


     write(name,'(A,A)') trim(adjustl(read_element)), "_evolved_4body"
     open(101,file=trim(adjustl(name)))
     print*, trim(adjustl(name)), "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
     write(name,'(A,A,A)') "Devolved/", trim(adjustl(read_element)), "_evolved_4body_gauss"
     open(102,file=trim(adjustl(name)))





     do while(1.eq.1)
        read(101,*,IOSTAT=stat) read_in
        !print*, read_in
        IF(IS_IOSTAT_END(stat)) exit
        
        do i=1, nbin         
           dist_height=1/(1+gaussian(2,i))
           gaussian(2,i)=gaussian(2,i)+dist_height*read_in(2)*(1/(sigma*(2*3.14159)**0.5))*&
                exp(-0.5*((gaussian(1,i)-read_in(1))/sigma)**2)
        end do
     end do
     norma_vector=maxval(gaussian,2)
     do i=1, nbin
        write(102,*) gaussian(1,i), gaussian(2,i)/norma_vector(2)
     end do
     close(102)
     close(101)
  end do

  


  close(103)

end subroutine bond_evolution

subroutine addposcar(option_generate_files_only,option_filepath,stage, prev_structures_overwrite)
  type(unitcell), dimension(:), allocatable :: formula
  character(1024) :: tmp, name,option_filepath
  character(3), dimension(:), allocatable :: elnames
  double precision :: cellmultiplier
  double precision, dimension(:), allocatable :: bondcutoff
  integer :: q,l,k,j,i,structno, ecount, eltot, leng,structures, stage
  integer :: prev_structures, option_generate_files_only, prev_structures_overwrite
  type (atom), dimension(:,:), allocatable :: tmplist,atomlist, alistrep
  integer, dimension(:), allocatable :: elno
  !! Wipes the randomly generated formula
  structures=1
  structno=1

  if (option_generate_files_only.eq.0) then
     prev_structures=structurecounter("pos")
     call touchpos()
     call touchposdir(structures,prev_structures)
  end if

  print*, structures, prev_structures, "!!"
  allocate(formula(structno))

  write(6,*) "Please enter the filename you wish to add to the database"
  !read(*, *) name
  if (option_generate_files_only.eq.0) then
     write(name,'(A)') "POSCAR"
  else 
     print*, "here"
     write(name,'(A,A)') trim(adjustl(option_filepath)),"/POSCAR"
  end if
  print*, name
  open(50, file=name)
  print*, "step 1"
  read(50, '(A)') tmp
  read(50, '(F16.0)') cellmultiplier


  do i=1, 3
     read(50,*) formula(1)%cell(i,1),formula(1)%cell(i,2),formula(1)%cell(i,3)
     formula(1)%cell(i,1)=formula(1)%cell(i,1)*cellmultiplier
     formula(1)%cell(i,2)=formula(1)%cell(i,2)*cellmultiplier
     formula(1)%cell(i,3)=formula(1)%cell(i,3)*cellmultiplier
     print*, formula(1)%cell(i,1),formula(1)%cell(i,2),formula(1)%cell(i,3)

  end do




  if (option_generate_files_only.eq.0) then

     write(name,'(A11,I0.3,A7)')"pos/POSCAR_",structures+prev_structures,"/POSCAR"
     open(structures+10000, file=name,status="new")
     write(structures+10000,*) "Test"
     write(structures+10000,*) 1.0
     do i=1, 3
        write(structures+10000,*) formula(1)%cell(i,1), &
             &formula(1)%cell(i,2), formula(1)%cell(i,3)
     end do
  end if
  read(50,'(A)') tmp


  ecount=0
  eltot=0
  allocate(tmplist(structno,1000))
  allocate(elnames(1024))
  do i=1,len(tmp)-1
     if(i.eq.1) then 
        if((scan(tmp(i:i+1)," ").eq.0).or.&
             &((scan(tmp(i:i+1)," ").eq.2))) then
           eltot=eltot+1
           !print*, tmp(i:i+1)
           if(scan(tmp(i:i+1)," ").eq.0) then
              elnames(eltot)=tmp(i:i+1)
           end if
           if((scan(tmp(i:i+1)," ").eq.2).and.(scan(tmp(i-1:i)," ").eq.1)) then
              elnames(eltot)=tmp(i:i)
           end if
        else
        end if


     else
        if((scan(tmp(i:i+1)," ").eq.0).or.&
             &((scan(tmp(i:i+1)," ").eq.2).and.(scan(tmp(i-1:i)," ")&
             &.eq.1))) then
           eltot=eltot+1
           !print*, tmp(i:i+1)
           if(scan(tmp(i:i+1)," ").eq.0) then
              elnames(eltot)=tmp(i:i+1)
           end if
           if((scan(tmp(i:i+1)," ").eq.2).and.(scan(tmp(i-1:i)," ").eq.1)) then
              elnames(eltot)=tmp(i:i)
           end if
        else
        end if
     end if
  end do
  print*, elnames(:), eltot
  allocate(bondcutoff(eltot))
  allocate(elno(eltot))
  !read(50,'(4X)', advance='no')
  k=0
  read(50,*) elno
  print*, elno(:)
  read(50, *) tmp
  k=0
  do i=1, eltot
     do j=1, elno(i)
        k=k+1

        tmplist(structures,k)%name=elnames(i)
     end do
  end do
  print*, k
  l=0
  do i=1, eltot
     do j=1, elno(i)
        l=l+1
        read(50, *)&
             &tmplist(structures,l)%position(1),tmplist(structures,l)%position(2),&
             &tmplist(structures,l)%position(3)
        print*, "Warning! Direct POSCAR import ONLY"
        print*, tmplist(structures,l)%position(1),tmplist(structures,l)%position(2),&
             &tmplist(structures,l)%position(3)
        wait(5)
        do q=1,3
           tmplist(structures,L)%position(q)=&
                &formula(structures)%cell(1,q)*tmplist(structures,L)%position(1)+&
                &formula(structures)%cell(2,q)*tmplist(structures,L)%position(2)+&
                &formula(structures)%cell(3,q)*tmplist(structures,L)%position(3)
           print*, tmplist(structures,L)%position(q)
        end do
        !print*, atomlist(structures,j)%position(:)
        !print*, atomlist(1,j)%position(1),atomlist(1,j)%position(2),atomlist(1,j)%position(3\

     end do
  end do
  allocate(atomlist(1,k))
  do i=1, k
     atomlist(structures,i)%position=tmplist(structures,i)%position
     atomlist(structures,i)%name=tmplist(structures,i)%name
     !print*, tmplist(structures,i)%position(:)

  end do
  deallocate(tmplist)
  if (option_generate_files_only.eq.0) then
     call poswrite(formula(structures)%cell,atomlist,k,1,1,prev_structures)
  end if
  allocate(alistrep(1,k*27))
  do i=1, k
     call atomrepeater(structures,atomlist(structures,i)%position,atomlist(structures,i)%name,&
          &alistrep,&
          &formula,i,k)
  end do
!!! The problem is here. Firstly, you need to pass STAGE to this function so it can be correctly passed to the following functions 
!!! Secondly, prev_structures shouldn't be calculated but also passed from the previous function, as loop number from calling loop
!!!
  if (prev_structures_overwrite.eq.0) then 
     prev_structures=structurecounter("bon")
  else  
     prev_structures=prev_structures_overwrite-1
  end if

  !do i=1, k
  !   call generatebondfiles(stage,structures,prev_structures,atomlist,alistrep,eltot,elno,i,elnames)
  !end do
  bondcutoff(1)=1.0
  print*, "MANUALLY EDIT"
  stop
  if (prev_structures_overwrite.eq.0) prev_structures=structurecounter("bad")
  !do i=1, k
  !   call generateanglefiles(stage,structures,prev_structures,atomlist,alistrep,eltot,elno,i,bondcutoff)
  !end do

  !do i=1, k
  !   call generate4files(stage,structures,prev_structures,atomlist,alistrep,eltot,elno,i,bondcutoff)
  !end do

  !do i=1, k 
  !   call generate4files(stage,structures,prev_structures,atomlist,alistrep,eltot,elno,i,bondcutoff)
  !end do

  write(tmp,'(A11,I0.3)')"pos/POSCAR_",structures+prev_structures
  call Incarwrite(trim(adjustl(tmp)),500, 20*k)
  call Jobwrite(tmp,1,1,1)
  call potwrite(tmp, elnames, eltot)




end subroutine addposcar


subroutine touchposdir(structures,prev_structures)
  character(1024) :: name, tmp, command  
  logical :: dir_e
  integer :: structures, prev_structures 
  write(name,'(A11,I0.3,A7)')"pos/POSCAR_",structures+prev_structures,"/POSCAR"
  !!Calculates the new structure number, and writes it to tmp
  write(tmp,'(A11,I0.3)')"pos/POSCAR_",structures+prev_structures
  !!Checks if a directory to contain that file exsts already (It should never exist, could add warning)
  inquire(file=tmp, exist=dir_e)
  if(dir_e) then
  else
     !!Writes a command to create said directory
     write(command,'(A17,I0.3)')"mkdir pos/POSCAR_",structures+prev_structures
     CALL execute_command_line(trim(adjustl(command)))
  end if
end subroutine touchposdir

subroutine angledistribution(atoma,resul)
  integer :: i,j,k,length
  character(3) :: atoma 
  character(1024) :: filename, buffer 
  logical :: dir_e 
  double precision, dimension(:), allocatable :: resul

  if(allocated(resul)) deallocate(resul)
  write(filename,*) "angledist", trim(adjustl(atoma)),".in"
  filename=trim(adjustL(filename))
  inquire(file=trim(adjustl(filename)), exist=dir_e)
  if(dir_e) then 
     !print*, "Found distribution for ", atoma
  else 
     print*, "STOP. YOU VIOLATED THE LAW. YOU MUST HAVE DISTRIBUTIONS FOR", atoma
     stop 
  end if
  open(50,file=filename, status="old")
  j=0
  read(50,"(1024A)") buffer
  close(50)
  k=0
  length=0
  loop1: do i=1, 1024 
     if(i.lt.k) cycle
     if((scan(buffer(i:i)," ").eq.1)) cycle  
     loop2: do j=0, 1024-i
        if((scan(buffer(i:i+j)," ").gt.0)) then 
           k=i+j
           length=length+1
           exit loop2
        end if

     end do loop2
  end do loop1
  allocate(resul(length))
  length=0
  j=0
  k=0
  loop1p: do i=1, 1024
     if(i.lt.k) cycle
     if((scan(buffer(i:i)," ").eq.1)) cycle
     loop2p: do j=0, 1024-i
        if((scan(buffer(i:i+j)," ").gt.0)) then
           length=length+1
           read(buffer(i:i+j-1),*) resul(length)

           k=i+j
           exit loop2p
        end if

     end do loop2p
  end do loop1p

end subroutine angledistribution



recursive subroutine invar(a,b,c)   
  implicit none 
  integer :: a,i, tmpvar, iostat, ierr
  character(1024) :: buffer, command
  integer, dimension(:), allocatable :: b
  character(1024), dimension(:), allocatable, intent(out) :: c 


  open(71, file="Infile.txt") 
  rewind(71)
  do i=1, a 
     read(71,*) buffer
  end do
  close(71) 
  ! print*, buffer
  write(command,'(4A,X,A,X,A,X,A)') "(echo",' "',trim(adjustl(buffer)),'" | ','sed -e "s/^.*=//g")>tmp.txt'
  CALL execute_command_line(command)
  !print*, command
  open(72,file="tmp.txt")

  if(a.le.3) then
     allocate(b(1))
     read(72,*) b
     close(72)
  else if((a.eq.4)) then 
     close(72)
     CALL invar(2,b,c)

     open(71, file="Infile.txt")
     rewind(71)
     do i=1, a
        read(71,'(A,X,A,X,A,X,A,X,A,X,A,X,A,X,A)') buffer
     end do
     close(71)
     ! print*, buffer
     write(command,'(4A,X,A,X,A,X,A)') "(echo",' "',trim(adjustl(buffer)),'" | ','sed -e "s/^.*=//g")>tmp.txt'
     CALL execute_command_line(command)
     ! print*, command
     open(72,file="tmp.txt")

     tmpvar=b(1)
     deallocate(b)
     allocate(c(tmpvar))
     read(72,*) c
     close(72)

  else if((a.eq.5)) then
     close(72)
     CALL invar(2,b,c)

     open(71, file="Infile.txt")
     rewind(71)
     do i=1, a
        read(71,'(A,X,A,X,A,X,A,X,A,X,A,X,A,X,A)') buffer
     end do
     close(71)
     ! print*, buffer
     write(command,'(4A,X,A,X,A,X,A)') "(echo",' "',trim(adjustl(buffer)),'" | ','sed -e "s/^.*=//g")>tmp.txt'
     CALL execute_command_line(command)
     ! print*, command
     open(72,file="tmp.txt")

     tmpvar=b(1)
     deallocate(b)
     allocate(b(tmpvar))
     read(72,*) b
     !print*, b
     close(72)
  else if (a.eq.13) then 
     open(71, file="Infile.txt")
     rewind(71)
     do i=1, a
        read(71,'(A,X,A,X,A,X,A,X,A,X,A,X,A,X,A)') buffer
     end do
     close(71)
     !print*, buffer
     write(command,'(4A,X,A,X,A,X,A)') "(echo",' "',trim(adjustl(buffer)),'" | ','sed -e "s/^.*=//g")>tmp.txt'
     CALL execute_command_line(command)
     !print*, command
     open(72,file="tmp.txt")


     allocate(c(1))
     read(72,*) c
     print*, c
     close(72)
  else if(a.eq.14) then    
     close(72) 
     write(command,*) "sed -n 's|buildmap_resolution=||p'<Infile.txt>tmp.txt"
     CALL execute_command_line(command)
     print*, command
     open(72,file="tmp.txt")
     allocate(b(3))
     read(72,*) b
     close(72)
  else if(a.eq.16) then 
     close(72)
     write(command,*) "sed -n 's|VoidPseudoScanratio=||p'<Infile.txt>tmp.txt"
     Call execute_command_line(command)
     print*, command
     open(72,file="tmp.txt")
     allocate(b(3))
     read(72,*) b
     close(72)
     
  else
     allocate(b(1))
     read(72,*, iostat=iostat) b
     close(72) 
  end if


  write(command,*) "rm tmp.txt"
  CALL execute_command_line(command)
end subroutine invar
!!!-------------------------------------------------------------!!!
!!!This function provides the cross product of two vectors      !!!
!!!-------------------------------------------------------------!!!

function cross(a,b) result(axb)
  implicit none
  integer,parameter :: wp=selected_real_kind(15, 307) !double precision
  double precision,dimension(3) :: axb
  double precision,dimension(3),intent(in) :: a
  double precision,dimension(3),intent(in) :: b 
  axb(1) = a(2)*b(3) - a(3)*b(2)
  axb(2) = a(3)*b(1) - a(1)*b(3)
  axb(3) = a(1)*b(2) - a(2)*b(1)
end function cross

!!!---------------------------------------------------------------------------!!!
!!!This function writes a list of atomic positions to the POSCAR in cartesian !!!
!!!---------------------------------------------------------------------------!!!

subroutine  poswrite(box,atomlist,len, structures, structno, prev_structures)

  integer :: i,j, len, reason, structures, structno, prev_structures
  type(atom), dimension(:,:) :: atomlist
  double precision, dimension(3,3) :: box
  do i=1, len 
     if(i.eq.len) then 
        write(structures+10000,'(5X,A2)', IOStat=reason, advance="no")&
             &atomlist(structures,i)%name
     else 
        !if(i.eq.(len-1)) cycle
        if (atomlist(structures,i)%name.eq.atomlist(structures,i+1)%name) cycle
        write(structures+10000,'(5X,A2)', IOStat=reason, advance="no") &
             &atomlist(structures,i)%name      
     end if
  end do

  write(structures+10000, '(2X)')
  !write(structures+10000,*) "HERE"
  j=1
  do i=1, len 
     if(i.eq.len) then 
        write(structures+10000,'(1X,I3)', IOStat=reason, advance="no") j
     else
        if (atomlist(structures,i)%name.eq.atomlist(structures,i+1)%name) then 
           j=j+1
           cycle
        else 
           write(structures+10000,'(1X,I3)', IOStat=reason, advance="no") j
           j=1

        end if
     end if
  end do


  write(structures+10000,'(/,A)') "Cartesian"

  do i=1, len 
     write(structures+10000,'(1X,F0.16,1X,F0.16,1X,F0.16)') atomlist(structures,i)%position(:)
  end do

end subroutine poswrite

!!!-------------------------------------------------!!!
!!!This function returns normvol as the cell volume !!!
!!!-------------------------------------------------!!!

function sphereoverlap(r,rp,b, pi) result(volume) 
  double precision :: r,rp,volume, b, step1, step2, pi, step3
  step1=pi*(r+rp-b)**2
  step2=(b**2)+(2*b*rp)-(3*(rp**2))+(2*b*r)+(6*rp*r)-3*(r**2)
  step3=1.0/(12.0*b)
  !  print*, "Sub.f90 ->" ,step1*step2*step3
  volume=step1*step2*step3
  if(volume.lt.0) volume=0
end function sphereoverlap

function  cellvol(x) result(normvol)
  double precision, dimension(3,3) :: x
  double precision, dimension(3) :: a,b,c, tmp
  integer :: i
  double precision :: normvol
  a=x(:,1)
  b=x(:,2)
  c=x(:,3)
  tmp=cross(b,c)
  normvol=0
  do i=1,3
     normvol=normvol+a(i)*tmp(i)
  end do
  !print*, "Cell volume is", abs(normvol)
end function cellvol


!!!Bondlength!!!
function bondlength(A,B) result(C)
  double precision, dimension(3), intent(in) :: A,B
  double precision :: C
  C=(((A(1)-B(1))**2)+((A(2)-B(2))**2)+((A(3)-B(3))**2))**0.5
end function bondlength

function structurecounter(dir_t)  result(n)
  logical :: file_exists
  character(1024) :: name, dir_name, dir_append
  character(3) :: dir_t
  integer :: n, gap_test
  n=1
  mainloop : do while(n.ne.0)
     if(dir_t.eq."don") dir_name="/DON_"     
     if(dir_t.eq."pos") dir_name="/POSCAR_"
     if(dir_t.eq."don") dir_append="/DON"
     if(dir_t.eq."pos") dir_append="/POSCAR"
     if(dir_t.eq."bon") dir_name="/BON_"
     if(dir_t.eq."bon") dir_append=""
     if(dir_t.eq."bad") dir_name="/BAD_"
     if(dir_t.eq."bad") dir_append=""
     write(name,'(A3,A,I0.3,A7)') dir_t,trim(adjustl(dir_name)), n,dir_append
     name=trim(adjustl(name))

     INQUIRE(FILE=name, EXIST=file_exists)
     if(file_exists) then 
        n=n+1
        cycle 
     else 
        do gap_test=1, 1000
           write(name,'(A3,A,I0.3,A7)') dir_t,trim(adjustl(dir_name)), n+gap_test,dir_append
           INQUIRE(FILE=name, EXIST=file_exists)
           if(file_exists) then
              n=n+gap_test
              cycle mainloop 
           end if
        end do
        n=n-1
        





        exit
     end if
  end do mainloop
end function structurecounter

function angle_between_vectors(A,B) result  (theta) 
  double precision :: theta
  double precision, dimension(3) :: A,B
  !if(isnan(A))
  theta=Acos((dot_product(A,B))/(norm2(A)*norm2(B)))
end function angle_between_vectors


function bondangle(A,B,C) result(theta)
  double precision :: theta, x
  double precision, dimension(3) :: A,B,C, bond1, bond2
  integer :: i

  theta=Acos((dot_product(A-B,C-B))/(norm2(A-B)*norm2(C-B)))
  if(isnan(theta)) then 
     theta=0
  end if
end function bondangle

function fourbody(A,B,C,D) result(theta)
  double precision :: theta, x
  double precision, dimension(3) :: A,B,C,D, bond1, bond2, AB, AC, AD
  !print*, ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
  !print*, A, B, C, D
  AB(:)=B(:)-A(:)
  AC(:)=C(:)-B(:)
  AD(:)=D(:)-B(:)
  !print*, AB, AC, AD
  !print*, "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"

  theta=angle_between_vectors(cross(AB,AC),AD)
  !print*, theta*180/3.141592
  
  !print*, cross(AB,AC)!
  !Theta set to 1000 so contribution function can recognise the special case
  if(isnan(theta)) theta=1000
  !if(theta.ne.0) print*, theta
end function fourbody






subroutine Incarwrite(filepath,nstep,bandno) 
  character(1024) ::  name
  integer :: nstep, bandno
  character(*) :: filepath

  write(name,'(A,A6)') trim(filepath), "/INCAR"
  print*, "The name of the new file will be", trim(adjustl(name))

  !print*, trim(adjustl(name))
  open(unit=11,file=trim(adjustl(name)))!"pos/POSCAR_001")!trim(name))

  write(11, *)"SYSTEM = RSS"
  write(11, *)"### sys ###"
  write(11, *)"GGA = PE "
  write(11, *)"ISYM = 0"
  write(11, *)"ENCUT = 500"
  write(11, *)"ENAUG = 500"
  write(11, *)"PREC = Accurate"
  write(11, *)"ISTART = 1"
  write(11, *)"ICHARG = 2"
  write(11, *)"LWAVE = .FALSE."
  write(11, *)"LAECHG =.FALSE."
  write(11, *)"LMAXMIX = 6"
  write(11, *)"LASPH = .TRUE."
  write(11, *)"LMIXTAU = .TRUE."
  write(11, *)"LREAL = .FALSE."
  write(11, *)"NWRITE = 3"
  write(11, *)"ALGO = N"
  write(11, *)"LHAR = .FALSE."
  write(11, *)"LVTOT = .FALSE."
  write(11, *)

  write(11, *)"### vdw ###"
  write(11, *)"#IVDW = 11"
  write(11, *)"#VDW_RADIUS = 80"
  write(11, *)"#VDW_CNRADIUS = 50.0"
  write(11, *)"#VDW_S8 = 0.722"
  write(11, *)"#VDW_SR = 1.217"
  write(11, *) 

  write(11, *)"### elc ###"
  write(11, *)"#AMIX = 0.6"
  write(name,'(I0.3)') nstep
  write(11, *)"NELM = ", adjustL(trim(name))
  write(11, *)"NELMIN = 5"
  write(11, *)"NELMDL = -5"
  write(11, '(1X,A,1X,I0.3)') "NBANDS =", bandno
  write(11, *)"EDIFF = 10d-8"
  write(11, *)"ISMEAR = 0"

!!! THIS SHOULD BE CHANGED BASED ON INTUITION

  write(11, *)"SIGMA = 0.2"
  write(11, *) 

  write(11, *)"### Mag ###"
  write(11, *)"#ISPIN = 2"
  write(11, *)"#MAGMOM = 1 -1 1 1 -1 -1"
  write(11, *) 

  write(11, *)"### ncl ###"
  write(11, *)"#LNONCOLLINEAR = .TRUE."
  write(11, *)"#LSORBIT = .TRUE."
  write(11, *)"#GGA_COMPAT = .FALSE."
  write(11, *) 

  write(11, *)"### MPI ###"
  write(11, *)"#NCORE = 24"
  write(11, *)"#KPAR = 2"
  write(11, *) 

  write(11, *)"### Rlx ###"
  write(11, *)"#LMAXPAW =-1"
  write(11, *)"ADDGRID = .TRUE."
  write(11, *)"#POTIM = 0.1"
  write(11, *)"#NFREE = 15"
  write(11, *)"#NSW = 150"
  write(11, *)"#ISIF = 2"
  write(11, *)"#IBRION = 1"
  write(11, *)"#EDIFFG = -0.001"
  write(11, *) 

  write(11, *)"### dos ###"
  write(11, *)"#EMIN = -3.0"
  write(11, *)"#EMAX =  6.0"
  write(11, *)"#NEDOS = 5000"
  write(11, *)"#LORBIT = 11"

  close(11)
end subroutine Incarwrite

subroutine Jobwrite(filepath, a,b,c) 

  integer :: a,b,c
  character(1024) ::  name
  character(20) :: filepath
  name=" "
  !print*, filepath
  write(name,*) "cp job_vasp_isca.in ", trim(filepath)
  call execute_command_line(name)
  write(name,'(A,A8)') trim(filepath), "/KPOINTS"
  open(unit=11,file=trim(adjustl(name)), status='new')!"pos/POSCAR_001")!trim(name))
  write(11, '(A7)') "KPOINTS"
  write(11, '(A1)') "0" 
  write(11, '(A1)') "G" 
  write(11, '(I0,X,I0,X,I0)') a, b, c
  write(11, '(I0,X,I0,X,I0)') 0,0,0
  close(11) 
end subroutine Jobwrite





subroutine potwrite(filepath, elnames, eltot) 

  character(3), dimension(:), allocatable :: elnames
  character(1024) ::  name
  character(20) :: filepath
  integer :: eltot, i
  name=" "
  !print*, filepath
  !print*, eltot
  do i=1, eltot 

     !print*, i
     if(i.eq.1) then
        !print*, "I reached here"
        write(name, '(A3,1X,A6,A2,1X,A1,A,A)') "cat", "potcar", elnames(i),">",trim(adjustl(filepath)),"/POTCAR"
        !print*, name
        call execute_command_line(name)
     else       
        write(name, '(A,1X,A,A,1X,A,A,1X,A,1X,A,A)') "cat ",trim(adjustl(filepath)),"/POTCAR", "potcar", elnames(i)," >>",&
             & trim(adjustl(filepath)),"/POTCAR1"   

        call execute_command_line(name)
        write(name,'(A,1X,A,A8,1X,A,A)') "mv ",trim(adjustl(filepath)),"/POTCAR1", trim(adjustl(filepath)),"/POTCAR"
        call execute_command_line(name)
     end if

  end do
end subroutine potwrite

subroutine initialisehost(leng)
  type(unitcell), dimension(:), allocatable :: formula
  character(1024) :: tmp, name, location
  character(1024), dimension(:), allocatable :: elnames,tmpelnames
  double precision :: cellmultiplier,meanvol
  integer :: l,k,j,i,structno, ecount,addedelements, eltot, leng,structures, prev_structures
  type (atom), dimension(:,:), allocatable :: tmplist,atomlist,tmplist2
  integer, dimension(:), allocatable :: stochio, tmpstochio, ts2, tmpstochiotot
  integer, dimension(:), allocatable :: tmpdig
  character(1024), dimension(:), allocatable :: tmpels


  call invar(13,tmpdig,tmpels)
  write(name, '(A,A)') "POSCAR_",trim(adjustl(tmpels(1)))
  open(50, file=trim(adjustl(name)))
  do i=1, 5
     read(50, * )
  end do
  deallocate(tmpels)
  read(50,'(A)') tmp
  eltot=0
  allocate(tmpelnames(100))
  allocate(tmplist(1,1000))
  do i=1,len(tmp)-1


     if(i.eq.1) then 
        if((scan(tmp(i:i+1)," ").eq.0).or.&
             &((scan(tmp(i:i+1)," ").eq.2))) then
           eltot=eltot+1
           if(scan(tmp(i:i+1)," ").eq.0) then
              tmpelnames(eltot)=tmp(i:i+1)
           end if

           if((scan(tmp(i:i+1)," ").eq.2)) then
              tmpelnames(eltot)=tmp(i:i)
           end if
        end if
     else
        if((scan(tmp(i:i+1)," ").eq.0).or.&
             &((scan(tmp(i:i+1)," ").eq.2).and.(scan(tmp(i-1:i)," ")&
             &.eq.1))) then

           eltot=eltot+1
           !print*, tmp(i:i+1)
           if(scan(tmp(i:i+1)," ").eq.0) then
              tmpelnames(eltot)=tmp(i:i+1)
           end if
           if((scan(tmp(i:i+1)," ").eq.2).and.(scan(tmp(i-1:i)," ").eq.1)) then
              tmpelnames(eltot)=tmp(i:i)
           end if
        else              
        end if
     end if
  end do
  allocate(tmpstochio(eltot))
  read(50,*) tmpstochio
  do i=1, eltot
     leng=leng+tmpstochio(i)
  end do
  close(50)
  
  
  !100 restart : do i=1, eltot-1
  !   if(tmpelnames(i).eq.tmpelnames(i+1)) then 
  !      deallocate(ts2)
  !      allocate(ts2(eltot-1))
  !      do j=1, eltot-1 
  !         ts2(j)=tmpstochio(j+1) 
  !         
  !      end do
  !      deallocate(tmpstochio)
  !      allocate(tmpstochio(eltot-1))
  !      tmpstochio=ts2
  !      eltot=eltot-1
  !      addedelements=
  !      GO TO 100  
  !     
  !      
  !   end if
  !end do restart
  
  !print*, eltot, tmpstochio, tmpelnames 
  

  !print*, "The new atom total is", leng
end subroutine initialisehost

subroutine generatebondfiles(stage,structures,prev_structures,array,repeatedarray,eltot,stochio,atomnumber,structure_elements)
  use atomtype
  implicit none
  character(1024) :: command,name,tmp, location_string, debug_string
  integer :: l,i,el,structures,prev_structures,stage, eltot, atomnumber, tmpint,m,j,k
  type (atom), dimension(:,:), allocatable :: array, repeatedarray
  integer, dimension(:), allocatable :: stochio
  logical :: dir_e
  double precision, dimension(:), allocatable :: list, list_tmp
  character(3), dimension(:), allocatable :: structure_elements

  if (stage.eq.1) location_string="bon"
  if (stage.eq.2) location_string="bon2"
  if (stage.eq.3) location_string="bon3"

  inquire(file=adjustl(trim(location_string)),exist=dir_e)
  if(dir_e) then
  else
     write(name,*) "mkdir ", trim(adjustl(location_string))
     call execute_command_line(name)
  end if
  !print*, structurecounter("bon"), "!!!!!!!!!!!!"
  !prev_structures=structurecounter("bon")
  write(tmp,'(A,A,I0.3)') trim(adjustl(location_string)),"/BON_",structures+prev_structures
  inquire(file=tmp, exist=dir_e)
  if(dir_e) then
  else
     write(command,'(A,A,A,I0.3)') "mkdir ",trim(adjustl(location_string)),"/BON_",structures+prev_structures
     call execute_command_line(command)
  end if

  i=1
  tmpint=atomnumber
  !print*, atomnumber,"-------------------------------"
  do while (i.ge.1)

     if((tmpint-stochio(i)).gt.0) then
        !print*, (tmpint-stochio(i)), tmpint, stochio(i)
        !print*, array(structures,atomnumber)%name, array(structures,tmpint-stochio(i))%name
        if(array(structures,atomnumber)%name.ne.(array(structures,tmpint-stochio(i))%name)) then 
           tmpint=tmpint-stochio(i)
        end if
        i=i+1
        !print*, tmpint+stochio(i),tmpint, i
        if(i.ge.size(stochio)) i=0
     else
        i=0
     end if
  end do
  write(name, '(A,A,I0.3,A1,A,A,I0.3)') trim(adjustl(location_string)),"/BON_",structures+prev_structures,"/",&
       &trim(adjustl(array(structures,atomnumber)%name)),"_",tmpint
  open(101, file=name)
  m=0
  k=0
  do el=1, eltot
     m=0
     allocate(list(27*stochio(el)))
     
     do i=1, stochio(el)
        do L=1, 27
           k=k+1
           if((L.eq.14).and.(bondlength(repeatedarray(structures,k)%position,&
                &array(structures,atomnumber)%position)).eq.0) cycle
           if(bondlength(repeatedarray(structures,k)%position,&
                &array(structures,atomnumber)%position).gt.5) cycle
              


           m=m+1
           list(m)=bondlength(repeatedarray(structures,k)%position,&
                &array(structures,atomnumber)%position)
           !write(101,*) list(m), m, "L"
        end do
     end do
     do i=1, m
        do j=1, m
           if(i.eq.j) cycle
           if(abs(list(i)-list(j)).lt.0.1) then
              list(j)=0
           end if
        end do
     end do
     write(101,*) structure_elements(el)
     
     do i=1, m
        if(list(i).lt.0.001) cycle

        write(101,*) list(i), i
    
     end do
     deallocate(list)


  end do
end subroutine generatebondfiles

subroutine generate4files(stage,structures,prev_structures,array,repeatedarray,&
     &eltot,stochio,atomnumber,bondlengthcutoff,structure_elements)
  character(1024) :: command,name,tmp,location_string
  integer :: stage,prev_structures,structures,eltot,atomnumber,tmpint,l,k,n,x,i,j,m,p,ii,jj,kk,LL, dim, tool, counter
  integer :: ip, bondint_1, bondint_2, bondint_3, bondint_4
  type (atom), dimension(:,:), allocatable :: array, repeatedarray
  integer, dimension(:), allocatable :: stochio
  double precision :: res
  double precision, dimension(:,:), allocatable :: bondlengthcutoff
  double precision, dimension(3) :: ab,ac,ad, t
  logical :: dir_e
  double precision, dimension(:), allocatable :: list, list_tmp
  character(3), dimension(:), allocatable :: structure_elements


  if (stage.eq.1) location_string="4body"
  if (stage.eq.2) location_string="4body2"
  if (stage.eq.3) location_string="4body3"
  inquire(file=location_string,exist=dir_e)
  if(dir_e) then
  else
     write(name,*) "mkdir ", trim(adjustl(location_string))
     call execute_command_line(name)
  end if
  write(tmp,'(A,A,I0.3)') trim(adjustl(location_string)),"/4BOD_",structures+prev_structures
  inquire(file=tmp, exist=dir_e)
  if(dir_e) then
  else
     write(command,'(A,A,A,I0.3)')"mkdir ", trim(adjustl(location_string)),"/4BOD_",structures+prev_structures
     call execute_command_line(command)
  end if
  i=1
  tmpint=atomnumber
  do while (i.ge.1)
     
     if((tmpint-stochio(i)).gt.0) then
        if(array(structures,atomnumber)%name.ne.(array(structures,tmpint-stochio(i))%name)) then
           tmpint=tmpint-stochio(i)
        end if
        
        i=i+1
        if(i.gt.size(stochio)) i=0
     else
        i=0
     end if
  end do

  write(name, '(A,A,I0.3,A1,A,A,I0.3)') trim(adjustl(location_string)),"/4BOD_",structures+prev_structures,"/",&
       &trim(adjustl(array(structures,atomnumber)%name)),"_",tmpint
  open(101, file=name)
  counter=0
  allocate(list(1))
  x=0
  
  do ip=1, size(structure_elements)
     if(trim(adjustl(structure_elements(ip))).ne.&
          &trim(adjustl(repeatedarray(structures,atomnumber)%name))) then
        cycle
     else
        bondint_1=ip
     end if
  end do
  




  do i=1, eltot
     do j=1, stochio(i)
        do p=1, 27







           x=x+1



           do ip=1, size(structure_elements)
              if(trim(adjustl(structure_elements(ip))).ne.&
                   &trim(adjustl(repeatedarray(structures,x)%name))) then
                 cycle
              else
                 bondint_2=ip
              end if
           end do



           if(x.eq.atomnumber) cycle
           if(bondlength(array(structures,atomnumber)%position,&
                &repeatedarray(structures,x)%position).gt.bondlengthcutoff(bondint_1,bondint_2)) cycle

           if(bondlength(array(structures,atomnumber)%position,&
                &repeatedarray(structures,x)%position).lt.0.1) cycle
           if((P.eq.14).and.(bondlength(repeatedarray(structures,x)%position,&
                &array(structures,atomnumber)%position)).eq.0) cycle




           m=0
           do k=1, eltot
              do L=1, stochio(k)
                 do n=1, 27
                    !print*, i,j,p,k,l,n
                    !Change from atomnumber to x
                    m=m+1


                    do ip=1, size(structure_elements)
                       if(trim(adjustl(structure_elements(ip))).ne.&
                            &trim(adjustl(repeatedarray(structures,m)%name))) then
                          cycle
                       else
                          bondint_3=ip
                       end if
                    end do
                    

                    if(m.eq.x) cycle
                    if(m.eq.atomnumber) cycle
                    if(bondlength(repeatedarray(structures,x)%position,&
                         &repeatedarray(structures,m)%position).lt.0.01) cycle
                    if(bondlength(repeatedarray(structures,x)%position,&
                         &repeatedarray(structures,m)%position).gt.bondlengthcutoff(bondint_2,bondint_3)) cycle
                    if((N.eq.14).and.(bondlength(repeatedarray(structures,m)%position,&
                         &repeatedarray(structures,x)%position)).eq.0) cycle
  
                    LL=0
                    do ii=1, eltot
                       do jj=1, stochio(ii)
                          do kk=1, 27
                             LL=LL+1
                             if(LL.eq.x) cycle
                             if(LL.eq.m) cycle
                             if(LL.eq.atomnumber) cycle
                             !change from atomnumber to x
                             
                             do ip=1, size(structure_elements)
                                if(trim(adjustl(structure_elements(ip))).ne.&
                                     &trim(adjustl(repeatedarray(structures,LL)%name))) then
                                   cycle
                                else
                                   bondint_4=ip
                                end if
                             end do





                             if(bondlength(repeatedarray(structures,x)%position,&
                                  &repeatedarray(structures,LL)%position).gt.bondlengthcutoff(bondint_2,bondint_4)) cycle
                             if((kk.eq.14).and.(bondlength(repeatedarray(structures,LL)%position,&
                                  &repeatedarray(structures,x)%position)).eq.0) cycle
                             ab=array(structures,atomnumber)%position-repeatedarray(structures,x)%position
                             ac=repeatedarray(structures,x)%position-repeatedarray(structures,m)%position
                             ad=repeatedarray(structures,x)%position-repeatedarray(structures,LL)%position
                             res=fourbody(array(structures,atomnumber)%position,&
                                  &repeatedarray(structures,x)%position,&
                                  &repeatedarray(structures,m)%position,&
                                  &repeatedarray(structures,LL)%position)
                             
                             
                             if(isnan(res)) cycle
                             if(res.eq.1000) cycle
                             if(res.eq.0) cycle
                             counter=counter+1
                             !print*, counter, size(list)
                             if(counter.ne.1) then
                                allocate(list_tmp(counter-1))
                                do tool=1, counter-1
                                   list_tmp(tool)=list(tool)
                                end do
                                deallocate(list)
                                allocate(list(counter))
                                do tool=1, counter-1
                                   list(tool)=list_tmp(tool)
                                end do
                                list(counter)=res
                                deallocate(list_tmp)
                             else 
                                list(1)=res
                             end if
                             !*bondlength(array(structures,atomnumber)%position,&
                             !                                 &repeatedarray(structures,LL)%position)
                             
                          end do
                       end do
                    end do


                    !write(101,*) res

                 end do

              end do
           end do
        end do
     end do
  end do
  !do i=1, size(list)
     !print*, list(i), i
  !end do
  
  do i=1, counter
     do j=1, counter
        if(i.eq.j) cycle 
        if(abs(list(i)-list(j)).lt.0.01) then 
           list(j)=0
        end if



        
     end do
  end do

     do i=1, counter
        !print*, list(i)
        if(list(i).lt.0.001) cycle
        
        if(list(i).gt.3.141592654/2) then 
           write(101,*) list(i)
        else
           write(101,*) list(i)
        end if
     end do
     

  close(101)
end subroutine generate4files




subroutine generateanglefiles(stage,structures,prev_structures,array,repeatedarray,eltot,stochio,&
     &atomnumber,bondlengthcutoff,structure_elements)
  character(1024) :: command,name,tmp,location_string, name2
  integer :: stage,prev_structures,structures,eltot,atomnumber,tmpint,l,k,n,x,i,j,m,p,ip, bondint_1&
       &,bondint_2,bondint_3
  type (atom), dimension(:,:), allocatable :: array, repeatedarray
  integer, dimension(:), allocatable :: stochio
  double precision, dimension(:,:), allocatable :: bondlengthcutoff
  logical :: dir_e
  character(3), dimension(:), allocatable :: structure_elements


  if (stage.eq.1) location_string="bad"
  if (stage.eq.2) location_string="bad2"
  if (stage.eq.3) location_string="bad3"
inquire(file=location_string,exist=dir_e)
if(dir_e) then
else
   write(name,*) "mkdir ", trim(adjustl(location_string))
   call execute_command_line(name)
end if
write(tmp,'(A,A,I0.3)') trim(adjustl(location_string)),"/BAD_",structures+prev_structures
inquire(file=tmp, exist=dir_e)
if(dir_e) then
else
   write(command,'(A,A,A,I0.3)')"mkdir ", trim(adjustl(location_string)),"/BAD_",structures+prev_structures
   call execute_command_line(command)
end if
i=1
tmpint=atomnumber
do while (i.ge.1)
   print*, tmpint, stochio(i)
   if((tmpint-stochio(i)).gt.0) then
      if(array(structures,atomnumber)%name.ne.(array(structures,tmpint-stochio(i))%name)) then
         tmpint=tmpint-stochio(i)
      end if
      i=i+1
      if(i.gt.size(stochio)) i=0
   else
      i=0
   end if
end do
print*, tmpint, trim(adjustl(array(structures,atomnumber)%name))
write(name, '(A,A,I0.3,A1,A,A,I0.3)') trim(adjustl(location_string)),"/BAD_",structures+prev_structures,"/",&
     &trim(adjustl(array(structures,atomnumber)%name)),"_",tmpint
open(101, file=name)




do ip=1, size(structure_elements)
   if(structure_elements(ip).ne.array(structures,atomnumber)%name) then
      cycle
   else
      bondint_1=ip
   end if
end do

x=0
do i=1, eltot
   do j=1, stochio(i)
      do p=1, 27
         do ip=1, size(structure_elements)
            if(trim(adjustl(structure_elements(ip))).ne.&
                 &trim(adjustl(repeatedarray(structures,x+1)%name))) then
               cycle
            else
               bondint_2=ip
            end if
         end do



         m=0
         x=x+1
         if(x.eq.atomnumber) cycle
         if(bondlength(array(structures,atomnumber)%position,&
              &repeatedarray(structures,x)%position).gt.bondlengthcutoff(bondint_1,bondint_2)) then 
            cycle
         end if
         if(bondlength(array(structures,atomnumber)%position,&
              &repeatedarray(structures,x)%position).lt.0.01) cycle

         do k=1, eltot
            do L=1, stochio(k)
               do n=1, 27
                  m=m+1
                  !print*, i,j,p,k,l,n
                  if(m.eq.x) cycle
                  if(m.eq.atomnumber) cycle
                  
                  do ip=1, size(structure_elements)
                     if(trim(adjustl(structure_elements(ip))).ne.&
                          &trim(adjustl(repeatedarray(structures,m)%name))) then
                         cycle
                     else
                        bondint_3=ip
                     end if
                  end do



                  if(bondlength(array(structures,atomnumber)%position,&
                       &repeatedarray(structures,m)%position).lt.0.01) cycle
                  if(bondlength(array(structures,atomnumber)%position,&
                       &repeatedarray(structures,m)%position).gt.bondlengthcutoff(bondint_1,bondint_3)) then 
                     cycle
                  end if
                  !if((bondangle(repeatedarray(structures,m)%position,&
                  !     &array(structures,atomnumber)%position,&
                  !     &repeatedarray(structures,x)%position).gt.0.1)&
                  !     &.AND.(bondangle(repeatedarray(structures,m)%position,&
                  !     &array(structures,atomnumber)%position,&
                  !     &repeatedarray(structures,x)%position).lt.3.1)) then
                     
                  if(bondangle(repeatedarray(structures,m)%position,&
                          &array(structures,atomnumber)%position,&
                          &repeatedarray(structures,x)%position).ge.3.141592654/2) then 
                     
                     write(101,*) bondangle(repeatedarray(structures,m)%position,&
                          &array(structures,atomnumber)%position,&
                          &repeatedarray(structures,x)%position)
                     
                     write(name2, '(A,A,I0.3,A1,A,A,A,A,A,A,I0.3)') &
                          &trim(adjustl(location_string)),"/BAD_",structures+prev_structures,"/",&
                          &trim(adjustl(structure_elements(bondint_2))),"_",&
                          &trim(adjustl(array(structures,atomnumber)%name)),"_",&
                          &trim(adjustl(structure_elements(bondint_3)))
                     open(102,file=name2, position="append")
                     write(102,*) 3.141592654 - bondangle(repeatedarray(structures,m)%position,&
                          &array(structures,atomnumber)%position,&
                          &repeatedarray(structures,x)%position)

                  else

                  write(101,*) bondangle(repeatedarray(structures,m)%position,&
                          &array(structures,atomnumber)%position,&
                          &repeatedarray(structures,x)%position)
                     
                     write(name2, '(A,A,I0.3,A1,A,A,A,A,A,A,I0.3)') &
                          &trim(adjustl(location_string)),"/BAD_",structures+prev_structures,"/",&
                          &trim(adjustl(structure_elements(bondint_2))),"_",&
                          &trim(adjustl(array(structures,atomnumber)%name)),"_",&
                          &trim(adjustl(structure_elements(bondint_3)))
                     open(102,file=name2, position="append")
                     write(102,*) bondangle(repeatedarray(structures,m)%position,&
                          &array(structures,atomnumber)%position,&
                          &repeatedarray(structures,x)%position)
                  end if


                    ! print*, bondlengthcutoff(bondint_1,bondint_2), structure_elements(bondint_1), structure_elements(bondint_2)
                    ! print*, bondlengthcutoff(bondint_1,bondint_3), structure_elements(bondint_1), structure_elements(bondint_3)


                     !end if
                  !if((bondangle(repeatedarray(structures,m)%position,&
                  !     &array(structures,atomnumber)%position,&
                  !     &repeatedarray(structures,x)%position).gt.3.0)) then 
                  !   print*, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
                  !   print*, repeatedarray(structures,m)%position
                  !   print*, array(structures,atomnumber)%position
                  !   print*, repeatedarray(structures,x)%position
                  !   print*, bondangle(repeatedarray(structures,m)%position,&
                  !     &array(structures,atomnumber)%position,&
                  !     &repeatedarray(structures,x)%position)
                  !   print*, array(structures,atomnumber)%position-repeatedarray(structures,m)%position
                  !   print*, array(structures,atomnumber)%position-repeatedarray(structures,x)%position
                  !   stop 
                  !end if
               end do
               
            end do
         end do
      end do
   end do
end do
print*, "ANGLEFINISHED"
close(101)
end subroutine generateanglefiles


subroutine atomrepeater(structures,position,element,array,unit,atomnumber,length)
use atomtype
implicit none
type (atom), dimension(:,:), allocatable :: array
double precision, dimension(3) :: position
integer :: j,atomnumber,length,x,y,z,structures,m
type(unitcell), dimension(:), allocatable :: unit
character(3) :: element
! print*, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
!print*,unit(structures)%cell, array(structures,1)%position
m=((atomnumber-1)*27)+1
do x=-1,1
   do y=-1,1
      do z=-1,1
         do j=1, 3
            !print*, array(structures,m)%position(j),position(j)
            array(structures,m)%position(j)=position(j)+&
                 &(x*unit(structures)%cell(1,j))+&
                 &(y*unit(structures)%cell(2,j))+&
                 &(z*unit(structures)%cell(3,j))


            !print*,   array(structures,m)%position(j)
         end do
         array(structures,m)%name=element


         m=m+1
      end do
   end do
end do
end subroutine atomrepeater

subroutine atomprojector(position,array,unit,atomnumber,structures)
use atomtype
implicit none
type (atom), dimension(:,:), allocatable :: array
double precision, dimension(3) :: position
integer :: j,atomnumber,length,x,y,z,structures,m
type(unitcell), dimension(:), allocatable :: unit
! print*, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
!print*,unit(structures)%cell, array(structures,1)%position
!m=((atomnumber-1)*27)+1
m=1



do x=-1,1
   do y=-1,1
      do z=-1,1
         if((x.eq.0).and.(y.eq.0).and.(z.eq.0)) cycle
         do j=1, 3
            

            array(1,m)%position(j)=position(j)+&
                 &(x*unit(structures)%cell(j,1))+&
                 &(y*unit(structures)%cell(j,2))+&
                 &(z*unit(structures)%cell(j,3))
           ! print*,   array(structures,m)%position(j)
         end do
         !print*, array(structures,m)%position(:)
         !stop
         m=m+1
      end do
   end do
end do

end subroutine atomprojector







end module help
