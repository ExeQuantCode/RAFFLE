module help
  use constants, only: real12, pi
  use geom, only: &
       get_bondlength, get_bondangle, get_dihedral_angle, &
       atomprojector
  use inputs, only: c_cut, c_min, filename_host
  use atomtype
  use vasp_file_handler, only: unitcell, structurecounter, &
       poscar_read, addposcar, &
       touchpos, touchposdir, &
       potwrite, jobwrite, incarwrite, poswrite, &
       atomrepeater
  use file_generator
  use contributions
  implicit none


  TYPE densitymatrix 
     real(real12), dimension(3) :: position 
     real(real12) :: density 
     integer :: checked
  end type densitymatrix


  TYPE nearmatrix
     real(real12), dimension(3) :: position 
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
    real(real12), dimension(2) :: Maximum_value 
    real(real12), dimension(:), allocatable :: read_in
    real(real12), dimension(:,:), allocatable :: Gaussian, tmp, comp_a, comp_b
    real(real12) :: Gaus_in, bondcut, sigma, similarity_index
    call execute_command_line("./similarity.sh")

!!!Inputs 
    bondcut=10  
    sigma=0.1
!!!
    error_calculator=0
    open(11,file="similarity_combinations.txt")
    read(11,*) loop_structures

    write(*,*) loop_structures
    allocate(database(loop_structures))
    allocate(total_atoms(loop_structures))
    do i=1, loop_structures-error_calculator

       read(11,*, END=99) assigned_index

       write(*,*) i, assigned_index, loop_structures-error_calculator

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
99  close(11)
    assigned_index=loop_structures
    error_calculator=assigned_index-i+1
    open(12,file="test")
    write(*,*) loop_structures-error_calculator
    do i=1, loop_structures-error_calculator
       do j=1, size(database(i)%stoichio,1)
          !write(*,*) i,j
          !write(*,*) database(i)%names(j), database(i)%stoichio(j)
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
    write(*,*) "HERE"
    do i=1, loop_structures-error_calculator
       do j=1, size(database(i)%stoichio,1)
          do k=1, database(i)%stoichio(j)
             do x=1, size(database(i)%stoichio,1)
                write(name,'(A,I0.3,A,A,A,I0.3,A,A,A,A,A)')"bon/BON_", i,"/",trim(adjustl(database(i)%names(j))),&
                     &"_",k,trim(adjustl(database(i)%names(j))),"_",trim(adjustl(database(i)%names(x))),"_bond_distribution"
                open(12,file=trim(adjustl(name)))

                do while(1.eq.1) 
                   read(12,*,IOSTAT=stat) Gaus_in
                   !write(*,*) Gaus_in
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
                write(*,*) trim(adjustl(name))

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
                !write(*,*) name 
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

                            !write(*,*) name
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
                            write(*,*) similarity_index, i,j,k,x,q,w,e,r
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
       &results_matrix,eltot, elnames,placed,num_VOID, append_matrix) 

    implicit none
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







subroutine generate_bondangle_distribution (tmpvector,atomlist,alistrep,structures,bondcutoff, &
 &atom_number_new, new_position_needed&
 &,agausssamp,distribution,eltot, sigma1, sigma2,elnames, predicted_positions)
character(3), dimension(:), allocatable :: elnames
real(real12), dimension(3) :: tmpvector
integer :: y, atom_number_new, structures, eltot,i,j,k,L, nlines 
logical :: new_position_needed, file_exists
type (atom), dimension(:,:), allocatable :: atomlist, alistrep, predicted_positions
type (nearmatrix), dimension(:), allocatable :: nearneighbourmatrix
real(real12) :: agausssamp, bondcutoff, bondangle_store, distribution, read_angle, sigma1, sigma2, step_significance
real(real12) :: tmpvalue, total_significance
character(1024) :: determined_anglefile
real(real12), allocatable, dimension(:) :: initial_angle_peaks

distribution=0
!if(allocated(nearneighbourmatrix)) deallocate(nearneighbourmatrix) 
allocate(nearneighbourmatrix(((atom_number_new-1)*27)+26))
k=0

do y=1, ((atom_number_new-1)*27)+26
 if(y.le.((atom_number_new-1)*27)) then 
    if(get_bondlength(tmpvector,alistrep(structures,y)%position).lt.bondcutoff) then
       k=k+1
       nearneighbourmatrix(k)%position=alistrep(structures,y)%position
    end if
 else
    if(get_bondlength(tmpvector,predicted_positions(1,y-(atom_number_new-1)*27)%position).lt.bondcutoff) then
       k=k+1
       nearneighbourmatrix(k)%position=predicted_positions(1,y-(atom_number_new-1)*27)%position
    end if
 end if
end do


!!! CHANGE alistrep etc to nearneighbourmatrix%position. SHould work.
total_significance=0
if(new_position_needed.neqv..TRUE.) then 
 do L=1, eltot
    !write(*,*) elnames(L), structures, atom_number_previous, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    !write(*,*) atomlist(structures,atom_number_previous+1)%name, "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"
    if(trim(adjustl(elnames(L))).ne.trim(adjustl(atomlist(structures,atom_number_new)%name))) cycle

    do j=1, k 
       do y=1, k 
          if (j.eq.y) cycle 
          bondangle_store=get_bondangle(nearneighbourmatrix(j)%position,tmpvector,&
               &nearneighbourmatrix(k)%position)
          if (get_bondlength(tmpvector,nearneighbourmatrix(k)%position).lt.1.0) then 
             new_position_needed=.TRUE.
             cycle
          end if
          if (get_bondlength(tmpvector,nearneighbourmatrix(k)%position).lt.1.0) then
             new_position_needed=.TRUE.
             cycle
          end if

          if(new_position_needed.eqv..TRUE.) cycle 

          call angledistribution(elnames(L),initial_angle_peaks)
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
10           CLOSE (11)
          end if
       end do
    end do
 end do
 distribution=abs(distribution/total_significance) 
end if


end subroutine generate_bondangle_distribution




subroutine generate_bondlength_distribution (sigma1,sigma2,atom_number_proposed,distribution,atomlist,alistrep,&
 &tmpvector,elrad,eltot,structures,elnames,new_position_needed, predicted_positions)
real(real12), dimension(:,:,:), allocatable :: elrad
integer :: nlines, atom_number_proposed, atom_number_previous, atom_total_repeated,j,k,l,eltot,structures
type (atom), dimension(:,:), allocatable :: atomlist, alistrep, predicted_positions
real(real12), dimension(3) :: tmpvector
real(real12) :: distribution, total_significance, step_significance, bondlength_store, sigma1,sigma2, read_bondlength 
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
          bondlength_store=get_bondlength(tmpvector,alistrep(structures,j)%position)
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
                !           write(*,*) distribution, nlines, total_significance 
                total_significance=total_significance+step_significance
                distribution=distribution+step_significance*(1.0/2.0)*(&
                     &abs(erf((bondlength_store-read_bondlength+0.5*sigma2)/((sigma1)*2**dble(1.0/2.0)))-&
                     &erf((bondlength_store-read_bondlength-0.5*sigma2)/((sigma1)*2**dble(1.0/2.0)))))
                nlines=nlines+1
             END DO
10           CLOSE (11)
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
       bondlength_store=get_bondlength(tmpvector,predicted_positions(1,j)%position)
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
             !           write(*,*) distribution, nlines, total_significance
             total_significance=total_significance+step_significance
             distribution=distribution+step_significance*(1.0/2.0)*(&
                  &abs(erf((bondlength_store-read_bondlength+0.5*sigma2)/((sigma1)*2**dble(1.0/2.0)))-&
                  &erf((bondlength_store-read_bondlength-0.5*sigma2)/((sigma1)*2**dble(1.0/2.0)))))
             nlines=nlines+1
          END DO
11        CLOSE (11)
       end if
    end do
 end do
 distribution=abs(distribution/total_significance)

 exit
end do

end subroutine generate_bondlength_distribution







subroutine regenerate_distribution_files (prev_structures)
implicit none
character(1024) :: command,name,tmp, location_string
integer :: l,i,structures,prev_structures,stage, eltot, atomnumber, tmpint,m,j,k, structure_loop
type (atom), dimension(:,:), allocatable :: array, repeatedarray
integer, dimension(:), allocatable :: stochio
logical :: dir_e, file_e, empty



character(1024) :: coordinate_type, header
character(LEN=:), allocatable :: read_in
character(1024) :: pathway, path_test
character(3), dimension(:), allocatable :: structure_elements
integer, dimension(:), allocatable :: structure_stochiometry
integer :: loop, status, size, stat
integer, parameter :: line_buf_len= 1024*4
character(LEN=line_buf_len) :: InS
type(unitcell), dimension(:), allocatable :: unit_cell
real(real12), dimension(:,:), allocatable :: atomic_positions
real(real12) :: structure_factor
real(real12), dimension(:,:), allocatable :: bondcutoff
logical ::  OK, set

allocate(unit_cell(1))
write(pathway,*) "rm -r bon bad 4body bon2 bad2 4body2 bon3 bad3 4body3"
call execute_command_line(pathway)
strucloop:do structure_loop=1, prev_structures
   !DEBUG
   !strucloop:do structure_loop=303, 304 
 write(pathway,'(A,I0.3,A)') "pos/POSCAR_",structure_loop,"/KPOINTS" 
 !write(*,*) pathway
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
 write(*,*) command
 call execute_command_line(command)
 write(path_test,'(A,I0.3,A)') "pos/POSCAR_",structure_loop,"/full.txt"

 inquire(file=trim(adjustl(path_test)),exist=dir_e)

 if(dir_e.eqv..FALSE.) then
    write(pathway,'(A,I0.3,A)') "pos/POSCAR_",structure_loop,"/POSCAR"
    inquire(file=trim(adjustl(pathway)),exist=dir_e)
    if(dir_e.eqv..FALSE.) cycle strucloop
 end if


 write(*,*) pathway, structure_loop 


 call poscar_read(pathway,unit_cell(1),structure_elements,structure_stochiometry, structure_factor,coordinate_type&
      &,atomic_positions, header)
 !write(*,*) pathway
 allocate(array(1,sum(structure_stochiometry)))
 allocate(repeatedarray(1,27*sum(structure_stochiometry)))
 k=0
 !write(*,*) size(array)
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




 write(*,*) "MANUALLY CHANGE"
 do i=1, sum(structure_stochiometry)
    call atomrepeater(1,array(1,i)%position,array(1,i)%name,&
         &repeatedarray,unit_cell,i,sum(structure_stochiometry)) 
 end do





 write(*,*) "YOU SHOULD ADD IN HERE A FUNCTION WHICH CERATES DYNAMIC BONDCUTOFF VALUES FOR DIFFERENT ELEMTNS."
 write(*,*) structure_stochiometry, sum(structure_stochiometry)
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
 write(*,*) pathway
 allocate(array(1,sum(structure_stochiometry)))
 allocate(repeatedarray(1,27*sum(structure_stochiometry)))
 k=0
 write(*,*) size(array)
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

 write(*,*) "YOU SHOULD ADD IN HERE A FUNCTION WHICH CERATES DYNAMIC BONDCUTOFF VALUES FOR DIFFERENT ELEMTNS."
 do i=1, sum(structure_stochiometry) 
    write(*,*) prev_structures,structure_loop,prev_structures+structure_loop
    call generatebondfiles(2,1,structure_loop-1,array,repeatedarray,size(structure_stochiometry),&
         &structure_stochiometry,i,structure_elements)
    write(*,*) "11"
    call generate4files(2,1,structure_loop-1,array,repeatedarray,size(structure_stochiometry),&
         &structure_stochiometry,i,bondcutoff,structure_elements)
    write(*,*) "22"
    call generateanglefiles(2,1,structure_loop-1,array,repeatedarray,size(structure_stochiometry),&
         &structure_stochiometry,i,bondcutoff,structure_elements)
    write(*,*) "33"
 end do
 deallocate(array)
 deallocate(repeatedarray)
 deallocate(bondcutoff)
end do strucloop2
strucloop3:do structure_loop=1, prev_structures
   ! strucloop3:do structure_loop=303, 304

 write(pathway,'(A,I0.3,A)') "pos/POSCAR_",structure_loop,"/RELAX/RELAX2/CONTCAR" 
 write(*,*) pathway

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
 write(*,*) pathway
 allocate(array(1,sum(structure_stochiometry)))
 allocate(repeatedarray(1,27*sum(structure_stochiometry)))
 k=0
 write(*,*) size(array)
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
 write(*,*) bondcutoff
 stop

 do i=1, sum(structure_stochiometry)
    call atomrepeater(1,array(1,i)%position,array(1,i)%name,&
         &repeatedarray,unit_cell,i,sum(structure_stochiometry)) 
 end do

 write(*,*) "YOU SHOULD ADD IN HERE A FUNCTION WHICH CERATES DYNAMIC BONDCUTOFF VALUES FOR DIFFERENT ELEMTNS."
 do i=1, sum(structure_stochiometry) 
    write(*,*) prev_structures,structure_loop,prev_structures+structure_loop
    call generatebondfiles(3,1,structure_loop-1,array,repeatedarray,size(structure_stochiometry),&
         &structure_stochiometry,i,structure_elements)
    write(*,*) "11"
    call generate4files(3,1,structure_loop-1,array,repeatedarray,size(structure_stochiometry),&
         &structure_stochiometry,i,bondcutoff,structure_elements)
    write(*,*) "22"
    call generateanglefiles(3,1,structure_loop-1,array,repeatedarray,size(structure_stochiometry),&
         &structure_stochiometry,i,bondcutoff,structure_elements)
    write(*,*) "33"
 end do
 deallocate(array)
 deallocate(repeatedarray)
 deallocate(bondcutoff)
end do strucloop3

end subroutine regenerate_distribution_files





subroutine bond_evolution(mode)
character(1024) :: name, read_element_pairing, read_element
integer :: prev_structures, mode,i, nbin, stat, exitst, exitst2, exitst3
logical :: dir_e
real(real12), dimension(:,:), allocatable :: gaussian
real(real12), dimension(2) :: read_in, norma_vector
real(real12) :: sigma, bondcut, dist_height

if(mode.eq.1) then 
 write(name,'(A)') "rm tmp_energies.txt"
 call execute_command_line(name,wait=.TRUE.) 
 write(name,'(A)') "prevstructures.txt" 
 open(13,file=name,status="old") 
 read(13, *) prev_structures
 write(*,*) prev_structures, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
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
 write(*,*) "#################################################"
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
 write(*,*) trim(adjustl(name)), "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
 write(name,'(A,A,A)') "Devolved/", trim(adjustl(read_element_pairing)), "_evolved_bondlength_gauss"
 open(102,file=trim(adjustl(name)))

 do while(1.eq.1)
    read(101,*,IOSTAT=stat) read_in
    !write(*,*) read_in
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
 write(*,*) trim(adjustl(name)), "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
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
 write(*,*) trim(adjustl(name)), "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
 write(name,'(A,A,A)') "Devolved/", trim(adjustl(read_element)), "_evolved_4body_gauss"
 open(102,file=trim(adjustl(name)))





 do while(1.eq.1)
    read(101,*,IOSTAT=stat) read_in
    !write(*,*) read_in
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



subroutine angledistribution(atoma,resul)
integer :: i,j,k,length
character(3) :: atoma 
character(1024) :: filename, buffer 
logical :: dir_e 
real(real12), dimension(:), allocatable :: resul

if(allocated(resul)) deallocate(resul)
write(filename,*) "angledist", trim(adjustl(atoma)),".in"
filename=trim(adjustL(filename))
inquire(file=trim(adjustl(filename)), exist=dir_e)
if(dir_e) then 
   !write(*,*) "Found distribution for ", atoma
else 
 write(*,*) "STOP. YOU VIOLATED THE LAW. YOU MUST HAVE DISTRIBUTIONS FOR", atoma
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
! write(*,*) buffer
write(command,'(4A,X,A,X,A,X,A)') "(echo",' "',trim(adjustl(buffer)),'" | ','sed -e "s/^.*=//g")>tmp.txt'
call execute_command_line(command)
!write(*,*) command
open(72,file="tmp.txt")

if(a.le.3) then
 allocate(b(1))
 read(72,*) b
 close(72)
else if((a.eq.4)) then 
 close(72)
 call invar(2,b,c)

 open(71, file="Infile.txt")
 rewind(71)
 do i=1, a
    read(71,'(A,X,A,X,A,X,A,X,A,X,A,X,A,X,A)') buffer
 end do
 close(71)
 ! write(*,*) buffer
 write(command,'(4A,X,A,X,A,X,A)') "(echo",' "',trim(adjustl(buffer)),'" | ','sed -e "s/^.*=//g")>tmp.txt'
 call execute_command_line(command)
 ! write(*,*) command
 open(72,file="tmp.txt")

 tmpvar=b(1)
 deallocate(b)
 allocate(c(tmpvar))
 read(72,*) c
 close(72)

else if((a.eq.5)) then
 close(72)
 call invar(2,b,c)

 open(71, file="Infile.txt")
 rewind(71)
 do i=1, a
    read(71,'(A,X,A,X,A,X,A,X,A,X,A,X,A,X,A)') buffer
 end do
 close(71)
 ! write(*,*) buffer
 write(command,'(4A,X,A,X,A,X,A)') "(echo",' "',trim(adjustl(buffer)),'" | ','sed -e "s/^.*=//g")>tmp.txt'
 call execute_command_line(command)
 ! write(*,*) command
 open(72,file="tmp.txt")

 tmpvar=b(1)
 deallocate(b)
 allocate(b(tmpvar))
 read(72,*) b
 !write(*,*) b
 close(72)
else if (a.eq.13) then 
 open(71, file="Infile.txt")
 rewind(71)
 do i=1, a
    read(71,'(A,X,A,X,A,X,A,X,A,X,A,X,A,X,A)') buffer
 end do
 close(71)
 !write(*,*) buffer
 write(command,'(4A,X,A,X,A,X,A)') "(echo",' "',trim(adjustl(buffer)),'" | ','sed -e "s/^.*=//g")>tmp.txt'
 call execute_command_line(command)
 !write(*,*) command
 open(72,file="tmp.txt")


 allocate(c(1))
 read(72,*) c
 write(*,*) c
 close(72)
else if(a.eq.14) then    
 close(72) 
 write(command,*) "sed -n 's|buildmap_resolution=||p'<Infile.txt>tmp.txt"
 call execute_command_line(command)
 write(*,*) command
 open(72,file="tmp.txt")
 allocate(b(3))
 read(72,*) b
 close(72)
else if(a.eq.16) then 
 close(72)
 write(command,*) "sed -n 's|VoidPseudoScanratio=||p'<Infile.txt>tmp.txt"
 call execute_command_line(command)
 write(*,*) command
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
call execute_command_line(command)
end subroutine invar





subroutine initialisehost(leng)
type(unitcell), dimension(:), allocatable :: formula
character(1024) :: tmp, name, location
character(1024), dimension(:), allocatable :: elnames,tmpelnames
real(real12) :: cellmultiplier,meanvol
integer :: l,k,j,i,structno, ecount,addedelements, eltot, leng,structures, prev_structures
type (atom), dimension(:,:), allocatable :: tmplist,atomlist,tmplist2
integer, dimension(:), allocatable :: stochio, tmpstochio, ts2, tmpstochiotot


open(50, file=trim(adjustl(filename_host)))
do i=1, 5
 read(50, * )
end do
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
       !write(*,*) tmp(i:i+1)
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
end subroutine initialisehost







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



end module help