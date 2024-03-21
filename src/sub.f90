module help
  use constants, only: real12, pi
  use geom, only: &
       get_bondlength, get_bondangle, get_dihedral_angle
  use atomtype
  use vasp_file_handler, only: unitcell
  use file_generator
  use contributions
  implicit none




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

end module help