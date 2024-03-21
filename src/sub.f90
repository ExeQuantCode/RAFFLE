module evolve
  use constants, only: real12
  use file_generator
  use contributions
  implicit none


  private

  public :: bond_evolution


contains

  subroutine bond_evolution(mode)
    implicit none
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

end module evolve