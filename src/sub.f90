module evolve
  use constants, only: real12
  use misc, only: grep
  use file_generator
  use contributions
  implicit none


  private

  public :: bond_evolution


contains

!!!#############################################################################
!!! 
!!!#############################################################################
  subroutine bond_evolution()
    implicit none
    character(1024) :: name, read_element_pairing, read_element
    integer :: prev_structures, i, nbin, stat, exitst, exitst2, exitst3
    logical :: dir_e
    real(real12), dimension(:,:), allocatable :: gaussian
    real(real12), dimension(2) :: read_in, norma_vector
    real(real12) :: sigma, bondcut, dist_height, energy
    character(50) :: buffer1, buffer2
    logical :: success

    integer :: previous_structures_unit, xml_unit


    
    !! read in the number of previous structures
    write(name,'(A)') "prevstructures.txt" 
    open(previous_structures_unit,file=name,status="old") 
    read(previous_structures_unit, *) prev_structures
    close(previous_structures_unit)
    write(*,*) prev_structures, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"

    !! read in the energies of the previous structures
    call execute_command_line("rm tmp_energies.txt",wait=.TRUE.)
    do i=1, prev_structures
       write(name,'(A,I0.3,A)') "pos/POSCAR_",i,"/vasprun.xml"
       open(newunit=xml_unit, file=trim(adjustl(name)), status="old")
       call grep(xml_unit,'   <i name="e_fr_energy">',lline=.true., success=success)
       if(.not.success) cycle
       backspace(xml_unit)
       read(xml_unit,*) buffer1, buffer2, energy
       close(xml_unit)
      
       !!! STORE THE ENERGY IN AN ARRAY
       ! either energy[], or store it alongside the basis, probably the latter
       ! probably new structure format of crystal
       ! where crystal contains lattice, basis, and energy
       !!! DO SOMETHING ABOUT NESTED RELAXATIONS 
 
    end do

    !!!! JUST WORK OUT THE GVECTORS HERE !!!!
    !! we have other codes for this
    !! read in the POSCAR (eventually from the xml file)
    !! calculate the neighbour tables





    write(*,*) "#################################################"
    call execute_command_line("./bond_evolution.sh tmp_energies.txt prevstructures.txt",wait=.TRUE.)
    !  The arguements in the next function are bond cap, bins and element 1 and 2. These should be changed from carbon et in this file when finished 
    !  The sigma here should be related to the sampling of the placement map; too small compared to map, and expanding regions to study in more 
    !  detail will be difficult
    bondcut=5
    nbin=1000
    sigma=0.1

    !dist_height=1/(sigma*(sqrt(2*3.141592654)))
    allocate(gaussian(2,nbin))
    call execute_command_line("mkdir -p Devolved; rm Devolved/*evolved_*_gauss; ./bond_database_builder.sh")
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
!!!#############################################################################

end module evolve