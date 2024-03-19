module contributions
  use constants, only: real12
  use misc, only: jump


contains

!!!#############################################################################
!!! evalulate the contribution of bondlength to the energy
!!!#############################################################################
  subroutine evaluate_contribution (element_1, element_2, &
       bondlength, amplitude)
    implicit none
    character(*), intent(in) :: element_1, element_2
    real(real12), intent(in) :: bondlength
    real(real12), intent(out) :: amplitude
    
    integer :: i, ierror, unit
    real(real12) :: step_size, fraction
    character(1024) :: filename
    real(real12), dimension(2) :: amplitude_list = 0._real12
    real(real12), dimension(2) :: bondlength_list = 0._real12

    
    !! Open the file containing evolved bondlengths
    write(filename,'("Devolved/",A,"_",A,"_evolved_bondlength_gauss")') &
         trim(adjustl(element_1)), trim(adjustl(element_2))
    open(newunit=unit, file=trim(adjustl(filename)))
    amplitude = 0._real12
    
    !! read the step size
    read(unit,*,iostat=ierror) step_size
    if(is_iostat_end(ierror)) then
       close(unit); return
    end if
    call jump(unit, floor(bondlength/step_size) - 1)
    do i = 1, 2
       read(unit,*,iostat=ierror) bondlength_list(i), amplitude_list(i)
       write(*,*) i, bondlength_list(i), amplitude_list(i)
       if(is_iostat_end(ierror)) then
          close(unit)
          return
       end if
    end do

    !!Assuming the line between points can be modelled as linear
    fraction = (bondlength - bondlength_list(1)) / step_size
 
    amplitude = amplitude_list(1) + &
       (amplitude_list(2) - amplitude_list(1)) * fraction

    !! THIS IS A HACK TO MAKE THE BOND LENGTH CONTRIBUTION DECAY
    !!Arbitrary could be added for when bondlength contribution starts to decay,
    !! could be changed to be bond specific (should be)!
    if (bondlength.lt.1.6_real12) then 
       amplitude = amplitude 
    else if (bondlength.lt.3._real12) then  
       amplitude = amplitude
    else 
       amplitude = amplitude / 26._real12
    end if
 
    close(unit)

  end subroutine evaluate_contribution
!!!#############################################################################
 
 
!!!#############################################################################
!!! evalulate the contribution of bondangle to the energy
!!!#############################################################################
 subroutine evaluate_angle_contribution (element_1,&
  &bondangle_target,return_slot)
 character(*) :: element_1
 character(1024) :: filename
 integer :: stat, i
 real(real12), dimension(3) :: value_read, angles_read
 real(real12) :: a,b,c,d, bondangle_target, return_slot
 
 integer :: unit
 
 write(filename,'(A,A,A,A,A)') "Devolved/",trim(adjustl(element_1))&
    &,"_evolved_angles_gauss"
 open(newunit=unit, file=trim(adjustl(filename)))
 angles_read=0
 value_read=0
 angles_read(2)=bondangle_target
 
 
 master :do while(1.eq.1)
 
  angles_read(1)=angles_read(3)
  value_read(1)=value_read(3)
 
  read(unit,*,IOSTAT=stat) angles_read(3), value_read(3)
 
  IF(IS_IOSTAT_END(stat)) then
     return_slot=0
     close(unit)
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
     close(unit)
     if(return_slot.gt.0) then
        !write(*,*) value_read(1), value_read(3), d,b, return_slot
     end if
 
     exit master
  end IF
 
 
  IF(IS_IOSTAT_END(stat)) exit
 end do master
 end subroutine evaluate_angle_contribution
!!!#############################################################################
 

!!!#############################################################################
!!! evalulate the contribution of 4body to the energy
!!!#############################################################################
 subroutine evaluate_4body_contribution (element_1,&
  &bondangle_target,return_slot)
 character(*) :: element_1
 character(1024) :: filename
 integer :: stat, i
 real(real12), dimension(3) :: value_read
 real(real12), dimension(3) :: angles_read
 real(real12) :: a,b,c,d, bondangle_target, return_slot
 
 integer :: unit

 write(filename,'(A,A,A,A,A)') "Devolved/",trim(adjustl(element_1))&
    &,"_evolved_4body_gauss"
 open(newunit=unit, file=trim(adjustl(filename)))
 angles_read=0
 value_read=0
 angles_read(2)=bondangle_target
 if(bondangle_target.eq.1000) then 
  return_slot=1
  close(unit)
 else
  master :do while(1.eq.1)
 
     angles_read(1)=angles_read(3)
     value_read(1)=value_read(3)
 
     read(unit,*,IOSTAT=stat) angles_read(3), value_read(3)
     IF(IS_IOSTAT_END(stat)) then
        return_slot=0
        close(unit)
        exit
     end IF
     IF((angles_read(2).lt.angles_read(3)).and.&
          &(angles_read(2).gt.angles_read(1))) then
        !!Assuming the line between points can be modelled as linear
        a=angles_read(1)
        b=angles_read(2)
        c=angles_read(3)
 
        d=(b-a)/(c-a)
        !write(*,*) a, b, c
 
        if(value_read(3).lt.value_read(1)) d=1-d
 
        !value_read(2)=value_read(1)+dble((value_read(3)-value_read(1))*d)
        value_read(2)=(value_read(3)-value_read(1))/(angles_read(3)-angles_read(1))*(angles_read(2)-angles_read(1))+value_read(1)
 
 
        !if(value_read(2).lt.0.01) value_read(2)=0 
 
        !!Arbitrary could be added for when bondlength contribution starts to decay,
        !! could be changed to be bond specific (should be)!
 
 
        return_slot=value_read(2)
        !write(*,*) "Angle read in, between plane and AD vector: value determined for that angle"
        if(return_slot.gt.0) then 
           !            write(*,*) value_read(1), value_read(3), d,b, return_slot
        end if
        close(unit)
        exit master
     end IF
     IF(IS_IOSTAT_END(stat)) exit
  end do master
 end if
 end subroutine evaluate_4body_contribution
!!!#############################################################################

end module contributions