module contributions
  use constants, only: real12


contains

!!!#############################################################################
!!! evalulate the contribution of bondlength to the energy
!!!#############################################################################
subroutine evaluate_contribution (element_a,element_b,&
  &bondlength_target,return_slot)
 character(*), intent(in) :: element_a, element_b
 real(real12), intent(in) :: bondlength_target
 real(real12), intent(out) :: return_slot
 character(1024) :: name 
 integer :: atom_number_a, atom_number_b, stat, i
 real(real12), dimension(3) :: value_read
 real(real12), dimension(3) :: bondlength_read  
 real(real12) :: a,b,c,d
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
!!!#############################################################################
 
 
!!!#############################################################################
!!! evalulate the contribution of bondangle to the energy
!!!#############################################################################
 subroutine evaluate_angle_contribution (element_a,&
  &bondangle_target,return_slot)
 character(*) :: element_a
 character(1024) :: name
 integer :: atom_number_a, atom_number_b, stat, i
 real(real12), dimension(3) :: value_read
 real(real12), dimension(3) :: angles_read
 real(real12) :: a,b,c,d, bondangle_target, return_slot
 
 
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
 subroutine evaluate_4body_contribution (element_a,&
  &bondangle_target,return_slot)
 character(*) :: element_a
 character(1024) :: name
 integer :: atom_number_a, atom_number_b, stat, i
 real(real12), dimension(3) :: value_read
 real(real12), dimension(3) :: angles_read
 real(real12) :: a,b,c,d, bondangle_target, return_slot
 
 
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
        close(199)
        exit master
     end IF
     IF(IS_IOSTAT_END(stat)) exit
  end do master
 end if
 end subroutine evaluate_4body_contribution
!!!#############################################################################

end module contributions