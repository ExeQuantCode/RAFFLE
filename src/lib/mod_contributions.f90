module contributions
  use constants, only: real12
  use misc, only: jump
  implicit none

  private

  public :: get_2body_contribution
  public :: get_3body_contribution
  public :: get_4body_contribution


contains

!!!#############################################################################
!!! evalulate the contribution of bondlength to the energy
!!!#############################################################################
  function get_2body_contribution (element_1, element_2, bondlength) &
       result(amplitude)
    implicit none
    character(*), intent(in) :: element_1, element_2
    real(real12), intent(in) :: bondlength
    real(real12) :: amplitude
    
    integer :: i, ierror, unit
    real(real12) :: step_size, fraction
    character(1024) :: filename
    real(real12), dimension(2) :: amplitude_list = 0._real12
    real(real12), dimension(2) :: bondlength_list = 0._real12


    amplitude = 0._real12
    !! Open the file containing evolved bondlengths
    write(filename,'("Devolved/",A,"_",A,"_evolved_bondlength_gauss")') &
         trim(adjustl(element_1)), trim(adjustl(element_2))
    open(newunit=unit, file=trim(adjustl(filename)))
    
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
    close(unit)

    !!Assuming the line between points can be modelled as linear
    fraction = (bondlength - bondlength_list(1)) / step_size
    amplitude = amplitude_list(1) + &
         (amplitude_list(2) - amplitude_list(1)) * fraction

    !! THIS IS A HACK TO MAKE THE BOND LENGTH CONTRIBUTION DECAY
    !! Could be added for when bondlength contribution starts to decay,
    !! Could be changed to be bond specific (should be)!
    if (bondlength.lt.1.6_real12) then 
       amplitude = amplitude 
    else if (bondlength.lt.3._real12) then  
       amplitude = amplitude
    else 
       amplitude = amplitude / 26._real12
    end if
 

  end function get_2body_contribution
!!!#############################################################################
 
 
!!!#############################################################################
!!! evalulate the contribution of bondangle to the energy
!!!#############################################################################
  function get_3body_contribution (element_1, bondangle) &
       result(amplitude)
    implicit none
    real(real12), intent(in) :: bondangle
    character(*), intent(in) :: element_1
    real(real12) :: amplitude

    integer :: i, ierror, unit
    real(real12) :: step_size, fraction
    character(1024) :: filename
    real(real12), dimension(2) :: amplitude_list = 0._real12
    real(real12), dimension(2) :: bondangle_list = 0._real12


    amplitude = 0._real12
    write(filename,'(A,A,A,A,A)') "Devolved/",trim(adjustl(element_1))&
       &,"_evolved_angles_gauss"
    open(newunit=unit, file=trim(adjustl(filename)))

    !! read the step size
    read(unit,*,iostat=ierror) step_size
    if(is_iostat_end(ierror)) then
       close(unit); return
    end if
    call jump(unit, floor(bondangle/step_size) - 1)
    do i = 1, 2
       read(unit,*,iostat=ierror) bondangle_list(i), amplitude_list(i)
       write(*,*) i, bondangle_list(i), amplitude_list(i)
       if(is_iostat_end(ierror)) then
          close(unit)
          return
       end if
    end do
    close(unit)

    !!Assuming the line between points can be modelled as linear
    fraction = (bondangle - bondangle_list(1)) / step_size
    amplitude = amplitude_list(1) + &
         (amplitude_list(2) - amplitude_list(1)) * fraction

  end function get_3body_contribution
!!!#############################################################################
 

!!!#############################################################################
!!! evalulate the contribution of 4body to the energy
!!!#############################################################################
  function get_4body_contribution(element_1, bondangle) &
       result(amplitude)
    implicit none
    real(real12), intent(in) :: bondangle
    character(*), intent(in) :: element_1
    real(real12) :: amplitude

    integer :: i, ierror, unit
    real(real12) :: step_size, fraction
    character(1024) :: filename
    real(real12), dimension(2) :: amplitude_list = 0._real12
    real(real12), dimension(2) :: bondangle_list = 0._real12
 

    amplitude = 0._real12
    write(filename,'(A,A,A,A,A)') "Devolved/",trim(adjustl(element_1))&
       &,"_evolved_4body_gauss"
    open(newunit=unit, file=trim(adjustl(filename)))
    !!!##########################
    !!! WHAT IS THIS HERE FOR???
    if(bondangle.eq.1000) then 
       amplitude = 1._real12
       close(unit)
       return
    end if
    !!!##########################
 
    !! read the step size
    read(unit,*,iostat=ierror) step_size
    if(is_iostat_end(ierror)) then
       close(unit); return
    end if
    call jump(unit, floor(bondangle/step_size) - 1)
    do i = 1, 2
       read(unit,*,iostat=ierror) bondangle_list(i), amplitude_list(i)
       write(*,*) i, bondangle_list(i), amplitude_list(i)
       if(is_iostat_end(ierror)) then
          close(unit)
          return
       end if
    end do
    close(unit)
  
    !!Assuming the line between points can be modelled as linear
    fraction = (bondangle - bondangle_list(1)) / step_size
    amplitude = amplitude_list(1) + &
         (amplitude_list(2) - amplitude_list(1)) * fraction

  end function get_4body_contribution
!!!#############################################################################

end module contributions