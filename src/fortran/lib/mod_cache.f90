module raffle__cache
  use raffle__constants, only: real32
  implicit none

  private
  public :: store_probability_density, retrieve_probability_density

  real(real32), allocatable, dimension(:,:), save :: cached_probability_density

contains

  subroutine store_probability_density(probability_density)
    implicit none
    real(real32), intent(in) :: probability_density(:,:)
    if (allocated(cached_probability_density)) &
         deallocate(cached_probability_density)
    allocate(cached_probability_density, source = probability_density)

  end subroutine store_probability_density

  function retrieve_probability_density() result(probability_density)
    implicit none
    real(real32), allocatable :: probability_density(:,:)
    if(.not.allocated(cached_probability_density)) then
       write(0,*) "Probability density not allocated. Returning zero array."
       probability_density = 0._real32
    else
       allocate(probability_density, source = cached_probability_density)
    end if
  end function retrieve_probability_density

end module raffle__cache
