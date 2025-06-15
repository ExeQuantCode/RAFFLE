program test_cache
  use raffle__io_utils
  use raffle__cache
  use raffle__constants, only: real32
  implicit none

  logical :: success = .true.

  call test_retrieve_unallocated(success)
  call test_store_retrieve_probability_density(success)

  !-----------------------------------------------------------------------------
  ! check for any failed tests
  !-----------------------------------------------------------------------------
  write(*,*) "----------------------------------------"
  if(success)then
     write(*,*) 'test_cache passed all tests'
  else
     write(0,*) 'test_cache failed one or more tests'
     stop 1
  end if

contains

  subroutine test_store_retrieve_probability_density(success)
    implicit none
    logical, intent(inout) :: success
    real(real32), allocatable :: test_data(:,:), retrieved_data(:,:)
    integer :: i, j

    ! Create test data - a 3x3 matrix with values
    allocate(test_data(3,3))
    do i = 1, 3
      do j = 1, 3
        test_data(i,j) = real(i*10 + j, real32)
      end do
    end do

    ! Store the test data in the cache
    call store_probability_density(test_data)

    ! Retrieve the data from the cache
    retrieved_data = retrieve_probability_density()

    ! Check dimensions match
    if (.not. all(shape(retrieved_data) == shape(test_data))) then
      write(0,*) 'test_store_retrieve_probability_density: array shapes do not match'
      success = .false.
      return
    end if

    ! Check all values match
    do i = 1, 3
      do j = 1, 3
        if (abs(retrieved_data(i,j) - test_data(i,j)) .gt. 1.e-6_real32) then
          write(0,*) &
               'test_store_retrieve_probability_density: values do not match at', &
               i, j, retrieved_data(i,j), test_data(i,j)
          success = .false.
          return
        end if
      end do
    end do

    ! Clean up
    if (allocated(test_data)) deallocate(test_data)
    if (allocated(retrieved_data)) deallocate(retrieved_data)

    write(*,*) 'test_store_retrieve_probability_density: PASSED'
  end subroutine test_store_retrieve_probability_density

  subroutine test_retrieve_unallocated(success)
    implicit none
    logical, intent(inout) :: success
    real(real32), allocatable :: retrieved_data(:,:)

    ! Try to retrieve data when cache is not allocated
    retrieved_data = retrieve_probability_density()

    ! Verify returned array is a scalar with value 0
    if (size(retrieved_data) .ne. 1) then
      write(0,*) &
           'test_retrieve_unallocated: returned array should be scalar but has size', &
           size(retrieved_data)
      success = .false.
      return
    end if

    if (abs(retrieved_data(1,1)) .gt. 1.e-6_real32) then
      write(0,*) &
           'test_retrieve_unallocated: returned value should be 0 but is', &
           retrieved_data(1,1)
      success = .false.
      return
    end if

    ! Clean up
    if (allocated(retrieved_data)) deallocate(retrieved_data)

    write(*,*) 'test_retrieve_unallocated: PASSED'
  end subroutine test_retrieve_unallocated

end program test_cache
