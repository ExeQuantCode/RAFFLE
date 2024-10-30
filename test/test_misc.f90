program test_misc
  use raffle__io_utils
  use raffle__misc
  use raffle__constants, only: real32
  implicit none

  logical :: success = .true.

  test_error_handling = .true.

  call test_sort_str(success)
  call test_sort_str_order(success)
  call test_isort1D(success)
  call test_rsort1D(success)
  call test_iset(success)
  call test_rset(success)
  call test_cset(success)
  call test_ishuffle(success)
  call test_rshuffle(success)
  call test_icount(success)
  call test_to_upper(success)
  call test_to_lower(success)
  call test_strip_null(success)
  call test_grep(success)
  call test_jump(success)
  call test_rswap(success)
  call test_rswap_vec(success)


  !-----------------------------------------------------------------------------
  ! check for any failed tests
  !-----------------------------------------------------------------------------
  write(*,*) "----------------------------------------"
  if(success)then
     write(*,*) 'test_misc passed all tests'
  else
     write(0,*) 'test_misc failed one or more tests'
     stop 1
  end if

contains

  subroutine test_sort_str(success)
    implicit none
    logical, intent(inout) :: success
    character(len=20), dimension(5) :: list
    character(len=20), dimension(5) :: expected_list
    list = [ &
         'banana    ', 'apple     ', 'cherry    ', 'date      ', 'elderberry' &
    ]
    expected_list = [ &
         'apple     ', 'banana    ', 'cherry    ', 'date      ', 'elderberry' &
    ]
    call sort_str(list)
    call assert( &
         all(list .eq. expected_list), &
         'test_sort_str failed', success &
    )

    list = [ &
         'Banana    ', 'cherry    ', 'banana    ', 'date      ', 'elderberry' &
    ]
    expected_list = [ &
         'Banana    ', 'banana    ', 'cherry    ', 'date      ', 'elderberry' &
    ]
    call sort_str(list, lcase=.true.)
    call assert( &
         all(list .eq. expected_list), &
         'test_sort_str failed with ignore case', success &
    )
  end subroutine test_sort_str

  subroutine test_sort_str_order(success)
    implicit none
    logical, intent(inout) :: success
    character(len=20), dimension(5) :: list
    integer, dimension(5) :: expected_order = (/ 2, 1, 3, 4, 5 /)
    integer, dimension(:), allocatable :: order
    list = [ &
         'banana    ', 'apple     ', 'cherry    ', 'date      ', 'elderberry' &
    ]
    allocate(order, source=sort_str_order(list))
    call assert( &
         all(order .eq. expected_order), &
         'test_sort_str_order failed', success &
    )

    list = [ &
         'banana    ', 'cherry    ', 'Banana    ', 'date      ', 'elderberry' &
    ]
    expected_order = [ 1, 3, 2, 4, 5 ]
    order = sort_str_order(list, lcase=.true.)
    call assert( &
         all(order .eq. expected_order), &
         'test_sort_str_order failed with ignore case', success &
    )
  end subroutine test_sort_str_order

  subroutine test_isort1D(success)
    implicit none
    logical, intent(inout) :: success
    integer, dimension(5) :: arr = [5, 3, 4, 1, 2]
    integer, dimension(5) :: expected_arr = [1, 2, 3, 4, 5]
    call sort1D(arr)
    call assert( &
         all(arr .eq. expected_arr), &
         'test_isort1D failed', success &
    )
    expected_arr = [5, 4, 3, 2, 1]
    call sort1D(arr, reverse=.true.)
    call assert( &
         all(arr .eq. expected_arr), &
         'test_isort1D failed with reverse', success &
    )
  end subroutine test_isort1D

  subroutine test_rsort1D(success)
    implicit none
    logical, intent(inout) :: success
    real(real32), dimension(5) :: arr = &
         [5._real32, 3._real32, 4._real32, 1._real32, 2._real32]
    real(real32), dimension(5) :: expected_arr = &
         [1._real32, 2._real32, 3._real32, 4._real32, 5._real32]
    call sort1D(arr)
    call assert( &
         all( abs(arr - expected_arr) .lt. 1.E-6), &
         'test_rsort1D failed', success &
    )
    expected_arr = [5._real32, 4._real32, 3._real32, 2._real32, 1._real32]
    call sort1D(arr, reverse=.true.)
    call assert( &
         all(arr .eq. expected_arr), &
         'test_rsort1D failed with reverse', success &
    )
  end subroutine test_rsort1D

  subroutine test_iset(success)
    implicit none
    logical, intent(inout) :: success
    integer, dimension(:), allocatable :: arr
    integer, dimension(:), allocatable :: expected_arr
    allocate(arr(6))
    arr = [1, 2, 2, 3, 3, 3]
    allocate(expected_arr(3))
    expected_arr = [1, 2, 3]
    call set(arr)
    call assert( &
         all(arr .eq. expected_arr), &
         'test_iset failed', success &
    )
  end subroutine test_iset

  subroutine test_rset(success)
    implicit none
    logical, intent(inout) :: success
    real(real32), dimension(:), allocatable :: arr
    real(real32), dimension(:), allocatable :: expected_arr
    allocate(arr(6))
    arr = [1._real32, 2._real32, 2._real32, 3._real32, 3._real32, 3._real32]
    allocate(expected_arr(3))
    expected_arr = [1._real32, 2._real32, 3._real32]
    call set(arr)
    call assert( &
         all( abs(arr - expected_arr) .lt. 1.E-6), &
         'test_rset failed', success &
    )
    arr = [1._real32, 2._real32, 2.00001_real32, 3._real32, 3._real32]
    expected_arr = [1._real32, 2._real32, 2.00001_real32, 3._real32]
    call set(arr, tol=1.E-6)
    call assert( &
         all( abs(arr - expected_arr) .lt. 1.E-6), &
         'test_rset failed with lower tolerance', success &
    )
  end subroutine test_rset

  subroutine test_cset(success)
    implicit none
    logical, intent(inout) :: success
    character(len=20), dimension(:), allocatable :: arr
    character(len=20), dimension(:), allocatable :: expected_arr
    allocate(arr(6))
    arr(:) = [ 'apple ', 'banana', 'banana', 'cherry', 'cherry', 'cherry' ]
    allocate(expected_arr(3))
    expected_arr(:) = [ 'apple ', 'banana', 'cherry' ]
    call set(arr)
    call assert( &
         all(arr .eq. expected_arr), &
         'test_cset failed', success &
    )
    deallocate(arr)
    allocate(arr(6))
    arr = [ 'apple ', 'Banana', 'banana', 'cherry', 'cherry', 'cherry' ]
    expected_arr(:) = [ 'apple ', 'banana', 'cherry' ]
    call set(arr, lcase=.true.)
    call assert( &
         all(arr .eq. expected_arr), &
         'test_cset failed with ignore case', success &
    )
    deallocate(arr)
    allocate(arr(6))
    arr = [ 'apple ', 'Banana', 'banana', 'cherry', 'cherry', 'cherry' ]
    deallocate(expected_arr)
    allocate(expected_arr(6))
    expected_arr(:4) = &
         [ 'Banana', 'apple ', 'banana', 'cherry' ]
    expected_arr(5:) = ''
    call set(arr, lkeep_size=.true.)
    call assert( &
         size(arr) .eq. 6, &
         'test_cset failed to keep size', success &
    )
    call assert( &
         all(arr(:4) .eq. expected_arr(:4)), &
         'test_cset failed with keep_size', success &
    )
  end subroutine test_cset

  subroutine test_ishuffle(success)
    implicit none
    logical, intent(inout) :: success
    integer  :: i
    logical :: ltmp1
    integer :: arr(1,5)
    integer :: original_arr(1,5)

    arr(1,:) = [1, 2, 3, 4, 5]
    original_arr(1,:) = arr(1,:)
    call shuffle(arr, dim=2, seed=0)
    ltmp1 = .true.
    do i = 1, size(arr,dim=2)
       if(all(abs(arr(1,i) - original_arr(1,:)).gt.0)) then
          ltmp1 = .false.
          exit
       end if
    end do
    call assert(ltmp1, "ishuffle failed", success)
  end subroutine test_ishuffle

  subroutine test_rshuffle(success)
    implicit none
    logical, intent(inout) :: success
    integer  :: i
    logical :: ltmp1
    real(real32) :: arr(1,5)
    real(real32) :: original_arr(1,5)

    arr(1,:) = [1._real32, 2._real32, 3._real32, 4._real32, 5._real32]
    original_arr(1,:) = arr(1,:)
    call shuffle(arr, dim=2, seed=0)
    ltmp1 = .true.
    do i = 1, size(arr,dim=2)
       if(all(abs(arr(1,i) - original_arr(1,:)).gt.1.E-6)) then
          ltmp1 = .false.
          exit
       end if
    end do
    call assert(ltmp1, "rshuffle failed", success)
  end subroutine test_rshuffle

  subroutine test_icount(success)
    implicit none
    logical, intent(inout) :: success
    character(len=20) :: line = "apple,banana,cherry"
    integer :: count

    count = icount(line, ",")
    call assert(count .eq. 3, "icount failed", success)
  end subroutine test_icount

  subroutine test_to_upper(success)
    implicit none
    logical, intent(inout) :: success
    character(len=10) :: str
    character(len=10) :: expected_str

    str = "hello"
    expected_str = "HELLO"
    str = to_upper(str)
    call assert(trim(str) .eq. trim(expected_str), "to_upper failed", success)
  end subroutine test_to_upper

  subroutine test_to_lower(success)
    implicit none
    logical, intent(inout) :: success
    character(len=10) :: str
    character(len=10) :: expected_str

    str = "HELLO"
    expected_str = "hello"
    str = to_lower(str)
    call assert(trim(str) .eq. trim(expected_str), "to_lower failed", success)
  end subroutine test_to_lower

  subroutine test_strip_null(success)
    implicit none
    logical, intent(inout) :: success
    character(len=16) :: str
    character(len=16) :: expected_str

    str = "hello"//char(0)//"world"
    expected_str = "hello"
    str = strip_null(str)
    call assert(trim(str) .eq. trim(expected_str), "strip_null failed", success)
  end subroutine test_strip_null

  subroutine test_grep(success)
    implicit none
    logical, intent(inout) :: success
    integer :: unit
    logical :: success_tmp

    ! Create a temporary scratch file for testing
    open(newunit=unit, status='scratch', action='readwrite')
    write(unit, '(A)') 'This is a test line.'
    write(unit, '(A)') 'Another line with test pattern.'
    write(unit, '(A)') 'Yet another line.'
    rewind(unit)

    ! Test case 1: Pattern found in the middle of the line
    call grep(unit, 'test pattern', success=success_tmp)
    call assert(success_tmp, 'Pattern not found', success)
    rewind(unit)

    ! Test case 2: Pattern not found
    call grep(unit, 'nonexistent pattern', lstart=.true., success=success_tmp)
    call assert(.not. success_tmp, 'Nonexistent pattern found', success)
    rewind(unit)

    ! Test case 3: Pattern at the start of the line
    call grep(unit, 'This is', lline=.true., success=success_tmp)
    call assert(success_tmp, 'Pattern at start of line not found', success)
    close(unit, status='delete')

  end subroutine test_grep

  subroutine test_jump(success)
    implicit none
    logical, intent(inout) :: success
    integer :: i, j, unit

    open(newunit=unit, status='scratch', action='readwrite')
    do i = 1, 10
       write(unit, '(I0)') i
    end do

    ! Test case 1: Jump to the end of the file
    do i = 1, 10
       rewind(unit)
       call jump(unit, i)
       backspace(unit)
       read(unit, *) j
       call assert(j .eq. i, 'Jump failed', success)
    end do

  end subroutine test_jump

  subroutine test_rswap(success)
    implicit none
    logical, intent(inout) :: success
    real(real32) :: a = 1._real32
    real(real32) :: b = 2._real32
    real(real32) :: expected_a = 2._real32
    real(real32) :: expected_b = 1._real32

    call swap(a, b)
    call assert( &
         abs(a - expected_a).lt. 1.E-6 .and. &
         abs(b - expected_b).lt. 1.E-6, &
         "rswap failed", success &
    )
    
  end subroutine test_rswap

  subroutine test_rswap_vec(success)
    implicit none
    logical, intent(inout) :: success
    real(real32), dimension(2) :: a = [1._real32, 2._real32]
    real(real32), dimension(2) :: b = [3._real32, 4._real32]

    call swap(a, b)
    call assert( &
         all( abs(a - [3._real32, 4._real32]) .lt. 1.E-6_real32 ) .and. &
         all( abs(b - [1._real32, 2._real32]) .lt. 1.E-6_real32 ), &
         "rswap_vec failed", success &
    )

  end subroutine test_rswap_vec

!###############################################################################

  subroutine assert(condition, message, success)
    implicit none
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    logical, intent(inout) :: success
    if (.not. condition) then
       write(0,*) "Test failed: ", message
       success = .false.
    end if
  end subroutine assert

end program test_misc
