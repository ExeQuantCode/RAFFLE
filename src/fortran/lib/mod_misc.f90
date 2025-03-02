module raffle__misc
  !! Module contains various miscellaneous functions and subroutines.
  use raffle__constants, only: real32
  use raffle__io_utils, only: stop_program
  implicit none


  private

  public :: sort1D, sort2D, sort_str, sort_str_order
  public :: set, swap
  public :: shuffle
  public :: icount, grep, flagmaker
  public :: jump, file_check, touch, to_upper, to_lower
  public :: strip_null


  interface sort1D
     !! Sort a 1D array from min to max.
     procedure isort1D,rsort1D
  end interface sort1D

  interface set
     !! Reduce an array to its unique elements.
     procedure iset,rset, cset
  end interface set

  interface swap
     !! Swap two elements.
     procedure iswap, rswap, rswap_vec, cswap
  end interface swap

  interface shuffle
     !! Shuffle an array.
     procedure ishuffle, rshuffle
  end interface shuffle



contains

!###############################################################################
  subroutine sort_str(list, lcase)
    !! Sort a list of strings.
    implicit none

    ! Arguments
    character(*), dimension(:), intent(inout) :: list
    !! List of strings to be sorted.
    logical, optional, intent(in) :: lcase
    !! Optional. Boolean whether case insensitive sorting is required.
    !! Default is .false.

    ! Local variables
    integer :: i,loc
    !! Loop index.
    integer :: charlen
    !! Length of the strings.
    logical :: lcase_
    !! Boolean whether case insensitive sorting is required.
    character(:), allocatable, dimension(:) :: tlist
    !! Temporary list for case insensitive sorting.

    charlen = len(list(1))
    if(present(lcase))then
       lcase_ = lcase
    else
       lcase_ = .false.
    end if
    if(lcase_)then
       allocate(character(len=charlen) :: tlist(size(list)))
       tlist = list
       do i = 1, size(tlist)
          list(i) = to_upper(list(i))
       end do
    end if
    do i = 1, size(list)
       loc = minloc(list(i:),dim=1)
       if(loc.eq.1) cycle
       if(lcase_) call cswap(tlist(i),tlist(loc+i-1))
       call cswap(list(i),list(loc+i-1))
    end do
    if(lcase_) list=tlist
    
  end subroutine sort_str
!###############################################################################


!###############################################################################
  function sort_str_order(list,lcase) result(order)
    !! Sort a list of strings and return the order.
    implicit none

    ! Arguments
    character(*), dimension(:), intent(inout) :: list
    !! List of strings to be sorted.
    logical, optional, intent(in) :: lcase
    !! Optional. Boolean whether case insensitive sorting is required.
    !! Default is .false.

    ! Local variables
    integer :: i,loc
    !! Loop index.
    integer :: charlen
    !! Length of the strings.
    logical :: lcase_
    !! Boolean whether case insensitive sorting is required.
    character(:), allocatable, dimension(:) :: tlist
    !! Temporary list for case insensitive sorting.

    integer, allocatable, dimension(:) :: torder,order
    !! Order of the sorted list.

    charlen = len(list(1))
    lcase_ = .false.
    if(present(lcase))then
       if(lcase)then
          lcase_ = lcase
          allocate(character(len=charlen) :: tlist(size(list)))
          tlist = list
          do i = 1, size(tlist)
             list(i) = to_upper(list(i))
          end do
       end if
    end if

    allocate(torder(size(list)))
    do i = 1, size(list)
       torder(i) = i
    end do
    
    do i = 1, size(list)
       loc = minloc(list(i:),dim=1)
       if(loc.eq.1) cycle
       if(lcase_) call cswap(tlist(i),tlist(loc+i-1))
       call cswap(list(i),list(loc+i-1))
       call iswap(torder(i),torder(loc+i-1))
    end do
    
    allocate(order(size(list)))
    do i = 1, size(list)
       order(i) = findloc(torder,i,dim=1)
    end do
    
    if(lcase_) list=tlist
    
    return
  end function sort_str_order
!###############################################################################


!###############################################################################
  subroutine isort1D(arr1,arr2,reverse)
    !! Sort a 1D integer array from min to max.
    implicit none

    ! Arguments
    integer, dimension(:), intent(inout) :: arr1
    !! Array to be sorted.
    integer, dimension(:),intent(inout),optional :: arr2
    !! Optional. Second array to be sorted.
    logical, optional, intent(in) :: reverse
    !! Optional. Boolean whether to sort in reverse order.

    ! Local variables
    integer :: i,dim,loc
    !! Loop index.
    integer :: ibuff
    !! Buffer for swapping elements.
    logical :: reverse_
    !! Boolean whether to sort in reverse order.

    if(present(reverse))then
       reverse_=reverse
    else
       reverse_=.false.
    end if

    dim=size(arr1,dim=1)
    do i = 1, dim
       if(reverse_)then
          loc=maxloc(arr1(i:dim),dim=1)+i-1          
       else
          loc=minloc(arr1(i:dim),dim=1)+i-1
       end if
       ibuff=arr1(i)
       arr1(i)=arr1(loc)
       arr1(loc)=ibuff

       if(present(arr2)) then
          ibuff=arr2(i)
          arr2(i)=arr2(loc)
          arr2(loc)=ibuff
       end if
    end do

    return
  end subroutine isort1D
!###############################################################################


!###############################################################################
  subroutine rsort1D(arr1,arr2,reverse)
    !! Sort a 1D real array from min to max.
    implicit none

    ! Arguments
    real(real32), dimension(:), intent(inout) :: arr1
    !! Array to be sorted.
    integer, dimension(:),intent(inout),optional :: arr2
    !! Optional. Second array to be sorted.
    logical, optional, intent(in) :: reverse
    !! Optional. Boolean whether to sort in reverse order.

    ! Local variables
    integer :: i,dim,loc
    !! Loop index.
    integer :: ibuff
    !! Buffer for swapping elements.
    real(real32) :: rbuff
    !! Buffer for swapping elements.
    logical :: reverse_
    !! Boolean whether to sort in reverse order.

    if(present(reverse))then
       reverse_=reverse
    else
       reverse_=.false.
    end if

    dim=size(arr1,dim=1)
    do i = 1, dim
       select case(reverse_)
       case(.true.)
          loc=maxloc(arr1(i:dim),dim=1)+i-1          
       case default
          loc=minloc(arr1(i:dim),dim=1)+i-1
       end select
       rbuff     = arr1(i)
       arr1(i)   = arr1(loc)
       arr1(loc) = rbuff

       if(present(arr2)) then
          ibuff     = arr2(i)
          arr2(i)   = arr2(loc)
          arr2(loc) = ibuff
       end if
    end do

    return
  end subroutine rsort1D
!###############################################################################


!###############################################################################
  pure recursive subroutine quicksort(arr, low, high)
    !! Sort a 1D real array from min to max.
    !!
    !! This is a recursive implementation of the quicksort algorithm.
    implicit none

    ! Arguments
    real(real32), dimension(:), intent(inout) :: arr
    !! Array to be sorted.
    integer, intent(in) :: low, high
    !! Lower and upper bounds of the array to be sorted.

    ! Local variables
    integer :: i, j
    !! Loop indices.
    real(real32) :: pivot, temp
    !! Pivot element and temporary buffer.

    if (low .lt. high) then
       pivot = arr((low + high) / 2)
       i = low
       j = high
       do
          do while (arr(i) .lt. pivot .and. i .lt. high)
             i = i + 1
          end do
          do while (arr(j) .gt. pivot .and. j .gt. low)
             j = j - 1
          end do
          if (i .le. j) then
             temp = arr(i)
             arr(i) = arr(j)
             arr(j) = temp
             i = i + 1
             j = j - 1
          end if
          ! Exit the loop when indices cross
          if (i .gt. j) exit
       end do
       ! Recursively apply quicksort to both partitions
       if (low .lt. j) call quicksort(arr, low, j)
       if (i .lt. high) call quicksort(arr, i, high)
    end if
  end subroutine quicksort
!###############################################################################


!###############################################################################
  subroutine sort2D(arr,dim)
    !! Sort a 2D array along the first column.
    implicit none

    ! Arguments
    integer, intent(in) :: dim
    !! Dimension to sort along.
    real(real32), dimension(dim,3), intent(inout) :: arr
    !! Array to be sorted.

    ! Local variables
    integer :: i,j,loc,istart
    !! Loop indices.
    integer, dimension(3) :: a123
    !! Array to store the order of sorting.
    real(real32), dimension(3) :: buff
    !! Buffer for swapping elements.

    a123(:)=(/1,2,3/)
    istart=1
    do j = 1, 3
       do i = j, dim
          loc = minloc( &
               abs(arr(i:dim,a123(1))), &
               dim = 1, &
               mask = (abs(arr(i:dim,a123(1))).gt.1.D-5) &
          ) + i - 1
          buff(:) = arr(i,:)
          arr(i,:) = arr(loc,:)
          arr(loc,:) = buff(:)
       end do

       scndrow: do i = j, dim
          if(abs(arr(j,a123(1))).ne.abs(arr(i,a123(1)))) exit scndrow
          loc = minloc( &
               abs( arr(i:dim,a123(2)) ) + abs( arr(i:dim,a123(3) ) ), &
               dim = 1, &
               mask = ( &
                    abs( arr(j,a123(1)) ) .eq. abs( arr(i:dim,a123(1)) ) &
               ) &
          ) + i - 1
          buff(:)=arr(i,:)
          arr(i,:)=arr(loc,:)
          arr(loc,:)=buff(:)
       end do scndrow

       a123=cshift(a123,1)
    end do

    return
  end subroutine sort2D
!###############################################################################


!###############################################################################
  subroutine iset(arr)
    !! Reduce an integer array to its unique elements.
    implicit none

    ! Arguments
    integer, dimension(:), allocatable, intent(inout) :: arr
    !! Array to be reduced.

    ! Local variables
    integer :: i,n
    !! Loop index.
    integer, dimension(:), allocatable :: tmp_arr
    !! Temporary array for storing unique elements.
    

    call sort1D(arr)
    allocate(tmp_arr(size(arr)))

    tmp_arr(1) = arr(1)
    n=1
    do i = 2, size(arr)
       if(arr(i)==tmp_arr(n)) cycle
       n = n + 1
       tmp_arr(n) = arr(i)
    end do
    deallocate(arr); allocate(arr(n))
    arr(:n) = tmp_arr(:n)
    !call move_alloc(tmp_arr, arr)
    
  end subroutine iset
!###############################################################################


!###############################################################################
  subroutine rset(arr, tol, count_list)
    !! Reduce a real array to its unique elements.
    implicit none

    ! Arguments
    real(real32), dimension(:), allocatable, intent(inout) :: arr
    !! Array to be reduced.
    real(real32), intent(in), optional :: tol
    !! Tolerance for comparing real numbers.
    integer, dimension(:), allocatable, intent(out), optional :: count_list
    !! List of counts for each unique element.

    ! Local variables
    integer :: i,n
    !! Loop index.
    real(real32) :: tol_
    !! Tolerance for comparing real numbers.
    real(real32), dimension(:), allocatable :: tmp_arr
    !! Temporary array for storing unique elements.
    integer, dimension(:), allocatable :: count_list_
    !! List of counts for each unique element.
    

    if(present(tol))then
       tol_ = tol
    else
       tol_ = 1.E-4_real32
    end if
    
    call quicksort(arr, 1, size(arr))
    allocate(tmp_arr(size(arr)))
    allocate(count_list_(size(arr)), source = 1)

    tmp_arr(1) = arr(1)
    n=1
    do i = 2, size(arr)
       if(abs(arr(i)-tmp_arr(n)).lt.tol_)then
          count_list_(i) = count_list_(i) + 1
          cycle
       end if
       n = n + 1
       tmp_arr(n) = arr(i)
    end do
    deallocate(arr); allocate(arr(n))
    arr(:n) = tmp_arr(:n)
    if(present(count_list)) count_list = count_list_(:n)
    
  end subroutine rset
!###############################################################################


!###############################################################################
  subroutine cset(arr,lcase,lkeep_size)
    !! Reduce a character array to its unique elements.
    !!
    !! This subroutine reduces a character array to its unique elements.
    !! i.e. each string in the array is compared with the rest of the strings
    !! in the array and if a match is found, the string is removed.
    !! This results in only the unique strings being preserved.
    implicit none

    ! Arguments
    character(*), allocatable, dimension(:), intent(inout) :: arr
    !! Array to be reduced.
    logical, intent(in), optional :: lcase
    !! Optional. Boolean whether to perform case insensitive comparison.
    logical, intent(in), optional :: lkeep_size
    !! Optional. Boolean whether to keep the original size of the array.

    ! Local variables
    integer :: i, n
    !! Loop index.
    logical :: lkeep_size_
    !! Boolean whether to keep the original size of the array.
    character(len=:), allocatable, dimension(:) :: tmp_arr
    !! Temporary array for storing unique elements.
    logical :: lcase_
    !! Boolean whether to perform case insensitive comparison.


    if(present(lcase))then
       lcase_ = lcase
    else
       lcase_ = .false.
    end if
    call sort_str(arr,lcase_)
    
    allocate(character(len=len(arr(1))) :: tmp_arr(size(arr)))
    tmp_arr(1) = arr(1)
    n=1
    
    do i = 2, size(arr)
       if(lcase_) arr(i) = to_lower(arr(i))
       if(trim(arr(i)).eq.trim(tmp_arr(n))) cycle
       n = n + 1
       tmp_arr(n) = arr(i)
    end do
    if(present(lkeep_size))then
       lkeep_size_=lkeep_size
    else
       lkeep_size_=.false.
    end if

    if(lkeep_size_)then
       call move_alloc(tmp_arr,arr)
    else
       deallocate(arr)
       allocate(arr(n))
       arr(:n) = tmp_arr(:n)
    end if

  end subroutine cset
!###############################################################################


!###############################################################################
  subroutine iswap(value1,value2)
    !! Swap two integers.
    implicit none

    ! Arguments
    integer, intent(inout) :: value1, value2
    !! Integers to be swapped.

    ! Local variables
    integer :: itmp1
    !! Temporary buffer for swapping elements.

    itmp1  = value1
    value1 = value2
    value2 = itmp1
  end subroutine iswap
!###############################################################################


!###############################################################################
  subroutine rswap(value1,value2)
    !! Swap two reals.
    implicit none

    ! Arguments
    real(real32), intent(inout) :: value1, value2
    !! Reals to be swapped.

    ! Local variables
    real(real32) :: rtmp1
    !! Temporary buffer for swapping elements.

    rtmp1  = value1
    value1 = value2
    value2 = rtmp1
  end subroutine rswap
!###############################################################################


!###############################################################################
  subroutine cswap(c1,c2)
    !! Swap two character strings.
    implicit none

    ! Arguments
    character(*), intent(inout) :: c1, c2
    !! Strings to be swapped.

    ! Local variables
    character(len=:), allocatable :: ctmp
    !! Temporary buffer for swapping elements.

    ctmp=c1
    c1=c2
    c2=ctmp
  end subroutine cswap
!###############################################################################


!###############################################################################
  subroutine rswap_vec(vec1,vec2)
    !! Swap two real vectors.
    implicit none

    ! Arguments
    real(real32),dimension(:), intent(inout) :: vec1, vec2
    !! Vectors to be swapped.

    ! Local variables
    real(real32),allocatable,dimension(:)::tvec
    !! Temporary buffer for swapping elements.

    allocate(tvec(size(vec1)))
    tvec=vec1(:)
    vec1(:)=vec2(:)
    vec2(:)=tvec
  end subroutine rswap_vec
!###############################################################################


!###############################################################################
  subroutine ishuffle(arr,dim,seed)
    !! Shuffle a 2D integer array.
    implicit none

    ! Arguments
    integer, dimension(:,:), intent(inout) :: arr
    !! Array to be shuffled.
    integer, intent(in) :: dim
    !! Dimension to shuffle along.
    integer, intent(in), optional :: seed
    !! Seed for random number generator.

    ! Local variables
    integer :: iseed
    !! Seed for random number generator.
    integer :: i, j, k, n_data, iother, istart
    !! Loop indices.
    integer :: i1s,i2s,i1e,i2e,j1s,j2s,j1e,j2e
    !! Indices for swapping elements.
    real(real32) :: r
    !! Random number for shuffling.
    integer, allocatable, dimension(:,:) :: tlist
    !! Temporary list for swapping elements.


    if(present(seed)) iseed = seed

    call random_seed(iseed)
    n_data = size(arr,dim=dim)
    if(dim.eq.1)then
       iother = 2
       i2s=1;i2e=size(arr,dim=iother)
       j2s=1;j2e=size(arr,dim=iother)
       allocate(tlist(1,size(arr,dim=iother)))
    else
       iother = 1
       i1s=1;i1e=size(arr,dim=iother)
       j1s=1;j1e=size(arr,dim=iother)
       allocate(tlist(size(arr,dim=iother),1))
    end if
    istart=1
    do k = 1, 2
       do i = 1, n_data
          call random_number(r)
          j = istart + floor((n_data+1-istart)*r)
          if(dim.eq.1)then
             i1s=i;i1e=i
             j1s=j;j1e=j
          else
             i2s=i;i2e=i
             j2s=j;j2e=j
          end if
          tlist(:,:) = arr(i1s:i1e,i2s:i2e)
          arr(i1s:i1e,i2s:i2e) = arr(j1s:j1e,j2s:j2e)
          arr(j1s:j1e,j2s:j2e) = tlist(:,:)
       end do
    end do

  end subroutine ishuffle
!###############################################################################


!###############################################################################
  subroutine rshuffle(arr,dim,seed)
    !! Shuffle a 2D real array.
    implicit none

    ! Arguments
    real(real32), dimension(:,:), intent(inout) :: arr
    !! Array to be shuffled.
    integer, intent(in) :: dim
    !! Dimension to shuffle along.
    integer, intent(in), optional :: seed
    !! Seed for random number generator.

    ! Local variables
    integer :: iseed
    !! Seed for random number generator.
    integer :: i, j, k, n_data, iother, istart
    !! Loop indices.
    integer :: i1s,i2s,i1e,i2e,j1s,j2s,j1e,j2e
    !! Indices for swapping elements.
    real(real32) :: r
    !! Random number for shuffling.
    real(real32), allocatable, dimension(:,:) :: tlist
    !! Temporary list for swapping elements.


    if(present(seed)) iseed = seed

    call random_seed(iseed)
    n_data = size(arr,dim=dim)
    if(dim.eq.1)then
       iother = 2
       i2s=1;i2e=size(arr,dim=iother)
       j2s=1;j2e=size(arr,dim=iother)
    else
       iother = 1
       i1s=1;i1e=size(arr,dim=iother)
       j1s=1;j1e=size(arr,dim=iother)
    end if
    istart=1
    allocate(tlist(1,size(arr,dim=iother)))
    do k = 1, 2
       do i = 1, n_data
          call random_number(r)
          j = istart + floor((n_data+1-istart)*r)
          if(dim.eq.1)then
             i1s=i;i1e=i
             j1s=j;j1e=j
          else
             i2s=i;i2e=i
             j2s=j;j2e=j
          end if
          tlist(1:1,:) = arr(i1s:i1e,i2s:i2e)
          arr(i1s:i1e,i2s:i2e) = arr(j1s:j1e,j2s:j2e)
          arr(j1s:j1e,j2s:j2e) = tlist(1:1,:)
       end do
    end do

  end subroutine rshuffle
!###############################################################################


!###############################################################################
  integer function icount(line,fs)
    !! Count the number of fields separated by specified delimiter.
    !!
    !! This function counts the number of fields separated by a specified
    !! delimiter in a string. The default delimiter is a space.
    implicit none

    ! Arguments
    character(*) :: line
    !! String to be counted.
    character(*), intent(in), optional :: fs
    !! Optional. Delimiter (aka field separator).

    ! Local variables
    integer :: k
    !! Loop index.
    integer :: items, pos, length
    !! Number of fields and position in the string.
    character(len=:), allocatable :: fs_
    !! Delimiter (aka field separator).


    items=0
    pos=1
    length=1
    if(present(fs)) length=len(trim(fs))
    allocate(character(len=length) :: fs_)
    if(present(fs)) then
       fs_=trim(fs)
    else
       fs_=" "
    end if

    loop: do
       k=verify(line(pos:),fs_)
       if (k.eq.0) exit loop
       items=items+1
       pos=k+pos-1
       k=scan(line(pos:),fs_)
       if (k.eq.0) exit loop
       pos=k+pos-1
    end do loop
    icount=items

  end function icount
!###############################################################################


!###############################################################################
  subroutine grep(unit,input,lstart,lline,success)
    !! Search a file for a pattern.
    !!
    !! This subroutine searches a file for a pattern. It can search for the
    !! first line that contains the pattern or for the first line that starts
    !! with the pattern.
    implicit none

    ! Arguments
    integer :: unit
    !! Unit number of the file.
    character(*) :: input
    !! Pattern to search for.
    logical, intent(in), optional :: lstart
    !! Optional. Boolean whether to rewind file.
    logical, intent(in), optional :: lline
    !! Optional. Boolean whether the pattern is at the start of the line.
    logical, intent(out), optional :: success
    !! Optional. Boolean whether the pattern is found.

    ! Local variables
    integer :: iostat
    !! I/O status.
    character(1024) :: buffer
    !! Buffer for reading lines.
    logical :: lline_
    !! Boolean whether the pattern is at the start of the line.
    logical :: success_
    !! Boolean whether the pattern is found.


    lline_ = .false.
    success_ = .false.
    if(present(lstart))then
       if(lstart) rewind(unit)
    else
       rewind(unit)
    end if

    if(present(lline)) lline_ = lline
    if(lline_)then
       wholeloop: do
          read(unit,'(A100)',iostat=iostat) buffer
          if(is_iostat_end(iostat))then
             exit wholeloop
          elseif(iostat.ne.0)then
             call stop_program('I/O stat error encounted when reading file')
          end if
          if(index(trim(buffer),trim(input)).eq.1)then
             success_ = .true.
             exit wholeloop
          end if
       end do wholeloop
    else
       greploop: do
          read(unit,'(A100)',iostat=iostat) buffer
          if(is_iostat_end(iostat))then
             exit greploop
          elseif(iostat.ne.0)then
             call stop_program('I/O stat error encounted when reading file')
          end if
          if(index(trim(buffer),trim(input)).ne.0)then
             success_ = .true.
             exit greploop
          end if
       end do greploop
    end if

    if(present(success)) success = success_
  end subroutine grep
!###############################################################################


!###############################################################################
  subroutine flagmaker(buffer,flag,i,skip,empty)
    !! Assign variables of flags from get_command_argument.
    implicit none

    ! Arguments
    character(*), intent(inout) :: buffer
    !! Buffer to be assigned a flag.
    character(*), intent(in) :: flag
    !! Flag to look for.
    integer :: i
    !! Index of command argument.
    logical :: skip
    !! Boolean whether to skip the next argument.
    logical, intent(out) :: empty
    !! Boolean whether the buffer is empty.


    skip = .false.
    empty = .false.
    if(len(trim(buffer)).eq.len(trim(flag))) then
       call get_command_argument(i+1,buffer)
       if(scan(buffer,'-').eq.1.or.buffer.eq.'') then
          buffer=""
          empty=.true.
       else
          skip=.true.
       end if
    else
       buffer=buffer(len(trim(flag))+1:)
    end if

  end subroutine flagmaker
!###############################################################################


!###############################################################################
  subroutine jump(unit,linenum)
    !! Go to a specific line in a file.
    implicit none

    ! Arguments
    integer :: unit
    !! Unit number of the file.
    integer :: linenum
    !! Line number to jump to.

    ! Local variables
    integer :: i
    !! Loop index.


    rewind(unit)
    do i = 1, linenum, 1
       read(unit,*)
    end do

  end subroutine jump
!###############################################################################


!###############################################################################
  subroutine file_check(unit,filename,action)
    !! Check if a file exists and open it.
    implicit none

    ! Arguments
    integer, intent(inout) :: unit
    !! Unit number of the file.
    character(*), intent(inout) :: filename
    !! Name of the file.
    character(len=20), optional, intent(in) :: action
    !! Optional. Action to be taken on the file.

    ! Local variables
    integer :: i
    !! Loop index.
    integer :: iostat
    !! I/O status.
    character(20) :: action_
    !! Action to be taken on the file.
    logical :: filefound
    !! Boolean whether the file is found.


    action_="READWRITE"
    if(present(action)) action_=action
    action_=to_upper(action_)
    do i = 1, 5
       inquire(file=trim(filename),exist=filefound)
       if(.not.filefound) then
          write(6,'("File name ",A," not found.")')&
               "'"//trim(filename)//"'"
          write(6,'("Supply another filename: ")')
          read(*,*) filename
       else
          write(6,'("Using file ",A)')  &
               "'"//trim(filename)//"'"
          exit
       end if
       if(i.ge.4) then
          stop "Nope"
       end if
    end do
    if(trim(adjustl(action_)).eq.'NONE')then
       write(6,*) "File found, but not opened."
    else
       open(newunit=unit,file=trim(filename),&
            action=trim(action_),iostat=iostat)
    end if

  end subroutine file_check
!###############################################################################


!###############################################################################
  subroutine touch(file)
    !! Create a directory if it does not exist.
    implicit none

    ! Arguments
    character(*), intent(in) :: file
    !! Directory to be created.

    ! Local variables
    logical :: exists
    !! Boolean whether the directory exists.
  
    inquire(file=file, exist=exists)
    if(.not.exists) call execute_command_line("mkdir -p "//file)
  end subroutine touch
!###############################################################################


!###############################################################################
  function to_upper(buffer) result(upper)
    !! Convert a string to upper case.
    implicit none

    ! Arguments
    character(*), intent(in) :: buffer
    !! String to be converted to upper case.
    character(len=:),allocatable :: upper
    !! Upper case string.

    ! Local variables
    integer :: i,j
    !! Loop index.


    allocate(character(len=len(buffer)) :: upper)
    do i = 1, len(buffer)
       j=iachar(buffer(i:i))
       if(j.ge.iachar("a").and.j.le.iachar("z"))then
          upper(i:i)=achar(j-32)
       else
          upper(i:i)=buffer(i:i)
       end if
    end do

  end function to_upper
!###############################################################################


!###############################################################################
  function to_lower(buffer) result(lower)
    !! Convert a string to lower case.
    implicit none

    ! Arguments
    character(*), intent(in) :: buffer
    !! String to be converted to lower case.
    character(len=:), allocatable :: lower
    !! Lower case string.

    ! Local variables
    integer :: i,j
    !! Loop index.


    allocate(character(len=len(buffer)) :: lower)
    do i = 1, len(buffer)
       j=iachar(buffer(i:i))
       if(j.ge.iachar("A").and.j.le.iachar("Z"))then
          lower(i:i)=achar(j+32)
       else
          lower(i:i)=buffer(i:i)
       end if
    end do

  end function to_lower
!###############################################################################


!###############################################################################
  function strip_null(buffer) result(stripped)
    !! Strip null characters from a string.
    !!
    !! This is meant for handling strings passed from Python, which gain
    !! null characters at the end. The procedure finds the first null
    !! character and truncates the string at that point.
    !! Null characters are represented by ASCII code 0.
    implicit none

    ! Arguments
    character(*), intent(in) :: buffer
    !! String to be stripped.
    character(len=len(buffer)) :: stripped
    !! Stripped string.

    ! Local variables
    integer :: i
    !! Loop index.

    stripped = ""
    do i = 1, len(buffer)
       if(iachar(buffer(i:i)).ne.0)then
          stripped(i:i)=buffer(i:i)
       else
          exit
       end if
    end do

  end function strip_null
!###############################################################################

end module raffle__misc
