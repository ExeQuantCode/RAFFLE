!!!#############################################################################
!!! Code written by Ned Thaddeus Taylor and Francis Huw Davies
!!! Code part of the ARTEMIS group (Hepplestone research group).
!!! Think Hepplestone, think HRG.
!!!#############################################################################
!!! module contains various miscellaneous functions and subroutines.
!!! module includes the following functions and subroutines:
!!! sort1D           (sort 1st col of array by size. Opt:sort 2nd array wrt 1st)
!!! sort2D           (sort 1st two columns of an array by size)
!!! sort_str         (sort a list of strings)
!!! set              (return the sorted set of unique elements)
!!! swap             (swap two variables around)
!!! shuffle          (randomly shuffle a 2D array along one dimension)
!!!##################
!!! Icount           (counts words on line)
!!! grep             (finds 1st line containing the pattern)
!!! flagmaker        (read flag inputs supplied and stores variable if present)
!!! jump             (moves file to specified line number)
!!! file_check       (checks whether file exists and prompts user otherwise)
!!! touch            (creates a file if it doesn't exist)
!!! to_upper         (converts all characters in string to upper case)
!!! to_lower         (converts all characters in string to lower case)
!!! strip_null       (removes null characters from a string)
!!!#############################################################################
module misc_raffle
  use constants, only: real12
  implicit none


  private

  public :: sort1D, sort2D, sort_str, sort_str_order
  public :: set
  public :: shuffle
  public :: Icount, grep, flagmaker
  public :: jump, file_check, touch, to_upper, to_lower
  public :: strip_null


  interface sort1D
     procedure isort1D,rsort1D
  end interface sort1D

  interface set
     procedure iset,rset, cset
  end interface set

  interface swap
     procedure iswap, rswap, rswap_vec, cswap
  end interface swap

  interface shuffle
     procedure ishuffle, rshuffle
  end interface shuffle


!!!updated 2024/07/19


contains

!!!#####################################################
!!! sorts a character list
!!!#####################################################
  subroutine sort_str(list,lcase)
    implicit none
    integer :: i,loc
    integer :: charlen
    logical :: ludef_case
    character(*), dimension(:), intent(inout) :: list
    character(:), allocatable, dimension(:) :: tlist
    logical, optional, intent(in) :: lcase !default is false

    charlen = len(list(1))
    if(present(lcase))then
       if(lcase)then
          ludef_case = lcase
          allocate(character(len=charlen) :: tlist(size(list)))
          tlist = list
          do i=1,size(tlist)
             list(i) = to_upper(list(i))
          end do
       end if
    else
       ludef_case = .false.
    end if
    do i=1,size(list)
       loc = minloc(list(i:),dim=1)
       if(loc.eq.1) cycle
       if(ludef_case) call cswap(tlist(i),tlist(loc+i-1))
       call cswap(list(i),list(loc+i-1))
    end do
    if(ludef_case) list=tlist
    
    return
  end subroutine sort_str
!!!-----------------------------------------------------
!!!-----------------------------------------------------
  function sort_str_order(list,lcase) result(order)
    implicit none
    integer :: i,loc
    integer :: charlen
    logical :: ludef_case
    character(*), dimension(:), intent(inout) :: list
    character(:), allocatable, dimension(:) :: tlist
    logical, optional, intent(in) :: lcase !default is false

    integer, allocatable, dimension(:) :: torder,order

    charlen = len(list(1))
    if(present(lcase))then
       if(lcase)then
          ludef_case = lcase
          allocate(character(len=charlen) :: tlist(size(list)))
          tlist = list
          do i=1,size(tlist)
             list(i) = to_upper(list(i))
          end do
       end if
    else
       ludef_case = .false.
    end if

    allocate(torder(size(list)))
    do i=1,size(list)
       torder(i) = i
    end do
    
    do i=1,size(list)
       loc = minloc(list(i:),dim=1)
       if(loc.eq.1) cycle
       if(ludef_case) call cswap(tlist(i),tlist(loc+i-1))
       call cswap(list(i),list(loc+i-1))
       call iswap(torder(i),torder(loc+i-1))
    end do
    
    allocate(order(size(list)))
    do i=1,size(list)
       order(i) = findloc(torder,i,dim=1)
    end do
    
    if(ludef_case) list=tlist
    
    return
  end function sort_str_order
!!!#####################################################


!!!#####################################################
!!! sorts two arrays from min to max
!!! sorts the optional second array wrt the first array
!!!#####################################################
  subroutine isort1D(arr1,arr2,reverse)
    implicit none
    integer :: i,dim,loc
    integer :: ibuff
    logical :: reverse_
    integer, dimension(:) :: arr1
    integer, dimension(:),intent(inout),optional :: arr2
    logical, optional, intent(in) :: reverse

    if(present(reverse))then
       reverse_=reverse
    else
       reverse_=.false.
    end if

    dim=size(arr1,dim=1)
    do i=1,dim
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
!!!-----------------------------------------------------
!!!-----------------------------------------------------
  subroutine rsort1D(arr1,arr2,reverse)
    implicit none
    integer :: i,dim,loc
    real(real12) :: rbuff
    logical :: reverse_
    real(real12), dimension(:) :: arr1
    integer, dimension(:),intent(inout),optional :: arr2
    logical, optional, intent(in) :: reverse

    if(present(reverse))then
       reverse_=reverse
    else
       reverse_=.false.
    end if

    dim=size(arr1,dim=1)
    do i=1,dim
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
          rbuff=arr2(i)
          arr2(i)=arr2(loc)
          arr2(loc)=rbuff
       end if
    end do

    return
  end subroutine rsort1D
!!!#####################################################


!!!#####################################################
!!! quicksort an array from min to max
!!!#####################################################
  pure recursive subroutine quicksort(arr, low, high)
    implicit none
    real(real12), dimension(:), intent(inout) :: arr
    integer, intent(in) :: low, high
    integer :: i, j
    real(real12) :: pivot, temp

    if (low .lt. high) then
        pivot = arr((low + high) / 2)
        i = low
        j = high
        do
            do while (arr(i) .lt. pivot)
                i = i + 1
            end do
            do while (arr(j) .gt. pivot)
                j = j - 1
            end do
            if (i .le. j) then
                temp = arr(i)
                arr(i) = arr(j)
                arr(j) = temp
                i = i + 1
                j = j - 1
            end if
        end do
        call quicksort(arr, low, j)
        call quicksort(arr, i, high)
    end if
  end subroutine quicksort
!!!#####################################################


!!!#####################################################
!!! sort an array from min to max
!!!#####################################################
  subroutine sort2D(arr,dim)
    implicit none
    integer :: i,j,dim,loc,istart
    integer, dimension(3) :: a123
    real(real12), dimension(3) :: buff
    real(real12), dimension(dim,3) :: arr

    a123(:)=(/1,2,3/)
    istart=1
    do j=1,3
       do i=j,dim
          loc=minloc(abs(arr(i:dim,a123(1))),dim=1,mask=(abs(arr(i:dim,a123(1))).gt.1.D-5))+i-1
          buff(:)=arr(i,:)
          arr(i,:)=arr(loc,:)
          arr(loc,:)=buff(:)
       end do

       scndrow: do i=j,dim
          if(abs(arr(j,a123(1))).ne.abs(arr(i,a123(1)))) exit scndrow
          loc=minloc(abs(arr(i:dim,a123(2)))+abs(arr(i:dim,a123(3))),dim=1,&
               mask=(abs(arr(j,a123(1))).eq.abs(arr(i:dim,a123(1)))))+i-1
          buff(:)=arr(i,:)
          arr(i,:)=arr(loc,:)
          arr(loc,:)=buff(:)
       end do scndrow

       a123=cshift(a123,1)
    end do

    return
  end subroutine sort2D
!!!#####################################################


!!!#####################################################
!!! return the sorted set of unique elements
!!!#####################################################
  subroutine iset(arr)
    implicit none
    integer :: i,n
    integer, allocatable, dimension(:) :: tmp_arr
    
    integer, allocatable, dimension(:) :: arr

    call sort1D(arr)
    allocate(tmp_arr(size(arr)))

    tmp_arr(1) = arr(1)
    n=1
    do i=2,size(arr)
       if(arr(i)==tmp_arr(n)) cycle
       n = n + 1
       tmp_arr(n) = arr(i)
    end do
    deallocate(arr); allocate(arr(n))
    arr(:n) = tmp_arr(:n)
    !call move_alloc(tmp_arr, arr)
    
  end subroutine iset
!!!-----------------------------------------------------
!!!-----------------------------------------------------
  subroutine rset(arr, tol, count_list)
    implicit none
    integer :: i,n
    real(real12) :: tol_
    real(real12), allocatable, dimension(:) :: tmp_arr
    integer, allocatable, dimension(:) :: count_list_
    
    real(real12), allocatable, dimension(:) :: arr
    real(real12), optional :: tol
    integer, dimension(:), allocatable, optional :: count_list

    if(present(tol))then
       tol_ = tol
    else
       tol_ = 1.E-4_real12
    end if
    
    call quicksort(arr, 1, size(arr))
    allocate(tmp_arr(size(arr)))
    allocate(count_list_(size(arr)), source = 1)

    tmp_arr(1) = arr(1)
    n=1
    do i=2,size(arr)
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
!!!-----------------------------------------------------
!!!-----------------------------------------------------
  subroutine cset(arr,lcase,lkeep_size)
    implicit none
    integer :: i, n
    logical :: ludef_keep_size
    character(len=:), allocatable, dimension(:) :: tmp_arr
    character(*), allocatable, dimension(:) :: arr
    logical, intent(in), optional :: lcase, lkeep_size

    if(present(lcase))then
       call sort_str(arr,lcase)
    else
       call sort_str(arr)
    end if
    
    allocate(character(len=len(arr(1))) :: tmp_arr(size(arr)))
    tmp_arr(1) = arr(1)
    n=1
    
    do i=2,size(arr)
       if(trim(arr(i)).eq.trim(tmp_arr(n))) cycle
       n = n + 1
       tmp_arr(n) = arr(i)
    end do
    if(present(lkeep_size))then
       ludef_keep_size=lkeep_size
    else
       ludef_keep_size=.false.
    end if

    if(ludef_keep_size)then
       call move_alloc(tmp_arr,arr)!!!CONSISTENCY WITH OTHER SET FORMS
    else
       deallocate(arr)
       allocate(arr(n))
       arr(:n) = tmp_arr(:n)
    end if

  end subroutine cset
!!!#####################################################


!!!#####################################################
!!! sort an array over specified column
!!!#####################################################
!!! Have it optionally take in an integer vector that ...
!!! ... lists the order of imporance of columns
  subroutine sort_col(arr1,col,reverse)
    implicit none
    integer :: i,dim,loc
    logical :: reverse_
    real(real12), allocatable, dimension(:) :: dbuff
    real(real12), dimension(:,:) :: arr1

    integer, intent(in) :: col
    logical, optional, intent(in) :: reverse


    if(present(reverse))then
       reverse_=reverse
    else
       reverse_=.false.
    end if

    allocate(dbuff(size(arr1,dim=2)))

    dim=size(arr1,dim=1)
    do i=1,dim
       if(reverse_)then
          loc=maxloc(arr1(i:dim,col),dim=1)+i-1          
       else
          loc=minloc(arr1(i:dim,col),dim=1)+i-1
       end if
       dbuff=arr1(i,:)
       arr1(i,:)=arr1(loc,:)
       arr1(loc,:)=dbuff

    end do

    return
  end subroutine sort_col
!!!#####################################################


!!!#####################################################
!!! swap two ints
!!!#####################################################
  subroutine iswap(i1,i2)
    implicit none
    integer :: i1,i2,itmp

    itmp=i1
    i1=i2
    i2=itmp
  end subroutine iswap
!!!-----------------------------------------------------
!!!-----------------------------------------------------
  subroutine rswap(d1,d2)
    implicit none
    real(real12) :: d1,d2,dtmp

    dtmp=d1
    d1=d2
    d2=dtmp
  end subroutine rswap
!!!-----------------------------------------------------
!!!-----------------------------------------------------
  subroutine cswap(c1,c2)
    implicit none
    character(*) :: c1,c2
    character(len=:), allocatable :: ctmp

    ctmp=c1
    c1=c2
    c2=ctmp
  end subroutine cswap
!!!-----------------------------------------------------
!!!-----------------------------------------------------
  subroutine rswap_vec(vec1,vec2)
    implicit none
    real(real12),dimension(:)::vec1,vec2
    real(real12),allocatable,dimension(:)::tvec

    allocate(tvec(size(vec1)))
    tvec=vec1(:)
    vec1(:)=vec2(:)
    vec2(:)=tvec
  end subroutine rswap_vec
!!!#####################################################


!!!#####################################################
!!! shuffle an array along one dimension
!!!#####################################################
  subroutine ishuffle(arr,dim,seed)
   implicit none
   integer :: iseed,istart
   integer :: i,j,k,n_data,iother
   integer :: i1s,i2s,i1e,i2e,j1s,j2s,j1e,j2e
   real(real12) :: r
   integer, allocatable, dimension(:,:) :: tlist

   integer, intent(in) :: dim
   integer, dimension(:,:), intent(inout) :: arr

   integer, optional, intent(in) :: seed

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
   do k=1,2
      do i=1,n_data
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

 end subroutine ishuffle
!!!-----------------------------------------------------
!!!-----------------------------------------------------
  subroutine rshuffle(arr,dim,seed)
    implicit none
    integer :: iseed,istart
    integer :: i,j,k,n_data,iother
    integer :: i1s,i2s,i1e,i2e,j1s,j2s,j1e,j2e
    real(real12) :: r
    real(real12), allocatable, dimension(:,:) :: tlist

    integer, intent(in) :: dim
    real(real12), dimension(:,:), intent(inout) :: arr

    integer, optional, intent(in) :: seed

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
    do k=1,2
       do i=1,n_data
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
!!!#####################################################


!!!#############################################################################
!!!#############################################################################
!!!  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
!!!#############################################################################
!!!#############################################################################


!!!#####################################################
!!! counts the number of words on a line
!!!#####################################################
  integer function Icount(full_line,tmpchar)
    character(*) :: full_line
    !ONLY WORKS WITH IFORT COMPILER
    !      character(1) :: fs
    character(len=:),allocatable :: fs
    character(*),optional :: tmpchar
    integer ::items,pos,k,length
    items=0
    pos=1

    length=1
    if(present(tmpchar)) length=len(trim(tmpchar))
    allocate(character(len=length) :: fs)
    if(present(tmpchar)) then
       fs=trim(tmpchar)
    else
       fs=" "
    end if

    loop: do
       k=verify(full_line(pos:),fs)
       if (k.eq.0) exit loop
       items=items+1
       pos=k+pos-1
       k=scan(full_line(pos:),fs)
       if (k.eq.0) exit loop
       pos=k+pos-1
    end do loop
    Icount=items
  end function Icount
!!!#####################################################


!!!#####################################################
!!! grep 
!!!#####################################################
!!! searches a file untill it finds the mattching patern
  subroutine grep(unit,input,lstart,lline,success)
    integer :: unit,Reason
    character(*) :: input
    character(1024) :: buffer
    logical :: lline_ = .false., success_ = .false.
    logical, optional, intent(out) :: success
    logical, optional, intent(in) :: lstart, lline
    !  character(1024), intent(out), optional :: linechar
    if(present(lstart))then
       if(lstart) rewind(unit)
    else
       rewind(unit)
    end if

    if(present(lline)) lline_ = lline
    if(lline_)then
       wholeloop: do
          read(unit,'(A100)',iostat=Reason) buffer
          if(is_iostat_end(Reason))then
             exit wholeloop
          elseif(Reason.ne.0)then
             write(0,*) "Error reading file"
             stop 1
          end if
          if(index(trim(buffer),trim(input)).eq.1)then
             success_ = .true.
             exit wholeloop
          end if
       end do wholeloop
    else
       greploop: do
          read(unit,'(A100)',iostat=Reason) buffer
          if(is_iostat_end(Reason))then
             exit greploop
          elseif(Reason.ne.0)then
             write(0,*) "Error reading file"
             stop 1
          end if
          if(index(trim(buffer),trim(input)).ne.0)then
            success_ = .true.
            exit greploop
         end if
       end do greploop
    end if

    if(present(success)) success = success_
  end subroutine grep
!!!#####################################################


!!!#####################################################
!!! Assigns variables of flags from getarg
!!!#####################################################
!!! SHOULD MAKE THIS A FUNCTION INSTEAD !!!
  subroutine flagmaker(buffer,flag,i,skip,empty)
    integer :: i
    logical :: skip,empty
    character(*) :: flag,buffer

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

    return
  end subroutine flagmaker
!!!#####################################################


!!!#####################################################
!!! Jumps UNIT to input line number
!!!#####################################################
  subroutine jump(unit,linenum)
    integer :: unit, linenum, move
    rewind(unit)
    do move = 1, linenum, 1
       read(unit,*)
    end do
    return
  end subroutine jump
!!!#####################################################


!!!#####################################################
!!! File checker
!!!#####################################################
  subroutine file_check(UNIT,FILENAME,ACTION)
    implicit none
    integer :: i,UNIT,Reason
    character(len=*) :: FILENAME
    character(20) :: udef_action
    character(20), optional :: ACTION
    logical :: filefound

    udef_action="READWRITE"
    if(present(ACTION)) udef_action=ACTION
    udef_action=to_upper(udef_action)
    do i=1,5
       inquire(file=trim(FILENAME),exist=filefound)
       if(.not.filefound) then
          write(6,'("File name ",A," not found.")')&
               "'"//trim(FILENAME)//"'"
          write(6,'("Supply another filename: ")')
          read(*,*) FILENAME
       else
          write(6,'("Using file ",A)')  &
               "'"//trim(FILENAME)//"'"
          exit
       end if
       if(i.ge.4) then
          stop "Nope"
       end if
    end do
    if(trim(adjustl(udef_action)).eq.'NONE')then
       write(6,*) "File found, but not opened."
    else
       open(newunit=UNIT,file=trim(FILENAME),&
            action=trim(udef_action),iostat=Reason)
    end if


    return
  end subroutine file_check
!!!#####################################################

!!!#####################################################
!!! create a file if it doesn't exist
!!!#####################################################
  subroutine touch(file)
    implicit none
    character(*), intent(in) :: file
    logical :: exists
  
    inquire(file=file, exist=exists)
    if(.not.exists) call execute_command_line("mkdir -p "//file)
  end subroutine touch
!!!#####################################################

!!!#####################################################
!!! converts all characters in string to upper case
!!!#####################################################
  function to_upper(buffer) result(upper)
    implicit none
    integer :: i,j
    character(*) :: buffer
    character(len=:),allocatable :: upper


    allocate(character(len=len(buffer)) :: upper)
    do i=1,len(buffer)
       j=iachar(buffer(i:i))
       if(j.ge.iachar("a").and.j.le.iachar("z"))then
          upper(i:i)=achar(j-32)
       else
          upper(i:i)=buffer(i:i)
       end if
    end do

    return
  end function to_upper
!!!#####################################################


!!!#####################################################
!!! converts all characters in string to lower case
!!!#####################################################
  function to_lower(buffer) result(lower)
    implicit none
    integer :: i,j
    character(*) :: buffer
    character(len=:),allocatable :: lower


    allocate(character(len=len(buffer)) :: lower)
    do i=1,len(buffer)
       j=iachar(buffer(i:i))
       if(j.ge.iachar("A").and.j.le.iachar("Z"))then
          lower(i:i)=achar(j+32)
       else
          lower(i:i)=buffer(i:i)
       end if
    end do

    return
  end function to_lower
!!!#####################################################


!!!#####################################################
!!! strip null characters from string
!!!#####################################################
  function strip_null(buffer) result(stripped)
    implicit none
    integer :: i
    character(*) :: buffer
    character(len=len(buffer)) :: stripped

    stripped = ""
    do i=1,len(buffer)
       if(iachar(buffer(i:i)).ne.0)then
          stripped(i:i)=buffer(i:i)
       else
          exit
       end if
    end do

    return
  end function strip_null
  !!!#####################################################

end module misc_raffle
