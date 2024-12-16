program test_tools_infile
  use raffle__constants, only: real32
  use raffle__tools_infile
  implicit none


  integer :: ival = 0
  logical :: ltmp1
  real(real32) :: rtmp1
  character(256) :: stmp1
  character(256) :: line

  logical :: success = .true.


  line = "APPLES = string"
  call assign_val(line, stmp1, ival, keyword="APPLES")
  if( trim(stmp1) .ne. "string" .or. ival .ne. 1 )then
     write(0,*) "assign_val failed for string"
     success = .false.
  end if

  line = "ORANGES = 1"
  call assign_val(line, ltmp1, ival, keyword="ORANGES")
  if( .not. ltmp1 .or. ival .ne. 2 )then
     write(0,*) "assign_val failed for logical"
     success = .false.
  end if
  line = "ORANGES = 0"
  call assign_val(line, ltmp1, ival, keyword="ORANGES")
  if( ltmp1 .or. ival .ne. 3 )then
     write(0,*) "assign_val failed for logical"
     success = .false.
  end if
  line = "ORANGES = T"
  call assign_val(line, ltmp1, ival, keyword="ORANGES")
  if( .not. ltmp1 .or. ival .ne. 4 )then
     write(0,*) "assign_val failed for logical"
     success = .false.
  end if
  line = "ORANGES = F"
  call assign_val(line, ltmp1, ival, keyword="ORANGES")
  if( ltmp1 .or. ival .ne. 5 )then
     write(0,*) "assign_val failed for logical"
     success = .false.
  end if
  line = "ORANGES = t"
  call assign_val(line, ltmp1, ival, keyword="ORANGES")
  if( .not. ltmp1 .or. ival .ne. 6 )then
     write(0,*) "assign_val failed for logical"
     success = .false.
  end if
  line = "ORANGES = f"
  call assign_val(line, ltmp1, ival, keyword="ORANGES")
  if( ltmp1 .or. ival .ne. 7 )then
     write(0,*) "assign_val failed for logical"
     success = .false.
  end if

  line = "BANANAS = 1.0 # comment"
  ! ival = line number here
  call rm_comments(line, ival)
  if( trim(line) .ne. "BANANAS = 1.0")then
     write(0,*) "rm_comments failed"
     write(0,'("\",A,"\")') trim(line)
     success = .false.
  end if


  !-----------------------------------------------------------------------------
  ! check for any failed tests
  !-----------------------------------------------------------------------------
  write(*,*) "----------------------------------------"
  if(success)then
     write(*,*) 'test_tools_infile passed all tests'
  else
     write(0,*) 'test_tools_infile failed one or more tests'
     stop 1
  end if



end program test_tools_infile