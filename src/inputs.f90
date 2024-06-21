!!!#############################################################################
!!! Code written by Ned Thaddeus Taylor
!!! Code for RAFFLE
!!! Code part of the ARTEMIS group
!!!#############################################################################
!!! module defines all global variables
!!!#############################################################################
module inputs
  use misc_raffle, only: file_check,flagmaker, icount, to_lower
  use constants, only: real12, ierror, pi
  implicit none
  

  private

  public :: verbose
  public :: vdW, volvar
  public :: bins, vps_ratio
  public :: seed
  public :: num_structures, num_species, task
  public :: stoichiometry_list, element_list
  public :: filename_host
  public :: database_format, database_list
  public :: cutoff_min_list, cutoff_max_list, width_list, sigma_list

  public :: set_global_vars


  logical :: lseed

  integer :: verbose = 0
  integer :: seed !random seed
  integer :: num_structures ! number of structures to generate
  integer :: num_species   ! total number of species to add
  integer :: task ! task setting (defines the RAFFLE task)
  integer, allocatable, dimension(:) :: stoichiometry_list ! stoichiometry of species to add
  character(3), allocatable, dimension(:) :: element_list !species names to add

  integer :: vdW, volvar

  real(real12), dimension(3) :: cutoff_min_list, cutoff_max_list
  real(real12), dimension(3) :: width_list, sigma_list

  integer, dimension(3) :: bins
  integer, dimension(3) :: vps_ratio

  character(1024), dimension(:), allocatable :: database_list ! list of directories containing input database
  character(1024) :: database_format !format of input file (POSCAR, XYZ, etc.
  character(1024) :: filename_host !host structure filename


!!!updated  2023/06/16


contains
!!!#############################################################################
  subroutine set_global_vars()
    implicit none
    integer :: Reason
    integer :: i,j
    integer :: nseed
    character(1024) :: pattern,buffer,flag,input_file
    logical :: skip,empty,filefound
    integer, dimension(:), allocatable :: seed_arr 


!!!-----------------------------------------------------------------------------
!!! initialises variables
!!!-----------------------------------------------------------------------------
    input_file = ""
    seed = 1

!!!-----------------------------------------------------------------------------
!!! Reads flags and assigns to variables
!!!-----------------------------------------------------------------------------
    flagloop: do i=0,command_argument_count()-1
       empty=.false.
       if (skip) then
          skip=.false.
          cycle flagloop
       end if
       call get_command_argument(i,buffer)
       buffer=trim(buffer)
!!!------------------------------------------------------------------------
!!! FILE AND DIRECTORY FLAGS
!!!------------------------------------------------------------------------
       if(index(buffer,'-f').eq.1)then
          flag="-f"
          call flagmaker(buffer,flag,i,skip,empty)
          if(.not.empty)then
             read(buffer,'(A)') input_file
          else
             write(6,'("ERROR: No input filename supplied, but the flag ''-f'' was used")')
             infilename_do: do j=1,3
                write(6,'("Please supply an input filename:")')
                read(5,'(A)') input_file
                if(trim(input_file).ne.'')then
                   write(6,'("Input filename supplied")')
                   exit infilename_do
                else
                   write(6,'(1X,"Not a valid filename")')
                end if
                if(j.eq.3)then
                   write(0,*) "ERROR: No valid input filename supplied"
                   stop "Exiting..."
                end if
             end do infilename_do
          end if
!!!------------------------------------------------------------------------
!!! VERBOSE PRINTS
!!!------------------------------------------------------------------------
       elseif(index(buffer,'-v').eq.1)then
          flag="-v"
          call flagmaker(buffer,flag,i,skip,empty)
          if(.not.empty) read(buffer,*) ierror
       elseif(index(buffer,'-h').eq.1)then
          write(6,'("Flags:")')
          write(6,'(2X,"-h              : Prints the help for each flag.")')
          write(6,'(2X,"-v              : Verbose printing.")')
          write(6,'("-----------------FILE-NAME-FLAGS-----------------")')
          write(6,'(2X,"-f<STR>         : Input structure file name (Default = (empty)).")')
          stop
       end if
    end do flagloop


!!!-----------------------------------------------------------------------------
!!! check if input file was specified and read if true
!!!-----------------------------------------------------------------------------
    if(trim(input_file).ne."")then
       call read_input_file(input_file)
    end if


!!!-----------------------------------------------------------------------------
!!! initialise random seed
!!!-----------------------------------------------------------------------------
    call random_seed(size=nseed)
    allocate(seed_arr(nseed))
    if(lseed)then
       seed_arr = seed
    else
       call system_clock(count=seed)
       seed_arr = seed + 37 * (/ (i - 1, i = 1, nseed) /)
    end if
    call random_seed(put=seed_arr)
    deallocate(seed_arr)

    return
  end subroutine set_global_vars
!!!#############################################################################



!!!#############################################################################
!!! read input file to get variables
!!!#############################################################################
  subroutine read_input_file(file_name)
    implicit none
    character(*), intent(in) :: file_name

    integer :: i
    integer :: Reason,unit
    character(1) :: fs
    character(1024) :: stoichiometry, elements, database
    real(real12), dimension(3) :: width, sigma
    character(50), dimension(3) :: cutoff_min, cutoff_max


!!!-----------------------------------------------------------------------------
!!! set up namelists for input file
!!!-----------------------------------------------------------------------------
    namelist /setup/        task, filename_host, seed, vps_ratio, bins, &
                            database_format, database
    namelist /structure/    num_structures,num_species,elements,stoichiometry
    namelist /volume/       vdW, volvar
    namelist /distribution/ cutoff_min, cutoff_max, width, sigma


!!!-----------------------------------------------------------------------------
!!! check input file exists and open
!!!-----------------------------------------------------------------------------
    unit=20
    call file_check(unit,file_name)


    cutoff_min = "-1.0"
    cutoff_max = "-1.0"
    width = -1._real12
    sigma = -1._real12
    database_format = "vasprun.xml"
!!!-----------------------------------------------------------------------------
!!! read namelists from input file
!!!-----------------------------------------------------------------------------
    read(unit,NML=setup,iostat=Reason)
    if(Reason.ne.0)then
       write(0,*) "THERE WAS AN ERROR IN READING SETUP"
    end if
    read(unit,NML=structure,iostat=Reason)
    if(.not.is_iostat_end(Reason).and.Reason.ne.0)then
       stop "THERE WAS AN ERROR IN READING STRUCTURE SETTINGS"
    end if
    read(unit,NML=volume,iostat=Reason)
    if(.not.is_iostat_end(Reason).and.Reason.ne.0)then
       stop "THERE WAS AN ERROR IN READING VOLUME SETTINGS"
    end if
    read(unit,NML=distribution,iostat=Reason)
    if(.not.is_iostat_end(Reason).and.Reason.ne.0)then
       stop "THERE WAS AN ERROR IN READING DISTRIBUTION SETTINGS"
    end if

    if(trim(database).ne."")then
       allocate(database_list(icount(database)))
       read(database,*) database_list
    end if


    if(trim(stoichiometry).ne."")then
       allocate(stoichiometry_list(num_species))
       read(stoichiometry,*) stoichiometry_list
    end if

    if(trim(elements).ne."")then
       allocate(element_list(num_species))
       read(elements,*) element_list
    end if

    
    do i = 1, 3
       cutoff_min_list(i) = read_value_from_string(cutoff_min(i))
       cutoff_max_list(i) = read_value_from_string(cutoff_max(i))   
       write(*,*) "Cutoff: ",cutoff_min_list(i),cutoff_max_list(i)    
    end do

    width_list = width
    sigma_list = sigma

!!!-----------------------------------------------------------------------------
!!! close input file
!!!-----------------------------------------------------------------------------
    close(unit)
    write(*,*) "Input file read successfully."

    return
  end subroutine read_input_file
!!!#############################################################################

  function read_value_from_string(string) result(output)
    implicit none
    character(*), intent(in) :: string
    real(real12) :: output

    integer :: k, pos
    real(real12) :: variable, power
    character(:), allocatable :: string_
    character(12) :: numeric_set = "0123456789.-"

    pos = 1
    output = 1._real12
    variable = 0._real12
    power = 1._real12
    string_ = trim(to_lower(string))
    loop: do
       !! read until first non-numeric character
       !! read string up to k - 1 to variable (multiply)
       k = verify(string_(pos:len_trim(string_)),numeric_set)
       if (k.eq.0)then
          read(string_(pos:),*) variable
          output = output * variable ** power
          exit loop
       elseif(k.gt.1)then
          read(string_(pos:pos+k-2),*) variable
          output = output * variable ** power
       end if

       pos = pos + k - 1
       !! identify what the next character is (*, /, pi)
       !! if *, then change power factor to 1._real12
       !! if /, then change power factor to -1._real12
       !! if pi, then change power factor to 1._real12 and variable = pi
       !! if blank space, move pos to next non-space character and cycle
       !! if end of string, exit loop

       if (string_(pos:pos).eq."*")then
          power = 1._real12
          pos = pos + 1
       elseif(string_(pos:pos).eq."/")then
          power = -1._real12
          pos = pos + 1
       elseif(string_(pos:pos+1).eq."pi")then
          power = 1._real12
          output = output * pi ** power
          pos = pos + 2
       end if
       if(pos.gt.len_trim(string_)) exit loop
       pos = pos + verify(string_(pos:), " ") - 1
       if(pos.gt.len_trim(string_)) exit loop

    end do loop


    return
  end function read_value_from_string


end module inputs
