!!!#############################################################################
!!! Code written by Ned Thaddeus Taylor
!!! Code for RAFFLE
!!! Code part of the ARTEMIS group
!!!#############################################################################
!!! module defines all global variables
!!!#############################################################################
module inputs
  use misc, only: file_check,flagmaker, icount
  use constants, only: real12, ierror, pi
  implicit none
  logical :: lseed

  integer :: seed !random seed
  integer :: num_structures ! number of structures to generate
  integer :: num_species   ! total number of species to add
  integer :: task ! task setting (defines the RAFFLE task)
  integer, allocatable, dimension(:) :: stoichiometry_list ! stoichiometry of species to add
  character(3), allocatable, dimension(:) :: element_list !species names to add

  integer :: vdW    
  integer :: volvar 
  integer :: minbond
  integer :: maxbond

  real(real12), dimension(3) :: cutoff_min_list, cutoff_max_list
  real(real12), dimension(3) :: width_list, sigma_list

  integer, dimension(3) :: bins
  integer, dimension(3) :: vps_ratio

  character(1024), dimension(:), allocatable :: database_list ! list of directories containing input database
  character(1024) :: database_format !format of input file (POSCAR, XYZ, etc.
  character(1024) :: filename_host !host structure filename


!!!!task:
!!!! 0) Run RSS
!!!! 1) Regenerate DIst Files (WIP)
!!!! 2) Run HOST_RSS
!!!! 3) Test
!!!! 4) Sphere_Overlap
!!!! 5) Bondangle_test 
!!!! 6) Run evo (Should be run after any set created)
!!!! 7) Add new poscar  
!!!! 8) Run evo, but don't regen energies or evolve distributions (only reformat gaussians) 
!!!! 9) Run evo, don't get energies but do evolve distributions
  
  

  private

  public :: vdW, volvar, minbond, maxbond
  public :: bins, vps_ratio
  public :: seed
  public :: num_structures, num_species, task
  public :: stoichiometry_list, element_list
  public :: filename_host
  public :: database_format, database_list
  public :: cutoff_min_list, cutoff_max_list, width_list, sigma_list

  public :: set_global_vars


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
    integer :: Reason,unit
    character(1) :: fs
    character(1024) :: stoichiometry, elements, database
    real(real12), dimension(3) :: cutoff_min, cutoff_max, width, sigma

    character(*), intent(in) :: file_name

!!!-----------------------------------------------------------------------------
!!! set up namelists for input file
!!!-----------------------------------------------------------------------------
    namelist /setup/        task, filename_host, seed, vps_ratio, bins, &
                            database_format, database
    namelist /structure/    num_structures,num_species,elements,stoichiometry
    namelist /volume/       vdW, volvar, minbond, maxbond
    namelist /distribution/ cutoff_min, cutoff_max, width, sigma


!!!-----------------------------------------------------------------------------
!!! check input file exists and open
!!!-----------------------------------------------------------------------------
    unit=20
    call file_check(unit,file_name)


    cutoff_min = -1._real12
    cutoff_max = -1._real12
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

    cutoff_min_list = cutoff_min
    cutoff_max_list = cutoff_max
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

end module inputs
