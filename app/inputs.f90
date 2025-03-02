module inputs
  !! Module for reading input files and setting global variables.
  !!
  !! This module reads the input file and sets the global variables for the
  !! program. NOTE: this module is not intended to be used in a library, it is
  !! only included in the executable program.
  use raffle__misc, only: file_check,flagmaker, icount, to_lower
  use raffle__generator, only: stoichiometry_type
  use raffle__constants, only: real32, pi
  use raffle__io_utils
  implicit none


  private

  public :: grid, grid_spacing, method_ratio
  public :: seed
  public :: num_structures, task
  public :: stoich
  public :: filename_host
  public :: database_format, database_list
  public :: cutoff_min_list, cutoff_max_list, width_list, sigma_list
  public :: output_dir
  public :: element_symbols, element_energies
  public :: bond_pairs, pair_radii

  public :: set_global_vars
  public :: verbose


  logical :: lseed
  !! Logical for random seed.

  integer :: verbose
  !! Verbose printing.
  integer :: seed
  !! Random seed.
  integer :: num_structures
  !! Number of structures to generate.
  integer :: task
  !! Task to perform (0 = run, 1 = continue).
  type(stoichiometry_type), dimension(:), allocatable :: stoich
  !! Stoichiometry of species to add

  real(real32), dimension(3) :: cutoff_min_list, cutoff_max_list
  !! Cutoff values for distribution functions.
  real(real32), dimension(3) :: width_list, sigma_list
  !! Width and sigma values for distribution functions.

  real(real32), dimension(:), allocatable :: element_energies
  !! Element reference energies.
  real(real32), dimension(:), allocatable :: pair_radii
  !! Element pair radii.
  character(3), dimension(:), allocatable :: element_symbols
  !! Element symbols.
  character(3), dimension(:,:), allocatable :: bond_pairs
  !! Element symbol pairs for pair_radii.

  integer, dimension(3) :: grid = [0, 0, 0]
  !! Grid dimensions.
  real(real32) :: grid_spacing = 0._real32
  !! Grid spacing.
  real(real32), dimension(5) :: method_ratio = &
       [1._real32, 0.1_real32, 0.5_real32, 0.5_real32, 1._real32]
  !! Placement method ratios.

  character(1024), dimension(:), allocatable :: database_list
  !! List of directories containing input database.
  character(1024) :: database_format
  !! Format of input file (POSCAR, XYZ, etc).
  character(1024) :: filename_host
  !! Host structure filename.
  character(1024) :: output_dir
  !! Output directory.


contains

!###############################################################################
  subroutine set_global_vars()
    !! Set global variables for the program.
    !!
    !! This subroutine reads the flags from the command line and sets the global
    !! variables for the program.
    implicit none

    ! Local variables
    integer :: i,j
    !! Loop indices.
    integer :: nseed
    !! Number of random seed values.
    character(1024) :: buffer, flag, input_file
    !! Character variables.
    logical :: skip, empty
    !! Logical variables.
    integer, dimension(:), allocatable :: seed_arr
    !! Random seed array.


    !---------------------------------------------------------------------------
    ! initialises variables
    !---------------------------------------------------------------------------
    input_file = ""
    seed = 1

    !---------------------------------------------------------------------------
    ! read flags and assign to variables
    !---------------------------------------------------------------------------
    flagloop: do i = 0, command_argument_count()
       empty=.false.
       if (skip) then
          skip=.false.
          cycle flagloop
       end if
       call get_command_argument(i,buffer)
       buffer=trim(buffer)
       !------------------------------------------------------------------------
       ! flag options
       !------------------------------------------------------------------------
       if(index(buffer,'-f').eq.1)then
          flag="-f"
          call flagmaker(buffer,flag,i,skip,empty)
          if(.not.empty)then
             read(buffer,'(A)') input_file
          else
             write(0,'("&
                  &ERROR: No input filename supplied, &
                  &but the flag ''-f'' was used&
             &")')
             infilename_do: do j = 1, 3
                write(6,'("Please supply an input filename:")')
                read(5,'(A)') input_file
                if(trim(input_file).ne.'')then
                   write(6,'("Input filename supplied")')
                   exit infilename_do
                else
                   write(0,'(1X,"Not a valid filename")')
                end if
                if(j.eq.3)then
                   call stop_program("No valid input filename supplied")
                   return
                end if
             end do infilename_do
          end if
       elseif(index(buffer,'-v').eq.1)then
          flag="-v"
          call print_build_info()
          stop
       elseif(index(buffer,'--verbose').eq.1)then
          flag="--verbose"
          call flagmaker(buffer,flag,i,skip,empty)
          if(.not.empty) read(buffer,*) verbose
       elseif(index(buffer,'-h').eq.1)then
          write(6,'("Flags:")')
          write(6,'(2X,"-h              : Prints the help for each flag.")')
          write(6,'(2X,"-v              : Verbose printing.")')
          write(6,'("-----------------FILE-NAME-FLAGS-----------------")')
          write(6,'(2X,"-f<STR>         : Input structure file name (&
               &Default = (empty)&
          &).")')
          stop
       end if
    end do flagloop


    !---------------------------------------------------------------------------
    ! check if input file was specified and read if true
    !---------------------------------------------------------------------------
    if(trim(input_file).ne."")then
       call read_input_file(input_file)
    end if


    !---------------------------------------------------------------------------
    ! initialise random seed
    !---------------------------------------------------------------------------
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
!###############################################################################


!###############################################################################
  subroutine read_input_file(file_name)
    !! Read input file to get variables.
    !!
    !! This subroutine reads a namelist input file to get the variables for the
    !! program.
    implicit none

    ! Arguments
    character(*), intent(inout) :: file_name
    !! Input file name.

    ! Local variables
    integer :: i
    !! Loop index.
    integer :: num_species, num_bonds
    !! Number of species and bonds.
    integer :: iostat, unit, l_pos, r_pos
    !! I/O status, unit, left and right positions.
    character(32) :: pair
    !! Pair of elements.
    character(1024) :: stoichiometry, database, buffer, energies, &
         bond_radii
    !! Strings buffers to hold input values (usually derived types).
    real(real32), dimension(3) :: width, sigma
    !! Width and sigma values for distribution functions.
    real(real32) :: void, rand, walk, grow, min
    !! Placement method ratios.
    character(50), dimension(3) :: cutoff_min, cutoff_max
    !! Cutoff values for distribution functions.


    !---------------------------------------------------------------------------
    ! set up namelists for input file
    !---------------------------------------------------------------------------
    namelist /setup/ task, filename_host, seed, grid, &
         grid_spacing, &
         database_format, database, verbose, output_dir
    namelist /placement_method/ void, rand, walk, grow, min
    namelist /structure/    num_structures,stoichiometry
    namelist /distribution/ cutoff_min, cutoff_max, width, sigma
    namelist /element_info/ energies, bond_radii


    !---------------------------------------------------------------------------
    ! check input file exists and open
    !---------------------------------------------------------------------------
    unit=20
    call file_check(unit,file_name)


    !---------------------------------------------------------------------------
    ! initialise variables
    !---------------------------------------------------------------------------
    output_dir = "iteration1"
    cutoff_min = "-1.0"
    cutoff_max = "-1.0"
    width = -1._real32
    sigma = -1._real32
    database_format = "vasprun.xml"
    void = 0._real32; rand = 0._real32
    walk = 0._real32; grow = 0._real32
    min = 0._real32


    !---------------------------------------------------------------------------
    ! read namelists from input file
    !---------------------------------------------------------------------------
    read(unit,NML=setup,iostat=iostat)
    if(iostat.ne.0)then
       write(0,*) "THERE WAS AN ERROR IN READING SETUP"
    end if
    read(unit,NML=placement_method,iostat=iostat)
    if(.not.is_iostat_end(iostat).and.iostat.ne.0)then
       call stop_program( &
            "THERE WAS AN ERROR IN READING PLACEMENT_METHOD SETTINGS" &
       )
       return
    end if
    read(unit,NML=structure,iostat=iostat)
    if(.not.is_iostat_end(iostat).and.iostat.ne.0)then
       call stop_program( &
            "THERE WAS AN ERROR IN READING STRUCTURE SETTINGS" &
       )
       return
    end if
    read(unit,NML=distribution,iostat=iostat)
    if(.not.is_iostat_end(iostat).and.iostat.ne.0)then
       call stop_program( &
            "THERE WAS AN ERROR IN READING DISTRIBUTION SETTINGS" &
       )
       return
    end if
    read(unit,NML=element_info,iostat=iostat)
    if(.not.is_iostat_end(iostat).and.iostat.ne.0)then
       call stop_program( &
            "THERE WAS AN ERROR IN READING ELEMENT_INFO SETTINGS" &
       )
       return
    end if


    !---------------------------------------------------------------------------
    ! handle variables
    !---------------------------------------------------------------------------
    if(trim(database).ne."")then
       allocate(database_list(icount(database)))
       read(database,*) database_list
       l_pos = 0
       do i = 1, size(database_list)
          read(database(l_pos+1:),'(A)') buffer
          r_pos = scan(buffer,",")
          if(r_pos.eq.0) r_pos = len_trim(buffer)
          read(buffer(:r_pos-1),'(A)') database_list(i)
          l_pos = scan(database(l_pos+1:),",") + l_pos
       end do
    end if

    method_ratio = [void, rand, walk, grow, min]
    if(all(abs(method_ratio).lt.1.E-6))then
       method_ratio = &
            [1._real32, 0.1_real32, 0.5_real32, 0.5_real32, 1._real32]
    end if

    if(trim(stoichiometry).ne."")then
       num_species = icount(stoichiometry,",")
       allocate(stoich(num_species))
       l_pos = scan(stoichiometry,"{")
       r_pos = scan(stoichiometry,"}", back=.true.)
       do i = 1, num_species
          read(stoichiometry(l_pos+1:r_pos-1),'(A)') buffer
          read(buffer(:scan(buffer,":")-1),*) stoich(i)%element
          read(buffer(scan(buffer,":")+1:),*) stoich(i)%num
          l_pos = scan(stoichiometry(l_pos+1:),",") + l_pos
       end do
    else
       call stop_program("No stoichiometry specified")
       return
    end if


    if(trim(energies).ne."")then
       num_species = icount(energies,",")
       allocate(element_symbols(num_species))
       allocate(element_energies(num_species))
       l_pos = scan(energies,"{")
       r_pos = scan(energies,"}", back=.true.)
       do i = 1, num_species
          read(energies(l_pos+1:r_pos-1),'(A)') buffer
          read(buffer(:scan(buffer,":")-1),*) element_symbols(i)
          read(buffer(scan(buffer,":")+1:),*) element_energies(i)
          l_pos = scan(energies(l_pos+1:),",") + l_pos
       end do
    else
       call stop_program("No element energies specified")
       return
    end if


    if(trim(bond_radii).ne."")then
       num_bonds = icount(bond_radii,",")
       allocate(bond_pairs(num_bonds,2))
       allocate(pair_radii(num_bonds))
       l_pos = scan(energies,"{")
       r_pos = scan(energies,"}", back=.true.)
       do i = 1, num_bonds
          read(bond_radii(l_pos+1:r_pos-1),'(A)') buffer
          read(buffer(:scan(buffer,":")-1),*) pair
          if(index(pair,"-").ne.0)then
             read(pair(:index(pair,"-")-1),*) bond_pairs(i,1)
             read(pair(index(pair,"-")+1:),*) bond_pairs(i,2)
          else
             read(pair(:),*) bond_pairs(i,1)
             bond_pairs(i,2) = bond_pairs(i,1)
          end if
          read(buffer(scan(buffer,":")+1:),*) pair_radii(i)
          l_pos = scan(energies(l_pos+1:),",") + l_pos
       end do
    end if


    do i = 1, 3
       cutoff_min_list(i) = read_value_from_string(cutoff_min(i))
       cutoff_max_list(i) = read_value_from_string(cutoff_max(i))
       write(*,*) "Cutoff: ",cutoff_min_list(i),cutoff_max_list(i)
    end do

    width_list = width
    sigma_list = sigma

    !---------------------------------------------------------------------------
    ! close input file
    !---------------------------------------------------------------------------
    close(unit)
    write(*,*) "Input file read successfully."

    return
  end subroutine read_input_file
!###############################################################################


!###############################################################################
  function read_value_from_string(string) result(output)
    !! Read a formatted numeric from a string.
    implicit none

    ! Arguments
    character(*), intent(in) :: string
    !! Input string.
    real(real32) :: output
    !! Output value.

    ! Local variables
    integer :: k, pos
    !! Loop index, position.
    real(real32) :: variable, power
    !! Variable and power values.
    character(:), allocatable :: string_
    !! Copy of input string.
    character(12) :: numeric_set = "0123456789.-"
    !! Numeric set.

    pos = 1
    output = 1._real32
    variable = 0._real32
    power = 1._real32
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
       !! if *, then change power factor to 1._real32
       !! if /, then change power factor to -1._real32
       !! if pi, then change power factor to 1._real32 and variable = pi
       !! if blank space, move pos to next non-space character and cycle
       !! if end of string, exit loop

       if (string_(pos:pos).eq."*")then
          power = 1._real32
          pos = pos + 1
       elseif(string_(pos:pos).eq."/")then
          power = -1._real32
          pos = pos + 1
       elseif(string_(pos:pos+1).eq."pi")then
          power = 1._real32
          output = output * pi ** power
          pos = pos + 2
       end if
       if(pos.gt.len_trim(string_)) exit loop
       pos = pos + verify(string_(pos:), " ") - 1
       if(pos.gt.len_trim(string_)) exit loop

    end do loop

  end function read_value_from_string
!###############################################################################

end module inputs
