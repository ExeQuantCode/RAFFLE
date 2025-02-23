program raffle_program
  !! Main program for the interface-based random structure search.
  !!
  !! This program generates random structures based on the host structure and
  !! the distribution functions of the elements and bonds. The structures are
  !! generated using a random structure generator, which is based on the
  !! the RAFFLE method.
  use raffle__constants, only: real32
  use raffle__io_utils, only: stop_program
  use raffle__misc, only: touch
  use inputs
  use read_structures, only: get_gdfs_from_data
  use raffle, only: raffle_generator_type, distribs_container_type
  use raffle__geom_rw, only: geom_read, geom_write, basis_type
  implicit none

  ! Local variables
  integer :: i, unit, itmp1
  !! Loop index, unit number, temporary integer
  character(1024) :: buffer
  !! Buffer for strings
  character(:), allocatable :: next_dir
  !! Next directory name
  type(basis_type) :: host

  real(real32), dimension(:), allocatable :: tmp_energies
  !! Temporary array for element energies
  character(len=3), dimension(:), allocatable :: tmp_symbols
  !! Temporary array for element symbols

  type(raffle_generator_type) :: generator
  !! Random structure generator


  !-----------------------------------------------------------------------------
  ! read input file
  !-----------------------------------------------------------------------------
  call set_global_vars()


  !-----------------------------------------------------------------------------
  ! check the task and run the appropriate case
  !-----------------------------------------------------------------------------
  ! 0) Run structure search
  ! 1) Continue structure search
  select case(task)
  case(0)
     allocate(character(len=len_trim(output_dir)+1) :: next_dir)
     write(next_dir, '(A,"0")') trim(output_dir)
     write(*,'("Running host-based random structure search for ",I0,&
          &" structures")' &
     ) num_structures
  case(1)
     call get_next_directory(trim(output_dir), itmp1, next_dir)
     write(*,'("Running iteration ",I0,&
          &" of host-based random structure search for ",I0," structures")' &
     ) itmp1 + 1, num_structures
  case default
     call stop_program("Invalid option")
  end select


  !-----------------------------------------------------------------------------
  ! set the element energies and bond radii, if they are provided
  !-----------------------------------------------------------------------------
  if(allocated(element_symbols).and.allocated(element_energies)) then
     call generator%distributions%set_element_energies( &
          element_symbols, &
          element_energies &
     )
  end if

  if(allocated(bond_pairs).and.allocated(pair_radii)) then
     call generator%distributions%set_bond_radii( &
          bond_pairs, &
          pair_radii &
     )
  end if


  !-----------------------------------------------------------------------------
  ! read structures from the database and generate gvectors
  !-----------------------------------------------------------------------------
  generator%distributions = get_gdfs_from_data( &
       input_dir    = database_list, &
       file_format  = database_format, &
       distribs_container_template = distribs_container_type(&
            width = width_list, &
            sigma = sigma_list, &
            cutoff_min = cutoff_min_list, &
            cutoff_max = cutoff_max_list ) &
  )

  call generator%distributions%write_2body(file="2body.txt")
  call generator%distributions%write_3body(file="3body.txt")
  call generator%distributions%write_4body(file="4body.txt")

  call generator%distributions%get_element_energies( &
       tmp_symbols, &
       tmp_energies &
  )
  do i = 1, size(tmp_symbols)
     write(*,*) "Element ", tmp_symbols(i), " energy: ", tmp_energies(i)
  end do


  !-----------------------------------------------------------------------------
  ! set the host structure
  !-----------------------------------------------------------------------------
  open(newunit=unit, file=filename_host, status='old')
  call geom_read(unit, host)
  call generator%set_host(host)
  close(unit)
  if(grid_spacing.gt.1.E-6.and.all(grid.ne.0))then
     call stop_program('Cannot specify grid spacing and grid at the same time')
  elseif(grid_spacing.gt.1.E-6)then
     call generator%set_grid(grid_spacing = grid_spacing)
  elseif(all(grid.ne.0))then
     call generator%set_grid(grid = grid)
  end if


  !-----------------------------------------------------------------------------
  ! generate random structures
  !-----------------------------------------------------------------------------
  write(*,*) "Generating structures"
  call generator%generate( num_structures, &
       stoich, &
       method_ratio )
  write(*,*) "Structures have been successfully generated"


  !-----------------------------------------------------------------------------
  ! save generated structures
  !-----------------------------------------------------------------------------
  do i = 1, generator%num_structures
     write(buffer,'(A,"/struc",I0.3)') trim(next_dir),i
     call touch(buffer)
     open(newunit = unit, file=trim(buffer)//"/POSCAR")
     call geom_write(unit, generator%structures(i))
     close(unit)
  end do

contains

!###############################################################################
  subroutine get_next_directory(prefix, num, next_dir)
    !! Get the next directory name
    implicit none

    ! Arguments
    character(len=*), intent(in) :: prefix
    !! Prefix of the directory name
    integer, intent(out) :: num
    !! Number of directories with the prefix
    character(len=:), allocatable, intent(out) :: next_dir
    !! Next directory name

    ! Local variables
    integer :: i, ierror, itmp1, max_num
    !! Loop index, error code, temporary integer
    integer :: unit
    !! Unit number for the pipe
    character(len=1024) :: pattern, line, ctmp1
    !! Pattern to match directories
    character(len=1024), dimension(:), allocatable :: dir_list
    !! List of directories

    ! Initialise variables
    num = 0
    max_num = 0
    pattern = trim(prefix) // "[0-9]*"

    ! Command to list directories matching the pattern
    call execute_command_line('ls -d ' // trim(pattern) // ' 2>/dev/null')

    ! Open a pipe to execute the command
    open( newunit=unit, file='iter_count.txt', &
         status='old', action='read', iostat=ierror)
    if (ierror .ne. 0) then
       call stop_program('Opening pipe to list directories')
       return
    end if

    ! Read the output of the command and get number of directories
    allocate(dir_list(0))  ! Allocate a large enough array
    num = 0
    do
       read(unit, '(A)', iostat=ierror) line
       if (ierror .ne. 0) exit
       dir_list = [ dir_list, line ]
       num = num + 1
    end do
    close(unit)

    ! Count the number of matching directories and find the highest number
    do i = 1, num
       ctmp1 = dir_list(i)
       ctmp1 = ctmp1(index(ctmp1,prefix)+len_trim(prefix):)
       read(ctmp1, *, iostat=ierror) itmp1
       if (ierror .eq. 0) then
          if (itmp1 .gt. max_num) max_num = itmp1
       end if
    end do

    ! Determine the next directory name
    allocate(character(len=len_trim(prefix)+ceiling(max_num/10.0)) :: next_dir)
    write(next_dir, '(A,I0)') prefix, max_num + 1

    ! Deallocate the array
    deallocate(dir_list)
  end subroutine get_next_directory
!###############################################################################

end program raffle_program