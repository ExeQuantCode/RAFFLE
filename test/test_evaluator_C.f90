program test_evaluator_C
  use raffle__io_utils
  use raffle__constants, only: real32, pi
  use raffle__misc_linalg, only: modu
  use raffle__geom_rw, only: basis_type, geom_write
  use raffle__geom_extd, only: extended_basis_type
  use raffle__evaluator, only: evaluate_point
  use raffle__generator, only: raffle_generator_type
  use raffle__viability, only: get_gridpoints_and_viability
  implicit none


  integer :: unit
  integer :: i, is, ia, ja, num_points
  integer :: best_loc
  real(real32) :: max_bondlength
  type(extended_basis_type) :: basis_host
  logical :: ltmp1
  type(basis_type), dimension(1) :: database
  character(3), dimension(1) :: element_symbols
  real(real32), dimension(1) :: element_energies
  real(real32), dimension(3) :: tolerance
  integer, dimension(:,:), allocatable :: atom_ignore_list

  integer :: iostat
  logical :: viability_printing
  character(len=256) :: arg, arg_prev, viability_printing_file, fmt

  real(real32), dimension(:,:), allocatable :: gridpoints, viability_grid

  type(raffle_generator_type) :: generator

  logical :: success = .true.

  test_error_handling = .true.


  !-----------------------------------------------------------------------------
  ! check for input argument flags
  !-----------------------------------------------------------------------------
  viability_printing = .false.
  viability_printing_file = 'viability_C.dat'
  if( command_argument_count() .ne. 0 ) then
     i = 1
     do
        call get_command_argument(i, arg, status=iostat)
        if( iostat .ne. 0 ) exit
        if(index(arg,'-').eq.1) then
           if( arg == '-h' .or. arg == '--help' ) then
              ! print description of unit test and associated flags
              write(*,*) "This unit test evaluates the evaluator module using &
                   &the BaTiO3 structure."
              write(*,*) "Flags:"
              write(*,*) "-h, --help: Print this help message"
              write(*,*) "-p, --print [filename]: Print the gridpoints and &
                   &their viability values to a file. If no filename is &
                   &given. Default filename = 'viability_C.dat'."
              stop 0
           elseif( index(arg,'-p').eq.1 .or. index(arg,'--print').eq.1 ) then
              viability_printing = .true.
              if( index(arg,'-p').eq.1 .and. trim(arg).ne.'-p' )then
                 viability_printing_file = trim(adjustl(arg(3:)))
              elseif( index(arg,'--print').eq.1 .and. &
                   trim(arg).ne.'--print' )then
                 viability_printing_file = trim(adjustl(arg(8:)))
              end if
           else
              write(0,*) "Unknown flag: ", arg
              stop 1
           end if
        else
           call get_command_argument(i-1, arg_prev, status=iostat)
           if( index(arg,'-p').eq.1 .or. index(arg,'--print').eq.1 ) then
              viability_printing_file = trim(adjustl(arg))
           else
              write(0,*) "Unknown argument: ", arg
              stop 1
           end if
        end if
        i = i + 1
     end do
  end if


  max_bondlength = 6._real32
  !-----------------------------------------------------------------------------
  ! set up database
  !-----------------------------------------------------------------------------
  database(1)%nspec = 1
  database(1)%natom = 8
  allocate(database(1)%spec(database(1)%nspec))
  database(1)%spec(1)%num = 8
  database(1)%spec(1)%name = 'C'
  allocate(database(1)%spec(1)%atom(database(1)%spec(1)%num, 3))
  database(1)%spec(1)%atom(1, :3) = [0.0, 0.0, 0.0]
  database(1)%spec(1)%atom(2, :3) = [0.5, 0.5, 0.0]
  database(1)%spec(1)%atom(3, :3) = [0.5, 0.0, 0.5]
  database(1)%spec(1)%atom(4, :3) = [0.0, 0.5, 0.5]
  database(1)%spec(1)%atom(5, :3) = [0.25, 0.25, 0.25]
  database(1)%spec(1)%atom(6, :3) = [0.75, 0.75, 0.25]
  database(1)%spec(1)%atom(7, :3) = [0.75, 0.25, 0.75]
  database(1)%spec(1)%atom(8, :3) = [0.25, 0.75, 0.75]

  database(1)%lat(1,:) = [3.5607451090903233, 0.0, 0.0]
  database(1)%lat(2,:) = [0.0, 3.5607451090903233, 0.0]
  database(1)%lat(3,:) = [0.0, 0.0, 3.5607451090903233]
  database(1)%energy = -72.213492


  !-----------------------------------------------------------------------------
  ! set up element energies
  !-----------------------------------------------------------------------------
  element_symbols(1) = 'C'
  element_energies(1) = -9.0266865
  call generator%distributions%set_element_energies( &
       element_symbols, &
       element_energies &
  )


  !-----------------------------------------------------------------------------
  ! set up host structure
  !-----------------------------------------------------------------------------
  basis_host%sysname = 'diamond'
  basis_host%nspec = 1
  allocate(basis_host%spec(basis_host%nspec))
  basis_host%spec(1)%num = 16
  basis_host%spec(1)%name = 'C'
  basis_host%natom = sum(basis_host%spec(:)%num)
  allocate(basis_host%spec(1)%atom(basis_host%spec(1)%num, 3))
  basis_host%spec(1)%atom(1, :3) = [0.0, 0.0, 0.0]
  basis_host%spec(1)%atom(2, :3) = [0.5, 0.5, 0.0]
  basis_host%spec(1)%atom(3, :3) = [0.5, 0.0, 0.25]
  basis_host%spec(1)%atom(4, :3) = [0.0, 0.5, 0.25]
  basis_host%spec(1)%atom(5, :3) = [0.25, 0.25, 0.125]
  basis_host%spec(1)%atom(6, :3) = [0.75, 0.75, 0.125]
  basis_host%spec(1)%atom(7, :3) = [0.75, 0.25, 0.375]
  basis_host%spec(1)%atom(8, :3) = [0.25, 0.75, 0.375]
  basis_host%spec(1)%atom(9, :3) = [0.0, 0.0, 0.5]
  basis_host%spec(1)%atom(10, :3) = [0.5, 0.5, 0.5]
  basis_host%spec(1)%atom(11, :3) = [0.75, 0.25, 0.875]
  basis_host%spec(1)%atom(12, :3) = [0.25, 0.75, 0.875]
  basis_host%spec(1)%atom(13, :3) = [0.0, 0.5, 0.75]
  basis_host%spec(1)%atom(14, :3) = [0.5, 0.0, 0.75]
  basis_host%spec(1)%atom(15, :3) = [0.75, 0.75, 0.625]
  basis_host%spec(1)%atom(16, :3) = [0.25, 0.25, 0.625]
  basis_host%lat(1,:) = [3.560745109, 0.0, 0.0]
  basis_host%lat(2,:) = [0.0, 3.560745109, 0.0]
  basis_host%lat(3,:) = [0.0, 0.0, 7.121490218]

  ! set up atom_ignore_list
  allocate(atom_ignore_list(2,8))
  do i = 1, size(atom_ignore_list,2)
     atom_ignore_list(1,i) = 1
     atom_ignore_list(2,i) = &
          basis_host%spec(1)%num - size(atom_ignore_list,2) + i
  end do

  call basis_host%create_images( &
       max_bondlength = max_bondlength, &
       atom_ignore_list = atom_ignore_list &
  )


  generator%distributions%kBT = 0.2
  call generator%set_host(basis_host)
  call generator%set_grid( grid_spacing = 0.2, grid_offset = [0.0, 0.0, 0.0] )
  generator%distributions%radius_distance_tol = [1.5, 2.5, 3.0, 5.0]
  call generator%distributions%set_width([0.025, pi/200.0, pi/200.0])


  !-----------------------------------------------------------------------------
  ! set up distribution functions
  !-----------------------------------------------------------------------------
  call generator%distributions%create( &
       basis_list = database, &
       deallocate_systems = .true. &
  )
  call generator%distributions%set_element_map( &
       [ "C  "] &
  )
  call generator%distributions%host_system%set_element_map( &
       generator%distributions%element_info &
  )


  !-----------------------------------------------------------------------------
  ! set up gridpoints
  !-----------------------------------------------------------------------------
  num_points = 0
  gridpoints = get_gridpoints_and_viability( &
       generator%distributions, &
       generator%grid, &
       generator%bounds, &
       basis_host, &
       [ 1 ], &
       [ generator%distributions%bond_info(:)%radius_covalent ], &
       atom_ignore_list, &
       grid_offset = generator%grid_offset &
  )
  do i = 1, 3
     tolerance(i) = 1._real32 / real(generator%grid(i),real32) / 2._real32
  end do


  !-----------------------------------------------------------------------------
  ! print viability data to file
  !-----------------------------------------------------------------------------
  if(viability_printing)then
     write(*,*) "Printing viability data to file: ", &
          trim(viability_printing_file)
     open(newunit=unit, file=viability_printing_file)
     write(unit,'("#grid",3(1X,I0),3(1X,F0.3))') &
          generator%grid, generator%grid_offset
     write(unit,'("#lat",3(1X,F0.3))') &
          modu(basis_host%lat(1,:)), &
          modu(basis_host%lat(2,:)), &
          modu(basis_host%lat(3,:))
     write(fmt,'("(""#species"",",I0,"(1X,A3))")') basis_host%nspec
     write(unit,fmt) basis_host%spec(:)%name
     do is = 1, basis_host%nspec
        atom_loop: do ia = 1, basis_host%spec(is)%num
           do i = 1, size(atom_ignore_list,2)
              if( all(atom_ignore_list(:,i).eq.[is,ia]) ) cycle atom_loop
           end do
           write(unit,*) basis_host%spec(is)%atom(ia,:3)
        end do atom_loop
     end do
     write(unit,*)
     do is = 1, basis_host%nspec
        do i = 1, size(gridpoints,dim=2)
           write(unit,*) gridpoints(1:3,i), gridpoints(4+is,i), is
        end do
     end do
     close(unit)
     stop 0
  end if


  !-----------------------------------------------------------------------------
  ! call evaluator
  !-----------------------------------------------------------------------------
  allocate(viability_grid(basis_host%nspec,size(gridpoints,2)))
  do ia = 1, size(atom_ignore_list,2)
     viability_grid(:,:) = 0._real32
     do i = 1, size(gridpoints,dim=2)
        viability_grid(1,i) = evaluate_point( generator%distributions, &
             gridpoints(1:3,i), atom_ignore_list(1,ia), basis_host, &
             atom_ignore_list(:,ia:), &
             [ generator%distributions%bond_info(:)%radius_covalent ] &
        )
     end do
     best_loc = maxloc(viability_grid(atom_ignore_list(1,ia),:),dim=1)
     ! Check point is correct
     ltmp1 = .false.
     do ja = ia, size(atom_ignore_list,2), 1
        if( &
             all( &
                  abs( &
                       gridpoints(1:3,best_loc) - &
                       basis_host%spec(1)%atom(atom_ignore_list(2,ja),:3) &
                  ) .lt. tolerance + 1.E-6_real32 &
             ) &
        ) ltmp1 = .true.
     end do
     call assert( &
          ltmp1, &
          "Incorrect gridpoint found.", &
          success &
     )
     call basis_host%update_images( &
          max_bondlength = max_bondlength, &
          is = 1, ia = atom_ignore_list(2,ia) &
     )
  end do


  !-----------------------------------------------------------------------------
  ! check for any failed tests
  !-----------------------------------------------------------------------------
  write(*,*) "----------------------------------------"
  if(success)then
     write(*,*) 'test_evaluator passed all tests'
  else
     write(0,*) 'test_evaluator failed one or more tests'
     stop 1
  end if

contains

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

end program test_evaluator_C
