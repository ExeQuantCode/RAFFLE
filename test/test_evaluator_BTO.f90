program test_evaluator_BTO
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
  character(3), dimension(3) :: element_symbols
  real(real32), dimension(3) :: element_energies
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
  viability_printing_file = 'viability_BTO.dat'
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
                   &given. Default filename = 'viability_BTO.dat'."
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
  database(1)%nspec = 3
  database(1)%natom = 5
  allocate(database(1)%spec(database(1)%nspec))
  database(1)%spec(1)%num = 1
  database(1)%spec(1)%name = 'Ba'
  allocate(database(1)%spec(1)%atom(database(1)%spec(1)%num, 3))
  database(1)%spec(1)%atom(1, :3) = [0.0, 0.0, 0.0]
  database(1)%spec(2)%num = 1
  database(1)%spec(2)%name = 'Ti'
  allocate(database(1)%spec(2)%atom(database(1)%spec(2)%num, 3))
  database(1)%spec(2)%atom(1, :3) = [0.5, 0.5, 0.5]
  database(1)%spec(3)%num = 3
  database(1)%spec(3)%name = 'O'
  allocate(database(1)%spec(3)%atom(database(1)%spec(3)%num, 3))
  database(1)%spec(3)%atom(1, :3) = [0.5, 0.5, 0.0]
  database(1)%spec(3)%atom(2, :3) = [0.5, 0.0, 0.5]
  database(1)%spec(3)%atom(3, :3) = [0.0, 0.5, 0.5]

  database(1)%lat(1,:) = [4.01, 0.0, 0.0]
  database(1)%lat(2,:) = [0.0, 4.01, 0.0]
  database(1)%lat(3,:) = [0.0, 0.0, 4.00]
  database(1)%energy = -41.0


  !-----------------------------------------------------------------------------
  ! set up element energies
  !-----------------------------------------------------------------------------
  element_symbols(1) = 'Ba'
  element_energies(1) = -1.8816533088684082
  element_symbols(2) = 'Ti'
  element_energies(2) = -7.784755706787109
  element_symbols(3) = 'O'
  element_energies(3) = -4.818619251251221
  call generator%distributions%set_element_energies( &
       element_symbols, &
       element_energies &
  )
  !   call generator%distributions%set_bond_radii( &
  !       reshape(['Ba ','Ti '],shape=[1,2]), &
  !       [1.5] &
  !   )
  !   call generator%distributions%set_bond_radii( &
  !       reshape(['Ti ','O  '],shape=[1,2]), &
  !       [0.7] &
  !   )


  !-----------------------------------------------------------------------------
  ! set up host structure
  !-----------------------------------------------------------------------------
  basis_host%sysname = 'BaTiO3'
  basis_host%nspec = 3
  basis_host%natom = 10
  allocate(basis_host%spec(basis_host%nspec))
  basis_host%spec(1)%num = 2
  basis_host%spec(1)%name = 'Ba'
  allocate(basis_host%spec(1)%atom(basis_host%spec(1)%num, 3))
  basis_host%spec(1)%atom(1, :3) = [0.0, 0.0, 0.0]
  basis_host%spec(1)%atom(2, :3) = [0.0, 0.0, 0.5]
  basis_host%spec(2)%num = 2
  basis_host%spec(2)%name = 'Ti'
  allocate(basis_host%spec(2)%atom(basis_host%spec(2)%num, 3))
  basis_host%spec(2)%atom(1, :3) = [0.5, 0.5, 0.25]
  basis_host%spec(2)%atom(2, :3) = [0.5, 0.5, 0.75]
  basis_host%spec(3)%num = 6
  basis_host%spec(3)%name = 'O'
  allocate(basis_host%spec(3)%atom(basis_host%spec(3)%num, 3))
  basis_host%spec(3)%atom(1, :3) = [0.5, 0.5, 0.0]
  basis_host%spec(3)%atom(2, :3) = [0.5, 0.0, 0.25]
  basis_host%spec(3)%atom(3, :3) = [0.0, 0.5, 0.25]
  basis_host%spec(3)%atom(4, :3) = [0.5, 0.5, 0.5]
  basis_host%spec(3)%atom(5, :3) = [0.5, 0.0, 0.75]
  basis_host%spec(3)%atom(6, :3) = [0.0, 0.5, 0.75]
  basis_host%lat(1,:) = [4.01, 0.0, 0.0]
  basis_host%lat(2,:) = [0.0, 4.01, 0.0]
  basis_host%lat(3,:) = [0.0, 0.0, 8.00]

  ! set up atom_ignore_list
  allocate(atom_ignore_list(2,5))
  atom_ignore_list(:,1) = [1, 2]
  atom_ignore_list(:,2) = [3, 4]
  atom_ignore_list(:,3) = [3, 5]
  atom_ignore_list(:,4) = [3, 6]
  atom_ignore_list(:,5) = [2, 2]

  call basis_host%create_images( &
       max_bondlength = max_bondlength, &
       atom_ignore_list = atom_ignore_list &
  )



  generator%distributions%kbt = 0.2
  call generator%set_host(basis_host)
  call generator%set_grid( grid_spacing = 0.1, grid_offset = [0.0, 0.0, 0.0] )
  generator%distributions%radius_distance_tol = [1.5, 2.5, 3.0, 6.0]
  call generator%distributions%set_width([0.025, pi/200.0, pi/200.0])


  !-----------------------------------------------------------------------------
  ! set up distribution functions
  !-----------------------------------------------------------------------------
  call generator%distributions%create( &
       basis_list = database, &
       deallocate_systems = .true. &
  )
  call generator%distributions%set_element_map( &
       [ "Ba ", "Ti ", "O  "] &
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
       [ 1, 2, 3 ], &
       [ generator%distributions%bond_info(:)%radius_covalent ], &
       atom_ignore_list, &
       grid_offset = generator%grid_offset &
  )
  do i = 1, 3
     tolerance(i) = 1._real32 / real(generator%grid(i),real32) !/ 2._real32
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
     do is = 1, basis_host%nspec
        do i = 1, size(gridpoints,dim=2)
           viability_grid(is,i) = evaluate_point( generator%distributions, &
                gridpoints(1:3,i), is, &
                basis_host, &
                atom_ignore_list(:,ia:), &
                [ generator%distributions%bond_info(:)%radius_covalent ] &
           )
        end do
     end do
     best_loc = maxloc(viability_grid(atom_ignore_list(1,ia),:),dim=1)
     ltmp1 = .false.
     if(any(viability_grid(atom_ignore_list(1,ia),:) .gt. 0._real32) )then
        do ja = ia, size(atom_ignore_list,2), 1
           if( atom_ignore_list(1,ja) .ne. atom_ignore_list(1,ia) ) cycle
           if( &
                all( &
                     abs( &
                          gridpoints(1:3,best_loc) - &
                          basis_host%spec(atom_ignore_list(1,ja))%atom( &
                               atom_ignore_list(2,ja),:3 &
                          ) &
                     ) .lt. tolerance + 1.E-6_real32 &
                ) &
           ) ltmp1 = .true.
        end do
     end if
     call assert( &
          ltmp1, &
          "Incorrect gridpoint found.", &
          success &
     )
     if( .not. ltmp1 ) then
        write(*,*) "species: ", atom_ignore_list(1,ia)
        write(*,*) "Best location: ", best_loc
        write(*,*) "viability: ", &
             viability_grid(atom_ignore_list(1,ia),best_loc)
        write(*,*) "Gridpoint: ", gridpoints(1:3,best_loc)
     end if
     call basis_host%update_images( &
          max_bondlength = max_bondlength, &
          is = atom_ignore_list(1,ia), ia = atom_ignore_list(2,ia) &
     )
  end do
  close(unit)


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

end program test_evaluator_BTO
