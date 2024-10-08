program test_evaluator
  use error_handling
  use constants, only: real12
  use rw_geom, only: basis_type, geom_write
  use extended_geom, only: extended_basis_type
  use evaluator, only: evaluate_point
  use generator, only: raffle_generator_type
  use add_atom, only: get_gridpoints_and_viability
  implicit none


  integer :: unit
  integer :: i, ia, num_points
  integer :: best_loc
  real(real12) :: max_bondlength
  type(extended_basis_type) :: basis_host
  type(basis_type), dimension(1) :: database
  character(3), dimension(1) :: element_symbols
  real(real12), dimension(1) :: element_energies
  real(real12), dimension(3) :: tolerance
  integer, dimension(:,:), allocatable :: atom_ignore_list

  real(real12), dimension(:), allocatable :: suitability_grid
  real(real12), dimension(:,:), allocatable :: gridpoints

  type(raffle_generator_type) :: generator

  logical :: success = .true.

  test_error_handling = .true.


  max_bondlength = 6._real12
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
  allocate(atom_ignore_list(8,2))
  do i = 1, size(atom_ignore_list,1)
     atom_ignore_list(i,1) = 1
     atom_ignore_list(i,2) = &
          basis_host%spec(1)%num - size(atom_ignore_list,1) + i
  end do

  call basis_host%create_images( &
       max_bondlength = max_bondlength, &
       atom_ignore_list = atom_ignore_list &
  )


  generator%distributions%kBT = 0.2
  call generator%host%copy(basis_host)
  call generator%set_grid( grid_spacing = 0.2, grid_offset = [0.0, 0.0, 0.0] )
  generator%distributions%radius_distance_tol = [1.5, 2.5, 3.0, 6.0]


  !-----------------------------------------------------------------------------
  ! set up distribution functions
  !-----------------------------------------------------------------------------
  call generator%distributions%create( &
       basis_list = database, &
       deallocate_systems = .true. &
  )


  !-----------------------------------------------------------------------------
  ! set up gridpoints
  !-----------------------------------------------------------------------------
  num_points = 0
  gridpoints = get_gridpoints_and_viability( &
       generator%distributions, &
       generator%grid, &
       basis_host, &
       [ 1 ], &
       [ generator%distributions%bond_info(:)%radius_covalent ], &
       atom_ignore_list, &
       grid_offset = generator%grid_offset &
  )
  do i = 1, 3
     tolerance(i) = 1._real12 / real(generator%grid(i),real12) / 2._real12
  end do


  !-----------------------------------------------------------------------------
  ! call evaluator
  !-----------------------------------------------------------------------------
  allocate(suitability_grid(size(gridpoints,2)))
  do ia = 1, size(atom_ignore_list,1)
     suitability_grid(:) = 0._real12
     do i = 1, size(gridpoints,dim=2)
        suitability_grid(i) = evaluate_point( generator%distributions, &
             gridpoints(1:3,i), atom_ignore_list(ia,1), basis_host, &
             atom_ignore_list(ia:,:), &
             [ generator%distributions%bond_info(:)%radius_covalent ] &
        )
     end do
     best_loc = maxloc(suitability_grid,dim=1)
     ! Check point is correct
     call assert( &
          all( &
               abs( &
                    gridpoints(1:3,best_loc) - &
                    basis_host%spec(1)%atom(atom_ignore_list(ia,2),:3) &
               ) .lt. tolerance + 1.E-6_real12 &
          ), &
          "Incorrect gridpoint found.", &
          success &
     )
     call basis_host%update_images( &
          max_bondlength = max_bondlength, &
          is = 1, ia = atom_ignore_list(ia,2) &
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

end program test_evaluator