program test_generator
  use error_handling
  use constants, only: real12
  use rw_geom, only: basis_type
  use generator, only: raffle_generator_type, stoichiometry_type
  implicit none

  integer :: i
  type(raffle_generator_type) :: generator
  class(raffle_generator_type), allocatable :: generator_var
  type(basis_type) :: basis_host, basis_host_expected
  type(basis_type), dimension(1) :: database
  character(3), dimension(1) :: element_symbols
  real(real12), dimension(1) :: element_energies
  real(real12), dimension(3) :: tolerance

  logical :: success = .true.

  test_error_handling = .true.


  !-----------------------------------------------------------------------------
  ! test generator setup
  !-----------------------------------------------------------------------------
  generator_var = raffle_generator_type( &
       width = [ 1.0, 1.0, 1.0 ], &
       sigma = [ 0.1, 0.1, 0.1 ], &
       cutoff_min = [ 0.0, 0.0, 0.0 ], &
       cutoff_max = [ 1.0, 1.0, 1.0 ] &
  )
  call assert( &
       all( abs( generator_var%distributions%width - 1.0 ) .lt. 1.E-6 ), &
       'Generator failed to set width', &
       success &
  )
  call assert( &
       all( abs( generator_var%distributions%sigma - 0.1 ) .lt. 1.E-6 ), &
       'Generator failed to set sigma', &
       success &
  )
  call assert( &
       all( abs( generator_var%distributions%cutoff_min - 0.0 ) .lt. 1.E-6 ), &
       'Generator failed to set cutoff_min', &
       success &
  )
  call assert( &
       all( abs( generator_var%distributions%cutoff_max - 1.0 ) .lt. 1.E-6 ), &
       'Generator failed to set cutoff_max', &
       success &
  )
  write(*,*) "Testing generator missing host handling"
  call generator_var%generate( num_structures = 1, &
       stoichiometry = [ stoichiometry_type(element='C  ', num = 8) ], &
       method_probab = [0.0, 0.0, 1.0], &
       seed = 0, &
       verbose = 1 &
  )
  call assert( &
       generator_var%num_structures .eq. 0, &
       'Generator failed to handle missing host structure', &
       success &
  )
  write(*,*) "Handled missing host"
  deallocate(generator_var)


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
  basis_host%spec(1)%num = 8
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
  basis_host%lat(1,:) = [3.560745109, 0.0, 0.0]
  basis_host%lat(2,:) = [0.0, 3.560745109, 0.0]
  basis_host%lat(3,:) = [0.0, 0.0, 7.121490218]


  !-----------------------------------------------------------------------------
  ! test generator setup with host
  !-----------------------------------------------------------------------------
  generator_var = raffle_generator_type(host = basis_host)
  call assert(compare_bas(generator_var%host, basis_host), &
       'Generator failed to set host structure', &
       success &
  )
  deallocate(generator_var)


  !-----------------------------------------------------------------------------
  ! set up expected generated host structure
  !-----------------------------------------------------------------------------
  basis_host_expected%sysname = 'diamond+inserts'
  basis_host_expected%nspec = 1
  allocate(basis_host_expected%spec(basis_host_expected%nspec))
  basis_host_expected%spec(1)%num = 16
  basis_host_expected%spec(1)%name = 'C'
  basis_host_expected%natom = sum(basis_host_expected%spec(:)%num)
  allocate(basis_host_expected%spec(1)%atom(basis_host_expected%spec(1)%num, 3))
  basis_host_expected%spec(1)%atom(1, :3) = [0.0, 0.0, 0.0]
  basis_host_expected%spec(1)%atom(2, :3) = [0.5, 0.5, 0.0]
  basis_host_expected%spec(1)%atom(3, :3) = [0.5, 0.0, 0.25]
  basis_host_expected%spec(1)%atom(4, :3) = [0.0, 0.5, 0.25]
  basis_host_expected%spec(1)%atom(5, :3) = [0.25, 0.25, 0.125]
  basis_host_expected%spec(1)%atom(6, :3) = [0.75, 0.75, 0.125]
  basis_host_expected%spec(1)%atom(7, :3) = [0.75, 0.25, 0.375]
  basis_host_expected%spec(1)%atom(8, :3) = [0.25, 0.75, 0.375]
  basis_host_expected%spec(1)%atom(9, :3) = [0.75, 0.25, 0.875]
  basis_host_expected%spec(1)%atom(10, :3) = [0.75, 0.75, 0.625]
  basis_host_expected%spec(1)%atom(11, :3) = [0.5, 0.0, 0.75]
  basis_host_expected%spec(1)%atom(12, :3) = [0.25, 0.25, 0.625]
  basis_host_expected%spec(1)%atom(13, :3) = [0.25, 0.75, 0.875]
  basis_host_expected%spec(1)%atom(14, :3) = [0.5, 0.5, 0.5]
  basis_host_expected%spec(1)%atom(15, :3) = [0.0, 0.5, 0.75]
  basis_host_expected%spec(1)%atom(16, :3) = [0.0, 0.0, 0.5]
  basis_host_expected%lat(1,:) = [3.560745109, 0.0, 0.0]
  basis_host_expected%lat(2,:) = [0.0, 3.560745109, 0.0]
  basis_host_expected%lat(3,:) = [0.0, 0.0, 7.121490218]


  !-----------------------------------------------------------------------------
  ! set up host
  !-----------------------------------------------------------------------------
  call generator%set_host( basis_host )
  call assert(compare_bas(generator%host, basis_host), &
       'Generator failed to set host structure', &
       success &
  )

  !-----------------------------------------------------------------------------
  ! test grid setup
  !-----------------------------------------------------------------------------
  write(*,*) "Testing grid error handling"
  call generator%set_grid( grid = [ 1, 1, 1 ], grid_spacing = 0.2 )
  call assert( &
       all( generator%grid .eq. 0 ), &
       'Generator failed to handle invalid grid', &
       success &
  )
  write(*,*) "Handled grid setup error"
  call generator%set_grid( grid = [ 1, 1, 1 ] )
  call assert( &
       all( generator%grid .eq. 1 ), &
       'Generator failed to handle grid', &
       success &
  )

  ! check grid resetting
  call generator%reset_grid()
  call assert( &
       all( generator%grid .eq. 0 ), &
       'Generator failed to reset grid', &
       success &
  )

  ! check grid spacing
  call generator%set_grid( grid_spacing = 0.2 )
  call assert( &
       all( generator%grid .eq. [ 18, 18, 36 ] ), &
       'Generator failed to handle grid_spacing', &
       success &
  )

  ! check grid_offset
  call generator%set_grid( grid_spacing = 0.2, grid_offset = [0.1, 0.2, 0.3] )
  call assert( &
       all( &
            abs( &
                 generator%grid_offset - &
                 [0.1, 0.2, 0.3] &
            ) .lt. 1.E-6_real12 &
       ), &
       'Generator failed to handle grid_offset', &
       success &
  )

  call generator%reset_grid()


  !-----------------------------------------------------------------------------
  ! set up generator
  !-----------------------------------------------------------------------------
  generator%distributions%kbt = 0.2
  call generator%set_grid( grid_spacing = 0.2, grid_offset = [0.0, 0.0, 0.0] )
  generator%distributions%radius_distance_tol = [1.5, 2.5, 3.0, 6.0]
  do i = 1, 3
     tolerance(i) = 1._real12 / real(generator%grid(i),real12) / 2._real12
  end do


  !-----------------------------------------------------------------------------
  ! set up distribution functions
  !-----------------------------------------------------------------------------
  call generator%distributions%create( &
       basis_list = database, &
       deallocate_systems = .true. &
  )

  !-----------------------------------------------------------------------------
  ! test generator random seed and no printing
  !-----------------------------------------------------------------------------
  call generator%generate( num_structures = 0, &
       stoichiometry = [ stoichiometry_type(element='C  ', num = 8) ], &
       method_probab = [0.0, 0.0, 1.0] &
  )

  !-----------------------------------------------------------------------------
  ! generate random structures
  !-----------------------------------------------------------------------------
  call generator%generate( num_structures = 1, &
       stoichiometry = [ stoichiometry_type(element='C  ', num = 8) ], &
       method_probab = [0.0, 0.0, 1.0], &
       seed = 0, &
       verbose = 1 &
  )


  !-----------------------------------------------------------------------------
  ! compare the generated host structure
  !-----------------------------------------------------------------------------
  call assert(compare_bas(generator%structures(1), basis_host_expected), &
       'Generated host structure does not match expected host structure', &
       success &
  )


  !-----------------------------------------------------------------------------
  ! check adding a new structure
  !-----------------------------------------------------------------------------
  call generator%generate( num_structures = 2, &
       stoichiometry = [ stoichiometry_type(element='C  ', num = 12) ], &
       method_probab = [0.0, 0.0, 1.0], &
       seed = 0, &
       verbose = 1 &
  )
  call assert( &
       generator%num_structures .eq. 3, &
       'Generator failed to add a new structure', &
       success &
  )
  call assert( &
       size(generator%structures) .eq. 3, &
       'Generator failed to add a new structure', &
       success &
  )


  !-----------------------------------------------------------------------------
  ! check for any failed tests
  !-----------------------------------------------------------------------------
  write(*,*) "----------------------------------------"
  if(success)then
     write(*,*) 'test_generator passed all tests'
  else
     write(0,*) 'test_generator failed one or more tests'
     stop 1
  end if

contains

  function compare_bas(bas1, bas2) result(output)
    type(basis_type), intent(in) :: bas1, bas2
    logical :: output

    integer :: is, ia
    output = .true.

  ! Compare the geometries
    if(any(abs(bas1%lat - bas2%lat).gt.1.E-6)) then
      output = .false.
   end if
   if(bas1%sysname .ne. bas2%sysname) then
      output = .false.
   end if
   if(bas1%natom .ne. bas2%natom) then
      output = .false.
   end if

   do is = 1, bas1%nspec
      do ia = 1, bas1%spec(is)%num
         if( &
              any( &
                   abs( &
                        bas1%spec(is)%atom(ia,:3) - &
                        bas2%spec(is)%atom(ia,:3) &
                   ) .ge. 2._real12 * tolerance + 1.E-6_real12 &
              ) &
         ) then
            write(0,*) 'Generator failed to produce expected atom: ', is, ia
            write(0,*) bas1%spec(is)%atom(ia,:3), bas2%spec(is)%atom(ia,:3)
         end if
      end do
   end do

  end function compare_bas

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

end program test_generator