program test_evaluator
  use constants, only: real12
  use rw_geom, only: basis_type, geom_write
  use extended_geom, only: extended_basis_type
  use evaluator, only: evaluate_point
  use generator, only: raffle_generator_type
  use add_atom, only: get_viable_gridpoints
  implicit none


  integer :: unit
  integer :: i, j, k, num_points
  integer :: best_loc
  type(extended_basis_type) :: basis_host
  type(basis_type), dimension(1) :: database
  character(3), dimension(1) :: element_symbols
  real(real12), dimension(1) :: element_energies
  integer, dimension(1,2) :: atom_ignore_list

  real(real12), dimension(:), allocatable :: tmp_vals
  real(real12), dimension(:), allocatable :: suitability_grid
  real(real12), dimension(:,:), allocatable :: gridpoints

  character(len=3), dimension(:), allocatable :: tmp_symbols
  character(len=3), dimension(:,:), allocatable :: tmp_pairs

  type(raffle_generator_type) :: generator

  
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



!   ! graphite cell
!   database(1)%nspec = 1
!   database(1)%natom = 4
! !   allocate(database(1)%spec(database(1)%nspec))
!   database(1)%spec(1)%num = 4
!   database(1)%spec(1)%name = 'C'
!   deallocate(database(1)%spec(1)%atom)
!   allocate(database(1)%spec(1)%atom(database(1)%spec(1)%num, 3))
!   database(1)%spec(1)%atom(1, :3) = [0.0, 0.0, 0.25]
!   database(1)%spec(1)%atom(2, :3) = [0.0, 0.0, 0.75]
!   database(1)%spec(1)%atom(3, :3) = [1.0/3.0, 2.0/3.0, 0.25]
!   database(1)%spec(1)%atom(4, :3) = [2.0/3.0, 1.0/3.0, 0.75]

!   database(1)%lat(1,:) = [1.2336456308015413, -2.1367369110836267, 0.0]
!   database(1)%lat(2,:) = [1.2336456308015413,  2.1367369110836267, 0.0]
!   database(1)%lat(3,:) = [0.0, 0.0, 7.8030730000000004]
!   database(1)%energy = -36.86795585

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
  basis_host%spec(1)%num = 11
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
!   basis_host%spec(1)%atom(11, :3) = [0.722, 0.722, 0.6388]
!   basis_host%spec(1)%atom(12, :3) = [0.222, 0.222, 0.6111]
!   basis_host%spec(1)%atom(13, :3) = [-0.056, 0.5, 0.75]
!   basis_host%spec(1)%atom(14, :3) = [0.5, -0.056, 0.75]
  basis_host%spec(1)%atom(basis_host%spec(1)%num, :3) = [0.0, 0.0, 0.0]
  atom_ignore_list(1,1) = 1
  atom_ignore_list(1,2) = basis_host%spec(1)%num

  basis_host%lat(1,:) = [3.560745109, 0.0, 0.0]
  basis_host%lat(2,:) = [0.0, 3.560745109, 0.0]
  basis_host%lat(3,:) = [0.0, 0.0, 7.121490218]


!   basis_host%sysname = 'graphite'
!   basis_host%nspec = 1
!   basis_host%natom = 7
!   allocate(basis_host%spec(database(1)%nspec))
!   basis_host%spec(1)%num = 7
!   basis_host%spec(1)%name = 'C'
!   allocate(basis_host%spec(1)%atom(basis_host%spec(1)%num, 3))
!   basis_host%spec(1)%atom(1, :3) = [0.000000000, 0.000000000, 0.125000000]
!   basis_host%spec(1)%atom(2, :3) = [0.000000000, 0.000000000, 0.375000000]
!   basis_host%spec(1)%atom(3, :3) = [0.000000000, 0.000000000, 0.875000000]
!   basis_host%spec(1)%atom(4, :3) = [0.333333333, 0.666666667, 0.125000000]
!   basis_host%spec(1)%atom(5, :3) = [0.666666667, 0.333333333, 0.375000000]
!   basis_host%spec(1)%atom(6, :3) = [0.666666667, 0.333333333, 0.875000000]
! !   basis_host%spec(1)%atom(7, :3) = [0.333333333, 0.666666667, 0.625000000]
!   basis_host%spec(1)%atom(7, :3) = [0.0, 0.0, 0.0]

!   basis_host%lat(1,:) = [1.233645631, -2.136736911,  0.000000000]
!   basis_host%lat(2,:) = [1.233645631,  2.136736911,  0.000000000]
!   basis_host%lat(3,:) = [0.000000000,  0.000000000, 15.606146000]
!     atom_ignore_list(1,1) = 1
!     atom_ignore_list(1,2) = 8
  call basis_host%create_images( &
       max_bondlength = 6._real12, &
       atom_ignore_list = atom_ignore_list &
  )
  open(newunit=unit, file="POSCAR_host_graphite", status="replace")
  call geom_write(unit, basis_host)
  close(unit)


  generator%distributions%kbt = 0.2
  call generator%host%copy(basis_host)
  call generator%set_grid( grid_spacing = 0.2, grid_offset = [0.0, 0.0, 0.0] )
  generator%distributions%radius_distance_tol = [1.5, 2.5, 3.0, 6.0]


  !-----------------------------------------------------------------------------
  ! set up distribution functions
  !-----------------------------------------------------------------------------
  call generator%distributions%create( &
       basis_list = database, &
       deallocate_systems = .false. &
  )


  !-----------------------------------------------------------------------------
  ! set up gridpoints
  !-----------------------------------------------------------------------------
  num_points = 0
!   allocate(gridpoints(3, product(grid)))
!   do i = 0, grid(1) - 1, 1
!      do j = 0, grid(2) - 1, 1
!         do k = 0, grid(3) - 1, 1
!            num_points = num_points + 1
!            gridpoints(:, num_points) = real([i + 0.5, j + 0.5, k + 0.5],real12) / real(grid, real12)
!         end do
!      end do
!   end do
  gridpoints = get_viable_gridpoints( generator%grid, &
       basis_host, &
       [ generator%distributions%bond_info(:)%radius_covalent ], &
       atom_ignore_list, &
       lowtol = generator%distributions%radius_distance_tol(1), &
       grid_offset = generator%grid_offset &
  )
  write(*,*) "Number of gridpoints: ", size(gridpoints, dim=2)


  !-----------------------------------------------------------------------------
  ! write distribution functions to file
  !-----------------------------------------------------------------------------
  call generator%distributions%write_2body(file="2body.txt")
  call generator%distributions%write_3body(file="3body.txt")
  call generator%distributions%write_4body(file="4body.txt")


  !-----------------------------------------------------------------------------
  ! check element energies
  !-----------------------------------------------------------------------------
  call generator%distributions%get_element_energies( &
       tmp_symbols, &
       tmp_vals &
  )
  do i = 1, size(tmp_symbols)
     write(*,*) "Element: ", tmp_symbols(i), " energy: ", tmp_vals(i)
  end do
  call generator%distributions%get_bond_radii( &
       tmp_pairs, &
       tmp_vals &
  )
  do i = 1, size(tmp_symbols)
     write(*,*) "Element1: ", tmp_pairs(i,1), "Element2: ", tmp_pairs(i,2), " radii: ", tmp_vals(i)
  end do



  !-----------------------------------------------------------------------------
  ! call evaluator
  !-----------------------------------------------------------------------------
  open(newunit=unit, file="evaluator.txt", status="replace")
  do i = 1, basis_host%spec(1)%num - 1
     write(unit, *) matmul(basis_host%spec(1)%atom(i, :), basis_host%lat)
  end do
!   do i = 1, basis_host%image_spec(1)%num
!      write(unit, *) matmul(basis_host%image_spec(1)%atom(i, :), basis_host%lat)
!   end do
  write(unit, *)

  write(*,*) "LOOK", basis_host%num_images
  write(*,*) "Number of images: ", basis_host%image_spec(1)%num
  allocate(suitability_grid(size(gridpoints,2)))
  do i = 1, size(gridpoints,dim=2)
     suitability_grid(i) = evaluate_point( generator%distributions, &
          gridpoints(:,i), basis_host, &
          atom_ignore_list, &
          [ generator%distributions%bond_info(:)%radius_covalent ] &
     )
     write(unit, *) matmul(gridpoints(:,i),basis_host%lat), suitability_grid(i)
  end do
  close(unit)
  write(*,*) "Size of suitability grid: ", size(suitability_grid)
  write(*,*) "Size of gridpoints: ", size(gridpoints,2)
  write(*,*) "Max suitability: ", maxval(suitability_grid)
  write(*,*) "Max location: ", maxloc(suitability_grid)
  best_loc = maxloc(suitability_grid,dim=1)
  write(*,*) gridpoints(:,best_loc)
  basis_host%spec(1)%atom(basis_host%spec(1)%num, 1:3) = &
       gridpoints(1:3,best_loc)
  call basis_host%update_images( &
       max_bondlength = generator%distributions%cutoff_max(1), &
       is = atom_ignore_list(1,1), &
       ia = atom_ignore_list(1,2) &
  )
  write(*,*) "Number of images after update: ", basis_host%image_spec(1)%num

  ! !-----------------------------------------------------------------------------
  ! ! generate random structures
  ! !-----------------------------------------------------------------------------
  ! write(*,*) "Generating structures"
  ! call generator%generate( num_structures, &
  !      stoich, &
  !      method_probab )
  ! write(*,*) "Structures have been successfully generated"


end program test_evaluator