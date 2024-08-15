program test_evaluator
  use constants, only: real12
  use rw_geom, only: basis_type
  use extended_geom, only: extended_basis_type
  use evaluator, only: evaluate_point_multiplier
  use generator, only: raffle_generator_type
  ! use add_atom, only: get_viable_gridpoints
  implicit none


  integer :: unit
  integer :: i, j, k, num_points
  type(extended_basis_type) :: basis_host
  type(basis_type), dimension(1) :: database
  integer, dimension(3) :: grid
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

  database(1)%lat(1,:) = [3.5668, 0.0, 0.0]
  database(1)%lat(2,:) = [0.0, 3.5668, 0.0]
  database(1)%lat(3,:) = [0.0, 0.0, 3.5668]
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
  basis_host%nspec = 1
  basis_host%natom = 9
  allocate(basis_host%spec(basis_host%nspec))
  basis_host%spec(1)%num = 9
  basis_host%spec(1)%name = 'C'
  allocate(basis_host%spec(1)%atom(basis_host%spec(1)%num, 3))
  basis_host%spec(1)%atom(1, :3) = [0.0, 0.0, 0.0]
  basis_host%spec(1)%atom(2, :3) = [0.5, 0.5, 0.0]
  basis_host%spec(1)%atom(3, :3) = [0.5, 0.0, 0.25]
  basis_host%spec(1)%atom(4, :3) = [0.0, 0.5, 0.25]
  basis_host%spec(1)%atom(5, :3) = [0.25, 0.25, 0.125]
  basis_host%spec(1)%atom(6, :3) = [0.75, 0.75, 0.125]
  basis_host%spec(1)%atom(7, :3) = [0.75, 0.25, 0.375]
  basis_host%spec(1)%atom(8, :3) = [0.25, 0.75, 0.375]
  basis_host%spec(1)%atom(9, :3) = [0.0, 0.0, 0.0]
  atom_ignore_list(1,1) = 1
  atom_ignore_list(1,2) = 9

  basis_host%lat(1,:) = [3.5668, 0.0, 0.0]
  basis_host%lat(2,:) = [0.0, 3.5668, 0.0]
  basis_host%lat(3,:) = [0.0, 0.0, 7.1336]
  call basis_host%create_images( &
       max_bondlength = 6._real12, &
       atom_ignore_list = atom_ignore_list &
  )


  grid = [10, 10, 20]
  call generator%host%copy(basis_host)
  generator%bins = grid

  num_points = 0
  allocate(gridpoints(3, product(grid)))
  do i = 0, grid(1) - 1, 1
     do j = 0, grid(2) - 1, 1
        do k = 0, grid(3) - 1, 1
           num_points = num_points + 1
           gridpoints(:, num_points) = real([i, j, k],real12) / real(grid, real12)
        end do
     end do
  end do
  ! gridpoints = get_viable_gridpoints( grid, &
  !      basis_host, &
  !      [ this%distributions%bond_info(:)%radius_covalent ], &
  !      atom_ignore_list &
  ! )


  !-----------------------------------------------------------------------------
  ! set up distribution functions
  !-----------------------------------------------------------------------------
  call generator%distributions%create( &
       basis_list = database, &
       deallocate_systems = .false. &
  )

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
  ! do i = 1, basis_host%image_spec(1)%num
  !    write(unit, *) matmul(basis_host%image_spec(1)%atom(i, :), basis_host%lat)
  ! end do
  write(unit, *)

  allocate(suitability_grid(product(grid)))
  do i = 1, size(gridpoints,dim=2)
     suitability_grid(i) = evaluate_point_multiplier( generator%distributions, &
          gridpoints(:,i), basis_host, &
          atom_ignore_list, &
          [ generator%distributions%bond_info(:)%radius_covalent ], &
          uptol=3._real12, lowtol=1.5_real12)
     write(unit, *) matmul(gridpoints(:,i),basis_host%lat), suitability_grid(i)
  end do
  close(unit)

  ! !-----------------------------------------------------------------------------
  ! ! generate random structures
  ! !-----------------------------------------------------------------------------
  ! write(*,*) "Generating structures"
  ! call generator%generate( num_structures, &
  !      stoich, &
  !      method_probab )
  ! write(*,*) "Structures have been successfully generated"


end program test_evaluator