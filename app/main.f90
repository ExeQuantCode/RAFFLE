program raffle_program
  use constants, only: real12
  use misc_raffle, only: touch
  use inputs
  use read_structures, only: get_evolved_gvectors_from_data
  use raffle, only: raffle_generator_type, gvector_container_type
  use rw_geom, only: geom_read, geom_write
  implicit none

  integer :: i, unit
  character(1024) :: buffer

  ! type(gvector_container_type) :: gvector_container
  real(real12), dimension(3) :: method_probab
  real(real12), dimension(:), allocatable :: tmp_energies
  character(len=3), dimension(:), allocatable :: tmp_symbols

  type(raffle_generator_type) :: generator



!!!-----------------------------------------------------------------------------
!!! read input file
!!!-----------------------------------------------------------------------------
  call set_global_vars()


!!!-----------------------------------------------------------------------------
!!! check the task and run the appropriate case
!!!-----------------------------------------------------------------------------
!!! OLD TASKS !!!
!!! 0) Run RSS
!!! 1) Regenerate DIst Files (WIP)
!!! 2) Run HOST_RSS
!!! 3) Test
!!! 4) Sphere_Overlap
!!! 5) Bondangle_test !!! THIS LITERALLY JUST TESTS THAT THE BONDANGLE METHOD WORKS! DO NOT USE! !!!
!!! 6) Run evo (Should be run after any set created)
!!! 7) Add new poscar  
!!! 8) Run evo, but don't regen energies or evolve distributions (only reformat gaussians) 
!!! 9) Run evo, don't get energies but do evolve distributions
  select case(task)
  case(0)
     write(*,*) "NOTHING WAS EVER SET UP FOR CASE 0"
  case(1)
     write(*,*) "Regenerating Distribution Files"
     write(*,*) "DEPRECATED"
     stop 0
  case(2)
     write(*,*) "Running HOST_RSS"
  case default
     write(*,*) "Invalid option"
     stop 1
  end select



!!!-----------------------------------------------------------------------------
!!! set the element energies
!!!-----------------------------------------------------------------------------
  call generator%distributions%set_element_energies( &
       element_symbols, &
       element_energies &
  )


!!!-----------------------------------------------------------------------------
!!! read structures from the database and generate gvectors
!!!-----------------------------------------------------------------------------
  generator%distributions = get_evolved_gvectors_from_data( &
       input_dir    = database_list, &
       file_format  = database_format, &
       gvector_container_template = gvector_container_type(&
            width = width_list, &
            sigma = sigma_list, &
            cutoff_min = cutoff_min_list, &
            cutoff_max = cutoff_max_list ) )

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


!!!-----------------------------------------------------------------------------
!!! calculate the probability of each placement method
!!!-----------------------------------------------------------------------------
  method_probab(1) = vps_ratio(1)/real(sum(vps_ratio),real12)
  method_probab(2) = method_probab(1) + &
       vps_ratio(2)/real(sum(vps_ratio),real12)
  method_probab(3) = method_probab(2) + &
       vps_ratio(3)/real(sum(vps_ratio),real12)
  write(*,*) "Method probabilities (void, scan, pseudorandom-walk): ", &
       method_probab


!!!-----------------------------------------------------------------------------
!!! set the host structure
!!!-----------------------------------------------------------------------------
  open(newunit=unit, file=filename_host, status='old')
  call geom_read(unit, generator%host)
  close(unit)
  generator%bins = bins


!!!-----------------------------------------------------------------------------
!!! generate random structures
!!!-----------------------------------------------------------------------------
  write(*,*) "Generating structures"
  call generator%generate( num_structures, &
       stoich, &
       method_probab )
  write(*,*) "Structures have been successfully generated"


!!!-----------------------------------------------------------------------------
!!! save generated structures
!!!-----------------------------------------------------------------------------
  do i = 1, generator%num_structures
     write(buffer,'(A,"/struc",I0.3)') trim(output_dir),i
     call touch(buffer)
     open(newunit = unit, file=trim(buffer)//"/POSCAR")
     call geom_write(unit, generator%structures(i))
     close(unit)
  end do

end program raffle_program