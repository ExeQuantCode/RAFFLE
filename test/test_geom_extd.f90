program test_geom_extd
  !! Test program for the module geom_extd.
  use raffle__io_utils
  use raffle__constants, only: real32
  use raffle__geom_extd
  implicit none

  type(extended_basis_type) :: basis_diamond

  logical :: success = .true.

  test_error_handling = .true.

  ! diamond cell
  basis_diamond%nspec = 1
  basis_diamond%natom = 8
  allocate(basis_diamond%spec(basis_diamond%nspec))
  basis_diamond%spec(1)%num = 8
  basis_diamond%spec(1)%name = 'C'
  allocate(basis_diamond%spec(1)%atom(basis_diamond%spec(1)%num, 3))
  basis_diamond%spec(1)%atom(1, :3) = [0.0, 0.0, 0.0]
  basis_diamond%spec(1)%atom(2, :3) = [0.5, 0.5, 0.0]
  basis_diamond%spec(1)%atom(3, :3) = [0.5, 0.0, 0.5]
  basis_diamond%spec(1)%atom(4, :3) = [0.0, 0.5, 0.5]
  basis_diamond%spec(1)%atom(5, :3) = [0.25, 0.25, 0.25]
  basis_diamond%spec(1)%atom(6, :3) = [0.75, 0.75, 0.25]
  basis_diamond%spec(1)%atom(7, :3) = [0.75, 0.25, 0.75]
  basis_diamond%spec(1)%atom(8, :3) = [0.25, 0.75, 0.75]

  basis_diamond%lat(1,:) = [3.5607451090903233, 0.0, 0.0]
  basis_diamond%lat(2,:) = [0.0, 3.5607451090903233, 0.0]
  basis_diamond%lat(3,:) = [0.0, 0.0, 3.5607451090903233]
  basis_diamond%energy = -72.213492

  call test_create_images(basis_diamond, success)
  call test_update_iamges(basis_diamond, success)

  !-----------------------------------------------------------------------------
  ! check for any failed tests
  !-----------------------------------------------------------------------------
  write(*,*) "----------------------------------------"
  if(success)then
     write(*,*) 'test_geom_extd passed all tests'
  else
     write(0,*) 'test_geom_extd failed one or more tests'
     stop 1
  end if

contains

  subroutine test_create_images(basis, success)
    implicit none
    type(extended_basis_type), intent(in) :: basis
    logical, intent(inout) :: success

    type(extended_basis_type) :: basis_copy

    call basis_copy%copy(basis)

    ! Create images
    call basis_copy%create_images( max_bondlength = 0._real32 )

    ! Check if the number of images is correct
    call assert( &
         basis_copy%num_images .eq. 10, &
         'Number of images is not correct', &
         success &
    )
    call assert( &
         basis_copy%image_spec(1)%num .eq. 10, &
         'Number of images in first species is not correct', &
         success &
    )
    call assert( &
         all( &
              basis_copy%image_spec(1)%atom(:, :) - &
              1._real32 .lt. 1.E-6_real32 &
         ), &
         'Atoms outside of max bondlength', &
         success &
    )

  end subroutine test_create_images

  subroutine test_update_iamges(basis, success)
    implicit none
    type(extended_basis_type), intent(in) :: basis
    logical, intent(inout) :: success

    type(extended_basis_type) :: basis_copy
    integer, dimension(2,1) :: atom_ignore_list

    call basis_copy%copy(basis)
    atom_ignore_list(:,1) = [1, 1]

    ! Create images
    call basis_copy%create_images( &
         max_bondlength = 0._real32, &
         atom_ignore_list = atom_ignore_list &
    )

    ! Check if the number of images is correct
    call assert( &
         basis_copy%num_images .eq. 3, &
         'Number of images is not correct', &
         success &
    )
    call assert( &
         basis_copy%image_spec(1)%num .eq. 3, &
         'Number of images in first species is not correct', &
         success &
    )

    ! Update images
    call basis_copy%update_images( max_bondlength = 0._real32, is = 1, ia = 1 )

    ! Check if the number of images is correct
    call assert( &
         basis_copy%num_images .eq. 10, &
         'Number of images is not correct', &
         success &
    )
    call assert( &
         basis_copy%image_spec(1)%num .eq. 10, &
         'Number of images in first species is not correct', &
         success &
    )
    call assert( &
         all( &
              basis_copy%image_spec(1)%atom(:, :) - &
              1._real32 .lt. 1.E-6_real32 &
         ), &
         'Atoms outside of max bondlength', &
         success &
    )

  end subroutine test_update_iamges

!###############################################################################

  function compare_bas(bas1, bas2) result(output)
    type(extended_basis_type), intent(in) :: bas1, bas2
    logical :: output
    output = .true.

    ! Compare the geometries
    if(any(abs(bas1%lat - bas2%lat).gt.1.E-6)) then
       write(0,*) 'Geometry read/write failed, lattice check failed'
       output = .false.
    end if
    if(bas1%sysname .ne. bas2%sysname) then
       write(0,*) 'Geometry read/write failed, system name check failed'
       write(0,*) bas1%sysname, bas2%sysname
       output = .false.
    end if
    if(bas1%natom .ne. bas2%natom) then
       write(0,*) 'Geometry read/write failed, number of atoms check failed'
       write(0,*) bas1%natom, bas2%natom
       output = .false.
    end if
    if(abs(bas1%energy - bas2%energy).gt.1.E-6) then
       write(0,*) 'Geometry read/write failed, energy check failed'
       write(0,*) bas1%energy, bas2%energy
       output = .false.
    end if

  end function compare_bas

  subroutine uninitialise_bas(bas)
    implicit none
    type(extended_basis_type), intent(inout) :: bas

    bas%natom = 0
    bas%nspec = 0
    bas%lat = 0.E0
    bas%energy = 0.E0
    bas%sysname = ""
    bas%lcart = .false.
    bas%pbc = .true.
    deallocate(bas%spec)
    
  end subroutine uninitialise_bas

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

end program test_geom_extd