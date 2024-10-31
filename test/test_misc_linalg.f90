program test_misc_linalg
  use raffle__io_utils
  use raffle__misc_linalg
  use raffle__constants, only: real32, pi
  implicit none

  logical :: success = .true.

  test_error_handling = .true.


  call test_modu(success)
  call test_cross(success)
  call test_get_distance(success)
  call test_get_angle_from_vectors(success)
  call test_get_angle_from_points(success)
  call test_get_dihedral_angle_from_vectors(success)
  call test_get_dihedral_angle_from_points(success)
  call test_inverse_3x3(success)


  !-----------------------------------------------------------------------------
  ! check for any failed tests
  !-----------------------------------------------------------------------------
  write(*,*) "----------------------------------------"
  if(success)then
     write(*,*) 'test_misc_linalg passed all tests'
  else
     write(0,*) 'test_misc_linalg failed one or more tests'
     stop 1
  end if

contains

  subroutine test_modu(success)
    logical, intent(inout) :: success
    real(real32), dimension(3) :: vector
    real(real32) :: result
    vector = [3.0_real32, 4.0_real32, 0.0_real32]
    result = modu(vector)
    call assert_almost_equal_scalar( &
         result, 5.0_real32, 1.E-6_real32, &
         "modu", success &
    )
  end subroutine test_modu

  subroutine test_cross(success)
    logical, intent(inout) :: success
    real(real32), dimension(3) :: a, b, result
    a = [1.0_real32, 0.0_real32, 0.0_real32]
    b = [0.0_real32, 1.0_real32, 0.0_real32]
    result = cross(a, b)
    call assert_almost_equal_vector( &
         result, [0.0_real32, 0.0_real32, 1.0_real32], 1.E-6_real32, &
         "cross", success &
    )
  end subroutine test_cross

  subroutine test_get_distance(success)
    logical, intent(inout) :: success
    real(real32), dimension(3) :: point1, point2
    real(real32) :: result
    point1 = [1.0_real32, 2.0_real32, 3.0_real32]
    point2 = [4.0_real32, 6.0_real32, 8.0_real32]
    result = get_distance(point1, point2)
    call assert_almost_equal_scalar( &
         result, 7.0710678118654755_real32, 1.E-6_real32, &
         "get_angle_from_vectors", success &
    )
  end subroutine test_get_distance

  subroutine test_get_angle_from_vectors(success)
    logical, intent(inout) :: success
    real(real32), dimension(3) :: vector1, vector2
    real(real32) :: result
    vector1 = [1.0_real32, 0.0_real32, 0.0_real32]
    vector2 = [0.0_real32, 1.0_real32, 0.0_real32]
    result = get_angle(vector1, vector2)
    call assert_almost_equal_scalar( &
         result, pi/2.0_real32, 1.E-6_real32, &
         "get_angle_from_vectors", success &
    )
  end subroutine test_get_angle_from_vectors

  subroutine test_get_angle_from_points(success)
    logical, intent(inout) :: success
    real(real32), dimension(3) :: point1, point2, point3
    real(real32) :: result
    point1 = [1.0_real32, 0.0_real32, 0.0_real32]
    point2 = [0.0_real32, 0.0_real32, 0.0_real32]
    point3 = [0.0_real32, 1.0_real32, 0.0_real32]
    result = get_angle(point1, point2, point3)
    call assert_almost_equal_scalar( &
         result, pi/2.0_real32, 1.E-6_real32, &
         "get_angle_from_points", success &
    )
  end subroutine test_get_angle_from_points

  subroutine test_get_dihedral_angle_from_vectors(success)
    logical, intent(inout) :: success
    real(real32), dimension(3) :: vector1, vector2, vector3
    real(real32) :: result
    vector1 = [1.0_real32, 0.0_real32, 0.0_real32]
    vector2 = [0.0_real32, 1.0_real32, 0.0_real32]
    vector3 = [1.0_real32, 0.0_real32, 0.0_real32]
    result = get_dihedral_angle(vector1, vector2, vector3)
    call assert_almost_equal_scalar( &
         result, pi/2.0_real32, 1.E-6_real32, &
         "get_dihedral_angle_from_vectors", success &
    )
  end subroutine test_get_dihedral_angle_from_vectors

  subroutine test_get_dihedral_angle_from_points(success)
    logical, intent(inout) :: success
    real(real32), dimension(3) :: point1, point2, point3, point4
    real(real32) :: result
    point1 = [1.0_real32, 0.0_real32, 0.0_real32]
    point2 = [0.0_real32, 0.0_real32, 0.0_real32]
    point3 = [0.0_real32, 1.0_real32, 0.0_real32]
    point4 = [1.0_real32, 0.0_real32, .0_real32]
    result = get_dihedral_angle(point1, point2, point3, point4)
    call assert_almost_equal_scalar( &
         result, pi/2.0_real32, 1.E-6_real32, &
         "get_dihedral_angle_from_points", success &
    )
  end subroutine test_get_dihedral_angle_from_points

  subroutine test_inverse_3x3(success)
    logical, intent(inout) :: success
    real(real32), dimension(3,3) :: matrix, result, expected
    matrix = reshape([ &
         4.0_real32, 3.0_real32, 0.0_real32, &
         3.0_real32, 2.0_real32, 1.0_real32, &
         0.0_real32, 1.0_real32, 1.0_real32 &
    ], [3,3])
    expected = reshape([ &
         -1.0_real32, 3.0_real32, -3.0_real32, &
         3.0_real32, -4.0_real32, 4.0_real32, &
         -3.0_real32, 4.0_real32, 1.0_real32 &
    ], [3,3])
    expected = expected / 5.0_real32
    result = inverse_3x3(matrix)
    call assert_almost_equal_matrix( &
         result, expected, 1.E-6_real32, "inverse_3x3", success &
    )
  end subroutine test_inverse_3x3

  subroutine assert_almost_equal_scalar(actual, expected, tol, message, success)
    real(real32), intent(in) :: actual
    real(real32), intent(in) :: expected
    character(len=*), intent(in) :: message
    logical, intent(inout) :: success
    real(real32), intent(in) :: tol

    if( abs(actual - expected) .gt. tol ) then
       write(0,*) "Test failed: ", message
       success = .false.
    end if
  end subroutine assert_almost_equal_scalar

  subroutine assert_almost_equal_vector(actual, expected, tol, message, success)
    real(real32), dimension(:), intent(in) :: actual
    real(real32), dimension(..), intent(in) :: expected
    character(len=*), intent(in) :: message
    logical, intent(inout) :: success
    real(real32), intent(in) :: tol

    select rank(expected)
    rank(0)
       if( any( abs(actual - expected) .gt. tol ) ) then
          write(0,*) "Test failed: ", message
          success = .false.
       end if
    rank(1)
       if( any( abs(actual - expected) .gt. tol ) ) then
          write(0,*) "Test failed: ", message
          success = .false.
       end if
    end select
  end subroutine assert_almost_equal_vector

  subroutine assert_almost_equal_matrix(actual, expected, tol, message, success)
    real(real32), dimension(:,:), intent(in) :: actual
    real(real32), dimension(..), intent(in) :: expected
    character(len=*), intent(in) :: message
    logical, intent(inout) :: success
    real(real32), intent(in) :: tol

    select rank(expected)
    rank(0)
       if( any( abs(actual - expected) .gt. tol ) ) then
          write(0,*) "Test failed: ", message
          success = .false.
       end if
    rank(2)
       if( any( abs(actual - expected) .gt. tol ) ) then
          write(0,*) "Test failed: ", message
          success = .false.
       end if
    end select
  end subroutine assert_almost_equal_matrix

end program test_misc_linalg