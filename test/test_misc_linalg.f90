program test_misc_linalg
  use error_handling
  use misc_linalg
  use constants, only: real12, pi
  implicit none

  logical :: success = .true.

  test_error_handling = .true.


  call test_uvec(success)
  call test_modu(success)
  call test_cross(success)
  call test_get_distance(success)
  call test_get_angle_from_vectors(success)
  call test_get_angle_from_points(success)
  call test_get_dihedral_angle_from_vectors(success)
  call test_get_dihedral_angle_from_points(success)
  call test_get_area(success)
  call test_get_vol(success)
  call test_LUinv(success)


  !-----------------------------------------------------------------------------
  ! check for any failed tests
  !-----------------------------------------------------------------------------
  write(*,*) "----------------------------------------"
  if(success)then
     write(*,*) 'test_misc_maths passed all tests'
  else
     write(0,*) 'test_misc_maths failed one or more tests'
     stop 1
  end if

contains

  subroutine test_uvec(success)
    logical, intent(inout) :: success
    real(real12), dimension(3) :: vector, result
    vector = [3.0_real12, 4.0_real12, 0.0_real12]
    result = uvec(vector)
    call assert_almost_equal_vector( &
         result, [0.6_real12, 0.8_real12, 0.0_real12], 1.E-6_real12, &
         "uvec", success &
    )
  end subroutine test_uvec

  subroutine test_modu(success)
    logical, intent(inout) :: success
    real(real12), dimension(3) :: vector
    real(real12) :: result
    vector = [3.0_real12, 4.0_real12, 0.0_real12]
    result = modu(vector)
    call assert_almost_equal_scalar( &
         result, 5.0_real12, 1.E-6_real12, &
         "modu", success &
    )
  end subroutine test_modu

  subroutine test_cross(success)
    logical, intent(inout) :: success
    real(real12), dimension(3) :: a, b, result
    a = [1.0_real12, 0.0_real12, 0.0_real12]
    b = [0.0_real12, 1.0_real12, 0.0_real12]
    result = cross(a, b)
    call assert_almost_equal_vector( &
         result, [0.0_real12, 0.0_real12, 1.0_real12], 1.E-6_real12, &
         "cross", success &
    )
  end subroutine test_cross

  subroutine test_get_distance(success)
    logical, intent(inout) :: success
    real(real12), dimension(3) :: point1, point2
    real(real12) :: result
    point1 = [1.0_real12, 2.0_real12, 3.0_real12]
    point2 = [4.0_real12, 6.0_real12, 8.0_real12]
    result = get_distance(point1, point2)
    call assert_almost_equal_scalar( &
         result, 7.0710678118654755_real12, 1.E-6_real12, &
         "get_angle_from_vectors", success &
    )
  end subroutine test_get_distance

  subroutine test_get_angle_from_vectors(success)
    logical, intent(inout) :: success
    real(real12), dimension(3) :: vector1, vector2
    real(real12) :: result
    vector1 = [1.0_real12, 0.0_real12, 0.0_real12]
    vector2 = [0.0_real12, 1.0_real12, 0.0_real12]
    result = get_angle(vector1, vector2)
    call assert_almost_equal_scalar( &
         result, pi/2.0_real12, 1.E-6_real12, &
         "get_angle_from_vectors", success &
    )
  end subroutine test_get_angle_from_vectors

  subroutine test_get_angle_from_points(success)
    logical, intent(inout) :: success
    real(real12), dimension(3) :: point1, point2, point3
    real(real12) :: result
    point1 = [1.0_real12, 0.0_real12, 0.0_real12]
    point2 = [0.0_real12, 0.0_real12, 0.0_real12]
    point3 = [0.0_real12, 1.0_real12, 0.0_real12]
    result = get_angle(point1, point2, point3)
    call assert_almost_equal_scalar( &
         result, pi/2.0_real12, 1.E-6_real12, &
         "get_angle_from_points", success &
    )
  end subroutine test_get_angle_from_points

  subroutine test_get_dihedral_angle_from_vectors(success)
    logical, intent(inout) :: success
    real(real12), dimension(3) :: vector1, vector2, vector3
    real(real12) :: result
    vector1 = [1.0_real12, 0.0_real12, 0.0_real12]
    vector2 = [0.0_real12, 1.0_real12, 0.0_real12]
    vector3 = [1.0_real12, 0.0_real12, 0.0_real12]
    result = get_dihedral_angle(vector1, vector2, vector3)
    write(*,*) "dihedral from vectors", result
    call assert_almost_equal_scalar( &
         result, pi/2.0_real12, 1.E-6_real12, &
         "get_dihedral_angle_from_vectors", success &
    )
  end subroutine test_get_dihedral_angle_from_vectors

  subroutine test_get_dihedral_angle_from_points(success)
    logical, intent(inout) :: success
    real(real12), dimension(3) :: point1, point2, point3, point4
    real(real12) :: result
    point1 = [1.0_real12, 0.0_real12, 0.0_real12]
    point2 = [0.0_real12, 0.0_real12, 0.0_real12]
    point3 = [0.0_real12, 1.0_real12, 0.0_real12]
    point4 = [1.0_real12, 0.0_real12, .0_real12]
    result = get_dihedral_angle(point1, point2, point3, point4)
    call assert_almost_equal_scalar( &
         result, pi/2.0_real12, 1.E-6_real12, &
         "get_dihedral_angle_from_points", success &
    )
  end subroutine test_get_dihedral_angle_from_points

  subroutine test_get_area(success)
    logical, intent(inout) :: success
    real(real12), dimension(3) :: a, b
    real(real12) :: result
    a = [1.0_real12, 0.0_real12, 0.0_real12]
    b = [0.0_real12, 1.0_real12, 0.0_real12]
    result = get_area(a, b)
    call assert_almost_equal_scalar( &
         result, 1.0_real12, 1.E-6_real12, "get_area", success &
    )
  end subroutine test_get_area

  subroutine test_get_vol(success)
    logical, intent(inout) :: success
    real(real12), dimension(3,3) :: matrix
    real(real12) :: result
    matrix = reshape([1.0_real12, 0.0_real12, 0.0_real12, &
                      0.0_real12, 1.0_real12, 0.0_real12, &
                      0.0_real12, 0.0_real12, 1.0_real12], [3,3])
    result = get_vol(matrix)
    call assert_almost_equal_scalar( &
         result, 1.0_real12, 1.E-6_real12, "get_vol", success &
    )
  end subroutine test_get_vol

  subroutine test_LUinv(success)
    logical, intent(inout) :: success
    real(real12), dimension(3,3) :: matrix, result, expected
    matrix = reshape([4.0_real12, 3.0_real12, 0.0_real12, &
                      3.0_real12, 2.0_real12, 1.0_real12, &
                      0.0_real12, 1.0_real12, 1.0_real12], [3,3])
    expected = reshape([-1.0_real12, 3.0_real12, -3.0_real12, &
                        3.0_real12, -4.0_real12, 4.0_real12, &
                        -3.0_real12, 4.0_real12, 1.0_real12], [3,3])
    expected = expected / 5.0_real12
    result = LUinv(matrix)
    call assert_almost_equal_matrix( &
         result, expected, 1.E-6_real12, "LUinv", success &
    )
  end subroutine test_LUinv

  subroutine assert_almost_equal_scalar(actual, expected, tol, message, success)
    real(real12), intent(in) :: actual
    real(real12), intent(in) :: expected
    character(len=*), intent(in) :: message
    logical, intent(inout) :: success
    real(real12), intent(in) :: tol

    if( abs(actual - expected) .gt. tol ) then
       write(0,*) "Test failed: ", message
       success = .false.
    end if
  end subroutine assert_almost_equal_scalar

  subroutine assert_almost_equal_vector(actual, expected, tol, message, success)
    real(real12), dimension(:), intent(in) :: actual
    real(real12), dimension(..), intent(in) :: expected
    character(len=*), intent(in) :: message
    logical, intent(inout) :: success
    real(real12), intent(in) :: tol

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
    real(real12), dimension(:,:), intent(in) :: actual
    real(real12), dimension(..), intent(in) :: expected
    character(len=*), intent(in) :: message
    logical, intent(inout) :: success
    real(real12), intent(in) :: tol

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