module raffle__dist_calcs
  !! Module containing distance calculators
  !!
  !! This module contains procedures to calculate the distance between atoms
  !! and other points in the system.
  use raffle__constants, only: pi,real32
  use raffle__geom_rw, only: basis_type
  use raffle__misc_linalg, only: modu, get_angle
  implicit none


  private

  public :: get_min_dist
  public :: get_min_dist_between_point_and_atom
  public :: get_min_dist_between_point_and_species
  public :: get_dist_between_point_and_atom


contains

!###############################################################################
  pure function get_min_dist(basis,loc,lignore_close,axis,labove,lreal,tol, &
       ignore_list) result(output)
    !! Return the minimum distance between a point and the nearest atom
    !! in a cell.
    !!
    !! This function returns the minimum distance between a point and the
    !! nearest atom in a periodic cell.
    implicit none

    ! Arguments
    logical, intent(in) :: lignore_close
    !! If true, ignore atoms that are really close to the point.
    class(basis_type), intent(in) :: basis
    !! The basis of the cell.
    real(real32), dimension(3), intent(in) :: loc
    !! The location of the point (in crystal coordinates).

    integer, intent(in), optional :: axis
    !! The axis along which to calculate the distance (if undefined, the
    !! distance is calculated in all directions).
    real(real32), intent(in), optional :: tol
    !! The tolerance for the distance.
    logical, intent(in), optional :: labove, lreal
    !! If true, return the real distance, otherwise return the vector.
    integer, dimension(:,:), intent(in), optional :: ignore_list
    !! List of atoms to ignore.
    real(real32), dimension(3) :: output
    !! The minimum distance between the point and the nearest atom.


    ! Local variables
    integer :: js, ja, i
    !! Loop counters.
    integer :: axis_
    !! Axis along which to calculate the distance.
    real(real32) :: dtmp1
    !! Temporary variables.
    real(real32) :: min_bond
    !! Minimum bond length.
    real(real32) :: tol_
    !! Tolerance for the distance.
    logical :: labove_, lreal_
    !! Booleans for above and real distance arguments
    real(real32), dimension(3) :: vdtmp1, vdtmp2
    !! Vectors for distance calculations.


    ! CORRECT tol TO ACCOUNT FOR LATTICE SIZE
    tol_ = 1.E-5_real32
    labove_ = .false.
    lreal_ = .true.
    axis_ = 0
    if(present(tol)) tol_ = tol

    if(present(labove)) labove_ = labove

    if(present(lreal)) lreal_ = lreal

    if(present(axis)) axis_=axis

    min_bond=huge(0._real32)
    output = 0._real32
    do js = 1, basis%nspec
       atmloop: do ja = 1, basis%spec(js)%num
          if(present(ignore_list))then
             do i = 1, size(ignore_list,2), 1
                if(all(ignore_list(:,i).eq.[js,ja])) cycle atmloop
             end do
          end if
          vdtmp1 = basis%spec(js)%atom(ja,:3) - loc
          if(lignore_close.and.modu(vdtmp1).lt.tol_) cycle atmloop
          if(axis_.gt.0)then
             if(abs(vdtmp1(axis_)).lt.tol_) cycle atmloop
             if(labove_)then
                vdtmp1(axis_) = 1._real32 + vdtmp1(axis_)
             else
                vdtmp1(axis_) = vdtmp1(axis_) - 1._real32
             end if
          else
             vdtmp1 = vdtmp1 - ceiling(vdtmp1 - 0.5_real32)
          end if
          vdtmp2 = matmul(vdtmp1,basis%lat)
          dtmp1 = modu(vdtmp2)
          if(dtmp1.lt.min_bond)then
             min_bond = dtmp1
             if(lreal_)then
                output = vdtmp2
             else
                output = vdtmp1
             end if
          end if
       end do atmloop
    end do

  end function get_min_dist
!###############################################################################


!###############################################################################
  pure function get_min_dist_between_point_and_atom(basis,loc,atom) result(dist)
    !! Return the minimum distance between a point and an atom in a cell.
    !!
    !! This function returns the minimum distance between a point and an atom
    !! in a periodic cell.
    implicit none

    ! Arguments
    class(basis_type), intent(in) :: basis
    !! The basis of the cell.
    integer, dimension(2), intent(in) :: atom
    !! The index of the atom in the cell (species, atom).
    real(real32), dimension(3), intent(in) :: loc
    !! The location of the point (in crystal coordinates).
    real(real32) :: dist
    !! The minimum distance between the point and the atom.

    ! Local variables
    real(real32), dimension(3) :: vec
    !! Vector between the point and the atom.

    vec = loc - basis%spec(atom(1))%atom(atom(2),:3)
    vec = vec - ceiling(vec - 0.5_real32)
    vec = matmul(vec,basis%lat)
    dist = modu(vec)

  end function get_min_dist_between_point_and_atom
!###############################################################################


!###############################################################################
  pure function get_min_dist_between_point_and_species( &
       basis, loc, species, ignore_list) result(dist)
    !! Return the minimum distance between a point and a species in a cell.
    !!
    !! This function returns the minimum distance between a point and any
    !! instance of the specified species in a periodic cell.
    implicit none

    ! Arguments
    class(basis_type), intent(in) :: basis
    !! The basis of the cell.
    integer, intent(in) :: species
    !! The index of the species in the cell.
    real(real32), dimension(3), intent(in) :: loc
    !! The location of the point (in crystal coordinates).
    integer, dimension(:,:), intent(in), optional :: ignore_list
    !! List of atoms to ignore.
    real(real32) :: dist
    !! The minimum distance between the point and the species.

    ! Local variables
    integer :: ia, i
    !! Loop indices.
    real(real32) :: rtmp1
    !! Temporary variable.
    real(real32), dimension(3) :: vec
    !! Vector between the point and the atom.


    dist = huge(0._real32)
    atom_loop: do ia = 1,basis%spec(species)%num
       if(present(ignore_list))then
          do i = 1, size(ignore_list,2), 1
             if(all(ignore_list(:,i).eq.[species,ia])) cycle atom_loop
          end do
       end if
       vec = loc - basis%spec(species)%atom(ia,:3)
       vec = vec - ceiling(vec - 0.5_real32)
       vec = matmul(vec, basis%lat)
       rtmp1 = modu(vec)
       if( rtmp1 .lt. dist ) dist = rtmp1
    end do atom_loop

  end function get_min_dist_between_point_and_species
!###############################################################################


!###############################################################################
  pure function get_dist_between_point_and_atom(basis,loc,atom) result(dist)
    !! Return the distance between a point and an atom in a cell.
    !!
    !! This function returns the distance between a point and an atom in a cell.
    implicit none

    ! Arguments
    class(basis_type), intent(in) :: basis
    !! The basis of the cell.
    integer, dimension(2), intent(in) :: atom
    !! The index of the atom in the cell (species, atom).
    real(real32), dimension(3), intent(in) :: loc
    !! The location of the point (in crystal coordinates).
    real(real32) :: dist
    !! The minimum distance between the point and the atom.

    ! Local variables
    real(real32), dimension(3) :: vec
    !! Vector between the point and the atom.

    vec = loc - basis%spec(atom(1))%atom(atom(2),:3)
    vec = matmul(vec,basis%lat)
    dist = modu(vec)

  end function get_dist_between_point_and_atom
!###############################################################################

end module raffle__dist_calcs
