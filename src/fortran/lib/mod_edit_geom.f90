module edit_geom
  !! Module to contain all geometry-manipulation related procedures
  !!
  !! This module contains procedures that are used to manipulate the geometry
  !! of the system. The geometry type used is defined in the rw_geom module.
  use raffle__constants, only: pi,real12
  use raffle__rw_geom, only: basis_type
  use raffle__misc_linalg, only: modu, get_angle
  implicit none


  private

  public :: get_min_dist
  public :: get_min_dist_between_point_and_atom
  public :: get_min_dist_between_point_and_species
  public :: get_dist_between_point_and_atom
  public :: basis_merge


contains

!###############################################################################
  function get_min_dist(basis,loc,lignore_close,axis,labove,lreal,tol, &
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
    real(real12), dimension(3), intent(in) :: loc
    !! The location of the point (in crystal coordinates).

    integer, intent(in), optional :: axis
    !! The axis along which to calculate the distance (if undefined, the
    !! distance is calculated in all directions).
    real(real12), intent(in), optional :: tol
    !! The tolerance for the distance.
    logical, intent(in), optional :: labove, lreal
    !! If true, return the real distance, otherwise return the vector.
    integer, dimension(:,:), intent(in), optional :: ignore_list
    !! List of atoms to ignore.
    real(real12), dimension(3) :: output
    !! The minimum distance between the point and the nearest atom.


    ! Local variables
    integer :: js, ja, i
    !! Loop counters.
    integer :: axis_ = 0
    !! Axis along which to calculate the distance.
    real(real12) :: dtmp1
    !! Temporary variables.
    real(real12) :: min_bond
    !! Minimum bond length.
    real(real12) :: tol_
    !! Tolerance for the distance.
    logical :: labove_, lreal_
    !! Booleans for above and real distance arguments
    real(real12), dimension(3) :: vdtmp1, vdtmp2
    !! Vectors for distance calculations.


    ! CORRECT tol TO ACCOUNT FOR LATTICE SIZE
    tol_ = 1.E-5_real12
    labove_ = .false.
    lreal_ = .true.
    if(present(tol)) tol_ = tol

    if(present(labove)) labove_ = labove

    if(present(lreal)) lreal_ = lreal

    if(present(axis)) axis_=axis

    min_bond=huge(0._real12)
    output = 0._real12
    do js = 1, basis%nspec
       atmloop: do ja=1,basis%spec(js)%num
          if(present(ignore_list))then
             do i = 1, size(ignore_list,1), 1
                if(all(ignore_list(i,:).eq.[js,ja])) cycle atmloop
             end do
          end if
          vdtmp1 = basis%spec(js)%atom(ja,:3) - loc
          if(lignore_close.and.modu(vdtmp1).lt.tol_) cycle atmloop
          if(axis_.gt.0)then
             if(abs(vdtmp1(axis_)).lt.tol_) cycle atmloop
             if(labove_)then
                vdtmp1(axis_) = 1._real12 + vdtmp1(axis_)
             else
                vdtmp1(axis_) = vdtmp1(axis_) - 1._real12
             end if
          else
             vdtmp1 = vdtmp1 - ceiling(vdtmp1 - 0.5_real12)
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
  pure function get_min_dist_between_point_and_atom(basis,loc,atom) &
       result(dist)
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
    real(real12), dimension(3), intent(in) :: loc
    !! The location of the point (in crystal coordinates).
    real(real12) :: dist
    !! The minimum distance between the point and the atom.

    ! Local variables
    real(real12), dimension(3) :: vec
    !! Vector between the point and the atom.

    vec = loc - basis%spec(atom(1))%atom(atom(2),:3)
    vec = vec - ceiling(vec - 0.5_real12)
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
    real(real12), dimension(3), intent(in) :: loc
    !! The location of the point (in crystal coordinates).
    integer, dimension(:,:), intent(in), optional :: ignore_list
    !! List of atoms to ignore.
    real(real12) :: dist
    !! The minimum distance between the point and the species.

    ! Local variables
    integer :: ia, i
    !! Loop indices.
    real(real12) :: rtmp1
    !! Temporary variable.
    real(real12), dimension(3) :: vec
    !! Vector between the point and the atom.


    dist = huge(0._real12)
    atom_loop: do ia = 1,basis%spec(species)%num
       if(present(ignore_list))then
          do i = 1, size(ignore_list,1), 1
             if(all(ignore_list(i,:).eq.[species,ia])) cycle atom_loop
          end do
       end if
       vec = loc - basis%spec(species)%atom(ia,:3)
       vec = vec - ceiling(vec - 0.5_real12)
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
    real(real12), dimension(3), intent(in) :: loc
    !! The location of the point (in crystal coordinates).
    real(real12) :: dist
    !! The minimum distance between the point and the atom.

    ! Local variables
    real(real12), dimension(3) :: vec
    !! Vector between the point and the atom.

    vec = loc - basis%spec(atom(1))%atom(atom(2),:3)
    vec = matmul(vec,basis%lat)
    dist = modu(vec)

  end function get_dist_between_point_and_atom
!###############################################################################


!###############################################################################
  function basis_merge(basis1,basis2,length,map1,map2) result(output)
    !! Merge two supplied bases
    !!
    !! Merge two bases assuming that the lattice is the same
    implicit none

    ! Arguments
    type(basis_type) :: output
    !! Output merged basis.
    class(basis_type), intent(in) :: basis1, basis2
    !! Input bases to merge.
    integer, intent(in), optional :: length
    !! Number of dimensions for atomic positions (default 3).
    integer, allocatable, dimension(:,:,:), optional, intent(inout) :: map1,map2
    !! Maps for atoms in the two bases.

    ! Local variables
    integer :: i, j, k, itmp, dim
    !! Loop counters.
    logical :: lmap
    !! Boolean for map presence.
    integer, allocatable, dimension(:) :: match
    !! Array to match species.
    integer, allocatable, dimension(:,:,:) :: new_map
    !! New map for merged basis.



    !---------------------------------------------------------------------------
    ! set up number of species
    !---------------------------------------------------------------------------
    dim=3
    if(present(length)) dim=length

    allocate(match(basis2%nspec))
    match=0
    output%nspec=basis1%nspec
    do i=1,basis2%nspec
       if(.not.any(basis2%spec(i)%name.eq.basis1%spec(:)%name))then
          output%nspec=output%nspec+1
       end if
    end do
    allocate(output%spec(output%nspec))
    output%spec(:basis1%nspec)%num=basis1%spec(:)%num
    output%spec(:basis1%nspec)%name=basis1%spec(:)%name


    write(output%sysname,'(A,"+",A)') &
         trim(basis1%sysname),trim(basis2%sysname)
    k=basis1%nspec
    spec1check: do i=1,basis2%nspec
       do j=1,basis1%nspec
          if(basis2%spec(i)%name.eq.basis1%spec(j)%name)then
             output%spec(j)%num=output%spec(j)%num+basis2%spec(i)%num
             match(i)=j
             cycle spec1check
          end if
       end do
       k=k+1
       match(i)=k
       output%spec(k)%num=basis2%spec(i)%num
       output%spec(k)%name=basis2%spec(i)%name
    end do spec1check


    !---------------------------------------------------------------------------
    ! if map is present, sets up new map
    !---------------------------------------------------------------------------
    lmap = .false.
    if_map: if(present(map1).and.present(map2))then
       if(all(map1.eq.-1)) exit if_map
       lmap = .true.
       allocate(new_map(&
            output%nspec,&
            maxval(output%spec(:)%num,dim=1),2))
       new_map = 0
    end if if_map


    !---------------------------------------------------------------------------
    ! set up atoms in merged basis
    !---------------------------------------------------------------------------
    do i=1,basis1%nspec
       allocate(output%spec(i)%atom(output%spec(i)%num,dim))
       output%spec(i)%atom(:,:)=0._real12
       output%spec(i)%atom(1:basis1%spec(i)%num,:3)=basis1%spec(i)%atom(:,:3)
       if(lmap) new_map(i,:basis1%spec(i)%num,:)=map1(i,:basis1%spec(i)%num,:)
    end do
    do i=1,basis2%nspec
       if(match(i).gt.basis1%nspec)then
          allocate(output%spec(match(i))%atom(output%spec(match(i))%num,dim))
          output%spec(match(i))%atom(:,:)=0._real12
          output%spec(match(i))%atom(:,:3)=basis2%spec(i)%atom(:,:3)
          if(lmap) new_map(match(i),:basis2%spec(i)%num,:) = &
               map2(i,:basis2%spec(i)%num,:)
       else
          itmp=basis1%spec(match(i))%num
          output%spec(match(i))%atom(itmp+1:basis2%spec(i)%num+itmp,:3) = &
               basis2%spec(i)%atom(:,:3)   
          if(lmap) new_map(match(i),itmp+1:basis2%spec(i)%num+itmp,:) = &
               map2(i,:basis2%spec(i)%num,:)      
       end if
    end do
    output%natom=sum(output%spec(:)%num)


    if(lmap) call move_alloc(new_map,map1)

    return
  end function basis_merge
!###############################################################################

end module edit_geom
