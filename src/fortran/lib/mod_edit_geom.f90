module edit_geom
  !! Module to contain all geometry-manipulation related procedures
  !!
  !! This module contains procedures that are used to manipulate the geometry
  !! of the system. The geometry type used is defined in the rw_geom module.
  use constants, only: pi,real12
  use rw_geom, only: basis_type
  use misc_linalg, only: modu, get_angle
  implicit none


  private

  public :: get_min_dist
  public :: get_min_dist_between_point_and_atom
  public :: bas_merge


contains

!###############################################################################
  function get_min_dist(bas,loc,lignore_close,axis,labove,lreal,tol, &
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
    type(basis_type), intent(in) :: bas
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
    real(real12) :: tol_ = 1.E-5
    !! Tolerance for the distance.
    logical :: labove_ = .false., lreal_ = .true.
    !! Booleans for above and real distance arguments
    real(real12), dimension(3) :: vdtmp1, vdtmp2
    !! Vectors for distance calculations.


    ! CORRECT tol TO ACCOUNT FOR LATTICE SIZE
    if(present(tol)) tol_ = tol

    if(present(labove)) labove_ = labove

    if(present(lreal)) lreal_ = lreal

    if(present(axis)) axis_=axis

    min_bond=huge(0._real12)
    output = 0._real12
    do js = 1, bas%nspec
       atmloop: do ja=1,bas%spec(js)%num
          if(present(ignore_list))then
             do i = 1, size(ignore_list,1), 1
                if(all(ignore_list(i,:).eq.[js,ja])) cycle atmloop
             end do
          end if
          vdtmp1 = bas%spec(js)%atom(ja,:3) - loc
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
          vdtmp2 = matmul(vdtmp1,bas%lat)
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
  pure function get_min_dist_between_point_and_atom(bas,loc,atom) &
       result(dist)
    !! Return the minimum distance between a point and an atom in a cell.
    !!
    !! This function returns the minimum distance between a point and an atom
    !! in a periodic cell.
    implicit none

    ! Arguments
    type(basis_type), intent(in) :: bas
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

    vec = loc - bas%spec(atom(1))%atom(atom(2),:3)
    vec = vec - ceiling(vec - 0.5_real12)
    vec = matmul(vec,bas%lat)
    dist = modu(vec)

  end function get_min_dist_between_point_and_atom
!###############################################################################


!###############################################################################
  function bas_merge(bas1,bas2,length,map1,map2) result(mergbas)
    !! Merge two supplied bases
    !!
    !! Merge two bases assuming that the lattice is the same
    implicit none

    ! Arguments
    type(basis_type) :: mergbas
    !! Output merged basis.
    type(basis_type), intent(in) :: bas1, bas2
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

    allocate(match(bas2%nspec))
    match=0
    mergbas%nspec=bas1%nspec
    do i=1,bas2%nspec
       if(.not.any(bas2%spec(i)%name.eq.bas1%spec(:)%name))then
          mergbas%nspec=mergbas%nspec+1
       end if
    end do
    allocate(mergbas%spec(mergbas%nspec))
    mergbas%spec(:bas1%nspec)%num=bas1%spec(:)%num
    mergbas%spec(:bas1%nspec)%name=bas1%spec(:)%name


    write(mergbas%sysname,'(A,"+",A)') &
         trim(bas1%sysname),trim(bas2%sysname)
    k=bas1%nspec
    spec1check: do i=1,bas2%nspec
       do j=1,bas1%nspec
          if(bas2%spec(i)%name.eq.bas1%spec(j)%name)then
             mergbas%spec(j)%num=mergbas%spec(j)%num+bas2%spec(i)%num
             match(i)=j
             cycle spec1check
          end if
       end do
       k=k+1
       match(i)=k
       mergbas%spec(k)%num=bas2%spec(i)%num
       mergbas%spec(k)%name=bas2%spec(i)%name
    end do spec1check


    !---------------------------------------------------------------------------
    ! if map is present, sets up new map
    !---------------------------------------------------------------------------
    lmap = .false.
    if_map: if(present(map1).and.present(map2))then
       if(all(map1.eq.-1)) exit if_map
       lmap = .true.
       allocate(new_map(&
            mergbas%nspec,&
            maxval(mergbas%spec(:)%num,dim=1),2))
       new_map = 0
    end if if_map


    !---------------------------------------------------------------------------
    ! set up atoms in merged basis
    !---------------------------------------------------------------------------
    do i=1,bas1%nspec
       allocate(mergbas%spec(i)%atom(mergbas%spec(i)%num,dim))
       mergbas%spec(i)%atom(:,:)=0._real12
       mergbas%spec(i)%atom(1:bas1%spec(i)%num,:3)=bas1%spec(i)%atom(:,:3)
       if(lmap) new_map(i,:bas1%spec(i)%num,:)=map1(i,:bas1%spec(i)%num,:)
    end do
    do i=1,bas2%nspec
       if(match(i).gt.bas1%nspec)then
          allocate(mergbas%spec(match(i))%atom(mergbas%spec(match(i))%num,dim))
          mergbas%spec(match(i))%atom(:,:)=0._real12
          mergbas%spec(match(i))%atom(:,:3)=bas2%spec(i)%atom(:,:3)
          if(lmap) new_map(match(i),:bas2%spec(i)%num,:) = &
               map2(i,:bas2%spec(i)%num,:)
       else
          itmp=bas1%spec(match(i))%num
          mergbas%spec(match(i))%atom(itmp+1:bas2%spec(i)%num+itmp,:3) = &
               bas2%spec(i)%atom(:,:3)   
          if(lmap) new_map(match(i),itmp+1:bas2%spec(i)%num+itmp,:) = &
               map2(i,:bas2%spec(i)%num,:)      
       end if
    end do
    mergbas%natom=sum(mergbas%spec(:)%num)


    if(lmap) call move_alloc(new_map,map1)

    return
  end function bas_merge
!###############################################################################

end module edit_geom
