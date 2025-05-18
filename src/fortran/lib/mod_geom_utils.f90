module raffle__geom_utils
  !! Module to contain all geometry-manipulation related procedures
  !!
  !! This module contains procedures that are used to manipulate the geometry
  !! of the system. The geometry type used is defined in the geom_rw module.
  use raffle__constants, only: pi,real32
  use raffle__geom_rw, only: basis_type
  use raffle__misc_linalg, only: modu, get_angle
  implicit none


  private

  public :: basis_merge


contains

!###############################################################################
  function basis_merge(basis1,basis2,length,map1,map2, mask1, mask2) &
       result(output)
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
    logical, intent(in), optional :: mask1, mask2
    !! Mask for atoms in the two bases.

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
    do i = 1, basis2%nspec
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
    spec1check: do i = 1, basis2%nspec
       do j = 1, basis1%nspec
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
    do i = 1, basis1%nspec
       allocate(output%spec(i)%atom_mask(output%spec(i)%num))
       allocate(output%spec(i)%atom_idx(output%spec(i)%num))
       allocate(output%spec(i)%atom(output%spec(i)%num,dim))
       output%spec(i)%atom_mask(1:basis1%spec(i)%num) = basis1%spec(i)%atom_mask
       output%spec(i)%atom(:,:)=0._real32
       output%spec(i)%atom(1:basis1%spec(i)%num,:3) = basis1%spec(i)%atom(:,:3)
       output%spec(i)%atom_idx(1:basis1%spec(i)%num) = basis1%spec(i)%atom_idx
       if(lmap) new_map(i,:basis1%spec(i)%num,:) = map1(i,:basis1%spec(i)%num,:)
       if(present(mask1)) output%spec(i)%atom_mask(1:basis1%spec(i)%num) = mask1
    end do
    do i = 1, basis2%nspec
       if(match(i).gt.basis1%nspec)then
          allocate(output%spec(match(i))%atom_mask(output%spec(match(i))%num))
          output%spec(match(i))%atom_mask(:) = basis2%spec(i)%atom_mask(:)
          allocate(output%spec(match(i))%atom_idx(output%spec(match(i))%num))
          output%spec(match(i))%atom_idx(:) = &
               basis2%spec(i)%atom_idx(:) + basis1%natom
          allocate(output%spec(match(i))%atom(output%spec(match(i))%num,dim))
          output%spec(match(i))%atom(:,:)=0._real32
          output%spec(match(i))%atom(:,:3)=basis2%spec(i)%atom(:,:3)
          if(lmap) new_map(match(i),:basis2%spec(i)%num,:) = &
               map2(i,:basis2%spec(i)%num,:)
          if(present(mask2)) output%spec(match(i))%atom_mask(:) = mask2
       else
          itmp=basis1%spec(match(i))%num
          output%spec(match(i))%atom_mask(itmp+1:basis2%spec(i)%num+itmp) = &
               basis2%spec(i)%atom_mask(:)
          output%spec(match(i))%atom_idx(itmp+1:basis2%spec(i)%num+itmp) = &
               basis2%spec(i)%atom_idx(:) + basis1%natom
          output%spec(match(i))%atom(itmp+1:basis2%spec(i)%num+itmp,:3) = &
               basis2%spec(i)%atom(:,:3)   
          if(lmap) new_map(match(i),itmp+1:basis2%spec(i)%num+itmp,:) = &
               map2(i,:basis2%spec(i)%num,:)    
          if(present(mask2)) &
               output%spec(match(i))%atom_mask( &
                    itmp+1:basis2%spec(i)%num+itmp &
               ) = mask2  
       end if
    end do
    output%natom=sum(output%spec(:)%num)

    if(lmap) call move_alloc(new_map,map1)

  end function basis_merge
!###############################################################################

end module raffle__geom_utils
