!!!#############################################################################
!!! Code written by Ned Thaddeus Taylor and Francis Huw Davies
!!! Code part of the ARTEMIS group (Hepplestone research group).
!!! Think Hepplestone, think HRG.
!!!#############################################################################
!!!module contains lattice- and basis-related functions and subroutines.
!!!module includes the following functions and subroutines:
!!! min_dist           (min distance between a point in a cell and nearest atom)
!!! min_dist_iterative (get min length of vector in a lattice, iterative method)
!!! get_surface_normal (return the vector normal to the surface of the plane ...
!!!                     ... constructed by the other two vectors)
!!! get_atom_height    (get the value of the atom along that axis)
!!! get_min_bulk_bond  (get the minimum bulk bond value)
!!! get_min_dist       (min distance between a point in a cell and nearest atom)
!!! get_centre_atom    (returns the location of the centre atom)
!!! get_closest_atom   (returns the closest atom to a height/3D coordinates)
!!! get_largest_gap    (returns the largest gap along a specified axis)
!!! get_shortest_bond  (returns the shortest bond for a specified atom)
!!! get_min_dist_between_two_atoms
!!! get_min_dist_between_point_and_atom
!!! centre_of_geom     (prints centre of geom of a molecule
!!! centre_of_mass     (prints centre of mass of a molecule)
!!!##################
!!! get_neighbour_fingerprint    (returns a fingerprint of the neighbours)
!!! get_nearest_neighbours_atom  (returns the nearest neighbours of an atom)
!!! get_nearest_neighbours_basis (returns the nearest neighbours of a basis)
!!!##################
!!! get_gvector_2body  (returns the g-vector for a 2-body interaction)
!!! get_gvector_3body  (returns the g-vector for a 3-body interaction)
!!! get_gvector_4body  (returns the g-vector for a 4-body interaction)
!!!##################
!!! shifter          (shifts the basis along the cell by an amount)
!!! shift_region     (shifts the basis within a region along the cell by amount)
!!! vacuumer         (adds a vacuum gap to the location specified)
!!! set_vacuum       (set the vacuum gap in a region to a specified value)
!!! vacater          (vacate a set of atoms from a basis)
!!! transformer      (applies a transformation matrix to a lattice and basis)
!!!##################
!!! MATNORM          (normalises a 3x3 matrix)
!!! ortho_axis       (makes specified axis perpendicular to plane of other two)
!!! change_basis     (convert basis into direct coords wrt another lattice)
!!! region_rot       (rotates a region specified along an axis about that axis)
!!! normalise_basis  (convert basis coordinates to be within val-> val-1)
!!!##################
!!! primitive_lat    (reorientates the lattice to the primitive lattice)
!!! reducer
!!! mkNiggli_lat
!!! reduced_check
!!! planecutter      (generates transformation mat to obtain miller plane)
!!! bas_merge        (merges two supplied bases)
!!! bas_lat_merge    (merges two supplied bases and lattices)
!!! split_bas        (split a basis into multiple sub-bases by region)
!!!##################
!!! get_bulk         (determines a bulk cell within a slab)
!!! get_wyckoff      (returns an array of the similar atoms)
!!!#############################################################################
module edit_geom
  use constants, only: pi,real12
  use rw_geom, only: bas_type,geom_write,convert_bas,clone_bas
  use misc, only: swap, sort1d
  use misc_linalg, only: cross,outer_product,cross_matrix,uvec,modu,&
       get_vol,det,inverse,inverse_3x3,LUinv,reduce_vec_gcd,get_vec_multiple,&
       proj,GramSchmidt,LLL_reduce, get_angle
  implicit none

  type wyck_atom_type
     integer, allocatable, dimension(:) :: atom
  end type wyck_atom_type
  type wyck_spec_type
     type(wyck_atom_type), allocatable, dimension(:) :: spec
  end type wyck_spec_type
  type bond_type
     real(real12) :: length
     integer, dimension(2,2) :: atoms
  end type bond_type

  !! %isa = integer list of spec and atom
  type neighbour_atom_type
     integer :: num
     integer, allocatable, dimension(:,:) :: isa !! set 2nd dimension to 2
     real(real12), allocatable, dimension(:) :: sep
  end type neighbour_atom_type

  type neighbour_spec_type
     type(neighbour_atom_type), allocatable, dimension(:) :: atom
  end type neighbour_spec_type

  
  interface get_closest_atom
     procedure get_closest_atom_1D,get_closest_atom_3D
  end interface get_closest_atom

  private

  public :: MATNORM
  public :: min_dist, get_atom_height, get_min_bulk_bond, get_min_bond
  public :: get_min_dist, get_shortest_bond
  public :: get_min_dist_between_point_and_atom, get_min_dist_between_two_atoms
  public :: shifter, shift_region, vacuumer, set_vacuum, ortho_axis
  public :: transformer, change_basis, region_rot, normalise_basis
  public :: centre_of_geom, centre_of_mass, primitive_lat, reducer
  public :: mkNiggli_lat, reduced_check, planecutter, bas_merge
  public :: bas_lat_merge, split_bas, get_bulk, get_centre_atom
  public :: get_wyckoff
  public :: get_closest_atom
  public :: wyck_atom_type, wyck_spec_type, bond_type

!!!updated 2024/03/20


contains
!!!#############################################################################
!!! Normalises a 3x3 matrix to the form:
!!! a 0 0
!!! b c 0
!!! d e f
!!! NEW NAME: lat_low
!!!#############################################################################
  pure function MATNORM(lat) result(nlat)
    implicit none
    real(real12), dimension(3,3), intent(in) :: lat
    real(real12), dimension(3,3) :: nlat
    nlat(1,1)=sqrt(lat(1,1)**2+lat(1,2)**2+lat(1,3)**2)
    nlat(1,2)=0.0
    nlat(1,3)=0.0

    nlat(2,1)=(lat(1,1)*lat(2,1)+lat(1,2)*lat(2,2)+lat(1,3)*lat(2,3))/nlat(1,1)
    nlat(2,2)=&
         sqrt((lat(1,2)*lat(2,3)-lat(1,3)*lat(2,2))**2+&
         (lat(1,3)*lat(2,1)-lat(1,1)*lat(2,3))**2+&
         (lat(1,1)*lat(2,2)-lat(1,2)*lat(2,1))**2)/nlat(1,1)
    nlat(2,3)=0.0

    nlat(3,1)=(lat(1,1)*lat(3,1)+lat(1,2)*lat(3,2)+lat(1,3)*lat(3,3))/nlat(1,1)
    nlat(3,2)=(&
         lat(2,1)*lat(3,1)+&
         lat(2,2)*lat(3,2)+&
         lat(2,3)*lat(3,3)-&
         nlat(2,1)*nlat(3,1))/nlat(2,2)
    nlat(3,3)=sqrt(&
         lat(3,1)**2+lat(3,2)**2+&
         lat(3,3)**2-nlat(3,1)**2-nlat(3,2)**2)
  end function MATNORM
!!!#############################################################################


!!!#############################################################################
!!! Finds distance between a location in a cell and ...
!!! ... the nearest atom to that point either above ...
!!! ... or below
!!!#############################################################################
  pure function min_dist(bas,axis,loc,above)
    implicit none
    integer, intent(in) :: axis
    type(bas_type), intent(in) :: bas
    logical, intent(in), optional :: above

    integer :: is
    real(real12) :: min_dist,pos
    real(real12), intent(in) :: loc
    logical :: labove


    pos=loc
    labove=.false.
    if(present(above)) labove=above
    aboveloop: if(labove)then
       min_dist=huge(0._real12)
       if(all( (/ (bas%spec(is)%atom(:,axis),is=1,bas%nspec) /).lt.pos))&
            pos=pos-1._real12
    else
       min_dist=-huge(0._real12)
       if(all( (/ (bas%spec(is)%atom(:,axis),is=1,bas%nspec) /).gt.pos))&
            pos=pos-1._real12
    end if aboveloop


    do is=1,bas%nspec
       if(.not.labove.and.maxval(bas%spec(is)%atom(:,axis)-pos,&
            mask=(bas%spec(is)%atom(:,axis)-pos.le.0._real12)).gt.min_dist) then
          min_dist=maxval(bas%spec(is)%atom(:,axis)-pos,&
               mask=(bas%spec(is)%atom(:,axis)-pos.le.0._real12))
       elseif(labove.and.minval(bas%spec(is)%atom(:,axis)-pos,&
            mask=(bas%spec(is)%atom(:,axis)-pos.ge.0._real12)).lt.min_dist) then
          min_dist=minval(bas%spec(is)%atom(:,axis)-pos,&
               mask=(bas%spec(is)%atom(:,axis)-pos.ge.0._real12))
       end if
    end do

  end function min_dist
!!!#############################################################################


!!!#############################################################################
!!! Iterative min distance
!!!#############################################################################
  pure function min_dist_iterative(lat,vec) result(dist)
    implicit none
    integer :: i,j,k
    real(real12) :: dtmp1,dist
    real(real12) :: a_mag2,b_mag2,c_mag2
    real(real12) :: ab_cosab,ac_cosac,bc_cosbc
    integer, dimension(3) :: expansion

    real(real12), dimension(3), intent(in) :: vec
    real(real12), dimension(3,3), intent(in) :: lat


   a_mag2 = dot_product(lat(1,:),lat(1,:))
   b_mag2 = dot_product(lat(2,:),lat(2,:))
   c_mag2 = dot_product(lat(3,:),lat(3,:))
   ab_cosab = dot_product(lat(1,:),lat(2,:))
   ac_cosac = dot_product(lat(1,:),lat(3,:))
   bc_cosbc = dot_product(lat(2,:),lat(3,:))

   dist = huge(1._real12)
   expansion = [0,0,0]
    do i=-10,10,1
       do j=-10,10,1
          do k=-10,10,1
             dtmp1 = sqrt(&
                  a_mag2 * (vec(1) + dble(i))**2._real12 + &
                  b_mag2 * (vec(2) + dble(j))**2._real12 + &
                  c_mag2 * (vec(3) + dble(k))**2._real12 + &
                  2._real12 * ab_cosab * (vec(1) + i) * (vec(2) + j) + &
                  2._real12 * ac_cosac * (vec(1) + i) * (vec(3) + k) + &
                  2._real12 * bc_cosbc * (vec(2) + j) * (vec(3) + k)&
                  )
             if(dtmp1.lt.dist)then
                dist = dtmp1
                expansion = [i,j,k]
             end if
          end do
       end do
    end do

  end function min_dist_iterative
!!!#############################################################################


!!!#############################################################################
!!! Return the surface normal vector
!!!#############################################################################
  function get_surface_normal(lat,axis) result(normal)
    implicit none
    real(real12) :: component
    integer, dimension(3) :: order=(/1,2,3/)
    real(real12), dimension(3) :: normal

    integer, intent(in) :: axis
    real(real12), dimension(3,3), intent(in) :: lat

    order = cshift(order,3-axis)
    normal = cross(lat(order(1),:),lat(order(2),:))
    component = dot_product(lat(3,:),normal) / modu(normal)**2._real12
    normal = normal * component

    return
  end function get_surface_normal
!!!#############################################################################


!!!#############################################################################
!!! Get the value of the atom along that axis
!!!#############################################################################
  function get_atom_height(bas,atom,axis) result(val)
    implicit none
    integer :: i,axis,atom,sum_atom
    real(real12) :: val
    type(bas_type) :: bas

    val=0._real12
    sum_atom=0
    do i=1,bas%nspec
       if(atom.le.sum_atom+bas%spec(i)%num)then
          val=bas%spec(i)%atom(atom-sum_atom,axis)
          return
       end if
       sum_atom=sum_atom+bas%spec(i)%num
    end do


  end function get_atom_height
!!!#############################################################################


!!!#############################################################################
!!! returns minimum bond within bulk
!!!#############################################################################
  function get_min_bulk_bond(lat,bas) result(min_bond)
    implicit none
    integer :: is,ia,js,ja
    real(real12) :: dtmp1,min_bond
    type(bas_type) :: bas
    real(real12), dimension(3) :: vdtmp1
    real(real12), dimension(3,3) :: lat


    min_bond=huge(0._real12)
    do is=1,bas%nspec
       do ia=1,bas%spec(is)%num

          do js=1,bas%nspec
             atmloop: do ja=1,bas%spec(js)%num
                if(is.eq.js.and.ia.eq.ja) cycle atmloop
                vdtmp1 = bas%spec(js)%atom(ja,:3) - bas%spec(is)%atom(ia,:3)
                vdtmp1 = &
                     vdtmp1(1)*lat(1,:3) + &
                     vdtmp1(2)*lat(2,:3) + &
                     vdtmp1(3)*lat(3,:3)
                dtmp1 = modu(vdtmp1)
                if(dtmp1.lt.min_bond) min_bond = dtmp1
             end do atmloop
          end do

       end do
    end do


  end function get_min_bulk_bond
!!!#############################################################################


!!!#############################################################################
!!! returns minimum bond for a specified atom
!!!#############################################################################
  function get_min_bond(lat,bas,is,ia,axis,labove,tol) result(vsave)
    implicit none
    integer :: js,ja
    integer :: axis_
    real(real12) :: dtmp1,min_bond,tol_
    logical :: labove_
    real(real12), dimension(3) :: vdtmp1, vsave

    integer, intent(in) :: is,ia
    type(bas_type), intent(in) :: bas
    real(real12), dimension(3,3), intent(in) :: lat

    integer, intent(in), optional :: axis
    real(real12), intent(in), optional :: tol
    logical, intent(in), optional :: labove

    if(present(tol))then
       tol_ = tol
    else
       tol_ = 1.E-5
    end if

    if(present(labove))then
       labove_=labove
    else
       labove_=.false.
    end if

    if(present(axis))then
       axis_=axis
    else
       axis_=0
    end if

    min_bond=huge(0._real12)
    
    do js=1,bas%nspec
       atmloop: do ja=1,bas%spec(js)%num
          if(is.eq.js.and.ia.eq.ja) cycle atmloop
          vdtmp1 = bas%spec(js)%atom(ja,:3) - bas%spec(is)%atom(ia,:3)
          if(axis_.gt.0)then
             if(abs(vdtmp1(axis_)).lt.tol_) cycle atmloop
             if(labove_)then
                vdtmp1(axis_) = 1._real12 + vdtmp1(axis_)
             else
                vdtmp1(axis_) = vdtmp1(axis_) - 1._real12
             end if
          end if
          vdtmp1 = &
               vdtmp1(1)*lat(1,:3) + &
               vdtmp1(2)*lat(2,:3) + &
               vdtmp1(3)*lat(3,:3)
          dtmp1 = modu(vdtmp1)
          if(dtmp1.lt.min_bond)then
             min_bond = dtmp1
             vsave = vdtmp1
          end if
       end do atmloop
    end do


  end function get_min_bond
!!!#############################################################################


!!!#############################################################################
!!! returns minimum bond for a specified atom
!!!#############################################################################
  function get_min_dist(lat,bas,loc,lignore_close,axis,labove,lreal,tol, &
       ignore_list) result(output)
    implicit none
    integer :: js,ja,i
    integer :: axis_ = 0
    real(real12) :: dtmp1,min_bond
    real(real12) :: tol_ = 1.E-5
    logical :: labove_ = .false.,lreal_ = .true.
    real(real12), dimension(3) :: vdtmp1,vdtmp2,output

    logical, intent(in) :: lignore_close
    type(bas_type), intent(in) :: bas
    real(real12), dimension(3), intent(in) :: loc
    real(real12), dimension(3,3), intent(in) :: lat

    integer, intent(in), optional :: axis
    real(real12), intent(in), optional :: tol
    logical, intent(in), optional :: labove, lreal
    integer, dimension(:,:), intent(in), optional :: ignore_list

    !! CORRECT tol TO ACCOUNT FOR LATTICE SIZE
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
          vdtmp2 = matmul(vdtmp1,lat)
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
!!!#############################################################################


!!!#############################################################################
!!! get the shortest distance between two atoms in a periodic cell
!!!#############################################################################
  pure function get_min_dist_between_two_atoms(lat,bas,atom_1,atom_2) &
       result(dist)
    implicit none
    type(bas_type), intent(in) :: bas
    integer, dimension(2), intent(in) :: atom_1,atom_2
    real(real12), dimension(3,3), intent(in) :: lat
    real(real12) :: dist
 
    real(real12), dimension(3) :: vec
 
    vec = bas%spec(atom_2(1))%atom(atom_2(2),:3) - &
          bas%spec(atom_1(1))%atom(atom_1(2),:3)
    vec = vec - ceiling(vec - 0.5_real12)
    vec = matmul(vec,lat)
    dist = modu(vec)
 
  end function get_min_dist_between_two_atoms
!!!#############################################################################


!!!#############################################################################
!!! get the shortest distance between a point and an atom in a periodic cell
!!!#############################################################################
  pure function get_min_dist_between_point_and_atom(lat,bas,loc,atom) &
       result(dist)
    implicit none
    type(bas_type), intent(in) :: bas
    integer, dimension(2), intent(in) :: atom
    real(real12), dimension(3), intent(in) :: loc
    real(real12), dimension(3,3), intent(in) :: lat
    real(real12) :: dist

    real(real12), dimension(3) :: vec

    vec = loc - bas%spec(atom(1))%atom(atom(2),:3)
    vec = vec - ceiling(vec - 0.5_real12)
    vec = matmul(vec,lat)
    dist = modu(vec)

  end function get_min_dist_between_point_and_atom
!!!#############################################################################


!!!#############################################################################
!!! identify the shortest bond in the crystal, takes in crystal basis
!!!#############################################################################
  function get_shortest_bond(lat,bas) result(bond)
    implicit none
    integer :: is,js,ia,ja,ja_start
    real(real12) :: dist,min_bond
    type(bas_type), intent(in) :: bas
    type(bond_type) :: bond
    real(real12), dimension(3) :: vec
    integer, dimension(2,2) :: atoms
    real(real12), dimension(3,3) :: lat
    
    min_bond = 100._real12
    atoms = 0
    do is=1,bas%nspec
       do js=is,bas%nspec
          do ia=1,bas%spec(is)%num
             if(is.eq.js)then
                ja_start = ia+1
             else
                ja_start = 1
             end if
             do ja=ja_start,bas%spec(js)%num
                vec = bas%spec(is)%atom(ia,:3) - bas%spec(js)%atom(ja,:3)
                vec = vec - ceiling(vec - 0.5_real12)
                vec = matmul(vec,lat)
                dist = modu(vec)
                if(dist.lt.min_bond)then
                   min_bond = dist
                   atoms(1,:) = (/is, ia/)
                   atoms(2,:) = (/js, ja/)
                end if
             end do
          end do
       end do
    end do
    bond%length = min_bond
    bond%atoms = atoms


  end function get_shortest_bond
!!!#############################################################################


!!!##########################################################################!!!
!!!##########################################################################!!!
!!! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  !!!
!!!##########################################################################!!!
!!!##########################################################################!!!


!!!#############################################################################
!!! get the neighbour fingerprint
!!!#############################################################################
  function get_neighbour_fingerprint(loc,lat,bas,cutoff,nbins,width) &
      result(fingerprint)
    implicit none
    real(real12), dimension(3), intent(in) :: loc
    real(real12), dimension(3,3), intent(in) :: lat
    type(bas_type), intent(in) :: bas

    integer, intent(in), optional :: nbins
    real(real12), intent(in), optional :: cutoff, width

    integer, allocatable, dimension(:,:) :: fingerprint

    integer :: i,j,k
    integer :: is,ia
    integer :: amax,bmax,cmax
    integer :: nbins_, bin
    real(real12) :: cutoff_,width_
    real(real12) :: dtmp1
    real(real12), dimension(3) :: vtmp1,diff


    if(present(cutoff))then
      cutoff_ = cutoff
    else
      cutoff_ = 6._real12
    end if
    if(present(width))then
      width_ = width
    else
      width_ = 0.1_real12
    end if
    if(present(nbins))then
      nbins_ = nbins
      width_ = cutoff_/nbins_
    else
      nbins_ = cutoff_/width_
    end if

    allocate(fingerprint(nbins_,bas%nspec))
    fingerprint = 0._real12

    !! this is not perfect
    !! won't work for extremely acute/obtuse angle cells
    !! (due to diagonal path being shorter than individual lattice vectors)
    amax = ceiling(cutoff_/modu(lat(1,:)))
    bmax = ceiling(cutoff_/modu(lat(2,:)))
    cmax = ceiling(cutoff_/modu(lat(3,:)))

    spec_loop: do is=1,bas%nspec
      atom_loop: do ia=1,bas%spec(is)%num
        diff = loc -  bas%spec(is)%atom(ia,:3)
        diff = diff - ceiling(diff - 0.5_real12)
        do i=-amax,amax+1,1
          vtmp1(1) = diff(1) + real(i)
          do j=-bmax,bmax+1,1
            vtmp1(2) = diff(2) + real(j)
            do k=-cmax,cmax+1,1
              vtmp1(3) = diff(3) + real(k)
              dtmp1 = modu(matmul(vtmp1,lat))
              if(dtmp1.lt.cutoff_ + width_/2._real12)then
                bin = nint(nbins_ * ( dtmp1 + width_/2._real12 ) / cutoff_)
                if(bin.gt.nbins_) cycle
                fingerprint(bin,is) = fingerprint(bin,is) + 1
              end if
            end do
          end do
        end do
      end do atom_loop
    end do spec_loop
   
 end function get_neighbour_fingerprint
!!!#############################################################################


!!!#############################################################################
!!! get the list of nearest neighbours
!!!#############################################################################
  function get_nearest_neighbours_atom(lat,bas,is,ia,max,cutoff) result(neighbour)
    implicit none
    integer :: i,j,k
    integer :: js,ja,dim,nneigh,natom,loc
    integer :: amax,bmax,cmax
    real(real12) :: dtmp1,tol
    real(real12), dimension(3) :: vtmp1,diff
    integer, allocatable, dimension(:,:) :: atom_list
    real(real12), allocatable, dimension(:) :: sep_list

    integer, intent(in) :: is,ia
    real(real12), dimension(3,3), intent(in) :: lat
    type(bas_type), intent(in) :: bas
    integer, intent(in), optional :: max
    real(real12), intent(in), optional :: cutoff

    type(neighbour_atom_type) :: neighbour

    if(present(cutoff))then
       tol = cutoff
    else
       tol = 6._real12
    end if
    !! define these based on tol val
    !! ... i.e. modu(lat(1,:)) compare to tol
    amax = 1
    bmax = 1
    cmax = 1
    nneigh = 0

    natom = bas%natom*(2*amax+1)*(2*bmax+1)*(2*cmax+1)
    allocate(atom_list(natom,2))
    allocate(sep_list(natom))

    spec_loop: do js=1,bas%nspec
       atom_loop: do ja=1,bas%spec(js)%num
          if(is.eq.js.and.ia.eq.ja) cycle atom_loop
          diff = bas%spec(is)%atom(ia,:3) -  bas%spec(js)%atom(ja,:3)
          diff = diff - ceiling(diff - 0.5_real12)
          do i=-amax,amax,1
             vtmp1(1) = diff(1) + real(i)
             do j=-bmax,bmax,1
                vtmp1(2) = diff(2) + real(j)
                do k=-cmax,cmax,1
                   vtmp1(3) = diff(3) + real(k)
                   dtmp1 = modu(matmul(vtmp1,lat))
                   if(dtmp1.le.tol)then
                      nneigh = nneigh + 1
                      atom_list(nneigh,:2) = [js,ja]
                      sep_list(nneigh) = dtmp1
                   end if
                end do
             end do
          end do
       end do atom_loop
    end do spec_loop
    
    dim = nneigh
    if(present(max))then
       if(max.gt.0) dim = min(nneigh,max)
    end if

    neighbour%num = dim
    allocate(neighbour%isa(dim,2))
    allocate(neighbour%sep(dim))
    do i=1,dim
       loc = minloc(sep_list(1:nneigh),dim=1)
       
       neighbour%sep(i) = sep_list(loc)
       neighbour%isa(i,:2) = atom_list(loc,:2)
       sep_list(loc) = huge(1._real12)
    end do

    deallocate(atom_list)
    deallocate(sep_list)


  end function get_nearest_neighbours_atom
!!!#############################################################################


!!!#############################################################################
!!! get the list of nearest neighbours
!!!#############################################################################
  function get_nearest_neighbours_basis(lat,bas,max,cutoff) result(neighbour_table)
    implicit none
    integer :: is,ia
    integer :: nmax
    real(real12) :: tol
    
    real(real12), dimension(3,3), intent(in) :: lat
    type(bas_type), intent(in) :: bas
    integer, intent(in), optional :: max
    real(real12), intent(in), optional :: cutoff

    type(neighbour_spec_type), allocatable, dimension(:) :: neighbour_table


    if(present(max))then
       nmax = max
    else
       nmax = -1
    end if

    if(present(cutoff))then
       tol = cutoff
    else
       tol = 6._real12
    end if


    allocate(neighbour_table(bas%nspec))

    do is=1,bas%nspec
       allocate(neighbour_table(is)%atom(bas%spec(is)%num))
       do ia=1,bas%spec(is)%num
          neighbour_table(is)%atom(ia) = &
               get_nearest_neighbours_atom(lat,bas,is,ia,nmax,tol)
       end do
    end do


  end function get_nearest_neighbours_basis
!!!#############################################################################


!!!##########################################################################!!!
!!!##########################################################################!!!
!!! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  !!!
!!!##########################################################################!!!
!!!##########################################################################!!!


!!!#############################################################################
!!! 2-body gvector radial descriptor
!!!#############################################################################
!!! Implementation follows original Behler and Parrinello paper
  pure function get_gvector_2body(lat, bas, species_1, species_2, &
       cutoff, nbins, width) &
       result(gvector)
    implicit none
    real(real12), dimension(3,3), intent(in) :: lat
    type(bas_type), intent(in) :: bas

    integer, intent(in), optional :: nbins
    real(real12), intent(in), optional :: cutoff, width
    character(*), intent(in), optional :: species_1, species_2

    real(real12), allocatable, dimension(:) :: gvector_tmp
    real(real12), allocatable, dimension(:,:) :: gvector

    integer :: i, j, k, b
    integer :: is, js, ia, ja
    integer :: amax, bmax, cmax
    integer :: species_1_num, species_2_num
    integer :: nbins_, bin
    real(real12) :: cutoff_, width_, fc, scale
    real(real12) :: rtmp1, rtmp2, eta
    real(real12), dimension(3) :: vtmp1, diff

    type :: bond_type
       integer, dimension(2) :: species
       real(real12) :: value
    end type bond_type

    type(bond_type), dimension(:), allocatable :: bond_list

    if(present(cutoff))then
       cutoff_ = cutoff
    else
       cutoff_ = 6._real12
    end if
    if(present(width))then
       width_ = width
    else
       width_ = 0.1_real12
    end if
    if(present(nbins))then
       nbins_ = nbins
       width_ = cutoff_/nbins_
    else
       nbins_ = cutoff_/width_
    end if

    eta = 1._real12 / ( 2._real12 * width_**2._real12 )

    if(present(species_1))then
       do i = 1, bas%nspec
          if(trim(bas%spec(i)%name).eq.trim(species_1))then
             species_1_num = i
             exit
          end if
       end do
       allocate(gvector(nbins_,species_1_num:species_1_num), source=0._real12)
    else
       species_1_num = 0
       allocate(gvector(nbins_,bas%nspec), source=0._real12)
    end if
    if(present(species_2))then
       do i = 1, bas%nspec
          if(trim(bas%spec(i)%name).eq.trim(species_2))then
             species_2_num = i
             exit
          end if
       end do
    else
       species_2_num = 0
    end if

    !! this is not perfect
    !! won't work for extremely acute/obtuse angle cells
    !! (due to diagonal path being shorter than individual lattice vectors)
    amax = ceiling(cutoff_/modu(lat(1,:)))
    bmax = ceiling(cutoff_/modu(lat(2,:)))
    cmax = ceiling(cutoff_/modu(lat(3,:)))

    allocate(bond_list(0)) !if doesn't work, allocate a dummy bond first
    spec_loop1: do is=1,bas%nspec
       if(present(species_1).and.species_1_num.ne.is) cycle spec_loop1
       atom_loop1: do ia=1,bas%spec(is)%num
          spec_loop2: do js=is,bas%nspec
             if(present(species_2).and.species_2_num.ne.js) cycle spec_loop2
             atom_loop2: do ja=1,bas%spec(is)%num
                if(is.eq.js.and.ja.lt.ia) cycle atom_loop2
                diff = bas%spec(is)%atom(ia,:3) -  bas%spec(js)%atom(ja,:3)
                diff = diff - ceiling(diff - 0.5_real12)
                do i=-amax,amax+1,1
                   vtmp1(1) = diff(1) + real(i)
                   do j=-bmax,bmax+1,1
                      vtmp1(2) = diff(2) + real(j)
                      do k=-cmax,cmax+1,1
                         vtmp1(3) = diff(3) + real(k)
                         rtmp1 = modu(matmul(vtmp1,lat))
                         if(rtmp1.lt.cutoff_ + width_/2._real12)then
                            bond_list = [ bond_list, bond_type([is,js],rtmp1) ]
                         end if
                      end do
                   end do
                end do
             end do atom_loop2
          end do spec_loop2
       end do atom_loop1
    end do spec_loop1

    allocate(gvector_tmp(nbins_))
    do i = 1, size(bond_list)
       if(abs(bond_list(i)%value).lt.1.E-3) cycle
       rtmp1 = bond_list(i)%value
       scale = 1._real12
       !! this is done before the lopp to help catch others just outside the cutoff
       bin = nint(nbins_ * ( rtmp1 + width_/2._real12 ) /cutoff_)
       do j = i+1, size(bond_list)
          !! don't need to look at reverse of species, as the ordering will ...
          !! always be enforced by the loop above, i.e. is <= js
          if( bond_list(i)%species(1).ne.bond_list(j)%species(1) .or. &
               bond_list(i)%species(2).ne.bond_list(j)%species(2) ) cycle
          if(abs(bond_list(j)%value-rtmp1).lt.1.E-3)then !!MAKE THIS LINKED TO WIDTH?
             bond_list(j)%value = 0._real12
             scale = scale + 1._real12
          end if          
       end do
       if(bin.gt.nbins_) cycle
   
       fc = 0.5_real12 * cos( pi * rtmp1 / cutoff_ ) + 0.5_real12
       fc = fc * scale
   
       gvector_tmp = 0._real12
       !! do forward loop to add gaussian for larger distances
       do b = bin, nbins_, 1
          rtmp2 = eta * ( rtmp1 - width_ * real(b) ) ** 2._real12
          if(rtmp2.gt.16._real12) cycle
          gvector_tmp(b) = gvector_tmp(b) + exp( -rtmp2 ) * fc
       end do
   
       !! do backward loop to add gaussian for smaller distances
       do b = bin-1, 1, -1
          rtmp2 = eta * ( rtmp1 - width_ * real(b) ) ** 2._real12
          if(rtmp2.gt.16._real12) cycle
          gvector_tmp(b) = gvector_tmp(b) + exp( -rtmp2 ) * fc
       end do
      
       gvector(:,bond_list(i)%species(1)) = &
            gvector(:,bond_list(i)%species(1)) + gvector_tmp
       gvector(:,bond_list(i)%species(2)) = &
            gvector(:,bond_list(i)%species(2)) + gvector_tmp
    end do

  end function get_gvector_2body
!!!-----------------------------------------------------------------------------
  pure function get_gvector_2body_alt(lat, bas, species_1, species_2, &
       cutoff, nbins, width) &
       result(gvector)
    implicit none
    real(real12), dimension(3,3), intent(in) :: lat
    type(bas_type), intent(in) :: bas

    integer, intent(in), optional :: nbins
    real(real12), intent(in), optional :: cutoff, width
    character(*), intent(in) :: species_1, species_2

    real(real12), allocatable, dimension(:) :: gvector

    integer :: i, j, k, b
    integer :: ia, ja
    integer :: amax, bmax, cmax
    integer :: species_1_num, species_2_num
    integer :: nbins_, bin
    real(real12) :: cutoff_, width_, fc, scale
    real(real12) :: rtmp1, rtmp2, eta
    real(real12), dimension(3) :: vtmp1, diff

    real(real12), dimension(:), allocatable :: bond_list

    if(present(cutoff))then
       cutoff_ = cutoff
    else
       cutoff_ = 6._real12
    end if
    if(present(width))then
       width_ = width
    else
       width_ = 0.1_real12
    end if
    if(present(nbins))then
       nbins_ = nbins
       width_ = cutoff_/nbins_
    else
       nbins_ = cutoff_/width_
    end if

    eta = 1._real12 / ( 2._real12 * width_**2._real12 )

    do i = 1, bas%nspec
       if(trim(bas%spec(i)%name).eq.trim(species_1))then
          species_1_num = i
          exit
       end if
    end do
    do i = 1, bas%nspec
       if(trim(bas%spec(i)%name).eq.trim(species_2))then
          species_2_num = i
          exit
       end if
    end do
    allocate(gvector(nbins_), source=0._real12)

    !! this is not perfect
    !! won't work for extremely acute/obtuse angle cells
    !! (due to diagonal path being shorter than individual lattice vectors)
    amax = ceiling(cutoff_/modu(lat(1,:)))
    bmax = ceiling(cutoff_/modu(lat(2,:)))
    cmax = ceiling(cutoff_/modu(lat(3,:)))

    allocate(bond_list(0)) !if doesn't work, allocate a dummy bond first
    atom_loop1: do ia=1,bas%spec(species_1_num)%num
       atom_loop2: do ja=1,bas%spec(species_2_num)%num
          if(species_1_num.eq.species_2_num.and.ja.lt.ia) cycle atom_loop2
          diff = bas%spec(species_1_num)%atom(ia,:3) -  &
               bas%spec(species_2_num)%atom(ja,:3)
          diff = diff - ceiling(diff - 0.5_real12)
          do i=-amax,amax+1,1
             vtmp1(1) = diff(1) + real(i)
             do j=-bmax,bmax+1,1
                vtmp1(2) = diff(2) + real(j)
                do k=-cmax,cmax+1,1
                   vtmp1(3) = diff(3) + real(k)
                   rtmp1 = modu(matmul(vtmp1,lat))
                   if(rtmp1.lt.cutoff_ + width_/2._real12)then
                      bond_list = [ bond_list, rtmp1 ]
                   end if
                end do
             end do
          end do
       end do atom_loop2
    end do atom_loop1

    do i = 1, size(bond_list)
       if(abs(bond_list(i)).lt.1.E-3) cycle
       rtmp1 = bond_list(i)
       scale = 1._real12
       !! this is done before the lopp to help catch others just outside the cutoff
       bin = nint(nbins_ * ( rtmp1 + width_/2._real12 ) /cutoff_)
       do j = i+1, size(bond_list)
          if(abs(bond_list(j)-rtmp1).lt.1.E-3)then !!MAKE THIS LINKED TO WIDTH?
             bond_list(j) = 0._real12
             scale = scale + 1._real12
          end if          
       end do
       if(bin.gt.nbins_) cycle
   
       fc = 0.5_real12 * cos( pi * rtmp1 / cutoff_ ) + 0.5_real12
       fc = fc * scale
   
       !! do forward loop to add gaussian for larger distances
       do b = bin, nbins_, 1
          rtmp2 = eta * ( rtmp1 - width_ * real(b) ) ** 2._real12
          if(rtmp2.gt.16._real12) cycle
          gvector(b) = gvector(b) + exp( -rtmp2 ) * fc
       end do
   
       !! do backward loop to add gaussian for smaller distances
       do b = bin-1, 1, -1
          rtmp2 = eta * ( rtmp1 - width_ * real(b) ) ** 2._real12
          if(rtmp2.gt.16._real12) cycle
          gvector(b) = gvector(b) + exp( -rtmp2 ) * fc
       end do
    end do

  end function get_gvector_2body_alt
!!!#############################################################################


!!!#############################################################################
!!! 3 body angular distribution function
!!!#############################################################################
!!! implemented as is done by Bircher et al. in their 2021 paper
!!! https://doi.org/10.1088/2632-2153/abf817
  pure function get_gvector_3body(lat, bas, species_1, x_min, x_max, &
       theta_min, theta_max, &
       nbins, width) &
       result(gvector)
    implicit none
    real(real12), dimension(3,3), intent(in) :: lat
    type(bas_type), intent(in) :: bas

    integer, intent(in), optional :: nbins
    real(real12), intent(in) :: x_min, x_max, theta_min, theta_max
    real(real12), intent(in), optional :: width
    character(*), intent(in) :: species_1

    real(real12), allocatable, dimension(:) :: gvector

    integer :: i, j, k, b
    integer :: is, js, ia, ja
    integer :: amax, bmax, cmax
    integer :: nbins_, bin
    integer :: species_1_num
    real(real12) :: cutoff, width_, fc
    real(real12) :: rtmp1, rtmp2, eta, x0, dx, dtheta, theta0, theta
    real(real12) :: fp1, fp2, fp3
    real(real12), dimension(3) :: vtmp1, vector, diff

    type :: bond_type
       real(real12), dimension(3) :: vector
    end type bond_type

    type(bond_type), dimension(:), allocatable :: bond_list


    if(present(width))then
       width_ = width
    else
       width_ = 0.1_real12
    end if
    if(present(nbins))then
       nbins_ = nbins
       width_ = (theta_max - theta_min)/nbins_
    else
       nbins_ = (theta_max - theta_min)/width_
     end if
     cutoff = x_max - x_min

    eta = 1._real12 / ( 2._real12 * width_**2._real12 )

    do i = 1, bas%nspec
       if(trim(bas%spec(i)%name).eq.trim(species_1))then
          species_1_num = i
          exit
       end if
    end do
    allocate(gvector(nbins_), source=0._real12)

    x0 = 0.5_real12 * ( x_max + x_min )
    dx = 0.5_real12 * ( x_max - x_min )
    !! default theta_min = 0
    !! default theta_max = pi
    theta0 = 0.5_real12 * ( theta_min + theta_max )
    dtheta = 0.5_real12 * ( theta_max - theta_min )

    !! this is not perfect
    !! won't work for extremely acute/obtuse angle cells
    !! (due to diagonal path being shorter than individual lattice vectors)
    amax = ceiling(cutoff/modu(lat(1,:)))
    bmax = ceiling(cutoff/modu(lat(2,:)))
    cmax = ceiling(cutoff/modu(lat(3,:)))

    atom_loop1: do ia=1,bas%spec(is)%num
       allocate(bond_list(0))
       spec_loop2: do js=is,bas%nspec
          atom_loop2: do ja=1,bas%spec(is)%num
             if(is.eq.js.and.ja.lt.ia) cycle atom_loop2
             diff = bas%spec(is)%atom(ia,:3) -  bas%spec(js)%atom(ja,:3)
             diff = diff - ceiling(diff - 0.5_real12)
             do i=-amax,amax+1,1
                vtmp1(1) = diff(1) + real(i)
                do j=-bmax,bmax+1,1
                   vtmp1(2) = diff(2) + real(j)
                   do k=-cmax,cmax+1,1
                      vtmp1(3) = diff(3) + real(k)
                      vector = matmul(vtmp1,lat)
                      rtmp1 = modu(vector)
                      if(rtmp1.lt.x_min.or.rtmp1.gt.x_max) cycle
                      bond_list = [ bond_list, bond_type(vector) ]
                   end do
                end do
             end do
          end do atom_loop2
       end do spec_loop2
       do i = 1, size(bond_list)
          fp1 = get_3body_cutoff_function( modu(bond_list(i)%vector), x0, dx )
          do j = i+1, size(bond_list), 1
             theta = get_angle( bond_list(i)%vector, bond_list(j)%vector)
   
             fp2 = get_3body_cutoff_function( modu(bond_list(j)%vector), x0, dx )
             fp3 = get_3body_cutoff_function(theta, theta0, dtheta)
             

             fc = fp1 * fp2 * fp3
   
             if(fc.lt.1.E-3) cycle

             !! do forward loop to add gaussian for larger distances
             do b = bin, nbins_, 1
                rtmp2 = eta * ( rtmp1 - width_ * real(b) ) ** 2._real12
                if(rtmp2.gt.16._real12) cycle
                gvector(b) = gvector(b) + exp( -rtmp2 ) * fc
             end do

             !! do backward loop to add gaussian for smaller distances
             do b = bin-1, 1, -1
                rtmp2 = eta * ( rtmp1 - width_ * real(b) ) ** 2._real12
                if(rtmp2.gt.16._real12) cycle
                gvector(b) = gvector(b) + exp( -rtmp2 ) * fc
             end do

          end do
       end do
       deallocate(bond_list)
    end do atom_loop1

  end function get_gvector_3body
!!!#############################################################################
  pure function get_3body_cutoff_function(x, x0, dx) result(output)
    implicit none
    real(real12), intent(in) :: x, x0, dx
    real(real12) :: output
    
    if(x.ge.x0.and.x.le.x0+dx)then
       output = get_fp2_function( ( x - x0 ) / ( 2._real12 * dx ) )
    elseif(x.ge.x0-dx.and.x.le.x0)then
       output = get_fp2_function( ( x0 - x ) / ( 2._real12 * dx ) )
    else
       output = 0._real12
    end if

  end function get_3body_cutoff_function
  pure function get_fp2_function(x) result(output)
    implicit none
    real(real12), intent(in) :: x
    real(real12) :: output

    output = x ** 3._real12 * &
             ( x * ( 15._real12 - 6._real12 * x ) - 10._real12 ) + 1._real12

  end function get_fp2_function




!!!##########################################################################!!!
!!!##########################################################################!!!
!!! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  !!!
!!!##########################################################################!!!
!!!##########################################################################!!!


!!!#############################################################################
!!! Shifts the basis along a, b or c by amount 'shift'
!!!#############################################################################
  subroutine shifter(bas,axis,shift,ltmp)
    implicit none
    integer :: i,j,k,axis
    real(real12) :: shift
    type(bas_type) :: bas
    logical, optional ::ltmp
    logical :: lrenorm

    k=axis
    lrenorm=.false.
    if(present(ltmp)) lrenorm=ltmp

    do i=1,bas%nspec
       do j=1,bas%spec(i)%num
          bas%spec(i)%atom(j,k)=bas%spec(i)%atom(j,k) + shift
          if(lrenorm) bas%spec(i)%atom(j,k)=bas%spec(i)%atom(j,k) - &
               floor(bas%spec(i)%atom(j,k))
       end do
    end do

  end subroutine shifter
!!!#############################################################################


!!!#############################################################################
!!! Shifts basis in a region by an amount
!!!#############################################################################
  subroutine shift_region(bas,region_axis,region_lw,region_up,shift_axis,shift,renorm)
    implicit none
    integer :: is,ia,shift_axis,region_axis
    real(real12) :: shift,region_lw,region_up
    type(bas_type) :: bas
    logical, optional ::renorm
    logical :: lrenorm

    lrenorm=.false.
    if(present(renorm)) lrenorm=renorm

    do is=1,bas%nspec
       do ia=1,bas%spec(is)%num
          if(bas%spec(is)%atom(ia,region_axis).ge.region_lw.and.&
               bas%spec(is)%atom(ia,region_axis).le.region_up)then
             bas%spec(is)%atom(ia,shift_axis) = &
                  bas%spec(is)%atom(ia,shift_axis) + shift
             if(lrenorm) bas%spec(is)%atom(ia,shift_axis) = &
                  bas%spec(is)%atom(ia,shift_axis) - &
                  floor(bas%spec(is)%atom(ia,shift_axis))
          end if
       end do
    end do

  end subroutine shift_region
!!!#############################################################################


!!!#############################################################################
!!! Adjusts the amount of vacuum at a location ...
!!! ... within a cell and adjusts the basis accordingly
!!!#############################################################################
  subroutine vacuumer(lat,bas,axis,loc,add,otol)
    implicit none
    integer :: is,ia
    integer, intent(in) :: axis
    real(real12) :: tol,tloc
    real(real12) :: cur_vac,inc,diff,mag_old,mag_new
    real(real12), intent(in) :: add,loc
    type(bas_type) :: bas
    real(real12), optional :: otol
    real(real12),dimension(3,3) :: lat


    tol=1.E-5
    inc=add
    if(present(otol)) tol=otol
    cur_vac=min_dist(bas,axis,loc,.true.)-min_dist(bas,axis,loc,.false.)
    cur_vac=cur_vac*modu(lat(axis,:))
    diff=cur_vac+inc
    if(diff.lt.0._real12)then
       write(0,*) "WARNING! Removing vacuum entirely"
    end if

    mag_old=modu(lat(axis,:))
    mag_new=(mag_old+inc)/mag_old
    lat(axis,:)=lat(axis,:)*mag_new
    inc=inc/modu(lat(axis,:))
    tol=tol/mag_old
    tloc=loc/mag_new+tol


    do is=1,bas%nspec
       do ia=1,bas%spec(is)%num
          bas%spec(is)%atom(ia,axis)=bas%spec(is)%atom(ia,axis)/mag_new
          if(bas%spec(is)%atom(ia,axis).gt.tloc) then
             bas%spec(is)%atom(ia,axis)=bas%spec(is)%atom(ia,axis)+inc
          end if
       end do
    end do


  end subroutine vacuumer
!!!#############################################################################


!!!#############################################################################
!!! Adjusts the amount of vacuum at a location ...
!!! ... within a cell and adjusts the basis accordingly
!!!#############################################################################
  subroutine set_vacuum(lat,bas,axis,loc,vac,otol)
    implicit none
    integer :: is,ia
    integer, intent(in) :: axis
    real(real12) :: tol,tloc
    real(real12) :: cur_vac,diff,mag_old,mag_new
    real(real12), intent(in) :: vac,loc
    type(bas_type) :: bas
    real(real12), optional :: otol
    real(real12),dimension(3,3) :: lat


    tol=0._real12
    if(present(otol)) tol=otol
    if(vac.lt.0._real12)then
       write(0,*) "WARNING! Removing vacuum entirely"
    end if
    cur_vac=min_dist(bas,axis,loc,.true.)-min_dist(bas,axis,loc,.false.)
    cur_vac=cur_vac*modu(lat(axis,:))
    diff=vac-cur_vac

    mag_old=modu(lat(axis,:))
    mag_new=(mag_old+diff)/mag_old
    lat(axis,:)=lat(axis,:)*mag_new
    diff=diff/modu(lat(axis,:))
    tol=tol/mag_old
    tloc=loc/mag_new+tol



    do is=1,bas%nspec
       do ia=1,bas%spec(is)%num
          bas%spec(is)%atom(ia,axis)=bas%spec(is)%atom(ia,axis)/mag_new
          if(bas%spec(is)%atom(ia,axis).gt.tloc) then
             bas%spec(is)%atom(ia,axis)=bas%spec(is)%atom(ia,axis)+diff
          end if
       end do
    end do


  end subroutine set_vacuum
!!!#############################################################################


!!!#############################################################################
!!! Takes a lattice and makes the define axis orthogonal to the other two
!!! WARNING! THIS IS FOR SLAB STRUCTURES! IT REMOVES PERIODICITY ALONG THAT AXIS
!!!#############################################################################
  subroutine ortho_axis(lat,bas,axis)
    implicit none
    integer :: axis
    real(real12) :: ortho_comp
    type(bas_type) :: bas
    integer, dimension(3) :: order
    real(real12), dimension(3) :: ortho_vec
    real(real12), dimension(3,3) :: invlat,lat


    bas=convert_bas(bas,transpose(lat))
    order=(/1,2,3/)
    order=cshift(order,3-axis)

    ortho_vec=cross(lat(order(1),:),lat(order(2),:))
    ortho_comp=dot_product(lat(3,:),ortho_vec)/modu(ortho_vec)**2._real12
    ortho_vec=ortho_vec*ortho_comp

    lat(3,:)=ortho_vec
    invlat=inverse_3x3(lat)
    bas=convert_bas(bas,transpose(invlat))


    return
  end subroutine ortho_axis
!!!#############################################################################


!!!#############################################################################
!!! Applies a transformation matrix to a lattice ...
!!! ... and extends the basis where needed
!!!#############################################################################
  subroutine transformer(lat,bas,tfmat,map)
    implicit none
    integer :: i,j,k,l,m,n,is,ia
    integer :: satom,dim
    real(real12) :: tol,vol_inc
    logical :: lmap
    type(bas_type) :: bas,sbas
    integer, dimension(3) :: latmin,latmax
    real(real12), dimension(3):: translvec,tolvec
    integer, allocatable, dimension(:) :: tmp_map_atom
    integer, allocatable, dimension(:,:,:) :: new_map
    real(real12), allocatable, dimension(:,:) :: tmpbas
    real(real12), dimension(3,3) :: lat,slat,tfmat,invmat

    integer, allocatable, dimension(:,:,:), optional, intent(inout) :: map

    vol_inc = abs(det(lat))
    if(vol_inc.lt.0.5_real12)then
       write(0,'(1X,"ERROR: Internal error in transformer function")')
       write(0,'(2X,"transformer in mod_edit_geom.f90 been supplied a&
            & lattice with almost zero determinant")')
       write(0,'(2X,"determinant = ",F0.9)') vol_inc
       write(0,'(3(1X,F7.2))') lat
       stop
    end if
    call normalise_basis(bas,1._real12,lfloor=.true.,lround=.false.)
    vol_inc=abs(det(tfmat))
    slat=matmul(tfmat,lat)
    invmat=inverse_3x3(tfmat)
    translvec=0._real12
    dim=size(bas%spec(1)%atom(1,:))
    
    
    !!--------------------------------------------------------------------------
    !! If map is present, sets up new map
    !!--------------------------------------------------------------------------
    lmap = .false.
    if_map: if(present(map))then
       if(all(map.eq.-1))then
          exit if_map
       end if
       lmap = .true.
       allocate(new_map(&
            bas%nspec,&
            ceiling(vol_inc)*maxval(bas%spec(:)%num,dim=1),2))
       new_map=0
       if(all(map.eq.0))then
          do is=1,bas%nspec
             map(is,:bas%spec(is)%num,1) = is
             do ia=1,bas%spec(is)%num
                map(is,ia,2) = ia
             end do
          end do
       end if
    end if if_map
    
    
    !!--------------------------------------------------------------------------
    !! Convert tolerance from Å to a fraction of each direction
    !!--------------------------------------------------------------------------
    tol=1.E-3 !! in Å
    do i=1,3
       tolvec(i)=tol/modu(slat(i,:))
    end do
    if(vol_inc.lt.minval(tolvec))then
       write(0,'(1X,"ERROR: Internal error in transformer function")')
       write(0,'(2X,"transformer in mod_edit_geom.f90 been supplied a&
            & transformation matrix with almost zero determinant")')
       write(0,'(2X,"determinant = ",F0.9)') vol_inc
       write(0,'(3(1X,F7.2))') tfmat
       stop
    end if
    
    
    !!--------------------------------------------------------------------------
    !! extends basis from min to sum(:) for each column
    !!--------------------------------------------------------------------------
    !! latmin is how far down you have to move to get it back to 0.
    !! hence, latmin=-nint(max*vol)
    !! latmax is how far you still have left to go.
    !! hence, latmax=vol-nint(max*vol)
    !! 
    !! This is why the sum works.
    !! Distance above origin using a is the sum of all positive a's in tfmat
    !! Distance below origin using a is the sum of all negative a's in tfmat
    !!
    !!     ____b
    !!    /   /
    !!  a/   /
    !!  /___/
    !! o 
    !!   
    !!    /\
    !!  a/  \b
    !!  /   /
    !! o\  /
    !!   \/
    !!----------------------------------
    !latmin(i)=(minval(invmat(i,:))-ceiling(minval(invmat(i,:))))*vol
    !latmax(i)=(maxval(invmat(i,:))-floor(minval(invmat(i,:))))*vol
    !latmin(i)=(min(minval(invmat(i,:)),0._real12)-ceiling(minval(invmat(i,:))))*vol
    !latmax(i)=(ceiling(maxval(invmat(i,:3)))-maxval(invmat(i,:3)) )*vol
    do i=1,3
       latmin(i)=floor(sum(tfmat(:3,i),mask=tfmat(:3,i).lt.0._real12))-1
       latmax(i)=ceiling(sum(tfmat(:3,i),mask=tfmat(:3,i).gt.0._real12))+1
    end do
    
    
    !!--------------------------------------------------------------------------
    !! transform the basis
    !!--------------------------------------------------------------------------
    do i=1,bas%nspec
       do j=1,bas%spec(i)%num
          bas%spec(i)%atom(j,:3)=matmul(bas%spec(i)%atom(j,:3),invmat)
       end do
    end do
    
    
    !!--------------------------------------------------------------------------
    !! generates atoms to fill the supercell
    !!--------------------------------------------------------------------------
    allocate(sbas%spec(bas%nspec))
    sbas%sysname=bas%sysname
    sbas%nspec=0
    sbas%natom=0
    spec_loop1: do is=1,bas%nspec
       if(allocated(tmpbas)) deallocate(tmpbas)
       allocate(tmpbas(bas%spec(is)%num*(&
            (abs(latmax(3))+abs(latmin(3))+1)*&
            (abs(latmax(2))+abs(latmin(2))+1)*&
            (abs(latmax(1))+abs(latmin(1))+1)),3))
       satom=0
       if(lmap)then
          allocate(tmp_map_atom(ceiling(vol_inc)*bas%spec(is)%num))
       end if
       do ia=1,bas%spec(is)%num
          do n=latmin(3),latmax(3)!,1
             translvec(3)=dble(n)
             do m=latmin(2),latmax(2)!,1
                translvec(2)=dble(m)
                inloop: do l=latmin(1),latmax(1)!,1
                   translvec(1)=dble(l)
                   tmpbas(satom+1,:3) = &
                        bas%spec(is)%atom(ia,:3) + matmul(translvec,invmat)
                   !!tmpbas(satom+1,:3)=&
                   !!     matmul((bas%spec(is)%atom(ia,:3)+translvec),invmat)
                   !where(abs(tmpbas(satom+1,:3)-nint(tmpbas(satom+1,k))).lt.tol)
                   !   tmpbas(satom+1,:3)=nint(tmpbas(satom+1,:3))
                   !end where
                   !if(any(tmpbas(satom+1,:).ge.1._real12).or.&
                   !     any(tmpbas(satom+1,:).lt.0._real12)) cycle
                   !if(any(tmpbas(satom+1,:).ge.1._real12+tol).or.&
                   !     any(tmpbas(satom+1,:).lt.0._real12-tol)) cycle
                   if(any(tmpbas(satom+1,:).ge.1._real12-tol).or.&
                        any(tmpbas(satom+1,:).lt.0._real12-tol)) cycle inloop !??? cycle inloop or spec_loop1?
                   tmpbas(satom+1,:3) = tmpbas(satom+1,:3) - &
                        dble(floor(tmpbas(satom+1,:3)))
                   do k=1,satom
                      if(all(mod(abs(tmpbas(satom+1,:3)-tmpbas(k,:3)),1._real12).le.&
                           tol)) cycle inloop
                   end do
                   if(lmap) tmp_map_atom(satom+1)=map(is,ia,2)
                   satom=satom+1
                end do inloop
             end do
          end do
       end do
       if(satom.eq.0)then
          if(lmap) deallocate(tmp_map_atom)
          cycle spec_loop1
       end if
       sbas%nspec=sbas%nspec+1
       sbas%spec(sbas%nspec)%num=satom
       sbas%natom=sbas%natom+satom
       sbas%spec(sbas%nspec)%name=bas%spec(is)%name
       allocate(sbas%spec(sbas%nspec)%atom(satom,dim))
       sbas%spec(sbas%nspec)%atom(1:satom,:3)=tmpbas(1:satom,:3)
       if(dim.eq.4) sbas%spec(sbas%nspec)%atom(1:satom,4)=1._real12
       deallocate(tmpbas)
       deallocate(bas%spec(is)%atom)
       if(lmap)then
          new_map(sbas%nspec,:satom,1) = is
          new_map(sbas%nspec,:satom,2) = tmp_map_atom(:satom)
          deallocate(tmp_map_atom)
       end if
    end do spec_loop1
    
    
    !!--------------------------------------------------------------------------
    !! check to see if successfully generated correct number of atoms
    !!--------------------------------------------------------------------------
    if(all(abs(tfmat-nint(tfmat)).lt.tol))then
       if(nint(bas%natom*vol_inc).ne.sbas%natom)then
          write(0,'(1X,"ERROR: Internal error in transformer function")')
          write(0,'(2X,"Transformer in mod_edit_geom.f90 has failed to &
               &generate enough atoms when extending the cell")')
          write(0,'(2X,"Generated ",I0," atoms, whilst expecting ",I0," atoms")') &
               sbas%natom,nint(bas%natom*vol_inc)
          write(0,*) bas%natom,nint(vol_inc)
          write(0,'(3(1X,F7.2))') tfmat
          open(60,file="broken_cell.vasp")
          call geom_write(60,slat,sbas)
          close(60)
          stop
       end if
    end if
    
    
    !!--------------------------------------------------------------------------
    !! saves new lattice and basis to original set
    !!--------------------------------------------------------------------------
    lat=slat
    deallocate(bas%spec)
    allocate(bas%spec(sbas%nspec))
    bas%sysname=sbas%sysname
    bas%nspec=sbas%nspec
    bas%natom=sbas%natom
    do i=1,sbas%nspec
       allocate(bas%spec(i)%atom(sbas%spec(i)%num,dim))
       bas%spec(i)=sbas%spec(i)
    end do
    
    
    !!--------------------------------------------------------------------------
    !! sets up the new map, if map supplied
    !!--------------------------------------------------------------------------
    if(lmap)then
       deallocate(map)
       call move_alloc(new_map,map)
    end if



  end subroutine transformer
!!!#############################################################################


!!!#############################################################################
!!! Convert basis from direct coords in one lattice ...
!!! ... into direct coords wrt another lattice
!!!#############################################################################
  function change_basis(vec,old_lat,new_lat)
    implicit none
    real(real12), dimension(3) :: change_basis,vec
    real(real12), dimension(3,3), intent(in) :: old_lat,new_lat
    real(real12), dimension(3,3) :: inew_lat
    inew_lat=inverse_3x3(new_lat)
    change_basis=matmul(transpose(inew_lat),matmul(old_lat,vec))
  end function change_basis
!!!#############################################################################


!!!#############################################################################
!!! rotates a region along an axis about that axis
!!!#############################################################################
  subroutine region_rot(bas,lat,angle,axis,bound1,bound2,tvec)
    implicit none
    integer :: axis,i,j
    real(real12) :: angle,bound1,bound2
    real(real12), dimension(3) :: u,centre
    real(real12), dimension(3,3) :: rotmat,ident,lat,invlat
    type(bas_type) :: bas
    real(real12), optional, dimension(3) :: tvec

    centre=(/0.5,0.5,0.0/)
    if(present(tvec)) centre=tvec
    ident=0._real12
    do i=1,3
       ident(i,i)=1._real12
    end do

!!! DEFINE ROTMAT BEFORE THIS
    u=0._real12
    u(axis)=-1._real12
    rotmat=&
         (cos(angle)*ident)+&
         (sin(angle))*cross_matrix(u)+&
         (1-cos(angle))*outer_product(u,u)


!!! Transform the rotation matrix into direct space
    invlat=LUinv(lat)
    rotmat=matmul(lat,rotmat)
    rotmat=matmul(rotmat,invlat)


!!! Rotate the basis within the bounds
    do i=1,bas%nspec
       do j=1,bas%spec(i)%num
          if(bas%spec(i)%atom(j,axis).lt.bound1.or.&
               bas%spec(i)%atom(j,axis).gt.bound2) cycle
          bas%spec(i)%atom(j,:3)=&
               matmul(rotmat,bas%spec(i)%atom(j,:3)-centre)+centre
       end do
    end do


    return
  end subroutine region_rot
!!!#############################################################################


!!!#############################################################################
!!! convert basis coordinates to be within +val -> val-1
!!!#############################################################################
!!! MAKE AN OPTIONAL SETTING TO MOVE CLOSES ATOM TO ZERO TO ACTUAL ZERO
  subroutine normalise_basis(bas,dtmp,lfloor,lround)
    implicit none
    integer :: is,ia,j
    real(real12) :: ceil,flr,dround
    real(real12), optional :: dtmp
    type(bas_type) :: bas
    logical :: lfloor1,lround1
    logical, optional :: lfloor,lround


    ceil=1._real12
    lfloor1=.false.
    if(present(dtmp)) ceil=dtmp
    if(present(lfloor)) lfloor1=lfloor
    flr=ceil-1._real12
    lround1=.false.
    dround=1.E-8
    if(present(lround)) lround1=lround

    do is=1,bas%nspec
       do ia=1,bas%spec(is)%num
          do j=1,3
             if(lfloor1)then
                bas%spec(is)%atom(ia,j)=bas%spec(is)%atom(ia,j)&
                     -floor(bas%spec(is)%atom(ia,j)-flr)
             else
                bas%spec(is)%atom(ia,j)=bas%spec(is)%atom(ia,j)&
                     -ceiling(bas%spec(is)%atom(ia,j)-ceil)
             end if
             if(lround1)then
                if(abs(bas%spec(is)%atom(ia,j)-ceil).lt.dround.or.&
                     abs(bas%spec(is)%atom(ia,j)).lt.dround) &
                     bas%spec(is)%atom(ia,j)=flr
             end if
          end do
       end do
    end do


    return
  end subroutine normalise_basis
!!!#############################################################################


!!!#############################################################################
!!! finds the centre of geometry of the supplied basis
!!!#############################################################################
  function centre_of_geom(bas) result(centre)
    implicit none
    integer :: is,ia,j
    real(real12), dimension(3) :: centre
    type(bas_type) :: bas

    centre=0._real12
    do is=1,bas%nspec
       do ia=1,bas%spec(is)%num
          do j=1,3
             centre(j)=centre(j)+bas%spec(is)%atom(ia,j)
          end do
       end do
    end do

    centre=centre/bas%natom

    return
  end function centre_of_geom
!!!#############################################################################


!!!#############################################################################
!!! finds the centre of mass of the supplied basis
!!!#############################################################################
  function centre_of_mass(bas) result(centre)
    implicit none
    integer :: is,ia,j
    real(real12) :: tot_mass
    real(real12), dimension(3) :: centre
    type(bas_type) :: bas

    centre=0._real12
    tot_mass=0._real12
    do is=1,bas%nspec
       tot_mass=tot_mass+bas%spec(is)%mass*bas%spec(is)%num
       do ia=1,bas%spec(is)%num
          do j=1,3
             centre(j)=centre(j)+bas%spec(is)%atom(ia,j)*bas%spec(is)%mass
          end do
       end do
    end do

    centre=centre/tot_mass

    return
  end function centre_of_mass
!!!#############################################################################


!!!#############################################################################
!!! Reorientates lattice to the primitive lattice of its type
!!!#############################################################################
!!! NEED TO SET UP TO WORK FOR THE EXTRA SWAPPINGS OF A, B AND C
  function primitive_lat(inlat) result(plat)
    implicit none
    integer :: i,j
    real(real12) :: dtmp1
    real(real12), dimension(3) :: scal
    real(real12), dimension(3,3) :: lat,plat,tmat1,tmat2
    real(real12), dimension(3,3), intent(in) :: inlat
    real(real12), dimension(6,3,3) :: special


    !!---------------------------------------------------------------
    !! makes all lattice vectors unity
    !!---------------------------------------------------------------
    lat=inlat
    plat=lat
    do i=1,3
       scal(i)=modu(lat(i,:))
       lat(i,:)=lat(i,:)/scal(i)
    end do


    !!---------------------------------------------------------------
    !! sets up the special set of primitive lattices
    !!---------------------------------------------------------------
    special(1,:,:) = transpose( reshape( (/&
         1._real12, 0._real12, 0._real12,&
         0._real12, 1._real12, 0._real12,&
         0._real12, 0._real12, 1._real12/), shape(lat) ) )
    special(2,:,:) = transpose( reshape( (/&
         1._real12, 0._real12, 0._real12,&
         -0.5_real12, sqrt(3._real12)/2._real12, 0._real12,&
         0._real12, 0._real12, 1._real12/), shape(lat) ) )
    special(3,:,:) = transpose( reshape( (/&
         0._real12, 1._real12, 1._real12,&
         1._real12, 0._real12, 1._real12,&
         1._real12, 1._real12, 0._real12/), shape(lat) ) )
    special(3,:,:) = special(3,:,:)/sqrt(2._real12)
    special(4,:,:) = transpose( reshape( (/&
         -1._real12,  1._real12,  1._real12,&
         1._real12, -1._real12,  1._real12,&
         1._real12,  1._real12, -1._real12/), shape(lat) ) )
    special(4,:,:) = special(4,:,:)/sqrt(3._real12)
    !! the following are obtuse versions of the acute
    special(5,:,:) = transpose( reshape( (/&
         1._real12, 0._real12, 0._real12,&
         0.5_real12, sqrt(3._real12)/2._real12, 0._real12,&
         0._real12, 0._real12, 1._real12/), shape(lat) ) )
    special(6,:,:) = transpose( reshape( (/&
         1._real12,    0.25_real12,  0.25_real12,&
         0.25_real12,  1._real12,    0.25_real12,&
         0.25_real12,  0.25_real12,  1._real12/), shape(lat) ) )
    special(6,:,:) = special(6,:,:)/sqrt(3._real12)


    !!---------------------------------------------------------------
    !! cycles special set to find primitive lattice of supplied lat
    !!---------------------------------------------------------------
    tmat1=matmul(lat,transpose(lat))
    checkloop: do i=1,6
       !tfmat=matmul(lat,inverse_3x3(special(i,:,:)))
       !tfmat=matmul(tfmat,transpose(tfmat))
       tmat2=matmul(special(i,:,:),transpose(special(i,:,:)))
       dtmp1=tmat2(1,1)/tmat1(1,1)
       !if(all(abs(tfmat-nint(tfmat)).lt.1.E-8))then
       if(all(abs(tmat1*dtmp1-tmat2).lt.real(1.E-6,real12)))then
          do j=1,3
             plat(j,:)=scal(j)*special(i,j,:)
          end do
          exit checkloop
       end if
    end do checkloop


  end function primitive_lat
!!!#############################################################################


!!!#############################################################################
!!! Uses Buerger's algorithm to reduce cell.
!!!#############################################################################
  subroutine reducer(lat,bas,tmptype,ltmp)
    implicit none
    integer :: cell_type
    integer :: i,j,k,count,limit
    real(real12), dimension(3,3) :: lat,newlat,transmat,S,tmp_mat
    real(real12) :: tiny,pi,pi2
    logical :: verb,lreduced
    integer, optional :: tmptype
    logical, optional :: ltmp
    type(bas_type) :: bas



!!!-----------------------------------------------------------------------------
!!! set up inital variable values
!!!-----------------------------------------------------------------------------
    verb=.false.
    if(present(ltmp)) verb=ltmp
    cell_type=2
    if(present(tmptype)) cell_type=tmptype
    S=0._real12
    count=0
    limit=100
    lreduced=.false.
    tiny=1E-5*(get_vol(lat))**(1.E0/3.E0)
    pi=4._real12*atan(1._real12)
    pi2=2._real12*atan(1._real12)
    transmat=0._real12
    do i=1,3
       transmat(i,i)=1._real12
    end do
    newlat=lat


!!!-----------------------------------------------------------------------------
!!! performs checks on the other main conditions defined by Niggli
!!!-----------------------------------------------------------------------------
    find_reduced: do while(.not.lreduced)
       count=count+1
       call mkNiggli_lat(lat,newlat,transmat,S)
       lreduced=reduced_check(newlat,cell_type,S)
       if(lreduced) exit
       if(verb) then
          write(67,*)
          write(67,*) count
          write(67,*) "###############"
          write(67,*) (transmat(i,:),i=1,3)
          write(67,*)
          write(67,*) (newlat(i,:),i=1,3)
       end if
       if(count.gt.limit) then
          write(0,'("FAILED to find the reduced cell within ",I0," steps")') count
          exit
       end if


       !! A1 & A2 
       do i=1,2
          j=i+1
          if(S(i,i)-S(j,j).gt.tiny) then
             call swap(transmat(i,:),transmat(j,:))
             transmat=-transmat
             if(i.eq.2) cycle find_reduced
             call mkNiggli_lat(lat,newlat,transmat,S)
          end if
       end do


       !! A3
       i=1;j=1;k=1
       if(S(2,3).lt.0) i=-1
       if(S(1,3).lt.0) j=-1
       if(S(1,2).lt.0) k=-1
       if(i*j*k.gt.0) then
          tmp_mat=reshape((/i,0,0,  0,j,0,  0,0,k/),shape(tmp_mat))
          transmat=matmul(transpose(tmp_mat),transmat)
          call mkNiggli_lat(lat,newlat,transmat,S)
       end if


       !! A4
       i=1;j=1;k=1
       if(S(1,2).lt.0) i=-1
       if(S(1,3).lt.0) j=-1
       if(S(1,2).lt.0) k=-1
       if(i*j*k.gt.0) then
          tmp_mat=reshape((/i,0,0,  0,j,0,  0,0,k/),shape(tmp_mat))
          transmat=matmul(transpose(tmp_mat),transmat)
          call mkNiggli_lat(lat,newlat,transmat,S)
       end if


       tmp_mat=reshape((/1,0,0,  0,1,0,  0,0,1/),shape(tmp_mat))
       !! A5
       if(abs(2*S(2,3)).gt.S(2,2)+tiny.or.&
            (abs(2*S(2,3)-S(2,2)).le.tiny.and.2*S(1,3).lt.S(1,2)).or.&
            (abs(2*S(2,3)+S(2,2)).le.tiny.and.S(1,2).lt.0._real12))then
          tmp_mat(2,3)=((-1)**(cell_type+1))*floor((2*S(2,3)+S(2,2))/(2*S(2,2)))
          transmat=matmul(transpose(tmp_mat),transmat)
          cycle find_reduced
          !       elseif(cell_type.eq.1.and.S(2,3).lt.0._real12)then
          !          tmp_mat(2,3)=1._real12
          !          transmat=matmul(transpose(tmp_mat),transmat)
          !          cycle find_reduced
       end if


       !! A6
       if(abs(2*S(1,3)).gt.S(1,1)+tiny.or.&
            (abs(2*S(1,3)-S(1,1)).le.tiny.and.2*S(2,3).lt.S(1,2)).or.&
            (abs(2*S(1,3)+S(1,1)).le.tiny.and.S(1,2).lt.0._real12))then
          tmp_mat(1,3)=((-1)**(cell_type+1))*floor((2*S(1,3)+S(1,1))/(2*S(1,1)))
          transmat=matmul(transpose(tmp_mat),transmat)
          cycle find_reduced
          !       elseif(cell_type.eq.1.and.S(1,3).lt.0._real12)then
          !          tmp_mat(1,3)=1._real12
          !          transmat=matmul(transpose(tmp_mat),transmat)
          !          cycle find_reduced
       end if


       !! A7
       if(abs(2*S(1,2)).gt.S(1,1)+tiny.or.&
            (abs(2*S(1,2)-S(1,1)).le.tiny.and.2*S(2,3).lt.S(1,3)).or.&
            (abs(2*S(1,2)+S(1,1)).le.tiny.and.S(1,3).lt.0._real12))then
          tmp_mat(1,2)=((-1)**(cell_type+1))*floor((2*S(1,2)+S(1,1))/(2*S(1,1)))
          transmat=matmul(transpose(tmp_mat),transmat)
          cycle find_reduced
          !       elseif(cell_type.eq.1.and.S(1,2).lt.0._real12)then
          !          tmp_mat(1,2)=1._real12
          !          transmat=matmul(transpose(tmp_mat),transmat)
          !          cycle find_reduced
       end if


       !! A8
       if(cell_type.eq.2.and.&
            2*(S(2,3)+S(1,3)+S(1,2))+S(1,1)+S(2,2).lt.tiny.or.&
            (abs(2*(S(2,3)+S(1,3)+S(1,2))+S(1,1)+S(2,2)).le.tiny.and.&
            2*(S(1,1)+2*S(1,3))+2*S(1,2).gt.tiny))then
          tmp_mat(1,3)=((-1)**(cell_type+1))*floor( &
               ( 2*(S(2,3) + S(1,3) + S(1,2)) + S(1,1) + S(2,2)  )/&
               (2*(2*S(1,2) + S(1,1) + S(2,2) ))  )
          tmp_mat(2,3)=tmp_mat(1,3)
          transmat=matmul(transpose(tmp_mat),transmat)
          cycle find_reduced     
       end if

       lreduced=.true.
    end do find_reduced


    if(abs(det(transmat)+1._real12).le.tiny)then
       tmp_mat=reshape((/-1,0,0,  0,-1,0,  0,0,-1/),shape(tmp_mat))
       transmat=matmul(transpose(tmp_mat),transmat)
    end if
    call mkNiggli_lat(lat,newlat,transmat,S)
    lreduced=reduced_check(newlat,cell_type,S,"n")
    if(verb) then
       write(67,*) lreduced
       write(67,*) (transmat(i,:),i=1,3)
    end if


!!!-----------------------------------------------------------------------------
!!! Renormalises the lattice and basis into the new lattice
!!!-----------------------------------------------------------------------------
    lat=newlat
    do i=1,bas%nspec
       do j=1,bas%spec(i)%num
          bas%spec(i)%atom(j,:3)=&
               matmul(bas%spec(i)%atom(j,:3),inverse_3x3(transmat))
          bas%spec(i)%atom(j,:3)=&
               bas%spec(i)%atom(j,:3)-floor(bas%spec(i)%atom(j,:3))
       end do
    end do


    return
  end subroutine reducer
!!!#############################################################################


!!!#############################################################################
!!! Subroutine to set up the required dot products of the lattice
!!!#############################################################################
!!! a = lat(:,1),  b=lat(:,2),  c=lat(:,3)
!!! S(1,1) = a.a,   S(2,2) = b.b,   S(3,3) = c.c
!!! S(1,2) = a.b,   S(1,3) = a.c,   S(2,3) = b.c
  subroutine mkNiggli_lat(lat,newlat,transmat,S)
    implicit none
    real(real12), dimension(3,3) :: lat,newlat,transmat,S
    real(real12), dimension(3) :: a,b,c


    newlat=matmul(transmat,lat)
    a=newlat(1,:);b=newlat(2,:);c=newlat(3,:)

    S(1,1)=dot_product(a,a)
    S(2,2)=dot_product(b,b)
    S(3,3)=dot_product(c,c)
    S(2,3)=dot_product(b,c)
    S(1,3)=dot_product(a,c)
    S(1,2)=dot_product(a,b)

    return
  end subroutine mkNiggli_lat
!!!#############################################################################


!!!#############################################################################
!!! Function to check whether cell satisfies all the main ...
  ! ... Niggli conditions (1928)
!!!#############################################################################
!!! tiny = tolerance to satisfy conditions
!!! lat = lattice being checked
!!! a = lat(:,1),  b=lat(:,2),  c=lat(:,3)
!!! S(1,1) = a.a,   S(2,2) = b.b,   S(3,3) = c.c
!!! S(1,2) = a.b,   S(1,3) = a.c,   S(2,3) = b.c
!!! Type I  = Sij (i!=j) are all positive (angles <90)
!!! Type II = Sij (i!=j) are all negative or any zero (angles >=90)
!!! Cell is reduced if, and only if, all conditions are ...
!!! ... satisfied (Niggli 1928)
  function reduced_check(lat,cell_type,S,tchar) result(check)
    implicit none
    integer :: cell_type
    real(real12) :: tiny,alpha,beta,gamma,pi2
    real(real12), dimension(3) :: a,b,c
    real(real12), dimension(3,3) :: lat,S
    character(1) :: quiet
    character(1), optional :: tchar
    logical :: check


    quiet="q"
    if(present(tchar)) quiet=tchar
    if(quiet.ne."y".and.quiet.ne."q") quiet="n"

    pi2 = 2._real12*atan(1._real12)
    check=.false.
    tiny=1E-3

    a=lat(1,:);b=lat(2,:);c=lat(3,:)
    S(1,1)=dot_product(a,a)
    S(2,2)=dot_product(b,b)
    S(3,3)=dot_product(c,c)
    S(2,3)=dot_product(b,c)
    S(1,3)=dot_product(a,c)
    S(1,2)=dot_product(a,b)

    alpha=acos(S(2,3)/sqrt(S(2,2)*S(3,3)))
    beta=acos(S(1,3)/sqrt(S(1,1)*S(3,3)))
    gamma=acos(S(1,2)/sqrt(S(1,1)*S(2,2)))

    if(S(2,2)-S(3,3).lt.tiny.and.S(1,1)-S(2,2).lt.tiny) then
       check=.true.
    else
       check=.false.
       return
    end if
    if(cell_type.eq.1.and.&
         alpha.le.pi2.and.beta.le.pi2.and.gamma.le.pi2.and.&
         S(1,2)-0.5_real12*S(1,1).lt.tiny.and.&
         S(1,3)-0.5_real12*S(1,1).lt.tiny.and.&
         S(2,3)-0.5_real12*S(2,2).lt.tiny) then !Type I
       check=.true.
       if(quiet.eq."n") write(0,*) "Found Type I reduced Niggli cell"
    elseif(cell_type.eq.2.and.&
         alpha.ge.pi2-tiny.and.beta.ge.pi2-tiny.and.gamma.ge.pi2-tiny.and.&
         abs(S(1,2))-0.5_real12*S(1,1).lt.tiny.and.&
         abs(S(1,3))-0.5_real12*S(1,1).lt.tiny.and.&
         abs(S(2,3))-0.5_real12*S(2,2).lt.tiny.and.&
         (abs(S(2,3))+abs(S(1,3))+abs(S(1,2)))-0.5_real12*(S(1,1)+S(2,2)).lt.tiny) then !Type II
       if(abs(S(1,2))-0.5_real12*S(1,1).le.tiny.and.S(1,3).gt.tiny) return
       if(abs(S(1,3))-0.5_real12*S(1,1).le.tiny.and.S(1,2).gt.tiny) return
       if(abs(S(2,3))-0.5_real12*S(2,2).le.tiny.and.S(1,2).gt.tiny) return
       if((abs(S(2,3))+abs(S(1,3))+abs(S(1,2)))-0.5_real12*(S(1,1)+S(2,2)).gt.tiny.and.&
            S(1,1)-(2._real12*abs(S(1,3))+abs(S(1,2))).gt.tiny) return
       check=.true.
       if(quiet.eq."n") write(0,*) "Found Type II reduced Niggli cell"
    else
       check=.false.
    end if

    return
  end function reduced_check
!!!#############################################################################


!!!#############################################################################
!!! planecutter
!!!#############################################################################
  function planecutter(inlat,invec) result(tfmat)
    implicit none
    integer :: i,j,itmp1
    real(real12) :: tol
    integer, dimension(3) :: order
    real(real12), dimension(3) :: vec,tvec1
    real(real12), dimension(3,3) :: lat,b,tfmat,invlat,reclat
    real(real12), dimension(3), intent(in) :: invec
    real(real12), dimension(3,3), intent(in) :: inlat



!!!-----------------------------------------------------------------------------
!!! Initialise variables and matrices
!!!-----------------------------------------------------------------------------
    tol=1.E-4
    vec=invec
    lat=inlat
    invlat=inverse(lat)
    reclat=transpose(invlat)*2._real12*pi
    vec=reduce_vec_gcd(vec)
    order=(/1,2,3/)


!!!-----------------------------------------------------------------------------
!!! Align the normal vector such that all non-zero values are left of all zeros
!!!-----------------------------------------------------------------------------
    do i=1,2
       if(vec(i).eq.0)then
          if(all(vec(i:).eq.0._real12)) exit
          itmp1=maxloc(vec(i+1:),mask=vec(i+1:).ne.0,dim=1)+i
          call swap(order(i),order(itmp1))
          call swap(vec(i),vec(itmp1))
          call swap(lat(:,i),lat(:,itmp1))
          call swap(lat(i,:),lat(itmp1,:))
          call swap(reclat(:,i),reclat(:,itmp1))
          call swap(reclat(i,:),reclat(itmp1,:))
       end if
    end do
    !vec=matmul(vec,reclat)


!!!-----------------------------------------------------------------------------
!!! Perform Lenstra-Lenstra-Lovász reduction
!!!-----------------------------------------------------------------------------
    b(1,:) = (/-vec(2),vec(1),0._real12/)
    b(2,:) = (/-vec(3),0._real12,vec(1)/)
    b(3,:) = vec
    tfmat = b
    b(:2,:) = LLL_reduce(b(:2,:))


!!!-----------------------------------------------------------------------------
!!! Checking whether b1 and b2 are still perpendicular to b3 and have size ...
!!! ... greater than zero
!!!-----------------------------------------------------------------------------
    if(dot_product(b(1,:),b(3,:)).gt.tol)then
       write(0,'("ERROR: Internatl error in planecutter")')
       write(0,'(2X,"Error in planecutter subroutine in mod_edit_geom.f90")')
       write(0,'(2X,"b1 not perpendicular to b3")')
       write(0,'(2X,"b1 = ",3(1X,F0.3))') b(1,:)
       write(0,'(2X,"b3 = ",3(1X,F0.3))') b(3,:)
       write(0,'(2X,"b1·b3 = ",F0.3)') dot_product(b(1,:),b(3,:))
       write(0,'("Inform developers of this issue")')
       write(0,'("Stopping...")')
       stop
    elseif(dot_product(b(2,:),b(3,:)).gt.tol)then
       write(0,'("ERROR: Internatl error in planecutter")')
       write(0,'(2X,"Error in planecutter subroutine in mod_edit_geom.f90")')
       write(0,'(2X,"b2 not perpendicular to b3")')
       write(0,'(2X,"b2 = ",3(1X,F6.2))') b(2,:)
       write(0,'(2X,"b3 = ",3(1X,F6.2))') b(3,:)
       write(0,'(2X,"b1·b3 = ",F6.3)') dot_product(b(2,:),b(3,:))
       write(0,'("Inform developers of this issue")')
       write(0,'("Stopping...")')
       stop
    elseif(dot_product(b(1,:),b(1,:)).lt.tol)then
       write(0,'("ERROR: Internatl error in planecutter")')
       write(0,'(2X,"Error in planecutter subroutine in mod_edit_geom.f90")')
       write(0,'(2X,"b1 has zero size")')
       write(0,'(2X,"b1 = ",3(1X,F6.2))') b(1,:)
       write(0,'("Inform developers of this issue")')
       write(0,'("Stopping...")')
       stop
    elseif(dot_product(b(2,:),b(2,:)).lt.tol)then
       write(0,'("ERROR: Internatl error in planecutter")')
       write(0,'(2X,"Error in planecutter subroutine in mod_edit_geom.f90")')
       write(0,'(2X,"b2 has zero size")')
       write(0,'(2X,"b2 = ",3(1X,F6.2))') b(2,:)
       write(0,'("Inform developers of this issue")')
       write(0,'("Stopping...")')
       stop
    end if

    !b = matmul(b,lat)
    

!!!-----------------------------------------------------------------------------
!!! Fix normal vector and lattice
!!!-----------------------------------------------------------------------------
    do i=1,3
       if(i.eq.order(i)) cycle
       call swap(lat(i,:),lat(order(i),:))
       call swap(lat(:,i),lat(:,order(i)))
       call swap(b(:,i),b(:,order(i)))
       call swap(order(order(i)),order(i))
    end do


!!!-----------------------------------------------------------------------------
!!! Convert the new lattice to direct coordinates
!!! Make it such that it is a fully integerised transformation matrix
!!!-----------------------------------------------------------------------------
    !b=matmul(b,invlat)
    where(abs(b(:,:)).lt.tol)
       b(:,:)=0._real12
    end where
    !write(0,'(3(2X,F9.3))') (b(j,:),j=1,3)
    !write(0,*) 
    reduce_loop: do i=1,3
       b(i,:)=reduce_vec_gcd(b(i,:))
       if(any(abs(b(i,:)-nint(b(i,:))).gt.tol))then
          write(0,'("Issue with plane ",3(1X,I0))') nint(invec)
          write(0,*) vec
          write(0,'("row ",I0," of the following matrix")') i
          write(0,'(3(2X,F9.3))') (b(j,:),j=1,3)
          write(0,'(1X,"ERROR: Internal error in planecutter function")')
          write(0,'(2X,"Planecutter in mod_edit_geom.f90 is unable to find a&
               & perpendicular plane")')
          b=0._real12
          exit
       end if
    end do reduce_loop
    if(det(b).lt.0._real12)then
       tvec1=b(2,:)
       b(2,:)=b(1,:)
       b(1,:)=tvec1
    end if
    if(abs(det(b)).lt.tol)then
       write(0,'(1X,"ERROR: Internal error in planecutter function")')
       write(0,'(2X,"Planecutter in mod_edit_geom.f90 has generated a 0&
            & determinant matrix")')
       write(0,'(3(2X,F9.3))') (b(j,:),j=1,3)
       b=0._real12
       !stop
    end if
    tfmat=b


    return
  end function planecutter
!!!#############################################################################


!!!#############################################################################
!!! merges two supplied bases
!!!#############################################################################
!!! Assumes the same lattice for each
  function bas_merge(bas1,bas2,length,map1,map2) result(mergbas)
    implicit none
    integer :: i,j,k,itmp,dim
    logical :: lmap
    integer, allocatable, dimension(:) :: match
    integer, allocatable, dimension(:,:,:) :: new_map

    type(bas_type) :: mergbas
    type(bas_type), intent(in) :: bas1,bas2
    integer, intent(in), optional :: length
    integer, allocatable, dimension(:,:,:), optional, intent(inout) :: map1,map2


    !!--------------------------------------------------------------------------
    !! Set up number of species
    !!--------------------------------------------------------------------------
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


    !!--------------------------------------------------------------------------
    !! If map is present, sets up new map
    !!--------------------------------------------------------------------------
    lmap = .false.
    if_map: if(present(map1).and.present(map2))then
       if(all(map1.eq.-1)) exit if_map
       lmap = .true.
       allocate(new_map(&
            mergbas%nspec,&
            maxval(mergbas%spec(:)%num,dim=1),2))
       new_map = 0
    end if if_map


    !!--------------------------------------------------------------------------
    !! Set up atoms in merged basis
    !!--------------------------------------------------------------------------
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
!!!#############################################################################


!!!#############################################################################
!!! merges two supplied bases and lattices
!!! Does so by stitching one onto the top of the other
!!!#############################################################################
  subroutine bas_lat_merge(merglat,mergbas,inlat1,inlat2,inbas1,inbas2,axis,inoffset,map1,map2)
    implicit none
    integer :: i,k,axis
    real(real12) :: c1_ratio,c2_ratio,add,loc,zgap
    type(bas_type) :: mergbas,bas1,bas2
    type(bas_type), intent(in) :: inbas1,inbas2
    integer, dimension(3) :: order
    real(real12), dimension(3) :: unit_vec,offset
    real(real12), dimension(3), intent(in) :: inoffset
    real(real12), dimension(3,3) :: merglat,lat1,lat2
    real(real12), dimension(3,3), intent(in) :: inlat1,inlat2

    integer, allocatable, dimension(:,:,:), optional, intent(inout) :: map1,map2

    offset=inoffset
    if(allocated(mergbas%spec))then
       do i=1,mergbas%nspec
          if(allocated(mergbas%spec(i)%atom)) deallocate(mergbas%spec(i)%atom)
       end do
       deallocate(mergbas%spec)
    end if

    call clone_bas(inbas1,bas1,inlat1,lat1)
    call clone_bas(inbas2,bas2,inlat2,lat2)
!!!-----------------------------------------------------------------------------
!!! Shifts cells to 
!!!-----------------------------------------------------------------------------
    loc=0._real12
    lat1=MATNORM(lat1)
    add=-min_dist(bas1,axis,loc,.true.)
    call shifter(bas1,axis,add,.true.)

    add=-min_dist(bas2,axis,loc,.true.)
    lat2=MATNORM(lat2)
    call shifter(bas2,axis,add,.true.)


!!!-----------------------------------------------------------------------------
!!! reduces vacuum between materials to desired sizes
!!!-----------------------------------------------------------------------------
    loc=1._real12
    call set_vacuum(lat1,bas1,axis,loc,offset(axis))
    call set_vacuum(lat2,bas2,axis,loc,offset(axis))

    order=(/1,2,3/)
    order=cshift(order,3-axis)
    do k=1,2
       offset(order(k))=offset(order(k))/modu(lat1(order(k),:))
    end do
    unit_vec=uvec(lat1(order(3),:))
    zgap=offset(order(3))/unit_vec(order(3))
    !!NOT SET UP OFFSET FEATURE MADE ABOVE!!! MIGHT BE FIXED NOW! NEED TO TEST

    !loc=1._real12
    !add=zgap+min_dist(bas1,axis,loc)*modu(lat1(axis,:))
    !call vacuumer(lat1,bas1,axis,loc,add)

    !add=zgap+min_dist(bas2,axis,loc)*modu(lat2(axis,:))
    !call vacuumer(lat2,bas2,axis,loc,add)


!!!-----------------------------------------------------------------------------
!!! makes supercell
!!!-----------------------------------------------------------------------------
    merglat(order(1),:)=lat1(order(1),:)
    merglat(order(2),:)=lat1(order(2),:)
    unit_vec=uvec(lat1(axis,:))
    !  slat(axis,:)=lat1(axis,:) + ( modu(lat2(axis,:)) + zgap/unit_vec(axis) )*unit_vec
    merglat(axis,:)=lat1(axis,:) + modu(lat2(axis,:))*unit_vec
    c1_ratio=modu(lat1(axis,:))/modu(merglat(axis,:))
    c2_ratio=modu(lat2(axis,:))/modu(merglat(axis,:))


!!!-----------------------------------------------------------------------------
!!! merge list of atomic types and respective numbers for both structures
!!!-----------------------------------------------------------------------------
    do i=1,bas1%nspec
       bas1%spec(i)%atom(:,axis)=bas1%spec(i)%atom(:,axis)*c1_ratio
    end do
    do i=1,bas2%nspec
       bas2%spec(i)%atom(:,axis)=bas2%spec(i)%atom(:,axis)*c2_ratio + c1_ratio
       do k=1,2
          bas2%spec(i)%atom(:,order(k))=bas2%spec(i)%atom(:,order(k))+offset(order(k))
       end do
    end do

    if(present(map1).and.present(map2))then
       mergbas=bas_merge(bas1,bas2,map1=map1,map2=map2)
    else
       mergbas=bas_merge(bas1,bas2)
    end if
    call normalise_basis(mergbas,1._real12,.true.)


    return
  end subroutine bas_lat_merge
!!!#############################################################################


!!!#############################################################################
!!! splits basis into an array of bases
!!!#############################################################################
  function split_bas(inbas,loc_vec,axis,lall_same_nspec,map1,map2) result(bas_arr)
    implicit none
    integer :: i,is,ia,itmp1,nregions,axis,nspec
    logical :: lsame
    logical :: lmap,lmove
    type(bas_type) :: tbas
    real(real12), allocatable, dimension(:,:) :: dloc_vec
    logical, optional :: lall_same_nspec

    type(bas_type),intent(in) :: inbas
    real(real12), dimension(:,:), intent(in) :: loc_vec
    type(bas_type), allocatable, dimension(:) :: bas_arr

    type map_type   
       integer, allocatable, dimension(:,:,:) :: spec       
    end type map_type
    type(map_type), dimension(2) :: map
    integer, allocatable, dimension(:,:,:), optional, intent(inout) :: map1,map2


    !!--------------------------------------------------------------------------
    !! If map1 is present, sets up map2
    !!--------------------------------------------------------------------------
    if(present(map1).and.present(map2))then
       lmap=.true.
    else
       lmap=.false.
    end if


    !!--------------------------------------------------------------------------
    !! Sets up region locations
    !!--------------------------------------------------------------------------
    lsame=.true.
    if(present(lall_same_nspec)) lsame=lall_same_nspec
    nregions=size(loc_vec(:,1),dim=1)
    allocate(dloc_vec(nregions,2))
    dloc_vec(:,:)=loc_vec(:,:)-floor(loc_vec(:,:))
    where(dloc_vec(:,2).lt.dloc_vec(:,1))
       dloc_vec(:,2)=dloc_vec(:,2)+1._real12
    end where
    allocate(bas_arr(nregions))

    
    !!--------------------------------------------------------------------------
    !! Loops over regions and assigns atoms to them
    !!--------------------------------------------------------------------------
    regionloop1: do i=1,nregions
       bas_arr(i)%natom = 0
       bas_arr(i)%nspec = inbas%nspec
       write(bas_arr(i)%sysname,'(A,"_region_",I0)') trim(inbas%sysname),i
       allocate(bas_arr(i)%spec(inbas%nspec))
       if(lmap) allocate(map(i)%spec(bas_arr(i)%nspec,maxval(inbas%spec(:)%num),2))

       specloop1: do is=1,inbas%nspec
          bas_arr(i)%spec(is)%name=inbas%spec(is)%name
          bas_arr(i)%spec(is)%num = &
               count(inbas%spec(is)%atom(:,axis) - &
               floor(inbas%spec(is)%atom(:,axis)-dloc_vec(i,1))&
               .lt.dloc_vec(i,2))
          bas_arr(i)%natom = bas_arr(i)%natom + bas_arr(i)%spec(is)%num
          allocate(bas_arr(i)%spec(is)%atom(bas_arr(i)%spec(is)%num,3))
          itmp1=0
          atomloop1: do ia=1,inbas%spec(is)%num
             if(inbas%spec(is)%atom(ia,axis) - &
                  floor(inbas%spec(is)%atom(ia,axis)-dloc_vec(i,1))&
                  .lt.dloc_vec(i,2))then
                itmp1=itmp1+1
                bas_arr(i)%spec(is)%atom(itmp1,:3) = inbas%spec(is)%atom(ia,:3)
                if(lmap) map(i)%spec(is,itmp1,:) = map1(is,ia,:)
             end if

          end do atomloop1
       end do specloop1
    end do regionloop1


    !!--------------------------------------------------------------------------
    !! Revmoes null species from regions, if specified
    !!--------------------------------------------------------------------------
    if(.not.lsame)then
       do i=1,nregions
          nspec=0
          lmove=.false.
          tbas%nspec=count(bas_arr(i)%spec(:)%num.gt.0)
          tbas%natom=bas_arr(i)%natom
          tbas%sysname=bas_arr(i)%sysname
          allocate(tbas%spec(tbas%nspec))

          if(lmap.and.i.eq.1)then
             if(allocated(map1)) deallocate(map1)
             allocate(map1(tbas%nspec,size(map(1)%spec(1,:,1),dim=1),2))
             map1=0
          elseif(lmap.and.i.eq.2)then
             if(allocated(map2)) deallocate(map2)
             allocate(map2(tbas%nspec,size(map(2)%spec(1,:,1),dim=1),2))
             map2=0
          end if

          specloop2: do is=1,bas_arr(i)%nspec
             if(bas_arr(i)%spec(is)%num.eq.0)then
                lmove=.true.
                cycle specloop2
             else
                tbas%spec(nspec+1)%num=bas_arr(i)%spec(is)%num
                tbas%spec(nspec+1)%name=bas_arr(i)%spec(is)%name
                if(lmap.and.i.eq.1) map1(nspec+1,:,:) = map(i)%spec(is,:,:)
                if(lmap.and.i.eq.2) map2(nspec+1,:,:) = map(i)%spec(is,:,:)
                call move_alloc(bas_arr(i)%spec(is)%atom,tbas%spec(nspec+1)%atom)
                lmove=.false.                
             end if
             nspec=nspec+1
          end do specloop2
          call clone_bas(tbas,bas_arr(i))
          deallocate(tbas%spec)
       end do
    end if

    


  end function split_bas
!!!#############################################################################


!!!#############################################################################
!!! returns the bulk basis and lattice of 
!!!#############################################################################
  subroutine get_bulk(lat,bas,axis,bulk_lat,bulk_bas)
    implicit none
    integer :: is,ia,ja,len,itmp1
    integer :: minspecloc,minatomloc,nxtatomloc
    real(real12), dimension(3) :: transvec
    real(real12), dimension(2,2) :: regions
    real(real12), dimension(3,3) :: tf
    logical, allocatable, dimension(:) :: atom_mask
    type(bas_type), allocatable, dimension(:) :: splitbas

    integer, intent(in) :: axis
    type(bas_type), intent(in) :: bas
    real(real12), dimension(3,3), intent(in):: lat
    type(bas_type), intent(out) :: bulk_bas
    real(real12), dimension(3,3), intent(out) :: bulk_lat


    minspecloc = minloc(bas%spec(:)%num,mask=bas%spec(:)%num.ne.0,dim=1)
    if(bas%spec(minspecloc)%num.eq.1)then
       write(0,'("ERROR: Internal error in get_bulk")')
       write(0,'(2X,"get_bulk subroutine in mod_edit_geom.f90 unable cannot &
            &find enough atoms to reproduce a bulk from")')
       stop
    end if
    minatomloc = minloc(bas%spec(minspecloc)%atom(:,axis),dim=1)



    allocate(atom_mask(bas%spec(minspecloc)%num))
    atom_mask = .true.



    itmp1 = minatomloc
    where(bas%spec(minspecloc)%atom(:,axis).le.&
         bas%spec(minspecloc)%atom(itmp1,axis))
       atom_mask(:) = .false.
    end where
    region_loop1: do

       atom_mask(itmp1) = .false.
       where(bas%spec(minspecloc)%atom(:,axis).lt.&
            bas%spec(minspecloc)%atom(itmp1,axis))
          atom_mask(:) = .false.
       end where
       nxtatomloc = minloc(&
            abs(&
            bas%spec(minspecloc)%atom(:,axis)-&
            bas%spec(minspecloc)%atom(itmp1,axis)),&
            mask = atom_mask,&
            dim=1)
       write(0,*) minatomloc,nxtatomloc
       
       transvec = bas%spec(minspecloc)%atom(nxtatomloc,:) - &
            bas%spec(minspecloc)%atom(minatomloc,:)

       regions(1,1:2) = [ &
            bas%spec(minspecloc)%atom(minatomloc,axis), &
            bas%spec(minspecloc)%atom(nxtatomloc,axis) ]
       regions(2,1:2) = [ &
            bas%spec(minspecloc)%atom(nxtatomloc,axis), &
            bas%spec(minspecloc)%atom(minatomloc,axis) ]
       splitbas = split_bas(bas,regions,axis)


       spec_loop1: do is=1,bas%nspec
          atom_loop1: do ia=1,splitbas(1)%spec(is)%num

             atom_loop2: do ja=1,splitbas(2)%spec(is)%num

                if( all( abs( ( splitbas(1)%spec(is)%atom(ia,:3) + transvec ) - &
                     splitbas(2)%spec(is)%atom(ja,:3) ).lt.1.E-5 ) )then
                   write(0,*) ia,ja
                   cycle atom_loop1


                end if

             end do atom_loop2
             itmp1 = nxtatomloc
             cycle region_loop1

          end do atom_loop1
       end do spec_loop1
       exit region_loop1

    end do region_loop1


    len=size(bas%spec(1)%atom(1,:))
    allocate(bulk_bas%spec(bas%nspec))
    bulk_bas%nspec=bas%nspec
    bulk_bas%sysname=bas%sysname


    bulk_lat = lat
    bulk_lat(axis,:) = matmul(transvec,lat)


    tf=matmul(inverse(bulk_lat),lat)
    write(0,*) tf
    call clone_bas(splitbas(1),bulk_bas)
    bulk_bas = convert_bas(splitbas(1),tf)





  end subroutine get_bulk
!!!#############################################################################


!!!#############################################################################
!!! returns the atom closest to the centre of the region for a species
!!!#############################################################################
  function get_centre_atom(bas,spec,axis,lw,up) result(iatom)
    implicit none
    integer :: ia
    integer :: iatom
    real(real12) :: dtmp1,dtmp2,centre
    real(real12) :: dlw,dup
    integer, intent(in) :: spec,axis
    real(real12), intent(in) :: lw,up
    type(bas_type), intent(in) :: bas


    iatom=0
    dtmp1 = 1._real12
    if(lw.gt.up)then
       dlw = lw
       dup = 1._real12 + up
    else
       dlw = lw
       dup = up
    end if
    centre = (dlw + dup)/2._real12
    do ia=1,bas%spec(spec)%num
       dtmp2=bas%spec(spec)%atom(ia,axis)&
            -ceiling(bas%spec(spec)%atom(ia,axis)-dup)
       if(abs(dtmp2-centre).lt.dtmp1)then
          dtmp1=abs(dtmp2-centre)
          iatom=ia
       end if
    end do


  end function get_centre_atom
!!!#############################################################################


!!!#############################################################################
!!! returns the atom closest to the location
!!!#############################################################################
  function get_closest_atom_1D(bas,axis,loc,species,above,below) result(atom)
    implicit none
    integer :: is,ia
    integer :: is_start,is_end
    real(real12) :: dtmp1,dtmp2
    logical :: labove,lbelow
    integer, intent(in) :: axis
    real(real12), intent(in) :: loc
    integer, dimension(2) :: atom
    type(bas_type), intent(in) :: bas

    integer, optional, intent(in) :: species
    logical, optional, intent(in) :: above,below

    
    atom=[0,0]
    dtmp1 = 1._real12
    if(present(species))then
       is_start=species
       is_end=species
    else
       is_start=1
       is_end=bas%nspec
    end if
    labove=.false.
    lbelow=.false.
    if(present(above)) labove=above
    if(present(below)) lbelow=below

    do is=is_start,is_end
       atom_loop1: do ia=1,bas%spec(is)%num
          if(labove.and.bas%spec(is)%atom(ia,axis).lt.loc)then
             cycle atom_loop1
          elseif(lbelow.and.bas%spec(is)%atom(ia,axis).gt.loc)then
             cycle atom_loop1
          end if
          dtmp2=bas%spec(is)%atom(ia,axis)&
               -ceiling(bas%spec(is)%atom(ia,axis)-(loc+0.5_real12))
          if(abs(dtmp2-loc).lt.dtmp1)then
             dtmp1=abs(dtmp2-loc)
             atom=[is,ia]
          end if
       end do atom_loop1
    end do
       



  end function get_closest_atom_1D
!!!-----------------------------------------------------
!!!-----------------------------------------------------
  function get_closest_atom_3D(lat,bas,loc,species) result(atom)
    implicit none
    integer :: is,ia
    integer :: is_start,is_end
    real(real12) :: dtmp1,dtmp2
    real(real12), dimension(3) :: vtmp1
    real(real12), dimension(3), intent(in) :: loc
    real(real12), dimension(3,3), intent(in) :: lat
    integer, dimension(2) :: atom
    type(bas_type), intent(in) :: bas

    integer, optional, intent(in) :: species

    
    atom=[0,0]
    dtmp1 = 1._real12
    if(present(species))then
       is_start=species
       is_end=species
    else
       is_start=1
       is_end=bas%nspec
    end if

    spec_loop1: do is=is_start,is_end
       atom_loop1: do ia=1,bas%spec(is)%num
          vtmp1 = bas%spec(is)%atom(ia,:) - loc
          vtmp1 = vtmp1 - ceiling(vtmp1 - 0.5_real12)
          vtmp1 = matmul(vtmp1,lat)
          dtmp2 = modu(vtmp1)
          if(dtmp2.lt.dtmp1)then
             dtmp1=dtmp2
             atom=[is,ia]
          end if
       end do atom_loop1
    end do spec_loop1
    

  end function get_closest_atom_3D
!!!#############################################################################


!!!#############################################################################
!!! returns the largest vacumm gap in the cell along a specified axis
!!!#############################################################################
  function get_largest_gap(lat,bas,axis,tol,return_lower) result(gap)
    implicit none
    integer :: i,init,iloc
    real(real12) :: dtmp1,max_sep,tol_
    logical :: l_lower
    real(real12), dimension(2) :: gap
    real(real12), allocatable, dimension(:) :: bas_list

    integer, intent(in) :: axis
    type(bas_type), intent(in) :: bas
    real(real12), dimension(3,3), intent(in) :: lat

    real(real12), optional, intent(in) :: tol
    logical, optional, intent(in) :: return_lower


!!!-----------------------------------------------------------------------------
!!! handle optional arguments and set defaults if not present
!!!-----------------------------------------------------------------------------
    if(present(return_lower))then
       l_lower=return_lower
    else
       l_lower=.true.
    end if

    if(present(tol))then
       tol_ = tol
    else
       tol_ = 1.E-5_real12
    end if

!!!-----------------------------------------------------------------------------
!!! Set up basis list that will order them wrt distance along 'axis'
!!!-----------------------------------------------------------------------------
    allocate(bas_list(bas%natom))
    init = 1
    do i=1,bas%nspec
       bas_list(init:init+bas%spec(i)%num-1) = bas%spec(i)%atom(:,axis)
       init = init + bas%spec(i)%num
    end do
    call sort1D(bas_list)

!!!-----------------------------------------------------------------------------
!!! cycle through basis to find largest gap
!!!-----------------------------------------------------------------------------
    max_sep = bas_list(1) - ( bas_list(bas%natom) - 1._real12 )
    if(l_lower)then
       iloc =  bas%natom
    else
       iloc = 1
    end if
    do i=1,bas%natom-1
       dtmp1 = bas_list(i+1) - bas_list(i)
       if(dtmp1.gt.max_sep)then
          max_sep = dtmp1
          if(l_lower)then
             iloc = i
          else
             iloc = i + 1
          end if
       end if
    end do

!!!-----------------------------------------------------------------------------
!!! set result variable with edge of gap and size of gap
!!!-----------------------------------------------------------------------------
    gap = [bas_list(iloc)+tol_,max_sep]


  end function get_largest_gap
!!!#############################################################################


!!!#############################################################################
!!! returns the wyckoff atom for each
!!!#############################################################################
  function get_wyckoff(bas,axis) result(wyckoff)
    implicit none
    integer :: is,ia,ja,itmp1,itmp2!ref_atom
    integer :: minspecloc,minatomloc,nxtatomloc
    real(real12) :: up_loc,lw_loc,up_loc2,lw_loc2
    real(real12), dimension(3) :: transvec,tmp_vec1,tmp_vec2,tmp_vec3,tvec
    logical, allocatable, dimension(:) :: atom_mask
    type(wyck_spec_type) :: wyckoff
    integer, intent(in) :: axis
    type(bas_type), intent(in) :: bas

    type l_bulk_type
       logical, allocatable, dimension(:) :: atom
    end type l_bulk_type
    type(l_bulk_type), allocatable, dimension(:) :: l_bulk_atoms

    
    
!!!-----------------------------------------------------------------------------
!!! Finds the species with the minimum number of atoms
!!! Finds upper and lower locations for "slab" and finds atom nearest to the ...
!!! ... centre of that region
!!!-----------------------------------------------------------------------------
    minspecloc = minloc(bas%spec(:)%num,mask=bas%spec(:)%num.ne.0,dim=1)
    minatomloc = minloc(bas%spec(minspecloc)%atom(:,axis),dim=1)
    nxtatomloc = maxloc(bas%spec(minspecloc)%atom(:,axis),dim=1)
    lw_loc = bas%spec(minspecloc)%atom(minatomloc,axis)
    up_loc = bas%spec(minspecloc)%atom(nxtatomloc,axis)
    minatomloc = &
         maxloc(bas%spec(minspecloc)%atom(:,axis)-(lw_loc+up_loc)/2._real12,dim=1,&
         mask=bas%spec(minspecloc)%atom(:,axis)-(lw_loc+up_loc)/2._real12.le.0._real12)
    allocate(atom_mask(bas%spec(minspecloc)%num))
    atom_mask = .true.



!!! INSTEAD OF STARTING FROM BOTTOM, START FROM CLOSEST BELOW MIDDLE AND CLOSEST ABOVE MIDDLE
!!! THEN WORK YOUR WAY OUT FROM THAT GOING 1 BELOW, THEN 1 ABOVE, etc.


!!!-----------------------------------------------------------------------------
!!! Set up lower atom location
!!!-----------------------------------------------------------------------------
    itmp1 = minatomloc
    lw_loc = bas%spec(minspecloc)%atom(minatomloc,axis)
    up_loc = 1._real12
    allocate(l_bulk_atoms(bas%nspec))
    do is=1,bas%nspec
       allocate(l_bulk_atoms(is)%atom(bas%spec(is)%num))
       l_bulk_atoms(is)%atom(:) = .false.
    end do

    
!!!-----------------------------------------------------------------------------
!!! Loops over atoms in cell until it finds a reproducible set to define ...
!!! ... as the bulk
!!!-----------------------------------------------------------------------------
    region_loop1: do
       !!-----------------------------------------------------------------------
       !! Mask of whether an atom has been checked for bulk limits or not
       !!-----------------------------------------------------------------------
       atom_mask(itmp1) = .false.
       where(bas%spec(minspecloc)%atom(:,axis).lt.&
            bas%spec(minspecloc)%atom(itmp1,axis))
          atom_mask(:) = .false.
       end where
       if(all(.not.atom_mask))then
          write(0,'("ERROR: Internal error in get_wyckoff")')
          write(0,'(2X,"Error in subroutine get_wyckoff in mod_edit_geom.f90")')
          write(0,'(2X,"No bulk found")')
          write(0,'(2X,"Exiting subroutine...")')
          return
       end if


       !!-----------------------------------------------------------------------
       !! Defines next atom to use as uper limit for bulk cell
       !!-----------------------------------------------------------------------
       nxtatomloc = minloc(&
            abs(&
            bas%spec(minspecloc)%atom(:,axis)-&
            bas%spec(minspecloc)%atom(itmp1,axis)),&
            mask = atom_mask,&
            dim=1)
       transvec = bas%spec(minspecloc)%atom(nxtatomloc,:) - &
            bas%spec(minspecloc)%atom(minatomloc,:)


       !!-----------------------------------------------------------------------
       !! Checks atoms within a region to see if they reproduce layer above
       !!-----------------------------------------------------------------------
       up_loc = lw_loc + transvec(axis)
       !if(lw_loc.eq.up_loc) up_loc = up_loc + 1.E-8  !! IS THIS NEEDED?
       if(lw_loc.gt.up_loc)then
          write(0,'("ERROR: Internal error in get_wyckoff")')
          write(0,'(2X,"Error in subroutine get_wyckoff in mod_edit_geom.f90")')
          write(0,'(2X,"Region size is negative")')
          write(0,'(2X,"Stopping...")')
          stop
       end if
       spec_loop1: do is=1,bas%nspec
          l_bulk_atoms(is)%atom(:)=.false.
          if(bas%spec(is)%num.eq.0) cycle spec_loop1
          atom_loop1: do ia=1,bas%spec(is)%num
             if(bas%spec(is)%atom(ia,axis).lt.lw_loc.or.&
                  bas%spec(is)%atom(ia,axis).ge.up_loc)then
                cycle atom_loop1
             else
                l_bulk_atoms(is)%atom(ia)=.true.
             end if
             atom_loop2: do ja=1,bas%spec(is)%num
                tmp_vec2 = ( bas%spec(is)%atom(ia,:3) + transvec ) - &
                     bas%spec(is)%atom(ja,:3)
                !! SAME ISSUE HERE AS BELOW
                !! NEED TO TAKE INTO ACCOUNT THAT THEY WORK IN UNISON
                tmp_vec2 = tmp_vec2 - ceiling( tmp_vec2 - 0.5_real12 )


                if( all( abs(tmp_vec2).lt.1.E-5 ) )then
                   cycle atom_loop1
                end if

             end do atom_loop2
             itmp1 = nxtatomloc
             cycle region_loop1

          end do atom_loop1
          if(all(.not.l_bulk_atoms(is)%atom(:)))then
             itmp1 = nxtatomloc
             cycle region_loop1
          end if
       end do spec_loop1


       !!-----------------------------------------------------------------------
       !! Checks the layer above
       !!-----------------------------------------------------------------------
       lw_loc2 = lw_loc + transvec(axis)
       up_loc2 = lw_loc2 + transvec(axis)
       spec_loop2: do is=1,bas%nspec
          if(bas%spec(is)%num.eq.0) cycle spec_loop2
          atom_loop3: do ia=1,bas%spec(is)%num
             if(bas%spec(is)%atom(ia,axis).lt.lw_loc2.or.&
                  bas%spec(is)%atom(ia,axis).ge.up_loc2) cycle atom_loop3
             tmp_vec1 = bas%spec(is)%atom(ia,:3) + transvec
             if( all(bas%spec(is)%atom(:,axis).lt.tmp_vec1(axis)-1.E-5) ) cycle atom_loop3
             atom_loop4: do ja=1,bas%spec(is)%num
                tmp_vec2 = tmp_vec1 - bas%spec(is)%atom(ja,:3)
                tmp_vec2 = tmp_vec2 - ceiling( tmp_vec2 - 0.5_real12 )
                if( all( abs(tmp_vec2).lt.1.E-5 ) )then
                   cycle atom_loop3
                end if
             end do atom_loop4
             itmp1 = nxtatomloc
             cycle region_loop1

          end do atom_loop3
          if(all(.not.l_bulk_atoms(is)%atom(:)))then
             itmp1 = nxtatomloc
             cycle region_loop1
          end if
       end do spec_loop2

       !!-----------------------------------------------------------------------
       !! If it gets to this point, then it has found a bulk cell and exits
       !!-----------------------------------------------------------------------
       exit region_loop1


    end do region_loop1



!!!-----------------------------------------------------------------------------
!!! Using the bulk definition, loop runs through checking which atom maps ...
!!! ... onto which through the bulk translation.
!!! Defines each atom's cell centre wyckoff atom
!!!-----------------------------------------------------------------------------
    allocate(wyckoff%spec(bas%nspec))
    do is=1,bas%nspec
       allocate(wyckoff%spec(is)%atom(bas%spec(is)%num))
       wyckoff%spec(is)%atom(:)=0
       atom_loop5: do ia=1,bas%spec(is)%num
          if(l_bulk_atoms(is)%atom(ia))then
             wyckoff%spec(is)%atom(ia) = ia
          end if
          !write(0,*) is,ia,l_bulk_atoms(is)%atom(ia)
          tmp_vec2 = bas%spec(is)%atom(ia,:3)
          atom_loop6: do ja=1,bas%spec(is)%num
             if(.not.l_bulk_atoms(is)%atom(ja)) cycle atom_loop6
             tmp_vec3 = tmp_vec2 - bas%spec(is)%atom(ja,:3)
             itmp1 = nint(tmp_vec3(axis)/transvec(axis))
             tvec = itmp1*transvec
             tvec = tvec - ceiling(tvec-1._real12)
             !tmp_vec3 = tmp_vec3/transvec
             !tmp_vec3 = reduce_vec_gcd(tmp_vec3)
             itmp2 = nint(get_vec_multiple(tvec,tmp_vec3))

             if(itmp1.eq.0) cycle atom_loop6
             tmp_vec3 = tmp_vec3 - tvec!itmp1*tvec
             tmp_vec3 = tmp_vec3 - ceiling(tmp_vec3 - 0.5_real12)
             !THIS IS WHERE WE NEED TO MAKE IT RIGHT
             !! FIND THE GCD AND DIVIDE
             if(all(abs(tmp_vec3).lt.1.E-5))then
                if(wyckoff%spec(is)%atom(ja).ne.0)then
                   wyckoff%spec(is)%atom(ia) = wyckoff%spec(is)%atom(ja)
                else
                   wyckoff%spec(is)%atom(ia) = ja
                end if
                cycle atom_loop5
             end if
          end do atom_loop6
       end do atom_loop5


       if(any(wyckoff%spec(is)%atom(:).eq.0))then
          write(0,'("ERROR: Internal error in get_wyckoff")')
          write(0,'(2X,"Error in subroutine get_wyckoff in mod_edit_geom.f90")')
          write(0,'(2X,"Not all wyckoff atoms found")')
          do ia=1,bas%spec(is)%num
             write(0,*) is,ia,wyckoff%spec(is)%atom(ia)
          end do
          write(0,'(2X,"Stopping...")')
          stop
       end if


    end do



  end function get_wyckoff
!!!#############################################################################

  
!!!#############################################################################
!!! shares strain between two lattices
!!!#############################################################################
  subroutine share_strain(lat1,lat2,bulk_mod1,bulk_mod2,axis,lcompensate)
    implicit none
    integer :: i
    integer :: axis_
    real(real12) :: area1,area2,delta1,delta2
    integer, dimension(3) :: abc=(/1,2,3/)
    real(real12), dimension(3) :: strain

    real(real12), intent(in) :: bulk_mod1,bulk_mod2
    real(real12), dimension(3,3), intent(inout) :: lat1,lat2

    integer, optional, intent(in) :: axis
    logical, optional, intent(in) :: lcompensate

    axis_=3
    if(present(axis)) axis_=axis
 
    abc=cshift(abc,3-axis_)
    area1 = modu(cross(lat1(abc(1),:),lat1(abc(2),:)))
    area2 = modu(cross(lat2(abc(1),:),lat2(abc(2),:)))
    delta1 = - (1._real12 - area2/area1)/(1._real12 + (area2/area1)*(bulk_mod1/bulk_mod2))
    delta2 = - (1._real12 - area1/area2)/(1._real12 + (area1/area2)*(bulk_mod2/bulk_mod1))
    write(0,*) "areas", area1,area2
    write(0,*) "deltas", delta1,delta2
    write(0,*) "modulus", bulk_mod1,bulk_mod2
    do i=1,3
       if(i.eq.axis_) cycle
       strain(:) = lat1(i,:)-lat2(i,:)
       lat1(i,:) = lat1(i,:) * (1._real12 + delta1)
       lat2(i,:) = lat1(i,:)
    end do
    
    if(present(lcompensate))then
       if(lcompensate)then
          lat1(abc(3),:) =  lat1(abc(3),:) * (1._real12 - delta1/(1._real12 + delta1))  
          lat2(abc(3),:) =  lat2(abc(3),:) * (1._real12 - delta2/(1._real12 + delta2))
       end if
    end if

  end subroutine share_strain
!!!#############################################################################

end module edit_geom
