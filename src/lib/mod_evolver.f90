module evolver
  use constants, only: real12, pi
  use misc, only: set
  use misc_linalg, only: get_angle, cross, modu
  use rw_geom, only: bas_type
  implicit none

  private


  public :: gvector_type, gvector_container_type
  

  type :: gvector_type
     integer :: num_atoms = 0
     real(real12) :: energy = 0.0_real12 !! should be formation energy
     character(len=3), dimension(:), allocatable :: species
     real(real12), dimension(:,:), allocatable :: df_2body
     real(real12), dimension(:,:), allocatable :: df_3body
     real(real12), dimension(:,:), allocatable :: df_4body
   contains
     procedure, pass(this) :: calculate
  end type gvector_type

!!! KEEP THE GVECTOR AS A DERIVED TYPE

  !! should gvector be for the entire prediction, or one for each system?
  !! if one for each system, do not contain nbins, width, as this would ...
  !! ... result in lots of duplicated data
  type :: gvector_container_type
     integer :: best_system = 0
     integer :: nbins_2body, nbins_3body, nbins_4body
     real(real12) :: width_2body, width_3body, width_4body
     real(real12) :: r_max = 6._real12, r_min = 0._real12
     real(real12) :: theta_max = pi, theta_min = 0._real12
     real(real12) :: phi_max = pi/2._real12, phi_min = 0._real12
     type(gvector_type) :: total !! name it best instead?
     type(gvector_type), dimension(:), allocatable :: system
   contains
     procedure, pass(this) :: evolve
  !   procedure :: read
  end type gvector_container_type

  !interface gvector_container_type
  !   module function init_gvector( &
  !        nbins_2body, nbins_3body, nbins_4body, &
  !        width_2body, width_3body, width_4body, &
  !        r_max, r_min, theta_max, theta_min, phi_max, phi_min &
  !        ) result(gvector)
  !   end function init_gvector
  !end interface gvector_container_type


contains

!!!#############################################################################
!!! 
!!!#############################################################################
  !!! TAKE gvector_type OR bas_type AND WORK FROM THERE
  !!! THEREFORE, ALSO ALLOWS FOR xml FILE
  subroutine evolve(this, system)
    implicit none
    class(gvector_container_type), intent(inout) :: this
    type(gvector_type), dimension(..), intent(in) :: system

    integer :: idx1, idx2
    integer :: i, j, is, js, num_structures_previous
    real(real12) :: weight
    real(real12), dimension(:), allocatable :: height
    integer, dimension(:,:), allocatable :: idx_list
    

    !!! SELECT TYPE OF THE SYSTEM
    !!! IF TYPE(gvector_type)
    !!! IF TYPE(bas_type)
    select rank(system)
    rank(0)
       this%system = [ this%system, system ]
       if( system%energy .lt. this%total%energy ) then
          this%total%energy = system%energy / system%num_atoms
          this%best_system = size(this%system)
       end if
    rank(1)
       num_structures_previous = size(this%system)
       this%system = [ this%system, system ]
       do i = 1, size(system)
          if( system(i)%energy .lt. this%total%energy ) then
             this%total%energy = system(i)%energy
             this%best_system = num_structures_previous + i
          end if
       end do
    rank default
       write(*,*) "ERROR: Invalid rank for system"
       write(*,*) "Expected rank 0 or 1, got ", rank(system)
       stop 1
    end select

    !! get list of species in dataset
    do i = 1, size(this%system)
        this%total%species = [ this%total%species, this%system(i)%species ]
    end do
    call set(this%total%species)

    !!! EASIER TO STORE THE LIST OF LENGTHS, AND ANGLES, OR THE INDIVIDUAL SYSTEM GVECTORS?
    this%total%df_2body = 0._real12
    this%total%df_3body = 0._real12
    this%total%df_4body = 0._real12
    do i = 1, size(this%system)
       
       j = 0
       allocate(idx_list(size(this%system(i)%species),size(this%system(i)%species)))
       do is = 1, size(this%system(i)%species)
          do js = is, size(this%system(i)%species), 1
             j = j + 1
             idx_list(is,js) = j
             idx_list(js,is) = j
          end do
       end do
       
       j = 0
       weight = exp( this%total%energy - &
                      this%system(i)%energy / this%system(i)%num_atoms )
       do is = 1, size(this%total%species)

          idx1 = findloc(this%system(i)%species, this%total%species(is), dim=1)
          if(idx1.eq.0) cycle

          height = 1._real12 / ( 1._real12 + this%total%df_3body(:,is) )
          this%total%df_3body(:,is) = this%total%df_3body(:,is) + &
               height * weight * this%system(i)%df_3body(:,idx1)
          
          height = 1._real12 / ( 1._real12 + this%total%df_4body(:,is) )
          this%total%df_4body(:,is) = this%total%df_4body(:,is) + &
               height * weight * this%system(i)%df_4body(:,idx1)
          
          do js = is, size(this%total%species), 1
             j = is * (js - 1) + js
             idx2 = findloc(this%system(i)%species, this%total%species(js), dim=1)
             height = 1._real12 / ( 1._real12 + this%total%df_2body(:,j) ) ** 2._real12
             this%total%df_2body(:,j) = this%total%df_2body(:,j) + &
                  height * weight * this%system(i)%df_2body(:,idx_list(idx1,idx2))
 
          end do
       end do
       deallocate(idx_list)
       !! Need to consider the fact that some structures may not have all ...
       !! ... the species, so their gvectors will be ordered differently
       !! Also, height should be calculated per pair/species, not per system
       !! So loop needed in here to go over every pair/species
       ! 
    end do

  end subroutine evolve
!!!#############################################################################

!!!#############################################################################
!!! 
!!!#############################################################################
  subroutine calculate(this, lattice, basis, &
       nbins, width, cutoff_min, cutoff_max)
    implicit none
    class(gvector_type), intent(inout) :: this
    real(real12), dimension(3,3), intent(in) :: lattice
    type(bas_type), intent(in) :: basis

    integer, dimension(3), intent(in), optional :: nbins
    real(real12), dimension(3), intent(in), optional :: width
    real(real12), dimension(3), intent(in), optional :: cutoff_min, cutoff_max

    integer, dimension(3) :: nbins_
    real(real12), dimension(3) :: width_ = [0.25_real12, pi/8._real12, pi/16._real12] !!! RANDOMLY CHOSEN DEFAULTS FOR NOW
    real(real12), dimension(3) :: cutoff_min_ = [0._real12, 0._real12, 0._real12]
    real(real12), dimension(3) :: cutoff_max_ = [6._real12, pi, pi/2._real12]

    integer :: bin, max_num_steps
    integer :: i, j, k, b
    integer :: is, js, ia, ja
    integer :: num_pairs
    integer :: amax, bmax, cmax
    real(real12) :: rtmp1, rtmp2, fc, weight, scale
    real(real12), dimension(3) :: eta, limit
    real(real12), dimension(3) :: vtmp1, vtmp2, vtmp3, diff
    real(real12), allocatable, dimension(:) :: gvector_tmp, angle

    integer, dimension(3,2) :: loop_limits
    integer, allocatable, dimension(:,:) :: idx

    type :: bond_type
       integer, dimension(2) :: species, atom !!! CONFIRM THIS ORDER IS KEPT
       logical :: skip = .false.
       real(real12), dimension(3) :: vector
    end type bond_type
    type(bond_type), dimension(:), allocatable :: bond_list


    if(present(cutoff_min)) cutoff_min_ = cutoff_min
    if(present(cutoff_max)) cutoff_max_ = cutoff_max
    if(present(width)) width_ = width
    if(present(nbins))then
       nbins_ = nbins
       width_ = (cutoff_max_ - cutoff_min_)/nbins_
    else
       nbins_ = (cutoff_max_ - cutoff_min_)/width_
    end if
    limit = cutoff_max_ - cutoff_min_


    !! combination calculator with repetition
    i = 0
    num_pairs = gamma(real(basis%nspec + 2, real12)) / &
                ( gamma(real(basis%nspec, real12)) * gamma( 3._real12 ) )
    allocate(idx(2,num_pairs))
    this%species = basis%spec(:)%name !!! CHECK LENGTH OF %name !!!
    do is = 1, basis%nspec
       do js = is, basis%nspec, 1
          i = i + 1
          idx(:,i) = [is, js]
       end do
    end do

    !! this is not perfect
    !! won't work for extremely acute/obtuse angle cells
    !! (due to diagonal path being shorter than individual lattice vectors)
    amax = ceiling(cutoff_max_(1)/modu(lattice(1,:)))
    bmax = ceiling(cutoff_max_(1)/modu(lattice(2,:)))
    cmax = ceiling(cutoff_max_(1)/modu(lattice(3,:)))

    allocate(bond_list(0)) !if doesn't work, allocate a dummy bond first
    spec_loop1: do is=1,basis%nspec
       atom_loop1: do ia=1,basis%spec(is)%num
          spec_loop2: do js=is,basis%nspec
             atom_loop2: do ja=1,basis%spec(is)%num
                if(is.eq.js.and.ja.lt.ia) cycle atom_loop2
                diff = basis%spec(is)%atom(ia,:3) -  basis%spec(js)%atom(ja,:3)
                diff = diff - ceiling(diff - 0.5_real12)
                do i=-amax,amax+1,1
                   vtmp1(1) = diff(1) + real(i, real12)
                   do j=-bmax,bmax+1,1
                      vtmp1(2) = diff(2) + real(j, real12)
                      do k=-cmax,cmax+1,1
                         vtmp1(3) = diff(3) + real(k, real12)
                         rtmp1 = modu(matmul(vtmp1,lattice))
                         if( rtmp1 .gt. cutoff_min_(1) - width_(1)/2._real12 .and. &
                             rtmp1 .lt. cutoff_max_(1) + width_(1)/2._real12 )then
                            bond_list = [ bond_list, bond_type([is,js], [ia,ja], .false., rtmp1) ]
                         end if
                      end do
                   end do
                end do
             end do atom_loop2
          end do spec_loop2
       end do atom_loop1
    end do spec_loop1


    !!--------------------------------------------------------------------------
    !! calculate the gaussian width
    !!--------------------------------------------------------------------------
    eta = 1._real12 / ( 2._real12 * width_**2._real12 )
    this%num_atoms = basis%natom


    !!--------------------------------------------------------------------------
    !! build the 2-body gvectors (radial distribution functions)
    !!--------------------------------------------------------------------------
    allocate(this%df_2body(nbins_(1),num_pairs), source = 0._real12)
    allocate(gvector_tmp(nbins_(1)),             source = 0._real12)
    do i = 1, size(bond_list)
       if(bond_list(i)%skip) cycle
       if(abs(modu(bond_list(i)%vector)).lt.1.E-3) cycle
       is = bond_list(i)%species(1)
       js = bond_list(i)%species(2)
       rtmp1 = modu(bond_list(i)%vector)
       !!-----------------------------------------------------------------------
       !! get number of equivalent bonds
       !!-----------------------------------------------------------------------
       scale = 1._real12
       do j = i+1, size(bond_list)
          !! don't need to look at reverse of species, as the ordering will ...
          !! always be enforced by the loop above, i.e. is <= js
          if( is .ne. bond_list(j)%species(1) .or. &
              js .ne. bond_list(j)%species(2) ) cycle
          if(abs(modu(bond_list(j)%vector)-rtmp1).lt.1.E-3)then !!MAKE THIS LINKED TO WIDTH?
             bond_list(j)%skip = .true.
             scale = scale + 1._real12
          end if          
       end do
       !!-----------------------------------------------------------------------
       bin = nint(nbins_(1) * ( ( rtmp1 - cutoff_min(1)) + &
            width_(1)/2._real12 ) / limit(1) )
       if(bin.gt.nbins_(1).or.bin.lt.1) cycle

       fc = 0.5_real12 * cos( pi * ( rtmp1 - cutoff_min(1) ) / limit(1) ) + &
            0.5_real12

       !!-----------------------------------------------------------------------
       !! calculate the gaussian for this bond
       !!-----------------------------------------------------------------------
       gvector_tmp = 0._real12
       max_num_steps = ceiling( sqrt(16._real12/eta(1)) / width )
       loop_limits(:,1) = [ bin, min(nbins_(1), bin + max_num_steps), 1 ]
       loop_limits(:,2) = [ bin - 1, max(1, bin - max_num_steps), -1 ]

       !! do forward and backward loops to add gaussian for larger distances
       do concurrent ( j = 1:2 )
          do concurrent ( b = loop_limits(1,j):loop_limits(2,j) )
             gvector_tmp(b) = gvector_tmp(b) + &
                  exp( -eta * ( rtmp1 - width * real(b) ) ** 2._real12 )
          end do
       end do
       get_pair_index_loop: do j = 1, num_pairs
          if( is .eq. idx(1,j) .and. js .eq. idx(2,j) )then
             k = j
             exit get_pair_index_loop
          end if
       end do get_pair_index_loop
       this%df_2body(:,k) = this%df_2body(:,k) + gvector_tmp * scale
    end do


    !!--------------------------------------------------------------------------
    !! build the 3-body gvectors (angular distribution functions)
    !!--------------------------------------------------------------------------
    deallocate(gvector_tmp)
    allocate(this%df_3body(nbins_(2), basis%nspec), source = 0._real12)
    allocate(gvector_tmp(nbins_(2)),                source = 0._real12)
    do is = 1, basis%nspec
       allocate(angle(0))
       !!-----------------------------------------------------------------------
       !! loop over all bonds to find the first bond
       !!-----------------------------------------------------------------------
       do i = 1, size(bond_list)
          if(abs(modu(bond_list(i)%vector)).lt.1.E-3) cycle
 
          if( is .eq. bond_list(i)%species(1) )then
             vtmp1 = bond_list(i)%vector
             ia = bond_list(i)%atom(1)
          elseif( is .eq. bond_list(i)%species(2) )then
             vtmp1 = -bond_list(i)%vector
             ia = bond_list(i)%atom(2)
          else
             cycle
          end if
        
          !!--------------------------------------------------------------------
          !! loop over all bonds to find the second bond
          !!--------------------------------------------------------------------
          do j = i + 1, size(bond_list)
             if(abs(modu(bond_list(j)%vector)).lt.1.E-3) cycle
             if( is .eq. bond_list(j)%species(1) .and. &
                 ia .eq. bond_list(j)%atom(1) )then
                vtmp2 = bond_list(j)%vector
             elseif( is .eq. bond_list(j)%species(2) .and. &
                     ia .eq. bond_list(j)%atom(2) )then
                vtmp2 = -bond_list(j)%vector
             else
                cycle
             end if
 
             !!-----------------------------------------------------------------
             !! calculate the angle between the two bonds
             !!-----------------------------------------------------------------
             angle = [ angle, get_angle( vtmp1, vtmp2 ) ]

          end do
       end do
       this%df_3body(:,is) = this%df_3body(:,is) + &
            get_angle_gvector(angle, nbins_(2), eta(2), width_(2), &
                              cutoff_min(2), limit(2))
       deallocate(angle)
    end do

  
    !!--------------------------------------------------------------------------
    !! build the 4-body gvectors (angular distribution functions)
    !!--------------------------------------------------------------------------
    deallocate(gvector_tmp)
    allocate(this%df_4body(nbins_(3),basis%nspec), source = 0._real12)
    allocate(gvector_tmp(nbins_(3)),               source = 0._real12)
    do is = 1, basis%nspec
       allocate(angle(0))
       !!-----------------------------------------------------------------------
       !! loop over all bonds to find the first bond
       !!-----------------------------------------------------------------------
       do i = 1, size(bond_list)
          if(abs(modu(bond_list(i)%vector)).lt.1.E-3) cycle
 
          if( is .eq. bond_list(i)%species(1) )then
             vtmp1 = bond_list(i)%vector
             ia = bond_list(i)%atom(1)
          elseif( is .eq. bond_list(i)%species(2) )then
             vtmp1 = -bond_list(i)%vector
             ia = bond_list(i)%atom(2)
          else
             cycle
          end if

          !!--------------------------------------------------------------------
          !! loop over all bonds to find the second bond
          !!--------------------------------------------------------------------
          do j = i + 1, size(bond_list)
             if(abs(modu(bond_list(j)%vector)).lt.1.E-3) cycle
             if( is .eq. bond_list(j)%species(1) .and. &
                 ia .eq. bond_list(j)%atom(1) )then
                vtmp2 = bond_list(j)%vector
             elseif( is .eq. bond_list(j)%species(2) .and. &
                     ia .eq. bond_list(j)%atom(2) )then
                vtmp2 = -bond_list(j)%vector
             else
                cycle
             end if
 
             !!-----------------------------------------------------------------
             !! loop over all bonds to find the third bond
             !!-----------------------------------------------------------------
             do k = j + 1, size(bond_list)
                if(abs(modu(bond_list(k)%vector)).lt.1.E-3) cycle
                if( is .eq. bond_list(k)%species(1) .and. &
                    ia .eq. bond_list(k)%atom(1) )then
                   vtmp3 = bond_list(k)%vector
                elseif( is .eq. bond_list(k)%species(2) .and. &
                        ia .eq. bond_list(k)%atom(2) )then
                   vtmp3 = -bond_list(j)%vector
                else
                   cycle
                end if
 
                !!--------------------------------------------------------------
                !! calculate the angle between the two bonds
                !!--------------------------------------------------------------
                angle = [ angle, get_angle( cross(vtmp1, vtmp2), vtmp3 ) ]

             end do
          end do
       end do
       this%df_4body(:,is) = this%df_4body(:,is) + &
            get_angle_gvector(angle, nbins_(3), eta(3), width_(3), &
                              cutoff_min(3), limit(3))
       deallocate(angle)
    end do

  end subroutine calculate
!!!#############################################################################


!!!#############################################################################
!!! get the gvector for angle distributions
!!!#############################################################################
  pure function get_angle_gvector(angle, nbins, eta, width, cutoff_min, limit) &
       result(gvector)
    implicit none
    integer, intent(in) :: nbins
    real(real12), intent(in) :: eta, width, cutoff_min, limit
    real(real12), dimension(:), intent(in) :: angle
    real(real12), dimension(nbins) :: gvector

    integer :: i, j, b, bin, max_num_steps
    real(real12) :: scale
    real(real12), dimension(nbins) :: gvector_tmp
    real(real12), dimension(:), allocatable :: angle_copy
    integer, dimension(3,2) :: loop_limits

    angle_copy = angle
    do i = 1, size(angle_copy)
       if( abs(angle_copy(i) + 1.E3_real12 ) .lt. 1.E-3 ) cycle
       scale = 1._real12
       do j = i + 1, size(angle_copy)
          if(abs(angle_copy(i)-angle_copy(j)) .lt. 1.E-3 )then
             angle_copy(j) = -1.E3_real12
             scale = scale + 1._real12
          end if
       end do
 
       bin = nint(nbins * ( angle_copy(i) - cutoff_min ) / limit )
       if(bin.gt.nbins.or.bin.lt.1) cycle
       !!------------------------------------------------------------------
       !! calculate the gaussian for this bond
       !!------------------------------------------------------------------
       gvector_tmp = 0._real12
       ! max_num_steps = ceiling( abs( angle_copy(i) - sqrt(16._real12/eta) ) / width )
       max_num_steps = ceiling( sqrt(16._real12/eta) / width )
       loop_limits(:,1) = [ bin, min(nbins, bin + max_num_steps), 1 ]
       loop_limits(:,2) = [ bin - 1, max(1, bin - max_num_steps), -1 ]

       do concurrent ( j = 1:2 )
          do concurrent ( b = loop_limits(1,j):loop_limits(2,j) )
             gvector_tmp(b) = gvector_tmp(b) + &
                  exp( -eta * ( angle_copy(i) - width * real(b) ) ** 2._real12 )
          end do
       end do
       gvector(:) = gvector(:) + gvector_tmp * scale
    end do

  end function get_angle_gvector
!!!#############################################################################

end module evolver