module evolver
  use constants, only: real12, pi
  use misc, only: set
  use misc_linalg, only: get_angle, cross, modu
  use rw_geom, only: bas_type
  use elements, only: element_type, load_elements, elements_database
  implicit none

  private


  public :: gvector_type, gvector_container_type
  
  
  type :: gvector_base_type
     real(real12), dimension(:,:), allocatable :: df_2body
     real(real12), dimension(:,:), allocatable :: df_3body
     real(real12), dimension(:,:), allocatable :: df_4body
  end type gvector_base_type

  type, extends(gvector_base_type) :: gvector_type
     integer :: num_atoms = 0
     real(real12) :: energy = 0.0_real12 !! should be formation energy
     integer, dimension(:), allocatable :: stoichiometry
     character(len=3), dimension(:), allocatable :: species
   contains
     procedure, pass(this) :: calculate
  end type gvector_type


  !! should gvector be for the entire prediction, or one for each system?
  !! if one for each system, do not contain nbins, width, as this would ...
  !! ... result in lots of duplicated data
  type :: gvector_container_type
     integer :: best_system = 0
     real(real12) :: best_energy = 0.0_real12
     integer, dimension(3) :: nbins
     real(real12), dimension(3) :: width = [ 0.25_real12, pi/8._real12, pi/16._real12 ]
     real(real12), dimension(3) :: cutoff_min = [ 0._real12, 0._real12, 0._real12 ]
     real(real12), dimension(3) :: cutoff_max = [ 6._real12, pi, pi/2._real12 ]
     type(gvector_base_type) :: total !! name it best instead?
     type(gvector_type), dimension(:), allocatable :: system
     type(element_type), dimension(:), allocatable :: species_info
   contains
     procedure, pass(this) :: add, add_basis
     procedure, pass(this) :: set_species_list
     procedure, pass(this) :: set_best_energy
     procedure, pass(this) :: evolve
  !   procedure :: read
  end type gvector_container_type

  !interface gvector_container_type
  !   module function init_gvector( &
  !        nbins, width, cutoff_min, cutoff_max &
  !        ) result(gvector)
  !   end function init_gvector
  !end interface gvector_container_type


contains

!!!#############################################################################
!!! add system (basis or gvector) to the container
!!!#############################################################################
  subroutine add(this, system, lattice)
    implicit none
    class(gvector_container_type), intent(inout) :: this
    class(*), dimension(..), intent(in) :: system
    real(real12), dimension(..), intent(in), optional :: lattice

    integer :: i, num_structures_previous
    character(128) :: buffer


    select rank(system)
    rank(0)
       select type(system)
       type is (gvector_type)
          this%system = [ this%system, system ]
       type is (bas_type)
          if(.not. present(lattice))then
             write(*,*) "ERROR: Lattice vectors not provided"
             stop 1
          end if
          select rank(lattice)
          rank(2)
             call this%add_basis(lattice(:,:), system)
          rank default
             write(*,*) "ERROR: Invalid rank for lattice"
             write(buffer,*) rank(lattice)
             write(*,*) "Expected rank 2, got ", trim(buffer)
             stop 1
          end select
       class default
          write(*,*) "ERROR: Invalid type for system"
          write(*,*) "Expected type gvector_type or bas_type"
          stop 1
       end select
    rank(1)
       num_structures_previous = size(this%system)
       select type(system)
       type is (gvector_type)
          this%system = [ this%system, system ]
       type is (bas_type)
          if(.not. present(lattice))then
             write(*,*) "ERROR: Lattice vectors not provided"
             stop 1
          end if
          select rank(lattice)
          rank(2)
             do i = 1, size(system)
                call this%add_basis(lattice(:,:), system(i))
             end do
          rank(3)
             do i = 1, size(system)
                call this%add_basis(lattice(:,:,i), system(i))
             end do
          rank default
             write(*,*) "ERROR: Invalid rank for lattice"
             write(buffer,*) rank(lattice)
             write(*,*) "Expected rank 2, got ", trim(buffer)
             stop 1
          end select
       class default
          write(*,*) "ERROR: Invalid type for system"
          write(*,*) "Expected type gvector_type or bas_type"
          stop 1
       end select
    rank default
       write(*,*) "ERROR: Invalid rank for system"
       write(buffer,*) rank(system)
       write(*,*) "Expected rank 0 or 1, got ", trim(buffer)
       stop 1
    end select
    call this%set_species_list()
    call this%set_best_energy()

  end subroutine add
!!!#############################################################################


!!!#############################################################################
!!! generate gvectors from basis and add to the container
!!!#############################################################################
  subroutine add_basis(this, lattice, basis)
    implicit none
    class(gvector_container_type), intent(inout) :: this
    type(bas_type), intent(in) :: basis
    real(real12), dimension(3,3), intent(in) :: lattice

    integer :: i, num_structures_previous
    type(gvector_type) :: system

    call system%calculate(lattice, basis, this%nbins, this%width, &
                     this%cutoff_min, this%cutoff_max)

    this%system = [ this%system, system ]
  end subroutine add_basis
!!!#############################################################################


!!!#############################################################################
!!! set the species list for the container
!!!#############################################################################
  subroutine set_species_list(this, elements_file)
    implicit none
    class(gvector_container_type), intent(inout) :: this
    character(*), intent(in), optional :: elements_file

    integer :: i, unit
    character(len=3), dimension(:), allocatable :: species_list


    !!--------------------------------------------------------------------------
    !! load the elements database
    !!--------------------------------------------------------------------------
    if(present(elements_file))then
       call load_elements(elements_file)
    elseif(.not.allocated(elements_database))then
       call load_elements()
    end if

    !!--------------------------------------------------------------------------
    !! get list of species in dataset
    !!--------------------------------------------------------------------------
    do i = 1, size(this%system)
       species_list = [ species_list, this%system(i)%species ]
    end do
    call set(species_list)
    allocate(this%species_info(size(species_list)))
    do i = 1, size(species_list)
       call this%species_info(i)%set(species_list(i))
    end do
    
  end subroutine set_species_list
!!!#############################################################################


!!!#############################################################################
!!! set the best energy and system
!!!#############################################################################
  subroutine set_best_energy(this)
    implicit none
    class(gvector_container_type), intent(inout) :: this

    integer :: i, is, idx
    real(real12) :: energy

    do i = 1, size(this%system)
       energy = this%system(i)%energy
       do is = 1, size(this%species_info)
          idx = findloc(this%system(i)%species, this%species_info(is)%name, dim=1)
          energy = energy - this%system(i)%stoichiometry(idx) * this%species_info(is)%energy
       end do
       energy = energy / this%system(i)%num_atoms
       if( energy .lt. this%best_energy ) then
          this%best_energy = energy
          this%best_system = i
       end if
    end do

  end subroutine set_best_energy
!!!#############################################################################


!!!#############################################################################
!!! generate the evolved gvectors
!!!#############################################################################
!!! EASIER TO STORE THE LIST OF LENGTHS, AND ANGLES, OR THE INDIVIDUAL SYSTEM GVECTORS?
  subroutine evolve(this, system)
    implicit none
    class(gvector_container_type), intent(inout) :: this
    type(gvector_type), dimension(..), intent(in), optional :: system

    integer :: idx1, idx2
    integer :: i, j, is, js, num_structures_previous
    real(real12) :: weight, energy
    real(real12), dimension(:), allocatable :: height
    integer, dimension(:,:), allocatable :: idx_list
    

    !!--------------------------------------------------------------------------
    !! if present, add the system to the container
    !!--------------------------------------------------------------------------
    if(present(system)) call this%add(system)


    !!--------------------------------------------------------------------------
    !! initialise the total gvectors
    !!--------------------------------------------------------------------------
    this%total%df_2body = 0._real12
    this%total%df_3body = 0._real12
    this%total%df_4body = 0._real12


    !!--------------------------------------------------------------------------
    !! loop over all systems to calculate the total gvectors
    !!--------------------------------------------------------------------------
    do i = 1, size(this%system)
       !!-----------------------------------------------------------------------
       !! get the list of 2-body species pairs the system
       !!-----------------------------------------------------------------------
       j = 0
       allocate(idx_list(size(this%system(i)%species),size(this%system(i)%species)))
       do is = 1, size(this%system(i)%species)
          do js = is, size(this%system(i)%species), 1
             j = j + 1
             idx_list(is,js) = j
             idx_list(js,is) = j
          end do
       end do
       

       !!-----------------------------------------------------------------------
       !! calculate the weight for the system
       !!-----------------------------------------------------------------------
       energy = this%system(i)%energy
       do is = 1, size(this%species_info)
          idx1 = findloc(this%system(i)%species, this%species_info(is)%name, dim=1)
          energy = energy - this%system(i)%stoichiometry(idx1) * this%species_info(is)%energy
       end do
       energy = energy / this%system(i)%num_atoms
       weight = exp( this%best_energy - energy )
       j = 0
       !!-----------------------------------------------------------------------
       !! loop over all species in the system to add the gvectors
       !!-----------------------------------------------------------------------
       do is = 1, size(this%species_info)

          idx1 = findloc(this%system(i)%species, this%species_info(is)%name, dim=1)
          if(idx1.eq.0) cycle

          height = 1._real12 / ( 1._real12 + this%total%df_3body(:,is) )
          this%total%df_3body(:,is) = this%total%df_3body(:,is) + &
               height * weight * this%system(i)%df_3body(:,idx1)
          
          height = 1._real12 / ( 1._real12 + this%total%df_4body(:,is) )
          this%total%df_4body(:,is) = this%total%df_4body(:,is) + &
               height * weight * this%system(i)%df_4body(:,idx1)
          
          do js = is, size(this%species_info), 1
             j = is * (js - 1) + js
             idx2 = findloc(this%system(i)%species, this%species_info(js)%name, dim=1)
             height = 1._real12 / ( 1._real12 + this%total%df_2body(:,j) ) ** 2._real12
             this%total%df_2body(:,j) = this%total%df_2body(:,j) + &
                  height * weight * this%system(i)%df_2body(:,idx_list(idx1,idx2))
 
          end do
       end do
       deallocate(idx_list)
    end do


    !!--------------------------------------------------------------------------
    !! deallocate the individual gvectors
    !!--------------------------------------------------------------------------
    deallocate(this%system)

  end subroutine evolve
!!!#############################################################################


!!!#############################################################################
!!! calculate the gvectors for a system from its lattice and basis
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


    !!--------------------------------------------------------------------------
    !! initialise optional variables
    !!--------------------------------------------------------------------------
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


    !!--------------------------------------------------------------------------
    !! get the number of pairs of species
    !!--------------------------------------------------------------------------
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


    !!--------------------------------------------------------------------------
    !! get the stoichiometry, energy, and number of atoms
    !!--------------------------------------------------------------------------
    this%stoichiometry = basis%spec(:)%num
    this%energy = basis%energy
    this%num_atoms = basis%natom


    !!--------------------------------------------------------------------------
    !! get the maximum number of lattice vectors to consider
    !!--------------------------------------------------------------------------
    !! this is not perfect
    !! won't work for extremely acute/obtuse angle cells
    !! (due to diagonal path being shorter than individual lattice vectors)
    amax = ceiling(cutoff_max_(1)/modu(lattice(1,:)))
    bmax = ceiling(cutoff_max_(1)/modu(lattice(2,:)))
    cmax = ceiling(cutoff_max_(1)/modu(lattice(3,:)))


    !!--------------------------------------------------------------------------
    !! build the bond list
    !!--------------------------------------------------------------------------
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
       max_num_steps = ceiling( sqrt(16._real12/eta(1)) / width_(1) )
       loop_limits(:,1) = [ bin, min(nbins_(1), bin + max_num_steps), 1 ]
       loop_limits(:,2) = [ bin - 1, max(1, bin - max_num_steps), -1 ]

       !! do forward and backward loops to add gaussian for larger distances
       do concurrent ( j = 1:2 )
          do concurrent ( b = loop_limits(1,j):loop_limits(2,j) )
             gvector_tmp(b) = gvector_tmp(b) + &
                  exp( -eta(1) * ( rtmp1 - width_(1) * real(b) ) ** 2._real12 )
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


    !!--------------------------------------------------------------------------
    !! calculate the gvector for a list of angles
    !!--------------------------------------------------------------------------
    angle_copy = angle
    do i = 1, size(angle_copy)
       if( abs(angle_copy(i) + 1.E3_real12 ) .lt. 1.E-3 ) cycle
       !!-----------------------------------------------------------------------
       !! remove duplicates
       !!-----------------------------------------------------------------------
       scale = 1._real12
       do j = i + 1, size(angle_copy)
          if(abs(angle_copy(i)-angle_copy(j)) .lt. 1.E-3 )then
             angle_copy(j) = -1.E3_real12
             scale = scale + 1._real12
          end if
       end do


       !!-----------------------------------------------------------------------
       !! get the bin closest to the angle
       !!-----------------------------------------------------------------------
       bin = nint(nbins * ( angle_copy(i) - cutoff_min ) / limit )
       if(bin.gt.nbins.or.bin.lt.1) cycle


       !!-----------------------------------------------------------------------
       !! calculate the gaussian for this bond
       !!-----------------------------------------------------------------------
       gvector_tmp = 0._real12
       ! max_num_steps = ceiling( abs( angle_copy(i) - sqrt(16._real12/eta) ) / width )
       max_num_steps = ceiling( sqrt(16._real12/eta) / width )
       loop_limits(:,1) = [ bin, min(nbins, bin + max_num_steps), 1 ]
       loop_limits(:,2) = [ bin - 1, max(1, bin - max_num_steps), -1 ]


       !!-----------------------------------------------------------------------
       !! do forward and backward loops to add gaussian from its centre
       !!-----------------------------------------------------------------------
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