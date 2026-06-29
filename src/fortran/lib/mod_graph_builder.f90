module raffle__graph_builder
  !! Module for building graph tensors and topology for RAFFLE.
  !!
  !! This module contains types and procedures for constructing multigraph
  !! topologies and graph tensors used in the RAFFLE neural network.
  use raffle__constants, only: real32
  implicit none

  private

  ! Public API
  public :: topology_type, graph_tensors_type

  !-----------------------------------------------------------------------------
  ! Type definitions
  !-----------------------------------------------------------------------------
  type :: topology_type
     !! Container for multigraph topology data.
     character(len=3), dimension(:), allocatable :: symbols
     !! Element symbols for each atom.
     integer, dimension(:), allocatable :: species_index
     !! Species index for each atom.
     real(real32), dimension(:), allocatable :: atomic_numbers
     !! Atomic numbers for each atom.
     real(real32), dimension(:), allocatable :: covalent_radii
     !! Covalent radii for each atom.

     ! Pair data
     integer, dimension(:,:), allocatable :: pair_image_shift
     !! Image shifts for each pair.
     integer, dimension(:), allocatable :: pair_target_species_index
     !! Target species index for each pair.
     integer, dimension(:,:), allocatable :: pair_index
     !! Pair indices (center, target).
     integer, dimension(:), allocatable :: pair_type_index
     !! Pair type index.
     real(real32), dimension(:), allocatable :: pair_cutoff_weight_3body
     !! Cutoff weights for 3-body.
     real(real32), dimension(:), allocatable :: pair_cutoff_weight_4body
     !! Cutoff weights for 4-body.

     ! Triplet data
     integer, dimension(:,:), allocatable :: triplet_index
     !! Triplet indices (center_atom, species, triplet_id).
     integer, dimension(:,:), allocatable :: triplet_pair_ids
     !! Pair IDs for each triplet.
     integer, dimension(:), allocatable :: triplet_center_index
     !! Center atom index for each triplet.
     integer, dimension(:), allocatable :: triplet_species_index
     !! Species index for each triplet.

     ! Quadruplet data
     integer, dimension(:,:), allocatable :: quadruplet_pair_ids
     !! Pair IDs for each quadruplet (pair_a, pair_b, pair_c).
     integer, dimension(:), allocatable :: quadruplet_species_index
     !! Species index for each quadruplet.

     ! Derived data
     integer :: num_atoms = 0
     !! Number of atoms.
     integer :: num_pairs = 0
     !! Number of pairs.
     integer :: num_triplets = 0
     !! Number of triplets.
     integer :: num_quadruplets = 0
     !! Number of quadruplets.
   contains
     procedure, pass(this) :: allocate_arrays
     procedure, pass(this) :: allocate_atoms
     procedure, pass(this) :: allocate_pairs
     procedure, pass(this) :: allocate_triplets
     procedure, pass(this) :: allocate_quadruplets
     procedure, pass(this) :: ensure_pair_capacity
     procedure, pass(this) :: resize_pairs
     procedure, pass(this) :: resize_triplets
     procedure, pass(this) :: resize_quadruplets
     procedure, pass(this) :: trim_arrays   ! final trimming to actual counts
     procedure, pass(this) :: finalize
  end type topology_type

  type :: graph_tensors_type
     !! Container for graph tensors used in neural network.
     real(real32), dimension(:,:), allocatable :: global_features
     !! Global features (lattice features).
     real(real32), dimension(:,:), allocatable :: atom_node_features
     !! Atom node features.
     real(real32), dimension(:,:), allocatable :: pair_node_features
     !! Pair node features.

     ! Edge data
     integer, dimension(:,:), allocatable :: atom_edge_index
     !! Atom edge indices.
     integer, dimension(:,:), allocatable :: pair_edge_index
     !! Pair edge indices.

     real(real32), dimension(:,:), allocatable :: atom_edge_attr
     !! Atom edge attributes.
     real(real32), dimension(:,:), allocatable :: pair_edge_attr
     !! Pair edge attributes.

     real(real32), dimension(:), allocatable :: atom_edge_weight
     !! Atom edge weights.
     real(real32), dimension(:), allocatable :: pair_edge_weight
     !! Pair edge weights.

     integer, allocatable :: hyperedge_index(:,:)
     real(real32), allocatable :: hyperedge_weight(:)
     real(real32), allocatable :: hyperedge_attr(:,:)
  end type graph_tensors_type


contains

!###############################################################################
  subroutine allocate_arrays(this, num_atoms, num_pairs, &
       num_triplets, num_quadruplets)
    !! Allocate arrays for the topology type.
    implicit none
    class(topology_type), intent(inout) :: this
    integer, intent(in), optional :: num_atoms, num_pairs
    integer, intent(in), optional :: num_triplets, num_quadruplets

    if(present(num_atoms))then
       call this%allocate_atoms(num_atoms)
    end if

    if(present(num_pairs))then
       call this%allocate_pairs(num_pairs)
    end if

    if(present(num_triplets))then
       call this%allocate_triplets(num_triplets)
    end if

    if(present(num_quadruplets))then
       call this%allocate_quadruplets(num_quadruplets)
    end if

  end subroutine allocate_arrays
!-------------------------------------------------------------------------------
  subroutine allocate_atoms(this, num_atoms)
    implicit none
    class(topology_type), intent(inout) :: this
    integer, intent(in) :: num_atoms
    this%num_atoms = num_atoms
    allocate(this%symbols(num_atoms), source="   ")
    allocate(this%species_index(num_atoms), source=0)
    allocate(this%atomic_numbers(num_atoms), source=0._real32)
    allocate(this%covalent_radii(num_atoms), source=0._real32)
  end subroutine allocate_atoms
!-------------------------------------------------------------------------------
  subroutine allocate_pairs(this, num_pairs)
    implicit none
    class(topology_type), intent(inout) :: this
    integer, intent(in) :: num_pairs
    this%num_pairs = num_pairs
    allocate(this%pair_image_shift(num_pairs, 3), source=0)
    allocate(this%pair_target_species_index(num_pairs), source=0)
    allocate(this%pair_index(num_pairs, 2), source=0)
    allocate(this%pair_type_index(num_pairs), source=0)
    allocate(this%pair_cutoff_weight_3body(num_pairs), source=0._real32)
    allocate(this%pair_cutoff_weight_4body(num_pairs), source=0._real32)
  end subroutine allocate_pairs
!-------------------------------------------------------------------------------
  subroutine allocate_triplets(this, num_triplets)
    implicit none
    class(topology_type), intent(inout) :: this
    integer, intent(in) :: num_triplets
    this%num_triplets = num_triplets
    allocate(this%triplet_index(num_triplets, 3), source=0)
    allocate(this%triplet_pair_ids(num_triplets, 2), source=0)
    allocate(this%triplet_center_index(num_triplets), source=0)
    allocate(this%triplet_species_index(num_triplets), source=0)
  end subroutine allocate_triplets
!-------------------------------------------------------------------------------
  subroutine allocate_quadruplets(this, num_quadruplets)
    implicit none
    class(topology_type), intent(inout) :: this
    integer, intent(in) :: num_quadruplets
    this%num_quadruplets = num_quadruplets
    allocate(this%quadruplet_pair_ids(num_quadruplets, 3), source=0)
    allocate(this%quadruplet_species_index(num_quadruplets), source=0)
  end subroutine allocate_quadruplets
!###############################################################################


!###############################################################################
  subroutine ensure_pair_capacity(this, required)
    implicit none
    class(topology_type), intent(inout) :: this
    integer, intent(in) :: required
    integer :: current, new_size
    if(.not. allocated(this%pair_index)) then
       ! First allocation: use a sensible initial size (e.g., 10*num_atoms or required)
       new_size = max(required, 10*this%num_atoms)
       call resize_pairs(this, new_size)
    else
       current = size(this%pair_index, 1)
       if(required > current) then
          new_size = max(required, nint(current * 1.5_real32))  ! 50% growth
          call resize_pairs(this, new_size)
       end if
    end if
  end subroutine ensure_pair_capacity
!###############################################################################


!###############################################################################
  subroutine resize_pairs(this, new_size)
    !! Resize all pair-related arrays to `new_size`.
    !! If arrays are not allocated, they are allocated.
    implicit none
    class(topology_type), intent(inout) :: this
    integer, intent(in) :: new_size
    integer :: old_size

    ! Temporary arrays (same rank/shape as originals)
    integer, allocatable :: tmp_pair_index(:, :)
    integer, allocatable :: tmp_pair_image_shift(:, :)
    integer, allocatable :: tmp_pair_target_species_index(:)
    integer, allocatable :: tmp_pair_type_index(:)
    real(real32), allocatable :: tmp_pair_cutoff_weight_3body(:)
    real(real32), allocatable :: tmp_pair_cutoff_weight_4body(:)

    ! --- Resize pair_index ---
    if(allocated(this%pair_index)) then
       old_size = size(this%pair_index, 1)
       call move_alloc(this%pair_index, tmp_pair_index)
       allocate(this%pair_index(new_size, 2), source=0)
       if(old_size > 0) then
          this%pair_index(1:min(old_size, new_size), :) = &
               tmp_pair_index(1:min(old_size, new_size), :)
       end if
       deallocate(tmp_pair_index)
    else
       allocate(this%pair_index(new_size, 2), source=0)
    end if

    ! --- Resize pair_image_shift ---
    if(allocated(this%pair_image_shift)) then
       old_size = size(this%pair_image_shift, 1)
       call move_alloc(this%pair_image_shift, tmp_pair_image_shift)
       allocate(this%pair_image_shift(new_size, 3), source=0)
       if(old_size > 0) then
          this%pair_image_shift(1:min(old_size, new_size), :) = &
               tmp_pair_image_shift(1:min(old_size, new_size), :)
       end if
       deallocate(tmp_pair_image_shift)
    else
       allocate(this%pair_image_shift(new_size, 3), source=0)
    end if

    ! --- Resize pair_target_species_index ---
    if(allocated(this%pair_target_species_index)) then
       old_size = size(this%pair_target_species_index)
       call move_alloc(this%pair_target_species_index, tmp_pair_target_species_index)
       allocate(this%pair_target_species_index(new_size), source=0)
       if(old_size > 0) then
          this%pair_target_species_index(1:min(old_size, new_size)) = &
               tmp_pair_target_species_index(1:min(old_size, new_size))
       end if
       deallocate(tmp_pair_target_species_index)
    else
       allocate(this%pair_target_species_index(new_size), source=0)
    end if

    ! --- Resize pair_type_index ---
    if(allocated(this%pair_type_index)) then
       old_size = size(this%pair_type_index)
       call move_alloc(this%pair_type_index, tmp_pair_type_index)
       allocate(this%pair_type_index(new_size), source=0)
       if(old_size > 0) then
          this%pair_type_index(1:min(old_size, new_size)) = &
               tmp_pair_type_index(1:min(old_size, new_size))
       end if
       deallocate(tmp_pair_type_index)
    else
       allocate(this%pair_type_index(new_size), source=0)
    end if

    ! --- Resize pair_cutoff_weight_3body ---
    if(allocated(this%pair_cutoff_weight_3body)) then
       old_size = size(this%pair_cutoff_weight_3body)
       call move_alloc(this%pair_cutoff_weight_3body, tmp_pair_cutoff_weight_3body)
       allocate(this%pair_cutoff_weight_3body(new_size), source=0._real32)
       if(old_size > 0) then
          this%pair_cutoff_weight_3body(1:min(old_size, new_size)) = &
               tmp_pair_cutoff_weight_3body(1:min(old_size, new_size))
       end if
       deallocate(tmp_pair_cutoff_weight_3body)
    else
       allocate(this%pair_cutoff_weight_3body(new_size), source=0._real32)
    end if

    ! --- Resize pair_cutoff_weight_4body ---
    if(allocated(this%pair_cutoff_weight_4body)) then
       old_size = size(this%pair_cutoff_weight_4body)
       call move_alloc(this%pair_cutoff_weight_4body, tmp_pair_cutoff_weight_4body)
       allocate(this%pair_cutoff_weight_4body(new_size), source=0._real32)
       if(old_size > 0) then
          this%pair_cutoff_weight_4body(1:min(old_size, new_size)) = &
               tmp_pair_cutoff_weight_4body(1:min(old_size, new_size))
       end if
       deallocate(tmp_pair_cutoff_weight_4body)
    else
       allocate(this%pair_cutoff_weight_4body(new_size), source=0._real32)
    end if

  end subroutine resize_pairs
!-------------------------------------------------------------------------------
  subroutine resize_triplets(this, new_size)
    implicit none
    class(topology_type), intent(inout) :: this
    integer, intent(in) :: new_size
    integer :: old_size
    integer, allocatable :: tmp_triplet_index(:, :)
    integer, allocatable :: tmp_triplet_pair_ids(:, :)
    integer, allocatable :: tmp_triplet_center_index(:)

    if(allocated(this%triplet_index)) then
       old_size = size(this%triplet_index, 1)
       call move_alloc(this%triplet_index, tmp_triplet_index)
       allocate(this%triplet_index(new_size, 3))
       if(old_size > 0) then
          this%triplet_index(1:min(old_size, new_size), :) = &
               tmp_triplet_index(1:min(old_size, new_size), :)
       end if
       deallocate(tmp_triplet_index)
    else
       allocate(this%triplet_index(new_size, 3))
    end if

    if(allocated(this%triplet_pair_ids)) then
       old_size = size(this%triplet_pair_ids, 1)
       call move_alloc(this%triplet_pair_ids, tmp_triplet_pair_ids)
       allocate(this%triplet_pair_ids(new_size, 2))
       if(old_size > 0) then
          this%triplet_pair_ids(1:min(old_size, new_size), :) = &
               tmp_triplet_pair_ids(1:min(old_size, new_size), :)
       end if
       deallocate(tmp_triplet_pair_ids)
    else
       allocate(this%triplet_pair_ids(new_size, 2))
    end if

    if(allocated(this%triplet_center_index)) then
       old_size = size(this%triplet_center_index)
       call move_alloc(this%triplet_center_index, tmp_triplet_center_index)
       allocate(this%triplet_center_index(new_size))
       if(old_size > 0) then
          this%triplet_center_index(1:min(old_size, new_size)) = &
               tmp_triplet_center_index(1:min(old_size, new_size))
       end if
       deallocate(tmp_triplet_center_index)
    else
       allocate(this%triplet_center_index(new_size))
    end if

    if(allocated(this%triplet_species_index)) then
       old_size = size(this%triplet_species_index)
       call move_alloc(this%triplet_species_index, tmp_triplet_center_index)
       allocate(this%triplet_species_index(new_size))
       if(old_size > 0) then
          this%triplet_species_index(1:min(old_size, new_size)) = &
               tmp_triplet_center_index(1:min(old_size, new_size))
       end if
       deallocate(tmp_triplet_center_index)
    else
       allocate(this%triplet_species_index(new_size))
    end if

  end subroutine resize_triplets
!-------------------------------------------------------------------------------
  subroutine resize_quadruplets(this, new_size)
    implicit none
    class(topology_type), intent(inout) :: this
    integer, intent(in) :: new_size
    integer :: old_size
    integer, allocatable :: tmp_quadruplet_pair_ids(:, :)
    integer, allocatable :: tmp_quadruplet_species_index(:)

    if(allocated(this%quadruplet_pair_ids)) then
       old_size = size(this%quadruplet_pair_ids, 1)
       call move_alloc(this%quadruplet_pair_ids, tmp_quadruplet_pair_ids)
       allocate(this%quadruplet_pair_ids(new_size, 3))
       if(old_size > 0) then
          this%quadruplet_pair_ids(1:min(old_size, new_size), :) = &
               tmp_quadruplet_pair_ids(1:min(old_size, new_size), :)
       end if
       deallocate(tmp_quadruplet_pair_ids)
    else
       allocate(this%quadruplet_pair_ids(new_size, 3))
    end if

    if(allocated(this%quadruplet_species_index)) then
       old_size = size(this%quadruplet_species_index)
       call move_alloc(this%quadruplet_species_index, tmp_quadruplet_species_index)
       allocate(this%quadruplet_species_index(new_size))
       if(old_size > 0) then
          this%quadruplet_species_index(1:min(old_size, new_size)) = &
               tmp_quadruplet_species_index(1:min(old_size, new_size))
       end if
       deallocate(tmp_quadruplet_species_index)
    else
       allocate(this%quadruplet_species_index(new_size))
    end if

  end subroutine resize_quadruplets
!###############################################################################


!###############################################################################
  subroutine trim_arrays(this)
    !! Trim **all** allocated arrays to the current `num_*` counters.
    !! Assumes that `num_pairs`, `num_triplets`, `num_quadruplets`
    !! hold the exact used counts.
    implicit none
    class(topology_type), intent(inout) :: this

    if(this%num_pairs .gt. 0) then
       call this%resize_pairs(this%num_pairs)
    end if
    if(this%num_triplets .gt. 0) then
       call this%resize_triplets(this%num_triplets)
    end if
    if(this%num_quadruplets .gt. 0) then
       call this%resize_quadruplets(this%num_quadruplets)
    end if
  end subroutine trim_arrays
!###############################################################################


!###############################################################################
  subroutine finalize(this)
    !! Finalize the topology type (trim arrays to actual size).
    implicit none
    class(topology_type), intent(inout) :: this

    call this%trim_arrays()

  end subroutine finalize
!###############################################################################

end module raffle__graph_builder
