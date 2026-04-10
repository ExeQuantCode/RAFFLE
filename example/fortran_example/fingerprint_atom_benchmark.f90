program fingerprint_atom_benchmark
  !! Benchmark: comparing full-structure vs. single-atom fingerprint generation.
  !!
  !! This example builds a diamond-carbon supercell of varying sizes and
  !! measures the wall-clock time for:
  !!   (a) computing the fingerprint for the whole structure, and
  !!   (b) computing the fingerprint for a single selected atom.
  !!
  !! Usage:
  !!   fpm run fingerprint_atom_benchmark --example --profile release
  !!
  !! Expected result: single-atom time ≈ (1 / natom) × whole-structure time.
  use raffle__constants, only: real32
  use raffle__geom_rw, only: basis_type
  use raffle__distribs, only: distribs_type
  use raffle__distribs_container, only: distribs_container_type
  implicit none

  integer, parameter :: dp = kind(1.0d0)
  !! Double-precision kind for timing variables.

  !-----------------------------------------------------------------------------
  ! Benchmark parameters
  !-----------------------------------------------------------------------------
  integer, parameter :: n_repeat = 5
  !! Number of repeated timing measurements to average.

  !-----------------------------------------------------------------------------
  ! Local variables
  !-----------------------------------------------------------------------------
  type(basis_type) :: basis_unit, basis_super
  type(distribs_container_type) :: distribs
  type(distribs_type) :: fp
  integer :: supercell_size, ir
  integer :: ix, iy, iz, ia
  integer :: natom_super
  real(real32) :: a0
  real(dp) :: t_start, t_end, t_full, t_atom
  integer :: itgt

  !-----------------------------------------------------------------------------
  ! Build the 8-atom diamond unit cell for carbon
  !-----------------------------------------------------------------------------
  a0 = 3.5607451_real32

  basis_unit%nspec = 1
  basis_unit%natom = 8
  basis_unit%energy = -72.213492_real32
  basis_unit%lat(1,:) = [ a0,  0._real32, 0._real32 ]
  basis_unit%lat(2,:) = [ 0._real32,  a0, 0._real32 ]
  basis_unit%lat(3,:) = [ 0._real32, 0._real32,  a0 ]

  allocate(basis_unit%spec(1))
  basis_unit%spec(1)%name = 'C'
  basis_unit%spec(1)%num  = 8
  allocate(basis_unit%spec(1)%atom(8, 3))
  basis_unit%spec(1)%atom(1,:) = [0.000_real32, 0.000_real32, 0.000_real32]
  basis_unit%spec(1)%atom(2,:) = [0.500_real32, 0.500_real32, 0.000_real32]
  basis_unit%spec(1)%atom(3,:) = [0.500_real32, 0.000_real32, 0.500_real32]
  basis_unit%spec(1)%atom(4,:) = [0.000_real32, 0.500_real32, 0.500_real32]
  basis_unit%spec(1)%atom(5,:) = [0.250_real32, 0.250_real32, 0.250_real32]
  basis_unit%spec(1)%atom(6,:) = [0.750_real32, 0.750_real32, 0.250_real32]
  basis_unit%spec(1)%atom(7,:) = [0.750_real32, 0.250_real32, 0.750_real32]
  basis_unit%spec(1)%atom(8,:) = [0.250_real32, 0.750_real32, 0.750_real32]

  !-----------------------------------------------------------------------------
  ! Loop over supercell sizes: 1x1x1, 2x2x2, 3x3x3
  !-----------------------------------------------------------------------------
  write(*,'(70("="))')
  write(*,'(A)') "  Fingerprint generation benchmark  (carbon diamond supercells)"
  write(*,'(70("="))')
  write(*,'(A5, A8, A18, A18, A14)') &
       "S×S×S", "N_atoms", "Full (ms/call)", "Atom (ms/call)", "Speedup"
  write(*,'(70("-"))')

  do supercell_size = 1, 3
    !---------------------------------------------------------------------------
    ! Build the S×S×S supercell
    !---------------------------------------------------------------------------
    call build_supercell(basis_unit, supercell_size, basis_super)
    natom_super = basis_super%natom

    !---------------------------------------------------------------------------
    ! Target atom for single-atom benchmark (middle of the list)
    !---------------------------------------------------------------------------
    itgt = natom_super / 2 + 1

    !---------------------------------------------------------------------------
    ! (a) Time the full-structure fingerprint
    !---------------------------------------------------------------------------
    distribs = distribs_container_type()
    t_full = 0.0_dp
    do ir = 1, n_repeat
      call cpu_time(t_start)
      fp = distribs%generate_fingerprint(basis_super)
      call cpu_time(t_end)
      t_full = t_full + real(t_end - t_start, dp)
    end do
    t_full = t_full / real(n_repeat, dp) * 1.0e3_dp   ! convert to ms

    !---------------------------------------------------------------------------
    ! (b) Time the single-atom fingerprint
    !---------------------------------------------------------------------------
    distribs = distribs_container_type()
    t_atom = 0.0_dp
    do ir = 1, n_repeat
      call cpu_time(t_start)
      fp = distribs%generate_fingerprint(basis_super, atom_index=itgt)
      call cpu_time(t_end)
      t_atom = t_atom + real(t_end - t_start, dp)
    end do
    t_atom = t_atom / real(n_repeat, dp) * 1.0e3_dp   ! convert to ms

    !---------------------------------------------------------------------------
    ! Report results
    !---------------------------------------------------------------------------
    write(*,'(I3,"×",I1,"×",I1, I8, F16.3, F18.3, F12.1, "×")') &
         supercell_size, supercell_size, supercell_size, &
         natom_super, t_full, t_atom, &
         t_full / max(t_atom, 1.0e-12_dp)

    ! clean up supercell for next iteration
    deallocate(basis_super%spec(1)%atom)
    deallocate(basis_super%spec)
  end do

  write(*,'(70("="))')
  write(*,*)
  write(*,'(A)') "Done."

contains

!###############################################################################
  subroutine build_supercell(unit_cell, s, super)
    !! Tile unit_cell into an s×s×s supercell.
    use raffle__geom_rw, only: basis_type
    use raffle__constants, only: real32
    implicit none

    type(basis_type), intent(in)  :: unit_cell
    !! Unit cell to tile.
    integer,          intent(in)  :: s
    !! Supercell multiplier (same in all three directions).
    type(basis_type), intent(out) :: super
    !! Resulting supercell.

    integer :: nu, ns, ix, iy, iz, ia, idx
    real(real32) :: inv_s

    nu = unit_cell%spec(1)%num
    ns = nu * s**3
    inv_s = 1.0_real32 / real(s, real32)

    super%nspec = 1
    super%natom = ns
    super%energy = unit_cell%energy * real(s**3, real32)
    super%lat(1,:) = unit_cell%lat(1,:) * real(s, real32)
    super%lat(2,:) = unit_cell%lat(2,:) * real(s, real32)
    super%lat(3,:) = unit_cell%lat(3,:) * real(s, real32)

    allocate(super%spec(1))
    super%spec(1)%name = unit_cell%spec(1)%name
    super%spec(1)%num  = ns
    allocate(super%spec(1)%atom(ns, 3))

    idx = 0
    do iz = 0, s-1
      do iy = 0, s-1
        do ix = 0, s-1
          do ia = 1, nu
            idx = idx + 1
            super%spec(1)%atom(idx, 1) = &
                 (unit_cell%spec(1)%atom(ia, 1) + real(ix, real32)) * inv_s
            super%spec(1)%atom(idx, 2) = &
                 (unit_cell%spec(1)%atom(ia, 2) + real(iy, real32)) * inv_s
            super%spec(1)%atom(idx, 3) = &
                 (unit_cell%spec(1)%atom(ia, 3) + real(iz, real32)) * inv_s
          end do
        end do
      end do
    end do

  end subroutine build_supercell
!###############################################################################

end program fingerprint_atom_benchmark
