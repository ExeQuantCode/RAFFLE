"""Fit a global target fingerprint by moving only one atom."""

import argparse
from pathlib import Path

import numpy as np
from ase.build import bulk
from ase.io import write

from raffle.fingerprint.differentiable import DifferentiableFingerprint, FingerprintRelaxation

if __package__:
    from .plot_comparison import plot_fingerprint_results
else:
    from plot_comparison import plot_fingerprint_results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--atom", type=int, default=0, help="Zero-based mobile atom index")
    parser.add_argument("--steps", type=int, default=80)
    parser.add_argument("--output-dir", type=Path, default=Path("."))
    args = parser.parse_args()
    atoms = bulk("C", "diamond", a=3.567, cubic=True)
    if not 0 <= args.atom < len(atoms) or args.steps < 1:
        parser.error("--atom must be a valid atom index and --steps must be positive")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    descriptor = DifferentiableFingerprint(["C"])
    target = descriptor.generate(atoms)
    initial = atoms.copy()
    initial.positions[args.atom, 0] += 0.3
    fixed = np.ones(len(initial), dtype=bool)
    fixed[args.atom] = False
    result = FingerprintRelaxation(
        descriptor,
        target,
        fixed_position_mask=fixed,
        loss_weights=(1.0, 1.0, 1.0),
        loss_reduction="density",
        radial_smoothing=0.25,
        angular_smoothing=0.2,
        position_learning_rate=0.01,
    ).run(initial, optimise_species=False, max_steps=args.steps,
          final_refine_steps=0, tolerance=0, verbose=True)

    movement = np.linalg.norm(result.structure.positions - initial.positions, axis=1)
    print(f"Mobile atom displacement: {movement[args.atom]:.6g} Angstrom")
    print(f"Maximum fixed-atom displacement: {movement[fixed].max():.6g} Angstrom")
    print(f"Initial and final smoothed loss: {result.loss_history[0]['total']:.6g}, {result.loss:.6g}")
    write(args.output_dir / "local_atom.extxyz", result.structure)
    write(args.output_dir / "local_atom.traj", result.trajectory)
    plot_fingerprint_results(target, result.fingerprint, args.output_dir / "local_atom.png")


if __name__ == "__main__":
    main()
