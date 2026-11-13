"""Fit graphite, FCC and BCC carbon fingerprints starting from diamond."""

import argparse
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.build import bulk
from ase.io import write

from raffle.fingerprint.differentiable import DifferentiableFingerprint, FingerprintRelaxation

if __package__:
    from .plot_comparison import plot_fingerprint_results
else:
    from plot_comparison import plot_fingerprint_results


def graphite():
    a, c = 2.46, 6.70
    return Atoms(
        "C4",
        scaled_positions=[[0, 0, 0], [1 / 3, 2 / 3, 0],
                          [0, 0, 0.5], [2 / 3, 1 / 3, 0.5]],
        cell=[[a, 0, 0], [-a / 2, np.sqrt(3) * a / 2, 0], [0, 0, c]],
        pbc=True,
    ).repeat((2, 1, 1))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--steps", type=int, default=100, help="Coarse steps per phase")
    parser.add_argument("--fine-steps", type=int, default=50, help="Fine steps per phase")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--output-dir", type=Path, default=Path("."))
    args = parser.parse_args()
    if args.steps < 1 or args.fine_steps < 0:
        parser.error("--steps must be positive and --fine-steps non-negative")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    descriptor = DifferentiableFingerprint(
        ["C"], parameters={"radius_distance_tol": [1.5, 2.8, 3.0, 6.0]}
    )
    volume_per_atom = 3.567 ** 3 / 8
    targets = {
        "graphite": graphite(),
        "fcc": bulk("C", "fcc", a=(4 * volume_per_atom) ** (1 / 3), cubic=True).repeat((2, 1, 1)),
        "bcc": bulk("C", "bcc", a=(2 * volume_per_atom) ** (1 / 3), cubic=True).repeat((2, 2, 1)),
    }
    rng = np.random.default_rng(args.seed)
    diamond = bulk("C", "diamond", a=3.567, cubic=True)
    diamond.positions += rng.normal(0.0, 0.15, diamond.positions.shape)
    for phase, target_atoms in targets.items():
        initial = diamond.copy()
        target = descriptor.generate(target_atoms)
        initial_fp = descriptor.generate(initial)
        initial_mse = [
            float((actual - expected).square().mean())
            for actual, expected in zip(
                (initial_fp.two_body, initial_fp.three_body, initial_fp.four_body),
                (target.two_body, target.three_body, target.four_body),
            )
        ]
        print(f"\ndiamond -> {phase}: fixed diamond cell, initial MSE = {initial_mse}", flush=True)
        coarse = FingerprintRelaxation(
            descriptor, target, loss_weights=(1.0, 1.0, 1.0),
            radial_smoothing=0.5, angular_smoothing=0.3,
            loss_reduction="density", position_learning_rate=0.02,
        ).run(initial, optimise_species=False, max_steps=args.steps,
              final_refine_steps=0, tolerance=0, verbose=True)
        result = coarse
        trajectory = coarse.trajectory
        if args.fine_steps:
            fine = FingerprintRelaxation(
                descriptor, target, position_learning_rate=0.01,
            ).run(coarse.structure, optimise_species=False, max_steps=args.fine_steps,
                  final_refine_steps=0, tolerance=0, verbose=True)
            coarse_mse = sum(
                float((actual - expected).square().mean().detach())
                for actual, expected in zip(
                    (coarse.fingerprint.two_body, coarse.fingerprint.three_body, coarse.fingerprint.four_body),
                    (target.two_body, target.three_body, target.four_body),
                )
            )
            if fine.loss <= coarse_mse:
                result = fine
                trajectory = trajectory + fine.trajectory[1:]
            else:
                print("Fine stage increased unsmoothed loss; retaining the coarse structure.")
        final_mse = [
            float((actual - expected).square().mean().detach())
            for actual, expected in zip(
                (result.fingerprint.two_body, result.fingerprint.three_body, result.fingerprint.four_body),
                (target.two_body, target.three_body, target.four_body),
            )
        ]
        print(f"diamond -> {phase}: final unsmoothed MSE = {final_mse}", flush=True)
        write(args.output_dir / f"diamond_to_{phase}.extxyz", result.structure)
        write(args.output_dir / f"diamond_to_{phase}.traj", trajectory)
        plot_fingerprint_results(
            target, result.fingerprint, args.output_dir / f"diamond_to_{phase}.png"
        )


if __name__ == "__main__":
    main()
