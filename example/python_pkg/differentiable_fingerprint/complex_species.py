"""Recover CaCu3Ti4O12 species on the supplied structure's known sites."""

import argparse
from pathlib import Path

import torch
import numpy as np
from ase.io import read, write

from raffle.fingerprint.differentiable import DifferentiableFingerprint, FingerprintRelaxation

if __package__:
    from .plot_comparison import plot_fingerprint_results
else:
    from plot_comparison import plot_fingerprint_results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--steps", type=int, default=200)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--parallel_workers", type=int, default=1)
    parser.add_argument("--relax-positions", action="store_true")
    parser.add_argument("--output-dir", type=Path, default=Path("."))
    args = parser.parse_args()
    if args.steps < 1 or args.threads < 1:
        parser.error("--steps and --threads must be positive")
    torch.set_num_threads(args.threads)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    reference = read(Path(__file__).with_name("POSCAR_CCTO"), format="vasp")
    all_species = reference.get_chemical_symbols()
    elements = sorted(set(reference.get_chemical_symbols()))
    elements += ["C"]
    descriptor = DifferentiableFingerprint(elements, body_orders=(2, 3, 4), parallel_workers=args.parallel_workers)
    target = descriptor.generate(reference)
    initial = reference.copy()
    initial.set_chemical_symbols(["C"] * len(initial))
    # randomly perturb positions to avoid starting from the target's known positions
    rng = np.random.default_rng(42)
    initial.positions += rng.normal(0.0, 0.25, initial.positions.shape)
    print("Starting from the target's known positions and cell; species are all C.")
    print(f"Exact target composition: {target.composition}")
    result1 = FingerprintRelaxation(
        descriptor,
        target,
        composition_constraint="exact",
        loss_weights=(1.0, 1.0, 1.0),
        loss_reduction="density",
        radial_smoothing=0.5,
        angular_smoothing=0.3,
        species_learning_rate=0.05,
        final_temperature=1.0,
        position_learning_rate=0.1,
        repulsion_weight=10.0,
    ).run(
        initial,
        optimise_positions=True,
        max_steps=args.steps,
        final_refine_steps=0,
        tolerance=0,
        verbose=True,
        # method="adam", exploration_temperature=0.01,
        # denoiser_steps=0, random_seed=42,  # for reproducibility
    )
    # result2 = FingerprintRelaxation(
    #     descriptor,
    #     target,
    #     composition_constraint="exact",
    #     loss_weights=(1.0, 1.0, 0.0),
    #     loss_reduction="density",
    #     radial_smoothing=0.5,
    #     angular_smoothing=0.3,
    #     species_learning_rate=0.05,
    #     final_temperature=1.0,
    #     position_learning_rate=0.1,
    #     repulsion_weight=10.0,
    # ).run(
    #     result1.structure,
    #     optimise_positions=True,
    #     final_refine_steps=0,
    #     max_steps=50,
    #     tolerance=0,
    #     verbose=True,
    #     method="adam",
    # )
    # result = result2
    result = result1
    correct = sum(
        actual == expected
        for actual, expected in zip(result.structure.get_chemical_symbols(), reference.get_chemical_symbols())
    )
    print(f"Initial and final soft loss: {result.loss_history[0]['total']:.6g}, "
          f"{result.loss_history[-1]['total']:.6g}")
    print(f"Final discrete loss: {result.loss:.6g}; exact site labels: {correct}/{len(reference)}")
    print(f"Final composition: {result.composition}")
    write(args.output_dir / "complex_species.extxyz", result.structure)
    write(args.output_dir / "complex_species.traj", result1.trajectory )#+ result2.trajectory)
    plot_fingerprint_results(target, result.fingerprint, args.output_dir / "complex_species.png")


if __name__ == "__main__":
    main()
