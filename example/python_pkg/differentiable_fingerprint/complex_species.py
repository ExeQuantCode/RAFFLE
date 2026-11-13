"""Recover CaCu3Ti4O12 species on the supplied structure's known sites."""

import argparse
from pathlib import Path

import torch
from ase.io import read, write

from raffle.fingerprint.differentiable import DifferentiableFingerprint, FingerprintRelaxation

if __package__:
    from .plot_comparison import plot_fingerprint_results
else:
    from plot_comparison import plot_fingerprint_results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--steps", type=int, default=60)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--relax-positions", action="store_true")
    parser.add_argument("--output-dir", type=Path, default=Path("."))
    args = parser.parse_args()
    if args.steps < 1 or args.threads < 1:
        parser.error("--steps and --threads must be positive")
    torch.set_num_threads(args.threads)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    reference = read(Path(__file__).with_name("CaTi4(CuO4)3.poscar"), format="vasp")
    elements = sorted(set(reference.get_chemical_symbols()))
    descriptor = DifferentiableFingerprint(elements, body_orders=(2, 3))
    target = descriptor.generate(reference)
    initial = reference.copy()
    initial.set_chemical_symbols(["Ca"] * len(initial))
    print("Starting from the target's known positions and cell; species are all Ca.")
    print(f"Exact target composition: {target.composition}")
    result = FingerprintRelaxation(
        descriptor,
        target,
        composition_constraint="exact",
        loss_weights=(1.0, 1.0, 0.0),
        loss_reduction="density",
        radial_smoothing=0.3,
        angular_smoothing=0.3,
        species_learning_rate=0.05,
        final_temperature=1.0,
    ).run(
        initial,
        optimise_positions=args.relax_positions,
        max_steps=args.steps,
        final_refine_steps=0,
        tolerance=0,
        verbose=True,
    )
    correct = sum(
        actual == expected
        for actual, expected in zip(result.structure.get_chemical_symbols(), reference.get_chemical_symbols())
    )
    print(f"Initial and final soft loss: {result.loss_history[0]['total']:.6g}, "
          f"{result.loss_history[-1]['total']:.6g}")
    print(f"Final discrete loss: {result.loss:.6g}; exact site labels: {correct}/{len(reference)}")
    print(f"Final composition: {result.composition}")
    write(args.output_dir / "complex_species.extxyz", result.structure)
    write(args.output_dir / "complex_species.traj", result.trajectory)
    plot_fingerprint_results(target, result.fingerprint, args.output_dir / "complex_species.png")


if __name__ == "__main__":
    main()
