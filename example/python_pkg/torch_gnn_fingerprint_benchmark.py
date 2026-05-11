"""Benchmark the PyTorch multigraph fingerprint surrogate on carbon data."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from ase.build import bulk
from ase.io import read, write

from raffle import symmetry_aware_rmsd
from raffle.torch_gnn_fingerprint import TorchGNNFingerprint
from torch_gnn_workflow_common import save_descriptor_comparison_report


PERTURBATION = np.array(
    [
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.00172792, 0.00410809, 0.00165219],
        [-0.00651579, 0.00452678, 0.00223187],
        [-0.00268477, 0.00290559, 0.00182286],
        [0.00147066, 0.00014211, 0.00273356],
    ],
    dtype=np.float32,
)
def build_augmented_structures(original, count: int, seed: int) -> list:
    rng = np.random.default_rng(seed)
    structures = []
    for _ in range(count):
        atoms = original.copy()
        atoms.set_positions(
            atoms.get_positions()
            + rng.normal(scale=0.04, size=atoms.positions.shape).astype(np.float32)
        )
        structures.append(atoms)
    return structures


def build_perturbed_structure(original):
    perturbed = original.copy()
    perturbed.set_positions(perturbed.get_positions() + PERTURBATION)
    fixed_atoms = np.zeros(len(perturbed), dtype=bool)
    fixed_atoms[:4] = True
    return perturbed, fixed_atoms


def run_workflow(
    repo_root: Path,
    output_dir: Path,
    carbon_count: int = 48,
    augmented_count: int = 16,
    num_epochs: int = 30,
    batch_size: int = 8,
    inverse_steps: int = 180,
    inverse_step_size: float = 5.0e-3,
    seed: int = 42,
) -> dict:
    carbon_xyz = repo_root / "example" / "data" / "carbon.xyz"
    carbon_structures = read(str(carbon_xyz), index=":")[:carbon_count]

    original = bulk("C", "diamond", a=3.567, cubic=True)
    original.pbc = True
    augmented_structures = build_augmented_structures(original, augmented_count, seed)
    perturbed, fixed_atoms = build_perturbed_structure(original)

    model = TorchGNNFingerprint(
        species_list=["C"],
        hidden_dim=80,
        num_message_layers=2,
        learning_rate=5.0e-4,
        lr_decay_rate=5.0e-3,
        smooth_cutoff_width=0.2,
        seed=seed,
    )

    reference_perturbed = model.compute_reference_fingerprint(perturbed)
    initial_prediction = model.predict(perturbed)
    initial_eval_mse = float(np.mean((initial_prediction - reference_perturbed) ** 2))

    training_history = model.fit(
        carbon_structures,
        num_epochs=num_epochs,
        batch_size=batch_size,
        augment_structures=augmented_structures,
        verbose=0,
    )

    trained_prediction = model.predict(perturbed)
    trained_eval_mse = float(np.mean((trained_prediction - reference_perturbed) ** 2))

    grad2, grad3, grad4 = model.compute_gradients(perturbed)
    target_fingerprint = model.compute_reference_fingerprint(original)
    initial_inverse_mse = float(np.mean((trained_prediction - target_fingerprint) ** 2))

    optimisation_history = []
    current_atoms = perturbed.copy()
    current_prediction = trained_prediction
    optimisation_history.append(
        {
            "step": 0,
            "fingerprint_mse": initial_inverse_mse,
            "rmsd": symmetry_aware_rmsd(original, current_atoms),
        }
    )
    for step in range(1, inverse_steps + 1):
        current_atoms = model.inverse_design(
            target_fingerprint=target_fingerprint,
            atoms=current_atoms,
            fixed_atoms=fixed_atoms,
            num_steps=1,
            step_size=inverse_step_size,
            verbose=0,
        )
        current_prediction = model.predict(current_atoms)
        optimisation_history.append(
            {
                "step": step,
                "fingerprint_mse": float(np.mean((current_prediction - target_fingerprint) ** 2)),
                "rmsd": symmetry_aware_rmsd(original, current_atoms),
            }
        )

    optimised = current_atoms
    final_prediction = current_prediction
    final_inverse_mse = float(np.mean((final_prediction - target_fingerprint) ** 2))
    final_rmsd = symmetry_aware_rmsd(original, optimised)

    output_dir.mkdir(parents=True, exist_ok=True)
    write(output_dir / "torch_gnn_diamond_original.xyz", original)
    write(output_dir / "torch_gnn_diamond_perturbed.xyz", perturbed)
    write(output_dir / "torch_gnn_diamond_optimised.xyz", optimised)
    descriptor_report = save_descriptor_comparison_report(
        model=model,
        target_fingerprint=target_fingerprint,
        final_atoms=optimised,
        structure_path=output_dir / "torch_gnn_diamond_optimised.xyz",
    )

    figure = plt.figure(figsize=(14, 4.5))
    ax1 = figure.add_subplot(1, 3, 1)
    ax1.plot(range(len(training_history)), training_history, marker="o")
    ax1.set_xlabel("Epoch")
    ax1.set_ylabel("Joint weighted MSE")
    ax1.set_title("Training convergence")

    ax2 = figure.add_subplot(1, 3, 2)
    ax2.plot(
        [entry["step"] for entry in optimisation_history],
        [entry["fingerprint_mse"] for entry in optimisation_history],
        marker="o",
        label="fingerprint MSE",
    )
    ax2.plot(
        [entry["step"] for entry in optimisation_history],
        [entry["rmsd"] for entry in optimisation_history],
        marker="s",
        label="RMSD",
    )
    ax2.set_xlabel("Inverse-design step")
    ax2.set_title("Inverse design")
    ax2.legend(loc="best")

    ax3 = figure.add_subplot(1, 3, 3, projection="3d")
    original_positions = original.get_positions()
    optimised_positions = optimised.get_positions()
    ax3.scatter(
        original_positions[:, 0],
        original_positions[:, 1],
        original_positions[:, 2],
        s=35,
        label="original",
    )
    ax3.scatter(
        optimised_positions[:, 0],
        optimised_positions[:, 1],
        optimised_positions[:, 2],
        s=35,
        marker="^",
        label="optimised",
    )
    ax3.set_title("Structure comparison")
    ax3.legend(loc="best")

    figure.tight_layout()
    figure_path = output_dir / "torch_gnn_fingerprint_benchmark.png"
    figure.savefig(figure_path, dpi=150)
    plt.close(figure)

    metrics = {
        "num_carbon_structures": len(carbon_structures),
        "num_augmented_structures": len(augmented_structures),
        "component_dims": {
            "2body": model.fingerprint_dim_2body,
            "3body": model.fingerprint_dim_3body,
            "4body": model.fingerprint_dim_4body,
        },
        "training_history": training_history,
        "initial_eval_mse": initial_eval_mse,
        "trained_eval_mse": trained_eval_mse,
        "gradient_shapes": {
            "2body": list(grad2.shape),
            "3body": list(grad3.shape),
            "4body": list(grad4.shape),
        },
        "gradient_norms": {
            "2body": float(np.linalg.norm(grad2)),
            "3body": float(np.linalg.norm(grad3)),
            "4body": float(np.linalg.norm(grad4)),
        },
        "initial_inverse_mse": initial_inverse_mse,
        "final_inverse_mse": final_inverse_mse,
        "final_rmsd": final_rmsd,
        "optimisation_history": optimisation_history,
        "output_files": {
            "plot": str(figure_path),
            "original": str(output_dir / "torch_gnn_diamond_original.xyz"),
            "perturbed": str(output_dir / "torch_gnn_diamond_perturbed.xyz"),
            "optimised": str(output_dir / "torch_gnn_diamond_optimised.xyz"),
            "optimised_descriptor_comparison": descriptor_report["report_file"],
            "optimised_descriptor_comparison_plot": descriptor_report["plot_file"],
        },
        "descriptor_comparisons": {"final": descriptor_report},
    }
    metrics_path = output_dir / "torch_gnn_fingerprint_benchmark_metrics.json"
    metrics_path.write_text(json.dumps(metrics, indent=2))
    return metrics


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--carbon-count", type=int, default=48)
    parser.add_argument("--augmented-count", type=int, default=16)
    parser.add_argument("--epochs", type=int, default=30)
    parser.add_argument("--batch-size", type=int, default=8)
    parser.add_argument("--inverse-steps", type=int, default=180)
    parser.add_argument("--inverse-step-size", type=float, default=5.0e-3)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("build") / "torch_gnn_fingerprint_benchmark",
    )
    parser.add_argument("--seed", type=int, default=42)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    repo_root = Path(__file__).resolve().parents[2]
    metrics = run_workflow(
        repo_root=repo_root,
        output_dir=args.output_dir,
        carbon_count=args.carbon_count,
        augmented_count=args.augmented_count,
        num_epochs=args.epochs,
        batch_size=args.batch_size,
        inverse_steps=args.inverse_steps,
        inverse_step_size=args.inverse_step_size,
        seed=args.seed,
    )
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
