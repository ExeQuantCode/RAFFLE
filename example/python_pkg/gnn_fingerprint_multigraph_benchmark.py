"""Benchmark the multigraph Fortran GNN fingerprint workflow on carbon data.

This script exercises the exact workflow required by plan_new.md:
  1. load example/data/carbon.xyz structures
  2. compute analytical RAFFLE fingerprints
  3. train the Fortran multigraph operator against those fingerprints
  4. evaluate on a perturbed diamond-carbon structure
  5. run inverse design with partial atom freezing
  6. visualise training and optimisation convergence

All model logic stays in the compiled Fortran implementation. Python is used
only for orchestration, data I/O, metrics, and visualisation.
"""

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

from raffle.gnn_fingerprint import GNNFingerprint


def dataset_mse(model: GNNFingerprint, structures) -> float:
    losses = []
    for atoms in structures:
        analytical = model.compute_fingerprint(atoms)
        predicted = model.predict(atoms)
        losses.append(float(np.mean((predicted - analytical) ** 2)))
    return float(np.mean(losses)) if losses else 0.0


def fingerprint_mse(model: GNNFingerprint, atoms, target: np.ndarray) -> float:
    predicted = model.predict(atoms)
    return float(np.mean((predicted - target) ** 2))


def build_perturbed_diamond() -> tuple:
    original = bulk("C", "diamond", a=3.567, cubic=True)
    original.pbc = True
    perturbed = original.copy()
    perturbation = np.array(
        [
            [0.020, -0.015, 0.000],
            [-0.010, 0.010, 0.015],
            [0.012, 0.000, -0.018],
            [-0.014, -0.010, 0.010],
            [0.018, 0.015, -0.010],
            [-0.020, 0.012, 0.000],
            [0.010, -0.018, 0.016],
            [0.000, 0.015, -0.012],
        ],
        dtype=np.float32,
    )
    perturbed.set_positions(perturbed.get_positions() + perturbation)
    return original, perturbed


def rmsd(reference, candidate) -> float:
    ref = reference.copy()
    cand = candidate.copy()
    ref.wrap()
    cand.wrap()
    delta = cand.get_positions() - ref.get_positions()
    return float(np.sqrt(np.mean(np.sum(delta ** 2, axis=1))))


def run_benchmark(
    repo_root: Path,
    epochs: int,
    inverse_steps: int,
    batch_size: int,
    learning_rate: float,
    output_dir: Path,
) -> dict:
    carbon_xyz = repo_root / "example" / "data" / "carbon.xyz"
    structures = read(str(carbon_xyz), index=":")
    if len(structures) == 0:
        raise RuntimeError(f"No structures found in {carbon_xyz}")

    model = GNNFingerprint(
        species_list=["C"],
        bond_cutoff=6.0,
        gnn_hidden_sizes=[32],
        learning_rate=learning_rate,
        lr_decay_rate=1.0e-2,
        num_time_steps=2,
        gnn_output_dim=16,
        max_degree=8,
        layer_type=1,
        n_rbf=12,
        kernel_hidden=32,
        seed=42,
    )

    analytical_fingerprints = [model.compute_fingerprint(atoms) for atoms in structures]
    analytical_component_dims = model.component_dims

    training_history = [dataset_mse(model, structures)]
    for _ in range(epochs):
        model.train(structures, num_epochs=1, batch_size=batch_size, verbose=0)
        training_history.append(dataset_mse(model, structures))

    original, perturbed = build_perturbed_diamond()
    target_fingerprint = model.compute_fingerprint(original)
    analytical_eval = model.compute_fingerprint(perturbed)
    predicted_before = model.predict(perturbed)
    eval_mse_before = float(np.mean((predicted_before - analytical_eval) ** 2))

    grad2, grad3, grad4 = model.compute_gradients(perturbed)
    gradient_norms = {
        "2body": float(np.linalg.norm(grad2)),
        "3body": float(np.linalg.norm(grad3)),
        "4body": float(np.linalg.norm(grad4)),
    }

    fixed_atoms = np.zeros(len(perturbed), dtype=bool)
    fixed_atoms[:4] = True

    optimisation_history = [
        {
            "step": 0,
            "fingerprint_mse": fingerprint_mse(model, perturbed, target_fingerprint),
            "rmsd": rmsd(original, perturbed),
        }
    ]
    optimised = perturbed.copy()
    for step in range(1, inverse_steps + 1):
        optimised = model.inverse_design(
            target_fingerprint=target_fingerprint,
            atoms=optimised,
            fixed_atoms=fixed_atoms,
            num_steps=1,
            step_size=2.0e-3,
            verbose=0,
            use_predict=True,
        )
        optimisation_history.append(
            {
                "step": step,
                "fingerprint_mse": fingerprint_mse(model, optimised, target_fingerprint),
                "rmsd": rmsd(original, optimised),
            }
        )

    predicted_after = model.predict(optimised)
    eval_mse_after = float(np.mean((predicted_after - target_fingerprint) ** 2))
    fp2, fp3, fp4 = model.compute_fingerprint_components(optimised)

    output_dir.mkdir(parents=True, exist_ok=True)
    write(output_dir / "diamond_original.xyz", original)
    write(output_dir / "diamond_perturbed.xyz", perturbed)
    write(output_dir / "diamond_optimised.xyz", optimised)

    fig = plt.figure(figsize=(14, 4.5))

    ax1 = fig.add_subplot(1, 3, 1)
    ax1.plot(range(len(training_history)), training_history, marker="o")
    ax1.set_xlabel("Epoch")
    ax1.set_ylabel("Dataset prediction MSE")
    ax1.set_title("Training convergence")

    ax2 = fig.add_subplot(1, 3, 2)
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
    ax2.set_title("Optimisation convergence")
    ax2.legend(loc="best")

    ax3 = fig.add_subplot(1, 3, 3, projection="3d")
    orig_pos = original.get_positions()
    opt_pos = optimised.get_positions()
    ax3.scatter(orig_pos[:, 0], orig_pos[:, 1], orig_pos[:, 2], s=35, label="original")
    ax3.scatter(
        opt_pos[:, 0],
        opt_pos[:, 1],
        opt_pos[:, 2],
        s=35,
        marker="^",
        label="optimised",
    )
    ax3.set_title("Structure comparison")
    ax3.legend(loc="best")

    fig.tight_layout()
    figure_path = output_dir / "gnn_multigraph_benchmark.png"
    fig.savefig(figure_path, dpi=150)
    plt.close(fig)

    metrics = {
        "num_structures": len(structures),
        "fingerprint_dim": model.fingerprint_dim,
        "component_dims": {
            "2body": analytical_component_dims[0],
            "3body": analytical_component_dims[1],
            "4body": analytical_component_dims[2],
        },
        "training_history": training_history,
        "eval_mse_before": eval_mse_before,
        "eval_mse_after": eval_mse_after,
        "inverse_design_history": optimisation_history,
        "gradient_norms": gradient_norms,
        "analytical_fingerprint_norm_mean": float(
            np.mean([np.linalg.norm(fp) for fp in analytical_fingerprints])
        ),
        "optimised_component_norms": {
            "2body": float(np.linalg.norm(fp2)),
            "3body": float(np.linalg.norm(fp3)),
            "4body": float(np.linalg.norm(fp4)),
        },
        "output_files": {
            "plot": str(figure_path),
            "original": str(output_dir / "diamond_original.xyz"),
            "perturbed": str(output_dir / "diamond_perturbed.xyz"),
            "optimised": str(output_dir / "diamond_optimised.xyz"),
        },
    }

    metrics_path = output_dir / "gnn_multigraph_benchmark_metrics.json"
    metrics_path.write_text(json.dumps(metrics, indent=2))
    return metrics


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--epochs", type=int, default=8)
    parser.add_argument("--inverse-steps", type=int, default=8)
    parser.add_argument("--batch-size", type=int, default=2)
    parser.add_argument("--learning-rate", type=float, default=5.0e-4)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("build") / "gnn_multigraph_benchmark",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    repo_root = Path(__file__).resolve().parents[2]
    metrics = run_benchmark(
        repo_root=repo_root,
        epochs=args.epochs,
        inverse_steps=args.inverse_steps,
        batch_size=args.batch_size,
        learning_rate=args.learning_rate,
        output_dir=args.output_dir,
    )
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
