"""
RAFFLE GNN Fingerprint Benchmark
================================

Benchmarks the Fortran-backed RAFFLE GNN fingerprint pipeline exposed through
the Python package. The example exercises training, inference, and inverse
design against the ATHENA-based Fortran implementation.

Metrics reported:
    - Training time
    - Inference time
    - Descriptor accuracy (MSE vs true RAFFLE fingerprint)
    - Inverse design convergence (loss reduction)

Usage:
        python gnn_fingerprint_benchmark.py

Requirements:
        - raffle (pip install -e .)
        - ase, numpy
"""
import time
import numpy as np
from pathlib import Path
from ase import Atoms
from ase.build import bulk
from ase.visualize import view


def assert_finite(label: str, values) -> np.ndarray:
    array = np.asarray(values, dtype=np.float64)
    if not np.all(np.isfinite(array)):
        raise RuntimeError(f"Non-finite values detected in {label}.")
    return array


def run_benchmark(
    use_mlip: bool,
    label: str,
    dataset_label: str,
    training_structures: list,
    carbon_diamond: Atoms,
    carbon_perturbed: Atoms,
    target_fp: np.ndarray,
    compute_fp,
    num_epochs: int = 20,
    batch_size: int = 8,
    inv_steps: int = 200,
    inv_step_size: float = 50.0,
    use_simple_fingerprint: bool = False,
):
    """Run full pipeline for one Fortran-backed architecture and return metrics."""
    from raffle import GNNFingerprint

    print(f"\n{'='*60}")
    print(f"  Architecture: {label}")
    print(f"  Dataset: {dataset_label} ({len(training_structures)} structures)")
    print(f"  use_mlip_layer = {use_mlip}")
    print(f"{'='*60}")

    gnn = GNNFingerprint(
        species_list=["C"],
        bond_cutoff=6.0,
        gnn_hidden_sizes=[16],
        learning_rate=0.001,
        use_mlip_layer=use_mlip,
    )

    # ---- Training ----
    print(
        f"\n  Training on {len(training_structures)} structures "
        f"for {num_epochs} epochs with batch_size={batch_size} ..."
    )
    t0 = time.perf_counter()
    loss_history = gnn.train(
        structures=training_structures,
        num_epochs=num_epochs,
        batch_size=batch_size,
        verbose=1,
        use_simple_fingerprint=use_simple_fingerprint,
    )
    loss_history = assert_finite(f"{label} loss history", loss_history)
    train_time = time.perf_counter() - t0
    print(f"  Training time: {train_time:.2f} s")
    print(f"  Loss estimate: {loss_history[0]:.6f} -> {loss_history[-1]:.6f}")
    if abs(loss_history[-1] - loss_history[0]) < 1e-9:
        raise RuntimeError(f"{label} training loss was flat.")

    # ---- Inference ----
    print("\n  Inference ...")
    t0 = time.perf_counter()
    n_infer = 10
    for _ in range(n_infer):
        predicted_fp = gnn.predict(carbon_diamond)
    predicted_fp = assert_finite(f"{label} predicted fingerprint", predicted_fp)
    infer_time = (time.perf_counter() - t0) / n_infer
    print(f"  Inference time (avg over {n_infer}): {infer_time*1e3:.2f} ms")

    true_fp = target_fp
    prediction_mse = float(np.mean((predicted_fp - true_fp) ** 2))
    print(f"  Prediction MSE: {prediction_mse:.8f}")

    grad_valid = None
    grad_norm = None
    has_nan = None

    # ---- Inverse design ----
    print(f"\n  Inverse design ({inv_steps} steps) ...")
    test_structure = carbon_perturbed.copy()
    n_atoms = len(test_structure)
    fixed_atoms = np.zeros(n_atoms, dtype=bool)

    t0 = time.perf_counter()
    optimised = gnn.inverse_design(
        target_fingerprint=target_fp,
        atoms=test_structure,
        fixed_atoms=fixed_atoms,
        num_steps=inv_steps,
        step_size=inv_step_size,
        verbose=1,
        use_simple_fingerprint=use_simple_fingerprint,
    )
    inv_time = time.perf_counter() - t0
    print(f"  Inverse design time: {inv_time:.2f} s")

    # Evaluate final descriptor distance
    final_fp = assert_finite(f"{label} optimised fingerprint", compute_fp(optimised))
    final_mse = float(np.mean((final_fp - target_fp) ** 2))
    initial_fp = assert_finite("initial perturbed fingerprint", compute_fp(carbon_perturbed))
    initial_mse = float(np.mean((initial_fp - target_fp) ** 2))
    print(f"  Initial descriptor MSE: {initial_mse:.8f}")
    print(f"  Final descriptor MSE:   {final_mse:.8f}")
    improvement = (initial_mse - final_mse) / max(initial_mse, 1e-12) * 100
    print(f"  Improvement: {improvement:.1f}%")

    view(optimised)
    return {
        "label": label,
        "dataset_label": dataset_label,
        "dataset_size": len(training_structures),
        "train_time": train_time,
        "infer_time": infer_time,
        "prediction_mse": prediction_mse,
        "loss_history": loss_history,
        "inv_time": inv_time,
        "initial_mse": initial_mse,
        "final_mse": final_mse,
        "improvement_pct": improvement,
        "grad_valid": grad_valid,
        "grad_norm": grad_norm,
        "has_nan": has_nan,
        "gnn": gnn,
        "optimised": optimised,
    }


def main():
    from raffle import GNNFingerprint

    print("=" * 60)
    print("RAFFLE GNN Benchmark: Fortran Duvenaud Layer")
    print("=" * 60)

    # -----------------------------------------------------------------
    # Prepare structures
    # -----------------------------------------------------------------
    carbon_diamond = bulk("C", "diamond", a=3.567) * (2, 1, 1)
    carbon_diamond.info["energy"] = -72.0

    carbon_perturbed = carbon_diamond.copy()
    rng = np.random.default_rng(42)
    positions = carbon_perturbed.get_positions()
    positions += rng.normal(0, 0.15, positions.shape)
    carbon_perturbed.set_positions(positions)
    carbon_perturbed.info["energy"] = -71.5

    use_simple_fingerprint = False

    view(carbon_perturbed)
    wait = input("\nPress Enter to continue with benchmarks...")

    # # Keep the Python example aligned with the stable native Fortran test setup.
    # training_structures = [carbon_diamond, carbon_perturbed]
    from ase.io import read
    data_path = Path(__file__).resolve().parents[1] / 'data' / 'carbon.xyz'
    training_structures = read(data_path, index=":")

    # Compute target fingerprint
    tmp_gnn = GNNFingerprint(species_list=["C"], bond_cutoff=6.0)
    if use_simple_fingerprint:
        compute_fp = tmp_gnn.compute_fingerprint_direct
    else:
        compute_fp = tmp_gnn.compute_fingerprint
    target_fp = compute_fp(carbon_diamond)

    # -----------------------------------------------------------------
    # Run benchmarks
    # -----------------------------------------------------------------
    dataset_cases = [
        ("small", training_structures[:8], 1, 8),
        ("large", training_structures[: min(len(training_structures), 64)], 8, 5),
        ("full", training_structures, 8, 3),
    ]

    results = []
    for dataset_label, dataset_structures, batch_size, num_epochs in dataset_cases:
        common_kw = dict(
            dataset_label=dataset_label,
            training_structures=dataset_structures,
            carbon_diamond=carbon_diamond,
            carbon_perturbed=carbon_perturbed,
            target_fp=target_fp,
            compute_fp=compute_fp,
            num_epochs=num_epochs,
            batch_size=batch_size,
            inv_steps=5,
            inv_step_size=0.005,
            use_simple_fingerprint=use_simple_fingerprint,
        )
        results.append(
            run_benchmark(
                use_mlip=False,
                label="Fortran Duvenaud GNN",
                **common_kw,
            )
        )
        results.append(
            run_benchmark(
                use_mlip=True,
                label="Fortran custom GNN",
                **common_kw,
            )
        )

    # -----------------------------------------------------------------
    # Summary table
    # -----------------------------------------------------------------
    print("\n")
    print("=" * 70)
    print("                      BENCHMARK SUMMARY")
    print("=" * 70)
    fmt = "  {:<16s} {:<24s} {:>10s} {:>12s} {:>12s} {:>12s}"
    print(fmt.format("Dataset", "Architecture", "Train (s)", "Infer (ms)", "Loss", "Pred MSE"))
    print("-" * 70)
    for result in results:
        print(
            fmt.format(
                result["dataset_label"],
                result["label"],
                f"{result['train_time']:.2f}",
                f"{result['infer_time'] * 1e3:.2f}",
                f"{result['loss_history'][-1]:.6f}",
                f"{result['prediction_mse']:.6f}",
            )
        )
    print("=" * 70)

    # -----------------------------------------------------------------
    # Plot comparison
    # -----------------------------------------------------------------
    try:
        import matplotlib.pyplot as plt

        plot_result = results[-1]
        fig, axes = plt.subplots(1, 2, figsize=(12, 4))

        # (a) Training loss
        ax = axes[0]
        ax.plot(plot_result["loss_history"], marker="o", label=plot_result["label"])
        ax.set_xlabel("Training checkpoint")
        ax.set_ylabel("Loss")
        ax.set_title("Pre/Post Training Loss")
        ax.set_yscale("log")
        ax.legend()
        ax.grid(True, alpha=0.3)

        # (b) Fingerprint comparison
        ax = axes[1]
        ax.plot(target_fp, label="Target", color="black", linewidth=1.5)
        pred_fp = plot_result["gnn"].predict(carbon_diamond)
        ax.plot(pred_fp, label="Prediction", linestyle="--", alpha=0.8)
        ax.set_xlabel("Fingerprint index")
        ax.set_ylabel("Value")
        ax.set_title("Predicted vs Target Fingerprint")
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.savefig("benchmark_comparison.png", dpi=150)
        print("\n  Plot saved to benchmark_comparison.png")
        if "agg" not in plt.get_backend().lower():
            plt.show()

    except ImportError:
        print("\n  matplotlib not available — skipping plots.")


if __name__ == "__main__":
    main()
