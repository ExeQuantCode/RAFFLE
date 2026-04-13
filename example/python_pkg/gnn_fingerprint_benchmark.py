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
from itertools import permutations
from ase import Atoms
from ase.build import bulk
from ase.visualize import view
from ase.geometry import find_mic


def assert_finite(label: str, values) -> np.ndarray:
    array = np.asarray(values, dtype=np.float64)
    if not np.all(np.isfinite(array)):
        raise RuntimeError(f"Non-finite values detected in {label}.")
    return array


def position_rmsd(atoms_a: Atoms, atoms_b: Atoms) -> float:
    """Compute per-atom RMSD using minimum image convention and best atom permutation."""
    n = len(atoms_a)
    best_rmsd = float("inf")
    for perm in permutations(range(n)):
        diff = atoms_a.get_positions()[list(perm)] - atoms_b.get_positions()
        _, mic_dist = find_mic(diff, atoms_a.cell)
        rmsd = float(np.sqrt(np.mean(mic_dist ** 2)))
        if rmsd < best_rmsd:
            best_rmsd = rmsd
    return best_rmsd


def max_position_error(atoms_a: Atoms, atoms_b: Atoms) -> float:
    """Compute max per-atom error using minimum image convention and best atom permutation."""
    n = len(atoms_a)
    best_max_err = float("inf")
    for perm in permutations(range(n)):
        diff = atoms_a.get_positions()[list(perm)] - atoms_b.get_positions()
        _, mic_dist = find_mic(diff, atoms_a.cell)
        max_err = float(np.max(mic_dist))
        if max_err < best_max_err:
            best_max_err = max_err
    return best_max_err


def run_benchmark(
    use_mlip: bool,
    label: str,
    dataset_label: str,
    training_structures: list,
    carbon_diamond: Atoms,
    carbon_perturbed: Atoms,
    target_fp: np.ndarray,
    num_epochs: int = 20,
    inv_steps: int = 200,
    inv_step_size: float = 50.0,
    use_simple_fingerprint: bool = False,
    layer_type: int = -1,
    learning_rate: float = 0.01,
    lr_decay_rate: float = 1e-2,
    seed: int = 42,
):
    """Run full pipeline for one Fortran-backed architecture and return metrics."""
    from raffle import GNNFingerprint

    print(f"\n{'='*60}")
    print(f"  Architecture: {label}")
    print(f"  Dataset: {dataset_label} ({len(training_structures)} structures)")
    if layer_type >= 0:
        print(f"  layer_type = {layer_type}")
    else:
        print(f"  use_mlip_layer = {use_mlip}")
    print(f"{'='*60}")

    gnn = GNNFingerprint(
        species_list=["C"],
        bond_cutoff=6.0,
        gnn_hidden_sizes=[16],
        learning_rate=learning_rate,
        lr_decay_rate=lr_decay_rate,
        use_mlip_layer=use_mlip,
        layer_type=layer_type,
        seed = seed,
    )

    # Use the GNN's own compute_fingerprint method for evaluation
    if use_simple_fingerprint:
        compute_fp = gnn.compute_fingerprint_direct
    else:
        compute_fp = gnn.compute_fingerprint

    # ---- Untrained inverse design (baseline) ----
    # Use untrained predict() for both target and loop so the comparison
    # shows the effect of training on the GNN landscape.
    untrained_pred_target = assert_finite(
        f"{label} untrained predict(diamond)", gnn.predict(carbon_diamond)
    )
    print(f"\n  Untrained inverse design ({inv_steps} steps, use_predict) ...")
    test_structure_untrained = carbon_perturbed.copy()
    n_atoms = len(test_structure_untrained)
    fixed_atoms = np.zeros(n_atoms, dtype=bool)

    untrained_optimised = gnn.inverse_design(
        target_fingerprint=untrained_pred_target,
        atoms=test_structure_untrained,
        fixed_atoms=fixed_atoms,
        num_steps=inv_steps,
        step_size=inv_step_size,
        verbose=0,
        use_simple_fingerprint=use_simple_fingerprint,
        use_predict=True,
    )
    untrained_rmsd = position_rmsd(untrained_optimised, carbon_diamond)
    untrained_max_err = max_position_error(untrained_optimised, carbon_diamond)
    print(f"  Untrained RMSD vs diamond: {untrained_rmsd:.4f} Å")
    print(f"  Untrained max error: {untrained_max_err:.4f} Å")

    # ---- Training ----
    print(
        f"\n  Training on {len(training_structures)} structures "
        f"for {num_epochs} epochs ..."
    )
    t0 = time.perf_counter()
    loss_history = gnn.train(
        structures=training_structures,
        num_epochs=num_epochs,
        verbose=0,
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

    # ---- Inverse design (trained, using GNN predict) ----
    trained_pred_target = assert_finite(
        f"{label} trained predict(diamond)", gnn.predict(carbon_diamond)
    )
    print(f"\n  Trained inverse design ({inv_steps} steps, use_predict) ...")
    test_structure = carbon_perturbed.copy()
    n_atoms = len(test_structure)
    fixed_atoms = np.zeros(n_atoms, dtype=bool)

    t0 = time.perf_counter()
    optimised = gnn.inverse_design(
        target_fingerprint=trained_pred_target,
        atoms=test_structure,
        fixed_atoms=fixed_atoms,
        num_steps=inv_steps,
        step_size=inv_step_size,
        verbose=1,
        use_simple_fingerprint=use_simple_fingerprint,
        use_predict=True,
    )
    inv_time = time.perf_counter() - t0
    print(f"  Inverse design time: {inv_time:.2f} s")

    # Evaluate final descriptor distance (using analytical descriptor for consistency)
    final_fp = assert_finite(f"{label} optimised fingerprint", compute_fp(optimised))
    final_mse = float(np.mean((final_fp - target_fp) ** 2))
    initial_fp = assert_finite("initial perturbed fingerprint", compute_fp(carbon_perturbed))
    initial_mse = float(np.mean((initial_fp - target_fp) ** 2))
    print(f"  Initial descriptor MSE: {initial_mse:.8f}")
    print(f"  Final descriptor MSE:   {final_mse:.8f}")
    improvement = (initial_mse - final_mse) / max(initial_mse, 1e-12) * 100
    print(f"  Improvement: {improvement:.1f}%")

    # Position comparison with diamond
    trained_rmsd = position_rmsd(optimised, carbon_diamond)
    trained_max_err = max_position_error(optimised, carbon_diamond)
    perturbed_rmsd = position_rmsd(carbon_perturbed, carbon_diamond)
    print(f"\n  Position comparison vs diamond:")
    print(f"    Perturbed RMSD:  {perturbed_rmsd:.4f} Å")
    print(f"    Untrained RMSD:  {untrained_rmsd:.4f} Å")
    print(f"    Trained RMSD:    {trained_rmsd:.4f} Å")
    print(f"    Trained max err: {trained_max_err:.4f} Å")
    diamond_reproduced = trained_max_err < 0.05
    trained_better = trained_rmsd < untrained_rmsd
    print(f"    Diamond reproduced (<0.05 Å): {'YES' if diamond_reproduced else 'NO'}")
    print(f"    Trained better than untrained: {'YES' if trained_better else 'NO'}")

    view(optimised)
    del gnn

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
        "optimised": optimised,
        "untrained_rmsd": untrained_rmsd,
        "trained_rmsd": trained_rmsd,
        "trained_max_err": trained_max_err,
        "diamond_reproduced": diamond_reproduced,
        "trained_better": trained_better,
    }


def main():
    from raffle import GNNFingerprint

    seed = 42
    # initialise random seed
    np.random.seed(seed)

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

    # view(carbon_perturbed)
    # wait = input("\nPress Enter to continue with benchmarks...")

    # # Keep the Python example aligned with the stable native Fortran test setup.
    # training_structures = [carbon_diamond, carbon_perturbed]
    from ase.io import read
    data_path = Path(__file__).resolve().parents[1] / 'data' / 'carbon.xyz'
    training_structures = read(data_path, index=":")

    # Compute target fingerprint
    tmp_gnn = GNNFingerprint(species_list=["C"], bond_cutoff=6.0)
    if use_simple_fingerprint:
        target_fp = tmp_gnn.compute_fingerprint_direct(carbon_diamond)
    else:
        target_fp = tmp_gnn.compute_fingerprint(carbon_diamond)
    del tmp_gnn

    # -----------------------------------------------------------------
    # Run benchmarks
    # -----------------------------------------------------------------
    layer_configs = [
        (0, "Duvenaud", 0.01),
        # (1, "RAFFLE MLIP", 0.01),
        # # (2, "SchNet", 0.001),
        # # (3, "DimeNet", 0.001),
        # (4, "Hybrid", 0.01),
    ]

    # Use a single dataset configuration for quick comparison
    dataset_structures = training_structures#[:min(len(training_structures), 32)]
    num_epochs = 50

    results = []
    for layer_id, layer_name, lr in layer_configs:
        common_kw = dict(
            dataset_label="large",
            training_structures=dataset_structures,
            carbon_diamond=carbon_diamond,
            carbon_perturbed=carbon_perturbed,
            target_fp=target_fp,
            num_epochs=num_epochs,
            inv_steps=500,
            inv_step_size=1.0,
            use_simple_fingerprint=use_simple_fingerprint,
            learning_rate=lr,
        )
        results.append(
            run_benchmark(
                use_mlip=False,
                label=f"Fortran {layer_name}",
                layer_type=layer_id,
                seed=seed,
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
    fmt = "  {:<16s} {:<24s} {:>10s} {:>12s} {:>10s} {:>10s} {:>8s} {:>8s}"
    print(fmt.format("Dataset", "Architecture", "Train (s)", "Pred MSE", "RMSD (Å)", "MaxErr(Å)", "Repro?", "Better?"))
    print("-" * 100)
    for result in results:
        print(
            fmt.format(
                result["dataset_label"],
                result["label"],
                f"{result['train_time']:.2f}",
                f"{result['prediction_mse']:.6f}",
                f"{result['trained_rmsd']:.4f}",
                f"{result['trained_max_err']:.4f}",
                "YES" if result["diamond_reproduced"] else "NO",
                "YES" if result["trained_better"] else "NO",
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

        # (b) Descriptor improvement per layer
        ax = axes[1]
        labels = [r["label"] for r in results]
        untrained = [r["untrained_rmsd"] for r in results]
        trained = [r["trained_rmsd"] for r in results]
        x = np.arange(len(labels))
        ax.bar(x - 0.2, untrained, 0.35, label="Untrained", alpha=0.7)
        ax.bar(x + 0.2, trained, 0.35, label="Trained", alpha=0.7)
        ax.set_xticks(x)
        ax.set_xticklabels([l.replace("Fortran ", "") for l in labels], rotation=30, ha="right")
        ax.set_ylabel("RMSD vs diamond (Å)")
        ax.set_title("Inverse Design: Untrained vs Trained")
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3, axis="y")

        plt.tight_layout()
        plt.savefig("benchmark_comparison.png", dpi=150)
        print("\n  Plot saved to benchmark_comparison.png")
        if "agg" not in plt.get_backend().lower():
            plt.show()

    except ImportError:
        print("\n  matplotlib not available — skipping plots.")


if __name__ == "__main__":
    main()
