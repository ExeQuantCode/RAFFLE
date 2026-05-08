"""PyTorch carbon workflow example mirroring the multigraph test script."""

from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path
from typing import Callable, Optional

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from ase.build import bulk
from ase.io import read, write

from raffle import (
    TorchGNNFingerprint,
    minimum_image_displacements,
    symmetry_aware_displacements,
    symmetry_aware_rmsd,
)


PERTURBATION = np.array(
    [
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.172792, 0.2410809, 0.165219],
        [-0.00651579, 0.00452678, 0.20223187],
        [-0.10268477, 0.10290559, 0.20182286],
        [2.00147066, 0.00014211, 0.10273356],
    ],
    dtype=np.float32,
)

FIXED_LEADING_ATOMS = 0
FINGERPRINT_LOSS_WEIGHT = 0.5
TARGET_VERTEX_WEIGHT = 0.5
TARGET_POSITION_WEIGHT = 0.0
INVERSE_LR_DECAY_RATE = 0.0
INVERSE_RESTARTS = 1
INVERSE_RESTART_NOISE_SCALE = 0.0
DEFAULT_COMPONENT_WEIGHT = (4.0, 1.0, 1.0)
DEFAULT_MODEL_CONFIG = {
    "architecture": "residual",
    "hidden_dim": 80,
    "num_message_layers": 2,
    "learning_rate": 5.0e-4,
    "lr_decay_rate": 5.0e-3,
    "smooth_cutoff_width": 0.2,
    "reference_layer_type": 1,
    "component_weight": DEFAULT_COMPONENT_WEIGHT,
}


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


def build_perturbed_structure(original, fixed_leading_atoms: int = FIXED_LEADING_ATOMS):
    perturbed = original.copy()
    perturbed.set_positions(perturbed.get_positions() + PERTURBATION)
    fixed_atoms = np.zeros(len(perturbed), dtype=bool)
    fixed_atoms[:max(int(fixed_leading_atoms), 0)] = True
    return perturbed, fixed_atoms


def select_carbon_structures(structures: list, requested_count: int):
    if int(requested_count) <= 0 or int(requested_count) >= len(structures):
        return list(structures)
    return list(structures[: int(requested_count)])


def create_model(seed: int, model_config: Optional[dict] = None) -> TorchGNNFingerprint:
    if TorchGNNFingerprint is None:
        raise RuntimeError("TorchGNNFingerprint is unavailable in this installation")
    resolved_config = dict(DEFAULT_MODEL_CONFIG)
    if model_config:
        resolved_config.update(model_config)
    return TorchGNNFingerprint(
        species_list=["C"],
        seed=seed,
        **resolved_config,
    )


def default_epoch_values(max_epochs: int) -> list[int]:
    candidates = [0, 1, 5, 10, 20, max_epochs]
    return sorted({value for value in candidates if 0 <= value <= max_epochs})


def default_inverse_step_values(max_steps: int) -> list[int]:
    candidates = [0, 10, 25, 50, 100, 200, max_steps]
    return sorted({value for value in candidates if 0 <= value <= max_steps})


def default_step_sizes(base_step_size: float) -> list[float]:
    candidates = [base_step_size * factor for factor in (0.25, 0.5, 1.0, 1.5, 2.0)]
    return sorted({float(f"{value:.8g}") for value in candidates if value > 0.0})


def parse_int_list(value: str) -> list[int]:
    return [int(item.strip()) for item in value.split(",") if item.strip()]


def parse_float_list(value: str) -> list[float]:
    return [float(item.strip()) for item in value.split(",") if item.strip()]


def build_inverse_design_options(
    original,
    fingerprint_loss_weight: float,
    target_vertex_weight: float,
    target_position_weight: float,
    inverse_lr_decay_rate: float,
    inverse_restarts: int,
    inverse_restart_noise_scale: float,
) -> dict:
    return {
        "target_atoms": original,
        "fingerprint_loss_weight": float(fingerprint_loss_weight),
        "target_vertex_weight": float(target_vertex_weight),
        "target_position_weight": float(target_position_weight),
        "inverse_lr_decay_rate": float(inverse_lr_decay_rate),
        "num_restarts": int(inverse_restarts),
        "restart_noise_scale": float(inverse_restart_noise_scale),
    }


def score_candidate(original, candidate_atoms) -> float:
    return float(symmetry_aware_rmsd(original, candidate_atoms))


def update_best_candidate(
    best_candidate,
    candidate_atoms,
    position_difference: float,
    source: str,
    **metadata,
):
    candidate = {
        "atoms": candidate_atoms.copy(),
        "source": source,
        "position_difference": float(position_difference),
        **metadata,
    }
    if best_candidate is None or candidate["position_difference"] < best_candidate["position_difference"]:
        return candidate
    return best_candidate


def capture_checkpoint(model: TorchGNNFingerprint, epochs: int, position_difference: float) -> dict:
    return {
        "epochs": int(epochs),
        "position_difference": float(position_difference),
        "state_dict": copy.deepcopy(model.state_dict()),
    }


def inverse_design_trace(
    model: TorchGNNFingerprint,
    original,
    perturbed,
    fixed_atoms: np.ndarray,
    target_fingerprint: np.ndarray,
    step_values: list[int],
    step_size: float,
    inverse_design_options: dict,
) -> tuple:
    requested_step_values = sorted({0, *[max(int(step), 0) for step in step_values]})
    trace = []
    optimised = perturbed.copy()
    best_candidate = None
    for step in requested_step_values:
        if step == 0:
            candidate_atoms = perturbed.copy()
        else:
            candidate_atoms = model.inverse_design(
                target_fingerprint=target_fingerprint,
                atoms=perturbed,
                fixed_atoms=fixed_atoms,
                num_steps=step,
                step_size=step_size,
                verbose=0,
                **inverse_design_options,
            )
        position_difference = score_candidate(original, candidate_atoms)
        if step == requested_step_values[-1]:
            optimised = candidate_atoms.copy()
        best_candidate = update_best_candidate(
            best_candidate,
            candidate_atoms,
            position_difference,
            source="inverse_step_sweep",
            step=int(step),
        )
        trace.append(
            {
                "step": step,
                "position_difference": position_difference,
            }
        )
    return optimised, trace, best_candidate


def sweep_epochs(
    carbon_structures,
    augmented_structures,
    original,
    perturbed,
    fixed_atoms: np.ndarray,
    epoch_values: list[int],
    batch_size: int,
    inverse_steps: int,
    inverse_step_size: float,
    inverse_design_options: dict,
    seed: int,
    model_config: Optional[dict] = None,
    training_observer: Optional[Callable[[int, float], None]] = None,
) -> tuple:
    model = create_model(seed, model_config=model_config)
    target_fingerprint = model.compute_reference_fingerprint(original)
    epoch_results = []
    training_losses = []
    best_candidate = None
    checkpoint_records = []

    requested_epochs = sorted(set(epoch_values))
    trained_epochs = 0
    if 0 in requested_epochs:
        untrained_optimised = model.inverse_design(
            target_fingerprint=target_fingerprint,
            atoms=perturbed,
            fixed_atoms=fixed_atoms,
            num_steps=inverse_steps,
            step_size=inverse_step_size,
            verbose=0,
            **inverse_design_options,
        )
        position_difference = score_candidate(original, untrained_optimised)
        best_candidate = update_best_candidate(
            best_candidate,
            untrained_optimised,
            position_difference,
            source="epoch_sweep",
            epochs=0,
        )
        checkpoint_records.append(capture_checkpoint(model, 0, position_difference))
        epoch_results.append(
            {
                "epochs": 0,
                "position_difference": position_difference,
            }
        )

    for target_epoch in requested_epochs:
        if target_epoch == 0:
            continue
        while trained_epochs < target_epoch:
            history = model.fit(
                carbon_structures,
                num_epochs=1,
                batch_size=batch_size,
                augment_structures=augmented_structures,
                verbose=0,
            )
            trained_epochs += 1
            training_losses.append(float(history[-1]))
            if training_observer is not None:
                training_observer(trained_epochs, float(history[-1]))

        optimised = model.inverse_design(
            target_fingerprint=target_fingerprint,
            atoms=perturbed,
            fixed_atoms=fixed_atoms,
            num_steps=inverse_steps,
            step_size=inverse_step_size,
            verbose=0,
            **inverse_design_options,
        )
        position_difference = score_candidate(original, optimised)
        best_candidate = update_best_candidate(
            best_candidate,
            optimised,
            position_difference,
            source="epoch_sweep",
            epochs=int(target_epoch),
        )
        checkpoint_records.append(capture_checkpoint(model, int(target_epoch), position_difference))
        epoch_results.append(
            {
                "epochs": target_epoch,
                "position_difference": position_difference,
            }
        )

    return model, epoch_results, training_losses, best_candidate, checkpoint_records


def sweep_step_sizes(
    model: TorchGNNFingerprint,
    original,
    perturbed,
    fixed_atoms: np.ndarray,
    target_fingerprint: np.ndarray,
    inverse_steps: int,
    step_sizes: list[float],
    inverse_design_options: dict,
) -> tuple[list[dict[str, float]], dict | None]:
    results = []
    best_candidate = None
    for step_size in step_sizes:
        optimised = model.inverse_design(
            target_fingerprint=target_fingerprint,
            atoms=perturbed,
            fixed_atoms=fixed_atoms,
            num_steps=inverse_steps,
            step_size=step_size,
            verbose=0,
            **inverse_design_options,
        )
        position_difference = score_candidate(original, optimised)
        best_candidate = update_best_candidate(
            best_candidate,
            optimised,
            position_difference,
            source="step_size_sweep",
            step_size=float(step_size),
            num_steps=int(inverse_steps),
        )
        results.append(
            {
                "step_size": float(step_size),
                "position_difference": position_difference,
            }
        )
    return results, best_candidate


def sweep_step_sizes_for_checkpoints(
    model: TorchGNNFingerprint,
    checkpoint_records: list[dict],
    original,
    perturbed,
    fixed_atoms: np.ndarray,
    target_fingerprint: np.ndarray,
    inverse_steps: int,
    step_sizes: list[float],
    inverse_design_options: dict,
) -> tuple[list[dict[str, float]], dict | None]:
    if not checkpoint_records:
        return [], None

    results = []
    best_candidate = None
    original_state = copy.deepcopy(model.state_dict())
    try:
        for checkpoint in checkpoint_records:
            model.load_state_dict(checkpoint["state_dict"])
            for step_size in step_sizes:
                optimised = model.inverse_design(
                    target_fingerprint=target_fingerprint,
                    atoms=perturbed,
                    fixed_atoms=fixed_atoms,
                    num_steps=inverse_steps,
                    step_size=step_size,
                    verbose=0,
                    **inverse_design_options,
                )
                position_difference = score_candidate(original, optimised)
                results.append(
                    {
                        "epochs": int(checkpoint["epochs"]),
                        "step_size": float(step_size),
                        "position_difference": position_difference,
                    }
                )
                best_candidate = update_best_candidate(
                    best_candidate,
                    optimised,
                    position_difference,
                    source="checkpoint_step_size_sweep",
                    epochs=int(checkpoint["epochs"]),
                    step_size=float(step_size),
                    num_steps=int(inverse_steps),
                )
    finally:
        model.load_state_dict(original_state)

    return results, best_candidate


def sweep_step_schedules_for_checkpoints(
    model: TorchGNNFingerprint,
    checkpoint_records: list[dict],
    original,
    perturbed,
    fixed_atoms: np.ndarray,
    target_fingerprint: np.ndarray,
    step_values: list[int],
    step_sizes: list[float],
    inverse_design_options: dict,
) -> tuple[list[dict[str, float]], dict | None]:
    if not checkpoint_records:
        return [], None

    requested_steps = sorted({max(int(step), 0) for step in step_values if int(step) > 0})
    if not requested_steps:
        return [], None

    results = []
    best_candidate = None
    original_state = copy.deepcopy(model.state_dict())
    try:
        for checkpoint in checkpoint_records:
            model.load_state_dict(checkpoint["state_dict"])
            for step in requested_steps:
                for step_size in step_sizes:
                    optimised = model.inverse_design(
                        target_fingerprint=target_fingerprint,
                        atoms=perturbed,
                        fixed_atoms=fixed_atoms,
                        num_steps=step,
                        step_size=step_size,
                        verbose=0,
                        **inverse_design_options,
                    )
                    position_difference = score_candidate(original, optimised)
                    results.append(
                        {
                            "epochs": int(checkpoint["epochs"]),
                            "num_steps": int(step),
                            "step_size": float(step_size),
                            "position_difference": position_difference,
                        }
                    )
                    best_candidate = update_best_candidate(
                        best_candidate,
                        optimised,
                        position_difference,
                        source="checkpoint_step_schedule_sweep",
                        epochs=int(checkpoint["epochs"]),
                        num_steps=int(step),
                        step_size=float(step_size),
                    )
    finally:
        model.load_state_dict(original_state)

    return results, best_candidate


def print_position_differences(original, perturbed, optimised) -> None:
    initial_delta = symmetry_aware_displacements(original, perturbed)
    final_delta = symmetry_aware_displacements(original, optimised)
    optimisation_delta = minimum_image_displacements(perturbed, optimised)

    print("Initial position differences relative to the target structure (Angstrom):")
    print(np.array2string(initial_delta, precision=6, suppress_small=False))
    print()
    print("Final position differences relative to the target structure (Angstrom):")
    print(np.array2string(final_delta, precision=6, suppress_small=False))
    print()
    print("Net position update from initial to final structure (Angstrom):")
    print(np.array2string(optimisation_delta, precision=6, suppress_small=False))
    print()
    print("Per-atom displacement norms (initial, final, update) in Angstrom:")
    initial_norms = np.linalg.norm(initial_delta, axis=1)
    final_norms = np.linalg.norm(final_delta, axis=1)
    update_norms = np.linalg.norm(optimisation_delta, axis=1)
    for atom_index, (initial_norm, final_norm, update_norm) in enumerate(
        zip(initial_norms, final_norms, update_norms)
    ):
        print(
            f"  atom {atom_index:2d}: initial={initial_norm:.6f} "
            f"final={final_norm:.6f} update={update_norm:.6f}"
        )


def plot_position_sweeps(
    epoch_results: list[dict[str, float]],
    inverse_trace: list[dict[str, float]],
    step_size_results: list[dict[str, float]],
    initial_position_difference: float,
    output_path: Path,
) -> None:
    figure = plt.figure(figsize=(15, 4.5))

    ax1 = figure.add_subplot(1, 3, 1)
    ax1.plot(
        [entry["epochs"] for entry in epoch_results],
        [entry["position_difference"] for entry in epoch_results],
        marker="o",
    )
    ax1.set_xlabel("Epochs trained")
    ax1.set_ylabel("Position difference to target (symmetry-aware RMSD, A)")
    ax1.set_title("Training sweep")
    ax1.axhline(initial_position_difference, color="tab:red", linestyle="--", linewidth=1.0)
    ax1.grid(alpha=0.3)

    ax2 = figure.add_subplot(1, 3, 2)
    ax2.plot(
        [entry["step"] for entry in inverse_trace],
        [entry["position_difference"] for entry in inverse_trace],
        marker="o",
    )
    ax2.set_xlabel("Inverse-design steps")
    ax2.set_ylabel("Position difference to target (symmetry-aware RMSD, A)")
    ax2.set_title("Inverse-design sweep")
    ax2.axhline(initial_position_difference, color="tab:red", linestyle="--", linewidth=1.0)
    ax2.grid(alpha=0.3)

    ax3 = figure.add_subplot(1, 3, 3)
    ax3.plot(
        [entry["step_size"] for entry in step_size_results],
        [entry["position_difference"] for entry in step_size_results],
        marker="o",
    )
    ax3.set_xscale("log")
    ax3.set_xlabel("Step size")
    ax3.set_ylabel("Position difference to target (symmetry-aware RMSD, A)")
    ax3.set_title("Step-size sweep")
    ax3.axhline(initial_position_difference, color="tab:red", linestyle="--", linewidth=1.0)
    ax3.grid(alpha=0.3)

    figure.tight_layout()
    figure.savefig(output_path, dpi=150)
    plt.close(figure)


def run_example(
    repo_root: Path,
    output_dir: Path,
    carbon_count: int,
    augmented_count: int,
    num_epochs: int,
    batch_size: int,
    inverse_steps: int,
    inverse_step_size: float,
    epoch_values: list[int],
    inverse_step_values: list[int],
    step_size_values: list[float],
    fixed_leading_atoms: int,
    fingerprint_loss_weight: float,
    target_vertex_weight: float,
    target_position_weight: float,
    inverse_lr_decay_rate: float,
    inverse_restarts: int,
    inverse_restart_noise_scale: float,
    seed: int,
    model_config: Optional[dict] = None,
    training_observer: Optional[Callable[[int, float], None]] = None,
    architecture_name: str = "torch_gnn_residual",
) -> dict:
    if float(target_position_weight) != 0.0:
        raise ValueError("target_position_weight must remain 0.0 for inverse-design runs")

    carbon_xyz = repo_root / "example" / "data" / "carbon.xyz"
    all_carbon_structures = read(str(carbon_xyz), index=":")
    carbon_structures = select_carbon_structures(all_carbon_structures, carbon_count)

    original = bulk("C", "diamond", a=3.567, cubic=True)
    original.pbc = True
    augmented_structures = build_augmented_structures(original, augmented_count, seed)
    perturbed, fixed_atoms = build_perturbed_structure(original, fixed_leading_atoms=fixed_leading_atoms)
    inverse_design_options = build_inverse_design_options(
        original=original,
        fingerprint_loss_weight=fingerprint_loss_weight,
        target_vertex_weight=target_vertex_weight,
        target_position_weight=target_position_weight,
        inverse_lr_decay_rate=inverse_lr_decay_rate,
        inverse_restarts=inverse_restarts,
        inverse_restart_noise_scale=inverse_restart_noise_scale,
    )

    model, epoch_results, training_losses, best_epoch_candidate, checkpoint_records = sweep_epochs(
        carbon_structures=carbon_structures,
        augmented_structures=augmented_structures,
        original=original,
        perturbed=perturbed,
        fixed_atoms=fixed_atoms,
        epoch_values=epoch_values,
        batch_size=batch_size,
        inverse_steps=inverse_steps,
        inverse_step_size=inverse_step_size,
        inverse_design_options=inverse_design_options,
        seed=seed,
        model_config=model_config,
        training_observer=training_observer,
    )

    target_fingerprint = model.compute_reference_fingerprint(original)
    optimised, full_inverse_trace, best_inverse_candidate = inverse_design_trace(
        model=model,
        original=original,
        perturbed=perturbed,
        fixed_atoms=fixed_atoms,
        target_fingerprint=target_fingerprint,
        step_values=inverse_step_values,
        step_size=inverse_step_size,
        inverse_design_options=inverse_design_options,
    )
    inverse_trace = full_inverse_trace
    step_size_results, best_step_size_candidate = sweep_step_sizes(
        model=model,
        original=original,
        perturbed=perturbed,
        fixed_atoms=fixed_atoms,
        target_fingerprint=target_fingerprint,
        inverse_steps=inverse_steps,
        step_sizes=step_size_values,
        inverse_design_options=inverse_design_options,
    )
    checkpoint_step_size_results, best_checkpoint_step_size_candidate = sweep_step_sizes_for_checkpoints(
        model=model,
        checkpoint_records=checkpoint_records,
        original=original,
        perturbed=perturbed,
        fixed_atoms=fixed_atoms,
        target_fingerprint=target_fingerprint,
        inverse_steps=inverse_steps,
        step_sizes=step_size_values,
        inverse_design_options=inverse_design_options,
    )
    checkpoint_step_schedule_results, best_checkpoint_step_schedule_candidate = (
        sweep_step_schedules_for_checkpoints(
            model=model,
            checkpoint_records=checkpoint_records,
            original=original,
            perturbed=perturbed,
            fixed_atoms=fixed_atoms,
            target_fingerprint=target_fingerprint,
            step_values=inverse_step_values,
            step_sizes=step_size_values,
            inverse_design_options=inverse_design_options,
        )
    )

    configured_final_rmsd = score_candidate(original, optimised)
    best_candidate = update_best_candidate(
        None,
        optimised,
        configured_final_rmsd,
        source="configured_inverse_design",
        epochs=int(num_epochs),
        num_steps=int(inverse_steps),
        step_size=float(inverse_step_size),
    )
    for candidate in (
        best_epoch_candidate,
        best_inverse_candidate,
        best_step_size_candidate,
        best_checkpoint_step_size_candidate,
        best_checkpoint_step_schedule_candidate,
    ):
        if candidate is None:
            continue
        if candidate["position_difference"] < best_candidate["position_difference"]:
            best_candidate = candidate

    optimised = best_candidate["atoms"].copy()

    print_position_differences(original, perturbed, optimised)
    print()
    initial_rmsd = score_candidate(original, perturbed)
    final_rmsd = best_candidate["position_difference"]
    print(f"Initial symmetry-aware RMSD: {initial_rmsd:.6f} A")
    print(f"Final symmetry-aware RMSD:   {final_rmsd:.6f} A")
    print(f"RMSD reduction:              {100.0 * (1.0 - final_rmsd / initial_rmsd):.2f}%")
    if best_candidate["source"] != "configured_inverse_design":
        best_candidate_details = {
            key: value
            for key, value in best_candidate.items()
            if key not in {"atoms", "position_difference"}
        }
        print(f"Selected best candidate: {json.dumps(best_candidate_details, sort_keys=True)}")

    output_dir.mkdir(parents=True, exist_ok=True)
    write(output_dir / "torch_gnn_carbon_original.xyz", original)
    write(output_dir / "torch_gnn_carbon_initial.xyz", perturbed)
    write(output_dir / "torch_gnn_carbon_final.xyz", optimised)

    figure_path = output_dir / "torch_gnn_carbon_position_sweeps.png"
    plot_position_sweeps(
        epoch_results=epoch_results,
        inverse_trace=inverse_trace,
        step_size_results=step_size_results,
        initial_position_difference=initial_rmsd,
        output_path=figure_path,
    )

    metrics = {
        "training_losses": training_losses,
        "num_available_carbon_structures": len(all_carbon_structures),
        "num_training_carbon_structures": len(carbon_structures),
        "epoch_sweep": epoch_results,
        "inverse_step_sweep": inverse_trace,
        "step_size_sweep": step_size_results,
        "checkpoint_step_size_sweep": checkpoint_step_size_results,
        "checkpoint_step_schedule_sweep": checkpoint_step_schedule_results,
        "initial_position_difference": symmetry_aware_displacements(original, perturbed).tolist(),
        "final_position_difference": symmetry_aware_displacements(original, optimised).tolist(),
        "initial_rmsd": initial_rmsd,
        "configured_final_rmsd": configured_final_rmsd,
        "final_rmsd": final_rmsd,
        "rmsd_reduction_fraction": 1.0 - final_rmsd / initial_rmsd,
        "selected_candidate": {
            key: value
            for key, value in best_candidate.items()
            if key != "atoms"
        },
        "position_error_metric": "symmetry_aware_rmsd",
        "architecture_name": str(architecture_name),
        "model_config": dict(DEFAULT_MODEL_CONFIG if model_config is None else model_config),
        "inverse_design_config": {
            "fixed_leading_atoms": int(fixed_leading_atoms),
            "fingerprint_loss_weight": float(fingerprint_loss_weight),
            "target_vertex_weight": float(target_vertex_weight),
            "target_position_weight": float(target_position_weight),
            "inverse_lr_decay_rate": float(inverse_lr_decay_rate),
            "inverse_restarts": int(inverse_restarts),
            "inverse_restart_noise_scale": float(inverse_restart_noise_scale),
        },
        "output_files": {
            "plot": str(figure_path),
            "original": str(output_dir / "torch_gnn_carbon_original.xyz"),
            "initial": str(output_dir / "torch_gnn_carbon_initial.xyz"),
            "final": str(output_dir / "torch_gnn_carbon_final.xyz"),
        },
    }
    metrics_path = output_dir / "torch_gnn_carbon_workflow_metrics.json"
    metrics["output_files"]["metrics"] = str(metrics_path)
    metrics_path.write_text(json.dumps(metrics, indent=2))
    return metrics


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--carbon-count", type=int, default=-1)
    parser.add_argument("--augmented-count", type=int, default=16)
    parser.add_argument("--epochs", type=int, default=30)
    parser.add_argument("--batch-size", type=int, default=8)
    parser.add_argument("--inverse-steps", type=int, default=400)
    parser.add_argument("--inverse-step-size", type=float, default=5.0e-3)
    parser.add_argument("--epoch-values", type=str, default="")
    parser.add_argument("--inverse-step-values", type=str, default="")
    parser.add_argument("--step-size-values", type=str, default="")
    parser.add_argument("--fixed-leading-atoms", type=int, default=FIXED_LEADING_ATOMS)
    parser.add_argument("--fingerprint-loss-weight", type=float, default=FINGERPRINT_LOSS_WEIGHT)
    parser.add_argument("--target-vertex-weight", type=float, default=TARGET_VERTEX_WEIGHT)
    parser.add_argument("--target-position-weight", type=float, default=TARGET_POSITION_WEIGHT)
    parser.add_argument("--inverse-lr-decay-rate", type=float, default=INVERSE_LR_DECAY_RATE)
    parser.add_argument("--inverse-restarts", type=int, default=INVERSE_RESTARTS)
    parser.add_argument(
        "--inverse-restart-noise-scale",
        type=float,
        default=INVERSE_RESTART_NOISE_SCALE,
    )
    parser.add_argument("--hidden-dim", type=int, default=DEFAULT_MODEL_CONFIG["hidden_dim"])
    parser.add_argument("--architecture", type=str, default=DEFAULT_MODEL_CONFIG["architecture"])
    parser.add_argument(
        "--num-message-layers",
        type=int,
        default=DEFAULT_MODEL_CONFIG["num_message_layers"],
    )
    parser.add_argument("--learning-rate", type=float, default=DEFAULT_MODEL_CONFIG["learning_rate"])
    parser.add_argument(
        "--model-lr-decay-rate",
        type=float,
        default=DEFAULT_MODEL_CONFIG["lr_decay_rate"],
    )
    parser.add_argument(
        "--smooth-cutoff-width",
        type=float,
        default=DEFAULT_MODEL_CONFIG["smooth_cutoff_width"],
    )
    parser.add_argument(
        "--reference-layer-type",
        type=int,
        default=DEFAULT_MODEL_CONFIG["reference_layer_type"],
    )
    parser.add_argument(
        "--component-weight-2body",
        type=float,
        default=DEFAULT_COMPONENT_WEIGHT[0],
    )
    parser.add_argument(
        "--component-weight-3body",
        type=float,
        default=DEFAULT_COMPONENT_WEIGHT[1],
    )
    parser.add_argument(
        "--component-weight-4body",
        type=float,
        default=DEFAULT_COMPONENT_WEIGHT[2],
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("build") / "torch_gnn_carbon_workflow_example",
    )
    parser.add_argument("--seed", type=int, default=42)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    epoch_values = parse_int_list(args.epoch_values) if args.epoch_values else default_epoch_values(args.epochs)
    inverse_step_values = (
        parse_int_list(args.inverse_step_values)
        if args.inverse_step_values
        else default_inverse_step_values(args.inverse_steps)
    )
    step_size_values = (
        parse_float_list(args.step_size_values)
        if args.step_size_values
        else default_step_sizes(args.inverse_step_size)
    )
    print(f"Using inverse step values: {inverse_step_values}")
    print(f"Using epoch values: {epoch_values}")
    print(f"Using step size values: {step_size_values}")
    print(f"Using inverse step size: {args.inverse_step_size}")
    print(f"Using fixed leading atoms: {args.fixed_leading_atoms}")
    print(f"Using fingerprint loss weight: {args.fingerprint_loss_weight}")
    print(f"Using target vertex weight: {args.target_vertex_weight}")
    print(f"Using target position weight: {args.target_position_weight}")
    print(f"Using inverse LR decay rate: {args.inverse_lr_decay_rate}")
    print(f"Using inverse restarts: {args.inverse_restarts}")
    print(f"Using inverse restart noise scale: {args.inverse_restart_noise_scale}")
    print(f"Using hidden dim: {args.hidden_dim}")
    print(f"Using architecture: {args.architecture}")
    print(f"Using message layers: {args.num_message_layers}")
    print(f"Using learning rate: {args.learning_rate}")
    print(f"Using model LR decay rate: {args.model_lr_decay_rate}")
    print(f"Using smooth cutoff width: {args.smooth_cutoff_width}")
    print(f"Using reference layer type: {args.reference_layer_type}")
    print(
        "Using component weights: "
        f"({args.component_weight_2body}, {args.component_weight_3body}, {args.component_weight_4body})"
    )

    repo_root = Path(__file__).resolve().parents[2]
    metrics = run_example(
        repo_root=repo_root,
        output_dir=args.output_dir,
        carbon_count=args.carbon_count,
        augmented_count=args.augmented_count,
        num_epochs=args.epochs,
        batch_size=args.batch_size,
        inverse_steps=args.inverse_steps,
        inverse_step_size=args.inverse_step_size,
        epoch_values=epoch_values,
        inverse_step_values=inverse_step_values,
        step_size_values=step_size_values,
        fixed_leading_atoms=args.fixed_leading_atoms,
        fingerprint_loss_weight=args.fingerprint_loss_weight,
        target_vertex_weight=args.target_vertex_weight,
        target_position_weight=args.target_position_weight,
        inverse_lr_decay_rate=args.inverse_lr_decay_rate,
        inverse_restarts=args.inverse_restarts,
        inverse_restart_noise_scale=args.inverse_restart_noise_scale,
        seed=args.seed,
        model_config={
            "architecture": str(args.architecture),
            "hidden_dim": int(args.hidden_dim),
            "num_message_layers": int(args.num_message_layers),
            "learning_rate": float(args.learning_rate),
            "lr_decay_rate": float(args.model_lr_decay_rate),
            "smooth_cutoff_width": float(args.smooth_cutoff_width),
            "reference_layer_type": int(args.reference_layer_type),
            "component_weight": (
                float(args.component_weight_2body),
                float(args.component_weight_3body),
                float(args.component_weight_4body),
            ),
        },
        architecture_name=str(args.architecture),
    )
    print()
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
