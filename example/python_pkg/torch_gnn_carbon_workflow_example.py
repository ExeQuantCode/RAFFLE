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
from torch_gnn_rollout import (
    PrioritizedReplayBuffer,
    ReplaySample,
    classify_rollout_step,
)
from torch_gnn_workflow_common import save_descriptor_comparison_report


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

FIXED_LEADING_ATOMS = 1
FINGERPRINT_LOSS_WEIGHT = 0.5
TARGET_VERTEX_WEIGHT = 0.0
TARGET_POSITION_WEIGHT = 0.0
INVERSE_LR_DECAY_RATE = 0.0
INVERSE_RESTARTS = 1
INVERSE_RESTART_NOISE_SCALE = 0.0
REPULSION_WEIGHT = 10.0
MINIMUM_DISTANCE_SCALE = 0.75
CELL_VIOLATION_WEIGHT = 0.0
COORDINATE_CLIP_VALUE = None
ROLLOUT_STAGES = 1
ROLLOUT_EPOCHS_PER_STAGE = 1
ROLLOUT_STEP_STRIDE = 25
REPLAY_BUFFER_CAPACITY = 64
REPLAY_SAMPLE_SIZE = 8
ROLLOUT_DRIFT_THRESHOLD = 5.0e-4
ROLLOUT_HIGH_ERROR_THRESHOLD = 1.0e-3
ROLLOUT_INSTABILITY_THRESHOLD = 1.0e-4
DEFAULT_REPLAY_CATEGORY_WEIGHTS = {
    "successful": 1.0,
    "failed": 2.0,
    "unstable": 2.0,
    "high_error": 2.0,
    "difficult": 2.0,
}
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


def parse_category_weights(value: str) -> dict[str, float]:
    resolved = dict(DEFAULT_REPLAY_CATEGORY_WEIGHTS)
    for item in value.split(","):
        entry = item.strip()
        if not entry:
            continue
        key, separator, raw_value = entry.partition("=")
        if separator != "=" or not key.strip() or not raw_value.strip():
            raise ValueError(
                "replay-category-weights entries must use category=value format"
            )
        resolved[key.strip()] = float(raw_value.strip())
    return resolved


def build_inverse_design_options(
    original,
    fingerprint_loss_weight: float,
    target_vertex_weight: float,
    target_position_weight: float,
    inverse_lr_decay_rate: float,
    inverse_restarts: int,
    inverse_restart_noise_scale: float,
    repulsion_weight: float,
    minimum_distance_scale: float,
    cell_violation_weight: float,
    coordinate_clip_value: float | None,
) -> dict:
    if float(target_vertex_weight) != 0.0:
        raise ValueError("target_vertex_weight must remain 0.0 for inverse-design runs")
    if float(target_position_weight) != 0.0:
        raise ValueError("target_position_weight must remain 0.0 for inverse-design runs")
    return {
        "target_atoms": None,
        "fingerprint_loss_weight": float(fingerprint_loss_weight),
        "target_vertex_weight": float(target_vertex_weight),
        "target_position_weight": float(target_position_weight),
        "inverse_lr_decay_rate": float(inverse_lr_decay_rate),
        "num_restarts": int(inverse_restarts),
        "restart_noise_scale": float(inverse_restart_noise_scale),
        "repulsion_weight": float(repulsion_weight),
        "minimum_distance_scale": float(minimum_distance_scale),
        "cell_violation_weight": float(cell_violation_weight),
        "coordinate_clip_value": (
            None if coordinate_clip_value is None else float(coordinate_clip_value)
        ),
    }


def minimum_training_carbon_count(total_count: int) -> int:
    return max(1, int(np.ceil(0.3 * max(int(total_count), 0))))


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


def build_rollout_step_observer(trajectory_records: list[dict]) -> Callable[[dict[str, object]], None]:
    def observer(step_record: dict[str, object]) -> None:
        trajectory_records.append(
            {
                "atoms": step_record["atoms"].copy(),
                "restart_index": int(step_record["restart_index"]),
                "num_restarts": int(step_record["num_restarts"]),
                "step": int(step_record["step"]),
                "num_steps": int(step_record["num_steps"]),
                "is_initial_state": bool(step_record["is_initial_state"]),
                "learning_rate": float(step_record.get("learning_rate", 0.0)),
                "total_loss": float(step_record.get("total_loss", 0.0)),
                "fingerprint_loss": float(step_record.get("fingerprint_loss", 0.0)),
                "repulsion_loss": float(step_record.get("repulsion_loss", 0.0)),
                "cell_violation_loss": float(step_record.get("cell_violation_loss", 0.0)),
            }
        )

    return observer


def save_inverse_design_path(
    output_dir: Path,
    trajectory_records: list[dict],
    *,
    prefix: str,
) -> dict[str, object]:
    if not trajectory_records:
        return {
            "traj_file": None,
            "step_structure_dir": None,
            "step_structure_files": [],
            "steps": [],
        }

    step_structure_dir = output_dir / f"{prefix}_steps"
    step_structure_dir.mkdir(parents=True, exist_ok=True)

    trajectory_frames = []
    step_structure_files: list[str] = []
    step_summaries: list[dict[str, object]] = []
    for record in trajectory_records:
        atoms_snapshot = record["atoms"].copy()
        restart_index = int(record["restart_index"]) + 1
        step = int(record["step"])
        num_steps = int(record["num_steps"])
        is_initial_state = bool(record["is_initial_state"])

        atoms_snapshot.info["inverse_restart_index"] = restart_index
        atoms_snapshot.info["inverse_num_restarts"] = int(record["num_restarts"])
        atoms_snapshot.info["inverse_step"] = step
        atoms_snapshot.info["inverse_num_steps"] = num_steps
        atoms_snapshot.info["inverse_is_initial_state"] = is_initial_state
        atoms_snapshot.info["inverse_learning_rate"] = float(record.get("learning_rate", 0.0))
        atoms_snapshot.info["inverse_total_loss"] = float(record.get("total_loss", 0.0))
        atoms_snapshot.info["inverse_fingerprint_loss"] = float(
            record.get("fingerprint_loss", 0.0)
        )
        atoms_snapshot.info["inverse_repulsion_loss"] = float(record.get("repulsion_loss", 0.0))
        atoms_snapshot.info["inverse_cell_violation_loss"] = float(
            record.get("cell_violation_loss", 0.0)
        )

        step_label = f"{prefix}_restart_{restart_index:02d}_step_{step:04d}"
        if is_initial_state:
            step_label += "_initial"
        elif step == num_steps:
            step_label += "_final"
        step_path = step_structure_dir / f"{step_label}.xyz"
        write(step_path, atoms_snapshot)

        trajectory_frames.append(atoms_snapshot)
        step_structure_files.append(str(step_path))
        step_summaries.append(
            {
                "restart_index": restart_index,
                "num_restarts": int(record["num_restarts"]),
                "step": step,
                "num_steps": num_steps,
                "is_initial_state": is_initial_state,
                "learning_rate": float(record.get("learning_rate", 0.0)),
                "total_loss": float(record.get("total_loss", 0.0)),
                "fingerprint_loss": float(record.get("fingerprint_loss", 0.0)),
                "repulsion_loss": float(record.get("repulsion_loss", 0.0)),
                "cell_violation_loss": float(record.get("cell_violation_loss", 0.0)),
                "structure_file": str(step_path),
            }
        )

    traj_path = output_dir / f"{prefix}_path.traj"
    write(traj_path, trajectory_frames)
    return {
        "traj_file": str(traj_path),
        "step_structure_dir": str(step_structure_dir),
        "step_structure_files": step_structure_files,
        "steps": step_summaries,
    }


def summarise_convergence(trace: list[dict[str, float]]) -> dict[str, float | int | bool]:
    if not trace:
        return {
            "best_step": 0,
            "best_position_difference": 0.0,
            "final_position_difference": 0.0,
            "best_fingerprint_mse": 0.0,
            "final_fingerprint_mse": 0.0,
            "tail_position_range": 0.0,
            "tail_fingerprint_range": 0.0,
            "fingerprint_nonincreasing_fraction": 1.0,
            "position_nonincreasing_fraction": 1.0,
            "converges_to_best_within_5pct": True,
        }

    position_values = np.asarray([entry["position_difference"] for entry in trace], dtype=np.float64)
    fingerprint_values = np.asarray([entry["fingerprint_mse"] for entry in trace], dtype=np.float64)
    step_values = np.asarray([entry["step"] for entry in trace], dtype=np.int64)
    best_index = int(np.argmin(position_values))
    tail_count = min(3, len(trace))
    tail_positions = position_values[-tail_count:]
    tail_fingerprint = fingerprint_values[-tail_count:]
    fingerprint_deltas = np.diff(fingerprint_values)
    position_deltas = np.diff(position_values)

    return {
        "best_step": int(step_values[best_index]),
        "best_position_difference": float(position_values[best_index]),
        "final_position_difference": float(position_values[-1]),
        "best_fingerprint_mse": float(np.min(fingerprint_values)),
        "final_fingerprint_mse": float(fingerprint_values[-1]),
        "tail_position_range": float(np.max(tail_positions) - np.min(tail_positions)),
        "tail_fingerprint_range": float(np.max(tail_fingerprint) - np.min(tail_fingerprint)),
        "fingerprint_nonincreasing_fraction": (
            1.0
            if len(fingerprint_deltas) == 0
            else float(np.mean(fingerprint_deltas <= 1.0e-12))
        ),
        "position_nonincreasing_fraction": (
            1.0 if len(position_deltas) == 0 else float(np.mean(position_deltas <= 1.0e-12))
        ),
        "converges_to_best_within_5pct": bool(
            position_values[-1] <= 1.05 * max(position_values[best_index], 1.0e-12)
        ),
    }


def run_rollout_retraining(
    model: TorchGNNFingerprint,
    carbon_structures,
    original,
    perturbed,
    fixed_atoms: np.ndarray,
    target_fingerprint: np.ndarray,
    batch_size: int,
    inverse_steps: int,
    inverse_step_size: float,
    inverse_design_options: dict,
    seed: int,
    rollout_stages: int,
    rollout_epochs_per_stage: int,
    rollout_step_stride: int,
    replay_buffer_capacity: int,
    replay_sample_size: int,
    replay_category_weights: dict[str, float],
    rollout_drift_threshold: float,
    rollout_high_error_threshold: float,
    rollout_instability_threshold: float,
    training_observer: Optional[Callable[[int, float], None]] = None,
    training_epoch_offset: int = 0,
) -> tuple[TorchGNNFingerprint, dict]:
    if int(rollout_stages) <= 0:
        return model, {
            "enabled": False,
            "num_rollout_epochs": 0,
            "replay_buffer": {
                "capacity": int(replay_buffer_capacity),
                "size": 0,
                "category_counts": {},
                "mean_priority": 0.0,
                "mean_fingerprint_drift_mse": 0.0,
                "mean_true_target_fingerprint_mse": 0.0,
            },
            "stages": [],
        }

    replay_buffer = PrioritizedReplayBuffer(capacity=replay_buffer_capacity, seed=seed)
    stage_summaries = []
    rollout_epoch_index = int(training_epoch_offset)

    for stage_index in range(int(rollout_stages)):
        trajectory_records: list[dict] = []
        step_observer = build_rollout_step_observer(trajectory_records)
        model.inverse_design(
            target_fingerprint=target_fingerprint,
            atoms=perturbed,
            fixed_atoms=fixed_atoms,
            num_steps=inverse_steps,
            step_size=inverse_step_size,
            verbose=0,
            step_observer=step_observer,
            **inverse_design_options,
        )

        retained_records = []
        last_step = max((record["step"] for record in trajectory_records), default=0)
        for record in trajectory_records:
            if record["is_initial_state"]:
                retained_records.append(record)
                continue
            if int(record["step"]) == int(last_step):
                retained_records.append(record)
                continue
            if int(rollout_step_stride) <= 1 or int(record["step"]) % int(rollout_step_stride) == 0:
                retained_records.append(record)

        new_samples = []
        previous_true_error = None
        category_counts: dict[str, int] = {}
        for record in retained_records:
            true_fingerprint = model.compute_reference_fingerprint(record["atoms"])
            surrogate_fingerprint = model.predict(record["atoms"])
            true_target_fingerprint_mse = float(
                np.mean((true_fingerprint - target_fingerprint) ** 2)
            )
            surrogate_target_fingerprint_mse = float(
                np.mean((surrogate_fingerprint - target_fingerprint) ** 2)
            )
            fingerprint_drift_mse = float(
                np.mean((surrogate_fingerprint - true_fingerprint) ** 2)
            )
            position_difference = float(score_candidate(original, record["atoms"]))
            categories = classify_rollout_step(
                true_target_fingerprint_mse=true_target_fingerprint_mse,
                previous_true_target_fingerprint_mse=previous_true_error,
                fingerprint_drift_mse=fingerprint_drift_mse,
                repulsion_loss=float(record["repulsion_loss"]),
                cell_violation_loss=float(record["cell_violation_loss"]),
                drift_threshold=rollout_drift_threshold,
                high_error_threshold=rollout_high_error_threshold,
                instability_threshold=rollout_instability_threshold,
            )
            for category in categories:
                category_counts[category] = category_counts.get(category, 0) + 1
            priority = (
                true_target_fingerprint_mse
                + surrogate_target_fingerprint_mse
                + fingerprint_drift_mse
                + float(record["repulsion_loss"])
                + float(record["cell_violation_loss"])
            )
            new_sample = ReplaySample(
                atoms=record["atoms"],
                priority=priority,
                categories=categories,
                stage_index=stage_index + 1,
                step=int(record["step"]),
                restart_index=int(record["restart_index"]),
                true_target_fingerprint_mse=true_target_fingerprint_mse,
                surrogate_target_fingerprint_mse=surrogate_target_fingerprint_mse,
                fingerprint_drift_mse=fingerprint_drift_mse,
                position_difference=position_difference,
                repulsion_loss=float(record["repulsion_loss"]),
                cell_violation_loss=float(record["cell_violation_loss"]),
                is_initial_state=bool(record["is_initial_state"]),
            )
            new_samples.append(new_sample)
            previous_true_error = true_target_fingerprint_mse

        replay_buffer.extend(new_samples)
        replay_batch = replay_buffer.sample(
            sample_size=replay_sample_size,
            category_weights=replay_category_weights,
        )

        stage_training_losses = []
        if replay_batch and int(rollout_epochs_per_stage) > 0:
            replay_structures = [sample.atoms for sample in replay_batch]
            for _ in range(int(rollout_epochs_per_stage)):
                history = model.fit(
                    carbon_structures,
                    num_epochs=1,
                    batch_size=batch_size,
                    augment_structures=replay_structures,
                    verbose=0,
                    reset_optimiser=False,
                    recalibrate_base=False,
                )
                epoch_loss = float(history[-1])
                stage_training_losses.append(epoch_loss)
                rollout_epoch_index += 1
                if training_observer is not None:
                    training_observer(rollout_epoch_index, epoch_loss)

        stage_true_errors = [sample.true_target_fingerprint_mse for sample in new_samples]
        stage_surrogate_errors = [sample.surrogate_target_fingerprint_mse for sample in new_samples]
        stage_drifts = [sample.fingerprint_drift_mse for sample in new_samples]
        stage_position_differences = [
            sample.position_difference for sample in new_samples if sample.position_difference is not None
        ]
        stage_summaries.append(
            {
                "stage_index": int(stage_index + 1),
                "num_stage_samples": int(len(new_samples)),
                "sampled_replay_count": int(len(replay_batch)),
                "buffer_size": int(len(replay_buffer)),
                "category_counts": category_counts,
                "initial_true_target_fingerprint_mse": (
                    None if not stage_true_errors else float(stage_true_errors[0])
                ),
                "best_true_target_fingerprint_mse": (
                    None if not stage_true_errors else float(np.min(stage_true_errors))
                ),
                "final_true_target_fingerprint_mse": (
                    None if not stage_true_errors else float(stage_true_errors[-1])
                ),
                "best_surrogate_target_fingerprint_mse": (
                    None if not stage_surrogate_errors else float(np.min(stage_surrogate_errors))
                ),
                "final_surrogate_target_fingerprint_mse": (
                    None if not stage_surrogate_errors else float(stage_surrogate_errors[-1])
                ),
                "mean_fingerprint_drift_mse": (
                    0.0 if not stage_drifts else float(np.mean(stage_drifts))
                ),
                "best_position_difference": (
                    None
                    if not stage_position_differences
                    else float(np.min(stage_position_differences))
                ),
                "final_position_difference": (
                    None
                    if not stage_position_differences
                    else float(stage_position_differences[-1])
                ),
                "training_losses": [float(loss) for loss in stage_training_losses],
                "trajectory": [
                    {
                        "restart_index": int(sample.restart_index),
                        "step": int(sample.step),
                        "is_initial_state": bool(sample.is_initial_state),
                        "categories": list(sample.categories),
                        "priority": float(sample.priority),
                        "true_target_fingerprint_mse": float(sample.true_target_fingerprint_mse),
                        "surrogate_target_fingerprint_mse": float(sample.surrogate_target_fingerprint_mse),
                        "fingerprint_drift_mse": float(sample.fingerprint_drift_mse),
                        "position_difference": (
                            None
                            if sample.position_difference is None
                            else float(sample.position_difference)
                        ),
                        "repulsion_loss": float(sample.repulsion_loss),
                        "cell_violation_loss": float(sample.cell_violation_loss),
                    }
                    for sample in new_samples
                ],
            }
        )

    return model, {
        "enabled": True,
        "num_rollout_epochs": int(max(rollout_epoch_index - int(training_epoch_offset), 0)),
        "replay_buffer": replay_buffer.stats(),
        "stages": stage_summaries,
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
        predicted_fingerprint = model.predict(candidate_atoms)
        fingerprint_mse = float(np.mean((predicted_fingerprint - target_fingerprint) ** 2))
        fingerprint_l2 = float(np.linalg.norm(predicted_fingerprint - target_fingerprint))
        structure_update_norm = float(
            np.linalg.norm(candidate_atoms.get_positions() - perturbed.get_positions())
        )
        if step == requested_step_values[-1]:
            optimised = candidate_atoms.copy()
        best_candidate = update_best_candidate(
            best_candidate,
            candidate_atoms,
            position_difference,
            source="inverse_step_sweep",
            step=int(step),
            fingerprint_mse=fingerprint_mse,
        )
        trace.append(
            {
                "step": step,
                "position_difference": position_difference,
                "fingerprint_mse": fingerprint_mse,
                "fingerprint_l2": fingerprint_l2,
                "structure_update_norm": structure_update_norm,
            }
        )
    return optimised, trace, best_candidate


def run_configured_inverse_design(
    model: TorchGNNFingerprint,
    perturbed,
    fixed_atoms: np.ndarray,
    target_fingerprint: np.ndarray,
    inverse_steps: int,
    inverse_step_size: float,
    inverse_design_options: dict,
) -> tuple[object, list[dict]]:
    trajectory_records: list[dict] = []
    step_observer = build_rollout_step_observer(trajectory_records)
    optimised = model.inverse_design(
        target_fingerprint=target_fingerprint,
        atoms=perturbed,
        fixed_atoms=fixed_atoms,
        num_steps=inverse_steps,
        step_size=inverse_step_size,
        verbose=0,
        step_observer=step_observer,
        **inverse_design_options,
    )
    return optimised, trajectory_records


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
    repulsion_weight: float,
    minimum_distance_scale: float,
    cell_violation_weight: float,
    coordinate_clip_value: float | None,
    rollout_stages: int,
    rollout_epochs_per_stage: int,
    rollout_step_stride: int,
    replay_buffer_capacity: int,
    replay_sample_size: int,
    replay_category_weights: dict[str, float],
    rollout_drift_threshold: float,
    rollout_high_error_threshold: float,
    rollout_instability_threshold: float,
    seed: int,
    model_config: Optional[dict] = None,
    training_observer: Optional[Callable[[int, float], None]] = None,
    architecture_name: str = "torch_gnn_residual",
    enable_checkpoint_step_size_sweep: bool = True,
    enable_checkpoint_step_schedule_sweep: bool = True,
) -> dict:
    if int(augmented_count) != 0:
        raise ValueError("augmented_count must remain 0 for inverse-design runs")
    if float(target_vertex_weight) != 0.0:
        raise ValueError("target_vertex_weight must remain 0.0 for inverse-design runs")
    if float(target_position_weight) != 0.0:
        raise ValueError("target_position_weight must remain 0.0 for inverse-design runs")

    carbon_xyz = repo_root / "example" / "data" / "carbon.xyz"
    all_carbon_structures = read(str(carbon_xyz), index=":")
    carbon_structures = select_carbon_structures(all_carbon_structures, carbon_count)
    minimum_carbon_count = minimum_training_carbon_count(len(all_carbon_structures))
    if len(carbon_structures) < minimum_carbon_count:
        raise ValueError(
            "carbon_count must select at least 30% of example/data/carbon.xyz "
            f"structures ({minimum_carbon_count} minimum, received {len(carbon_structures)})"
        )

    original = bulk("C", "diamond", a=3.567, cubic=True)
    original.pbc = True
    augmented_structures = []
    perturbed, fixed_atoms = build_perturbed_structure(original, fixed_leading_atoms=fixed_leading_atoms)
    inverse_design_options = build_inverse_design_options(
        original=original,
        fingerprint_loss_weight=fingerprint_loss_weight,
        target_vertex_weight=target_vertex_weight,
        target_position_weight=target_position_weight,
        inverse_lr_decay_rate=inverse_lr_decay_rate,
        inverse_restarts=inverse_restarts,
        inverse_restart_noise_scale=inverse_restart_noise_scale,
        repulsion_weight=repulsion_weight,
        minimum_distance_scale=minimum_distance_scale,
        cell_violation_weight=cell_violation_weight,
        coordinate_clip_value=coordinate_clip_value,
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
    model, rollout_metrics = run_rollout_retraining(
        model=model,
        carbon_structures=carbon_structures,
        original=original,
        perturbed=perturbed,
        fixed_atoms=fixed_atoms,
        target_fingerprint=target_fingerprint,
        batch_size=batch_size,
        inverse_steps=inverse_steps,
        inverse_step_size=inverse_step_size,
        inverse_design_options=inverse_design_options,
        seed=seed,
        rollout_stages=rollout_stages,
        rollout_epochs_per_stage=rollout_epochs_per_stage,
        rollout_step_stride=rollout_step_stride,
        replay_buffer_capacity=replay_buffer_capacity,
        replay_sample_size=replay_sample_size,
        replay_category_weights=replay_category_weights,
        rollout_drift_threshold=rollout_drift_threshold,
        rollout_high_error_threshold=rollout_high_error_threshold,
        rollout_instability_threshold=rollout_instability_threshold,
        training_observer=training_observer,
        training_epoch_offset=len(training_losses),
    )
    if rollout_metrics["enabled"]:
        final_rollout_position_difference = rollout_metrics["stages"][-1]["final_position_difference"]
        if final_rollout_position_difference is not None:
            checkpoint_records.append(
                capture_checkpoint(
                    model,
                    int(num_epochs + rollout_metrics["num_rollout_epochs"]),
                    float(final_rollout_position_difference),
                )
            )

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
    convergence_summary = summarise_convergence(inverse_trace)
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
    checkpoint_step_size_results = []
    best_checkpoint_step_size_candidate = None
    if enable_checkpoint_step_size_sweep:
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
    checkpoint_step_schedule_results = []
    best_checkpoint_step_schedule_candidate = None
    if enable_checkpoint_step_schedule_sweep:
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

    configured_optimised, configured_trajectory_records = run_configured_inverse_design(
        model=model,
        perturbed=perturbed,
        fixed_atoms=fixed_atoms,
        target_fingerprint=target_fingerprint,
        inverse_steps=inverse_steps,
        inverse_step_size=inverse_step_size,
        inverse_design_options=inverse_design_options,
    )
    configured_final_rmsd = score_candidate(original, configured_optimised)
    best_candidate = update_best_candidate(
        None,
        configured_optimised,
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
    configured_inverse_design_path = save_inverse_design_path(
        output_dir,
        configured_trajectory_records,
        prefix="torch_gnn_carbon_inverse_design",
    )
    descriptor_report = save_descriptor_comparison_report(
        model=model,
        target_fingerprint=target_fingerprint,
        final_atoms=optimised,
        structure_path=output_dir / "torch_gnn_carbon_final.xyz",
    )

    initial_prediction = model.predict(perturbed)
    final_prediction = model.predict(optimised)
    initial_fingerprint_mse = float(np.mean((initial_prediction - target_fingerprint) ** 2))
    final_fingerprint_mse = float(np.mean((final_prediction - target_fingerprint) ** 2))
    initial_fingerprint_l2 = float(np.linalg.norm(initial_prediction - target_fingerprint))
    final_fingerprint_l2 = float(np.linalg.norm(final_prediction - target_fingerprint))

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
        "convergence_summary": convergence_summary,
        "rollout": rollout_metrics,
        "checkpoint_step_size_sweep": checkpoint_step_size_results,
        "checkpoint_step_schedule_sweep": checkpoint_step_schedule_results,
        "initial_fingerprint_mse": initial_fingerprint_mse,
        "final_fingerprint_mse": final_fingerprint_mse,
        "initial_fingerprint_l2": initial_fingerprint_l2,
        "final_fingerprint_l2": final_fingerprint_l2,
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
            "repulsion_weight": float(repulsion_weight),
            "minimum_distance_scale": float(minimum_distance_scale),
            "cell_violation_weight": float(cell_violation_weight),
            "coordinate_clip_value": (
                None if coordinate_clip_value is None else float(coordinate_clip_value)
            ),
            "rollout_stages": int(rollout_stages),
            "rollout_epochs_per_stage": int(rollout_epochs_per_stage),
            "rollout_step_stride": int(rollout_step_stride),
            "replay_buffer_capacity": int(replay_buffer_capacity),
            "replay_sample_size": int(replay_sample_size),
            "replay_category_weights": {
                key: float(value) for key, value in replay_category_weights.items()
            },
            "rollout_drift_threshold": float(rollout_drift_threshold),
            "rollout_high_error_threshold": float(rollout_high_error_threshold),
            "rollout_instability_threshold": float(rollout_instability_threshold),
        },
        "configured_inverse_design_path": configured_inverse_design_path,
        "output_files": {
            "plot": str(figure_path),
            "original": str(output_dir / "torch_gnn_carbon_original.xyz"),
            "initial": str(output_dir / "torch_gnn_carbon_initial.xyz"),
            "final": str(output_dir / "torch_gnn_carbon_final.xyz"),
            "configured_inverse_design_traj": str(
                configured_inverse_design_path["traj_file"]
            ),
            "final_descriptor_comparison": descriptor_report["report_file"],
            "final_descriptor_comparison_plot": descriptor_report["plot_file"],
        },
        "descriptor_comparisons": {"final": descriptor_report},
    }
    metrics_path = output_dir / "torch_gnn_carbon_workflow_metrics.json"
    metrics["output_files"]["metrics"] = str(metrics_path)
    metrics_path.write_text(json.dumps(metrics, indent=2))
    return metrics


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--carbon-count", type=int, default=-1)
    parser.add_argument("--augmented-count", type=int, default=0)
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
    parser.add_argument("--repulsion-weight", type=float, default=REPULSION_WEIGHT)
    parser.add_argument("--minimum-distance-scale", type=float, default=MINIMUM_DISTANCE_SCALE)
    parser.add_argument("--cell-violation-weight", type=float, default=CELL_VIOLATION_WEIGHT)
    parser.add_argument("--coordinate-clip-value", type=float, default=COORDINATE_CLIP_VALUE)
    parser.add_argument("--rollout-stages", type=int, default=ROLLOUT_STAGES)
    parser.add_argument(
        "--rollout-epochs-per-stage",
        type=int,
        default=ROLLOUT_EPOCHS_PER_STAGE,
    )
    parser.add_argument("--rollout-step-stride", type=int, default=ROLLOUT_STEP_STRIDE)
    parser.add_argument("--replay-buffer-capacity", type=int, default=REPLAY_BUFFER_CAPACITY)
    parser.add_argument("--replay-sample-size", type=int, default=REPLAY_SAMPLE_SIZE)
    parser.add_argument(
        "--replay-category-weights",
        type=str,
        default=",".join(
            f"{key}={value}" for key, value in DEFAULT_REPLAY_CATEGORY_WEIGHTS.items()
        ),
    )
    parser.add_argument("--rollout-drift-threshold", type=float, default=ROLLOUT_DRIFT_THRESHOLD)
    parser.add_argument(
        "--rollout-high-error-threshold",
        type=float,
        default=ROLLOUT_HIGH_ERROR_THRESHOLD,
    )
    parser.add_argument(
        "--rollout-instability-threshold",
        type=float,
        default=ROLLOUT_INSTABILITY_THRESHOLD,
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
    parser.add_argument("--skip-checkpoint-step-size-sweep", action="store_true")
    parser.add_argument("--skip-checkpoint-step-schedule-sweep", action="store_true")
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
    print(f"Using repulsion weight: {args.repulsion_weight}")
    print(f"Using minimum distance scale: {args.minimum_distance_scale}")
    print(f"Using cell violation weight: {args.cell_violation_weight}")
    print(f"Using coordinate clip value: {args.coordinate_clip_value}")
    print(f"Using rollout stages: {args.rollout_stages}")
    print(f"Using rollout epochs per stage: {args.rollout_epochs_per_stage}")
    print(f"Using rollout step stride: {args.rollout_step_stride}")
    print(f"Using replay buffer capacity: {args.replay_buffer_capacity}")
    print(f"Using replay sample size: {args.replay_sample_size}")
    print(f"Using replay category weights: {args.replay_category_weights}")
    print(f"Using rollout drift threshold: {args.rollout_drift_threshold}")
    print(f"Using rollout high-error threshold: {args.rollout_high_error_threshold}")
    print(f"Using rollout instability threshold: {args.rollout_instability_threshold}")
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
        repulsion_weight=args.repulsion_weight,
        minimum_distance_scale=args.minimum_distance_scale,
        cell_violation_weight=args.cell_violation_weight,
        coordinate_clip_value=args.coordinate_clip_value,
        rollout_stages=args.rollout_stages,
        rollout_epochs_per_stage=args.rollout_epochs_per_stage,
        rollout_step_stride=args.rollout_step_stride,
        replay_buffer_capacity=args.replay_buffer_capacity,
        replay_sample_size=args.replay_sample_size,
        replay_category_weights=parse_category_weights(args.replay_category_weights),
        rollout_drift_threshold=args.rollout_drift_threshold,
        rollout_high_error_threshold=args.rollout_high_error_threshold,
        rollout_instability_threshold=args.rollout_instability_threshold,
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
        enable_checkpoint_step_size_sweep=not args.skip_checkpoint_step_size_sweep,
        enable_checkpoint_step_schedule_sweep=not args.skip_checkpoint_step_schedule_sweep,
    )
    print()
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
