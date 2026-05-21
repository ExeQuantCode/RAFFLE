"""Legacy combined carbon workflow example for TorchGNN inverse design.

The supported reproducible path is the split workflow in torch_gnn_train_model.py
and torch_gnn_inverse_design.py. This script remains useful for local carbon-only
experimentation and sweep prototyping.
"""

from __future__ import annotations

import argparse
import copy
import importlib.util
import json
from pathlib import Path
import sys
from typing import Callable, Optional

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from ase.build import bulk
from ase.io import read, write

def _load_local_torch_gnn_fingerprint() -> None:
    import raffle as raffle_package

    module_path = Path(__file__).resolve().parents[2] / "src" / "raffle" / "torch_gnn_fingerprint.py"
    spec = importlib.util.spec_from_file_location("raffle.torch_gnn_fingerprint", module_path)
    if spec is None or spec.loader is None:
        return
    module = importlib.util.module_from_spec(spec)
    sys.modules["raffle.torch_gnn_fingerprint"] = module
    spec.loader.exec_module(module)
    raffle_package.TorchGNNFingerprint = module.TorchGNNFingerprint


_load_local_torch_gnn_fingerprint()

from raffle import (
    TorchGNNFingerprint,
    minimum_image_displacements,
    structure_similarity_rmsd,
    symmetry_aware_displacements,
)

try:
    from raffle import InferenceEnsemble, InferenceEnsembleConfig
except ImportError:
    InferenceEnsemble = None
    InferenceEnsembleConfig = None
from torch_gnn_rollout import (
    PrioritizedReplayBuffer,
    ReplaySample,
    classify_rollout_step,
)
from torch_gnn_workflow_common import (
    save_descriptor_comparison_report,
    save_model_checkpoint,
    save_target_fingerprint,
    write_structure,
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

FIXED_LEADING_ATOMS = 1
FINGERPRINT_LOSS_WEIGHT = 0.5
TARGET_VERTEX_WEIGHT = 0.0
TARGET_POSITION_WEIGHT = 0.0
INVERSE_LR_DECAY_RATE = 0.0
# Deprecated compatibility defaults. Restart controls are ignored by the workflow.
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
INVERSE_STEP_LOG_INTERVAL = 10


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
    fingerprint_loss_weight: float,
    target_vertex_weight: float,
    target_position_weight: float,
    inverse_lr_decay_rate: float,
    repulsion_weight: float,
    minimum_distance_scale: float,
    cell_violation_weight: float,
    coordinate_clip_value: float | None,
    inverse_restarts: int | None = None,
    inverse_restart_noise_scale: float | None = None,
) -> dict:
    del inverse_restarts
    del inverse_restart_noise_scale
    if float(target_vertex_weight) != 0.0:
        raise ValueError("target_vertex_weight must remain 0.0 for inverse-design runs")
    if float(target_position_weight) != 0.0:
        raise ValueError("target_position_weight must remain 0.0 for inverse-design runs")
    return {
        "fingerprint_loss_weight": float(fingerprint_loss_weight),
        "target_vertex_weight": float(target_vertex_weight),
        "target_position_weight": float(target_position_weight),
        "inverse_lr_decay_rate": float(inverse_lr_decay_rate),
        "repulsion_weight": float(repulsion_weight),
        "minimum_distance_scale": float(minimum_distance_scale),
        "cell_violation_weight": float(cell_violation_weight),
        "coordinate_clip_value": (
            None if coordinate_clip_value is None else float(coordinate_clip_value)
        ),
    }


def build_ensemble_config(
    ensemble_enabled: bool = False,
    ensemble_num_trajectories: int = 16,
    ensemble_perturbation_scale: float = 0.01,
    ensemble_langevin_noise_scale: float = 0.005,
    ensemble_temperature: float = 1.0,
    ensemble_adaptive_scaling: bool = True,
    ensemble_trajectory_pruning: bool = True,
    ensemble_escape_detection: bool = True,
    ensemble_consensus_metric: str = "cluster",
    ensemble_aggregation: str = "mean_variance",
) -> dict:
    """Build ensemble configuration dictionary.

    Returns a dict that can be passed to InferenceEnsembleConfig or None if disabled.
    """
    if not ensemble_enabled:
        return None

    return {
        "enabled": True,
        "num_trajectories": int(ensemble_num_trajectories),
        "perturbation_scale": float(ensemble_perturbation_scale),
        "langevin_noise_scale": float(ensemble_langevin_noise_scale),
        "temperature": float(ensemble_temperature),
        "adaptive_scaling": bool(ensemble_adaptive_scaling),
        "trajectory_pruning": bool(ensemble_trajectory_pruning),
        "escape_detection": bool(ensemble_escape_detection),
        "consensus_metric": str(ensemble_consensus_metric),
        "aggregation": str(ensemble_aggregation),
    }


def minimum_training_carbon_count(total_count: int) -> int:
    return max(1, int(np.ceil(0.3 * max(int(total_count), 0))))


def score_candidate(original, candidate_atoms) -> float:
    return float(structure_similarity_rmsd(original, candidate_atoms))


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
        # Ensemble callbacks may emit summary-only records without atomic snapshots.
        # Skip those here and only retain trajectory-style records.
        if "atoms" not in step_record:
            return
        trajectory_records.append(
            {
                "atoms": step_record["atoms"].copy(),
                "step": int(step_record.get("step", 0)),
                "num_steps": int(step_record.get("num_steps", 0)),
                "is_initial_state": bool(step_record.get("is_initial_state", False)),
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
        step = int(record["step"])
        num_steps = int(record["num_steps"])
        is_initial_state = bool(record["is_initial_state"])

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
        if "position_difference" in record:
            atoms_snapshot.info["inverse_position_difference"] = float(
                record["position_difference"]
            )
        if "fingerprint_mse" in record:
            atoms_snapshot.info["inverse_fingerprint_mse"] = float(record["fingerprint_mse"])
        if "fingerprint_l2" in record:
            atoms_snapshot.info["inverse_fingerprint_l2"] = float(record["fingerprint_l2"])
        if "structure_update_norm" in record:
            atoms_snapshot.info["inverse_structure_update_norm"] = float(
                record["structure_update_norm"]
            )

        step_label = f"{prefix}_step_{step:04d}"
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
                "step": step,
                "num_steps": num_steps,
                "is_initial_state": is_initial_state,
                "learning_rate": float(record.get("learning_rate", 0.0)),
                "total_loss": float(record.get("total_loss", 0.0)),
                "fingerprint_loss": float(record.get("fingerprint_loss", 0.0)),
                "repulsion_loss": float(record.get("repulsion_loss", 0.0)),
                "cell_violation_loss": float(record.get("cell_violation_loss", 0.0)),
                "structure_file": str(step_path),
                "position_difference": record.get("position_difference"),
                "fingerprint_mse": record.get("fingerprint_mse"),
                "fingerprint_l2": record.get("fingerprint_l2"),
                "structure_update_norm": record.get("structure_update_norm"),
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
            "best_total_loss_step": 0,
            "best_total_loss": 0.0,
            "position_difference_at_best_total_loss": 0.0,
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
    total_loss_values = np.asarray([entry["total_loss"] for entry in trace], dtype=np.float64)
    step_values = np.asarray([entry["step"] for entry in trace], dtype=np.int64)
    best_index = int(np.argmin(position_values))
    best_total_loss_index = int(np.argmin(total_loss_values))
    tail_count = min(3, len(trace))
    tail_positions = position_values[-tail_count:]
    tail_fingerprint = fingerprint_values[-tail_count:]
    fingerprint_deltas = np.diff(fingerprint_values)
    position_deltas = np.diff(position_values)

    return {
        "best_step": int(step_values[best_index]),
        "best_position_difference": float(position_values[best_index]),
        "final_position_difference": float(position_values[-1]),
        "best_total_loss_step": int(step_values[best_total_loss_index]),
        "best_total_loss": float(total_loss_values[best_total_loss_index]),
        "position_difference_at_best_total_loss": float(
            position_values[best_total_loss_index]
        ),
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


def collect_rollout_stage_samples(
    model: TorchGNNFingerprint,
    original,
    perturbed,
    fixed_atoms: np.ndarray,
    target_fingerprint: np.ndarray,
    inverse_steps: int,
    inverse_step_size: float,
    inverse_design_options: dict,
    stage_index: int,
    rollout_step_stride: int,
    rollout_drift_threshold: float,
    rollout_high_error_threshold: float,
    rollout_instability_threshold: float,
    ) -> tuple[list[ReplaySample], dict]:
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
        true_target_fingerprint_mse = float(np.mean((true_fingerprint - target_fingerprint) ** 2))
        surrogate_target_fingerprint_mse = float(
            np.mean((surrogate_fingerprint - target_fingerprint) ** 2)
        )
        fingerprint_drift_mse = float(np.mean((surrogate_fingerprint - true_fingerprint) ** 2))
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
        new_samples.append(
            ReplaySample(
                atoms=record["atoms"],
                priority=priority,
                categories=categories,
                stage_index=int(stage_index),
                step=int(record["step"]),
                restart_index=0,
                true_target_fingerprint_mse=true_target_fingerprint_mse,
                surrogate_target_fingerprint_mse=surrogate_target_fingerprint_mse,
                fingerprint_drift_mse=fingerprint_drift_mse,
                position_difference=position_difference,
                repulsion_loss=float(record["repulsion_loss"]),
                cell_violation_loss=float(record["cell_violation_loss"]),
                is_initial_state=bool(record["is_initial_state"]),
            )
        )
        previous_true_error = true_target_fingerprint_mse

    stage_true_errors = [sample.true_target_fingerprint_mse for sample in new_samples]
    stage_surrogate_errors = [sample.surrogate_target_fingerprint_mse for sample in new_samples]
    stage_drifts = [sample.fingerprint_drift_mse for sample in new_samples]
    stage_position_differences = [
        sample.position_difference for sample in new_samples if sample.position_difference is not None
    ]
    return new_samples, {
        "stage_index": int(stage_index),
        "num_stage_samples": int(len(new_samples)),
        "sampled_replay_count": 0,
        "buffer_size": 0,
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
        "mean_fingerprint_drift_mse": 0.0 if not stage_drifts else float(np.mean(stage_drifts)),
        "best_position_difference": (
            None if not stage_position_differences else float(np.min(stage_position_differences))
        ),
        "final_position_difference": (
            None if not stage_position_differences else float(stage_position_differences[-1])
        ),
        "training_losses": [],
        "trajectory": [
            {
                "step": int(sample.step),
                "is_initial_state": bool(sample.is_initial_state),
                "categories": list(sample.categories),
                "priority": float(sample.priority),
                "true_target_fingerprint_mse": float(sample.true_target_fingerprint_mse),
                "surrogate_target_fingerprint_mse": float(sample.surrogate_target_fingerprint_mse),
                "fingerprint_drift_mse": float(sample.fingerprint_drift_mse),
                "position_difference": (
                    None if sample.position_difference is None else float(sample.position_difference)
                ),
                "repulsion_loss": float(sample.repulsion_loss),
                "cell_violation_loss": float(sample.cell_violation_loss),
            }
            for sample in new_samples
        ],
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
    del carbon_structures
    del original
    del perturbed
    del fixed_atoms
    del target_fingerprint
    del batch_size
    del inverse_steps
    del inverse_step_size
    del inverse_design_options
    del seed
    del rollout_stages
    del rollout_epochs_per_stage
    del rollout_step_stride
    del replay_sample_size
    del replay_category_weights
    del rollout_drift_threshold
    del rollout_high_error_threshold
    del rollout_instability_threshold
    del training_observer
    del training_epoch_offset
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


def resolve_effective_inverse_steps(
    inverse_steps: int,
    step_values: list[int] | None,
) -> int:
    candidates = [int(inverse_steps)]
    if step_values is not None:
        candidates.extend(int(step) for step in step_values)
    return max((max(step, 0) for step in candidates), default=0)


def select_logged_inverse_records(
    trajectory_records: list[dict],
    *,
    initial_atoms,
    final_atoms,
    num_steps: int,
    step_size: float,
    step_interval: int = INVERSE_STEP_LOG_INTERVAL,
) -> list[dict]:
    normalized_records = [
        {
            **record,
            "atoms": record["atoms"].copy(),
            "step": int(record.get("step", 0)),
            "num_steps": int(record.get("num_steps", num_steps)),
            "is_initial_state": bool(record.get("is_initial_state", False)),
            "learning_rate": float(record.get("learning_rate", step_size)),
            "total_loss": float(record.get("total_loss", 0.0)),
            "fingerprint_loss": float(record.get("fingerprint_loss", 0.0)),
            "repulsion_loss": float(record.get("repulsion_loss", 0.0)),
            "cell_violation_loss": float(record.get("cell_violation_loss", 0.0)),
        }
        for record in trajectory_records
    ]

    if not any(bool(record["is_initial_state"]) for record in normalized_records):
        normalized_records.insert(
            0,
            {
                "atoms": initial_atoms.copy(),
                "step": 0,
                "num_steps": int(num_steps),
                "is_initial_state": True,
                "learning_rate": float(step_size),
                "total_loss": 0.0,
                "fingerprint_loss": 0.0,
                "repulsion_loss": 0.0,
                "cell_violation_loss": 0.0,
            },
        )
    if not any(int(record["step"]) == int(num_steps) for record in normalized_records):
        normalized_records.append(
            {
                "atoms": final_atoms.copy(),
                "step": int(num_steps),
                "num_steps": int(num_steps),
                "is_initial_state": False,
                "learning_rate": float(step_size),
                "total_loss": 0.0,
                "fingerprint_loss": 0.0,
                "repulsion_loss": 0.0,
                "cell_violation_loss": 0.0,
            }
        )

    initial_record = None
    records_by_step: dict[int, dict] = {}
    for record in sorted(
        normalized_records,
        key=lambda value: (int(value["step"]), 0 if bool(value["is_initial_state"]) else 1),
    ):
        step = int(record["step"])
        if bool(record["is_initial_state"]):
            initial_record = record
            continue
        if (
            int(step_interval) <= 1
            or step == int(num_steps)
            or (step > 0 and step % int(step_interval) == 0)
        ):
            records_by_step[step] = record

    filtered_records: list[dict] = []
    if initial_record is not None:
        filtered_records.append(initial_record)
    filtered_records.extend(records_by_step[step] for step in sorted(records_by_step))
    return filtered_records


def annotate_inverse_rollout_records(
    model: TorchGNNFingerprint,
    trajectory_records: list[dict],
    original,
    perturbed,
    target_fingerprint: np.ndarray,
) -> tuple[list[dict], list[dict[str, float]]]:
    annotated_records: list[dict] = []
    trace: list[dict[str, float]] = []
    for record in trajectory_records:
        candidate_atoms = record["atoms"].copy()
        position_difference = score_candidate(original, candidate_atoms)
        predicted_fingerprint = model.predict(candidate_atoms)
        fingerprint_mse = float(np.mean((predicted_fingerprint - target_fingerprint) ** 2))
        fingerprint_l2 = float(np.linalg.norm(predicted_fingerprint - target_fingerprint))
        structure_update_norm = float(
            np.linalg.norm(candidate_atoms.get_positions() - perturbed.get_positions())
        )
        annotated_record = {
            **record,
            "atoms": candidate_atoms,
            "position_difference": position_difference,
            "fingerprint_mse": fingerprint_mse,
            "fingerprint_l2": fingerprint_l2,
            "structure_update_norm": structure_update_norm,
        }
        annotated_records.append(annotated_record)
        trace.append(
            {
                "step": int(annotated_record["step"]),
                "num_steps": int(annotated_record["num_steps"]),
                "is_initial_state": bool(annotated_record["is_initial_state"]),
                "position_difference": position_difference,
                "fingerprint_mse": fingerprint_mse,
                "fingerprint_l2": fingerprint_l2,
                "structure_update_norm": structure_update_norm,
                "learning_rate": float(annotated_record.get("learning_rate", 0.0)),
                "total_loss": float(annotated_record.get("total_loss", 0.0)),
                "fingerprint_loss": float(annotated_record.get("fingerprint_loss", 0.0)),
                "repulsion_loss": float(annotated_record.get("repulsion_loss", 0.0)),
                "cell_violation_loss": float(annotated_record.get("cell_violation_loss", 0.0)),
            }
        )
    return annotated_records, trace


def inverse_design_trace(
    model: TorchGNNFingerprint,
    original,
    perturbed,
    fixed_atoms: np.ndarray,
    target_fingerprint: np.ndarray,
    inverse_steps: int,
    step_values: list[int],
    step_size: float,
    inverse_design_options: dict,
    step_interval: int = INVERSE_STEP_LOG_INTERVAL,
) -> tuple:
    resolved_inverse_steps = resolve_effective_inverse_steps(inverse_steps, step_values)
    trajectory_records: list[dict] = []
    step_observer = build_rollout_step_observer(trajectory_records)
    optimised = model.inverse_design(
        target_fingerprint=target_fingerprint,
        atoms=perturbed,
        fixed_atoms=fixed_atoms,
        num_steps=resolved_inverse_steps,
        step_size=step_size,
        verbose=0,
        step_observer=step_observer,
        **inverse_design_options,
    )
    logged_records = select_logged_inverse_records(
        trajectory_records,
        initial_atoms=perturbed,
        final_atoms=optimised,
        num_steps=resolved_inverse_steps,
        step_size=step_size,
        step_interval=step_interval,
    )
    annotated_records, trace = annotate_inverse_rollout_records(
        model=model,
        trajectory_records=logged_records,
        original=original,
        perturbed=perturbed,
        target_fingerprint=target_fingerprint,
    )
    return optimised, annotated_records, trace


def compute_rollout_stage_start_epochs(total_epochs: int, rollout_stages: int) -> list[int]:
    if int(total_epochs) <= 0 or int(rollout_stages) <= 0:
        return []
    starts = []
    for stage_index in range(int(rollout_stages)):
        start_epoch = int(round(((stage_index + 1) * int(total_epochs)) / (int(rollout_stages) + 1)))
        start_epoch = min(max(start_epoch, 1), int(total_epochs))
        if start_epoch not in starts:
            starts.append(start_epoch)
    return starts


def sweep_epochs(
    carbon_structures,
    augmented_structures,
    original,
    perturbed,
    fixed_atoms: np.ndarray,
    num_epochs: int,
    batch_size: int,
    inverse_steps: int,
    inverse_step_size: float,
    inverse_design_options: dict,
    seed: int,
    model_config: Optional[dict] = None,
    training_observer: Optional[Callable[[int, float], None]] = None,
    rollout_stages: int = 0,
    rollout_epochs_per_stage: int = 0,
    rollout_step_stride: int = ROLLOUT_STEP_STRIDE,
    replay_buffer_capacity: int = REPLAY_BUFFER_CAPACITY,
    replay_sample_size: int = REPLAY_SAMPLE_SIZE,
    replay_category_weights: dict[str, float] | None = None,
    rollout_drift_threshold: float = ROLLOUT_DRIFT_THRESHOLD,
    rollout_high_error_threshold: float = ROLLOUT_HIGH_ERROR_THRESHOLD,
    rollout_instability_threshold: float = ROLLOUT_INSTABILITY_THRESHOLD,
) -> tuple:
    model = create_model(seed, model_config=model_config)
    learnable_parameter_count = sum(
        parameter.numel() for parameter in model.parameters() if parameter.requires_grad
    )
    target_fingerprint = model.compute_reference_fingerprint(original)
    replay_buffer = PrioritizedReplayBuffer(capacity=int(replay_buffer_capacity), seed=int(seed))
    training_losses = []
    epoch_results = []
    best_candidate = None
    checkpoint_records = []
    stage_summaries = []

    requested_epoch = max(int(num_epochs), 0)
    stage_start_epochs = compute_rollout_stage_start_epochs(requested_epoch, rollout_stages)
    stage_schedule = {
        epoch: stage_index + 1 for stage_index, epoch in enumerate(stage_start_epochs)
    }

    print(
        f"[train] learnable_parameters={learnable_parameter_count:,}"
    )

    print(
        f"[train] epochs={requested_epoch} "
        f"rollout_starts={stage_start_epochs if stage_start_epochs else 'disabled'}"
    )

    trained_epochs = 0
    active_stage_summary = None
    active_stage_epochs_remaining = 0

    if requested_epoch == 0:
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
        epoch_results.append({"epochs": 0, "position_difference": position_difference})
        print(f"[train] epoch checkpoint 0/0: rmsd={position_difference:.6f} A")

    while trained_epochs < requested_epoch:
        next_epoch = trained_epochs + 1
        if next_epoch in stage_schedule and int(rollout_epochs_per_stage) > 0:
            stage_index = stage_schedule[next_epoch]
            print(
                f"[train] rollout stage {stage_index}/{len(stage_schedule)} collecting samples "
                f"before epoch {next_epoch}"
            )
            new_samples, stage_summary = collect_rollout_stage_samples(
                model=model,
                original=original,
                perturbed=perturbed,
                fixed_atoms=fixed_atoms,
                target_fingerprint=target_fingerprint,
                inverse_steps=inverse_steps,
                inverse_step_size=inverse_step_size,
                inverse_design_options=inverse_design_options,
                stage_index=stage_index,
                rollout_step_stride=rollout_step_stride,
                rollout_drift_threshold=rollout_drift_threshold,
                rollout_high_error_threshold=rollout_high_error_threshold,
                rollout_instability_threshold=rollout_instability_threshold,
            )
            replay_buffer.extend(new_samples)
            stage_summary["buffer_size"] = int(len(replay_buffer))
            stage_summaries.append(stage_summary)
            active_stage_summary = stage_summary if new_samples else None
            active_stage_epochs_remaining = int(rollout_epochs_per_stage) if new_samples else 0
            print(
                f"[train] rollout stage {stage_index}: retained {len(new_samples)} samples, "
                f"buffer={len(replay_buffer)}"
            )

        augment_batch = list(augmented_structures)
        if active_stage_summary is not None and active_stage_epochs_remaining > 0:
            sampled_replay = replay_buffer.sample(
                sample_size=int(replay_sample_size),
                category_weights=(
                    dict(DEFAULT_REPLAY_CATEGORY_WEIGHTS)
                    if replay_category_weights is None
                    else dict(replay_category_weights)
                ),
            )
            augment_batch.extend(sample.atoms for sample in sampled_replay)
            active_stage_summary["sampled_replay_count"] += len(sampled_replay)

        history = model.fit(
            carbon_structures,
            num_epochs=1,
            batch_size=batch_size,
            augment_structures=augment_batch,
            verbose=0,
            reset_optimiser=(trained_epochs == 0),
            recalibrate_base=(trained_epochs == 0),
        )
        trained_epochs += 1
        epoch_loss = float(history[-1])
        training_losses.append(epoch_loss)
        if training_observer is not None:
            training_observer(trained_epochs, epoch_loss)
        if active_stage_summary is not None and active_stage_epochs_remaining > 0:
            active_stage_summary["training_losses"].append(epoch_loss)
            active_stage_epochs_remaining -= 1
            if active_stage_epochs_remaining == 0:
                active_stage_summary = None

        if trained_epochs == requested_epoch:
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
                epochs=int(trained_epochs),
            )
            checkpoint_records.append(capture_checkpoint(model, int(trained_epochs), position_difference))
            epoch_results.append(
                {"epochs": int(trained_epochs), "position_difference": position_difference}
            )
            print(
                f"[train] epoch checkpoint {trained_epochs}/{requested_epoch}: "
                f"rmsd={position_difference:.6f} A"
            )

    rollout_metrics = {
        "enabled": bool(stage_summaries),
        "num_rollout_epochs": int(sum(len(stage["training_losses"]) for stage in stage_summaries)),
        "replay_buffer": replay_buffer.stats(),
        "stages": stage_summaries,
    }
    return model, epoch_results, training_losses, best_candidate, checkpoint_records, rollout_metrics


def run_for_epochs(
    carbon_structures,
    augmented_structures,
    original,
    perturbed,
    fixed_atoms: np.ndarray,
    num_epochs: int,
    batch_size: int,
    inverse_steps: int,
    inverse_step_size: float,
    inverse_design_options: dict,
    seed: int,
    model_config: Optional[dict] = None,
    training_observer: Optional[Callable[[int, float], None]] = None,
) -> tuple:
    model, epoch_results, training_losses, best_candidate, checkpoint_records, _ = sweep_epochs(
        carbon_structures=carbon_structures,
        augmented_structures=augmented_structures,
        original=original,
        perturbed=perturbed,
        fixed_atoms=fixed_atoms,
        num_epochs=int(num_epochs),
        batch_size=batch_size,
        inverse_steps=inverse_steps,
        inverse_step_size=inverse_step_size,
        inverse_design_options=inverse_design_options,
        seed=seed,
        model_config=model_config,
        training_observer=training_observer,
    )
    return (
        model,
        epoch_results[-1] if epoch_results else {"epochs": int(num_epochs), "position_difference": 0.0},
        training_losses,
        best_candidate,
        checkpoint_records[-1] if checkpoint_records else capture_checkpoint(model, int(num_epochs), 0.0),
    )


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


def evaluate_checkpoints_with_step_size(
    model: TorchGNNFingerprint,
    checkpoint_records: list[dict],
    original,
    perturbed,
    fixed_atoms: np.ndarray,
    target_fingerprint: np.ndarray,
    inverse_steps: int,
    step_size: float,
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
                source="checkpoint_best_step_size",
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
    inverse_steps: int,
    inverse_step_size: float,
    inverse_design_options: dict,
    step_interval: int = INVERSE_STEP_LOG_INTERVAL,
) -> list[dict[str, float]]:
    if not checkpoint_records:
        return []

    results = []
    original_state = copy.deepcopy(model.state_dict())
    try:
        for checkpoint in checkpoint_records:
            model.load_state_dict(checkpoint["state_dict"])
            _, _, trace = inverse_design_trace(
                model=model,
                original=original,
                perturbed=perturbed,
                fixed_atoms=fixed_atoms,
                target_fingerprint=target_fingerprint,
                inverse_steps=inverse_steps,
                step_values=[int(inverse_steps)],
                step_size=inverse_step_size,
                inverse_design_options=inverse_design_options,
                step_interval=step_interval,
            )
            for entry in trace:
                results.append(
                    {
                        "epochs": int(checkpoint["epochs"]),
                        "num_steps": int(entry["step"]),
                        "step_size": float(inverse_step_size),
                        "position_difference": float(entry["position_difference"]),
                        "fingerprint_mse": float(entry["fingerprint_mse"]),
                        "fingerprint_l2": float(entry["fingerprint_l2"]),
                        "structure_update_norm": float(entry["structure_update_norm"]),
                        "total_loss": float(entry.get("total_loss", 0.0)),
                    }
                )
    finally:
        model.load_state_dict(original_state)

    return results


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
    figure = plt.figure(figsize=(10, 4.5))

    ax1 = figure.add_subplot(1, 2, 1)
    ax1.plot(
        [entry["step"] for entry in inverse_trace],
        [entry["position_difference"] for entry in inverse_trace],
        marker="o",
    )
    ax1.set_xlabel("Inverse-design steps")
    ax1.set_ylabel("Position difference to target (symmetry-aware RMSD, A)")
    ax1.set_title("Inverse-design rollout")
    ax1.axhline(initial_position_difference, color="tab:red", linestyle="--", linewidth=1.0)
    ax1.grid(alpha=0.3)

    ax2 = figure.add_subplot(1, 2, 2)
    ax2.plot(
        [entry["step_size"] for entry in step_size_results],
        [entry["position_difference"] for entry in step_size_results],
        marker="o",
    )
    ax2.set_xscale("log")
    ax2.set_xlabel("Step size")
    ax2.set_ylabel("Position difference to target (symmetry-aware RMSD, A)")
    ax2.set_title("Step-size sweep")
    ax2.axhline(initial_position_difference, color="tab:red", linestyle="--", linewidth=1.0)
    ax2.grid(alpha=0.3)

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
    inverse_step_values: list[int],
    step_size_values: list[float],
    fixed_leading_atoms: int,
    fingerprint_loss_weight: float,
    target_vertex_weight: float,
    target_position_weight: float,
    inverse_lr_decay_rate: float,
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
    inverse_restarts: int | None = None,
    inverse_restart_noise_scale: float | None = None,
    model_config: Optional[dict] = None,
    training_observer: Optional[Callable[[int, float], None]] = None,
    architecture_name: str = "torch_gnn_residual",
    enable_checkpoint_step_size_sweep: bool = True,
    enable_checkpoint_step_schedule_sweep: bool = True,
    ensemble_enabled: bool = False,
    ensemble_num_trajectories: int = 16,
    ensemble_perturbation_scale: float = 0.01,
    ensemble_langevin_noise_scale: float = 0.005,
    ensemble_temperature: float = 1.0,
    ensemble_adaptive_scaling: bool = True,
    ensemble_trajectory_pruning: bool = True,
    ensemble_escape_detection: bool = True,
    ensemble_consensus_metric: str = "cluster",
    ensemble_aggregation: str = "mean_variance",
) -> dict:
    if int(augmented_count) != 0:
        raise ValueError("augmented_count must remain 0 for inverse-design runs")
    if float(target_vertex_weight) != 0.0:
        raise ValueError("target_vertex_weight must remain 0.0 for inverse-design runs")
    if float(target_position_weight) != 0.0:
        raise ValueError("target_position_weight must remain 0.0 for inverse-design runs")

    print("Loading carbon structures...")
    carbon_xyz = repo_root / "example" / "data" / "carbon.xyz"
    all_carbon_structures = read(str(carbon_xyz), index=":")
    carbon_structures = select_carbon_structures(all_carbon_structures, carbon_count)
    minimum_carbon_count = minimum_training_carbon_count(len(all_carbon_structures))
    if len(carbon_structures) < minimum_carbon_count:
        raise ValueError(
            "carbon_count must select at least 30% of example/data/carbon.xyz "
            f"structures ({minimum_carbon_count} minimum, received {len(carbon_structures)})"
        )

    print("Setting up original and perturbed structures...")
    original = bulk("C", "diamond", a=3.567, cubic=True)
    original.pbc = True
    augmented_structures = []
    perturbed, fixed_atoms = build_perturbed_structure(original, fixed_leading_atoms=fixed_leading_atoms)
    inverse_design_options = build_inverse_design_options(
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

    print(f"[workflow] training {int(num_epochs)} epochs")
    model, epoch_results, training_losses, best_epoch_candidate, checkpoint_records, rollout_metrics = sweep_epochs(
        carbon_structures=carbon_structures,
        augmented_structures=augmented_structures,
        original=original,
        perturbed=perturbed,
        fixed_atoms=fixed_atoms,
        num_epochs=int(num_epochs),
        batch_size=batch_size,
        inverse_steps=inverse_steps,
        inverse_step_size=inverse_step_size,
        inverse_design_options=inverse_design_options,
        seed=seed,
        model_config=model_config,
        training_observer=training_observer,
        rollout_stages=rollout_stages,
        rollout_epochs_per_stage=rollout_epochs_per_stage,
        rollout_step_stride=rollout_step_stride,
        replay_buffer_capacity=replay_buffer_capacity,
        replay_sample_size=replay_sample_size,
        replay_category_weights=replay_category_weights,
        rollout_drift_threshold=rollout_drift_threshold,
        rollout_high_error_threshold=rollout_high_error_threshold,
        rollout_instability_threshold=rollout_instability_threshold,
    )

    target_fingerprint = model.compute_reference_fingerprint(original)

    # Initialize ensemble exploration if enabled
    ensemble_instance = None
    ensemble_stats_collector = []
    if ensemble_enabled:
        if InferenceEnsemble is None or InferenceEnsembleConfig is None:
            raise RuntimeError(
                "Ensemble inverse-design support is unavailable in this raffle installation"
            )
        ensemble_config = build_ensemble_config(
            ensemble_enabled=True,
            ensemble_num_trajectories=ensemble_num_trajectories,
            ensemble_perturbation_scale=ensemble_perturbation_scale,
            ensemble_langevin_noise_scale=ensemble_langevin_noise_scale,
            ensemble_temperature=ensemble_temperature,
            ensemble_adaptive_scaling=ensemble_adaptive_scaling,
            ensemble_trajectory_pruning=ensemble_trajectory_pruning,
            ensemble_escape_detection=ensemble_escape_detection,
            ensemble_consensus_metric=ensemble_consensus_metric,
            ensemble_aggregation=ensemble_aggregation,
        )
        ensemble_instance = InferenceEnsemble(
            model,
            InferenceEnsembleConfig(**ensemble_config)
        )

        # Create wrapper for inverse_design that uses ensemble
        original_inverse_design = model.inverse_design
        def ensemble_wrapped_inverse_design(**kwargs):
            num_steps = int(kwargs.get('num_steps', 100))
            step_size = float(kwargs.get('step_size', 1.0e-3))
            external_step_observer = kwargs.get('step_observer')
            source_atoms = kwargs.get('atoms')

            # Extract needed params for ensemble
            optimised, stats = ensemble_instance.run_ensemble_inverse_design(
                target_fingerprint=kwargs.get('target_fingerprint'),
                atoms=kwargs.get('atoms'),
                fixed_atoms=kwargs.get('fixed_atoms'),
                num_steps=num_steps,
                step_size=step_size,
                fingerprint_loss_weight=kwargs.get('fingerprint_loss_weight', 1.0),
                repulsion_weight=kwargs.get('repulsion_weight', 10.0),
                minimum_distance_scale=kwargs.get('minimum_distance_scale', 0.75),
                cell_violation_weight=kwargs.get('cell_violation_weight', 0.0),
                seed=int(np.random.randint(0, 2**31 - 1)),
                step_observer=None,
            )

            # Preserve legacy observer contract expected by rollout/traj consumers.
            if external_step_observer is not None:
                best_record = next(
                    (
                        record
                        for record in stats.trajectory_records
                        if int(record.trajectory_id) == int(stats.best_trajectory_id)
                    ),
                    None,
                )
                if best_record is not None and source_atoms is not None:
                    initial_atoms = source_atoms.copy()
                    initial_atoms.set_positions(best_record.initial_positions)
                    external_step_observer(
                        {
                            'atoms': initial_atoms,
                            'step': 0,
                            'num_steps': num_steps,
                            'is_initial_state': True,
                            'learning_rate': float(best_record.step_size or step_size),
                            'total_loss': 0.0,
                            'fingerprint_loss': 0.0,
                            'repulsion_loss': 0.0,
                            'cell_violation_loss': 0.0,
                        }
                    )
                    for step, positions in best_record.sampled_positions:
                        observed_atoms = source_atoms.copy()
                        observed_atoms.set_positions(positions)
                        loss_index = max(min(int(step) - 1, len(best_record.loss_history) - 1), 0)
                        total_loss = (
                            float(best_record.loss_history[loss_index])
                            if best_record.loss_history
                            else float(stats.best_loss)
                        )
                        external_step_observer(
                            {
                                'atoms': observed_atoms,
                                'step': int(step),
                                'num_steps': num_steps,
                                'is_initial_state': False,
                                'learning_rate': float(best_record.step_size or step_size),
                                'total_loss': total_loss,
                                'fingerprint_loss': total_loss,
                                'repulsion_loss': 0.0,
                                'cell_violation_loss': 0.0,
                            }
                        )
                else:
                    external_step_observer(
                        {
                            'atoms': optimised.copy(),
                            'step': num_steps,
                            'num_steps': num_steps,
                            'is_initial_state': False,
                            'learning_rate': step_size,
                            'total_loss': float(stats.best_loss),
                            'fingerprint_loss': float(stats.mean_loss),
                            'repulsion_loss': 0.0,
                            'cell_violation_loss': 0.0,
                        }
                    )

            # Store stats for logging
            ensemble_stats_collector.append(stats)
            return optimised

        # Monkey-patch the inverse_design method
        model.inverse_design = ensemble_wrapped_inverse_design

    resolved_inverse_steps = resolve_effective_inverse_steps(
        inverse_steps=inverse_steps,
        step_values=inverse_step_values,
    )
    print(
        f"[eval] running inverse-design rollout trace to {resolved_inverse_steps} steps "
        f"with logging every {INVERSE_STEP_LOG_INTERVAL} steps"
    )
    configured_optimised, configured_trajectory_records, inverse_trace = inverse_design_trace(
        model=model,
        original=original,
        perturbed=perturbed,
        fixed_atoms=fixed_atoms,
        target_fingerprint=target_fingerprint,
        inverse_steps=resolved_inverse_steps,
        step_values=inverse_step_values,
        step_size=inverse_step_size,
        inverse_design_options=inverse_design_options,
    )
    convergence_summary = summarise_convergence(inverse_trace)
    print("[eval] running inverse-design step-size sweep")
    step_size_results, best_step_size_candidate = sweep_step_sizes(
        model=model,
        original=original,
        perturbed=perturbed,
        fixed_atoms=fixed_atoms,
        target_fingerprint=target_fingerprint,
        inverse_steps=resolved_inverse_steps,
        step_sizes=step_size_values,
        inverse_design_options=inverse_design_options,
    )
    best_step_size = (
        float(best_step_size_candidate["step_size"])
        if best_step_size_candidate is not None
        else float(inverse_step_size)
    )
    checkpoint_step_size_results = []
    best_checkpoint_step_size_candidate = None
    if enable_checkpoint_step_size_sweep:
        print(
            f"[eval] collecting checkpoint results at best step size {best_step_size:.6g}"
        )
        checkpoint_step_size_results, best_checkpoint_step_size_candidate = evaluate_checkpoints_with_step_size(
            model=model,
            checkpoint_records=checkpoint_records,
            original=original,
            perturbed=perturbed,
            fixed_atoms=fixed_atoms,
            target_fingerprint=target_fingerprint,
            inverse_steps=resolved_inverse_steps,
            step_size=best_step_size,
            inverse_design_options=inverse_design_options,
        )
    checkpoint_step_schedule_results = []
    if enable_checkpoint_step_schedule_sweep:
        print("[eval] collecting checkpoint rollout traces")
        checkpoint_step_schedule_results = sweep_step_schedules_for_checkpoints(
            model=model,
            checkpoint_records=checkpoint_records,
            original=original,
            perturbed=perturbed,
            fixed_atoms=fixed_atoms,
            target_fingerprint=target_fingerprint,
            inverse_steps=resolved_inverse_steps,
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
        num_steps=int(resolved_inverse_steps),
        step_size=float(inverse_step_size),
    )
    for candidate in (
        best_epoch_candidate,
        best_step_size_candidate,
        best_checkpoint_step_size_candidate,
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
    checkpoint_path = output_dir / "torch_gnn_model_checkpoint.pt"
    fingerprint_path = output_dir / "torch_gnn_target_fingerprint.npy"
    reference_structure_path = output_dir / "torch_gnn_reference_structure.xyz"
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
        "position_error_metric": "structure_similarity_rmsd",
        "architecture_name": str(architecture_name),
        "model_config": dict(DEFAULT_MODEL_CONFIG if model_config is None else model_config),
        "inverse_design_config": {
            "inverse_steps": int(inverse_steps),
            "resolved_inverse_steps": int(resolved_inverse_steps),
            "inverse_step_size": float(inverse_step_size),
            "best_step_size": float(best_step_size),
            "inverse_step_log_interval": int(INVERSE_STEP_LOG_INTERVAL),
            "fixed_leading_atoms": int(fixed_leading_atoms),
            "fingerprint_loss_weight": float(fingerprint_loss_weight),
            "target_vertex_weight": float(target_vertex_weight),
            "target_position_weight": float(target_position_weight),
            "inverse_lr_decay_rate": float(inverse_lr_decay_rate),
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
            "checkpoint": str(checkpoint_path),
            "target_fingerprint": str(fingerprint_path),
            "reference_structure": str(reference_structure_path),
            "configured_inverse_design_traj": str(
                configured_inverse_design_path["traj_file"]
            ),
            "final_descriptor_comparison": descriptor_report["report_file"],
            "final_descriptor_comparison_plot": descriptor_report["plot_file"],
        },
        "descriptor_comparisons": {"final": descriptor_report},
        "ensemble_exploration": {
            "enabled": bool(ensemble_enabled),
            "num_trajectories": int(ensemble_num_trajectories) if ensemble_enabled else 0,
            "statistics_collected": len(ensemble_stats_collector),
            "final_stats": (
                {
                    "mean_loss": float(ensemble_stats_collector[-1].mean_loss) if ensemble_stats_collector else None,
                    "position_spread": float(ensemble_stats_collector[-1].position_spread) if ensemble_stats_collector else None,
                    "consensus_strength": float(ensemble_stats_collector[-1].consensus_strength) if ensemble_stats_collector else None,
                    "escape_count": int(ensemble_stats_collector[-1].escape_count) if ensemble_stats_collector else None,
                    "num_clusters": int(ensemble_stats_collector[-1].num_clusters) if ensemble_stats_collector else None,
                    "collapse_detected": bool(ensemble_stats_collector[-1].collapse_detected) if ensemble_stats_collector else None,
                }
                if ensemble_stats_collector
                else {}
            ),
        } if ensemble_enabled else {},
    }
    checkpoint_training_config = {
        "repo_root": str(repo_root),
        "seed": int(seed),
        "architecture_name": str(architecture_name),
        "workflow_config": {
            "carbon_count": int(carbon_count),
            "augmented_count": int(augmented_count),
            "epochs": int(num_epochs),
            "batch_size": int(batch_size),
            "inverse_steps": int(inverse_steps),
            "inverse_step_size": float(inverse_step_size),
            "fixed_leading_atoms": int(fixed_leading_atoms),
            "fingerprint_loss_weight": float(fingerprint_loss_weight),
            "target_vertex_weight": float(target_vertex_weight),
            "target_position_weight": float(target_position_weight),
            "inverse_lr_decay_rate": float(inverse_lr_decay_rate),
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
            "ignored_deprecated_config_fields": {
                key: value
                for key, value in {
                    "inverse_restarts": (
                        None if inverse_restarts is None else int(inverse_restarts)
                    ),
                    "inverse_restart_noise_scale": (
                        None
                        if inverse_restart_noise_scale is None
                        else float(inverse_restart_noise_scale)
                    ),
                }.items()
                if value is not None
            },
        },
        "num_available_carbon_structures": len(all_carbon_structures),
        "num_training_carbon_structures": len(carbon_structures),
        "convergence_summary": convergence_summary,
        "selected_candidate": metrics["selected_candidate"],
        "rollout": rollout_metrics,
    }
    save_model_checkpoint(
        model=model,
        checkpoint_path=checkpoint_path,
        species_list=["C"],
        model_config=metrics["model_config"],
        training_config=checkpoint_training_config,
        training_history=training_losses,
    )
    save_target_fingerprint(target_fingerprint, fingerprint_path)
    write_structure(reference_structure_path, original)
    metrics_path = output_dir / "torch_gnn_carbon_workflow_metrics.json"
    metrics["output_files"]["metrics"] = str(metrics_path)
    metrics_path.write_text(json.dumps(metrics, indent=2))
    return metrics


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--carbon-count", type=int, default=-1)
    parser.add_argument("--augmented-count", type=int, default=0)
    parser.add_argument("--epochs", type=int, default=30)
    parser.add_argument("--batch-size", type=int, default=8)
    parser.add_argument("--inverse-steps", type=int, default=400)
    parser.add_argument("--inverse-step-size", type=float, default=5.0e-3)
    parser.add_argument("--inverse-step-values", type=str, default="")
    parser.add_argument("--step-size-values", type=str, default="")
    parser.add_argument("--fixed-leading-atoms", type=int, default=FIXED_LEADING_ATOMS)
    parser.add_argument("--fingerprint-loss-weight", type=float, default=FINGERPRINT_LOSS_WEIGHT)
    parser.add_argument(
        "--target-vertex-weight",
        type=float,
        default=TARGET_VERTEX_WEIGHT,
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--target-position-weight",
        type=float,
        default=TARGET_POSITION_WEIGHT,
        help=argparse.SUPPRESS,
    )
    parser.add_argument("--inverse-lr-decay-rate", type=float, default=INVERSE_LR_DECAY_RATE)
    parser.add_argument(
        "--inverse-restarts",
        type=int,
        default=INVERSE_RESTARTS,
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--inverse-restart-noise-scale",
        type=float,
        default=INVERSE_RESTART_NOISE_SCALE,
        help=argparse.SUPPRESS,
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
    inverse_step_values = (
        parse_int_list(args.inverse_step_values)
        if args.inverse_step_values
        else default_inverse_step_values(args.inverse_steps)
    )
    resolved_inverse_steps = resolve_effective_inverse_steps(
        inverse_steps=args.inverse_steps,
        step_values=inverse_step_values,
    )
    step_size_values = (
        parse_float_list(args.step_size_values)
        if args.step_size_values
        else default_step_sizes(args.inverse_step_size)
    )
    print(
        "[config] training "
        f"architecture={args.architecture} epochs={args.epochs} batch_size={args.batch_size} "
        f"hidden_dim={args.hidden_dim} message_layers={args.num_message_layers}"
    )
    print(
        "[config] inverse design "
        f"steps={resolved_inverse_steps} step_size={args.inverse_step_size} "
        f"log_interval={INVERSE_STEP_LOG_INTERVAL} "
        f"fixed_leading_atoms={args.fixed_leading_atoms}"
    )
    print(
        "[config] rollout "
        f"stages={args.rollout_stages} epochs_per_stage={args.rollout_epochs_per_stage} "
        f"step_stride={args.rollout_step_stride} replay_capacity={args.replay_buffer_capacity} "
        f"replay_sample_size={args.replay_sample_size}"
    )
    print(
        "[config] evaluation "
        f"step_size_values={step_size_values}"
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
    print(
        json.dumps(
            {
                "initial_rmsd": metrics["initial_rmsd"],
                "final_rmsd": metrics["final_rmsd"],
                "rmsd_reduction_fraction": metrics["rmsd_reduction_fraction"],
                "metrics_file": metrics["output_files"]["metrics"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
