"""Shared helpers for the PyTorch GNN training and inverse-design examples."""

from __future__ import annotations

import copy
import json
from pathlib import Path
from typing import Any, Optional, Sequence

import numpy as np
import torch
from ase.build import bulk
from ase.io import read, write

from raffle import (
    TorchGNNFingerprint,
    minimum_image_displacements,
    symmetry_aware_displacements,
    symmetry_aware_rmsd,
)


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


def ensure_torch_gnn_available() -> None:
    if TorchGNNFingerprint is None:
        raise RuntimeError("TorchGNNFingerprint is unavailable in this installation")


def default_reference_structure():
    reference = bulk("C", "diamond", a=3.567, cubic=True)
    reference.pbc = True
    return reference


def read_json_config(config_path: Path) -> dict[str, Any]:
    payload = json.loads(config_path.read_text())
    if not isinstance(payload, dict):
        raise ValueError(f"Configuration file must contain a JSON object: {config_path}")
    return payload


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True))


def fingerprint_component_dimensions(model) -> dict[str, int]:
    if all(
        hasattr(model, attribute)
        for attribute in ("fingerprint_dim_2body", "fingerprint_dim_3body", "fingerprint_dim_4body")
    ):
        return {
            "2body": int(model.fingerprint_dim_2body),
            "3body": int(model.fingerprint_dim_3body),
            "4body": int(model.fingerprint_dim_4body),
        }
    if hasattr(model, "component_dims"):
        dims = tuple(int(value) for value in model.component_dims)
        if len(dims) != 3:
            raise ValueError("model.component_dims must contain exactly three entries")
        return {"2body": dims[0], "3body": dims[1], "4body": dims[2]}
    raise AttributeError("Model does not expose fingerprint component dimensions")


def fingerprint_component_ranges(model) -> list[tuple[str, int, int]]:
    start = 0
    ranges: list[tuple[str, int, int]] = []
    for label, width in fingerprint_component_dimensions(model).items():
        end = start + int(width)
        ranges.append((label, start, end))
        start = end
    return ranges


def split_fingerprint_components(model, fingerprint: np.ndarray) -> dict[str, np.ndarray]:
    array = np.asarray(fingerprint, dtype=np.float32).reshape(-1)
    components: dict[str, np.ndarray] = {}
    expected_size = 0
    for label, start, end in fingerprint_component_ranges(model):
        components[label] = array[start:end].astype(np.float32, copy=False)
        expected_size = end
    if array.size != expected_size:
        raise ValueError(
            "Fingerprint length does not match model component dimensions: "
            f"{array.size} != {expected_size}"
        )
    return components


def _resolve_true_fingerprint(model, atoms) -> np.ndarray:
    if hasattr(model, "compute_reference_fingerprint"):
        values = model.compute_reference_fingerprint(atoms)
    elif hasattr(model, "compute_fingerprint"):
        values = model.compute_fingerprint(atoms)
    else:
        raise AttributeError("Model does not expose a true RAFFLE fingerprint API")
    return np.asarray(values, dtype=np.float32).reshape(-1)


def _resolve_true_fingerprint_components(model, atoms) -> dict[str, np.ndarray]:
    if hasattr(model, "compute_reference_components"):
        values = model.compute_reference_components(atoms)
        return {
            "2body": np.asarray(values[0], dtype=np.float32).reshape(-1),
            "3body": np.asarray(values[1], dtype=np.float32).reshape(-1),
            "4body": np.asarray(values[2], dtype=np.float32).reshape(-1),
        }
    if hasattr(model, "compute_fingerprint_components"):
        values = model.compute_fingerprint_components(atoms)
        return {
            "2body": np.asarray(values[0], dtype=np.float32).reshape(-1),
            "3body": np.asarray(values[1], dtype=np.float32).reshape(-1),
            "4body": np.asarray(values[2], dtype=np.float32).reshape(-1),
        }
    return split_fingerprint_components(model, _resolve_true_fingerprint(model, atoms))


def _resolve_predicted_fingerprint(model, atoms) -> np.ndarray:
    if not hasattr(model, "predict"):
        raise AttributeError("Model does not expose a predicted fingerprint API")
    return np.asarray(model.predict(atoms), dtype=np.float32).reshape(-1)


def _resolve_predicted_fingerprint_components(model, atoms) -> dict[str, np.ndarray]:
    if hasattr(model, "predict_components"):
        values = model.predict_components(atoms)
        return {
            "2body": np.asarray(values[0], dtype=np.float32).reshape(-1),
            "3body": np.asarray(values[1], dtype=np.float32).reshape(-1),
            "4body": np.asarray(values[2], dtype=np.float32).reshape(-1),
        }
    return split_fingerprint_components(model, _resolve_predicted_fingerprint(model, atoms))


def fingerprint_error_summary(reference: np.ndarray, candidate: np.ndarray) -> dict[str, float]:
    reference_array = np.asarray(reference, dtype=np.float32).reshape(-1)
    candidate_array = np.asarray(candidate, dtype=np.float32).reshape(-1)
    if reference_array.shape != candidate_array.shape:
        raise ValueError(
            "Fingerprint arrays must have the same shape: "
            f"{reference_array.shape} != {candidate_array.shape}"
        )
    delta = candidate_array - reference_array
    mse = float(np.mean(delta ** 2)) if delta.size else 0.0
    return {
        "mae": float(np.mean(np.abs(delta))) if delta.size else 0.0,
        "rmse": float(np.sqrt(mse)),
        "mse": mse,
        "l2": float(np.linalg.norm(delta)),
        "max_abs_error": float(np.max(np.abs(delta))) if delta.size else 0.0,
    }


def descriptor_report_path(structure_path: Path) -> Path:
    return structure_path.with_name(f"{structure_path.stem}_descriptor_comparison.json")


def save_descriptor_comparison_report(
    model,
    target_fingerprint: np.ndarray,
    final_atoms,
    structure_path: Path,
    structure_label: str = "final",
) -> dict[str, Any]:
    resolved_structure_path = structure_path.resolve()
    target_global = np.asarray(target_fingerprint, dtype=np.float32).reshape(-1)
    predicted_global = _resolve_predicted_fingerprint(model, final_atoms)
    true_global = _resolve_true_fingerprint(model, final_atoms)
    if predicted_global.shape != target_global.shape:
        raise ValueError(
            "Predicted and target fingerprints must have the same shape: "
            f"{predicted_global.shape} != {target_global.shape}"
        )
    if true_global.shape != target_global.shape:
        raise ValueError(
            "True and target fingerprints must have the same shape: "
            f"{true_global.shape} != {target_global.shape}"
        )

    target_components = split_fingerprint_components(model, target_global)
    predicted_components = _resolve_predicted_fingerprint_components(model, final_atoms)
    true_components = _resolve_true_fingerprint_components(model, final_atoms)

    global_error_vectors = {
        "ml_predicted_minus_target": (predicted_global - target_global).astype(np.float32),
        "true_raffle_minus_target": (true_global - target_global).astype(np.float32),
        "ml_predicted_minus_true_raffle": (predicted_global - true_global).astype(np.float32),
    }
    component_reports: dict[str, Any] = {}
    component_summaries: dict[str, Any] = {}
    for component_name in ("2body", "3body", "4body"):
        target_component = target_components[component_name]
        predicted_component = predicted_components[component_name]
        true_component = true_components[component_name]
        component_metrics = {
            "ml_predicted_vs_target": fingerprint_error_summary(
                target_component,
                predicted_component,
            ),
            "true_raffle_vs_target": fingerprint_error_summary(
                target_component,
                true_component,
            ),
            "ml_predicted_vs_true_raffle": fingerprint_error_summary(
                true_component,
                predicted_component,
            ),
        }
        component_reports[component_name] = {
            "target_raffle": target_component.tolist(),
            "ml_predicted": predicted_component.tolist(),
            "true_raffle": true_component.tolist(),
            "error_vectors": {
                "ml_predicted_minus_target": (predicted_component - target_component).tolist(),
                "true_raffle_minus_target": (true_component - target_component).tolist(),
                "ml_predicted_minus_true_raffle": (predicted_component - true_component).tolist(),
            },
            "metrics": component_metrics,
        }
        component_summaries[component_name] = component_metrics

    report = {
        "structure_label": str(structure_label),
        "structure_file": str(resolved_structure_path),
        "fingerprint_length": int(target_global.size),
        "component_dimensions": fingerprint_component_dimensions(model),
        "global_fingerprints": {
            "target_raffle": target_global.tolist(),
            "ml_predicted": predicted_global.tolist(),
            "true_raffle": true_global.tolist(),
        },
        "global_error_vectors": {
            key: value.tolist() for key, value in global_error_vectors.items()
        },
        "global_metrics": {
            "ml_predicted_vs_target": fingerprint_error_summary(target_global, predicted_global),
            "true_raffle_vs_target": fingerprint_error_summary(target_global, true_global),
            "ml_predicted_vs_true_raffle": fingerprint_error_summary(true_global, predicted_global),
        },
        "components": component_reports,
        "component_summaries": component_summaries,
    }
    report_path = descriptor_report_path(resolved_structure_path)
    write_json(report_path, report)
    report["report_file"] = str(report_path)
    return report


def resolve_path(value: Optional[str | Path], base_dir: Path) -> Optional[Path]:
    if value is None:
        return None
    path = Path(value)
    if path.is_absolute():
        return path
    return (base_dir / path).resolve()


def build_augmented_structures(
    original,
    count: int,
    seed: int,
    noise_scale: float = 0.04,
) -> list:
    rng = np.random.default_rng(seed)
    structures = []
    for _ in range(max(int(count), 0)):
        atoms = original.copy()
        atoms.set_positions(
            atoms.get_positions()
            + rng.normal(scale=float(noise_scale), size=atoms.positions.shape).astype(np.float32)
        )
        structures.append(atoms)
    return structures


def select_structures(structures: list, requested_count: int) -> list:
    if int(requested_count) <= 0 or int(requested_count) >= len(structures):
        return list(structures)
    return list(structures[: int(requested_count)])


def infer_species_list(structures: Sequence) -> list[str]:
    species_list: list[str] = []
    seen: set[str] = set()
    for atoms in structures:
        for symbol in atoms.get_chemical_symbols():
            resolved = str(symbol).strip()
            if resolved in seen:
                continue
            seen.add(resolved)
            species_list.append(resolved)
    if not species_list:
        raise ValueError("Unable to infer species list from an empty structure collection")
    return species_list


def normalise_component_weight(value: Sequence[float]) -> tuple[float, float, float]:
    if len(value) != 3:
        raise ValueError("component_weight must contain exactly three entries")
    return tuple(float(entry) for entry in value)


def normalise_model_config(model_config: Optional[dict[str, Any]] = None) -> dict[str, Any]:
    resolved = dict(DEFAULT_MODEL_CONFIG)
    if model_config:
        resolved.update(model_config)
    resolved["component_weight"] = normalise_component_weight(
        resolved.get("component_weight", DEFAULT_COMPONENT_WEIGHT)
    )
    return resolved


def create_model(
    seed: int,
    species_list: Sequence[str],
    model_config: Optional[dict[str, Any]] = None,
) -> TorchGNNFingerprint:
    ensure_torch_gnn_available()
    resolved_config = normalise_model_config(model_config)
    return TorchGNNFingerprint(
        species_list=[str(symbol).strip() for symbol in species_list],
        seed=int(seed),
        **resolved_config,
    )


def load_structures(structure_path: Path, requested_count: int = 0) -> list:
    structures = read(str(structure_path), index=":")
    return select_structures(list(structures), requested_count)


def read_single_structure(structure_path: Path, index: int = 0):
    return read(str(structure_path), index=int(index))


def write_structure(structure_path: Path, atoms) -> None:
    write(structure_path, atoms, format="extxyz")


def save_target_fingerprint(target_fingerprint: np.ndarray, output_path: Path) -> None:
    np.save(output_path, np.asarray(target_fingerprint, dtype=np.float32))


def load_target_fingerprint(fingerprint_path: Path) -> np.ndarray:
    suffix = fingerprint_path.suffix.lower()
    if suffix == ".npy":
        target_fingerprint = np.load(fingerprint_path)
    elif suffix == ".json":
        target_fingerprint = np.asarray(json.loads(fingerprint_path.read_text()), dtype=np.float32)
    else:
        raise ValueError(
            "target fingerprint files must use .npy or .json format; "
            f"received {fingerprint_path}"
        )
    return np.asarray(target_fingerprint, dtype=np.float32).reshape(-1)


def save_model_checkpoint(
    model: TorchGNNFingerprint,
    checkpoint_path: Path,
    species_list: Sequence[str],
    model_config: dict[str, Any],
    training_config: dict[str, Any],
    training_history: Sequence[float],
) -> dict[str, Any]:
    payload = {
        "checkpoint_version": 1,
        "species_list": [str(symbol).strip() for symbol in species_list],
        "model_config": {
            **normalise_model_config(model_config),
            "component_weight": list(
                normalise_component_weight(model_config.get("component_weight", DEFAULT_COMPONENT_WEIGHT))
            ),
        },
        "training_config": copy.deepcopy(training_config),
        "training_history": [float(value) for value in training_history],
        "is_fitted": bool(model.is_fitted),
        "state_dict": model.state_dict(),
    }
    torch.save(payload, checkpoint_path)
    return payload


def load_model_from_checkpoint(
    checkpoint_path: Path,
    map_location: str = "cpu",
) -> tuple[TorchGNNFingerprint, dict[str, Any]]:
    ensure_torch_gnn_available()
    checkpoint = torch.load(checkpoint_path, map_location=map_location)
    if not isinstance(checkpoint, dict):
        raise ValueError(f"Unexpected checkpoint payload in {checkpoint_path}")
    species_list = checkpoint.get("species_list")
    model_config = checkpoint.get("model_config")
    if not species_list or not isinstance(model_config, dict):
        raise ValueError(f"Checkpoint is missing model metadata: {checkpoint_path}")
    model = create_model(
        seed=int(model_config.get("seed", checkpoint.get("training_config", {}).get("seed", 42))),
        species_list=species_list,
        model_config=model_config,
    )
    state_dict = checkpoint.get("state_dict")
    if state_dict is None:
        raise ValueError(f"Checkpoint is missing state_dict: {checkpoint_path}")
    model.load_state_dict(state_dict)
    model._is_fitted = bool(checkpoint.get("is_fitted", True))
    return model, checkpoint


def build_inverse_design_options(
    target_atoms=None,
    fingerprint_loss_weight: float = 1.0,
    target_vertex_weight: float = 0.0,
    target_position_weight: float = 0.0,
    inverse_lr_decay_rate: float = 0.0,
    inverse_restarts: int = 1,
    inverse_restart_noise_scale: float = 0.0,
    repulsion_weight: float = 10.0,
    minimum_distance_scale: float = 0.75,
    cell_violation_weight: float = 0.0,
    coordinate_clip_value: float | None = None,
) -> dict[str, Any]:
    if float(target_vertex_weight) != 0.0:
        raise ValueError("target_vertex_weight must remain 0.0 for plan-compliant inverse design")
    if float(target_position_weight) != 0.0:
        raise ValueError("target_position_weight must remain 0.0 for plan-compliant inverse design")
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


def fingerprint_mse(model: TorchGNNFingerprint, atoms, target_fingerprint: np.ndarray) -> float:
    prediction = model.predict(atoms)
    return float(np.mean((prediction - target_fingerprint) ** 2))


def score_candidate(target_atoms, candidate_atoms) -> float:
    return float(symmetry_aware_rmsd(target_atoms, candidate_atoms))


def compute_inverse_design_metrics(
    model: TorchGNNFingerprint,
    initial_atoms,
    optimised_atoms,
    target_fingerprint: np.ndarray,
    target_atoms=None,
) -> dict[str, Any]:
    initial_prediction = model.predict(initial_atoms)
    final_prediction = model.predict(optimised_atoms)
    initial_mse = float(np.mean((initial_prediction - target_fingerprint) ** 2))
    final_mse = float(np.mean((final_prediction - target_fingerprint) ** 2))
    metrics: dict[str, Any] = {
        "initial_fingerprint_mse": initial_mse,
        "final_fingerprint_mse": final_mse,
        "fingerprint_mse_reduction_fraction": (
            0.0 if initial_mse == 0.0 else 1.0 - final_mse / initial_mse
        ),
        "initial_fingerprint_l2": float(np.linalg.norm(initial_prediction - target_fingerprint)),
        "final_fingerprint_l2": float(np.linalg.norm(final_prediction - target_fingerprint)),
        "structure_update_norm": float(
            np.linalg.norm(optimised_atoms.get_positions() - initial_atoms.get_positions())
        ),
    }
    if target_atoms is not None:
        initial_delta = symmetry_aware_displacements(target_atoms, initial_atoms)
        final_delta = symmetry_aware_displacements(target_atoms, optimised_atoms)
        metrics.update(
            {
                "position_error_metric": "symmetry_aware_rmsd",
                "initial_rmsd": score_candidate(target_atoms, initial_atoms),
                "final_rmsd": score_candidate(target_atoms, optimised_atoms),
                "initial_position_difference": initial_delta.tolist(),
                "final_position_difference": final_delta.tolist(),
            }
        )
        if metrics["initial_rmsd"] != 0.0:
            metrics["rmsd_reduction_fraction"] = (
                1.0 - metrics["final_rmsd"] / metrics["initial_rmsd"]
            )
    return metrics


def write_inverse_design_log(log_path: Path, metrics: dict[str, Any]) -> None:
    lines = [
        f"Initial fingerprint MSE: {metrics['initial_fingerprint_mse']:.8e}",
        f"Final fingerprint MSE:   {metrics['final_fingerprint_mse']:.8e}",
        (
            "Fingerprint MSE reduction: "
            f"{100.0 * metrics['fingerprint_mse_reduction_fraction']:.2f}%"
        ),
        f"Initial fingerprint L2: {metrics['initial_fingerprint_l2']:.8e}",
        f"Final fingerprint L2:   {metrics['final_fingerprint_l2']:.8e}",
        f"Structure update norm:  {metrics['structure_update_norm']:.8e}",
    ]
    if "initial_rmsd" in metrics and "final_rmsd" in metrics:
        lines.extend(
            [
                f"Initial RMSD:          {metrics['initial_rmsd']:.8e}",
                f"Final RMSD:            {metrics['final_rmsd']:.8e}",
                f"RMSD reduction:        {100.0 * metrics['rmsd_reduction_fraction']:.2f}%",
            ]
        )
    log_path.write_text("\n".join(lines) + "\n")


def print_position_differences(target_atoms, initial_atoms, optimised_atoms) -> None:
    initial_delta = symmetry_aware_displacements(target_atoms, initial_atoms)
    final_delta = symmetry_aware_displacements(target_atoms, optimised_atoms)
    optimisation_delta = minimum_image_displacements(initial_atoms, optimised_atoms)

    print("Initial position differences relative to the target structure (Angstrom):")
    print(np.array2string(initial_delta, precision=6, suppress_small=False))
    print()
    print("Final position differences relative to the target structure (Angstrom):")
    print(np.array2string(final_delta, precision=6, suppress_small=False))
    print()
    print("Net position update from initial to final structure (Angstrom):")
    print(np.array2string(optimisation_delta, precision=6, suppress_small=False))
