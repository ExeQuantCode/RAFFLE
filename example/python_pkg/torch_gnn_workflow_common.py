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
) -> dict[str, Any]:
    if target_atoms is None and (
        float(target_vertex_weight) > 0.0 or float(target_position_weight) > 0.0
    ):
        raise ValueError(
            "target_vertex_weight and target_position_weight require a target structure"
        )
    return {
        "target_atoms": target_atoms,
        "fingerprint_loss_weight": float(fingerprint_loss_weight),
        "target_vertex_weight": float(target_vertex_weight),
        "target_position_weight": float(target_position_weight),
        "inverse_lr_decay_rate": float(inverse_lr_decay_rate),
        "num_restarts": int(inverse_restarts),
        "restart_noise_scale": float(inverse_restart_noise_scale),
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
