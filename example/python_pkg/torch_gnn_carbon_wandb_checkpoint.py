"""Replay a W&B carbon inverse-design run and export the trained surrogate checkpoint."""

from __future__ import annotations

import argparse
import ast
import json
from pathlib import Path
import sys
from typing import Any

import wandb


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from torch_gnn_carbon_wandb import (  # noqa: E402
    DEFAULT_ARCHITECTURE,
    WANDB_PROJECT,
    build_model_config,
    validate_plan_constraints,
)
from torch_gnn_carbon_workflow_example import (  # noqa: E402
    CELL_VIOLATION_WEIGHT,
    COORDINATE_CLIP_VALUE,
    DEFAULT_MODEL_CONFIG,
    DEFAULT_REPLAY_CATEGORY_WEIGHTS,
    FINGERPRINT_LOSS_WEIGHT,
    FIXED_LEADING_ATOMS,
    INVERSE_LR_DECAY_RATE,
    INVERSE_RESTART_NOISE_SCALE,
    INVERSE_RESTARTS,
    MINIMUM_DISTANCE_SCALE,
    REPLAY_BUFFER_CAPACITY,
    REPLAY_SAMPLE_SIZE,
    REPULSION_WEIGHT,
    ROLLOUT_DRIFT_THRESHOLD,
    ROLLOUT_EPOCHS_PER_STAGE,
    ROLLOUT_HIGH_ERROR_THRESHOLD,
    ROLLOUT_INSTABILITY_THRESHOLD,
    ROLLOUT_STAGES,
    ROLLOUT_STEP_STRIDE,
    TARGET_POSITION_WEIGHT,
    TARGET_VERTEX_WEIGHT,
    build_inverse_design_options,
    build_perturbed_structure,
    minimum_training_carbon_count,
    parse_category_weights,
    run_rollout_retraining,
    select_carbon_structures,
    sweep_epochs,
)
from torch_gnn_workflow_common import (  # noqa: E402
    default_reference_structure,
    load_structures,
    save_model_checkpoint,
    save_target_fingerprint,
    write_json,
    write_structure,
)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "run_id",
        type=str,
        help="W&B run id, project/run_id, or entity/project/run_id.",
    )
    parser.add_argument(
        "--entity",
        type=str,
        default=None,
        help="W&B entity when run_id is not a full entity/project/run_id path.",
    )
    parser.add_argument("--project", type=str, default=WANDB_PROJECT)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path.cwd(),
        help=(
            "Directory for the exported checkpoint and replay artifacts. "
            "Defaults to the current working directory."
        ),
    )
    parser.add_argument(
        "--epochs",
        type=int,
        default=None,
        help="Override the W&B-configured number of training epochs.",
    )
    return parser.parse_args(argv)


def _strip_internal_config_keys(config: dict[str, Any]) -> dict[str, Any]:
    return {
        key: value
        for key, value in config.items()
        if not str(key).startswith("_")
    }


def _resolve_optional_float(value: Any, default: float | None) -> float | None:
    if value is None:
        return default
    if isinstance(value, str) and value.strip().lower() in {"", "none", "null"}:
        return None
    return float(value)


def _resolve_replay_category_weights(value: Any) -> dict[str, float]:
    if value is None:
        return dict(DEFAULT_REPLAY_CATEGORY_WEIGHTS)
    if isinstance(value, dict):
        return {str(key): float(raw_value) for key, raw_value in value.items()}
    return parse_category_weights(str(value))


def _parse_debug_log_config(debug_log_path: Path) -> dict[str, Any]:
    for line in reversed(debug_log_path.read_text().splitlines()):
        prefix = "config: "
        if not line.startswith(prefix):
            continue
        payload = ast.literal_eval(line[len(prefix) :])
        if not isinstance(payload, dict):
            raise ValueError(f"Unexpected config payload in {debug_log_path}")
        return payload
    raise ValueError(f"No config payload found in {debug_log_path}")


def _resolve_local_run(run_id: str, project: str) -> dict[str, Any] | None:
    wandb_root = SCRIPT_DIR / "wandb"
    if not wandb_root.exists():
        return None
    matches = sorted(wandb_root.glob(f"run-*-{run_id}"))
    if not matches:
        return None
    debug_log_path = matches[-1] / "logs" / "debug.log"
    if not debug_log_path.exists():
        raise ValueError(
            f"Local W&B cache for run '{run_id}' is missing {debug_log_path.name}"
        )
    return {
        "source": "local-cache",
        "entity": None,
        "project": project,
        "run_id": run_id,
        "run_name": None,
        "run_path": str(matches[-1]),
        "url": None,
        "config": _strip_internal_config_keys(_parse_debug_log_config(debug_log_path)),
    }


def _resolve_remote_run_path(run_id: str, project: str, entity: str | None) -> str:
    parts = [part for part in run_id.split("/") if part]
    if len(parts) == 3:
        return "/".join(parts)
    if entity is None:
        raise ValueError(
            "Provide --entity or a full entity/project/run_id path when the run id is "
            "not available in the local W&B cache."
        )
    if len(parts) == 2:
        return f"{entity}/{parts[0]}/{parts[1]}"
    if len(parts) == 1:
        return f"{entity}/{project}/{parts[0]}"
    raise ValueError(f"Unsupported run identifier '{run_id}'")


def resolve_run_record(
    run_id: str,
    *,
    project: str,
    entity: str | None = None,
) -> dict[str, Any]:
    local_record = None
    if "/" not in run_id:
        local_record = _resolve_local_run(run_id, project)
    if local_record is not None:
        return local_record

    run_path = _resolve_remote_run_path(run_id, project, entity)
    api = wandb.Api(overrides={"project": project})
    run = api.run(run_path)
    return {
        "source": "wandb-api",
        "entity": str(getattr(run, "entity", "")) or entity,
        "project": str(getattr(run, "project", project)),
        "run_id": str(getattr(run, "id", run_id.rsplit("/", 1)[-1])),
        "run_name": getattr(run, "name", None),
        "run_path": run_path,
        "url": getattr(run, "url", None),
        "config": _strip_internal_config_keys(dict(getattr(run, "config", {}))),
    }


def build_replay_config(
    run_config: dict[str, Any],
    *,
    epochs_override: int | None = None,
) -> dict[str, Any]:
    model_source = {
        "architecture": str(run_config.get("architecture", DEFAULT_ARCHITECTURE)),
        "hidden_dim": int(run_config.get("hidden_dim", DEFAULT_MODEL_CONFIG["hidden_dim"])),
        "num_message_layers": int(
            run_config.get(
                "num_message_layers",
                DEFAULT_MODEL_CONFIG["num_message_layers"],
            )
        ),
        "learning_rate": float(
            run_config.get("learning_rate", DEFAULT_MODEL_CONFIG["learning_rate"])
        ),
        "model_lr_decay_rate": float(
            run_config.get("model_lr_decay_rate", DEFAULT_MODEL_CONFIG["lr_decay_rate"])
        ),
        "smooth_cutoff_width": float(
            run_config.get(
                "smooth_cutoff_width",
                DEFAULT_MODEL_CONFIG["smooth_cutoff_width"],
            )
        ),
        "reference_layer_type": int(
            run_config.get(
                "reference_layer_type",
                DEFAULT_MODEL_CONFIG["reference_layer_type"],
            )
        ),
        "component_weight_2body": float(
            run_config.get(
                "component_weight_2body",
                DEFAULT_MODEL_CONFIG["component_weight"][0],
            )
        ),
        "component_weight_3body": float(
            run_config.get(
                "component_weight_3body",
                DEFAULT_MODEL_CONFIG["component_weight"][1],
            )
        ),
        "component_weight_4body": float(
            run_config.get(
                "component_weight_4body",
                DEFAULT_MODEL_CONFIG["component_weight"][2],
            )
        ),
    }
    resolved = {
        "architecture": model_source["architecture"],
        "carbon_count": int(run_config.get("carbon_count", -1)),
        "augmented_count": int(run_config.get("augmented_count", 0)),
        "epochs": int(
            run_config.get("epochs", 50)
            if epochs_override is None
            else epochs_override
        ),
        "batch_size": int(run_config.get("batch_size", 8)),
        "inverse_steps": int(run_config.get("inverse_steps", 400)),
        "inverse_step_size": float(run_config.get("inverse_step_size", 5.0e-3)),
        "fixed_leading_atoms": int(
            run_config.get("fixed_leading_atoms", FIXED_LEADING_ATOMS)
        ),
        "fingerprint_loss_weight": float(
            run_config.get("fingerprint_loss_weight", FINGERPRINT_LOSS_WEIGHT)
        ),
        "target_vertex_weight": float(
            run_config.get("target_vertex_weight", TARGET_VERTEX_WEIGHT)
        ),
        "target_position_weight": float(
            run_config.get("target_position_weight", TARGET_POSITION_WEIGHT)
        ),
        "inverse_lr_decay_rate": float(
            run_config.get("inverse_lr_decay_rate", INVERSE_LR_DECAY_RATE)
        ),
        "inverse_restarts": int(run_config.get("inverse_restarts", INVERSE_RESTARTS)),
        "inverse_restart_noise_scale": float(
            run_config.get(
                "inverse_restart_noise_scale",
                INVERSE_RESTART_NOISE_SCALE,
            )
        ),
        "repulsion_weight": float(
            run_config.get("repulsion_weight", REPULSION_WEIGHT)
        ),
        "minimum_distance_scale": float(
            run_config.get("minimum_distance_scale", MINIMUM_DISTANCE_SCALE)
        ),
        "cell_violation_weight": float(
            run_config.get("cell_violation_weight", CELL_VIOLATION_WEIGHT)
        ),
        "coordinate_clip_value": _resolve_optional_float(
            run_config.get("coordinate_clip_value"),
            COORDINATE_CLIP_VALUE,
        ),
        "rollout_stages": int(run_config.get("rollout_stages", ROLLOUT_STAGES)),
        "rollout_epochs_per_stage": int(
            run_config.get(
                "rollout_epochs_per_stage",
                ROLLOUT_EPOCHS_PER_STAGE,
            )
        ),
        "rollout_step_stride": int(
            run_config.get("rollout_step_stride", ROLLOUT_STEP_STRIDE)
        ),
        "replay_buffer_capacity": int(
            run_config.get("replay_buffer_capacity", REPLAY_BUFFER_CAPACITY)
        ),
        "replay_sample_size": int(
            run_config.get("replay_sample_size", REPLAY_SAMPLE_SIZE)
        ),
        "replay_category_weights": _resolve_replay_category_weights(
            run_config.get("replay_category_weights")
        ),
        "rollout_drift_threshold": float(
            run_config.get("rollout_drift_threshold", ROLLOUT_DRIFT_THRESHOLD)
        ),
        "rollout_high_error_threshold": float(
            run_config.get(
                "rollout_high_error_threshold",
                ROLLOUT_HIGH_ERROR_THRESHOLD,
            )
        ),
        "rollout_instability_threshold": float(
            run_config.get(
                "rollout_instability_threshold",
                ROLLOUT_INSTABILITY_THRESHOLD,
            )
        ),
        "seed": int(run_config.get("seed", 42)),
        "model_config": build_model_config(model_source),
    }
    validate_plan_constraints(resolved)
    return resolved


def replay_training(
    config: dict[str, Any],
) -> dict[str, Any]:
    repo_root = Path(__file__).resolve().parents[2]
    carbon_xyz = repo_root / "example" / "data" / "carbon.xyz"
    all_carbon_structures = load_structures(carbon_xyz)
    carbon_structures = select_carbon_structures(
        all_carbon_structures,
        int(config["carbon_count"]),
    )
    minimum_carbon_count = minimum_training_carbon_count(len(all_carbon_structures))
    if len(carbon_structures) < minimum_carbon_count:
        raise ValueError(
            "carbon_count must select at least 30% of example/data/carbon.xyz "
            f"structures ({minimum_carbon_count} minimum, received {len(carbon_structures)})"
        )

    reference_structure = default_reference_structure()
    augmented_structures: list[Any] = []
    perturbed_structure, fixed_atoms = build_perturbed_structure(
        reference_structure,
        fixed_leading_atoms=int(config["fixed_leading_atoms"]),
    )
    inverse_design_options = build_inverse_design_options(
        original=reference_structure,
        fingerprint_loss_weight=float(config["fingerprint_loss_weight"]),
        target_vertex_weight=float(config["target_vertex_weight"]),
        target_position_weight=float(config["target_position_weight"]),
        inverse_lr_decay_rate=float(config["inverse_lr_decay_rate"]),
        inverse_restarts=int(config["inverse_restarts"]),
        inverse_restart_noise_scale=float(config["inverse_restart_noise_scale"]),
        repulsion_weight=float(config["repulsion_weight"]),
        minimum_distance_scale=float(config["minimum_distance_scale"]),
        cell_violation_weight=float(config["cell_violation_weight"]),
        coordinate_clip_value=config["coordinate_clip_value"],
    )

    training_history: list[float] = []

    def training_observer(epoch: int, loss: float) -> None:
        del epoch
        training_history.append(float(loss))

    model, _, _, _, _ = sweep_epochs(
        carbon_structures=carbon_structures,
        augmented_structures=augmented_structures,
        original=reference_structure,
        perturbed=perturbed_structure,
        fixed_atoms=fixed_atoms,
        epoch_values=[int(config["epochs"])],
        batch_size=int(config["batch_size"]),
        inverse_steps=int(config["inverse_steps"]),
        inverse_step_size=float(config["inverse_step_size"]),
        inverse_design_options=inverse_design_options,
        seed=int(config["seed"]),
        model_config=config["model_config"],
        training_observer=training_observer,
    )
    target_fingerprint = model.compute_reference_fingerprint(reference_structure)
    model, rollout_metrics = run_rollout_retraining(
        model=model,
        carbon_structures=carbon_structures,
        original=reference_structure,
        perturbed=perturbed_structure,
        fixed_atoms=fixed_atoms,
        target_fingerprint=target_fingerprint,
        batch_size=int(config["batch_size"]),
        inverse_steps=int(config["inverse_steps"]),
        inverse_step_size=float(config["inverse_step_size"]),
        inverse_design_options=inverse_design_options,
        seed=int(config["seed"]),
        rollout_stages=int(config["rollout_stages"]),
        rollout_epochs_per_stage=int(config["rollout_epochs_per_stage"]),
        rollout_step_stride=int(config["rollout_step_stride"]),
        replay_buffer_capacity=int(config["replay_buffer_capacity"]),
        replay_sample_size=int(config["replay_sample_size"]),
        replay_category_weights=dict(config["replay_category_weights"]),
        rollout_drift_threshold=float(config["rollout_drift_threshold"]),
        rollout_high_error_threshold=float(config["rollout_high_error_threshold"]),
        rollout_instability_threshold=float(config["rollout_instability_threshold"]),
        training_observer=training_observer,
        training_epoch_offset=len(training_history),
    )
    return {
        "repo_root": repo_root,
        "species_list": ["C"],
        "model": model,
        "reference_structure": reference_structure,
        "target_fingerprint": target_fingerprint,
        "training_history": training_history,
        "rollout": rollout_metrics,
        "num_available_carbon_structures": len(all_carbon_structures),
        "num_training_carbon_structures": len(carbon_structures),
    }


def export_checkpoint(
    *,
    run_record: dict[str, Any],
    config: dict[str, Any],
    replay_result: dict[str, Any],
    output_dir: Path,
) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)

    config_path = output_dir / "torch_gnn_wandb_replay_config.json"
    checkpoint_path = output_dir / "torch_gnn_model_checkpoint.pt"
    fingerprint_path = output_dir / "torch_gnn_target_fingerprint.npy"
    reference_structure_path = output_dir / "torch_gnn_reference_structure.xyz"
    metrics_path = output_dir / "torch_gnn_wandb_replay_metrics.json"

    replay_metadata = {
        "source": run_record["source"],
        "entity": run_record.get("entity"),
        "project": run_record["project"],
        "run_id": run_record["run_id"],
        "run_name": run_record.get("run_name"),
        "run_path": run_record.get("run_path"),
        "url": run_record.get("url"),
    }
    effective_config = {
        "repo_root": str(replay_result["repo_root"]),
        "source_wandb_run": replay_metadata,
        "workflow_config": {
            **config,
            "coordinate_clip_value": config["coordinate_clip_value"],
            "replay_category_weights": dict(config["replay_category_weights"]),
        },
        "species_list": list(replay_result["species_list"]),
        "rollout": replay_result["rollout"],
        "output_files": {
            "config": str(config_path),
            "checkpoint": str(checkpoint_path),
            "target_fingerprint": str(fingerprint_path),
            "reference_structure": str(reference_structure_path),
            "metrics": str(metrics_path),
        },
    }

    write_json(config_path, effective_config)
    save_model_checkpoint(
        model=replay_result["model"],
        checkpoint_path=checkpoint_path,
        species_list=replay_result["species_list"],
        model_config=config["model_config"],
        training_config=effective_config,
        training_history=replay_result["training_history"],
    )
    save_target_fingerprint(replay_result["target_fingerprint"], fingerprint_path)
    write_structure(reference_structure_path, replay_result["reference_structure"])

    metrics = {
        "source_wandb_run": replay_metadata,
        "num_available_carbon_structures": replay_result["num_available_carbon_structures"],
        "num_training_carbon_structures": replay_result["num_training_carbon_structures"],
        "training_history": [float(value) for value in replay_result["training_history"]],
        "initial_training_loss": (
            None
            if not replay_result["training_history"]
            else float(replay_result["training_history"][0])
        ),
        "final_training_loss": (
            None
            if not replay_result["training_history"]
            else float(replay_result["training_history"][-1])
        ),
        "species_list": list(replay_result["species_list"]),
        "model_config": config["model_config"],
        "training_config": {
            key: value
            for key, value in effective_config.items()
            if key not in {"output_files"}
        },
        "output_files": effective_config["output_files"],
    }
    write_json(metrics_path, metrics)
    return metrics


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    run_record = resolve_run_record(
        args.run_id,
        project=str(args.project),
        entity=args.entity,
    )
    config = build_replay_config(
        run_record["config"],
        epochs_override=args.epochs,
    )
    replay_result = replay_training(config)
    output_dir = args.output_dir.resolve()
    metrics = export_checkpoint(
        run_record=run_record,
        config=config,
        replay_result=replay_result,
        output_dir=output_dir,
    )
    print(
        json.dumps(
            {
                "run_id": run_record["run_id"],
                "project": run_record["project"],
                "output_dir": str(output_dir),
                "checkpoint": metrics["output_files"]["checkpoint"],
                "metrics": metrics["output_files"]["metrics"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
