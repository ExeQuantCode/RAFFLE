"""Replay a W&B carbon inverse-design run and export the trained surrogate checkpoint."""

from __future__ import annotations

import argparse
import ast
import copy
import json
import numpy as np
from pathlib import Path
import shutil
import sys
from typing import Any

import wandb


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from torch_gnn_carbon_wandb import (  # noqa: E402
    REPRODUCIBILITY_METADATA_KEY,
    RESOLVED_WORKFLOW_CONFIG_KEY,
    WANDB_PROJECT,
    parse_serialized_run_payload,
    validate_plan_constraints,
)
from torch_gnn_carbon_workflow_example import (  # noqa: E402
    build_inverse_design_options,
    minimum_training_carbon_count,
    select_carbon_structures,
    sweep_epochs,
)
from torch_gnn_workflow_common import (  # noqa: E402
    load_structures,
    load_model_from_checkpoint,
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
    parser.add_argument(
        "--replay-training",
        action="store_true",
        help=(
            "Replay training from the serialized resolved W&B run config instead of "
            "loading a persisted checkpoint."
        ),
    )
    return parser.parse_args(argv)


def _strip_internal_config_keys(config: dict[str, Any]) -> dict[str, Any]:
    return {
        key: value
        for key, value in config.items()
        if not str(key).startswith("_")
    }


def _restore_reference_structure(payload: dict[str, Any]):
    try:
        from ase import Atoms  # type: ignore
    except ImportError as exc:
        raise RuntimeError(
            "ASE is required to restore the serialized reference structure for replay"
        ) from exc
    return Atoms(
        symbols=list(payload["symbols"]),
        positions=payload["positions"],
        cell=payload["cell"],
        pbc=payload.get("pbc", True),
    )


def _restore_perturbed_structure(reference_structure, config: dict[str, Any]):
    perturbation = np.asarray(config["perturbation_matrix"], dtype=np.float32)
    reference_positions = np.asarray(reference_structure.get_positions(), dtype=np.float32)
    if perturbation.shape != reference_positions.shape:
        raise ValueError(
            "Serialized perturbation matrix shape does not match the serialized "
            f"reference structure: {perturbation.shape} != {reference_positions.shape}"
        )
    perturbed = reference_structure.copy()
    perturbed.set_positions(reference_positions + perturbation)
    fixed_atoms = np.zeros(len(perturbed), dtype=bool)
    fixed_atoms[: max(int(config["fixed_leading_atoms"]), 0)] = True
    return perturbed, fixed_atoms


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


def _parse_wandb_scalar(value: str) -> Any:
    text = value.strip()
    lowered = text.lower()
    if lowered == "true":
        return True
    if lowered == "false":
        return False
    if lowered in {"none", "null"}:
        return None
    try:
        return ast.literal_eval(text)
    except (SyntaxError, ValueError):
        return text


def _parse_wandb_files_config(config_path: Path) -> dict[str, Any]:
    resolved: dict[str, Any] = {}
    current_key: str | None = None

    for line in config_path.read_text().splitlines():
        if not line.strip():
            continue
        if not line.startswith(" "):
            current_key = line[:-1] if line.endswith(":") else None
            continue
        if current_key is None:
            continue
        if not line.startswith("    value:"):
            continue
        raw_value = line.split(":", 1)[1].strip()
        if raw_value:
            resolved[current_key] = _parse_wandb_scalar(raw_value)
        current_key = None

    if not resolved:
        raise ValueError(f"No persisted W&B config values found in {config_path}")
    return resolved


def _parse_json_file(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text())
    if not isinstance(payload, dict):
        raise ValueError(f"Expected JSON object in {path}")
    return payload


def _resolve_local_identity(debug_log_path: Path) -> dict[str, str | None]:
    identity: dict[str, str | None] = {"entity": None, "project": None}
    if not debug_log_path.exists():
        return identity

    marker = "finishing run "
    for line in reversed(debug_log_path.read_text().splitlines()):
        if marker not in line:
            continue
        run_path = line.split(marker, 1)[1].strip()
        parts = [part for part in run_path.split("/") if part]
        if len(parts) >= 3:
            identity["entity"] = parts[-3]
            identity["project"] = parts[-2]
        break
    return identity


def _resolve_local_run_name(run_dir: Path) -> str | None:
    output_log_path = run_dir / "files" / "output.log"
    if not output_log_path.exists():
        return None
    marker = '  "run_name": '
    for line in output_log_path.read_text().splitlines():
        if not line.startswith(marker):
            continue
        raw_value = line[len(marker) :].strip().rstrip(",")
        if raw_value in {"null", '""'}:
            return None
        try:
            resolved = json.loads(raw_value)
        except json.JSONDecodeError:
            return raw_value.strip('"') or None
        return None if resolved in {None, ""} else str(resolved)
    return None


def _resolve_local_output_dir(run_dir: Path) -> Path | None:
    summary_path = run_dir / "files" / "wandb-summary.json"
    if not summary_path.exists():
        return None
    output_dir = _parse_json_file(summary_path).get("output_dir")
    if output_dir in {None, ""}:
        return None
    resolved = Path(str(output_dir))
    if not resolved.is_absolute():
        resolved = (SCRIPT_DIR / resolved).resolve()
    return resolved


def _find_adjacent_checkpoint_assets(checkpoint_path: Path) -> dict[str, Path]:
    source_dir = checkpoint_path.parent
    assets: dict[str, Path] = {}
    for label, filename in (
        ("target_fingerprint", "torch_gnn_target_fingerprint.npy"),
        ("reference_structure", "torch_gnn_reference_structure.xyz"),
    ):
        candidate = source_dir / filename
        if candidate.exists():
            assets[label] = candidate
    return assets


def _resolve_wandb_run_path(run_record: dict[str, Any]) -> str | None:
    entity = run_record.get("entity")
    project = run_record.get("project")
    run_id = run_record.get("run_id")
    if entity and project and run_id:
        return f"{entity}/{project}/{run_id}"

    run_path = run_record.get("run_path")
    if not run_path:
        return None
    parts = [part for part in str(run_path).split("/") if part]
    if len(parts) == 3 and not Path(str(run_path)).is_absolute():
        return "/".join(parts)
    return None


def _serialise_checkpoint_source(checkpoint_source: dict[str, Any]) -> dict[str, Any]:
    serialised: dict[str, Any] = {}
    for key, value in checkpoint_source.items():
        if key == "assets":
            continue
        serialised[key] = str(value) if isinstance(value, Path) else value
    return serialised


def _resolve_remote_exact_checkpoint_source(
    run_record: dict[str, Any],
) -> dict[str, Any] | None:
    run_path = _resolve_wandb_run_path(run_record)
    if run_path is None:
        return None

    api = wandb.Api(overrides={"project": run_record["project"]})
    run = api.run(run_path)
    download_root = SCRIPT_DIR / ".wandb_artifact_cache" / str(run_record["run_id"])
    for artifact in run.logged_artifacts():
        if str(getattr(artifact, "type", "")) != "inverse-design-results":
            continue
        downloaded_dir = Path(artifact.download(root=str(download_root)))
        checkpoint_path = downloaded_dir / "torch_gnn_model_checkpoint.pt"
        if not checkpoint_path.exists():
            matches = list(downloaded_dir.rglob("torch_gnn_model_checkpoint.pt"))
            if not matches:
                continue
            checkpoint_path = matches[0]
        return {
            "mode": "wandb-artifact",
            "artifact_name": str(getattr(artifact, "name", "")) or None,
            "checkpoint_path": checkpoint_path,
            "assets": _find_adjacent_checkpoint_assets(checkpoint_path),
        }
    return None


def _resolve_local_run(run_id: str, project: str) -> dict[str, Any] | None:
    wandb_root = SCRIPT_DIR / "wandb"
    if not wandb_root.exists():
        return None
    matches = sorted(wandb_root.glob(f"run-*-{run_id}"))
    if not matches:
        return None
    run_dir = matches[-1]
    files_config_path = run_dir / "files" / "config.yaml"
    debug_log_path = run_dir / "logs" / "debug.log"
    local_identity = _resolve_local_identity(debug_log_path)
    if files_config_path.exists():
        config = _strip_internal_config_keys(
            _parse_wandb_files_config(files_config_path)
        )
    else:
        if not debug_log_path.exists():
            raise ValueError(
                f"Local W&B cache for run '{run_id}' is missing config.yaml and "
                f"{debug_log_path.name}"
            )
        config = _strip_internal_config_keys(_parse_debug_log_config(debug_log_path))
    if not debug_log_path.exists() and not files_config_path.exists():
        raise ValueError(
            f"Local W&B cache for run '{run_id}' is missing {debug_log_path.name}"
        )
    return {
        "source": "local-cache",
        "entity": local_identity["entity"],
        "project": local_identity["project"] or project,
        "run_id": run_id,
        "run_name": _resolve_local_run_name(run_dir),
        "run_path": str(run_dir),
        "url": None,
        "config": config,
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
    if entity is None and "/" not in run_id:
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
    if RESOLVED_WORKFLOW_CONFIG_KEY not in run_config:
        raise ValueError(
            "Run config does not contain the serialized resolved workflow payload "
            f"'{RESOLVED_WORKFLOW_CONFIG_KEY}'. This run predates config-only replay."
        )
    resolved = copy.deepcopy(
        parse_serialized_run_payload(
            run_config[RESOLVED_WORKFLOW_CONFIG_KEY],
            key=RESOLVED_WORKFLOW_CONFIG_KEY,
        )
    )
    if not resolved.get("spec_version"):
        raise ValueError(
            "Serialized resolved workflow payload is missing 'spec_version' and "
            "cannot be trusted for deterministic replay."
        )
    if epochs_override is not None:
        resolved["epochs"] = int(epochs_override)
    resolved["inverse_step_values"] = [
        int(value) for value in resolved.get("inverse_step_values", [])
    ]
    resolved["step_size_values"] = [
        float(value) for value in resolved.get("step_size_values", [])
    ]
    resolved["replay_category_weights"] = {
        str(key): float(value)
        for key, value in dict(resolved["replay_category_weights"]).items()
    }
    resolved["model_config"] = dict(resolved["model_config"])
    if "component_weight" in resolved["model_config"]:
        resolved["model_config"]["component_weight"] = tuple(
            float(value)
            for value in resolved["model_config"]["component_weight"]
        )
    if resolved.get("coordinate_clip_value") is not None:
        resolved["coordinate_clip_value"] = float(resolved["coordinate_clip_value"])
    metadata_value = run_config.get(REPRODUCIBILITY_METADATA_KEY)
    resolved["source_reproducibility_metadata"] = (
        None
        if metadata_value is None
        else parse_serialized_run_payload(
            metadata_value,
            key=REPRODUCIBILITY_METADATA_KEY,
        )
    )
    validate_plan_constraints(resolved)
    return resolved


def replay_training(
    config: dict[str, Any],
) -> dict[str, Any]:
    repo_root = Path(__file__).resolve().parents[2]
    carbon_dataset_path = Path(str(config["carbon_dataset_path"]))
    carbon_xyz = (
        carbon_dataset_path
        if carbon_dataset_path.is_absolute()
        else repo_root / carbon_dataset_path
    )
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

    reference_structure = _restore_reference_structure(
        dict(config["reference_structure"])
    )
    augmented_structures: list[Any] = []
    perturbed_structure, fixed_atoms = _restore_perturbed_structure(
        reference_structure,
        config,
    )
    inverse_design_options = build_inverse_design_options(
        fingerprint_loss_weight=float(config["fingerprint_loss_weight"]),
        target_vertex_weight=float(config["target_vertex_weight"]),
        target_position_weight=float(config["target_position_weight"]),
        inverse_lr_decay_rate=float(config["inverse_lr_decay_rate"]),
        repulsion_weight=float(config["repulsion_weight"]),
        minimum_distance_scale=float(config["minimum_distance_scale"]),
        cell_violation_weight=float(config["cell_violation_weight"]),
        coordinate_clip_value=config["coordinate_clip_value"],
    )

    training_history: list[float] = []

    def training_observer(epoch: int, loss: float) -> None:
        del epoch
        training_history.append(float(loss))

    model, _, _, _, _, rollout_metrics = sweep_epochs(
        carbon_structures=carbon_structures,
        augmented_structures=augmented_structures,
        original=reference_structure,
        perturbed=perturbed_structure,
        fixed_atoms=fixed_atoms,
        num_epochs=int(config["epochs"]),
        batch_size=int(config["batch_size"]),
        inverse_steps=int(config["inverse_steps"]),
        inverse_step_size=float(config["inverse_step_size"]),
        inverse_design_options=inverse_design_options,
        seed=int(config["seed"]),
        model_config=config["model_config"],
        training_observer=training_observer,
        rollout_stages=int(config["rollout_stages"]),
        rollout_epochs_per_stage=int(config["rollout_epochs_per_stage"]),
        rollout_step_stride=int(config["rollout_step_stride"]),
        replay_buffer_capacity=int(config["replay_buffer_capacity"]),
        replay_sample_size=int(config["replay_sample_size"]),
        replay_category_weights=dict(config["replay_category_weights"]),
        rollout_drift_threshold=float(config["rollout_drift_threshold"]),
        rollout_high_error_threshold=float(config["rollout_high_error_threshold"]),
        rollout_instability_threshold=float(config["rollout_instability_threshold"]),
    )
    target_fingerprint = model.compute_reference_fingerprint(reference_structure)
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


def resolve_exact_checkpoint_source(run_record: dict[str, Any]) -> dict[str, Any] | None:
    run_path = run_record.get("run_path")
    if run_path:
        run_dir = Path(str(run_path))
        candidates: list[Path] = []
        output_dir = _resolve_local_output_dir(run_dir)
        if output_dir is not None:
            candidates.append(output_dir / "torch_gnn_model_checkpoint.pt")
        candidates.append(run_dir / "files" / "torch_gnn_model_checkpoint.pt")

        for candidate in candidates:
            if not candidate.exists():
                continue
            return {
                "mode": "local-checkpoint",
                "checkpoint_path": candidate,
                "assets": _find_adjacent_checkpoint_assets(candidate),
            }
    return _resolve_remote_exact_checkpoint_source(run_record)


def export_exact_checkpoint(
    *,
    run_record: dict[str, Any],
    checkpoint_source: dict[str, Any],
    output_dir: Path,
) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)

    source_checkpoint_path = Path(str(checkpoint_source["checkpoint_path"]))
    config_path = output_dir / "torch_gnn_wandb_replay_config.json"
    checkpoint_path = output_dir / "torch_gnn_model_checkpoint.pt"
    fingerprint_path = output_dir / "torch_gnn_target_fingerprint.npy"
    reference_structure_path = output_dir / "torch_gnn_reference_structure.xyz"
    metrics_path = output_dir / "torch_gnn_wandb_replay_metrics.json"

    if source_checkpoint_path.resolve() != checkpoint_path.resolve():
        shutil.copy2(source_checkpoint_path, checkpoint_path)

    copied_output_files: dict[str, str | None] = {
        "config": str(config_path),
        "checkpoint": str(checkpoint_path),
        "target_fingerprint": None,
        "reference_structure": None,
        "metrics": str(metrics_path),
    }
    for label, source_path in checkpoint_source.get("assets", {}).items():
        destination = (
            fingerprint_path if label == "target_fingerprint" else reference_structure_path
        )
        if source_path.resolve() != destination.resolve():
            shutil.copy2(source_path, destination)
        copied_output_files[label] = str(destination)

    _, checkpoint_payload = load_model_from_checkpoint(checkpoint_path)
    training_history = [
        float(value) for value in checkpoint_payload.get("training_history", [])
    ]
    replay_metadata = {
        "source": run_record["source"],
        "entity": run_record.get("entity"),
        "project": run_record["project"],
        "run_id": run_record["run_id"],
        "run_name": run_record.get("run_name"),
        "run_path": run_record.get("run_path"),
        "url": run_record.get("url"),
    }
    exact_recovery_config = {
        "source_wandb_run": replay_metadata,
        "recovery_mode": "exact-checkpoint",
        "checkpoint_source": _serialise_checkpoint_source(checkpoint_source),
        "checkpoint_training_config": checkpoint_payload.get("training_config"),
        "model_config": checkpoint_payload.get("model_config"),
        "species_list": checkpoint_payload.get("species_list"),
        "output_files": copied_output_files,
    }
    write_json(config_path, exact_recovery_config)

    metrics = {
        "source_wandb_run": replay_metadata,
        "recovery_mode": "exact-checkpoint",
        "checkpoint_source": _serialise_checkpoint_source(checkpoint_source),
        "training_history": training_history,
        "initial_training_loss": (
            None if not training_history else float(training_history[0])
        ),
        "final_training_loss": (
            None if not training_history else float(training_history[-1])
        ),
        "species_list": checkpoint_payload.get("species_list"),
        "model_config": checkpoint_payload.get("model_config"),
        "training_config": checkpoint_payload.get("training_config"),
        "output_files": copied_output_files,
    }
    write_json(metrics_path, metrics)
    return metrics


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
    workflow_config = {
        key: value for key, value in config.items() if key != "reference_structure_weight"
    }
    effective_config = {
        "repo_root": str(replay_result["repo_root"]),
        "source_wandb_run": replay_metadata,
        "recovery_mode": "config-replay",
        "workflow_config": {
            **workflow_config,
            "coordinate_clip_value": workflow_config["coordinate_clip_value"],
            "replay_category_weights": dict(workflow_config["replay_category_weights"]),
        },
        "source_reproducibility_metadata": config.get("source_reproducibility_metadata"),
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
        "recovery_mode": "config-replay",
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
        "source_reproducibility_metadata": config.get("source_reproducibility_metadata"),
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
    output_dir = args.output_dir.resolve()
    if args.epochs is not None and not args.replay_training:
        raise ValueError(
            "--epochs only applies to approximate retraining. Pass --replay-training "
            "to opt into replaying training."
        )

    if args.replay_training:
        config = build_replay_config(
            run_record["config"],
            epochs_override=args.epochs,
        )
        replay_result = replay_training(config)
        metrics = export_checkpoint(
            run_record=run_record,
            config=config,
            replay_result=replay_result,
            output_dir=output_dir,
        )
    else:
        checkpoint_source = resolve_exact_checkpoint_source(run_record)
        if checkpoint_source is None:
            raise ValueError(
                "Exact checkpoint recovery is unavailable for this run because no "
                "torch_gnn_model_checkpoint.pt was persisted in the run outputs. "
                "Pass --replay-training for approximate retraining instead."
            )
        metrics = export_exact_checkpoint(
            run_record=run_record,
            checkpoint_source=checkpoint_source,
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
