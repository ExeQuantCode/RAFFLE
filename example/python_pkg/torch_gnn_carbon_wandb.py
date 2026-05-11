"""W&B online runner and sweep entrypoint for torch carbon inverse design."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import wandb

from torch_gnn_carbon_workflow_example import (
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
    TARGET_VERTEX_WEIGHT,
    default_epoch_values,
    default_inverse_step_values,
    default_step_sizes,
    parse_category_weights,
    parse_float_list,
    parse_int_list,
    run_example,
)


WANDB_PROJECT = "raffle-inverse-design-new"
DEFAULT_ARCHITECTURE = "torch_gnn_residual"
DEFAULT_TAGS = ["carbon", "diamond", "inverse-design"]


def parse_tags(value: str) -> list[str]:
    return [item.strip() for item in value.split(",") if item.strip()]


def build_model_config(config: dict) -> dict:
    return {
        "architecture": str(config["architecture"]),
        "hidden_dim": int(config["hidden_dim"]),
        "num_message_layers": int(config["num_message_layers"]),
        "learning_rate": float(config["learning_rate"]),
        "lr_decay_rate": float(config["model_lr_decay_rate"]),
        "smooth_cutoff_width": float(config["smooth_cutoff_width"]),
        "reference_layer_type": int(config["reference_layer_type"]),
        "component_weight": (
            float(config["component_weight_2body"]),
            float(config["component_weight_3body"]),
            float(config["component_weight_4body"]),
        ),
    }


def validate_plan_constraints(config: dict) -> None:
    if int(config["augmented_count"]) != 0:
        raise ValueError("augmented_count must remain 0 for plan_model benchmarks")
    if float(config["target_vertex_weight"]) != 0.0:
        raise ValueError("target_vertex_weight must remain 0.0 for plan_model benchmarks")
    if float(config["target_position_weight"]) != 0.0:
        raise ValueError("target_position_weight must remain 0.0 for plan_model benchmarks")


def build_sweep_config(args: argparse.Namespace) -> dict:
    parameters = {
        "carbon_count": {"values": [-1]},
        "hidden_dim": {"values": [64, 80, 96, 128]},
        "num_message_layers": {"values": [2, 3, 4]},
        "learning_rate": {"values": [1.0e-3, 5.0e-4, 2.5e-4]},
        "model_lr_decay_rate": {"values": [1.0e-2, 5.0e-3, 1.0e-3]},
        "smooth_cutoff_width": {"values": [0.1, 0.15, 0.2, 0.3]},
        "component_weight_2body": {"values": [1.0, 2.0, 4.0, 6.0]},
        "component_weight_3body": {"values": [0.0, 1.0, 2.0, 3.0]},
        "component_weight_4body": {"values": [0.0, 1.0, 2.0, 3.0]},
        "batch_size": {"values": [4, 8, 16]},
        "fingerprint_loss_weight": {"values": [1.0]},
        "target_vertex_weight": {"values": [0.0]},
        "target_position_weight": {"values": [0.0]},
        "inverse_steps": {"values": [100, 200, 300]},
        "inverse_step_size": {"values": [1.0e-3, 5.0e-3, 1.0e-2, 1.e-1]},
        "augmented_count": {"values": [0]},
        "inverse_restarts": {"values": [0]},
        "inverse_restart_noise_scale": {"values": [0.0, 0.005, 0.01]},
        "repulsion_weight": {"values": [5.0, 10.0, 20.0]},
        "minimum_distance_scale": {"values": [0.7, 0.75, 0.8]},
        "cell_violation_weight": {"values": [0.0, 0.01]},
        "coordinate_clip_value": {"values": [0.25, 0.5, 1.0]},
        "rollout_stages": {"values": [1, 2]},
        "rollout_epochs_per_stage": {"values": [1, 2]},
        "rollout_step_stride": {"values": [10, 25, 50]},
        "replay_buffer_capacity": {"values": [32, 64, 128]},
        "replay_sample_size": {"values": [4, 8, 16]},
        "seed": {"values": [11, 42, 101]},
    }
    if args.sweep_profile in {"vertex-focus", "fingerprint-focus"}:
        parameters = {
            "carbon_count": {"values": [-1]},
            "hidden_dim": {"values": [128, 192, 256]},
            "num_message_layers": {"values": [4, 6]},
            "learning_rate": {"values": [5.0e-4, 2.5e-4, 1.0e-4]},
            "model_lr_decay_rate": {"values": [1.0e-3, 5.0e-4]},
            "smooth_cutoff_width": {"values": [0.15, 0.2]},
            "component_weight_2body": {"values": [2.0, 4.0, 6.0, 8.0]},
            "component_weight_3body": {"values": [1.0, 2.0]},
            "component_weight_4body": {"values": [1.0, 2.0]},
            "batch_size": {"values": [4, 8]},
            "fingerprint_loss_weight": {"values": [0.5, 1.0, 2.0]},
            "target_vertex_weight": {"values": [0.0]},
            "target_position_weight": {"values": [0.0]},
            "inverse_steps": {"values": [400, 800, 1200]},
            "inverse_step_size": {"values": [5.0e-4, 1.0e-3, 2.5e-3]},
            "augmented_count": {"values": [0]},
            "inverse_restarts": {"values": [1, 2, 4]},
            "inverse_restart_noise_scale": {"values": [0.0, 0.005, 0.01]},
            "repulsion_weight": {"values": [5.0, 10.0, 20.0]},
            "minimum_distance_scale": {"values": [0.7, 0.75, 0.8]},
            "cell_violation_weight": {"values": [0.0, 0.01]},
            "coordinate_clip_value": {"values": [0.25, 0.5]},
            "rollout_stages": {"values": [1, 2]},
            "rollout_epochs_per_stage": {"values": [1, 2]},
            "rollout_step_stride": {"values": [10, 25]},
            "replay_buffer_capacity": {"values": [64, 128]},
            "replay_sample_size": {"values": [8, 16]},
            "seed": {"values": [11, 42, 101]},
        }
    if args.sweep_profile == "inverse-basin":
        parameters = {
            "carbon_count": {"values": [-1]},
            "hidden_dim": {"values": [96, 128, 160]},
            "num_message_layers": {"values": [3, 4, 5]},
            "learning_rate": {"values": [5.0e-4, 2.5e-4]},
            "model_lr_decay_rate": {"values": [5.0e-3, 1.0e-3]},
            "smooth_cutoff_width": {"values": [0.15, 0.2]},
            "component_weight_2body": {"values": [4.0, 5.0, 6.0]},
            "component_weight_3body": {"values": [1.0]},
            "component_weight_4body": {"values": [1.0]},
            "batch_size": {"values": [8]},
            "fingerprint_loss_weight": {"values": [0.25, 0.275, 0.3]},
            "target_vertex_weight": {"values": [0.0]},
            "target_position_weight": {"values": [0.0]},
            "inverse_steps": {"values": [20, 25, 30]},
            "inverse_step_size": {"values": [0.0205, 0.021, 0.0215]},
            "augmented_count": {"values": [0]},
            "inverse_restarts": {"values": [1]},
            "inverse_restart_noise_scale": {"values": [0.0]},
            "repulsion_weight": {"values": [5.0, 10.0, 15.0]},
            "minimum_distance_scale": {"values": [0.7, 0.75, 0.8]},
            "cell_violation_weight": {"values": [0.0, 0.01]},
            "coordinate_clip_value": {"values": [0.25, 0.5]},
            "rollout_stages": {"values": [1, 2]},
            "rollout_epochs_per_stage": {"values": [1, 2]},
            "rollout_step_stride": {"values": [5, 10, 25]},
            "replay_buffer_capacity": {"values": [32, 64]},
            "replay_sample_size": {"values": [4, 8]},
            "seed": {"values": [11]},
        }
    if args.sweep_profile == "plan-frontier":
        parameters = {
            "carbon_count": {"values": [-1]},
            "hidden_dim": {"values": [96, 128, 160, 192]},
            "num_message_layers": {"values": [3, 4, 5]},
            "learning_rate": {"values": [5.0e-4, 2.5e-4]},
            "model_lr_decay_rate": {"values": [5.0e-3, 1.0e-3]},
            "smooth_cutoff_width": {"values": [0.15, 0.2]},
            "component_weight_2body": {"values": [3.0, 4.0, 5.0, 6.0]},
            "component_weight_3body": {"values": [0.0, 1.0]},
            "component_weight_4body": {"values": [0.0, 1.0]},
            "batch_size": {"values": [4, 8]},
            "fingerprint_loss_weight": {"values": [0.25, 0.5, 1.0]},
            "target_vertex_weight": {"values": [0.0]},
            "target_position_weight": {"values": [0.0]},
            "inverse_steps": {"values": [100, 200, 400]},
            "inverse_step_size": {"values": [0.02, 0.05, 0.1]},
            "augmented_count": {"values": [0]},
            "inverse_restarts": {"values": [0, 1]},
            "inverse_restart_noise_scale": {"values": [0.0, 0.0025, 0.005]},
            "repulsion_weight": {"values": [5.0, 10.0, 20.0]},
            "minimum_distance_scale": {"values": [0.7, 0.75, 0.8]},
            "cell_violation_weight": {"values": [0.0, 0.01]},
            "coordinate_clip_value": {"values": [0.25, 0.5, 1.0]},
            "rollout_stages": {"values": [1, 2]},
            "rollout_epochs_per_stage": {"values": [1, 2]},
            "rollout_step_stride": {"values": [10, 25, 50]},
            "replay_buffer_capacity": {"values": [32, 64, 128]},
            "replay_sample_size": {"values": [4, 8, 16]},
            "seed": {"values": [11, 42, 101]},
        }
    if args.sweep_profile == "plan-attention-refine":
        parameters = {
            "carbon_count": {"values": [-1]},
            "epochs": {"values": [15]},
            "epoch_values": {"values": ["0,5,10,15"]},
            "hidden_dim": {"values": [96, 128, 160]},
            "num_message_layers": {"values": [3, 4, 5]},
            "learning_rate": {"values": [5.0e-4, 2.5e-4]},
            "model_lr_decay_rate": {"values": [5.0e-3, 1.0e-3]},
            "smooth_cutoff_width": {"values": [0.15, 0.2]},
            "component_weight_2body": {"values": [5.0, 6.0]},
            "component_weight_3body": {"values": [1.0]},
            "component_weight_4body": {"values": [1.0, 2.0]},
            "batch_size": {"values": [4, 8]},
            "fingerprint_loss_weight": {"values": [0.25, 0.275, 0.3, 0.5, 1.0]},
            "target_vertex_weight": {"values": [0.0]},
            "target_position_weight": {"values": [0.0]},
            "inverse_steps": {"values": [30]},
            "inverse_step_values": {"values": ["0,20,25,30"]},
            "inverse_step_size": {"values": [0.021]},
            "step_size_values": {"values": ["0.018,0.0205,0.021,0.0215,0.022"]},
            "augmented_count": {"values": [0]},
            "inverse_restarts": {"values": [0, 1]},
            "inverse_restart_noise_scale": {"values": [0.0, 0.0025]},
            "repulsion_weight": {"values": [5.0, 10.0, 15.0]},
            "minimum_distance_scale": {"values": [0.7, 0.75, 0.8]},
            "cell_violation_weight": {"values": [0.0, 0.01]},
            "coordinate_clip_value": {"values": [0.25, 0.5]},
            "rollout_stages": {"values": [1, 2]},
            "rollout_epochs_per_stage": {"values": [1, 2]},
            "rollout_step_stride": {"values": [5, 10, 25]},
            "replay_buffer_capacity": {"values": [32, 64]},
            "replay_sample_size": {"values": [4, 8]},
            "seed": {"values": [11, 42, 101]},
        }
    if args.sweep_profile == "architecture-frontier":
        parameters = {
            "architecture": {
                "values": [
                    "torch_gnn_residual",
                    "torch_gnn_attention",
                    "torch_gnn_graph_transformer",
                    "torch_gnn_graph_operator",
                    "torch_gnn_multkan",
                ]
            },
            "carbon_count": {"values": [-1]},
            "hidden_dim": {"values": [96, 128, 160]},
            "num_message_layers": {"values": [2, 3, 4]},
            "learning_rate": {"values": [5.0e-4, 2.5e-4]},
            "model_lr_decay_rate": {"values": [5.0e-3, 1.0e-3]},
            "smooth_cutoff_width": {"values": [0.15, 0.2]},
            "component_weight_2body": {"values": [3.0, 5.0]},
            "component_weight_3body": {"values": [0.0, 1.0]},
            "component_weight_4body": {"values": [0.0, 1.0]},
            "batch_size": {"values": [4, 8]},
            "fingerprint_loss_weight": {"values": [0.25, 0.5, 1.0]},
            "target_vertex_weight": {"values": [0.0]},
            "target_position_weight": {"values": [0.0]},
            "inverse_steps": {"values": [25, 50, 100]},
            "inverse_step_size": {"values": [0.01, 0.02, 0.05]},
            "augmented_count": {"values": [0]},
            "inverse_restarts": {"values": [0, 1]},
            "inverse_restart_noise_scale": {"values": [0.0, 0.0025]},
            "repulsion_weight": {"values": [5.0, 10.0, 20.0]},
            "minimum_distance_scale": {"values": [0.7, 0.75, 0.8]},
            "cell_violation_weight": {"values": [0.0, 0.01]},
            "coordinate_clip_value": {"values": [0.25, 0.5]},
            "rollout_stages": {"values": [1, 2]},
            "rollout_epochs_per_stage": {"values": [1, 2]},
            "rollout_step_stride": {"values": [5, 10, 25]},
            "replay_buffer_capacity": {"values": [32, 64]},
            "replay_sample_size": {"values": [4, 8]},
            "seed": {"values": [11, 42, 101]},
        }
    return {
        "name": (
            f"multi-architecture-{args.sweep_profile}-carbon-sweep"
            if args.sweep_profile == "architecture-frontier"
            else f"{args.architecture}-{args.sweep_profile}-carbon-sweep"
        ),
        "method": "bayes",
        "metric": {"name": "inverse/final_rmsd", "goal": "minimize"},
        "parameters": parameters,
    }


def resolve_existing_sweep(sweep_id: str) -> dict[str, str]:
    api = wandb.Api(overrides={"project": WANDB_PROJECT})
    try:
        sweep = api.sweep(str(sweep_id))
    except Exception as exc:
        raise ValueError(
            f"Sweep id '{sweep_id}' could not be resolved in W&B project "
            f"'{WANDB_PROJECT}': {exc}"
        ) from exc
    if str(sweep.project) != WANDB_PROJECT:
        raise ValueError(
            f"Sweep id '{sweep_id}' belongs to project '{sweep.project}', "
            f"not '{WANDB_PROJECT}'."
        )
    return {
        "entity": str(sweep.entity),
        "project": str(sweep.project),
        "sweep_id": str(sweep.id),
        "url": str(sweep.url),
    }


def launch_sweep_agent(
    *,
    sweep_id: str,
    config: dict,
    sweep_count: int,
    project: str,
    entity: str | None = None,
) -> None:
    def agent() -> None:
        execute_run(config=dict(config), sweep_run=True)

    agent_kwargs = {
        "function": agent,
        "project": project,
        "count": int(sweep_count),
    }
    if entity:
        agent_kwargs["entity"] = entity
    wandb.agent(sweep_id, **agent_kwargs)


def resolve_epoch_values(config: dict) -> list[int]:
    if config["epoch_values"]:
        return parse_int_list(str(config["epoch_values"]))
    return default_epoch_values(int(config["epochs"]))


def resolve_inverse_step_values(config: dict) -> list[int]:
    if config["inverse_step_values"]:
        return parse_int_list(str(config["inverse_step_values"]))
    return default_inverse_step_values(int(config["inverse_steps"]))


def resolve_step_size_values(config: dict) -> list[float]:
    if config["step_size_values"]:
        return parse_float_list(str(config["step_size_values"]))
    return default_step_sizes(float(config["inverse_step_size"]))


def resolve_tags(config: dict, sweep_run: bool) -> list[str]:
    tags = list(DEFAULT_TAGS)
    tags.extend(parse_tags(str(config.get("variant_tags", ""))))
    if sweep_run:
        tags.append("sweep")
    return sorted(set(tags))


def resolve_output_dir(base_output_dir: Path, architecture: str, run_id: str) -> Path:
    return base_output_dir / architecture / run_id


def log_series_tables(metrics: dict) -> None:
    payload = {
        "tables/epoch_sweep": wandb.Table(
            columns=["epochs", "position_difference"],
            data=[
                [entry["epochs"], entry["position_difference"]]
                for entry in metrics["epoch_sweep"]
            ],
        ),
        "tables/inverse_step_sweep": wandb.Table(
            columns=["step", "position_difference", "fingerprint_mse", "fingerprint_l2"],
            data=[
                [
                    entry["step"],
                    entry["position_difference"],
                    entry.get("fingerprint_mse"),
                    entry.get("fingerprint_l2"),
                ]
                for entry in metrics["inverse_step_sweep"]
            ],
        ),
        "tables/step_size_sweep": wandb.Table(
            columns=["step_size", "position_difference"],
            data=[
                [entry["step_size"], entry["position_difference"]]
                for entry in metrics["step_size_sweep"]
            ],
        ),
    }
    descriptor_report = metrics.get("descriptor_comparisons", {}).get("final")
    if descriptor_report:
        component_rows = []
        for component_name, comparisons in descriptor_report.get("component_summaries", {}).items():
            for comparison_name, summary in comparisons.items():
                component_rows.append(
                    [
                        component_name,
                        comparison_name,
                        summary["mae"],
                        summary["rmse"],
                        summary["mse"],
                        summary["l2"],
                        summary["max_abs_error"],
                    ]
                )
        if component_rows:
            payload["tables/final_descriptor_components"] = wandb.Table(
                columns=[
                    "component",
                    "comparison",
                    "mae",
                    "rmse",
                    "mse",
                    "l2",
                    "max_abs_error",
                ],
                data=component_rows,
            )

        value_rows = []
        component_dimensions = descriptor_report.get("component_dimensions", {})
        start = 0
        for component_name in ("2body", "3body", "4body"):
            width = int(component_dimensions.get(component_name, 0))
            end = start + width
            target_values = descriptor_report["components"][component_name]["target_raffle"]
            predicted_values = descriptor_report["components"][component_name]["ml_predicted"]
            true_values = descriptor_report["components"][component_name]["true_raffle"]
            for index in range(width):
                value_rows.append(
                    [
                        start + index,
                        component_name,
                        target_values[index],
                        predicted_values[index],
                        true_values[index],
                        predicted_values[index] - target_values[index],
                        true_values[index] - target_values[index],
                        predicted_values[index] - true_values[index],
                    ]
                )
            start = end
        if value_rows:
            payload["tables/final_descriptor_values"] = wandb.Table(
                columns=[
                    "fingerprint_index",
                    "component",
                    "target_raffle",
                    "ml_predicted",
                    "true_raffle",
                    "ml_predicted_minus_target",
                    "true_raffle_minus_target",
                    "ml_predicted_minus_true_raffle",
                ],
                data=value_rows,
            )
    if metrics.get("checkpoint_step_size_sweep"):
        payload["tables/checkpoint_step_size_sweep"] = wandb.Table(
            columns=["epochs", "step_size", "position_difference"],
            data=[
                [entry["epochs"], entry["step_size"], entry["position_difference"]]
                for entry in metrics["checkpoint_step_size_sweep"]
            ],
        )
    if metrics.get("checkpoint_step_schedule_sweep"):
        payload["tables/checkpoint_step_schedule_sweep"] = wandb.Table(
            columns=["epochs", "num_steps", "step_size", "position_difference"],
            data=[
                [
                    entry["epochs"],
                    entry["num_steps"],
                    entry["step_size"],
                    entry["position_difference"],
                ]
                for entry in metrics["checkpoint_step_schedule_sweep"]
            ],
        )
    rollout = metrics.get("rollout", {})
    if rollout.get("stages"):
        payload["tables/rollout_stages"] = wandb.Table(
            columns=[
                "stage_index",
                "num_stage_samples",
                "sampled_replay_count",
                "buffer_size",
                "best_true_target_fingerprint_mse",
                "final_true_target_fingerprint_mse",
                "best_position_difference",
                "final_position_difference",
                "mean_fingerprint_drift_mse",
            ],
            data=[
                [
                    entry["stage_index"],
                    entry["num_stage_samples"],
                    entry["sampled_replay_count"],
                    entry["buffer_size"],
                    entry["best_true_target_fingerprint_mse"],
                    entry["final_true_target_fingerprint_mse"],
                    entry["best_position_difference"],
                    entry["final_position_difference"],
                    entry["mean_fingerprint_drift_mse"],
                ]
                for entry in rollout["stages"]
            ],
        )
    wandb.log(payload)


def log_artifacts(run: wandb.sdk.wandb_run.Run, metrics: dict) -> None:
    wandb.log({"plots/position_sweeps": wandb.Image(metrics["output_files"]["plot"])})
    artifact = wandb.Artifact(
        name=f"{metrics['architecture_name']}-{run.id}-outputs",
        type="inverse-design-results",
    )
    for file_path in metrics["output_files"].values():
        artifact.add_file(file_path)
    run.log_artifact(artifact)


def execute_run(config: dict, sweep_run: bool = False) -> dict:
    validate_plan_constraints(config)
    architecture = str(config["architecture"])
    tags = resolve_tags(config, sweep_run=sweep_run)
    group = architecture
    with wandb.init(
        project=WANDB_PROJECT,
        group=group,
        tags=tags,
        job_type="sweep" if sweep_run else "benchmark",
        name=str(config["run_name"]) if config.get("run_name") else None,
        mode="online",
        config=config,
    ) as run:
        run_config = dict(run.config)
        repo_root = Path(__file__).resolve().parents[2]
        output_dir = resolve_output_dir(
            Path(str(run_config["output_dir"])),
            architecture=architecture,
            run_id=run.id,
        )

        def training_observer(epoch: int, loss: float) -> None:
            wandb.log({"train/epoch": epoch, "train/loss": loss}, step=epoch)

        metrics = run_example(
            repo_root=repo_root,
            output_dir=output_dir,
            carbon_count=int(run_config["carbon_count"]),
            augmented_count=int(run_config["augmented_count"]),
            num_epochs=int(run_config["epochs"]),
            batch_size=int(run_config["batch_size"]),
            inverse_steps=int(run_config["inverse_steps"]),
            inverse_step_size=float(run_config["inverse_step_size"]),
            epoch_values=resolve_epoch_values(run_config),
            inverse_step_values=resolve_inverse_step_values(run_config),
            step_size_values=resolve_step_size_values(run_config),
            fixed_leading_atoms=int(run_config["fixed_leading_atoms"]),
            fingerprint_loss_weight=float(run_config["fingerprint_loss_weight"]),
            target_vertex_weight=float(run_config["target_vertex_weight"]),
            target_position_weight=float(run_config["target_position_weight"]),
            inverse_lr_decay_rate=float(run_config["inverse_lr_decay_rate"]),
            inverse_restarts=int(run_config["inverse_restarts"]),
            inverse_restart_noise_scale=float(run_config["inverse_restart_noise_scale"]),
            repulsion_weight=float(run_config.get("repulsion_weight", REPULSION_WEIGHT)),
            minimum_distance_scale=float(
                run_config.get("minimum_distance_scale", MINIMUM_DISTANCE_SCALE)
            ),
            cell_violation_weight=float(
                run_config.get("cell_violation_weight", CELL_VIOLATION_WEIGHT)
            ),
            coordinate_clip_value=(
                None
                if run_config.get("coordinate_clip_value", COORDINATE_CLIP_VALUE) is None
                else float(run_config["coordinate_clip_value"])
            ),
            rollout_stages=int(run_config.get("rollout_stages", ROLLOUT_STAGES)),
            rollout_epochs_per_stage=int(
                run_config.get("rollout_epochs_per_stage", ROLLOUT_EPOCHS_PER_STAGE)
            ),
            rollout_step_stride=int(run_config.get("rollout_step_stride", ROLLOUT_STEP_STRIDE)),
            replay_buffer_capacity=int(
                run_config.get("replay_buffer_capacity", REPLAY_BUFFER_CAPACITY)
            ),
            replay_sample_size=int(run_config.get("replay_sample_size", REPLAY_SAMPLE_SIZE)),
            replay_category_weights=parse_category_weights(
                str(
                    run_config.get(
                        "replay_category_weights",
                        ",".join(
                            f"{key}={value}"
                            for key, value in DEFAULT_REPLAY_CATEGORY_WEIGHTS.items()
                        ),
                    )
                )
            ),
            rollout_drift_threshold=float(
                run_config.get("rollout_drift_threshold", ROLLOUT_DRIFT_THRESHOLD)
            ),
            rollout_high_error_threshold=float(
                run_config.get("rollout_high_error_threshold", ROLLOUT_HIGH_ERROR_THRESHOLD)
            ),
            rollout_instability_threshold=float(
                run_config.get("rollout_instability_threshold", ROLLOUT_INSTABILITY_THRESHOLD)
            ),
            seed=int(run_config["seed"]),
            model_config=build_model_config(run_config),
            training_observer=training_observer,
            architecture_name=architecture,
            enable_checkpoint_step_size_sweep=not bool(
                run_config.get("skip_checkpoint_step_size_sweep", False)
            ),
            enable_checkpoint_step_schedule_sweep=not bool(
                run_config.get("skip_checkpoint_step_schedule_sweep", False)
            ),
        )

        descriptor_report = metrics.get("descriptor_comparisons", {}).get("final")
        descriptor_payload = {}
        if descriptor_report:
            global_metrics = descriptor_report["global_metrics"]
            descriptor_payload.update(
                {
                    "descriptor/final_predicted_vs_target_mae": global_metrics[
                        "ml_predicted_vs_target"
                    ]["mae"],
                    "descriptor/final_predicted_vs_target_rmse": global_metrics[
                        "ml_predicted_vs_target"
                    ]["rmse"],
                    "descriptor/final_true_vs_target_mae": global_metrics[
                        "true_raffle_vs_target"
                    ]["mae"],
                    "descriptor/final_true_vs_target_rmse": global_metrics[
                        "true_raffle_vs_target"
                    ]["rmse"],
                    "descriptor/final_predicted_vs_true_mae": global_metrics[
                        "ml_predicted_vs_true_raffle"
                    ]["mae"],
                    "descriptor/final_predicted_vs_true_rmse": global_metrics[
                        "ml_predicted_vs_true_raffle"
                    ]["rmse"],
                }
            )
            for component_name, component_summary in descriptor_report.get(
                "component_summaries",
                {},
            ).items():
                descriptor_payload[
                    f"descriptor/{component_name}_predicted_vs_target_mae"
                ] = component_summary["ml_predicted_vs_target"]["mae"]
                descriptor_payload[
                    f"descriptor/{component_name}_predicted_vs_target_rmse"
                ] = component_summary["ml_predicted_vs_target"]["rmse"]
                descriptor_payload[
                    f"descriptor/{component_name}_true_vs_target_mae"
                ] = component_summary["true_raffle_vs_target"]["mae"]
                descriptor_payload[
                    f"descriptor/{component_name}_true_vs_target_rmse"
                ] = component_summary["true_raffle_vs_target"]["rmse"]

        wandb.log(
            {
                "inverse/initial_fingerprint_mse": metrics["initial_fingerprint_mse"],
                "inverse/final_fingerprint_mse": metrics["final_fingerprint_mse"],
                "inverse/initial_rmsd": metrics["initial_rmsd"],
                "inverse/final_rmsd": metrics["final_rmsd"],
                "inverse/rmsd_reduction_fraction": metrics["rmsd_reduction_fraction"],
                "convergence/best_step": metrics["convergence_summary"]["best_step"],
                "convergence/tail_position_range": metrics["convergence_summary"]["tail_position_range"],
                "convergence/tail_fingerprint_range": metrics["convergence_summary"]["tail_fingerprint_range"],
                "convergence/fingerprint_nonincreasing_fraction": metrics["convergence_summary"]["fingerprint_nonincreasing_fraction"],
                "convergence/position_nonincreasing_fraction": metrics["convergence_summary"]["position_nonincreasing_fraction"],
                "convergence/converges_to_best_within_5pct": float(
                    metrics["convergence_summary"]["converges_to_best_within_5pct"]
                ),
                "rollout/enabled": float(metrics["rollout"]["enabled"]),
                "rollout/buffer_size": metrics["rollout"]["replay_buffer"]["size"],
                "rollout/mean_fingerprint_drift_mse": metrics["rollout"]["replay_buffer"]["mean_fingerprint_drift_mse"],
                **descriptor_payload,
            }
        )
        log_series_tables(metrics)
        log_artifacts(run, metrics)
        run.summary.update(
            {
                "architecture_name": metrics["architecture_name"],
                "initial_rmsd": metrics["initial_rmsd"],
                "configured_final_rmsd": metrics.get("configured_final_rmsd"),
                "final_rmsd": metrics["final_rmsd"],
                "initial_fingerprint_mse": metrics["initial_fingerprint_mse"],
                "final_fingerprint_mse": metrics["final_fingerprint_mse"],
                "descriptor_comparison": (
                    None
                    if descriptor_report is None
                    else {
                        "global_metrics": descriptor_report["global_metrics"],
                        "report_file": descriptor_report["report_file"],
                        "structure_file": descriptor_report["structure_file"],
                    }
                ),
                "convergence_summary": metrics["convergence_summary"],
                "rollout": metrics["rollout"],
                "rmsd_reduction_fraction": metrics["rmsd_reduction_fraction"],
                "selected_candidate": metrics.get("selected_candidate"),
                "output_dir": str(output_dir),
            }
        )
        print(
            json.dumps(
                {
                    "run_id": run.id,
                    "run_name": run.name,
                    "project": WANDB_PROJECT,
                    "group": group,
                    "tags": tags,
                    "final_rmsd": metrics["final_rmsd"],
                    "output_dir": str(output_dir),
                },
                indent=2,
            )
        )
        return metrics


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--architecture", type=str, default=DEFAULT_ARCHITECTURE)
    parser.add_argument("--variant-tags", type=str, default="baseline")
    parser.add_argument("--run-name", type=str, default="")
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
    parser.add_argument("--target-position-weight", type=float, default=0.0)
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
        default=DEFAULT_MODEL_CONFIG["component_weight"][0],
    )
    parser.add_argument(
        "--component-weight-3body",
        type=float,
        default=DEFAULT_MODEL_CONFIG["component_weight"][1],
    )
    parser.add_argument(
        "--component-weight-4body",
        type=float,
        default=DEFAULT_MODEL_CONFIG["component_weight"][2],
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("build") / "wandb_inverse_design",
    )
    parser.add_argument("--skip-checkpoint-step-size-sweep", action="store_true")
    parser.add_argument("--skip-checkpoint-step-schedule-sweep", action="store_true")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--print-sweep-config", action="store_true")
    sweep_group = parser.add_mutually_exclusive_group()
    sweep_group.add_argument("--launch-sweep", action="store_true")
    sweep_group.add_argument(
        "--sweep-id",
        type=str,
        help="Attach to an existing W&B sweep in this project instead of creating a new one.",
    )
    parser.add_argument("--sweep-count", type=int, default=8)
    parser.add_argument("--sweep-profile", type=str, default="broad")
    return parser.parse_args(argv)


def main() -> None:
    args = parse_args()
    if args.print_sweep_config:
        print(json.dumps(build_sweep_config(args), indent=2))
        return

    config = vars(args).copy()
    config["output_dir"] = str(args.output_dir)
    config.pop("sweep_id", None)

    if args.sweep_id:
        try:
            existing_sweep = resolve_existing_sweep(args.sweep_id)
        except ValueError as exc:
            raise SystemExit(str(exc)) from exc
        print(json.dumps({"action": "continue", **existing_sweep}, indent=2))
        launch_sweep_agent(
            sweep_id=existing_sweep["sweep_id"],
            config=config,
            sweep_count=int(args.sweep_count),
            project=existing_sweep["project"],
            entity=existing_sweep["entity"],
        )
        return

    if args.launch_sweep:
        sweep_config = build_sweep_config(args)
        sweep_id = wandb.sweep(sweep=sweep_config, project=WANDB_PROJECT)
        print(json.dumps({"project": WANDB_PROJECT, "sweep_id": sweep_id}, indent=2))
        launch_sweep_agent(
            sweep_id=str(sweep_id),
            config=config,
            sweep_count=int(args.sweep_count),
            project=WANDB_PROJECT,
        )
        return

    execute_run(config=config, sweep_run=False)


if __name__ == "__main__":
    main()
