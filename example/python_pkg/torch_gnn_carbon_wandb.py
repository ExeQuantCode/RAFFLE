"""W&B online runner and sweep entrypoint for torch carbon inverse design."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import wandb

from torch_gnn_carbon_workflow_example import (
    DEFAULT_MODEL_CONFIG,
    FINGERPRINT_LOSS_WEIGHT,
    FIXED_LEADING_ATOMS,
    INVERSE_LR_DECAY_RATE,
    INVERSE_RESTART_NOISE_SCALE,
    INVERSE_RESTARTS,
    TARGET_VERTEX_WEIGHT,
    default_epoch_values,
    default_inverse_step_values,
    default_step_sizes,
    parse_float_list,
    parse_int_list,
    run_example,
)


WANDB_PROJECT = "raffle-inverse-design"
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


def build_sweep_config(args: argparse.Namespace) -> dict:
    parameters = {
        "hidden_dim": {"values": [64, 80, 96, 128]},
        "num_message_layers": {"values": [2, 3, 4]},
        "learning_rate": {"values": [1.0e-3, 5.0e-4, 2.5e-4]},
        "model_lr_decay_rate": {"values": [1.0e-2, 5.0e-3, 1.0e-3]},
        "smooth_cutoff_width": {"values": [0.1, 0.15, 0.2, 0.3]},
        "component_weight_2body": {"values": [2.0, 4.0, 6.0, 8.0]},
        "component_weight_3body": {"values": [1.0, 2.0]},
        "component_weight_4body": {"values": [1.0, 2.0]},
        "batch_size": {"values": [4, 8, 12]},
        "fingerprint_loss_weight": {"values": [0.25, 0.5, 0.75, 1.0]},
        "target_vertex_weight": {"values": [0.0, 0.25, 0.5, 0.75]},
        "inverse_steps": {"values": [200, 300, 400]},
        "inverse_step_size": {"values": [1.0e-3, 2.5e-3, 5.0e-3]},
        "augmented_count": {"values": [8, 16, 24]},
        "inverse_restarts": {"values": [1, 2, 4]},
        "inverse_restart_noise_scale": {"values": [0.0, 0.005, 0.01]},
        "seed": {"values": [11, 42, 101]},
    }
    if args.sweep_profile == "vertex-focus":
        parameters = {
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
            "target_vertex_weight": {"values": [0.5, 1.0, 2.0, 4.0]},
            "inverse_steps": {"values": [400, 800, 1200]},
            "inverse_step_size": {"values": [5.0e-4, 1.0e-3, 2.5e-3]},
            "augmented_count": {"values": [16, 24, 32]},
            "inverse_restarts": {"values": [1, 2, 4]},
            "inverse_restart_noise_scale": {"values": [0.0, 0.005, 0.01]},
            "seed": {"values": [11, 42, 101]},
        }
    if args.sweep_profile == "inverse-basin":
        parameters = {
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
            "target_vertex_weight": {"values": [0.5, 0.55, 0.6]},
            "inverse_steps": {"values": [20, 25, 30]},
            "inverse_step_size": {"values": [0.0205, 0.021, 0.0215]},
            "augmented_count": {"values": [16]},
            "inverse_restarts": {"values": [1]},
            "inverse_restart_noise_scale": {"values": [0.0]},
            "seed": {"values": [11]},
        }
    return {
        "name": f"{args.architecture}-{args.sweep_profile}-carbon-sweep",
        "method": "bayes",
        "metric": {"name": "inverse/final_rmsd", "goal": "minimize"},
        "parameters": parameters,
    }


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
            columns=["step", "position_difference"],
            data=[
                [entry["step"], entry["position_difference"]]
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
            seed=int(run_config["seed"]),
            model_config=build_model_config(run_config),
            training_observer=training_observer,
            architecture_name=architecture,
        )

        wandb.log(
            {
                "inverse/initial_rmsd": metrics["initial_rmsd"],
                "inverse/final_rmsd": metrics["final_rmsd"],
                "inverse/rmsd_reduction_fraction": metrics["rmsd_reduction_fraction"],
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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--architecture", type=str, default=DEFAULT_ARCHITECTURE)
    parser.add_argument("--variant-tags", type=str, default="baseline")
    parser.add_argument("--run-name", type=str, default="")
    parser.add_argument("--carbon-count", type=int, default=0)
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
    parser.add_argument("--target-position-weight", type=float, default=0.0)
    parser.add_argument("--inverse-lr-decay-rate", type=float, default=INVERSE_LR_DECAY_RATE)
    parser.add_argument("--inverse-restarts", type=int, default=INVERSE_RESTARTS)
    parser.add_argument(
        "--inverse-restart-noise-scale",
        type=float,
        default=INVERSE_RESTART_NOISE_SCALE,
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
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--print-sweep-config", action="store_true")
    parser.add_argument("--launch-sweep", action="store_true")
    parser.add_argument("--sweep-count", type=int, default=8)
    parser.add_argument("--sweep-profile", type=str, default="broad")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.print_sweep_config:
        print(json.dumps(build_sweep_config(args), indent=2))
        return

    config = vars(args).copy()
    config["output_dir"] = str(args.output_dir)

    if args.launch_sweep:
        sweep_config = build_sweep_config(args)
        sweep_id = wandb.sweep(sweep=sweep_config, project=WANDB_PROJECT)
        print(json.dumps({"project": WANDB_PROJECT, "sweep_id": sweep_id}, indent=2))

        def agent() -> None:
            execute_run(config=dict(config), sweep_run=True)

        wandb.agent(sweep_id, function=agent, project=WANDB_PROJECT, count=int(args.sweep_count))
        return

    execute_run(config=config, sweep_run=False)


if __name__ == "__main__":
    main()
