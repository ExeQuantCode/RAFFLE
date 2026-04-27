"""Run inverse design from a saved PyTorch GNN surrogate model checkpoint."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np
from ase.io import write


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from torch_gnn_workflow_common import (
    build_inverse_design_options,
    compute_inverse_design_metrics,
    load_model_from_checkpoint,
    load_target_fingerprint,
    print_position_differences,
    read_single_structure,
    write_inverse_design_log,
    write_json,
    write_structure,
)


def build_optimisation_step_observer(
    trajectory_frames: list,
    optimisation_history: list[dict[str, int | bool]],
):
    def observer(step_record: dict[str, object]) -> None:
        atoms_snapshot = step_record["atoms"].copy()
        atoms_snapshot.info["inverse_restart_index"] = int(step_record["restart_index"]) + 1
        atoms_snapshot.info["inverse_num_restarts"] = int(step_record["num_restarts"])
        atoms_snapshot.info["inverse_step"] = int(step_record["step"])
        atoms_snapshot.info["inverse_num_steps"] = int(step_record["num_steps"])
        atoms_snapshot.info["inverse_is_initial_state"] = bool(step_record["is_initial_state"])
        trajectory_frames.append(atoms_snapshot)
        optimisation_history.append(
            {
                "restart_index": int(step_record["restart_index"]) + 1,
                "num_restarts": int(step_record["num_restarts"]),
                "step": int(step_record["step"]),
                "num_steps": int(step_record["num_steps"]),
                "is_initial_state": bool(step_record["is_initial_state"]),
            }
        )

    return observer


def save_2body_fingerprint_comparison_plot(
    target_fingerprint_2body: np.ndarray,
    predicted_fingerprint_2body: np.ndarray,
    output_path: Path,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    indices = np.arange(target_fingerprint_2body.size)
    difference = predicted_fingerprint_2body - target_fingerprint_2body

    figure = plt.figure(figsize=(12, 6))
    ax1 = figure.add_subplot(2, 1, 1)
    ax1.plot(indices, target_fingerprint_2body, label="target 2-body fingerprint", linewidth=2.0)
    ax1.plot(
        indices,
        predicted_fingerprint_2body,
        label="optimised structure inferred 2-body fingerprint",
        linewidth=1.5,
    )
    ax1.set_ylabel("Fingerprint value")
    ax1.set_title("2-body fingerprint comparison")
    ax1.grid(alpha=0.3)
    ax1.legend(loc="best")

    ax2 = figure.add_subplot(2, 1, 2)
    ax2.plot(indices, difference, color="tab:red", linewidth=1.5)
    ax2.axhline(0.0, color="black", linewidth=1.0, linestyle="--")
    ax2.set_xlabel("2-body fingerprint index")
    ax2.set_ylabel("Predicted - target")
    ax2.grid(alpha=0.3)

    figure.tight_layout()
    figure.savefig(output_path, dpi=150)
    plt.close(figure)


def parse_args() -> argparse.Namespace:
    repo_root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-structure", type=Path, required=True)
    parser.add_argument("--input-structure-index", type=int, default=0)
    parser.add_argument("--target-fingerprint", type=Path, required=True)
    parser.add_argument("--model-checkpoint", type=Path, required=True)
    parser.add_argument("--target-structure", type=Path, default=None)
    parser.add_argument("--target-structure-index", type=int, default=0)
    parser.add_argument("--inverse-steps", type=int, default=400)
    parser.add_argument("--inverse-step-size", type=float, default=5.0e-3)
    parser.add_argument("--fixed-leading-atoms", type=int, default=0)
    parser.add_argument("--fingerprint-loss-weight", type=float, default=1.0)
    parser.add_argument("--target-vertex-weight", type=float, default=0.0)
    parser.add_argument("--target-position-weight", type=float, default=0.0)
    parser.add_argument("--inverse-lr-decay-rate", type=float, default=0.0)
    parser.add_argument("--inverse-restarts", type=int, default=1)
    parser.add_argument("--inverse-restart-noise-scale", type=float, default=0.0)
    parser.add_argument("--save-optimisation-traj", action="store_true")
    parser.add_argument("--plot-2body-fingerprint-comparison", action="store_true")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=repo_root / "build" / "torch_gnn_inverse_design",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    model, checkpoint = load_model_from_checkpoint(args.model_checkpoint.resolve())
    input_atoms = read_single_structure(args.input_structure.resolve(), index=args.input_structure_index)
    target_atoms = None
    if args.target_structure is not None:
        target_atoms = read_single_structure(
            args.target_structure.resolve(),
            index=args.target_structure_index,
        )
    fixed_atoms = np.zeros(len(input_atoms), dtype=bool)
    fixed_atoms[: max(int(args.fixed_leading_atoms), 0)] = True

    target_fingerprint = load_target_fingerprint(args.target_fingerprint.resolve())
    if target_fingerprint.size != int(model.fingerprint_dim):
        raise ValueError(
            "Target fingerprint length does not match model fingerprint dimension: "
            f"{target_fingerprint.size} != {model.fingerprint_dim}"
        )

    trajectory_frames: list = []
    optimisation_history: list[dict[str, int | bool]] = []
    step_observer = None
    if args.save_optimisation_traj:
        step_observer = build_optimisation_step_observer(trajectory_frames, optimisation_history)

    inverse_design_options = build_inverse_design_options(
        target_atoms=target_atoms,
        fingerprint_loss_weight=args.fingerprint_loss_weight,
        target_vertex_weight=args.target_vertex_weight,
        target_position_weight=args.target_position_weight,
        inverse_lr_decay_rate=args.inverse_lr_decay_rate,
        inverse_restarts=args.inverse_restarts,
        inverse_restart_noise_scale=args.inverse_restart_noise_scale,
    )
    optimised = model.inverse_design(
        target_fingerprint=target_fingerprint,
        atoms=input_atoms,
        fixed_atoms=fixed_atoms,
        num_steps=args.inverse_steps,
        step_size=args.inverse_step_size,
        verbose=1,
        step_observer=step_observer,
        **inverse_design_options,
    )

    metrics = compute_inverse_design_metrics(
        model=model,
        initial_atoms=input_atoms,
        optimised_atoms=optimised,
        target_fingerprint=target_fingerprint,
        target_atoms=target_atoms,
    )
    predicted_fingerprint_2body, _, _ = model.predict_components(optimised)
    target_fingerprint_2body = target_fingerprint[: model.fingerprint_dim_2body]
    metrics.update(
        {
            "model_checkpoint": str(args.model_checkpoint.resolve()),
            "target_fingerprint_file": str(args.target_fingerprint.resolve()),
            "input_structure_file": str(args.input_structure.resolve()),
            "final_fingerprint_mse_2body": float(
                np.mean((predicted_fingerprint_2body - target_fingerprint_2body) ** 2)
            ),
            "final_fingerprint_l2_2body": float(
                np.linalg.norm(predicted_fingerprint_2body - target_fingerprint_2body)
            ),
            "inverse_design_config": {
                "input_structure_index": int(args.input_structure_index),
                "target_structure": (
                    None if args.target_structure is None else str(args.target_structure.resolve())
                ),
                "target_structure_index": int(args.target_structure_index),
                "inverse_steps": int(args.inverse_steps),
                "inverse_step_size": float(args.inverse_step_size),
                "fixed_leading_atoms": int(args.fixed_leading_atoms),
                "fingerprint_loss_weight": float(args.fingerprint_loss_weight),
                "target_vertex_weight": float(args.target_vertex_weight),
                "target_position_weight": float(args.target_position_weight),
                "inverse_lr_decay_rate": float(args.inverse_lr_decay_rate),
                "inverse_restarts": int(args.inverse_restarts),
                "inverse_restart_noise_scale": float(args.inverse_restart_noise_scale),
            },
            "checkpoint_training_config": checkpoint.get("training_config", {}),
        }
    )

    structure_path = output_dir / "torch_gnn_inverse_design_final.xyz"
    metrics_path = output_dir / "torch_gnn_inverse_design_metrics.json"
    log_path = output_dir / "torch_gnn_inverse_design_metrics.log"
    write_structure(structure_path, optimised)
    metrics["output_files"] = {
        "optimised_structure": str(structure_path),
        "metrics": str(metrics_path),
        "log": str(log_path),
    }

    if args.save_optimisation_traj:
        traj_path = output_dir / "torch_gnn_inverse_design_path.traj"
        write(traj_path, trajectory_frames)
        metrics["optimisation_history"] = optimisation_history
        metrics["output_files"]["optimisation_traj"] = str(traj_path)

    if args.plot_2body_fingerprint_comparison:
        plot_path = output_dir / "torch_gnn_inverse_design_2body_fingerprint.png"
        save_2body_fingerprint_comparison_plot(
            target_fingerprint_2body=target_fingerprint_2body,
            predicted_fingerprint_2body=predicted_fingerprint_2body,
            output_path=plot_path,
        )
        metrics["target_fingerprint_2body"] = target_fingerprint_2body.tolist()
        metrics["predicted_fingerprint_2body"] = predicted_fingerprint_2body.tolist()
        metrics["output_files"]["two_body_fingerprint_plot"] = str(plot_path)

    write_json(metrics_path, metrics)
    write_inverse_design_log(log_path, metrics)

    if target_atoms is not None:
        print_position_differences(target_atoms, input_atoms, optimised)
        print()
    print(f"Saved optimised structure: {structure_path}")
    print(f"Saved metrics log: {log_path}")
    if args.save_optimisation_traj:
        print(f"Saved optimisation trajectory: {metrics['output_files']['optimisation_traj']}")
    if args.plot_2body_fingerprint_comparison:
        print(f"Saved 2-body fingerprint comparison plot: {metrics['output_files']['two_body_fingerprint_plot']}")
    # print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
