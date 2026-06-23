"""Run inverse design from a saved PyTorch GNN surrogate checkpoint.

This is the supported split-workflow entry point for replaying inverse design from
an input structure plus either a saved target fingerprint or a target structure.
It writes the final structure, descriptor comparison report, and machine-readable
metrics bundle into the requested output directory.
"""

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
    save_descriptor_comparison_report,
    structures_have_matching_atom_count,
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
        atoms_snapshot.info["inverse_step"] = int(step_record["step"])
        atoms_snapshot.info["inverse_num_steps"] = int(step_record["num_steps"])
        atoms_snapshot.info["inverse_is_initial_state"] = bool(step_record["is_initial_state"])
        trajectory_frames.append(atoms_snapshot)
        optimisation_history.append(
            {
                "step": int(step_record["step"]),
                "num_steps": int(step_record["num_steps"]),
                "is_initial_state": bool(step_record["is_initial_state"]),
            }
        )

    return observer


def save_optimisation_path(
    output_dir: Path,
    trajectory_frames: list,
    optimisation_history: list[dict[str, int | bool]],
) -> tuple[str | None, list[dict[str, int | bool]]]:
    if not trajectory_frames:
        return None, optimisation_history
    traj_path = output_dir / "torch_gnn_inverse_design_path.traj"
    write(traj_path, trajectory_frames)
    return str(traj_path), optimisation_history


def parse_fixed_atom_indices(value: str) -> list[int]:
    text = str(value).strip()
    if not text:
        return []
    return [int(item.strip()) for item in text.split(",") if item.strip()]


def build_fixed_atoms_mask(
    num_atoms: int,
    *,
    fixed_leading_atoms: int,
    fixed_atom_indices: list[int] | None = None,
) -> tuple[np.ndarray, list[int]]:
    resolved_indices = set(range(max(min(int(fixed_leading_atoms), int(num_atoms)), 0)))
    for atom_index in fixed_atom_indices or []:
        resolved_index = int(atom_index)
        if resolved_index < 0 or resolved_index >= int(num_atoms):
            raise ValueError(
                f"fixed atom index {resolved_index} is out of bounds for a structure "
                f"with {num_atoms} atoms"
            )
        resolved_indices.add(resolved_index)
    ordered_indices = sorted(resolved_indices)
    fixed_atoms = np.zeros(int(num_atoms), dtype=bool)
    if ordered_indices:
        fixed_atoms[np.asarray(ordered_indices, dtype=np.int64)] = True
    return fixed_atoms, ordered_indices


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


def resolve_target_fingerprint(
    model,
    *,
    target_fingerprint_path: Path | None,
    target_atoms,
    target_structure_path: Path | None,
    target_structure_index: int,
) -> tuple[np.ndarray, dict[str, int | str]]:
    if target_fingerprint_path is not None:
        resolved_path = target_fingerprint_path.resolve()
        return load_target_fingerprint(resolved_path), {
            "type": "file",
            "path": str(resolved_path),
        }
    if target_atoms is None or target_structure_path is None:
        raise ValueError(
            "Provide --target-fingerprint or --target-structure so the target descriptor "
            "can be resolved."
        )
    fingerprint = np.asarray(
        model.compute_reference_fingerprint(target_atoms),
        dtype=np.float32,
    ).reshape(-1)
    return fingerprint, {
        "type": "target_structure",
        "path": str(target_structure_path.resolve()),
        "index": int(target_structure_index),
    }


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    repo_root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--input-structure",
        "--reference-structure",
        dest="input_structure",
        type=Path,
        required=True,
        help="Starting structure to optimise during inverse design.",
    )
    parser.add_argument(
        "--input-structure-index",
        "--reference-structure-index",
        dest="input_structure_index",
        type=int,
        default=0,
        help="Frame index to read when --input-structure contains multiple structures.",
    )
    parser.add_argument(
        "--target-fingerprint",
        type=Path,
        default=None,
        help=(
            "Optional saved target fingerprint (.npy or .json). If omitted, the script "
            "computes the analytical target descriptor from --target-structure."
        ),
    )
    parser.add_argument(
        "--model-checkpoint",
        type=Path,
        required=True,
        help="Checkpoint written by torch_gnn_train_model.py or the W&B replay/export flow.",
    )
    parser.add_argument(
        "--target-structure",
        type=Path,
        default=None,
        help=(
            "Structure used for evaluation metrics and, when --target-fingerprint is not "
            "provided, for the analytical target descriptor."
        ),
    )
    parser.add_argument(
        "--target-structure-index",
        type=int,
        default=0,
        help="Frame index to read when --target-structure contains multiple structures.",
    )
    parser.add_argument(
        "--inverse-steps",
        type=int,
        default=400,
        help="Number of gradient-based inverse-design optimisation steps.",
    )
    parser.add_argument(
        "--inverse-step-size",
        type=float,
        default=5.0e-3,
        help="Initial inverse-design step size.",
    )
    parser.add_argument(
        "--fixed-leading-atoms",
        type=int,
        default=0,
        help="Freeze the first N atoms before applying any explicit --fixed-atoms list.",
    )
    parser.add_argument(
        "--fixed-atoms",
        "--fix-atoms",
        dest="fixed_atoms",
        type=parse_fixed_atom_indices,
        default=[],
        help=(
            "Comma-separated atom indices to keep fixed during inverse design. "
            "Combined with --fixed-leading-atoms."
        ),
    )
    parser.add_argument(
        "--fingerprint-loss-weight",
        type=float,
        default=1.0,
        help="Weight on the fingerprint-matching loss during inverse design.",
    )
    parser.add_argument("--target-vertex-weight", type=float, default=0.0, help=argparse.SUPPRESS)
    parser.add_argument("--target-position-weight", type=float, default=0.0, help=argparse.SUPPRESS)
    parser.add_argument(
        "--inverse-lr-decay-rate",
        type=float,
        default=0.0,
        help="Optional exponential decay applied to the inverse-design step size.",
    )
    parser.add_argument("--inverse-restarts", type=int, default=1, help=argparse.SUPPRESS)
    parser.add_argument(
        "--inverse-restart-noise-scale",
        type=float,
        default=0.0,
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--repulsion-weight",
        type=float,
        default=10.0,
        help="Penalty weight that discourages atoms from moving too close together.",
    )
    parser.add_argument(
        "--minimum-distance-scale",
        type=float,
        default=0.75,
        help="Minimum allowed distance as a fraction of the covalent-radius sum.",
    )
    parser.add_argument(
        "--cell-violation-weight",
        type=float,
        default=0.0,
        help="Penalty weight for moving atoms outside the fixed simulation cell.",
    )
    parser.add_argument(
        "--coordinate-clip-value",
        type=float,
        default=None,
        help="Optional absolute clamp applied to each optimisation coordinate update.",
    )
    parser.add_argument(
        "--save-optimisation-traj",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Save the full inverse-design path as a multi-frame .traj file.",
    )
    parser.add_argument(
        "--plot-2body-fingerprint-comparison",
        action="store_true",
        help="Save a target-vs-final 2-body fingerprint comparison plot.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=repo_root / "build" / "torch_gnn_inverse_design",
        help="Directory for the final structure, metrics bundle, and optional trajectory/plots.",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
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
    fixed_atoms, fixed_atom_indices = build_fixed_atoms_mask(
        len(input_atoms),
        fixed_leading_atoms=args.fixed_leading_atoms,
        fixed_atom_indices=args.fixed_atoms,
    )

    target_fingerprint, target_fingerprint_source = resolve_target_fingerprint(
        model,
        target_fingerprint_path=args.target_fingerprint,
        target_atoms=target_atoms,
        target_structure_path=args.target_structure,
        target_structure_index=args.target_structure_index,
    )
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
        fingerprint_loss_weight=args.fingerprint_loss_weight,
        inverse_lr_decay_rate=args.inverse_lr_decay_rate,
        repulsion_weight=args.repulsion_weight,
        minimum_distance_scale=args.minimum_distance_scale,
        cell_violation_weight=args.cell_violation_weight,
        coordinate_clip_value=args.coordinate_clip_value,
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
            "target_fingerprint_file": (
                target_fingerprint_source["path"]
                if target_fingerprint_source["type"] == "file"
                else None
            ),
            "target_fingerprint_source": target_fingerprint_source,
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
                "fixed_atoms": fixed_atom_indices,
                "fingerprint_loss_weight": float(args.fingerprint_loss_weight),
                "target_vertex_weight": float(args.target_vertex_weight),
                "target_position_weight": float(args.target_position_weight),
                "inverse_lr_decay_rate": float(args.inverse_lr_decay_rate),
                "repulsion_weight": float(args.repulsion_weight),
                "minimum_distance_scale": float(args.minimum_distance_scale),
                "cell_violation_weight": float(args.cell_violation_weight),
                "coordinate_clip_value": (
                    None if args.coordinate_clip_value is None else float(args.coordinate_clip_value)
                ),
            },
            "checkpoint_training_config": checkpoint.get("training_config", {}),
        }
    )

    structure_path = output_dir / "torch_gnn_inverse_design_final.xyz"
    metrics_path = output_dir / "torch_gnn_inverse_design_metrics.json"
    log_path = output_dir / "torch_gnn_inverse_design_metrics.log"
    write_structure(structure_path, optimised)
    descriptor_report = save_descriptor_comparison_report(
        model=model,
        target_fingerprint=target_fingerprint,
        final_atoms=optimised,
        structure_path=structure_path,
    )
    metrics["output_files"] = {
        "optimised_structure": str(structure_path),
        "optimised_structure_descriptor_comparison": descriptor_report["report_file"],
        "optimised_structure_descriptor_comparison_plot": descriptor_report["plot_file"],
        "metrics": str(metrics_path),
        "log": str(log_path),
    }
    metrics["descriptor_comparisons"] = {"final": descriptor_report}

    if args.save_optimisation_traj:
        traj_path, optimisation_path = save_optimisation_path(
            output_dir,
            trajectory_frames,
            optimisation_history,
        )
        metrics["optimisation_history"] = optimisation_path
        if traj_path is not None:
            metrics["output_files"]["optimisation_traj"] = traj_path

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

    if target_atoms is not None and structures_have_matching_atom_count(target_atoms, input_atoms):
        print_position_differences(target_atoms, input_atoms, optimised)
        print()
    print(f"Saved optimised structure: {structure_path}")
    print(
        "Saved descriptor comparison: "
        f"{metrics['output_files']['optimised_structure_descriptor_comparison']}"
    )
    print(
        "Saved descriptor comparison plot: "
        f"{metrics['output_files']['optimised_structure_descriptor_comparison_plot']}"
    )
    print(f"Saved metrics log: {log_path}")
    if args.save_optimisation_traj:
        print(f"Saved optimisation trajectory: {metrics['output_files']['optimisation_traj']}")
    if args.plot_2body_fingerprint_comparison:
        print(f"Saved 2-body fingerprint comparison plot: {metrics['output_files']['two_body_fingerprint_plot']}")
    # print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
