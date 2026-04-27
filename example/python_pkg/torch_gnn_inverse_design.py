"""Run inverse design from a saved PyTorch GNN surrogate model checkpoint."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np


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
        **inverse_design_options,
    )

    metrics = compute_inverse_design_metrics(
        model=model,
        initial_atoms=input_atoms,
        optimised_atoms=optimised,
        target_fingerprint=target_fingerprint,
        target_atoms=target_atoms,
    )
    metrics.update(
        {
            "model_checkpoint": str(args.model_checkpoint.resolve()),
            "target_fingerprint_file": str(args.target_fingerprint.resolve()),
            "input_structure_file": str(args.input_structure.resolve()),
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
    write_json(metrics_path, metrics)
    write_inverse_design_log(log_path, metrics)

    if target_atoms is not None:
        print_position_differences(target_atoms, input_atoms, optimised)
        print()
    print(f"Saved optimised structure: {structure_path}")
    print(f"Saved metrics log: {log_path}")
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
