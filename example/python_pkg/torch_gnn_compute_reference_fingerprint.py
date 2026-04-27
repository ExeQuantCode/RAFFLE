"""Compute and save the analytical RAFFLE fingerprint for a structure."""

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
    load_model_from_checkpoint,
    read_single_structure,
    save_target_fingerprint,
)


def parse_args() -> argparse.Namespace:
    repo_root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-structure", type=Path, required=True)
    parser.add_argument("--input-structure-index", type=int, default=0)
    parser.add_argument("--model-checkpoint", type=Path, required=True)
    parser.add_argument(
        "--output-path",
        type=Path,
        default=repo_root / "build" / "torch_gnn_reference_fingerprint" / "torch_gnn_target_fingerprint.npy",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output_path = args.output_path.resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)

    model, checkpoint = load_model_from_checkpoint(args.model_checkpoint.resolve())
    structure = read_single_structure(
        args.input_structure.resolve(),
        index=args.input_structure_index,
    )
    fingerprint = np.asarray(model.compute_reference_fingerprint(structure), dtype=np.float32).reshape(-1)
    save_target_fingerprint(fingerprint, output_path)

    summary = {
        "input_structure": str(args.input_structure.resolve()),
        "input_structure_index": int(args.input_structure_index),
        "model_checkpoint": str(args.model_checkpoint.resolve()),
        "fingerprint_length": int(fingerprint.size),
        "fingerprint_min": float(np.min(fingerprint)) if fingerprint.size > 0 else 0.0,
        "fingerprint_max": float(np.max(fingerprint)) if fingerprint.size > 0 else 0.0,
        "output_path": str(output_path),
        "checkpoint_training_config": checkpoint.get("training_config", {}),
    }
    print(f"Saved analytical RAFFLE fingerprint: {output_path}")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
