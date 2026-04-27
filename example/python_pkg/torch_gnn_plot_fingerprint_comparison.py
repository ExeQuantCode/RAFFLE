"""Plot analytical RAFFLE and model-predicted fingerprints for a structure."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from torch_gnn_workflow_common import load_model_from_checkpoint, read_single_structure


def parse_args() -> argparse.Namespace:
    repo_root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--structure", type=Path, required=True)
    parser.add_argument("--structure-index", type=int, default=0)
    parser.add_argument("--model-checkpoint", type=Path, required=True)
    parser.add_argument(
        "--output-path",
        type=Path,
        default=repo_root / "build" / "torch_gnn_fingerprint_comparison.png",
    )
    return parser.parse_args()


def _component_ranges(model) -> list[tuple[str, int, int]]:
    start = 0
    ranges = []
    for label, width in (
        ("2-body", int(model.fingerprint_dim_2body)),
        ("3-body", int(model.fingerprint_dim_3body)),
        ("4-body", int(model.fingerprint_dim_4body)),
    ):
        end = start + width
        ranges.append((label, start, end))
        start = end
    return ranges


def _component_metrics(reference: np.ndarray, predicted: np.ndarray, model) -> dict[str, dict[str, float]]:
    metrics: dict[str, dict[str, float]] = {}
    for label, start, end in _component_ranges(model):
        reference_slice = reference[start:end]
        predicted_slice = predicted[start:end]
        metrics[label] = {
            "mse": float(np.mean((predicted_slice - reference_slice) ** 2)),
            "mae": float(np.mean(np.abs(predicted_slice - reference_slice))),
            "l2": float(np.linalg.norm(predicted_slice - reference_slice)),
        }
    return metrics


def plot_fingerprint_comparison(
    reference_fingerprint: np.ndarray,
    predicted_fingerprint: np.ndarray,
    model,
    output_path: Path,
    title: str,
) -> None:
    indices = np.arange(reference_fingerprint.size)
    difference = predicted_fingerprint - reference_fingerprint
    total_mse = float(np.mean(difference ** 2))
    total_mae = float(np.mean(np.abs(difference)))

    figure = plt.figure(figsize=(14, 8))
    ax1 = figure.add_subplot(2, 1, 1)
    ax1.plot(indices, reference_fingerprint, label="RAFFLE analytical fingerprint", linewidth=2.0)
    ax1.plot(indices, predicted_fingerprint, label="Model predicted fingerprint", linewidth=1.5)
    ax1.set_ylabel("Fingerprint value")
    ax1.set_title(title)
    ax1.grid(alpha=0.3)
    ax1.legend(loc="best")

    ax2 = figure.add_subplot(2, 1, 2)
    ax2.plot(indices, difference, color="tab:red", linewidth=1.5)
    ax2.axhline(0.0, color="black", linewidth=1.0, linestyle="--")
    ax2.set_xlabel("Fingerprint index")
    ax2.set_ylabel("Predicted - analytical")
    ax2.grid(alpha=0.3)

    component_ranges = _component_ranges(model)
    for _, _, end in component_ranges[:-1]:
        ax1.axvline(end - 0.5, color="0.4", linestyle=":", linewidth=1.0)
        ax2.axvline(end - 0.5, color="0.4", linestyle=":", linewidth=1.0)
    for label, start, end in component_ranges:
        midpoint = 0.5 * (start + end - 1)
        ax1.text(
            midpoint,
            1.01,
            label,
            ha="center",
            va="bottom",
            transform=ax1.get_xaxis_transform(),
        )

    figure.suptitle(f"Fingerprint comparison: total MSE={total_mse:.6e}, total MAE={total_mae:.6e}")
    figure.tight_layout(rect=(0.0, 0.0, 1.0, 0.96))
    figure.savefig(output_path, dpi=150)
    plt.close(figure)


def main() -> None:
    args = parse_args()
    structure = read_single_structure(args.structure.resolve(), index=args.structure_index)
    model, checkpoint = load_model_from_checkpoint(args.model_checkpoint.resolve())

    reference_fingerprint = model.compute_reference_fingerprint(structure)
    predicted_fingerprint = model.predict(structure)
    if reference_fingerprint.shape != predicted_fingerprint.shape:
        raise ValueError(
            "Analytical and predicted fingerprints must have the same shape: "
            f"{reference_fingerprint.shape} != {predicted_fingerprint.shape}"
        )

    output_path = args.output_path.resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plot_fingerprint_comparison(
        reference_fingerprint=reference_fingerprint,
        predicted_fingerprint=predicted_fingerprint,
        model=model,
        output_path=output_path,
        title=f"Structure: {args.structure.name} (index {args.structure_index})",
    )

    metrics = {
        "structure_file": str(args.structure.resolve()),
        "structure_index": int(args.structure_index),
        "model_checkpoint": str(args.model_checkpoint.resolve()),
        "output_path": str(output_path),
        "fingerprint_length": int(reference_fingerprint.size),
        "total_mse": float(np.mean((predicted_fingerprint - reference_fingerprint) ** 2)),
        "total_mae": float(np.mean(np.abs(predicted_fingerprint - reference_fingerprint))),
        "component_metrics": _component_metrics(reference_fingerprint, predicted_fingerprint, model),
        "checkpoint_training_config": checkpoint.get("training_config", {}),
    }

    print(f"Saved fingerprint comparison plot: {output_path}")
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
