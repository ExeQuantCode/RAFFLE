"""Validate Python graph topology against analytical RAFFLE descriptors.

This script compares graph-level statistics and descriptor surrogate agreement
on identical structures using the TorchGNN multigraph model and the analytical
Fortran RAFFLE descriptor backend.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from ase.io import read

from raffle import TorchGNNFingerprint


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--structures",
        type=Path,
        required=True,
        help="Path to an XYZ file containing one or more structures.",
    )
    parser.add_argument(
        "--max-structures",
        type=int,
        default=8,
        help="Maximum number of structures to validate.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("build") / "torch_gnn_graph_parity_validation.json",
        help="Output JSON report path.",
    )
    parser.add_argument("--seed", type=int, default=42)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if TorchGNNFingerprint is None:
        raise RuntimeError("TorchGNNFingerprint is unavailable in this installation")

    structures = read(str(args.structures), index=":")
    if not isinstance(structures, list):
        structures = [structures]
    selected = structures[: max(int(args.max_structures), 1)]

    if not selected:
        raise ValueError("No structures were loaded for validation")

    first_symbols = sorted(set(selected[0].get_chemical_symbols()))
    model = TorchGNNFingerprint(
        species_list=first_symbols,
        hidden_dim=64,
        num_message_layers=2,
        seed=int(args.seed),
    )

    per_structure = []
    for index, atoms in enumerate(selected):
        graph_stats = model.graph_statistics(atoms)
        agreement = model.descriptor_surrogate_agreement(atoms)
        per_structure.append(
            {
                "index": int(index),
                "formula": str(atoms.get_chemical_formula()),
                "graph_stats": graph_stats,
                "descriptor_agreement": agreement,
            }
        )

    param_counts = model.parameter_counts()
    summary = {
        "structure_count": len(per_structure),
        "mean_num_2body_edges": float(np.mean([entry["graph_stats"]["num_2body_edges"] for entry in per_structure])),
        "mean_num_3body_edges": float(np.mean([entry["graph_stats"]["num_3body_edges"] for entry in per_structure])),
        "mean_num_4body_edges": float(np.mean([entry["graph_stats"]["num_4body_edges"] for entry in per_structure])),
        "mean_rmse_2body": float(np.mean([entry["descriptor_agreement"]["rmse_2body"] for entry in per_structure])),
        "mean_rmse_3body": float(np.mean([entry["descriptor_agreement"]["rmse_3body"] for entry in per_structure])),
        "mean_rmse_4body": float(np.mean([entry["descriptor_agreement"]["rmse_4body"] for entry in per_structure])),
        "parameter_counts": param_counts,
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "summary": summary,
        "structures": per_structure,
    }
    args.output.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    print("Graph parity validation complete")
    print(json.dumps(summary, indent=2))
    print(f"Report: {args.output}")


if __name__ == "__main__":
    main()
