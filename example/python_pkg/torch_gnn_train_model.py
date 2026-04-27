"""Train and save a single PyTorch GNN surrogate model for inverse design."""

from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path
import sys
from typing import Any


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from torch_gnn_workflow_common import (
    DEFAULT_COMPONENT_WEIGHT,
    DEFAULT_MODEL_CONFIG,
    build_augmented_structures,
    create_model,
    default_reference_structure,
    infer_species_list,
    load_structures,
    normalise_model_config,
    read_json_config,
    read_single_structure,
    resolve_path,
    save_model_checkpoint,
    save_target_fingerprint,
    write_json,
    write_structure,
)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=None)
    parser.add_argument("--training-structures", type=str, default=None)
    parser.add_argument("--training-structure-limit", type=int, default=None)
    parser.add_argument("--reference-structure", type=str, default=None)
    parser.add_argument("--reference-structure-index", type=int, default=None)
    parser.add_argument("--augmented-count", type=int, default=None)
    parser.add_argument("--augmentation-noise-scale", type=float, default=None)
    parser.add_argument("--epochs", type=int, default=None)
    parser.add_argument("--batch-size", type=int, default=None)
    parser.add_argument("--seed", type=int, default=None)
    parser.add_argument("--species-list", type=str, default=None)
    parser.add_argument("--hidden-dim", type=int, default=None)
    parser.add_argument("--architecture", type=str, default=None)
    parser.add_argument("--num-message-layers", type=int, default=None)
    parser.add_argument("--learning-rate", type=float, default=None)
    parser.add_argument("--model-lr-decay-rate", type=float, default=None)
    parser.add_argument("--smooth-cutoff-width", type=float, default=None)
    parser.add_argument("--reference-layer-type", type=int, default=None)
    parser.add_argument("--component-weight-2body", type=float, default=None)
    parser.add_argument("--component-weight-3body", type=float, default=None)
    parser.add_argument("--component-weight-4body", type=float, default=None)
    parser.add_argument("--output-dir", type=str, default=None)
    return parser


def _default_config(repo_root: Path) -> dict[str, Any]:
    return {
        "training_structures": str(repo_root / "example" / "data" / "carbon.xyz"),
        "training_structure_limit": 0,
        "reference_structure": None,
        "reference_structure_index": 0,
        "augmented_count": 16,
        "augmentation_noise_scale": 0.04,
        "epochs": 30,
        "batch_size": 8,
        "seed": 42,
        "species_list": None,
        "output_dir": str(repo_root / "build" / "torch_gnn_train_model"),
        "model_config": copy.deepcopy(DEFAULT_MODEL_CONFIG),
    }


def _merge_config(resolved: dict[str, Any], incoming: dict[str, Any]) -> None:
    for key, value in incoming.items():
        if value is None:
            continue
        if key == "model_config":
            if not isinstance(value, dict):
                raise ValueError("model_config must be a JSON object")
            resolved["model_config"].update(
                {nested_key: nested_value for nested_key, nested_value in value.items() if nested_value is not None}
            )
            continue
        resolved[key] = value


def _parse_species_list(value: Any) -> list[str]:
    if value is None:
        return []
    if isinstance(value, str):
        return [item.strip() for item in value.split(",") if item.strip()]
    if isinstance(value, (list, tuple)):
        return [str(item).strip() for item in value if str(item).strip()]
    raise ValueError("species_list must be either a comma-delimited string or a JSON array")


def parse_args() -> dict[str, Any]:
    args = _parser().parse_args()
    repo_root = Path(__file__).resolve().parents[2]
    resolved = _default_config(repo_root)
    config_base_dir = Path.cwd()

    if args.config is not None:
        config_path = args.config.resolve()
        config_base_dir = config_path.parent
        _merge_config(resolved, read_json_config(config_path))

    cli_overrides = {
        "training_structures": args.training_structures,
        "training_structure_limit": args.training_structure_limit,
        "reference_structure": args.reference_structure,
        "reference_structure_index": args.reference_structure_index,
        "augmented_count": args.augmented_count,
        "augmentation_noise_scale": args.augmentation_noise_scale,
        "epochs": args.epochs,
        "batch_size": args.batch_size,
        "seed": args.seed,
        "species_list": args.species_list,
        "output_dir": args.output_dir,
        "model_config": {
            "architecture": args.architecture,
            "hidden_dim": args.hidden_dim,
            "num_message_layers": args.num_message_layers,
            "learning_rate": args.learning_rate,
            "lr_decay_rate": args.model_lr_decay_rate,
            "smooth_cutoff_width": args.smooth_cutoff_width,
            "reference_layer_type": args.reference_layer_type,
        },
    }
    if (
        args.component_weight_2body is not None
        or args.component_weight_3body is not None
        or args.component_weight_4body is not None
    ):
        cli_overrides["model_config"]["component_weight"] = [
            DEFAULT_COMPONENT_WEIGHT[0] if args.component_weight_2body is None else args.component_weight_2body,
            DEFAULT_COMPONENT_WEIGHT[1] if args.component_weight_3body is None else args.component_weight_3body,
            DEFAULT_COMPONENT_WEIGHT[2] if args.component_weight_4body is None else args.component_weight_4body,
        ]
    _merge_config(resolved, cli_overrides)

    resolved["training_structures"] = str(
        resolve_path(resolved["training_structures"], config_base_dir)
    )
    resolved["output_dir"] = str(resolve_path(resolved["output_dir"], config_base_dir))
    resolved["reference_structure"] = (
        None
        if resolved["reference_structure"] is None
        else str(resolve_path(resolved["reference_structure"], config_base_dir))
    )
    resolved["model_config"] = {
        **normalise_model_config(resolved["model_config"]),
        "component_weight": list(
            normalise_model_config(resolved["model_config"])["component_weight"]
        ),
    }
    return resolved


def main() -> None:
    config = parse_args()
    repo_root = Path(__file__).resolve().parents[2]
    output_dir = Path(config["output_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)

    training_structures = load_structures(
        Path(config["training_structures"]),
        requested_count=int(config["training_structure_limit"]),
    )
    if not training_structures:
        raise ValueError("No training structures were loaded")

    species_list = _parse_species_list(config["species_list"])
    if not species_list:
        species_list = infer_species_list(training_structures)

    if config["reference_structure"] is None:
        if species_list != ["C"]:
            raise ValueError(
                "A reference structure file is required when species_list is not exactly ['C']"
            )
        reference_structure = default_reference_structure()
    else:
        reference_structure = read_single_structure(
            Path(config["reference_structure"]),
            index=int(config["reference_structure_index"]),
        )

    augmented_structures = build_augmented_structures(
        reference_structure,
        count=int(config["augmented_count"]),
        seed=int(config["seed"]),
        noise_scale=float(config["augmentation_noise_scale"]),
    )

    model = create_model(
        seed=int(config["seed"]),
        species_list=species_list,
        model_config=config["model_config"],
    )
    history = model.fit(
        training_structures,
        num_epochs=int(config["epochs"]),
        batch_size=int(config["batch_size"]),
        augment_structures=augmented_structures,
        verbose=1,
    )
    target_fingerprint = model.compute_reference_fingerprint(reference_structure)

    config_path = output_dir / "torch_gnn_train_config.json"
    checkpoint_path = output_dir / "torch_gnn_model_checkpoint.pt"
    fingerprint_path = output_dir / "torch_gnn_target_fingerprint.npy"
    reference_structure_path = output_dir / "torch_gnn_reference_structure.xyz"
    metrics_path = output_dir / "torch_gnn_training_metrics.json"

    effective_config = {
        **config,
        "repo_root": str(repo_root),
        "species_list": species_list,
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
        model=model,
        checkpoint_path=checkpoint_path,
        species_list=species_list,
        model_config=config["model_config"],
        training_config=effective_config,
        training_history=history,
    )
    save_target_fingerprint(target_fingerprint, fingerprint_path)
    write_structure(reference_structure_path, reference_structure)

    metrics = {
        "num_available_training_structures": len(load_structures(Path(config["training_structures"]))),
        "num_training_structures": len(training_structures),
        "num_augmented_structures": len(augmented_structures),
        "training_history": [float(value) for value in history],
        "initial_training_loss": float(history[0]),
        "final_training_loss": float(history[-1]),
        "species_list": species_list,
        "model_config": config["model_config"],
        "training_config": {
            key: value
            for key, value in effective_config.items()
            if key not in {"model_config", "output_files"}
        },
        "output_files": effective_config["output_files"],
    }
    write_json(metrics_path, metrics)

    print(f"Saved checkpoint: {checkpoint_path}")
    print(f"Saved target fingerprint: {fingerprint_path}")
    print(f"Saved config: {config_path}")
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
