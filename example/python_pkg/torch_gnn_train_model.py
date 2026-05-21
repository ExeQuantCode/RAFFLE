"""Train and save a single PyTorch GNN surrogate model for inverse design.

This is the recommended split-workflow training entry point. It writes a reusable
checkpoint together with the resolved training config, target fingerprint, and
reference structure needed for reproducible inverse-design runs.
"""

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
    build_augmented_dataset,
    create_model,
    infer_species_list,
    load_structures,
    normalise_perturbation_settings,
    perturb_structure,
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
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=None,
        help="Optional JSON config file. Explicit CLI flags override config values.",
    )
    parser.add_argument(
        "--training-structures",
        type=str,
        default=None,
        help="Structure file used for surrogate training.",
    )
    parser.add_argument(
        "--training-structure-limit",
        type=int,
        default=None,
        help="Limit the number of loaded training structures. Use 0 for all structures.",
    )
    parser.add_argument(
        "--reference-structure",
        type=str,
        default=None,
        help="Structure whose analytical fingerprint is exported with the checkpoint bundle.",
    )
    parser.add_argument(
        "--reference-structure-index",
        type=int,
        default=None,
        help="Frame index when --reference-structure contains multiple structures.",
    )
    parser.add_argument(
        "--augmented-count",
        type=int,
        default=None,
        help="Number of noisy reference-structure augmentations to add during training.",
    )
    parser.add_argument(
        "--augmentation-noise-scale",
        type=float,
        default=None,
        help="Maximum uniform displacement magnitude in angstrom applied to each augmented structure.",
    )
    parser.add_argument(
        "--minimum-interatomic-distance",
        type=float,
        default=None,
        help="Reject or resample perturbations whose minimum pair distance falls below this threshold.",
    )
    parser.add_argument(
        "--augmentation-max-resamples",
        type=int,
        default=None,
        help="Maximum perturbation resampling attempts per structure.",
    )
    parser.add_argument("--epochs", type=int, default=None, help="Number of training epochs.")
    parser.add_argument("--batch-size", type=int, default=None, help="Mini-batch size.")
    parser.add_argument("--seed", type=int, default=None, help="Random seed for training and augmentation.")
    parser.add_argument(
        "--species-list",
        type=str,
        default=None,
        help="Optional comma-separated species list. By default it is inferred from the data.",
    )
    parser.add_argument("--hidden-dim", type=int, default=None, help="Hidden width of each message-passing branch.")
    parser.add_argument("--architecture", type=str, default=None, help="Torch GNN branch architecture alias.")
    parser.add_argument(
        "--num-message-layers",
        type=int,
        default=None,
        help="Number of message-passing layers in each branch.",
    )
    parser.add_argument("--learning-rate", type=float, default=None, help="Initial optimiser learning rate.")
    parser.add_argument(
        "--model-lr-decay-rate",
        type=float,
        default=None,
        help="Exponential learning-rate decay applied once per epoch.",
    )
    parser.add_argument(
        "--smooth-cutoff-width",
        type=float,
        default=None,
        help="Smooth cutoff width used in the graph construction features.",
    )
    parser.add_argument(
        "--reference-layer-type",
        type=int,
        default=None,
        help="Fortran reference model layer type forwarded into TorchGNNFingerprint.",
    )
    parser.add_argument("--component-weight-2body", type=float, default=None, help="Relative loss weight for the 2-body component.")
    parser.add_argument("--component-weight-3body", type=float, default=None, help="Relative loss weight for the 3-body component.")
    parser.add_argument("--component-weight-4body", type=float, default=None, help="Relative loss weight for the 4-body component.")
    parser.add_argument(
        "--output-dir",
        type=str,
        default=None,
        help="Directory for the checkpoint, target fingerprint, resolved config, and metrics.",
    )
    return parser


def _default_config(repo_root: Path) -> dict[str, Any]:
    return {
        "training_structures": str(repo_root / "example" / "data" / "carbon.xyz"),
        "training_structure_limit": 0,
        "reference_structure": None,
        "reference_structure_index": 0,
        "augmented_count": 2,
        "augmentation_noise_scale": 1.0,
        "minimum_interatomic_distance": 0.8,
        "augmentation_max_resamples": 64,
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
        "minimum_interatomic_distance": args.minimum_interatomic_distance,
        "augmentation_max_resamples": args.augmentation_max_resamples,
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

    perturbation_settings = normalise_perturbation_settings(
        {
            "min_displacement": 0.0,
            "max_displacement": float(config["augmentation_noise_scale"]),
            "minimum_interatomic_distance": float(config["minimum_interatomic_distance"]),
            "max_resamples": int(config["augmentation_max_resamples"]),
        }
    )

    if int(config["augmented_count"]) <= 0:
        raise ValueError("augmented_count must be positive so training uses perturbed dataset structures only")

    if config["reference_structure"] is None:
        reference_source_index = int(int(config["seed"]) % len(training_structures))
        reference_source = training_structures[reference_source_index]
        reference_source_path = str(Path(config["training_structures"]))
    else:
        reference_source_index = int(config["reference_structure_index"])
        reference_source = read_single_structure(
            Path(config["reference_structure"]),
            index=reference_source_index,
        )
        reference_source_path = str(Path(config["reference_structure"]))

    reference_structure = perturb_structure(
        reference_source,
        seed=int(config["seed"]),
        perturbation_settings=perturbation_settings,
    )

    model = create_model(
        seed=int(config["seed"]),
        species_list=species_list,
        model_config=config["model_config"],
    )
    augmentation_rng = np.random.default_rng(int(config["seed"]) + 1)
    preview_structures = build_augmented_dataset(
        training_structures,
        variants_per_structure=int(config["augmented_count"]),
        rng=augmentation_rng,
        perturbation_settings=perturbation_settings,
    )
    if not preview_structures:
        raise ValueError("No perturbed training structures were generated")

    requested_epochs = max(int(config["epochs"]), 0)
    if requested_epochs == 0:
        history = model.fit(
            preview_structures,
            num_epochs=0,
            batch_size=int(config["batch_size"]),
            verbose=0,
            reset_optimiser=True,
            recalibrate_base=True,
        )
    else:
        history: list[float] = []
        epoch_structures = preview_structures
        for epoch_index in range(requested_epochs):
            epoch_history = model.fit(
                epoch_structures,
                num_epochs=1,
                batch_size=int(config["batch_size"]),
                verbose=1,
                reset_optimiser=(epoch_index == 0),
                recalibrate_base=(epoch_index == 0),
            )
            if not history:
                history.extend(float(value) for value in epoch_history)
            else:
                history.append(float(epoch_history[-1]))
            if epoch_index + 1 < requested_epochs:
                epoch_structures = build_augmented_dataset(
                    training_structures,
                    variants_per_structure=int(config["augmented_count"]),
                    rng=augmentation_rng,
                    perturbation_settings=perturbation_settings,
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
        "reference_structure_source": {
            "path": reference_source_path,
            "index": int(reference_source_index),
        },
        "perturbation_settings": dict(perturbation_settings),
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
        "num_augmented_structures": len(preview_structures),
        "training_history": [float(value) for value in history],
        "initial_training_loss": float(history[0]),
        "final_training_loss": float(history[-1]),
        "species_list": species_list,
        "model_config": config["model_config"],
        "perturbation_settings": dict(perturbation_settings),
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
