# %%
import os
import math
import numpy as np
import torch
import wandb
from pathlib import Path
from matplotlib import pyplot as plt
from ase.io import read, write
from ase.build import bulk
from raffle import TorchGNNFingerprint, structure_similarity_rmsd

# %%
# Load data
_script_dir = os.path.dirname(os.path.abspath(__file__))
_database_path = os.path.join(_script_dir, "..", "data", "carbon.xyz")
database = read(_database_path, index=":")

# Add silicon to the database
silicon = bulk("Si", "diamond", a=5.43)
database.append(silicon)

species_list = []
for atoms in database:
    species_list.extend(atoms.symbols.species())
species_list = sorted(set(species_list))
print(f"Species: {species_list}")

# %%
# Define sweep configuration
sweep_config = {
    "method": "bayes",  # or "grid", "random"
    "metric": {
        "name": "best_position_difference",
        "goal": "minimize"
    },
    "parameters": {
        "model.component_weights": {
            "values": [[1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 0], [1, 0, 1], [0, 1, 1], [2, 1, 1], [4, 1, 1], [1, 2, 1], [1, 4, 1], [1, 1, 2], [1, 1, 4]]
        },
        "model.hidden_dim_2body": {
            "values": [8, 16, 32, 64]
        },
        "model.hidden_dim_3body": {
            "values": [16, 32, 64, 128]
        },
        "model.hidden_dim_4body": {
            "values": [16, 32, 64, 128]
        },
        "model.num_message_layers_2body": {
            "values": [1, 2, 3]
        },
        "model.num_message_layers_3body": {
            "values": [1, 2, 3]
        },
        "model.num_message_layers_4body": {
            "values": [1, 2, 3]
        },
        "train.batch_size": {
            "values": [8, 16, 32]
        },
        "train.num_epochs": {
            "values": [50, 100, 200]
        },
        "train.learning_rate": {
            "values": [1e-3, 5e-3, 1e-2, 5e-2]
        },
        "train.lr_decay_rate": {
            "values": [0.001, 0.005, 0.01]
        },
        "seed": {
            "min": 1,
            "max": 1000
        }
    }
}

# %%
def plot_position_history(
    traces: dict,
    step_size_results: list[dict],
    initial_position_diff: float,
    output_path: Path,
) -> None:
    """
    Plot inverse design position difference traces and step size sweep summary.

    Args:
        traces: Dictionary mapping step_size -> list of trace records
        step_size_results: List of dicts with step_size, final_position_difference,
                          and best_position_difference
        initial_position_diff: Initial position difference
        output_path: Path to save the plot
    """
    figure = plt.figure(figsize=(10.5, 5))

    # Left subplot: Position difference traces for all step sizes
    ax1 = figure.add_subplot(1, 2, 1)
    colors = plt.cm.viridis(np.linspace(0, 1, len(traces)))

    for idx, (step_size_key, trace) in enumerate(sorted(traces.items())):
        steps = [entry["step"] for entry in trace if not entry["is_initial_state"]]
        pos_diffs = [entry["position_difference"] for entry in trace if not entry["is_initial_state"]]

        # Convert step_size_key to float for display
        step_size_val = float(step_size_key)
        ax1.plot(steps, pos_diffs, marker='o', linewidth=2, markersize=4,
                color=colors[idx], label=f'Step size = {step_size_val:.3f}')

    ax1.axhline(initial_position_diff, color='tab:red', linestyle='--', linewidth=1.0)
    ax1.set_xlabel("Inverse-design steps")
    ax1.set_ylabel("Position difference (RMSD, Å)")
    ax1.set_title("Inverse Design Structure Evolution", fontsize=14)
    ax1.grid(alpha=0.3)
    ax1.legend(loc='upper right', framealpha=0.7)

    # Right subplot: Step size sweep summary
    ax2 = figure.add_subplot(1, 2, 2)
    sorted_data = sorted(step_size_results, key=lambda x: x["step_size"])
    step_sizes = [d["step_size"] for d in sorted_data]
    final_pos_diffs = [d["final_position_difference"] for d in sorted_data]
    best_pos_diffs = [d["best_position_difference"] for d in sorted_data]

    ax2.plot(step_sizes, final_pos_diffs, marker='o', linewidth=2, linestyle='-.', markersize=6,
             label='Final position difference', color='tab:blue')
    ax2.plot(step_sizes, best_pos_diffs, marker='s', linewidth=2, linestyle=':', markersize=6,
             label='Best position difference', color='tab:green')
    ax2.axhline(initial_position_diff, color='tab:red', linestyle='--', linewidth=1.0,
                label=f'Initial: {initial_position_diff:.4f} Å')

    ax2.set_xscale('log')
    ax2.set_xlabel("Step size")
    ax2.set_ylabel("Position difference (RMSD, Å)")
    ax2.set_title("Inverse Design Structure Evolution for Different Step Sizes", fontsize=14)
    ax2.grid(alpha=0.3)
    ax2.legend()

    figure.tight_layout()
    figure.savefig(output_path, dpi=150)
    plt.close(figure)


def plot_fingerprint_comparison(target_fingerprint: np.ndarray, predicted_fingerprint: np.ndarray,
                                save_path: Path, title_suffix: str = ""):
    """Plot reference vs predicted fingerprint."""
    plt.figure(figsize=(10, 6))
    plt.plot(target_fingerprint, label="Target (Diamond Reference)", linewidth=2, alpha=0.8)
    plt.plot(predicted_fingerprint, label="Predicted by Model", linewidth=2, alpha=0.8)
    plt.xlabel("Fingerprint Index", fontsize=12)
    plt.ylabel("Fingerprint Value", fontsize=12)
    plt.title(f"Diamond Reference Fingerprint Comparison{title_suffix}", fontsize=14)
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(save_path, dpi=150)
    plt.close()


def get_atoms_structure(structure):
    """Helper function to extract Atoms object from various return types."""
    if hasattr(structure, 'get_positions'):  # Is already an Atoms object
        return structure
    elif isinstance(structure, (list, tuple)):
        # If it's a list, get the last element or the first Atoms object
        for item in reversed(structure):
            if hasattr(item, 'get_positions'):
                return item
        # If no Atoms found, return the first item
        return structure[0] if structure else None
    else:
        return structure


def run_inverse_design_for_step_size(model, target_fingerprint, perturbed_diamond,
                                     diamond, step_size, inverse_config):
    """Run inverse design for a specific step size and return results."""
    # Create copy of config with specific step size
    config = inverse_config.copy()
    config['step_size'] = step_size
    config['return_trajectory'] = True

    atoms_traj = model.inverse_design(
        target_fingerprint,
        perturbed_diamond,
        **config
    )

    # Handle different return types

    # Ensure atoms_relaxed is an Atoms object
    final_atoms = atoms_traj[-1]

    # Calculate position differences
    initial_pos_diff = structure_similarity_rmsd(diamond, perturbed_diamond)
    final_pos_diff = structure_similarity_rmsd(diamond, final_atoms)

    # Track best position difference from trajectory
    best_pos_diff = final_pos_diff
    trace_data = []

    for i, atoms in enumerate(atoms_traj):

        pos_diff = structure_similarity_rmsd(diamond, atoms)
        best_pos_diff = min(best_pos_diff, pos_diff)
        trace_data.append({
            "step": i,
            "is_initial_state": i == 0,
            "position_difference": pos_diff
        })

    return {
        "step_size": step_size,
        "atoms": final_atoms,
        "trajectory": atoms_traj,
        "trace_data": trace_data,
        "initial_position_difference": initial_pos_diff,
        "final_position_difference": final_pos_diff,
        "best_position_difference": best_pos_diff
    }


def train_and_evaluate(use_wandb=True, run_name=None, config=None):
    """Train a model and perform inverse design for multiple step sizes."""
    # Setup wandb if requested
    if use_wandb:
        wandb.init()
        config = wandb.config

    # If we have a config dict directly, use it
    if not use_wandb and config is None:
        # Use default sweep config values for testing
        config = {}
        for key, param in sweep_config["parameters"].items():
            if "values" in param:
                config[key] = param["values"][0]
            elif "value" in param:
                config[key] = param["value"]
            elif "min" in param and "max" in param:
                config[key] = (param["min"] + param["max"]) // 2

    # Extract configs
    if config is not None:
        # From wandb
        model_config = {k.replace('model.', ''): v for k, v in config.items() if k.startswith('model.')}
        train_config = {k.replace('train.', ''): v for k, v in config.items() if k.startswith('train.')}
        component_weights = model_config.pop('component_weights', [2, 2, 2])
        seed = config.get('seed', 42)
        num_epochs = train_config.get("num_epochs", 100)
        batch_size = train_config.get("batch_size", 16)
        learning_rate = train_config.get("learning_rate", 1e-2)
        lr_decay_rate = train_config.get("lr_decay_rate", 0.005)

    # Create model with seed from config
    model = TorchGNNFingerprint(
        species_list=species_list,
        component_weight=component_weights,
        seed=seed,
        **model_config
    )

    # Log parameter count and component weights
    num_params = sum(p.numel() for p in model.parameters() if p.requires_grad)
    if use_wandb:
        wandb.log({
            "parameter_count": num_params,
            "component_weight_2body": component_weights[0],
            "component_weight_3body": component_weights[1],
            "component_weight_4body": component_weights[2],
            "seed": seed
        })
    else:
        print(f"Parameter count: {num_params}")
        print(f"Component weights: {component_weights}")
        print(f"Seed: {seed}")

    # Train model using the fit method with wandb integration
    print(f"Training model with config: {model_config}")
    print(f"Component weights: {component_weights}")

    # Train with wandb logging only if requested
    history = model.fit(
        database,
        num_epochs=num_epochs,
        batch_size=batch_size,
        learning_rate=learning_rate,
        lr_decay_rate=lr_decay_rate,
        verbose=1,
        use_wandb=use_wandb,
        wandb_project="raffle_inverse_design_sweep" if use_wandb else None,
        wandb_run_name=run_name if use_wandb else None,
        wandb_config=model_config if use_wandb else None
    )

    # Log final training loss
    if use_wandb:
        wandb.log({"final_training_loss": history[-1]})

    # Save model
    model_path = Path(f"model_{run_name if run_name else 'test'}.pth")
    save_dict = {
        "state_dict": model.state_dict(),
        "config": {"species_list": species_list, **model_config}
    }
    torch.save(save_dict, model_path)
    if use_wandb:
        wandb.save(str(model_path))
    else:
        print(f"Model saved to {model_path}")

    # Load diamond structure for inverse design
    _diamond_path = os.path.join(_script_dir, "diamond.xyz")
    diamond = read(_diamond_path)

    # Create perturbed structure (replace one C with Si)
    _perturbed_path = os.path.join(_script_dir, "perturbed_diamond.xyz")
    perturbed_diamond = read(_perturbed_path)
    for atom in perturbed_diamond:
        if atom.symbol == "C":
            atom.symbol = "Si"
            break

    # Compute target fingerprint
    target_fingerprint = model._compute_reference_fingerprint(diamond)

    # Set up inverse design parameters
    inverse_config = {
        'update_topology_every_n_steps': 1,
        'optimize_species': True,
        'species_learning_rate': 0.5,
        'verbose': 1,
        'num_restarts': 1,
        'use_augmented_lagrangian': False,
        'num_steps': 100,
    }

    # Define step sizes to test (from 1e-2 to 1e0)
    step_sizes = [0.01, 0.05, 0.1, 0.5, 1.0]

    # Run inverse design for each step size
    print(f"Running inverse design for {len(step_sizes)} step sizes...")
    step_results = {}

    for step_size in step_sizes:
        print(f"  Step size: {step_size}")
        result = run_inverse_design_for_step_size(
            model, target_fingerprint, perturbed_diamond,
            diamond, step_size, inverse_config
        )
        step_results[str(step_size)] = result

        # Log individual step size results
        if result["atoms"] is not None and use_wandb:
            wandb.log({
                f"step_size_{step_size:.3f}_final_pos_diff": result["final_position_difference"],
                f"step_size_{step_size:.3f}_best_pos_diff": result["best_position_difference"],
            })

    # Find best overall result
    valid_results = [r for key, r in step_results.items() if r["atoms"] is not None]
    if valid_results:
        best_result = min(valid_results, key=lambda x: x["best_position_difference"])
        initial_position_diff = step_results["0.01"]["initial_position_difference"]
    else:
        # If all failed, return a large number
        if use_wandb:
            wandb.finish()
        return 999.0

    # Log best metrics
    if use_wandb:
        wandb.log({
            "initial_position_diff": initial_position_diff,
            "best_final_position_difference": best_result["final_position_difference"],
            "best_position_difference": best_result["best_position_difference"],
            "best_step_size": best_result["step_size"],
            "position_difference_improvement": initial_position_diff - best_result["best_position_difference"]
        })
    else:
        print(f"\nResults:")
        print(f"  Initial position difference: {initial_position_diff:.4f} Å")
        print(f"  Best final position difference: {best_result['final_position_difference']:.4f} Å")
        print(f"  Best position difference: {best_result['best_position_difference']:.4f} Å")
        print(f"  Best step size: {best_result['step_size']:.3f}")
        print(f"  Improvement: {initial_position_diff - best_result['best_position_difference']:.4f} Å")

    # Create and save plots
    run_id = run_name if run_name else "test"
    temp_dir = Path(f"temp_{run_id}")
    temp_dir.mkdir(exist_ok=True)

    # Plot 1: Position difference traces for all step sizes
    traces_dict = {}
    traces_dict["0.01"] = step_results["0.01"]["trace_data"]

    # Plot 2: Step size sweep summary
    step_size_data = []
    for key, result in step_results.items():
        if result["atoms"] is not None:
            step_size_data.append({
                "step_size": result["step_size"],
                "final_position_difference": result["final_position_difference"],
                "best_position_difference": result["best_position_difference"]
            })

    # Replace the two separate plotting calls with:
    if traces_dict and step_size_data:
        combined_plot_path = temp_dir / "position_history.png"
        plot_position_history(
            traces=traces_dict,
            step_size_results=step_size_data,
            initial_position_diff=initial_position_diff,
            output_path=combined_plot_path
        )
        if use_wandb:
            wandb.log({"position_history_plot": wandb.Image(str(combined_plot_path))})
        else:
            print(f"Position difference plot saved to {combined_plot_path}")

    # Plot 3: Fingerprint comparison for best result
    if best_result["atoms"] is not None:
        predicted_fingerprint = model.compute_fingerprint(best_result["atoms"])
        fingerprint_plot_path = temp_dir / "fingerprint_comparison.png"
        plot_fingerprint_comparison(
            target_fingerprint,
            predicted_fingerprint,
            fingerprint_plot_path,
            f" (best step size = {best_result['step_size']:.3f})"
        )
        if use_wandb:
            wandb.log({"fingerprint_comparison_plot": wandb.Image(str(fingerprint_plot_path))})
        else:
            print(f"Plot saved to {fingerprint_plot_path}")

    # Also log the best structure
    if best_result["atoms"] is not None:
        best_structure_path = temp_dir / "best_structure.xyz"
        write(best_structure_path, best_result["atoms"])
        if use_wandb:
            wandb.save(str(best_structure_path))
        else:
            print(f"Best structure saved to {best_structure_path}")

    if not use_wandb:
        print(f"Saving optimisation trajectory for best step size to {temp_dir / 'best_trajectory.xyz'}")
        write(temp_dir / "best_trajectory.xyz", best_result["trajectory"])

    # Log summary table of step size results
    if step_size_data and use_wandb:
        step_table = wandb.Table(
            columns=["step_size", "final_position_difference", "best_position_difference"]
        )
        for data in step_size_data:
            step_table.add_data(
                data["step_size"],
                data["final_position_difference"],
                data["best_position_difference"]
            )
        wandb.log({"step_size_summary": step_table})
    elif not use_wandb and step_size_data:
        print("\nStep size summary:")
        print("Step Size | Final Pos Diff | Best Pos Diff")
        print("-" * 45)
        for data in step_size_data:
            print(f"{data['step_size']:8.3f} | {data['final_position_difference']:14.4f} | {data['best_position_difference']:14.4f}")

    # Clean up
    if use_wandb:
      import shutil
      shutil.rmtree(temp_dir, ignore_errors=True)
    if model_path.exists() and not use_wandb:
        # Keep the model for test runs
        pass
    elif model_path.exists() and use_wandb:
        os.remove(model_path)

    # Finish wandb run
    if use_wandb:
        wandb.finish()

    return best_result["best_position_difference"]


def run_sweep(sweep_id=None, count=200, project="raffle_inverse_design_sweep"):
    """Run a wandb sweep (either new or continuing)."""
    if sweep_id is None:
        # Start a new sweep
        sweep_id = wandb.sweep(sweep_config, project=project)
        print(f"Created new sweep with ID: {sweep_id}")

    # Initialize the agent
    wandb.agent(sweep_id, function=lambda: train_and_evaluate(use_wandb=True), count=count)
    return sweep_id


def run_single_test(config_override=None):
    """Run a single test run without wandb."""
    # Build config from sweep defaults
    config_dict = {}
    for key, param in sweep_config["parameters"].items():
        if "values" in param:
            config_dict[key] = param["values"][0]
        elif "value" in param:
            config_dict[key] = param["value"]
        elif "min" in param and "max" in param:
            config_dict[key] = (param["min"] + param["max"]) // 2

    # Apply any overrides
    if config_override:
        for key, value in config_override.items():
            if key in config_dict:
                config_dict[key] = value
            else:
                print(f"Warning: Unknown parameter '{key}' in config override")

    print("Running single test with parameters:")
    for key, value in config_dict.items():
        print(f"  {key}: {value}")

    # Run training without wandb
    return train_and_evaluate(use_wandb=False, run_name="test", config=config_dict)


# %%
# Run sweep
if __name__ == "__main__":
    import argparse
    import sys

    parser = argparse.ArgumentParser(description='Run hyperparameter optimization')
    parser.add_argument('--mode', type=str, default='sweep',
                        choices=['sweep', 'test', 'continue'],
                        help='Run mode: sweep (new sweep), test (single test run), or continue (existing sweep)')
    parser.add_argument('--sweep_id', type=str, default=None,
                        help='Sweep ID to continue (required for continue mode)')
    parser.add_argument('--count', type=int, default=None,
                        help='Number of sweep runs to execute')
    parser.add_argument('--project', type=str, default='raffle_inverse_design_sweep',
                        help='Wandb project name')
    parser.add_argument('--config', type=str, default=None,
                        help='Path to JSON config file with hyperparameters (for test mode)')

    args = parser.parse_args()

    # Set default counts if not specified
    if args.count is None:
        if args.mode == 'sweep':
            args.count = 200
        elif args.mode == 'continue':
            args.count = 10
        else:
            args.count = 1

    if args.mode == 'sweep':
        print("Starting new wandb sweep...")
        run_sweep(sweep_id=None, count=args.count, project=args.project)
    elif args.mode == 'continue':
        if args.sweep_id is None:
            print("Error: --sweep_id is required for continue mode")
            sys.exit(1)
        print(f"Continuing wandb sweep: {args.sweep_id}")
        run_sweep(sweep_id=args.sweep_id, count=args.count, project=args.project)
    else:  # test mode
        if args.config:
            print(f"Running test with config file: {args.config}")
            import json
            with open(args.config, 'r') as f:
                config_override = json.load(f)
            run_single_test(config_override)
        else:
            print("Running single test with default parameters...")
            run_single_test(None)
