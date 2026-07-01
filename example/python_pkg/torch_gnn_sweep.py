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
def plot_position_difference_multi_step(traces: dict, initial_position_diff: float, save_path: Path):
    """
    Plot inverse design position difference over steps for multiple step sizes.

    Args:
        traces: Dictionary mapping step_size -> list of trace records
        initial_position_diff: Initial position difference
        save_path: Path to save the plot
    """
    plt.figure(figsize=(10, 7))

    # Use a color map for different step sizes
    colors = plt.cm.viridis(np.linspace(0, 1, len(traces)))

    for idx, (step_size, trace) in enumerate(sorted(traces.items())):
        steps = [entry["step"] for entry in trace if not entry["is_initial_state"]]
        pos_diffs = [entry["position_difference"] for entry in trace if not entry["is_initial_state"]]

        plt.plot(steps, pos_diffs, marker='o', linewidth=2, markersize=6,
                 color=colors[idx], label=f'Step size = {step_size:.3f}')

    plt.axhline(initial_position_diff, color='red', linestyle='--', linewidth=1.5,
                label=f'Initial: {initial_position_diff:.4f} Å')
    plt.xlabel("Inverse Design Steps", fontsize=12)
    plt.ylabel("Position Difference (RMSD, Å)", fontsize=12)
    plt.title("Inverse Design Structure Evolution for Different Step Sizes", fontsize=14)
    plt.grid(alpha=0.3)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(save_path, dpi=150)
    plt.close()


def plot_step_size_summary(step_size_data: list[dict], initial_position_diff: float,
                           save_path: Path):
    """Plot summary of final position difference vs step size."""
    plt.figure(figsize=(8, 6))

    # Sort by step size
    sorted_data = sorted(step_size_data, key=lambda x: x["step_size"])
    step_sizes = [d["step_size"] for d in sorted_data]
    final_pos_diffs = [d["final_position_difference"] for d in sorted_data]
    best_pos_diffs = [d["best_position_difference"] for d in sorted_data]

    plt.plot(step_sizes, final_pos_diffs, marker='o', linewidth=2, markersize=8,
             label='Final position difference', color='blue')
    plt.plot(step_sizes, best_pos_diffs, marker='s', linewidth=2, markersize=8,
             label='Best position difference', color='green')
    plt.axhline(initial_position_diff, color='red', linestyle='--', linewidth=1.5,
                label=f'Initial: {initial_position_diff:.4f} Å')

    plt.xscale('log')
    plt.xlabel("Step Size (log scale)", fontsize=12)
    plt.ylabel("Position Difference (RMSD, Å)", fontsize=12)
    plt.title("Step Size Sweep: Final Structure Quality", fontsize=14)
    plt.grid(alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(save_path, dpi=150)
    plt.close()


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

    result = model.inverse_design(
        target_fingerprint,
        perturbed_diamond,
        **config
    )

    # Handle different return types
    if isinstance(result, tuple) and len(result) == 2:
        atoms_relaxed, trajectory = result
    else:
        atoms_relaxed = result
        trajectory = None

    # Ensure atoms_relaxed is an Atoms object
    if atoms_relaxed is not None:
        atoms_relaxed = get_atoms_structure(atoms_relaxed)

    # Calculate position differences
    initial_pos_diff = structure_similarity_rmsd(diamond, perturbed_diamond)

    if atoms_relaxed is not None:
        final_pos_diff = structure_similarity_rmsd(diamond, atoms_relaxed)
    else:
        final_pos_diff = initial_pos_diff

    # Track best position difference from trajectory
    best_pos_diff = final_pos_diff
    trace_data = []

    if trajectory is not None:
        for record in trajectory:
            candidate_atoms = record["atoms"]
            # Ensure candidate_atoms is an Atoms object
            if hasattr(candidate_atoms, 'copy'):
                candidate_atoms = candidate_atoms.copy()
            else:
                candidate_atoms = get_atoms_structure(candidate_atoms)

            if candidate_atoms is not None:
                pos_diff = structure_similarity_rmsd(diamond, candidate_atoms)
                best_pos_diff = min(best_pos_diff, pos_diff)
                trace_data.append({
                    "step": record["step"],
                    "is_initial_state": record["step"] == 0,
                    "position_difference": pos_diff
                })

    return {
        "step_size": step_size,
        "atoms": atoms_relaxed,
        "trajectory": trajectory,
        "trace_data": trace_data,
        "initial_position_difference": initial_pos_diff,
        "final_position_difference": final_pos_diff,
        "best_position_difference": best_pos_diff
    }


def train_and_evaluate():
    """Train a model and perform inverse design for multiple step sizes."""
    # Setup wandb FIRST
    wandb.init()
    config = wandb.config

    # Extract configs
    model_config = {k.replace('model.', ''): v for k, v in config.items() if k.startswith('model.')}
    train_config = {k.replace('train.', ''): v for k, v in config.items() if k.startswith('train.')}

    # Get component weights from config or use default
    component_weights = model_config.pop('component_weights', [2, 2, 2])

    # Create model with seed from config
    seed = config.get('seed', 42)

    model = TorchGNNFingerprint(
        species_list=species_list,
        component_weight=component_weights,
        seed=seed,
        **model_config
    )

    # Log parameter count and component weights
    num_params = sum(p.numel() for p in model.parameters() if p.requires_grad)
    wandb.log({
        "parameter_count": num_params,
        "component_weight_2body": component_weights[0],
        "component_weight_3body": component_weights[1],
        "component_weight_4body": component_weights[2],
        "seed": seed
    })

    # Train model using the fit method with wandb integration
    print(f"Training model with config: {model_config}")
    print(f"Component weights: {component_weights}")

    # Train with wandb logging
    history = model.fit(
        database,
        num_epochs=train_config.get("num_epochs", 100),
        batch_size=train_config.get("batch_size", 16),
        learning_rate=train_config.get("learning_rate", 1e-2),
        lr_decay_rate=train_config.get("lr_decay_rate", 0.005),
        verbose=1,
        use_wandb=True,
        wandb_project="raffle_inverse_design_sweep",
        wandb_run_name=f"run_{wandb.run.id}",
        wandb_config=model_config
    )

    # Log final training loss
    wandb.log({"final_training_loss": history[-1]})

    # Save model
    model_path = Path(f"model_{wandb.run.id}.pth")
    save_dict = {
        "state_dict": model.state_dict(),
        "config": {"species_list": species_list, **model_config}
    }
    torch.save(save_dict, model_path)
    wandb.save(str(model_path))

    # Load diamond structure for inverse design
    _diamond_path = os.path.join(_script_dir, "diamond.xyz")
    diamond = read(_diamond_path)

    # Create perturbed structure (replace one C with Si)
    perturbed_diamond = diamond.copy()
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
        'num_steps': 200,
    }

    # Define step sizes to test (from 1e-2 to 1e0)
    step_sizes = [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0]

    # Run inverse design for each step size
    print(f"Running inverse design for {len(step_sizes)} step sizes...")
    step_results = []

    for step_size in step_sizes:
        print(f"  Step size: {step_size}")
        result = run_inverse_design_for_step_size(
            model, target_fingerprint, perturbed_diamond,
            diamond, step_size, inverse_config
        )
        step_results.append(result)

        # Log individual step size results
        if result["atoms"] is not None:
            wandb.log({
                f"step_size_{step_size:.3f}_final_pos_diff": result["final_position_difference"],
                f"step_size_{step_size:.3f}_best_pos_diff": result["best_position_difference"],
            })

    # Find best overall result
    valid_results = [r for r in step_results if r["atoms"] is not None]
    if valid_results:
        best_result = min(valid_results, key=lambda x: x["best_position_difference"])
        initial_position_diff = step_results[0]["initial_position_difference"]
    else:
        # If all failed, return a large number
        wandb.finish()
        return 999.0

    # Log best metrics
    wandb.log({
        "initial_position_difference": initial_position_diff,
        "best_final_position_difference": best_result["final_position_difference"],
        "best_position_difference": best_result["best_position_difference"],
        "best_step_size": best_result["step_size"],
        "position_difference_improvement": initial_position_diff - best_result["best_position_difference"]
    })

    # Create and save plots
    temp_dir = Path(f"temp_{wandb.run.id}")
    temp_dir.mkdir(exist_ok=True)

    # Plot 1: Position difference traces for all step sizes
    traces_dict = {}
    for result in step_results:
        if result["trace_data"]:
            traces_dict[result["step_size"]] = result["trace_data"]

    if traces_dict:
        multi_step_plot_path = temp_dir / "position_difference_multi_step.png"
        plot_position_difference_multi_step(traces_dict, initial_position_diff, multi_step_plot_path)
        wandb.log({"position_difference_multi_step_plot": wandb.Image(str(multi_step_plot_path))})

    # Plot 2: Step size sweep summary
    step_size_data = []
    for result in step_results:
        if result["atoms"] is not None:
            step_size_data.append({
                "step_size": result["step_size"],
                "final_position_difference": result["final_position_difference"],
                "best_position_difference": result["best_position_difference"]
            })

    if step_size_data:
        sweep_summary_path = temp_dir / "step_size_sweep_summary.png"
        plot_step_size_summary(step_size_data, initial_position_diff, sweep_summary_path)
        wandb.log({"step_size_sweep_summary": wandb.Image(str(sweep_summary_path))})

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
        wandb.log({"fingerprint_comparison_plot": wandb.Image(str(fingerprint_plot_path))})

    # Also log the best structure
    if best_result["atoms"] is not None:
        best_structure_path = temp_dir / "best_structure.xyz"
        write(best_structure_path, best_result["atoms"])
        wandb.save(str(best_structure_path))

    # Log summary table of step size results
    if step_size_data:
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

    # Clean up
    import shutil
    shutil.rmtree(temp_dir, ignore_errors=True)
    if model_path.exists():
        os.remove(model_path)

    # Finish wandb run
    wandb.finish()
    return best_result["best_position_difference"]


# %%
# Run sweep
if __name__ == "__main__":
    sweep_id = wandb.sweep(sweep_config, project="raffle_inverse_design_sweep")
    wandb.agent(sweep_id, function=train_and_evaluate, count=200)
