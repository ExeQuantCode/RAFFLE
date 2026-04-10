"""
RAFFLE Graph Neural Network Fingerprint Example
================================================

This example demonstrates the full workflow of the RAFFLE GNN fingerprint
module, which uses message-passing on molecular graphs:

  1. Loading/constructing atomic structures
  2. Converting structures to molecular graphs (atoms→vertices, bonds→edges)
  3. Forward inference to obtain descriptor fingerprints via GNN
  4. Inverse design to reconstruct/optimise a structure
  5. Atom mask functionality (fixed vs optimisable atoms)

Usage:
    python gnn_fingerprint_example.py

Requirements:
    - raffle (pip install .)
    - ase
    - numpy
"""
import numpy as np
from pathlib import Path
from ase import Atoms
from ase.build import bulk


def main():
    from raffle import GNNFingerprint

    print("=" * 60)
    print("RAFFLE Graph Neural Network Fingerprint Example")
    print("=" * 60)

    # -------------------------------------------------------------------------
    # Step 1: Create atomic structures
    # -------------------------------------------------------------------------
    print("\n--- Step 1: Creating atomic structures ---")

    carbon_diamond = bulk('C', 'diamond', a=3.567) * (2, 1, 1)
    carbon_diamond.info['energy'] = -72.0
    print(f"  Structure 1: C diamond, {len(carbon_diamond)} atoms")
    print(f"  Cell: {carbon_diamond.cell.lengths()}")

    carbon_perturbed = carbon_diamond.copy()
    positions = carbon_perturbed.get_positions()
    rng = np.random.default_rng(42)
    positions += rng.normal(0, 0.2, positions.shape)
    carbon_perturbed.set_positions(positions)
    carbon_perturbed.info['energy'] = -71.0
    print(f"  Structure 2: C diamond (perturbed), {len(carbon_perturbed)} atoms")

    # -------------------------------------------------------------------------
    # Step 2: Initialise the GNN fingerprint module
    # -------------------------------------------------------------------------
    print("\n--- Step 2: Initialising GNN fingerprint module ---")

    gnn = GNNFingerprint(
        species_list=['C'],
        bond_cutoff=6.0,
        gnn_hidden_sizes=[32, 32],
        learning_rate=0.001,
    )
    print(f"  Species: {gnn.species_list}")
    print(f"  Bond cutoff: {gnn.bond_cutoff} A")
    print(f"  Vertex features: {gnn.num_vertex_features}")

    use_simple_fingerprint = False
    compute_fp = (
        gnn.compute_fingerprint_direct if use_simple_fingerprint
        else gnn.compute_fingerprint
    )
    data_path = Path(__file__).resolve().parents[1] / 'data' / 'carbon.xyz'
    print(
        "  Fingerprint source: "
        + ("simple numpy radial fingerprint" if use_simple_fingerprint
           else "RAFFLE descriptor fingerprint")
    )

    # -------------------------------------------------------------------------
    # Step 3: Convert structures to molecular graphs
    # -------------------------------------------------------------------------
    print("\n--- Step 3: Converting structures to molecular graphs ---")

    vf1, adj1 = gnn.atoms_to_graph(carbon_diamond)
    vf2, adj2 = gnn.atoms_to_graph(carbon_perturbed)

    n_edges1 = int((np.sum(adj1) - len(carbon_diamond)) / 2)
    n_edges2 = int((np.sum(adj2) - len(carbon_perturbed)) / 2)
    print(f"  Graph 1: {vf1.shape[0]} vertices, {n_edges1} edges")
    print(f"    Vertex feature shape: {vf1.shape}")
    print(f"    Adjacency matrix shape: {adj1.shape}")
    print(f"  Graph 2: {vf2.shape[0]} vertices, {n_edges2} edges")

    # -------------------------------------------------------------------------
    # Step 4: Compute descriptor fingerprints
    # -------------------------------------------------------------------------
    print("\n--- Step 4: Computing RAFFLE descriptor fingerprints ---")

    fp1 = compute_fp(carbon_diamond)
    fp2 = compute_fp(carbon_perturbed)
    print(f"  Fingerprint dim: {gnn.fingerprint_dim}")
    print(f"  FP1 (diamond) max: {np.max(np.abs(fp1)):.6f}")
    print(f"  FP2 (perturbed) max: {np.max(np.abs(fp2)):.6f}")

    fp_distance = np.sqrt(np.sum((fp1 - fp2) ** 2))
    print(f"  L2 distance between FP1 and FP2: {fp_distance:.6f}")

    # -------------------------------------------------------------------------
    # Step 5: Train the GNN
    # -------------------------------------------------------------------------
    print("\n--- Step 5: Training GNN on structures ---")

    from ase.io import read
    max_training_structures = 32
    num_epochs = 200
    training_structures = read(data_path, index=":")[:max_training_structures]
    print(f"  Training structures used: {len(training_structures)}")
    print(f"  Training epochs: {num_epochs}")
    loss_history = gnn.train(
        structures=training_structures,
        num_epochs=num_epochs,
        verbose=1,
        use_simple_fingerprint=use_simple_fingerprint,
    )
    print(f"  Final training loss: {loss_history[-1]:.8f}")
    print(f"  GNN trained: {gnn.is_trained}")

    # -------------------------------------------------------------------------
    # Step 6: Forward inference via GNN
    # -------------------------------------------------------------------------
    print("\n--- Step 6: Forward inference (predict) ---")

    predicted_fp = gnn.predict(carbon_diamond)
    true_fp = fp1

    prediction_error = np.mean((predicted_fp - true_fp) ** 2)
    print(f"  Prediction MSE: {prediction_error:.8f}")
    print(f"  Predicted FP max: {np.max(np.abs(predicted_fp)):.6f}")
    print(f"  True FP max: {np.max(np.abs(true_fp)):.6f}")

    # -------------------------------------------------------------------------
    # Step 7: Inverse design with atom masking
    # -------------------------------------------------------------------------
    print("\n--- Step 7: Inverse design with atom mask ---")

    target_fp = fp1.copy()
    test_structure = carbon_perturbed.copy()
    original_positions = test_structure.get_positions().copy()
    structures = []
    structures.append(test_structure.copy())

    n_atoms = len(test_structure)
    fixed_atoms = np.zeros(n_atoms, dtype=bool)
    #fixed_atoms[0] = True
    print(f"  Total atoms: {n_atoms}")
    print(f"  Fixed atoms: {np.sum(fixed_atoms)} (atom indices: "
          f"{np.where(fixed_atoms)[0]})")
    print(f"  Movable atoms: {np.sum(~fixed_atoms)}")

    optimised = gnn.inverse_design(
        target_fingerprint=target_fp,
        atoms=test_structure,
        fixed_atoms=fixed_atoms,
        num_steps=100,
        step_size=0.005,
        verbose=1,
        use_simple_fingerprint=use_simple_fingerprint,
    )

    # Save the optimised structure for visualization
    from ase.io import write
    structures.append(optimised.copy())
    write("structures.traj", structures)

    final_positions = optimised.get_positions()
    if np.any(fixed_atoms):
        fixed_displacement = np.linalg.norm(
            final_positions[np.where(fixed_atoms)[0][0]]
            - original_positions[np.where(fixed_atoms)[0][0]]
        )
        fixed_constraint_ok = fixed_displacement < 1e-6
    else:
        fixed_displacement = 0.0
        fixed_constraint_ok = True
    movable_indices = np.where(~fixed_atoms)[0]
    if movable_indices.size > 0:
        movable_displacement = np.mean([
            np.linalg.norm(final_positions[i] - original_positions[i])
            for i in movable_indices
        ])
    else:
        movable_displacement = 0.0

    print(f"\n  Fixed atom displacement: {fixed_displacement:.8f} A")
    print(f"  Average movable atom displacement: "
          f"{movable_displacement:.8f} A")
    print(f"  Fixed atom correctly constrained: "
          f"{fixed_constraint_ok}")

    # -------------------------------------------------------------------------
    # Summary
    # -------------------------------------------------------------------------
    print("\n" + "=" * 60)
    print("GNN Fingerprint Example Complete")
    print("=" * 60)
    print(f"  Fingerprint dimension: {gnn.fingerprint_dim}")
    print(f"  Training loss (initial → final): "
          f"{loss_history[0]:.6f} → {loss_history[-1]:.6f}")
    print(f"  Prediction MSE: {prediction_error:.8f}")
    print(f"  Inverse design: fixed atom unmoved = "
            f"{fixed_constraint_ok}")
    print("=" * 60)


if __name__ == "__main__":
    main()
