"""
RAFFLE Neural Network Fingerprint Example
==========================================

This example demonstrates the full workflow of the RAFFLE neural network
fingerprint module:

  1. Loading/constructing atomic structures
  2. Forward inference to obtain descriptor fingerprints
  3. Inverse design to reconstruct/optimise a structure
  4. Atom mask functionality (fixed vs optimisable atoms)

Usage:
    python nn_fingerprint_example.py

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
    # Import RAFFLE NN fingerprint module
    from raffle import NNFingerprint

    print("=" * 60)
    print("RAFFLE Neural Network Fingerprint Example")
    print("=" * 60)

    # -------------------------------------------------------------------------
    # Step 1: Create atomic structures
    # -------------------------------------------------------------------------
    print("\n--- Step 1: Creating atomic structures ---")

    # Diamond cubic carbon (8 atoms in a 2x1x1 supercell)
    carbon_diamond = bulk('C', 'diamond', a=3.567) * (2, 1, 1)
    carbon_diamond.info['energy'] = -72.0
    print(f"  Structure 1: C diamond, {len(carbon_diamond)} atoms")
    print(f"  Cell: {carbon_diamond.cell.lengths()}")

    # Perturbed diamond structure
    carbon_perturbed = carbon_diamond.copy()
    positions = carbon_perturbed.get_positions()
    rng = np.random.default_rng(42)
    positions += rng.normal(0, 0.2, positions.shape)
    carbon_perturbed.set_positions(positions)
    carbon_perturbed.info['energy'] = -71.0
    print(f"  Structure 2: C diamond (perturbed), {len(carbon_perturbed)} atoms")

    # Silicon structure
    silicon = bulk('Si', 'diamond', a=5.431) * (2, 1, 1)
    silicon.info['energy'] = -40.0
    print(f"  Structure 3: Si diamond, {len(silicon)} atoms")

    # -------------------------------------------------------------------------
    # Step 2: Initialise the NN fingerprint module
    # -------------------------------------------------------------------------
    print("\n--- Step 2: Initialising NN fingerprint module ---")

    nn = NNFingerprint(
        species_list=['C'],
        max_atoms=8,
        hidden_layer_sizes=[64, 32],
        learning_rate=0.001,
    )
    print(f"  Species: {nn.species_list}")
    print(f"  Max atoms: {nn.max_atoms}")
    print(f"  Input dim: {nn.input_dim}")

    use_simple_fingerprint = False
    compute_fp = (
        nn.compute_fingerprint_direct if use_simple_fingerprint
        else nn.compute_fingerprint
    )
    data_path = Path(__file__).resolve().parents[1] / 'data' / 'carbon.xyz'
    print(
        "  Fingerprint source: "
        + ("simple numpy radial fingerprint" if use_simple_fingerprint
           else "RAFFLE descriptor fingerprint")
    )

    # -------------------------------------------------------------------------
    # Step 3: Compute descriptor fingerprints
    # -------------------------------------------------------------------------
    print("\n--- Step 3: Computing RAFFLE descriptor fingerprints ---")

    fp1 = compute_fp(carbon_diamond)
    fp2 = compute_fp(carbon_perturbed)
    print(f"  Fingerprint dim: {nn.fingerprint_dim}")
    print(f"  FP1 (diamond) max: {np.max(np.abs(fp1)):.6f}")
    print(f"  FP2 (perturbed) max: {np.max(np.abs(fp2)):.6f}")

    # Distance between fingerprints
    fp_distance = np.sqrt(np.sum((fp1 - fp2) ** 2))
    print(f"  L2 distance between FP1 and FP2: {fp_distance:.6f}")
    print(f"  (Different structures should have different fingerprints)")

    # -------------------------------------------------------------------------
    # Step 4: Train the neural network
    # -------------------------------------------------------------------------
    print("\n--- Step 4: Training neural network ---")

    training_structures = [carbon_diamond, carbon_perturbed]
    from ase.io import read
    training_structures = read(data_path, index=":")
    loss_history = nn.train(
        structures=training_structures,
        num_epochs=10000,
        verbose=1,
        use_simple_fingerprint=use_simple_fingerprint,
    )
    print(f"  Final training loss: {loss_history[-1]:.8f}")
    print(f"  Network trained: {nn.is_trained}")

    # -------------------------------------------------------------------------
    # Step 5: Forward inference
    # -------------------------------------------------------------------------
    print("\n--- Step 5: Forward inference (predict) ---")

    predicted_fp = nn.predict(
        carbon_diamond,
        use_simple_fingerprint=use_simple_fingerprint,
    )
    true_fp = fp1

    prediction_error = np.mean((predicted_fp - true_fp) ** 2)
    print(f"  Prediction MSE: {prediction_error:.8f}")
    print(f"  Predicted FP max: {np.max(np.abs(predicted_fp)):.6f}")
    print(f"  True FP max: {np.max(np.abs(true_fp)):.6f}")

    # -------------------------------------------------------------------------
    # Step 6: Inverse design with atom masking
    # -------------------------------------------------------------------------
    print("\n--- Step 6: Inverse design with atom mask ---")

    # Target: fingerprint of the ideal diamond structure
    target_fp = fp1.copy()

    # Start from perturbed structure
    test_structure = carbon_perturbed.copy()
    original_positions = test_structure.get_positions().copy()
    structures = []
    structures.append(test_structure.copy())

    # Define atom mask: fix first atom, allow rest to move
    n_atoms = len(test_structure)
    fixed_atoms = np.zeros(n_atoms, dtype=bool)
    fixed_atoms[0] = True  # Fix first atom
    print(f"  Total atoms: {n_atoms}")
    print(f"  Fixed atoms: {np.sum(fixed_atoms)} (atom indices: "
          f"{np.where(fixed_atoms)[0]})")
    print(f"  Movable atoms: {np.sum(~fixed_atoms)}")

    # Run inverse design
    optimised = nn.inverse_design(
        target_fingerprint=target_fp,
        atoms=test_structure,
        fixed_atoms=fixed_atoms,
        num_steps=1000,
        step_size=0.002,
        verbose=1,
        use_simple_fingerprint=use_simple_fingerprint,
    )

    # Save the optimised structure for visualization
    from ase.io import write
    structures.append(optimised.copy())
    write("structures.traj", structures)

    # Verify fixed atoms didn't move
    final_positions = optimised.get_positions()
    fixed_displacement = np.linalg.norm(
        final_positions[0] - original_positions[0]
    )
    movable_displacement = np.mean([
        np.linalg.norm(final_positions[i] - original_positions[i])
        for i in range(1, n_atoms)
    ])

    print(f"\n  Fixed atom displacement: {fixed_displacement:.8f} A")
    print(f"  Average movable atom displacement: "
          f"{movable_displacement:.8f} A")
    print(f"  Fixed atom correctly constrained: "
          f"{fixed_displacement < 1e-6}")

    # -------------------------------------------------------------------------
    # Summary
    # -------------------------------------------------------------------------
    print("\n" + "=" * 60)
    print("Summary")
    print("=" * 60)
    print(f"  Network architecture: "
          f"{nn.input_dim} -> {nn._hidden_sizes} -> {nn.fingerprint_dim}")
    print(f"  Training structures: {len(training_structures)}")
    print(f"  Final training loss: {loss_history[-1]:.8f}")
    print(f"  Forward prediction MSE: {prediction_error:.8f}")
    print(f"  Inverse design loss: "
            f"{np.mean((compute_fp(optimised) - target_fp)**2):.8e}")
    print(f"  Atom masking verified: {fixed_displacement < 1e-6}")
    print("\nAll operations completed successfully!")


if __name__ == '__main__':
    main()
