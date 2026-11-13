"""Jointly relax perturbed diamond positions and species assignments."""

import numpy as np
from ase.build import bulk
from ase.io import write

from raffle.fingerprint.differentiable import (
    DifferentiableFingerprint,
    FingerprintRelaxation,
)

if __package__:
    from .plot_comparison import plot_fingerprints
else:
    from plot_comparison import plot_fingerprints


target_atoms = bulk("C", "diamond", a=3.567, cubic=True)
descriptor = DifferentiableFingerprint(["C", "Si"])
target = descriptor.generate(target_atoms)
initial = target_atoms.copy()
initial.positions[0] += np.array([1.2, 0.0, 0.0])
initial[1].symbol = "Si"

coarse_relaxation = FingerprintRelaxation(
    descriptor,
    target,
    target_composition={"C": len(target_atoms), "Si": 0},
    composition_constraint="exact",
    loss_weights=(1000.0, 0.0, 0.0),
    radial_smoothing=0.5,
    position_learning_rate=0.02,
)
print("Coarse radial relaxation:")
coarse_result = coarse_relaxation.run(
    initial, max_steps=50, final_refine_steps=0, tolerance=0, verbose=True
)

relaxation = FingerprintRelaxation(
    descriptor,
    target,
    target_composition={"C": len(target_atoms), "Si": 0},
    composition_constraint="exact",
    position_learning_rate=0.01,
)
print("Fine fingerprint relaxation:")
result = relaxation.run(
    coarse_result.structure, max_steps=20, tolerance=0, verbose=True
)
write("diamond_inverse_design.extxyz", result.structure)

write("diamond.traj", coarse_result.trajectory + result.trajectory[1:])

# print(result.species_probabilities_history)

print("fine initial loss:", result.loss_history[0]["total"])
print("final loss:", result.loss)
print("composition:", result.composition)
print("steps:", result.n_steps)

plot_fingerprints(
    [value.detach().cpu().numpy() for value in (target.two_body, target.three_body, target.four_body)],
    [value.detach().cpu().numpy() for value in (
        result.fingerprint.two_body,
        result.fingerprint.three_body,
        result.fingerprint.four_body,
    )],
    {name: grid.detach().cpu().numpy() for name, grid in target.grids.items()},
    descriptor.elements,
    "Target",
    "Final",
    "diamond_fingerprint_comparison.png",
)
