"""Jointly relax perturbed diamond positions and species assignments."""

import numpy as np
from ase.build import bulk
from ase.io import write

from raffle.fingerprint.differentiable import (
    DifferentiableFingerprint,
    FingerprintRelaxation,
)


target_atoms = bulk("C", "diamond", a=3.567, cubic=True)
descriptor = DifferentiableFingerprint(["C", "Si"])
target = descriptor.generate(target_atoms)
initial = target_atoms.copy()
initial.positions[0] += np.array([1.2, 0.0, 0.0])
initial[1].symbol = "Si"

relaxation = FingerprintRelaxation(
    descriptor,
    target,
    target_composition={"C": len(target_atoms), "Si": 0},
    composition_constraint="exact",
)
result = relaxation.run(initial, max_steps=500)
write("diamond_inverse_design.extxyz", result.structure)

write("diamond.traj", result.trajectory)

print(result.species_probabilities_history)

print("initial loss:", result.loss_history[0]["total"])
print("final loss:", result.loss)
print("composition:", result.composition)
print("steps:", result.n_steps)
