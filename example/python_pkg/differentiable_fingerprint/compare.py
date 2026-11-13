"""Generate a differentiable fingerprint and compare with the Fortran wrapper."""

import numpy as np
import torch
from ase.build import bulk

from raffle.fingerprint.differentiable import DifferentiableFingerprint
from raffle.raffle import Raffle__Distribs_Container


atoms = bulk("C", "diamond", a=3.567, cubic=True)
parameters = Raffle__Distribs_Container.distribs_container_type()
descriptor = DifferentiableFingerprint(["C"], parameters=parameters)
reference = parameters.generate_fingerprint(atoms)
positions = torch.tensor(atoms.positions, dtype=torch.float64, requires_grad=True)
species = descriptor.probabilities_from_atoms(atoms).requires_grad_()
fingerprint = descriptor(positions, atoms.cell.array, species, atoms.pbc)

for name, expected, actual in zip(
    ("2-body", "3-body", "4-body"),
    reference,
    (fingerprint.two_body, fingerprint.three_body, fingerprint.four_body),
):
    difference = np.abs(expected - actual.detach().numpy())
    print(f"{name}: max_abs={difference.max():.6g}, mean_abs={difference.mean():.6g}")

objective = sum(
    component.square().sum()
    for component in (
        fingerprint.two_body,
        fingerprint.three_body,
        fingerprint.four_body,
    )
)
position_gradient, species_gradient = torch.autograd.grad(objective, (positions, species))
print("position-gradient norm:", position_gradient.norm().item())
print("species-gradient norm:", species_gradient.norm().item())
