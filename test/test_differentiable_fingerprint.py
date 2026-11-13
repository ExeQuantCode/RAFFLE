import math

import numpy as np
import pytest
import torch
from ase import Atoms
from ase.build import bulk

from raffle.fingerprint.differentiable import (
    DifferentiableFingerprint,
    FingerprintCalculator,
    FingerprintRelaxation,
)
from raffle.raffle import Geom_Rw


def test_periodic_fingerprint_and_autograd_gradients():
    atoms = bulk("C", "diamond", a=3.567, cubic=False)
    descriptor = DifferentiableFingerprint(["C", "Si"], skin=0.25)
    positions = torch.tensor(atoms.positions, dtype=torch.float64, requires_grad=True)
    probabilities = descriptor.probabilities_from_atoms(atoms).requires_grad_()

    fingerprint = descriptor(positions, atoms.cell.array, probabilities, atoms.pbc)
    scalar = sum(
        body.square().sum()
        for body in (
            fingerprint.two_body,
            fingerprint.three_body,
            fingerprint.four_body,
        )
    )
    position_gradient, species_gradient = torch.autograd.grad(
        scalar, (positions, probabilities)
    )

    assert fingerprint.two_body.shape == (221, 3)
    assert fingerprint.three_body.shape == (65, 2)
    assert fingerprint.four_body.shape == (65, 2)
    assert torch.isfinite(position_gradient).all()
    assert torch.isfinite(species_gradient).all()
    assert torch.linalg.vector_norm(position_gradient) > 0
    assert torch.linalg.vector_norm(species_gradient) > 0
    assert torch.allclose(
        fingerprint.two_body.sum(dim=0),
        torch.tensor([1.0, 0.0, 0.0], dtype=torch.float64),
        atol=1e-10,
    )


def test_position_and_species_logit_gradcheck():
    atoms = Atoms(
        symbols=["C", "Si", "C"],
        positions=[[0.0, 0.0, 0.0], [1.25, 0.1, 0.0], [0.2, 1.35, 0.45]],
        cell=[10.0, 10.0, 10.0],
        pbc=False,
    )
    descriptor = DifferentiableFingerprint(["C", "Si"])
    positions = torch.tensor(atoms.positions, dtype=torch.float64, requires_grad=True)
    logits = torch.tensor(
        [[0.2, -0.1], [-0.3, 0.4], [0.1, 0.3]],
        dtype=torch.float64,
        requires_grad=True,
    )

    def objective(position_values, species_logits):
        fingerprint = descriptor(
            position_values,
            atoms.cell.array,
            species_logits.softmax(dim=-1),
            atoms.pbc,
        )
        return fingerprint.two_body.square().sum() + fingerprint.three_body.square().sum()

    assert torch.autograd.gradcheck(
        objective, (positions, logits), eps=1e-6, atol=5e-4, rtol=5e-3
    )


def test_zero_skin_rebuilds_neighbor_topology_after_position_change():
    atoms = Atoms("CC", positions=[[0.0, 0.0, 0.0], [7.0, 0.0, 0.0]], cell=[20.0] * 3, pbc=False)
    descriptor = DifferentiableFingerprint(["C"], skin=0.0)
    first = descriptor.generate(atoms)
    atoms.positions[1, 0] = 5.9
    second = descriptor.generate(atoms)

    assert first.unnormalized["two_body"].sum() == 0.0
    assert second.unnormalized["two_body"].sum() > 0.0


def test_calculator_returns_negative_position_gradient_and_resets_species():
    target_atoms = bulk("C", "diamond", a=3.567, cubic=False)
    atoms = target_atoms.copy()
    atoms.positions[0, 0] += 0.08
    descriptor = DifferentiableFingerprint(["C", "Si"], skin=0.25)
    target = descriptor.generate(target_atoms)
    calculator = FingerprintCalculator(descriptor, target)
    atoms.calc = calculator

    initial_energy = atoms.get_potential_energy()
    forces = atoms.get_forces()
    assert forces.shape == atoms.positions.shape
    assert np.isfinite(forces).all()
    atoms.positions += 1e-4 * forces
    assert atoms.get_potential_energy() < initial_energy

    probabilities = descriptor.probabilities_from_atoms(atoms).numpy()
    probabilities[0] = [0.0, 1.0]
    calculator.set_species_probabilities(probabilities)
    assert atoms.get_potential_energy() >= 0.0


def test_exact_composition_and_fixed_atoms_are_preserved():
    target_atoms = Atoms(
        symbols=["C", "Si", "C"],
        positions=[[0.0, 0.0, 0.0], [1.8, 0.0, 0.0], [0.0, 1.8, 0.0]],
        cell=[8.0, 8.0, 8.0],
        pbc=False,
    )
    initial = target_atoms.copy()
    initial.set_chemical_symbols(["Si", "C", "C"])
    initial.positions[1, 0] += 0.2
    descriptor = DifferentiableFingerprint(["C", "Si"], skin=0.2)
    target = descriptor.generate(target_atoms)
    relaxation = FingerprintRelaxation(
        descriptor,
        target,
        target_composition={"C": 2, "Si": 1},
        composition_constraint=math.inf,
        fixed_position_mask=[True, False, False],
        fixed_species_mask=[False, True, False],
        position_learning_rate=0.001,
        species_learning_rate=0.01,
    )
    result = relaxation.run(
        initial, max_steps=3, trajectory_interval=1, final_refine_steps=0
    )

    assert result.composition == {"C": 2, "Si": 1}
    assert np.array_equal(result.structure.positions[0], initial.positions[0])
    assert result.structure[1].symbol == initial[1].symbol
    assert np.array_equal(result.trajectory[-1].positions, result.structure.positions)
    assert result.trajectory[-1].get_chemical_symbols() == result.structure.get_chemical_symbols()
    assert all(
        np.array_equal(frame.positions[0], initial.positions[0])
        and frame[1].symbol == initial[1].symbol
        for frame in result.trajectory
    )
    assert torch.allclose(
        result.species_probabilities.sum(dim=1),
        torch.ones(3, dtype=torch.float64),
        atol=1e-8,
    )


def test_invalid_exact_composition_is_rejected():
    atoms = bulk("C", "diamond", a=3.567, cubic=False)
    descriptor = DifferentiableFingerprint(["C", "Si"])
    target = descriptor.generate(atoms)
    relaxation = FingerprintRelaxation(
        descriptor,
        target,
        target_composition={"C": 1, "Si": 0},
        composition_constraint="exact",
    )
    with pytest.raises(ValueError, match="sum to atom count"):
        relaxation.run(atoms, max_steps=1, final_refine_steps=0)


def test_non_integer_exact_composition_is_rejected():
    atoms = Atoms("CC", positions=[[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]])
    descriptor = DifferentiableFingerprint(["C", "Si"])
    target = descriptor.generate(atoms)
    relaxation = FingerprintRelaxation(
        descriptor,
        target,
        target_composition={"C": 1.5, "Si": 0.5},
        composition_constraint="exact",
    )

    with pytest.raises(ValueError, match="non-negative integers"):
        relaxation.run(atoms, max_steps=1, final_refine_steps=0)


def test_fortran_fingerprint_parity():
    from raffle.raffle import Raffle__Distribs_Container

    atoms = bulk("C", "diamond", a=3.567, cubic=True)
    basis = Geom_Rw.basis(atoms)
    assert np.allclose(basis.toase().positions, atoms.positions)
    parameters = Raffle__Distribs_Container.distribs_container_type()
    reference = parameters.generate_fingerprint(atoms)
    descriptor = DifferentiableFingerprint(["C"], parameters=parameters, skin=0.25)
    actual = descriptor.generate(atoms)

    assert np.max(np.abs(reference[0] - actual.two_body.detach().numpy())) < 2e-6
    assert np.max(np.abs(reference[1] - actual.three_body.detach().numpy())) < 2e-6
    assert np.max(np.abs(reference[2] - actual.four_body.detach().numpy())) < 2e-3
