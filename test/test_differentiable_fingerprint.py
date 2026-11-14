import ast
import math
import os
from pathlib import Path
import runpy
import subprocess
import sys

import numpy as np
import pytest
import torch
from ase import Atoms
from ase.build import bulk
from ase.io import read, write
from ase.neighborlist import neighbor_list

from raffle.fingerprint.differentiable import (
    DifferentiableFingerprint,
    FingerprintCalculator,
    FingerprintRelaxation,
)
from raffle.fingerprint.differentiable.descriptor import _fingerprint_loss
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


def test_four_body_position_and_multispecies_logit_gradcheck():
    atoms = Atoms(
        "CSiOC",
        positions=[[0, 0, 0], [1.35, 0.1, 0.1], [0.2, 1.4, 0.3], [0.1, 0.3, 2.8]],
        cell=[12.0] * 3, pbc=False,
    )
    descriptor = DifferentiableFingerprint(["C", "O", "Si"], body_orders=(4,))
    positions = torch.tensor(atoms.positions, dtype=torch.float64, requires_grad=True)
    logits = torch.tensor(
        [[0.2, 0.1, -0.3], [-0.3, 0.4, 0.2], [0.1, -0.2, 0.3], [0.2, 0.3, -0.1]],
        dtype=torch.float64, requires_grad=True,
    )

    def objective(position_values, species_logits):
        fingerprint = descriptor(
            position_values, atoms.cell.array, species_logits.softmax(dim=-1), atoms.pbc
        )
        return fingerprint.four_body.square().sum()

    assert float(objective(positions, logits).detach()) > 0
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


def test_cached_neighbor_groups_match_fresh_descriptor_after_topology_changes():
    atoms = Atoms(
        "CCCCC", positions=[[9, 9, 9], [0, 0, 0], [1.3, 0.1, 0], [0.2, 1.3, 0.1], [0, 0.1, 2.7]],
        cell=[20.0] * 3, pbc=False,
    )
    cached = DifferentiableFingerprint(["C"], skin=0.2)
    for action in ("initial", "moved", "periodic"):
        if action == "moved":
            atoms.positions[0] = [1.1, 1.1, 1.1]
        elif action == "periodic":
            atoms.set_cell([6.0] * 3)
            atoms.pbc = True
        actual = cached.generate(atoms)
        expected = DifferentiableFingerprint(["C"], skin=0.2).generate(atoms)
        for name in ("two_body", "three_body", "four_body"):
            torch.testing.assert_close(getattr(actual, name), getattr(expected, name))


@pytest.mark.parametrize(
    "device",
    ["cpu"] + (["mps"] if torch.backends.mps.is_available() else [])
    + (["cuda"] if torch.cuda.is_available() else []),
)
def test_periodic_wrap_reuses_topology_and_matches_fresh_gradients(device):
    cell = np.array([[4.0, 0.0, 0.0], [1.3, 4.2, 0.0], [0.3, 0.8, 10.0]])
    atoms = Atoms(
        "CSiCC", scaled_positions=[[0.98, 0.15, 0.2], [0.1, 0.4, 0.2],
                                   [0.3, 0.15, 0.22], [0.65, 0.65, 0.24]],
        cell=cell, pbc=[True, True, False],
    )
    descriptor = DifferentiableFingerprint(["C", "Si"], skin=0.8)
    dtype = torch.float64 if device == "cpu" else torch.float32
    logits_data = torch.tensor(
        [[0.4, -0.1], [0.2, 0.3], [0.1, 0.5], [0.5, 0.1]],
        dtype=dtype, device=device,
    )
    descriptor(
        torch.tensor(atoms.positions, dtype=dtype, device=device),
        cell, logits_data.softmax(dim=1), atoms.pbc,
    )
    reference_positions = descriptor._topology.reference_positions

    def evaluate(current):
        positions = torch.tensor(atoms.positions, dtype=dtype, device=device, requires_grad=True)
        logits = logits_data.clone().requires_grad_()
        probabilities = logits.softmax(dim=1)
        fp = current(positions, cell, probabilities, atoms.pbc)
        values = (fp.two_body, fp.three_body, fp.four_body)
        loss = sum(value.square().sum() for value in values)
        return values, torch.autograd.grad(loss, (positions, probabilities, logits))

    for move in (0.04 * cell[0] - cell[0], 0.01 * cell[0], cell[0]):
        atoms.positions[0] += move
        actual, actual_gradients = evaluate(descriptor)
        expected, expected_gradients = evaluate(DifferentiableFingerprint(["C", "Si"], skin=0.8))
        assert descriptor._topology.reference_positions is reference_positions
        for lhs, rhs in zip((*actual, *actual_gradients), (*expected, *expected_gradients)):
            torch.testing.assert_close(
                lhs, rhs, rtol=1e-7 if device == "cpu" else 1e-4,
                atol=1e-9 if device == "cpu" else 3e-6,
            )
    center, neighbor, shifts = descriptor._topology.pairs(
        torch.tensor(atoms.positions), torch.tensor(cell), atoms.pbc
    )
    assert not ((center == neighbor) & ~shifts.ne(0).any(dim=1)).any()
    assert ((center == neighbor) & shifts.ne(0).any(dim=1)).any()
    assert ((center == 0) & (neighbor == 0)).sum() > 1


def test_neighbor_skin_covers_cutoff_crossings_and_invalidates_configuration():
    atoms = Atoms("CC", positions=[[0, 0, 0], [6.2, 0, 0]], cell=[20.0] * 3, pbc=False)
    descriptor = DifferentiableFingerprint(["C"], skin=0.5, body_orders=(2,))
    assert descriptor.generate(atoms).two_body.sum() == 0
    reference_positions = descriptor._topology.reference_positions
    atoms.positions[1, 0] = 5.97
    cached = descriptor.generate(atoms)
    assert descriptor._topology.reference_positions is reference_positions
    fresh = DifferentiableFingerprint(["C"], skin=0.5, body_orders=(2,)).generate(atoms)
    torch.testing.assert_close(cached.two_body, fresh.two_body)
    assert cached.two_body.sum() > 0
    atoms.positions[1, 0] = 6.1
    leaving = descriptor.generate(atoms)
    assert descriptor._topology.reference_positions is reference_positions
    assert leaving.two_body.sum() == 0
    atoms.positions[1, 0] = 5.9
    descriptor.generate(atoms)
    assert descriptor._topology.reference_positions is not reference_positions

    atoms.positions[1, 0] = 6.8
    descriptor.generate(atoms)
    reference_positions = descriptor._topology.reference_positions
    descriptor.cutoff_max[0] = 8.0
    changed = descriptor.generate(atoms)
    expected = DifferentiableFingerprint(
        ["C"], parameters={"cutoff_max": [8.0, math.pi, math.pi]},
        skin=0.5, body_orders=(2,),
    ).generate(atoms)
    assert descriptor._topology.reference_positions is not reference_positions
    torch.testing.assert_close(changed.two_body, expected.two_body)
    assert changed.two_body.sum() > 0


@pytest.mark.parametrize("change", ["cell", "pbc", "atom_count", "skin"])
def test_topology_rebuilds_when_structure_metadata_changes(change):
    atoms = Atoms(
        "CC", positions=[[0, 0, 0], [1.4, 0.1, 0]], cell=[8.0] * 3, pbc=False
    )
    descriptor = DifferentiableFingerprint(["C"], skin=0.5, body_orders=(2,))
    descriptor.generate(atoms)
    old_reference = descriptor._topology.reference_positions
    if change == "cell":
        atoms.cell[0, 0] += 1e-13
    elif change == "pbc":
        atoms.pbc[0] = True
    elif change == "atom_count":
        atoms += Atoms("C", positions=[[0, 1.4, 0]])
    else:
        descriptor._topology.skin = 0.8
    actual = descriptor.generate(atoms)
    expected = DifferentiableFingerprint(["C"], skin=0.5, body_orders=(2,)).generate(atoms)
    assert descriptor._topology.reference_positions is not old_reference
    torch.testing.assert_close(actual.two_body, expected.two_body)


@pytest.mark.parametrize("sigma", [0.1, 2.0])
def test_radial_projection_matches_dense_reference_at_grid_edges(sigma):
    descriptor = DifferentiableFingerprint(
        ["C", "Si"], parameters={"sigma": [sigma, 0.1, 0.1]}, body_orders=(2,)
    )
    positions = torch.tensor(
        [[0., 0., 0.], [0.52, 0., 0.], [5.9, 0.12, 0.]],
        dtype=torch.float64, requires_grad=True,
    )
    logits = torch.tensor(
        [[0.2, -0.3], [0.4, 0.1], [-0.2, 0.3]],
        dtype=torch.float64, requires_grad=True,
    )
    probabilities = logits.softmax(dim=1)
    fp = descriptor(positions, np.diag([20.] * 3), probabilities, [False] * 3)
    center = torch.tensor([0, 0, 1, 1, 2, 2])
    neighbor = torch.tensor([1, 2, 0, 2, 0, 1])
    distances = torch.linalg.vector_norm(positions[neighbor] - positions[center], dim=1)
    denominator = probabilities.new_zeros((3, 2)).index_add(0, center, probabilities[neighbor])
    first, second = probabilities[center], probabilities[neighbor]
    weights = torch.stack(
        (first[:, 0] * second[:, 0] / denominator[center, 0],
         first[:, 0] * second[:, 1] / denominator[center, 1]
         + first[:, 1] * second[:, 0] / denominator[center, 0],
         first[:, 1] * second[:, 1] / denominator[center, 1]),
        dim=1,
    )
    grid = fp.grids["two_body"]
    eta = 1.0 / (2.0 * sigma ** 2)
    support_radius = math.ceil(math.sqrt(16.0 / eta) / descriptor.width[0])
    nearest = torch.floor((distances - descriptor.cutoff_min[0]) / descriptor.width[0] + 0.5)
    support = (torch.arange(grid.numel())[None, :] - nearest[:, None]).abs() <= support_radius
    projected = torch.exp(-eta * (distances[:, None] - grid[None, :]).square()) * support
    taper = torch.ones_like(grid)
    taper = torch.where(
        grid > descriptor.cutoff_max[0] - 0.25,
        0.5 * (1 + torch.cos(math.pi * (grid - descriptor.cutoff_max[0] + 0.25) / 0.25)),
        taper,
    )
    taper = torch.where(
        grid < descriptor.cutoff_min[0] + 0.25,
        0.5 * (1 + torch.cos(math.pi * (grid - descriptor.cutoff_min[0] - 0.25) / 0.25)),
        taper,
    )
    expected = (
        projected.T @ weights * math.sqrt(eta / math.pi) * taper[:, None]
        / grid[:, None].square()
    )
    actual = fp.unnormalized["two_body"]
    torch.testing.assert_close(actual, expected, rtol=1e-12, atol=1e-12)
    factors = torch.linspace(0.1, 1.0, grid.numel(), dtype=grid.dtype)[:, None]
    actual_grad = torch.autograd.grad(
        (actual * factors).sum(), (positions, logits), retain_graph=True
    )
    expected_grad = torch.autograd.grad((expected * factors).sum(), (positions, logits))
    for lhs, rhs in zip(actual_grad, expected_grad):
        torch.testing.assert_close(lhs, rhs, rtol=1e-11, atol=1e-11)


def test_candidate_topology_covers_larger_allowed_species_radii():
    atoms = Atoms(
        "CCCC", positions=[[0, 0, 0], [2.6, 0, 0], [0, 2.6, 0], [0, 0, 4.4]],
        cell=[20.0] * 3, pbc=False,
    )
    descriptor = DifferentiableFingerprint(
        ["C", "Si"],
        bond_radii={("C", "C"): 0.7, ("C", "Si"): 1.4, ("Si", "Si"): 1.4},
    )
    positions = torch.tensor(atoms.positions, dtype=torch.float64, requires_grad=True)
    logits = torch.tensor([[0.5, 0.1]] * len(atoms), dtype=torch.float64, requires_grad=True)
    fp = descriptor(positions, atoms.cell.array, logits.softmax(dim=1), atoms.pbc)
    assert fp.three_body.sum() > 0
    assert fp.four_body.sum() > 0
    gradients = torch.autograd.grad(
        fp.three_body.square().sum() + fp.four_body.square().sum(), (positions, logits)
    )
    assert all(torch.isfinite(gradient).all() for gradient in gradients)
    assert all(torch.linalg.vector_norm(gradient) > 0 for gradient in gradients)


def test_repulsion_excludes_isolated_atom_self_pair():
    atoms = Atoms("C", positions=[[0, 0, 0]], cell=[10.0] * 3, pbc=False)
    descriptor = DifferentiableFingerprint(["C", "Si"], body_orders=(2,))
    relaxation = FingerprintRelaxation(
        descriptor, descriptor.generate(atoms), loss_weights=(1.0, 0.0, 0.0),
        repulsion_weight=1.0,
    )
    result = relaxation.run(atoms, max_steps=1, final_refine_steps=0)

    assert result.loss_history[0]["repulsion"] == 0.0
    assert np.isfinite(result.loss)
    assert torch.isfinite(result.species_probabilities).all()
    np.testing.assert_array_equal(result.structure.positions, atoms.positions)


def test_repulsion_keeps_distinct_atom_interactions():
    atoms = Atoms("CC", positions=[[0, 0, 0], [0.4, 0, 0]], cell=[10.0] * 3, pbc=False)
    descriptor = DifferentiableFingerprint(
        ["C"], body_orders=(2,), bond_radii={("C", "C"): 1.0}
    )
    relaxation = FingerprintRelaxation(
        descriptor, descriptor.generate(atoms), loss_weights=(1.0, 0.0, 0.0),
        repulsion_weight=1.0,
    )
    result = relaxation.run(
        atoms, optimise_species=False, max_steps=1, final_refine_steps=0
    )

    assert result.loss_history[0]["repulsion"] == pytest.approx((0.7 - 0.4) ** 2)
    assert result.structure.get_distance(0, 1) > atoms.get_distance(0, 1)


def test_repulsion_keeps_periodic_self_images():
    atoms = Atoms("C", positions=[[0, 0, 0]], cell=[0.4, 10.0, 10.0], pbc=[True, False, False])
    descriptor = DifferentiableFingerprint(
        ["C"], body_orders=(2,), bond_radii={("C", "C"): 1.0}
    )
    relaxation = FingerprintRelaxation(
        descriptor, descriptor.generate(atoms), loss_weights=(1.0, 0.0, 0.0),
        repulsion_weight=1.0,
    )
    result = relaxation.run(
        atoms, optimise_species=False, max_steps=1, final_refine_steps=0
    )

    distances = neighbor_list("d", atoms, descriptor.cutoff_max[0] + descriptor._topology.skin)
    expected = np.square(np.maximum(0.7 - distances, 0.0)).mean()
    assert expected > 0
    assert result.loss_history[0]["repulsion"] == pytest.approx(expected)
    assert np.isfinite(result.structure.positions).all()


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


@pytest.mark.parametrize("constraint", [0.0, "exact"])
def test_exact_composition_and_fixed_atoms_are_preserved(constraint):
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
        composition_constraint=constraint,
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


def test_composition_temperature_scales_fractional_mismatch():
    atoms = Atoms("CC", positions=[[0, 0, 0], [1.5, 0, 0]], cell=[10] * 3)
    initial = atoms.copy()
    initial[1].symbol = "Si"
    descriptor = DifferentiableFingerprint(["C", "Si"], body_orders=(2,))
    target = descriptor.generate(atoms)

    def run(constraint):
        relaxation = FingerprintRelaxation(
            descriptor, target, composition_constraint=constraint,
            loss_weights=(1.0, 0.0, 0.0),
        )
        return relaxation.run(initial, optimise_positions=False, max_steps=1, final_refine_steps=0)

    unrestricted = [run(value) for value in (None, "none", math.inf)]
    default = FingerprintRelaxation(
        descriptor, target, loss_weights=(1.0, 0.0, 0.0)
    ).run(initial, optimise_positions=False, max_steps=1, final_refine_steps=0)
    assert default.loss_history[0]["total"] == pytest.approx(unrestricted[0].loss_history[0]["total"])
    for result in unrestricted:
        assert result.structure.get_chemical_symbols() == ["C", "Si"]
        assert "composition" not in result.loss_history[0]
        assert result.loss_history[0]["total"] == pytest.approx(unrestricted[0].loss_history[0]["total"])

    for temperature in (0.1, 2.0):
        result = run(temperature)
        record = result.loss_history[0]
        mismatch = record["composition"]
        assert mismatch > 0
        assert record["total"] == pytest.approx(
            unrestricted[0].loss_history[0]["total"] + mismatch / temperature
        )


def test_finite_composition_temperature_biases_species_relaxation():
    atoms = Atoms("CC", positions=[[0, 0, 0], [1.5, 0, 0]], cell=[10] * 3)
    initial = atoms.copy()
    initial[1].symbol = "Si"
    descriptor = DifferentiableFingerprint(["C", "Si"], body_orders=(2,))
    target = descriptor.generate(atoms)
    results = []
    for constraint in (None, 0.1):
        relaxation = FingerprintRelaxation(
            descriptor, target, composition_constraint=constraint,
            loss_weights=(0.0, 0.0, 0.0), fixed_species_mask=[True, False],
        )
        results.append(relaxation.run(
            initial, optimise_positions=False, max_steps=40, final_refine_steps=0,
            tolerance=0,
        ))
    assert results[0].structure.get_chemical_symbols() == ["C", "Si"]
    assert results[1].structure.get_chemical_symbols() == ["C", "C"]


def test_finite_composition_temperature_allows_fixed_species_mismatch():
    atoms = Atoms("CC", positions=[[0, 0, 0], [1.5, 0, 0]], cell=[10] * 3)
    initial = atoms.copy()
    initial[1].symbol = "Si"
    descriptor = DifferentiableFingerprint(["C", "Si"], body_orders=(2,))
    relaxation = FingerprintRelaxation(
        descriptor, descriptor.generate(atoms), composition_constraint=0.1,
        fixed_species_mask=[False, True], loss_weights=(1.0, 0.0, 0.0),
    )
    result = relaxation.run(initial, optimise_positions=False, max_steps=2, final_refine_steps=0)
    assert result.structure[1].symbol == "Si"
    assert result.loss_history[-1]["composition"] > 0


@pytest.mark.parametrize("constraint", [-1.0, -math.inf, math.nan, "invalid", 1e-320])
def test_invalid_composition_temperature_is_rejected(constraint):
    atoms = Atoms("C", positions=[[0, 0, 0]], cell=[10] * 3)
    descriptor = DifferentiableFingerprint(["C"], body_orders=(2,))
    with pytest.raises(ValueError, match="composition_constraint"):
        FingerprintRelaxation(descriptor, descriptor.generate(atoms), composition_constraint=constraint)


@pytest.mark.parametrize("constraint", [0.0, 1.0])
def test_composition_temperature_requires_target_counts(constraint):
    atoms = Atoms("C", positions=[[0, 0, 0]], cell=[10] * 3)
    descriptor = DifferentiableFingerprint(["C"], body_orders=(2,))
    target = descriptor.generate(atoms)
    target.composition = None
    with pytest.raises(ValueError, match="target_composition"):
        FingerprintRelaxation(descriptor, target, composition_constraint=constraint)
    FingerprintRelaxation(descriptor, target, composition_constraint=None)


def test_composition_temperature_rejects_gradient_overflow():
    atoms = Atoms("CC", positions=[[0, 0, 0], [1.5, 0, 0]], cell=[10] * 3)
    initial = atoms.copy()
    initial[1].symbol = "Si"
    descriptor = DifferentiableFingerprint(["C", "Si"], body_orders=(2,))
    relaxation = FingerprintRelaxation(
        descriptor, descriptor.generate(atoms), composition_constraint=1e-200,
        loss_weights=(1.0, 0.0, 0.0),
    )
    with pytest.raises(FloatingPointError, match="species gradient"):
        relaxation.run(initial, optimise_positions=False, max_steps=1, final_refine_steps=0)


def test_trajectory_and_callback_reflect_evaluated_species_and_positions():
    atoms = bulk("C", "diamond", a=3.567, cubic=False)
    initial = atoms.copy()
    initial[1].symbol = "Si"
    initial.positions[0, 0] += 0.08
    descriptor = DifferentiableFingerprint(["C", "Si"])
    target = descriptor.generate(atoms)
    relaxation = FingerprintRelaxation(
        descriptor,
        target,
        target_composition={"C": 2, "Si": 0},
        composition_constraint="exact",
    )
    callbacks = []
    result = relaxation.run(
        initial, max_steps=3, final_refine_steps=0,
        callback=lambda step, frame, probabilities, record: callbacks.append(
            (frame, probabilities, record)
        ),
    )

    assert result.trajectory[0][1].symbol == "Si"
    for (frame, probabilities, record), saved, saved_probabilities in zip(
        callbacks, result.trajectory[1:-1], result.species_probabilities_history[1:-1]
    ):
        assert frame.get_chemical_symbols() == ["C", "C"]
        assert saved.get_chemical_symbols() == frame.get_chemical_symbols()
        assert np.array_equal(saved.positions, frame.positions)
        assert torch.equal(probabilities, saved_probabilities)
        fingerprint = descriptor(frame.positions, frame.cell.array, probabilities, frame.pbc)
        loss = sum(
            (actual - expected).square().mean()
            for actual, expected in zip(
                (fingerprint.two_body, fingerprint.three_body, fingerprint.four_body),
                (target.two_body, target.three_body, target.four_body),
            )
        )
        assert float(loss) == pytest.approx(record["total"], abs=1e-12)


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


def test_joint_relaxation_does_not_stop_on_small_absolute_loss_changes():
    atoms = bulk("C", "diamond", a=3.567, cubic=True)
    descriptor = DifferentiableFingerprint(["C", "Si"])
    target = descriptor.generate(atoms)
    initial = atoms.copy()
    initial.positions[0, 0] += 1.2
    initial[1].symbol = "Si"
    relaxation = FingerprintRelaxation(
        descriptor,
        target,
        target_composition={"C": 8, "Si": 0},
        composition_constraint="exact",
        position_learning_rate=0.01,
    )

    result = relaxation.run(initial, max_steps=40, final_refine_steps=0)

    assert result.n_steps == 40
    assert result.loss < 0.8 * result.loss_history[0]["total"]
    assert np.linalg.norm(result.structure.positions[0] - initial.positions[0]) > 0.1


def test_radial_smoothing_gives_useful_long_range_gradient():
    atoms = bulk("C", "diamond", a=3.567, cubic=True)
    descriptor = DifferentiableFingerprint(["C", "Si"])
    target = descriptor.generate(atoms)
    positions = torch.tensor(atoms.positions, dtype=torch.float64)
    positions[0, 0] += 1.2
    positions.requires_grad_()
    fingerprint = descriptor(
        positions, atoms.cell.array, descriptor.probabilities_from_atoms(atoms), atoms.pbc
    )
    unsmoothed, _ = _fingerprint_loss(fingerprint, target, (1.0, 0.0, 0.0))
    smoothed, _ = _fingerprint_loss(fingerprint, target, (1.0, 0.0, 0.0), 0.5)
    unsmoothed_gradient = torch.autograd.grad(unsmoothed, positions, retain_graph=True)[0]
    smoothed_gradient = torch.autograd.grad(smoothed, positions)[0]

    assert unsmoothed_gradient[0, 0] < 0
    assert smoothed_gradient[0, 0] > 0


def test_angular_smoothing_and_density_losses_are_differentiable():
    atoms = bulk("C", "diamond", a=3.567, cubic=False)
    descriptor = DifferentiableFingerprint(["C"])
    target = descriptor.generate(atoms)
    positions = torch.tensor(atoms.positions, dtype=torch.float64)
    positions[0, 0] += 0.1
    positions.requires_grad_()
    fingerprint = descriptor(
        positions, atoms.cell.array, descriptor.probabilities_from_atoms(atoms), atoms.pbc
    )
    _, mean_losses = _fingerprint_loss(fingerprint, target)
    _, density_losses = _fingerprint_loss(
        fingerprint, target, loss_reduction="density"
    )
    for name in ("two_body", "three_body", "four_body"):
        grid = fingerprint.grids[name]
        assert float(density_losses[name].detach()) == pytest.approx(
            float(mean_losses[name].detach()) * grid.numel() / float(grid[1] - grid[0])
        )
    for weights in ((0.0, 1.0, 0.0), (0.0, 0.0, 1.0)):
        loss, _ = _fingerprint_loss(
            fingerprint, target, weights, angular_smoothing=0.25
        )
        gradient = torch.autograd.grad(loss, positions, retain_graph=True)[0]
        assert torch.isfinite(gradient).all()
        assert torch.linalg.vector_norm(gradient) > 0


def test_enabled_body_orders_and_multielement_neighbors():
    atoms = Atoms(
        "OOOO", positions=[[0, 0, 0], [1.2, 0, 0], [0, 1.2, 0], [0, 0, 1.2]],
        cell=[4.2] * 3, pbc=True,
    )
    mixed = DifferentiableFingerprint(["Ca", "O"]).generate(atoms)
    pure = DifferentiableFingerprint(["O"]).generate(atoms)
    assert torch.allclose(mixed.three_body[:, 1], pure.three_body[:, 0], atol=1e-10)
    assert torch.allclose(mixed.four_body[:, 1], pure.four_body[:, 0], atol=1e-10)
    assert pure.three_body.max() > 0 and pure.four_body.max() > 0

    limited = DifferentiableFingerprint(["Ca", "O"], body_orders=(2, 3)).generate(atoms)
    assert limited.body_orders == (2, 3)
    assert torch.allclose(limited.three_body, mixed.three_body)
    with pytest.raises(ValueError, match="4-body loss requires"):
        _fingerprint_loss(limited, limited)


def test_nonperiodic_four_body_includes_in_cell_neighbors():
    atoms = Atoms(
        "CCCC", positions=[[0, 0, 0], [1.3, 0, 0], [0, 1.3, 0], [0, 0, 2.7]],
        cell=[8.0] * 3, pbc=False,
    )
    descriptor = DifferentiableFingerprint(["C"])
    positions = torch.tensor(atoms.positions, dtype=torch.float64, requires_grad=True)
    fingerprint = descriptor(
        positions, atoms.cell.array, descriptor.probabilities_from_atoms(atoms), atoms.pbc
    )
    assert fingerprint.four_body.max() > 0
    gradient = torch.autograd.grad(fingerprint.four_body.square().sum(), positions)[0]
    assert torch.isfinite(gradient).all()
    assert torch.linalg.vector_norm(gradient) > 0


def test_sinkhorn_enforces_rows_and_composition_with_sharp_logits():
    torch.manual_seed(4)
    logits = (torch.randn((40, 4), dtype=torch.float64) * 1.5).requires_grad_()
    counts = torch.tensor([2.0, 6.0, 24.0, 8.0], dtype=torch.float64)
    probabilities = FingerprintRelaxation._sinkhorn(logits, counts, 0.5)
    assert torch.allclose(probabilities.sum(dim=1), torch.ones(40, dtype=torch.float64), atol=1e-8)
    assert torch.allclose(probabilities.sum(dim=0), counts, atol=1e-5)
    probabilities.square().sum().backward()
    assert torch.isfinite(logits.grad).all()


def test_refinement_reduces_loss_and_verbose_reports_steps(capsys):
    atoms = bulk("C", "diamond", a=3.567, cubic=False)
    descriptor = DifferentiableFingerprint(["C"])
    target = descriptor.generate(atoms)
    initial = atoms.copy()
    initial.positions[0, 0] += 0.18
    relaxation = FingerprintRelaxation(
        descriptor, target, fixed_position_mask=[False, True], position_learning_rate=0.01
    )
    before = relaxation.run(
        initial, optimise_species=False, max_steps=1, final_refine_steps=0
    )
    assert capsys.readouterr().out == ""
    after = relaxation.run(
        initial, optimise_species=False, max_steps=1, final_refine_steps=20, verbose=True
    )
    output = capsys.readouterr().out

    assert "Phase" in output and "Fmax" in output and "SpeciesGrad" in output
    assert "Joint" in output and "Refine" in output and "MaxStep" in output
    assert after.loss < before.loss * 0.2
    assert np.linalg.norm(after.structure.positions[0] - before.structure.positions[0]) > 0.01
    assert np.array_equal(after.structure.positions[1], initial.positions[1])


@pytest.mark.parametrize("multielement", [False, True])
def test_plot_comparison_script_writes_image(tmp_path, multielement):
    pytest.importorskip("matplotlib")
    script = (
        Path(__file__).resolve().parents[1]
        / "example/python_pkg/differentiable_fingerprint/plot_comparison.py"
    )
    output = tmp_path / "comparison.png"
    structure_args = []
    if multielement:
        structure = tmp_path / "structure.extxyz"
        write(structure, bulk("SiC", "zincblende", a=4.36, cubic=True))
        structure_args = [str(structure)]
    subprocess.run(
        [sys.executable, str(script), *structure_args, "--output", str(output)],
        check=True,
        capture_output=True,
        text=True,
        env={**os.environ, "MPLBACKEND": "Agg"},
    )
    assert output.read_bytes().startswith(b"\x89PNG\r\n\x1a\n")


def test_inverse_design_plots_final_and_target_fingerprints(tmp_path, monkeypatch):
    pytest.importorskip("matplotlib")
    example_dir = (
        Path(__file__).resolve().parents[1]
        / "example/python_pkg/differentiable_fingerprint"
    )
    monkeypatch.setenv("MPLBACKEND", "Agg")
    monkeypatch.syspath_prepend(str(example_dir))
    monkeypatch.chdir(tmp_path)
    import plot_comparison

    original_plot = plot_comparison.plot_fingerprints
    captured = {}

    def capture_plot(reference, actual, grids, elements, reference_label, actual_label, output):
        captured.update(reference=reference, actual=actual, labels=(reference_label, actual_label))
        original_plot(reference, actual, grids, elements, reference_label, actual_label, output)

    monkeypatch.setattr(plot_comparison, "plot_fingerprints", capture_plot)
    namespace = runpy.run_path(str(example_dir / "inverse_design.py"), run_name="__main__")

    assert (tmp_path / "diamond_fingerprint_comparison.png").read_bytes().startswith(
        b"\x89PNG\r\n\x1a\n"
    )
    assert captured["labels"] == ("Target", "Final")
    descriptor = DifferentiableFingerprint(["C", "Si"])
    target = descriptor.generate(bulk("C", "diamond", a=3.567, cubic=True))
    final = descriptor.generate(read(tmp_path / "diamond_inverse_design.extxyz"))
    for plotted, expected in zip(captured["reference"], (target.two_body, target.three_body, target.four_body)):
        np.testing.assert_allclose(plotted, expected.detach().numpy(), atol=1e-10)
    for plotted, expected in zip(captured["actual"], (final.two_body, final.three_body, final.four_body)):
        np.testing.assert_allclose(plotted, expected.detach().numpy(), atol=1e-7)
    result = namespace["result"]
    assert result.loss < 0.01 * result.loss_history[0]["total"]
    assert np.max(np.abs(captured["actual"][0] - captured["reference"][0])) < 1e-4
    trajectory = read(tmp_path / "diamond.traj", index=":")
    assert trajectory[0][1].symbol == "Si"
    assert trajectory[1][1].symbol == "C"
    assert np.linalg.norm(trajectory[-1].positions[0] - trajectory[0].positions[0]) > 0.5


@pytest.mark.parametrize(
    ("script_name", "arguments"),
    [
        ("diamond_to_carbon_phases.py", ["--steps", "2", "--fine-steps", "0"]),
        ("complex_species.py", ["--steps", "2"]),
        ("local_atom.py", ["--steps", "3"]),
    ],
)
def test_additional_examples_run(tmp_path, script_name, arguments):
    pytest.importorskip("matplotlib")
    example_dir = (
        Path(__file__).resolve().parents[1]
        / "example/python_pkg/differentiable_fingerprint"
    )
    completed = subprocess.run(
        [sys.executable, str(example_dir / script_name), *arguments, "--output-dir", str(tmp_path)],
        check=True, capture_output=True, text=True,
        env={**os.environ, "MPLBACKEND": "Agg", "OMP_NUM_THREADS": "1"},
    )
    if script_name == "diamond_to_carbon_phases.py":
        initial = bulk("C", "diamond", a=3.567, cubic=True)
        initial.positions += np.random.default_rng(42).normal(0.0, 0.15, initial.positions.shape)
        volume_per_atom = initial.get_volume() / len(initial)
        graphite = Atoms(
            "C4", scaled_positions=[[0, 0, 0], [1 / 3, 2 / 3, 0],
                                     [0, 0, 0.5], [2 / 3, 1 / 3, 0.5]],
            cell=[[2.46, 0, 0], [-1.23, np.sqrt(3) * 1.23, 0], [0, 0, 6.70]],
            pbc=True,
        ).repeat((2, 1, 1))
        targets = {
            "graphite": graphite,
            "fcc": bulk("C", "fcc", a=(4 * volume_per_atom) ** (1 / 3), cubic=True).repeat((2, 1, 1)),
            "bcc": bulk("C", "bcc", a=(2 * volume_per_atom) ** (1 / 3), cubic=True).repeat((2, 2, 1)),
        }
        descriptor = DifferentiableFingerprint(
            ["C"], parameters={"radius_distance_tol": [1.5, 2.8, 3.0, 6.0]}
        )
        initial_fp = descriptor.generate(initial)
        target_fingerprints = {}
        for phase, target_atoms in targets.items():
            frames = read(tmp_path / f"diamond_to_{phase}.traj", index=":")
            np.testing.assert_allclose(frames[0].positions, initial.positions)
            np.testing.assert_allclose(frames[0].cell.array, initial.cell.array)
            assert np.array_equal(frames[0].cell.array, frames[-1].cell.array)
            assert (tmp_path / f"diamond_to_{phase}.png").is_file()
            line = next(
                line for line in completed.stdout.splitlines()
                if line.startswith(f"diamond -> {phase}: fixed diamond cell, initial MSE = ")
            )
            reported = ast.literal_eval(line.split("initial MSE = ", 1)[1])
            target_fp = descriptor.generate(target_atoms)
            target_fingerprints[phase] = target_fp
            assert target_fp.three_body.max() > 0
            assert target_fp.four_body.max() > 0
            expected = [
                float((actual - target).square().mean())
                for actual, target in zip(
                    (initial_fp.two_body, initial_fp.three_body, initial_fp.four_body),
                    (target_fp.two_body, target_fp.three_body, target_fp.four_body),
                )
            ]
            np.testing.assert_allclose(reported, expected, rtol=1e-6, atol=1e-10)
        assert not torch.allclose(
            target_fingerprints["fcc"].three_body, target_fingerprints["bcc"].three_body
        )
    elif script_name == "complex_species.py":
        frames = read(tmp_path / "complex_species.traj", index=":")
        assert set(frames[0].get_chemical_symbols()) == {"Ca"}
        assert frames[-1].symbols.formula.count() == {"Ca": 2, "Cu": 6, "O": 24, "Ti": 8}
        assert np.array_equal(frames[0].positions, frames[-1].positions)
    else:
        frames = read(tmp_path / "local_atom.traj", index=":")
        assert np.array_equal(frames[0].positions[1:], frames[-1].positions[1:])
        assert not np.array_equal(frames[0].positions[0], frames[-1].positions[0])


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
    assert np.max(np.abs(reference[2] - actual.four_body.detach().numpy())) < 1e-4


@pytest.mark.parametrize(
    "atoms",
    [
        Atoms(
            "CCCC", positions=[[0, 0, 0], [1.3, 0, 0], [0, 1.3, 0], [0, 0, 2.7]],
            cell=[8.0] * 3, pbc=False,
        ),
        bulk("SiC", "zincblende", a=4.36, cubic=True),
    ],
)
def test_four_body_fortran_parity_for_in_cell_neighbors(atoms):
    from raffle.raffle import Raffle__Distribs_Container

    parameters = Raffle__Distribs_Container.distribs_container_type()
    descriptor = DifferentiableFingerprint(sorted(set(atoms.get_chemical_symbols())), parameters)
    expected = parameters.generate_fingerprint(atoms)[2]
    actual = descriptor.generate(atoms).four_body.detach().numpy()
    assert np.max(np.abs(expected - actual)) < 2e-6
