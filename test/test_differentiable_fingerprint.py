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
    assert namespace["result"].loss < 1e-8
    assert np.max(np.abs(captured["actual"][0] - captured["reference"][0])) < 1e-5
    trajectory = read(tmp_path / "diamond.traj", index=":")
    assert trajectory[0][1].symbol == "Si"
    assert trajectory[1][1].symbol == "C"
    assert np.linalg.norm(trajectory[-1].positions[0] - trajectory[0].positions[0]) > 0.5


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
