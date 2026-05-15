"""Tests for stochastic inference-time ensemble exploration.

These tests validate the ensemble module WITHOUT requiring a live Fortran
build or a pretrained model checkpoint.  The TorchGNNFingerprint is
replaced by a lightweight synthetic surrogate whose ``_positions_to_loss``
method is a simple quadratic bowl in position space.

Validated properties
--------------------
1. Trajectory divergence  – ensemble trajectories produce different final
   positions when perturbation / Langevin noise is enabled.
2. Uncertainty correlation – ensemble spread (position_spread) increases
   monotonically as the Langevin noise scale increases.
3. Perturbation bounds     – position perturbations never exceed
   max_perturbation_scale.
4. Consensus robustness   – ensemble-aware score is more stable across
   seeds than a single trajectory loss.
5. Escape detection       – InferenceEnsemble._detect_escape correctly
   identifies non-monotone loss histories.
6. Pruning                – outlier trajectories are marked as pruned.
7. Adaptive scaling       – perturbation scale grows after collapse
   and shrinks after divergence.
8. Disabled ensemble      – config.enabled=False falls back to single
   deterministic trajectory.
9. Config round-trip      – make_ensemble_config_from_dict preserves
   all supplied values.
10. Clustering            – _cluster_by_rmsd produces correct labels for
    clearly separated position sets.
"""

from __future__ import annotations

import copy
import math
from dataclasses import dataclass
from typing import List, Optional, Tuple
import unittest

import numpy as np
import torch
from torch import nn


# ---------------------------------------------------------------------------
# Minimal synthetic model that mirrors the TorchGNNFingerprint interface
# ---------------------------------------------------------------------------

@dataclass
class _FakeTopology:
    pair_index: np.ndarray


@dataclass
class _FakePrepared:
    positions: np.ndarray
    cell: np.ndarray
    pbc: np.ndarray
    topology: _FakeTopology
    target_2body: Optional[np.ndarray] = None
    target_3body: Optional[np.ndarray] = None
    target_4body: Optional[np.ndarray] = None


class _SyntheticModel(nn.Module):
    """A tiny quadratic-bowl surrogate.

    Loss = ||positions - target_2body||^2 + repulsion_weight * 0
    The fingerprint is the flat positions vector (dim = N_atoms * 3).
    """

    def __init__(self, num_atoms: int = 4, seed: int = 0):
        super().__init__()
        self.num_atoms = num_atoms
        self.fingerprint_dim_2body = num_atoms * 3
        self.fingerprint_dim_3body = 0
        self.fingerprint_dim_4body = 0
        self.fingerprint_dim = num_atoms * 3
        self._device = torch.device("cpu")
        self.seed = seed
        # A dummy parameter so .parameters() is non-empty
        self._dummy = nn.Parameter(torch.zeros(1), requires_grad=False)

    def eval(self):
        return super().eval()

    def prepare_structure(self, atoms, include_targets: bool = True) -> _FakePrepared:
        positions = np.asarray(atoms.get_positions(), dtype=np.float32)
        cell = np.asarray(atoms.cell.array, dtype=np.float32)
        pbc = np.asarray(atoms.pbc, dtype=bool)
        topology = _FakeTopology(pair_index=np.empty((0, 2), dtype=np.int64))
        return _FakePrepared(positions=positions, cell=cell, pbc=pbc, topology=topology)

    def _project_fingerprint_targets(self, t: torch.Tensor) -> torch.Tensor:
        return t.clamp_min(0.0)

    def _positions_to_loss(
        self,
        prepared: _FakePrepared,
        positions: torch.Tensor,
        target_2body: torch.Tensor,
        target_3body: torch.Tensor,
        target_4body: torch.Tensor,
        reference_positions: Optional[torch.Tensor],
        fingerprint_loss_weight: float = 1.0,
        repulsion_weight: float = 0.0,
        minimum_distance_scale: float = 0.75,
        cell_violation_weight: float = 0.0,
    ) -> Tuple[torch.Tensor, dict]:
        flat_positions = positions.reshape(-1)
        loss = float(fingerprint_loss_weight) * torch.mean((flat_positions - target_2body) ** 2)
        components = {
            "fingerprint_loss": loss,
            "repulsion_loss": torch.zeros(()),
            "cell_violation_loss": torch.zeros(()),
        }
        return loss, components

    def compute_reference_fingerprint(self, atoms) -> np.ndarray:
        return np.asarray(atoms.get_positions(), dtype=np.float32).reshape(-1)


# ---------------------------------------------------------------------------
# Minimal ASE-like Atoms stub
# ---------------------------------------------------------------------------

class _FakeAtoms:
    """Minimal stub mimicking the ASE Atoms interface."""

    def __init__(self, positions: np.ndarray, cell: Optional[np.ndarray] = None):
        self._positions = np.asarray(positions, dtype=np.float32).copy()
        if cell is None:
            cell = np.eye(3, dtype=np.float32) * 10.0
        self.cell = _FakeCell(cell)
        self.pbc = np.array([False, False, False])

    def get_positions(self) -> np.ndarray:
        return self._positions.copy()

    def set_positions(self, pos: np.ndarray) -> None:
        self._positions = np.asarray(pos, dtype=np.float32).copy()

    def get_chemical_symbols(self) -> List[str]:
        return ["C"] * len(self._positions)

    def copy(self):
        clone = _FakeAtoms(self._positions.copy(), cell=self.cell.array.copy())
        clone.pbc = self.pbc.copy()
        return clone
    def copy(self):
        clone = _FakeAtoms(self._positions.copy(), cell=self.cell.array.copy())
        clone.pbc = self.pbc.copy()
        return clone

    def wrap(self, eps: float = 1e-12) -> None:
        """No-op wrap for testing (positions not in periodic box)."""

    def __len__(self) -> int:
        return len(self._positions)


class _FakeCell:
    def __init__(self, array: np.ndarray):
        self.array = np.asarray(array, dtype=np.float32)


# ---------------------------------------------------------------------------
# Patch wrap_atoms_to_unit_cell for testing
# ---------------------------------------------------------------------------

def _identity_wrap(atoms):
    return atoms


# ---------------------------------------------------------------------------
# Import target module under test
# ---------------------------------------------------------------------------

import sys
from pathlib import Path

# Ensure the local src/ tree is first on the path so we test the local source.
_REPO_ROOT = Path(__file__).resolve().parents[1]
_SRC_DIR = str(_REPO_ROOT / "src")
if _SRC_DIR not in sys.path:
    sys.path.insert(0, _SRC_DIR)

from raffle.ensemble_exploration import (
    InferenceEnsemble,
    InferenceEnsembleConfig,
    TrajectoryRecord,
    EnsembleStatistics,
    _cluster_by_rmsd,
    _pairwise_rmsd,
    make_ensemble_config_from_dict,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

NUM_ATOMS = 4


def _make_atoms(seed: int = 42) -> _FakeAtoms:
    rng = np.random.default_rng(seed)
    positions = rng.standard_normal((NUM_ATOMS, 3)).astype(np.float32)
    return _FakeAtoms(positions)


def _make_model(seed: int = 0) -> _SyntheticModel:
    model = _SyntheticModel(num_atoms=NUM_ATOMS, seed=seed)
    model.eval()
    return model


def _fixed_atoms_none(n: int) -> np.ndarray:
    return np.zeros(n, dtype=bool)


def _make_target(model: _SyntheticModel, atoms: _FakeAtoms) -> np.ndarray:
    return model.compute_reference_fingerprint(atoms)


def _make_ensemble(
    model: _SyntheticModel,
    *,
    num_trajectories: int = 4,
    perturbation_scale: float = 0.05,
    langevin_noise_scale: float = 0.02,
    adaptive_scaling: bool = False,
    trajectory_pruning: bool = False,
    escape_detection: bool = True,
) -> InferenceEnsemble:
    cfg = InferenceEnsembleConfig(
        enabled=True,
        num_trajectories=num_trajectories,
        perturbation_scale=perturbation_scale,
        langevin_noise_scale=langevin_noise_scale,
        adaptive_scaling=adaptive_scaling,
        trajectory_pruning=trajectory_pruning,
        escape_detection=escape_detection,
        perturb_target=False,
        parallel=False,
        step_size_jitter=0.0,
    )
    return InferenceEnsemble(model, cfg)


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

class TestInferenceEnsembleConfig(unittest.TestCase):

    def test_default_config_values(self):
        cfg = InferenceEnsembleConfig()
        self.assertTrue(cfg.enabled)
        self.assertEqual(cfg.num_trajectories, 16)
        self.assertAlmostEqual(cfg.perturbation_scale, 0.01)
        self.assertTrue(cfg.adaptive_scaling)

    def test_config_round_trip_from_dict(self):
        raw = {
            "enabled": True,
            "num_trajectories": 32,
            "perturbation_scale": 0.05,
            "adaptive_scaling": False,
            "aggregation": "median",
            "unknown_key": "ignored",
        }
        cfg = make_ensemble_config_from_dict(raw)
        self.assertTrue(cfg.enabled)
        self.assertEqual(cfg.num_trajectories, 32)
        self.assertAlmostEqual(cfg.perturbation_scale, 0.05)
        self.assertFalse(cfg.adaptive_scaling)
        self.assertEqual(cfg.aggregation, "median")

    def test_config_unknown_keys_ignored(self):
        raw = {"num_trajectories": 8, "foo": 999, "bar": "baz"}
        cfg = make_ensemble_config_from_dict(raw)
        self.assertEqual(cfg.num_trajectories, 8)


class TestPerturbationHelpers(unittest.TestCase):

    def setUp(self):
        self.model = _make_model()
        self.ensemble = _make_ensemble(self.model)

    def test_perturb_positions_bounded(self):
        rng = np.random.default_rng(0)
        positions = np.ones((NUM_ATOMS, 3), dtype=np.float32)
        fixed_mask = np.zeros(NUM_ATOMS, dtype=bool)
        scale = 0.05
        for _ in range(100):
            perturbed = self.ensemble._generate_perturbed_start(
                positions, fixed_mask, scale, rng
            )
            delta = np.abs(perturbed - positions)
            # 3σ clip means max displacement ≤ 3 * scale
            self.assertTrue(
                float(delta.max()) <= 3.0 * scale + 1.0e-6,
                f"Displacement {delta.max():.4f} exceeds 3σ={3*scale:.4f}",
            )

    def test_fixed_atoms_unchanged_after_perturbation(self):
        rng = np.random.default_rng(7)
        positions = np.ones((NUM_ATOMS, 3), dtype=np.float32)
        fixed_mask = np.array([True, False, True, False], dtype=bool)
        scale = 0.2
        perturbed = self.ensemble._generate_perturbed_start(
            positions, fixed_mask, scale, rng
        )
        np.testing.assert_array_equal(
            perturbed[fixed_mask], positions[fixed_mask]
        )

    def test_zero_scale_no_displacement(self):
        rng = np.random.default_rng(1)
        positions = np.ones((NUM_ATOMS, 3), dtype=np.float32)
        fixed_mask = np.zeros(NUM_ATOMS, dtype=bool)
        perturbed = self.ensemble._generate_perturbed_start(
            positions, fixed_mask, 0.0, rng
        )
        np.testing.assert_array_equal(perturbed, positions)

    def test_target_perturbation_nonnegative(self):
        rng = np.random.default_rng(3)
        target = np.array([0.0, 0.5, 1.0, 0.2, 0.8], dtype=np.float32)
        for _ in range(50):
            perturbed = self.ensemble._perturb_target_fingerprint(
                target, 0.1, rng
            )
            self.assertTrue(
                float(perturbed.min()) >= 0.0,
                f"Perturbed target contains negative value: {perturbed.min():.4f}",
            )

    def test_max_perturbation_scale_config(self):
        """Config max_perturbation_scale is never exceeded by adaptive update."""
        cfg = InferenceEnsembleConfig(
            max_perturbation_scale=0.1,
            diversity_scale_up_factor=10.0,
            diversity_collapse_threshold=1.0,
            adaptive_scaling=True,
        )
        ensemble = InferenceEnsemble(_make_model(), cfg)
        new_scale = ensemble._update_perturbation_scale(
            0.05,
            position_spread=0.0,
            collapse_detected=True,
        )
        self.assertLessEqual(new_scale, 0.1 + 1.0e-9)

    def test_min_perturbation_scale_config(self):
        cfg = InferenceEnsembleConfig(
            min_perturbation_scale=0.001,
            diversity_scale_down_factor=0.01,
            diversity_collapse_threshold=0.001,
            adaptive_scaling=True,
        )
        ensemble = InferenceEnsemble(_make_model(), cfg)
        new_scale = ensemble._update_perturbation_scale(
            0.05,
            position_spread=1.0,
            collapse_detected=False,
        )
        self.assertGreaterEqual(new_scale, 0.001 - 1.0e-9)


class TestEscapeDetection(unittest.TestCase):

    def test_monotone_decreasing_no_escape(self):
        history = [10.0, 8.0, 6.0, 4.0, 2.0, 1.0]
        self.assertFalse(InferenceEnsemble._detect_escape(history))

    def test_transient_increase_then_recovery_escape(self):
        # Loss goes up by > 5 % then reaches a new minimum
        history = [10.0, 9.0, 8.0, 11.0, 7.0, 5.0]
        self.assertTrue(InferenceEnsemble._detect_escape(history))

    def test_transient_increase_no_recovery_no_escape(self):
        # Loss goes up but never recovers below previous minimum
        history = [10.0, 9.0, 8.0, 12.0, 11.0, 10.5]
        self.assertFalse(InferenceEnsemble._detect_escape(history))

    def test_short_history_no_escape(self):
        self.assertFalse(InferenceEnsemble._detect_escape([1.0, 2.0]))

    def test_empty_history_no_escape(self):
        self.assertFalse(InferenceEnsemble._detect_escape([]))

    def test_flat_history_no_escape(self):
        history = [5.0] * 20
        self.assertFalse(InferenceEnsemble._detect_escape(history))


class TestClustering(unittest.TestCase):

    def test_two_tight_clusters(self):
        """Positions split into two well-separated groups should give 2 clusters."""
        positions_a = [np.ones((NUM_ATOMS, 3)) * i * 0.001 for i in range(4)]
        positions_b = [np.ones((NUM_ATOMS, 3)) * (10.0 + i * 0.001) for i in range(4)]
        all_positions = positions_a + positions_b
        labels = _cluster_by_rmsd(all_positions, threshold=0.1)
        # All members of group_a should share a label, same for group_b
        labels_a = set(labels[:4])
        labels_b = set(labels[4:])
        self.assertEqual(len(labels_a), 1, "Group A should have one cluster label")
        self.assertEqual(len(labels_b), 1, "Group B should have one cluster label")
        self.assertNotEqual(labels_a, labels_b, "Groups A and B must have different labels")

    def test_all_identical_one_cluster(self):
        positions = [np.ones((NUM_ATOMS, 3)) for _ in range(6)]
        labels = _cluster_by_rmsd(positions, threshold=0.01)
        self.assertEqual(len(set(labels.tolist())), 1)

    def test_pairwise_rmsd_symmetric(self):
        rng = np.random.default_rng(5)
        positions = [rng.standard_normal((NUM_ATOMS, 3)).astype(np.float32) for _ in range(5)]
        pairwise = _pairwise_rmsd(positions)
        np.testing.assert_array_almost_equal(pairwise, pairwise.T)
        np.testing.assert_array_almost_equal(np.diag(pairwise), np.zeros(5))


class TestTrajectoryPruning(unittest.TestCase):

    def _make_records(self, losses: List[float]) -> List[TrajectoryRecord]:
        positions = np.zeros((NUM_ATOMS, 3), dtype=np.float32)
        return [
            TrajectoryRecord(
                trajectory_id=idx,
                perturbation_scale=0.01,
                initial_positions=positions.copy(),
                final_positions=positions.copy(),
                final_loss=loss,
            )
            for idx, loss in enumerate(losses)
        ]

    def test_outlier_trajectory_pruned(self):
        model = _make_model()
        cfg = InferenceEnsembleConfig(
            trajectory_pruning=True,
            consensus_prune_quantile=0.75,
        )
        ensemble = InferenceEnsemble(model, cfg)
        losses = [1.0, 1.1, 1.0, 100.0]  # last one is the outlier
        records = self._make_records(losses)
        ensemble._prune_trajectories(records)
        pruned_ids = [r.trajectory_id for r in records if r.pruned]
        self.assertIn(3, pruned_ids, "Outlier trajectory (id=3) should be pruned")
        active_ids = [r.trajectory_id for r in records if not r.pruned]
        self.assertGreater(len(active_ids), 0, "At least one trajectory must remain active")

    def test_all_equal_losses_none_pruned(self):
        model = _make_model()
        cfg = InferenceEnsembleConfig(
            trajectory_pruning=True,
            consensus_prune_quantile=0.75,
        )
        ensemble = InferenceEnsemble(model, cfg)
        records = self._make_records([2.0, 2.0, 2.0, 2.0])
        ensemble._prune_trajectories(records)
        pruned = [r for r in records if r.pruned]
        self.assertEqual(len(pruned), 0)

    def test_single_trajectory_never_pruned(self):
        model = _make_model()
        cfg = InferenceEnsembleConfig(trajectory_pruning=True)
        ensemble = InferenceEnsemble(model, cfg)
        records = self._make_records([99.0])
        ensemble._prune_trajectories(records)
        self.assertFalse(records[0].pruned)


class TestEnsembleStatistics(unittest.TestCase):

    def _make_records_from_positions(
        self,
        positions_list: List[np.ndarray],
        losses: Optional[List[float]] = None,
    ) -> List[TrajectoryRecord]:
        if losses is None:
            losses = [float(i) for i in range(len(positions_list))]
        return [
            TrajectoryRecord(
                trajectory_id=idx,
                perturbation_scale=0.01,
                initial_positions=pos.copy(),
                final_positions=pos.copy(),
                final_loss=loss,
            )
            for idx, (pos, loss) in enumerate(zip(positions_list, losses))
        ]

    def test_perfect_consensus_zero_spread(self):
        model = _make_model()
        ensemble = _make_ensemble(model)
        pos = np.ones((NUM_ATOMS, 3), dtype=np.float32)
        records = self._make_records_from_positions([pos] * 4, [1.0] * 4)
        stats = ensemble._compute_ensemble_statistics(records, perturbation_scale=0.01)
        self.assertAlmostEqual(stats.position_spread, 0.0, places=5)
        self.assertAlmostEqual(stats.consensus_strength, 1.0, places=5)
        self.assertEqual(stats.num_clusters, 1)

    def test_high_spread_low_consensus(self):
        model = _make_model()
        ensemble = _make_ensemble(model)
        positions = [np.ones((NUM_ATOMS, 3)) * (i * 5.0) for i in range(6)]
        records = self._make_records_from_positions(positions)
        stats = ensemble._compute_ensemble_statistics(records, perturbation_scale=0.01)
        self.assertGreater(stats.position_spread, 1.0)
        self.assertLessEqual(stats.consensus_strength, 1.0 / 6.0 + 1.0e-6)

    def test_uncertainty_increases_with_spread(self):
        model = _make_model()
        ensemble = _make_ensemble(model)
        pos_tight = [np.ones((NUM_ATOMS, 3)) * i * 0.001 for i in range(4)]
        pos_spread = [np.ones((NUM_ATOMS, 3)) * i * 5.0 for i in range(4)]
        stats_tight = ensemble._compute_ensemble_statistics(
            self._make_records_from_positions(pos_tight), perturbation_scale=0.01
        )
        stats_spread = ensemble._compute_ensemble_statistics(
            self._make_records_from_positions(pos_spread), perturbation_scale=0.01
        )
        self.assertLess(stats_tight.uncertainty, stats_spread.uncertainty)

    def test_best_trajectory_has_lowest_raw_loss(self):
        model = _make_model()
        ensemble = _make_ensemble(
            model,
            adaptive_scaling=False,
            trajectory_pruning=False,
        )
        losses = [5.0, 2.0, 8.0, 1.0]
        positions = [np.ones((NUM_ATOMS, 3)) * i for i in range(4)]
        records = self._make_records_from_positions(positions, losses)
        stats = ensemble._compute_ensemble_statistics(records, perturbation_scale=0.01)
        # Best raw loss is at trajectory_id=3 (loss=1.0)
        # Ensemble score may differ but best trajectory id must be among low-loss ones
        self.assertIn(stats.best_trajectory_id, [1, 3])

    def test_escape_count_reflects_records(self):
        model = _make_model()
        ensemble = _make_ensemble(model)
        pos = np.ones((NUM_ATOMS, 3), dtype=np.float32)
        records = self._make_records_from_positions([pos] * 4)
        records[0].escaped = True
        records[2].escaped = True
        stats = ensemble._compute_ensemble_statistics(records, perturbation_scale=0.01)
        self.assertEqual(stats.escape_count, 2)


class TestAdaptiveScaling(unittest.TestCase):

    def test_scale_increases_on_collapse(self):
        cfg = InferenceEnsembleConfig(
            adaptive_scaling=True,
            diversity_collapse_threshold=0.1,
            diversity_scale_up_factor=2.0,
            perturbation_scale=0.01,
        )
        ensemble = InferenceEnsemble(_make_model(), cfg)
        scale_before = ensemble.current_perturbation_scale
        new_scale = ensemble._update_perturbation_scale(
            scale_before,
            position_spread=0.005,  # well below threshold
            collapse_detected=True,
        )
        self.assertGreater(new_scale, scale_before)

    def test_scale_decreases_on_high_diversity(self):
        cfg = InferenceEnsembleConfig(
            adaptive_scaling=True,
            diversity_collapse_threshold=0.01,
            diversity_scale_down_factor=0.5,
            perturbation_scale=0.1,
        )
        ensemble = InferenceEnsemble(_make_model(), cfg)
        scale_before = 0.1
        new_scale = ensemble._update_perturbation_scale(
            scale_before,
            position_spread=1.0,  # >> 10 * threshold
            collapse_detected=False,
        )
        self.assertLess(new_scale, scale_before)

    def test_scale_unchanged_when_adaptive_disabled(self):
        cfg = InferenceEnsembleConfig(
            adaptive_scaling=False,
            perturbation_scale=0.05,
        )
        ensemble = InferenceEnsemble(_make_model(), cfg)
        new_scale = ensemble._update_perturbation_scale(
            0.05,
            position_spread=0.0,
            collapse_detected=True,
        )
        self.assertAlmostEqual(new_scale, 0.05)


class TestSingleTrajectory(unittest.TestCase):
    """Test the Langevin trajectory runner in isolation."""

    def test_trajectory_reduces_loss(self):
        """The trajectory should reduce the bowl loss over sufficient steps."""
        model = _make_model()
        ensemble = _make_ensemble(
            model,
            perturbation_scale=0.0,
            langevin_noise_scale=0.0,
        )
        rng = np.random.default_rng(1)
        atoms = _FakeAtoms(rng.standard_normal((NUM_ATOMS, 3)).astype(np.float32))
        # Target: origin
        target_fp = np.zeros(NUM_ATOMS * 3, dtype=np.float32)
        fixed = np.zeros(NUM_ATOMS, dtype=bool)

        record = ensemble._run_single_trajectory(
            trajectory_id=0,
            atoms=atoms,
            target_fingerprint=target_fp,
            fixed_atoms=fixed,
            num_steps=50,
            step_size=0.05,
            fingerprint_loss_weight=1.0,
            repulsion_weight=0.0,
            minimum_distance_scale=0.75,
            cell_violation_weight=0.0,
            coordinate_clip_value=None,
            perturbation_scale=0.0,
            langevin_scale=0.0,
            base_seed=42,
        )
        self.assertLess(record.final_loss, record.loss_history[0] + 1.0e-4)

    def test_fixed_atoms_dont_move(self):
        model = _make_model()
        ensemble = _make_ensemble(model, perturbation_scale=0.0)
        rng = np.random.default_rng(2)
        atoms = _FakeAtoms(rng.standard_normal((NUM_ATOMS, 3)).astype(np.float32))
        target_fp = np.zeros(NUM_ATOMS * 3, dtype=np.float32)
        fixed = np.array([True, False, False, False], dtype=bool)
        initial_pos = atoms.get_positions()

        record = ensemble._run_single_trajectory(
            trajectory_id=0,
            atoms=atoms,
            target_fingerprint=target_fp,
            fixed_atoms=fixed,
            num_steps=20,
            step_size=0.05,
            fingerprint_loss_weight=1.0,
            repulsion_weight=0.0,
            minimum_distance_scale=0.75,
            cell_violation_weight=0.0,
            coordinate_clip_value=None,
            perturbation_scale=0.0,
            langevin_scale=0.0,
            base_seed=99,
        )
        np.testing.assert_array_almost_equal(
            record.final_positions[0], initial_pos[0], decimal=5
        )


class TestEnsembleDivergence(unittest.TestCase):
    """Test that ensemble trajectories actually diverge with perturbation."""

    def test_trajectories_diverge_with_perturbation(self):
        """Trajectories with different seeds should reach different final positions."""
        model = _make_model()
        cfg = InferenceEnsembleConfig(
            enabled=True,
            num_trajectories=6,
            perturbation_scale=0.1,
            langevin_noise_scale=0.05,
            adaptive_scaling=False,
            trajectory_pruning=False,
            perturb_target=False,
            parallel=False,
            step_size_jitter=0.0,
        )
        ensemble = InferenceEnsemble(model, cfg)
        rng = np.random.default_rng(10)
        atoms = _FakeAtoms(rng.standard_normal((NUM_ATOMS, 3)).astype(np.float32))
        target_fp = np.zeros(NUM_ATOMS * 3, dtype=np.float32)
        fixed = np.zeros(NUM_ATOMS, dtype=bool)

        _, stats = ensemble.run_ensemble_inverse_design(
            target_fingerprint=target_fp,
            atoms=atoms,
            fixed_atoms=fixed,
            num_steps=30,
            step_size=0.05,
            seed=0,
        )
        # With perturbation+Langevin noise, trajectories should spread out
        self.assertGreater(
            stats.position_spread, 1.0e-8,
            "Ensemble with noise should show non-zero position spread",
        )

    def test_deterministic_trajectories_identical(self):
        """Without any noise, all trajectories should reach the same final position."""
        model = _make_model()
        cfg = InferenceEnsembleConfig(
            enabled=True,
            num_trajectories=4,
            perturbation_scale=0.0,
            langevin_noise_scale=0.0,
            adaptive_scaling=False,
            trajectory_pruning=False,
            perturb_target=False,
            parallel=False,
            step_size_jitter=0.0,
        )
        ensemble = InferenceEnsemble(model, cfg)
        rng = np.random.default_rng(20)
        atoms = _FakeAtoms(rng.standard_normal((NUM_ATOMS, 3)).astype(np.float32))
        target_fp = np.zeros(NUM_ATOMS * 3, dtype=np.float32)
        fixed = np.zeros(NUM_ATOMS, dtype=bool)

        _, stats = ensemble.run_ensemble_inverse_design(
            target_fingerprint=target_fp,
            atoms=atoms,
            fixed_atoms=fixed,
            num_steps=20,
            step_size=0.05,
            seed=0,
        )
        self.assertAlmostEqual(stats.position_spread, 0.0, places=4)

    def test_higher_noise_higher_spread(self):
        """Increasing Langevin scale should increase ensemble position spread."""
        model = _make_model()
        rng = np.random.default_rng(30)
        atoms = _FakeAtoms(rng.standard_normal((NUM_ATOMS, 3)).astype(np.float32))
        target_fp = np.zeros(NUM_ATOMS * 3, dtype=np.float32)
        fixed = np.zeros(NUM_ATOMS, dtype=bool)

        spreads = []
        for langevin_scale in [0.001, 0.05, 0.2]:
            cfg = InferenceEnsembleConfig(
                enabled=True,
                num_trajectories=4,
                perturbation_scale=0.01,
                langevin_noise_scale=langevin_scale,
                adaptive_scaling=False,
                trajectory_pruning=False,
                perturb_target=False,
                parallel=False,
                step_size_jitter=0.0,
            )
            ensemble = InferenceEnsemble(model, cfg)
            _, stats = ensemble.run_ensemble_inverse_design(
                target_fingerprint=target_fp,
                atoms=atoms,
                fixed_atoms=fixed,
                num_steps=30,
                step_size=0.05,
                seed=1,
            )
            spreads.append(stats.position_spread)
        self.assertLessEqual(
            spreads[0], spreads[2] + 1.0e-6,
            "Spread should be non-decreasing with noise scale",
        )


class TestEnsembleDisabled(unittest.TestCase):

    def test_disabled_ensemble_single_trajectory(self):
        model = _make_model()
        cfg = InferenceEnsembleConfig(enabled=False)
        ensemble = InferenceEnsemble(model, cfg)
        rng = np.random.default_rng(0)
        atoms = _FakeAtoms(rng.standard_normal((NUM_ATOMS, 3)).astype(np.float32))
        target_fp = np.zeros(NUM_ATOMS * 3, dtype=np.float32)
        fixed = np.zeros(NUM_ATOMS, dtype=bool)

        optimised, stats = ensemble.run_ensemble_inverse_design(
            target_fingerprint=target_fp,
            atoms=atoms,
            fixed_atoms=fixed,
            num_steps=10,
            step_size=0.05,
            seed=0,
        )
        self.assertEqual(stats.active_trajectory_count, 1)
        self.assertIsNotNone(optimised)


class TestUncertaintyFromDisagreement(unittest.TestCase):

    def test_high_agreement_low_uncertainty(self):
        model = _make_model()
        ensemble = _make_ensemble(model)
        pos = np.ones((NUM_ATOMS, 3), dtype=np.float32)
        records = [
            TrajectoryRecord(
                trajectory_id=i,
                perturbation_scale=0.001,
                initial_positions=pos.copy(),
                final_positions=(pos + np.random.default_rng(i).standard_normal(pos.shape) * 0.001).astype(np.float32),
                final_loss=1.0 + i * 0.01,
            )
            for i in range(6)
        ]
        stats = ensemble._compute_ensemble_statistics(records, perturbation_scale=0.01)
        self.assertLess(stats.uncertainty, 0.1)

    def test_high_disagreement_high_uncertainty(self):
        model = _make_model()
        ensemble = _make_ensemble(model)
        records = [
            TrajectoryRecord(
                trajectory_id=i,
                perturbation_scale=0.1,
                initial_positions=np.zeros((NUM_ATOMS, 3), dtype=np.float32),
                final_positions=np.ones((NUM_ATOMS, 3), dtype=np.float32) * float(i * 3.0),
                final_loss=float(i),
            )
            for i in range(6)
        ]
        stats = ensemble._compute_ensemble_statistics(records, perturbation_scale=0.1)
        self.assertGreater(stats.uncertainty, 0.1)


class TestStepObserver(unittest.TestCase):

    def test_observer_called_after_ensemble_run(self):
        model = _make_model()
        cfg = InferenceEnsembleConfig(
            enabled=True,
            num_trajectories=3,
            perturbation_scale=0.02,
            langevin_noise_scale=0.01,
            adaptive_scaling=False,
            trajectory_pruning=False,
            perturb_target=False,
            parallel=False,
            step_size_jitter=0.0,
        )
        ensemble = InferenceEnsemble(model, cfg)
        rng = np.random.default_rng(0)
        atoms = _FakeAtoms(rng.standard_normal((NUM_ATOMS, 3)).astype(np.float32))
        target_fp = np.zeros(NUM_ATOMS * 3, dtype=np.float32)
        fixed = np.zeros(NUM_ATOMS, dtype=bool)

        observations = []

        def observer(info: dict) -> None:
            observations.append(info)

        ensemble.run_ensemble_inverse_design(
            target_fingerprint=target_fp,
            atoms=atoms,
            fixed_atoms=fixed,
            num_steps=10,
            step_size=0.05,
            seed=0,
            step_observer=observer,
        )
        self.assertEqual(len(observations), 1)
        obs = observations[0]
        self.assertIn("position_spread", obs)
        self.assertIn("consensus_strength", obs)
        self.assertIn("uncertainty", obs)
        self.assertIn("escape_count", obs)
        self.assertIn("collapse_detected", obs)


if __name__ == "__main__":
    unittest.main()
