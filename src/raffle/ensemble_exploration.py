"""Stochastic inference-time ensemble exploration for inverse design.

This module implements a "murmuration of thoughts" ensemble strategy:
multiple stochastic trajectories run simultaneously from the same frozen
model weights, with statistical aggregation to estimate uncertainty and
guide exploration of the structural search space.

The trained model weights remain completely frozen throughout.  All
stochasticity is injected at the inference/search level only via:

* Atom-displacement jitter on the starting structure.
* Small perturbations to the target fingerprint descriptor.
* Lattice parameter jitter (optional).
* Langevin noise added at every optimisation step (stochastic gradient
  dynamics with temperature-controlled exploration scale).
* Per-trajectory variation of hyperparameters (step size, noise scale).

Reference SDE update form:
    x_(t+1) = x_t - lr * grad_L(x_t) + sqrt(2 * lr * T) * eta

where eta ~ N(0, I), T is the temperature, and lr is the step size.
"""

from __future__ import annotations

import math
import copy
from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np
import torch
from torch import nn

from .structure_metrics import wrap_atoms_to_unit_cell


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _as_float_tensor(array: np.ndarray, device: torch.device) -> torch.Tensor:
    return torch.as_tensor(np.asarray(array, dtype=np.float32), dtype=torch.float32, device=device)


def _pairwise_rmsd(positions_list: List[np.ndarray]) -> np.ndarray:
    """Compute pairwise RMSDs between a list of position arrays.

    Parameters
    ----------
    positions_list:
        List of (N_atoms, 3) position arrays.

    Returns
    -------
    np.ndarray of shape (n, n) containing pairwise RMSD values.
    """
    n = len(positions_list)
    matrix = np.zeros((n, n), dtype=np.float64)
    for i in range(n):
        for j in range(i + 1, n):
            delta = positions_list[i] - positions_list[j]
            rmsd = float(np.sqrt(np.mean(delta ** 2)))
            matrix[i, j] = rmsd
            matrix[j, i] = rmsd
    return matrix


def _cluster_by_rmsd(
    positions_list: List[np.ndarray],
    threshold: float,
) -> np.ndarray:
    """Single-linkage clustering based on pairwise RMSD.

    Returns integer cluster label array of length n.
    """
    n = len(positions_list)
    labels = np.arange(n, dtype=np.int64)
    pairwise = _pairwise_rmsd(positions_list)
    for i in range(n):
        for j in range(i + 1, n):
            if pairwise[i, j] < float(threshold):
                old_label = labels[j]
                new_label = labels[i]
                labels[labels == old_label] = new_label
    # Remap to 0..K-1
    unique = np.unique(labels)
    remap = {old: new for new, old in enumerate(unique)}
    return np.array([remap[label] for label in labels], dtype=np.int64)


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

@dataclass
class InferenceEnsembleConfig:
    """Configuration for stochastic inference-time ensemble exploration.

    All perturbations act *only* during inference; no model weights are
    modified.

    Parameters
    ----------
    enabled:
        Master switch.  When False, ``InferenceEnsemble.run_ensemble_inverse_design``
        falls back to a single deterministic trajectory.
    num_trajectories:
        Number of parallel stochastic trajectories to maintain.
    perturbation_scale:
        Initial standard deviation of atom-displacement noise applied to
        the starting structure (in Angstroms).
    adaptive_scaling:
        When True, the perturbation scale is adjusted online based on
        trajectory diversity.
    aggregation:
        Strategy used to derive the *reported* fingerprint from the
        ensemble.  One of ``"mean_variance"``, ``"median"``,
        ``"cluster"``.
    consensus_metric:
        How to compute consensus strength.  One of ``"cluster"``,
        ``"variance"``, ``"entropy"``.
    escape_detection:
        Whether to detect and record local-minimum escape events.
    trajectory_pruning:
        Dynamically prune trajectories that are persistent outliers.
    parallel:
        When True, trajectories are evaluated using a thread pool.
        **Note**: PyTorch CPU forward/backward through shared frozen
        parameters is generally thread-safe as long as no optimiser
        updates model parameters.  Set to False if you observe race
        conditions on exotic hardware.
    perturb_positions:
        Add Gaussian displacement noise to starting atom positions.
    perturb_lattice:
        Add fractional strain jitter to the unit-cell matrix.
    perturb_target:
        Add small Gaussian noise to the target fingerprint descriptor
        per trajectory, simulating descriptor uncertainty.
    target_perturbation_scale:
        Standard deviation of per-trajectory target fingerprint noise,
        expressed as a fraction of the target fingerprint RMS norm.
    langevin_noise_scale:
        Base standard deviation of the Langevin noise injected at each
        optimisation step (Angstroms).  The effective scale is
        ``langevin_noise_scale * sqrt(temperature)``.
    temperature:
        Temperature parameter for Langevin dynamics.  Higher values
        increase exploration; lower values make trajectories more
        deterministic.
    min_perturbation_scale / max_perturbation_scale:
        Bounds on adaptive perturbation scaling.
    diversity_collapse_threshold:
        Position spread (Angstroms) below which the ensemble is
        considered collapsed and the perturbation scale is increased.
    diversity_scale_up_factor / diversity_scale_down_factor:
        Multiplicative factors for adaptive scale updates.
    consensus_prune_quantile:
        Trajectories with final loss above this quantile of the
        ensemble distribution are pruned (if ``trajectory_pruning``
        is True).
    variance_weight / consensus_weight / escape_weight:
        Coefficients in the ensemble-aware acquisition function::

            score = mean_loss
                    + variance_weight * loss_variance
                    - consensus_weight * consensus_strength
                    - escape_weight * escape_rate
    cluster_distance_threshold:
        RMSD threshold (Angstroms) used for consensus clustering.
    step_size_jitter:
        Relative jitter on the per-trajectory step size
        (uniform in ``[1 - jitter, 1 + jitter]``).
    """

    enabled: bool = True
    num_trajectories: int = 16
    perturbation_scale: float = 0.01
    adaptive_scaling: bool = True
    aggregation: str = "mean_variance"
    consensus_metric: str = "cluster"
    escape_detection: bool = True
    trajectory_pruning: bool = True
    parallel: bool = False

    # Perturbation channels
    perturb_positions: bool = True
    perturb_lattice: bool = False
    perturb_target: bool = True
    target_perturbation_scale: float = 0.005

    # Langevin dynamics
    langevin_noise_scale: float = 0.005
    temperature: float = 1.0

    # Adaptive scaling bounds
    min_perturbation_scale: float = 1.0e-4
    max_perturbation_scale: float = 0.3

    # Collapse / divergence thresholds
    diversity_collapse_threshold: float = 0.02
    diversity_scale_up_factor: float = 1.5
    diversity_scale_down_factor: float = 0.85

    # Pruning
    consensus_prune_quantile: float = 0.75

    # Acquisition weights
    variance_weight: float = 0.05
    consensus_weight: float = 0.1
    escape_weight: float = 0.05

    # Clustering
    cluster_distance_threshold: float = 0.1

    # Per-trajectory hyperparameter jitter
    step_size_jitter: float = 0.1


# ---------------------------------------------------------------------------
# Data containers
# ---------------------------------------------------------------------------

@dataclass
class TrajectoryRecord:
    """Full history of a single stochastic trajectory.

    Attributes
    ----------
    trajectory_id:
        Integer index in the ensemble (0-based).
    perturbation_scale:
        The position perturbation scale actually used at initialisation.
    initial_positions:
        Atom positions at the start of this trajectory, *after*
        perturbation, shape (N_atoms, 3).
    final_positions:
        Atom positions at the end of the trajectory, shape (N_atoms, 3).
    loss_history:
        Total loss value at each optimisation step.
    final_loss:
        Total loss at the last step.
    escaped:
        True if the trajectory escaped at least one local minimum (loss
        increased transiently before reaching a better basin).
    converged:
        True if the loss improvement in the final 20 % of steps is below
        1 % of the initial loss.
    fingerprint_history:
        (Sparse) list of concatenated fingerprint vectors sampled during
        the trajectory.  Not stored by default (empty list).
    pruned:
        True if this trajectory was pruned by the ensemble manager.
    step_size:
        The actual per-trajectory step size used.
    langevin_scale:
        The effective Langevin noise scale for this trajectory.
    perturbation_seed:
        RNG seed used for reproducibility.
    """

    trajectory_id: int
    perturbation_scale: float
    initial_positions: np.ndarray
    final_positions: np.ndarray
    loss_history: List[float] = field(default_factory=list)
    final_loss: float = float("inf")
    escaped: bool = False
    converged: bool = False
    fingerprint_history: List[np.ndarray] = field(default_factory=list)
    pruned: bool = False
    step_size: float = 0.0
    langevin_scale: float = 0.0
    perturbation_seed: int = 0


@dataclass
class EnsembleStatistics:
    """Aggregated statistics derived from a completed ensemble run.

    Attributes
    ----------
    mean_positions:
        Element-wise mean of final positions across active trajectories,
        shape (N_atoms, 3).
    median_positions:
        Element-wise median of final positions, shape (N_atoms, 3).
    position_variance:
        Element-wise variance of final positions, shape (N_atoms, 3).
    position_spread:
        Scalar diversity measure: mean RMSD between all pairs of
        trajectory endpoints (Angstroms).
    mean_loss / median_loss / loss_variance:
        Distribution of final trajectory losses.
    entropy:
        Shannon entropy of the normalised loss distribution (nats).
    pairwise_divergence:
        (n_active, n_active) matrix of pairwise endpoint RMSDs.
    cluster_labels:
        Integer cluster label per active trajectory.
    num_clusters:
        Number of distinct convergence clusters.
    consensus_strength:
        Fraction of trajectories in the largest cluster; 1 = perfect
        consensus, 0 = maximum diversity.
    uncertainty:
        Operational uncertainty estimate derived from ensemble spread.
        ``uncertainty = position_spread / (1 + consensus_strength)``.
    best_trajectory_id:
        Index of the trajectory with the lowest ensemble-aware score.
    best_positions:
        Atom positions from the best trajectory, shape (N_atoms, 3).
    best_loss:
        Raw final loss of the best trajectory.
    best_ensemble_score:
        Ensemble-aware score of the best trajectory (lower is better).
    escape_count:
        Number of trajectories that escaped a local minimum.
    collapse_detected:
        True if trajectory diversity fell below the collapse threshold.
    active_trajectory_count:
        Number of non-pruned trajectories.
    current_perturbation_scale:
        Perturbation scale *after* adaptive update.
    trajectory_records:
        The full list of ``TrajectoryRecord`` objects for further analysis.
    """

    mean_positions: np.ndarray
    median_positions: np.ndarray
    position_variance: np.ndarray
    position_spread: float
    mean_loss: float
    median_loss: float
    loss_variance: float
    entropy: float
    pairwise_divergence: np.ndarray
    cluster_labels: np.ndarray
    num_clusters: int
    consensus_strength: float
    uncertainty: float
    best_trajectory_id: int
    best_positions: np.ndarray
    best_loss: float
    best_ensemble_score: float
    escape_count: int
    collapse_detected: bool
    active_trajectory_count: int
    current_perturbation_scale: float
    trajectory_records: List[TrajectoryRecord] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Main class
# ---------------------------------------------------------------------------

class InferenceEnsemble:
    """Stochastic inference-time ensemble explorer for inverse design.

    Wraps a **frozen** ``TorchGNNFingerprint`` surrogate model and runs
    multiple stochastic Langevin trajectories simultaneously to explore
    the structural search space more robustly than a single deterministic
    optimisation path.

    Parameters
    ----------
    model:
        A fitted ``TorchGNNFingerprint`` instance.  **Weights are never
        modified.**
    config:
        ``InferenceEnsembleConfig`` controlling all stochastic parameters.
    """

    def __init__(self, model, config: Optional[InferenceEnsembleConfig] = None) -> None:
        if config is None:
            config = InferenceEnsembleConfig()
        self._model = model
        self._config = config
        self._current_perturbation_scale = float(config.perturbation_scale)
        # Ensure the model is in eval mode and gradients are not tracked for params.
        self._model.eval()

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    @property
    def config(self) -> InferenceEnsembleConfig:
        return self._config

    @property
    def current_perturbation_scale(self) -> float:
        return self._current_perturbation_scale

    def run_ensemble_inverse_design(
        self,
        target_fingerprint: np.ndarray,
        atoms,
        fixed_atoms: np.ndarray,
        *,
        num_steps: int = 200,
        step_size: float = 1.0e-2,
        fingerprint_loss_weight: float = 1.0,
        repulsion_weight: float = 10.0,
        minimum_distance_scale: float = 0.75,
        cell_violation_weight: float = 0.0,
        coordinate_clip_value: Optional[float] = None,
        seed: int = 0,
        step_observer: Optional[Callable[[dict], None]] = None,
    ) -> Tuple[object, EnsembleStatistics]:
        """Run ensemble inverse design and return the best structure.

        Parameters
        ----------
        target_fingerprint:
            1-D target descriptor array of length
            ``model.fingerprint_dim``.
        atoms:
            Starting structure (ASE ``Atoms``).
        fixed_atoms:
            Boolean mask of length ``len(atoms)``; True = atom is fixed.
        num_steps:
            Number of Langevin optimisation steps per trajectory.
        step_size:
            Base Adam step size (learning rate).
        fingerprint_loss_weight:
            Weight on the fingerprint MSE loss term.
        repulsion_weight:
            Weight on the soft pairwise repulsion penalty.
        minimum_distance_scale:
            Minimum allowed distance as fraction of sum of covalent radii.
        cell_violation_weight:
            Weight on the periodic-boundary violation penalty.
        coordinate_clip_value:
            If set, clip atom displacements from initial positions.
        seed:
            Base RNG seed (each trajectory gets ``seed + trajectory_id``).
        step_observer:
            Optional callback called at the end of each *ensemble step*
            with a summary dict.

        Returns
        -------
        optimised_atoms:
            The best structure found across all trajectories.
        stats:
            ``EnsembleStatistics`` containing ensemble-level diagnostics.
        """
        cfg = self._config
        if not cfg.enabled:
            # Fallback: single deterministic trajectory
            record = self._run_single_trajectory(
                trajectory_id=0,
                atoms=atoms,
                target_fingerprint=target_fingerprint,
                fixed_atoms=fixed_atoms,
                num_steps=num_steps,
                step_size=step_size,
                fingerprint_loss_weight=fingerprint_loss_weight,
                repulsion_weight=repulsion_weight,
                minimum_distance_scale=minimum_distance_scale,
                cell_violation_weight=cell_violation_weight,
                coordinate_clip_value=coordinate_clip_value,
                perturbation_scale=0.0,
                langevin_scale=0.0,
                base_seed=seed,
            )
            stats = self._compute_ensemble_statistics(
                [record],
                perturbation_scale=0.0,
            )
            optimised = atoms.copy()
            optimised.set_positions(record.final_positions)
            return wrap_atoms_to_unit_cell(optimised), stats

        num_traj = max(int(cfg.num_trajectories), 1)
        records: List[TrajectoryRecord] = []

        if cfg.parallel:
            from concurrent.futures import ThreadPoolExecutor
            futures = {}
            with ThreadPoolExecutor(max_workers=num_traj) as pool:
                for traj_id in range(num_traj):
                    fut = pool.submit(
                        self._run_single_trajectory,
                        trajectory_id=traj_id,
                        atoms=atoms,
                        target_fingerprint=target_fingerprint,
                        fixed_atoms=fixed_atoms,
                        num_steps=num_steps,
                        step_size=step_size,
                        fingerprint_loss_weight=fingerprint_loss_weight,
                        repulsion_weight=repulsion_weight,
                        minimum_distance_scale=minimum_distance_scale,
                        cell_violation_weight=cell_violation_weight,
                        coordinate_clip_value=coordinate_clip_value,
                        perturbation_scale=self._current_perturbation_scale,
                        langevin_scale=float(cfg.langevin_noise_scale),
                        base_seed=seed + traj_id,
                    )
                    futures[fut] = traj_id
                from concurrent.futures import as_completed
                for fut in as_completed(futures):
                    records.append(fut.result())
        else:
            for traj_id in range(num_traj):
                record = self._run_single_trajectory(
                    trajectory_id=traj_id,
                    atoms=atoms,
                    target_fingerprint=target_fingerprint,
                    fixed_atoms=fixed_atoms,
                    num_steps=num_steps,
                    step_size=step_size,
                    fingerprint_loss_weight=fingerprint_loss_weight,
                    repulsion_weight=repulsion_weight,
                    minimum_distance_scale=minimum_distance_scale,
                    cell_violation_weight=cell_violation_weight,
                    coordinate_clip_value=coordinate_clip_value,
                    perturbation_scale=self._current_perturbation_scale,
                    langevin_scale=float(cfg.langevin_noise_scale),
                    base_seed=seed + traj_id,
                )
                records.append(record)

        if cfg.trajectory_pruning:
            self._prune_trajectories(records)

        stats = self._compute_ensemble_statistics(
            records,
            perturbation_scale=self._current_perturbation_scale,
        )

        if cfg.adaptive_scaling:
            self._current_perturbation_scale = self._update_perturbation_scale(
                self._current_perturbation_scale,
                stats,
            )

        if step_observer is not None:
            step_observer(
                {
                    "ensemble_statistics": stats,
                    "active_trajectories": stats.active_trajectory_count,
                    "position_spread": stats.position_spread,
                    "consensus_strength": stats.consensus_strength,
                    "uncertainty": stats.uncertainty,
                    "best_loss": stats.best_loss,
                    "escape_count": stats.escape_count,
                    "collapse_detected": stats.collapse_detected,
                }
            )

        optimised = atoms.copy()
        optimised.set_positions(stats.best_positions)
        return wrap_atoms_to_unit_cell(optimised), stats

    # ------------------------------------------------------------------
    # Perturbation helpers
    # ------------------------------------------------------------------

    def _generate_perturbed_start(
        self,
        positions_initial: np.ndarray,
        fixed_mask: np.ndarray,
        scale: float,
        rng: np.random.Generator,
    ) -> np.ndarray:
        """Return position array with small Gaussian noise on movable atoms."""
        if scale <= 0.0:
            return positions_initial.copy()
        noise = rng.standard_normal(positions_initial.shape).astype(np.float32)
        # Clamp noise to ±3σ to avoid catastrophic displacements
        noise = np.clip(noise, -3.0, 3.0)
        perturbed = positions_initial + scale * noise
        if fixed_mask is not None:
            perturbed[fixed_mask] = positions_initial[fixed_mask]
        return perturbed.astype(np.float32)

    def _perturb_lattice(
        self,
        cell: np.ndarray,
        scale: float,
        rng: np.random.Generator,
    ) -> np.ndarray:
        """Apply small isotropic strain jitter to cell matrix."""
        if scale <= 0.0:
            return cell.copy()
        strain = 1.0 + scale * np.clip(rng.standard_normal(), -2.0, 2.0)
        return (cell * strain).astype(np.float32)

    def _perturb_target_fingerprint(
        self,
        target: np.ndarray,
        scale: float,
        rng: np.random.Generator,
    ) -> np.ndarray:
        """Add small descriptor noise to the target fingerprint."""
        if scale <= 0.0:
            return target.copy()
        target_rms = float(np.sqrt(np.mean(target ** 2)))
        abs_scale = scale * max(target_rms, 1.0e-6)
        noise = abs_scale * rng.standard_normal(target.shape).astype(np.float32)
        perturbed = (target + noise).clip(min=0.0)
        return perturbed.astype(np.float32)

    # ------------------------------------------------------------------
    # Single-trajectory Langevin optimisation
    # ------------------------------------------------------------------

    def _run_single_trajectory(
        self,
        *,
        trajectory_id: int,
        atoms,
        target_fingerprint: np.ndarray,
        fixed_atoms: np.ndarray,
        num_steps: int,
        step_size: float,
        fingerprint_loss_weight: float,
        repulsion_weight: float,
        minimum_distance_scale: float,
        cell_violation_weight: float,
        coordinate_clip_value: Optional[float],
        perturbation_scale: float,
        langevin_scale: float,
        base_seed: int,
    ) -> TrajectoryRecord:
        """Run one Langevin trajectory and return its record.

        This is the core stochastic dynamics engine.  Mirrors the Adam
        optimisation loop inside ``TorchGNNFingerprint.inverse_design``
        but injects Langevin noise at every step and supports per-
        trajectory hyperparameter variation.

        The **model weights are never modified** — we only optimise
        the position tensor.
        """
        cfg = self._config
        rng = np.random.default_rng(base_seed)
        torch_rng = torch.Generator(device="cpu")
        torch_rng.manual_seed(int(base_seed))

        # --- Per-trajectory step-size jitter ---
        jitter = 1.0
        if cfg.step_size_jitter > 0.0:
            jitter = 1.0 + cfg.step_size_jitter * float(rng.uniform(-1.0, 1.0))
        traj_step_size = float(step_size) * jitter

        # --- Effective Langevin noise ---
        traj_langevin = float(langevin_scale) * math.sqrt(max(float(cfg.temperature), 0.0))

        model = self._model
        device = model._device

        # Prepare structure topology (does not invoke Fortran backend)
        prepared = model.prepare_structure(atoms, include_targets=False)

        fixed_mask_np = np.asarray(fixed_atoms, dtype=bool)
        fixed_mask = torch.as_tensor(fixed_mask_np, dtype=torch.bool, device=device)
        movable_mask = ~fixed_mask

        # --- Initial position perturbation ---
        positions_np = np.asarray(atoms.get_positions(), dtype=np.float32)
        if cfg.perturb_positions and perturbation_scale > 0.0:
            positions_np = self._generate_perturbed_start(
                positions_np, fixed_mask_np, perturbation_scale, rng
            )
        if cfg.perturb_lattice and perturbation_scale > 0.0:
            prepared = self._apply_lattice_perturbation(prepared, perturbation_scale, rng)

        positions_initial = _as_float_tensor(positions_np, device)

        # --- Target fingerprint (optionally perturbed per trajectory) ---
        target_np = np.asarray(target_fingerprint, dtype=np.float32)
        if cfg.perturb_target and cfg.target_perturbation_scale > 0.0:
            target_np = self._perturb_target_fingerprint(
                target_np, cfg.target_perturbation_scale, rng
            )
        target = model._project_fingerprint_targets(
            _as_float_tensor(target_np, device)
        )
        target_2body = target[:model.fingerprint_dim_2body]
        offset = model.fingerprint_dim_2body
        target_3body = target[offset:offset + model.fingerprint_dim_3body]
        offset += model.fingerprint_dim_3body
        target_4body = target[offset:offset + model.fingerprint_dim_4body]

        # --- Optimisation loop with Langevin noise ---
        positions_param = nn.Parameter(positions_initial.clone())
        optimiser = torch.optim.Adam([positions_param], lr=traj_step_size)

        loss_history: List[float] = []
        initial_loss: float = float("inf")

        for step in range(int(num_steps)):
            optimiser.zero_grad()
            # Fixed atoms: use initial positions
            candidate_positions = torch.where(
                fixed_mask.unsqueeze(-1),
                positions_initial,
                positions_param,
            )
            total_loss, _ = model._positions_to_loss(
                prepared,
                candidate_positions,
                target_2body,
                target_3body,
                target_4body,
                reference_positions=None,
                fingerprint_loss_weight=fingerprint_loss_weight,
                repulsion_weight=repulsion_weight,
                minimum_distance_scale=minimum_distance_scale,
                cell_violation_weight=cell_violation_weight,
            )
            total_loss.backward()

            # Zero gradients on model parameters to avoid accumulation
            for p in model.parameters():
                if p.grad is not None:
                    p.grad.detach_()
                    p.grad.zero_()

            # Zero fixed-atom gradients
            if positions_param.grad is not None:
                positions_param.grad[fixed_mask] = 0.0

            torch.nn.utils.clip_grad_value_([positions_param], 1.0e-1)
            optimiser.step()

            # --- Langevin noise injection (SDE update) ---
            if traj_langevin > 0.0 and bool(movable_mask.any()):
                with torch.no_grad():
                    noise = torch.randn(
                        positions_param.data.shape,
                        generator=torch_rng,
                        dtype=positions_param.dtype,
                        device=device,
                    )
                    langevin_dt = math.sqrt(2.0 * traj_step_size)
                    positions_param.data[movable_mask] += (
                        traj_langevin * langevin_dt * noise[movable_mask]
                    )

            # Enforce fixed atoms and optional clip
            with torch.no_grad():
                positions_param.data[fixed_mask] = positions_initial[fixed_mask]
                if coordinate_clip_value is not None and bool(movable_mask.any()):
                    max_d = float(coordinate_clip_value)
                    delta = positions_param.data[movable_mask] - positions_initial[movable_mask]
                    positions_param.data[movable_mask] = (
                        positions_initial[movable_mask] + delta.clamp(-max_d, max_d)
                    )

            step_loss = float(total_loss.item())
            loss_history.append(step_loss)
            if step == 0:
                initial_loss = step_loss

        with torch.no_grad():
            final_positions_t = torch.where(
                fixed_mask.unsqueeze(-1),
                positions_initial,
                positions_param,
            )
            final_positions_np = final_positions_t.detach().cpu().numpy().astype(np.float32)
        final_loss = loss_history[-1] if loss_history else float("inf")

        # --- Escape detection ---
        escaped = False
        if cfg.escape_detection and len(loss_history) >= 4:
            escaped = self._detect_escape(loss_history)

        # --- Convergence detection ---
        converged = False
        if len(loss_history) >= 10:
            tail_length = max(len(loss_history) // 5, 2)
            tail = loss_history[-tail_length:]
            tail_improvement = abs(tail[0] - tail[-1])
            converged = tail_improvement < 0.01 * max(abs(initial_loss), 1.0e-12)

        return TrajectoryRecord(
            trajectory_id=trajectory_id,
            perturbation_scale=perturbation_scale,
            initial_positions=positions_np.copy(),
            final_positions=final_positions_np,
            loss_history=loss_history,
            final_loss=final_loss,
            escaped=escaped,
            converged=converged,
            step_size=traj_step_size,
            langevin_scale=traj_langevin,
            perturbation_seed=base_seed,
        )

    def _apply_lattice_perturbation(self, prepared, scale: float, rng: np.random.Generator):
        """Return a copy of PreparedStructure with a jittered cell."""
        import copy as copy_module
        new_prepared = copy_module.copy(prepared)
        new_prepared.cell = self._perturb_lattice(prepared.cell, scale, rng)
        return new_prepared

    # ------------------------------------------------------------------
    # Escape detection
    # ------------------------------------------------------------------

    @staticmethod
    def _detect_escape(loss_history: List[float]) -> bool:
        """True if the trajectory escaped at least one local minimum.

        Detection criterion: the loss increased by more than 5 % of the
        initial loss at some step (indicating the trajectory traversed a
        barrier) before reaching a final loss lower than the pre-increase
        value.
        """
        if len(loss_history) < 4:
            return False
        initial = loss_history[0]
        threshold = 0.05 * abs(initial) + 1.0e-12
        min_before = loss_history[0]
        for idx in range(1, len(loss_history)):
            current = loss_history[idx]
            if current > min_before + threshold:
                # Loss increased — check if it recovers
                recovery_min = min(loss_history[idx:])
                if recovery_min < min_before:
                    return True
            min_before = min(min_before, current)
        return False

    # ------------------------------------------------------------------
    # Trajectory pruning
    # ------------------------------------------------------------------

    def _prune_trajectories(self, records: List[TrajectoryRecord]) -> None:
        """Mark outlier trajectories as pruned in-place.

        Trajectories with final loss above the configured quantile of the
        ensemble loss distribution are pruned, provided that at least one
        active trajectory would remain.
        """
        active = [r for r in records if not r.pruned]
        if len(active) <= 1:
            return
        losses = np.array([r.final_loss for r in active], dtype=np.float64)
        threshold = float(np.quantile(losses, self._config.consensus_prune_quantile))
        kept = 0
        for record in active:
            if record.final_loss <= threshold:
                kept += 1
        if kept == 0:
            return
        for record in active:
            if record.final_loss > threshold:
                record.pruned = True

    # ------------------------------------------------------------------
    # Statistical aggregation
    # ------------------------------------------------------------------

    def _compute_ensemble_statistics(
        self,
        records: List[TrajectoryRecord],
        perturbation_scale: float,
    ) -> EnsembleStatistics:
        """Aggregate trajectory endpoints into ensemble statistics."""
        active = [r for r in records if not r.pruned]
        if not active:
            active = records  # fallback: use all records

        positions_list = [r.final_positions for r in active]
        losses = np.array([r.final_loss for r in active], dtype=np.float64)

        # Positional statistics
        positions_array = np.stack(positions_list, axis=0)  # (n, N_atoms, 3)
        mean_pos = float_array(positions_array.mean(axis=0))
        median_pos = float_array(np.median(positions_array, axis=0))
        pos_var = float_array(positions_array.var(axis=0))

        # Pairwise RMSD divergence matrix
        pairwise = _pairwise_rmsd(positions_list)
        n = len(active)
        if n > 1:
            off_diag = pairwise[np.triu_indices(n, k=1)]
            position_spread = float(off_diag.mean()) if off_diag.size > 0 else 0.0
        else:
            position_spread = 0.0

        # Loss statistics
        mean_loss = float(losses.mean())
        median_loss = float(np.median(losses))
        loss_var = float(losses.var()) if len(losses) > 1 else 0.0

        # Entropy of normalised loss distribution
        loss_range = float(losses.max() - losses.min())
        if loss_range > 1.0e-12:
            probs = (losses.max() - losses) / loss_range
            probs = probs / probs.sum().clip(min=1.0e-12)
            entropy = float(-np.sum(probs * np.log(probs.clip(min=1.0e-300))))
        else:
            entropy = 0.0

        # Clustering for consensus
        cluster_labels = _cluster_by_rmsd(
            positions_list,
            threshold=float(self._config.cluster_distance_threshold),
        )
        cluster_counts = np.bincount(cluster_labels)
        num_clusters = int(cluster_counts.size)
        consensus_strength = float(cluster_counts.max()) / n if n > 0 else 1.0

        # Uncertainty: high spread + low consensus = high uncertainty
        uncertainty = position_spread / max(1.0 + consensus_strength, 1.0e-12)

        # Collapse detection
        collapse_detected = position_spread < float(self._config.diversity_collapse_threshold)

        # Escape count
        escape_count = sum(1 for r in active if r.escaped)

        # Ensemble-aware score per trajectory
        escape_rate = float(escape_count) / max(n, 1)
        best_idx = 0
        best_score = float("inf")
        for traj_idx, (record, loss_val) in enumerate(zip(active, losses)):
            score = (
                float(loss_val)
                + float(self._config.variance_weight) * loss_var
                - float(self._config.consensus_weight) * consensus_strength
                - float(self._config.escape_weight) * escape_rate
            )
            if score < best_score:
                best_score = score
                best_idx = traj_idx

        best_record = active[best_idx]

        new_scale = self._update_perturbation_scale(perturbation_scale, None, position_spread, collapse_detected)

        return EnsembleStatistics(
            mean_positions=mean_pos,
            median_positions=median_pos,
            position_variance=pos_var,
            position_spread=position_spread,
            mean_loss=mean_loss,
            median_loss=median_loss,
            loss_variance=loss_var,
            entropy=entropy,
            pairwise_divergence=pairwise,
            cluster_labels=cluster_labels,
            num_clusters=num_clusters,
            consensus_strength=consensus_strength,
            uncertainty=uncertainty,
            best_trajectory_id=best_record.trajectory_id,
            best_positions=best_record.final_positions.copy(),
            best_loss=float(losses[best_idx]),
            best_ensemble_score=best_score,
            escape_count=escape_count,
            collapse_detected=collapse_detected,
            active_trajectory_count=len(active),
            current_perturbation_scale=new_scale,
            trajectory_records=list(records),
        )

    # ------------------------------------------------------------------
    # Adaptive scaling
    # ------------------------------------------------------------------

    def _update_perturbation_scale(
        self,
        current_scale: float,
        stats: Optional[EnsembleStatistics] = None,
        position_spread: Optional[float] = None,
        collapse_detected: Optional[bool] = None,
    ) -> float:
        """Return updated perturbation scale based on ensemble diversity.

        Called both internally (during ``_compute_ensemble_statistics``)
        and externally (from ``run_ensemble_inverse_design``).
        """
        cfg = self._config
        if not cfg.adaptive_scaling:
            return current_scale

        # Accept spread/collapse either from stats object or directly
        if stats is not None:
            position_spread = stats.position_spread
            collapse_detected = stats.collapse_detected

        if position_spread is None:
            return current_scale

        new_scale = current_scale
        if collapse_detected or position_spread < float(cfg.diversity_collapse_threshold):
            # Trajectories have collapsed — increase exploration scale
            new_scale = current_scale * float(cfg.diversity_scale_up_factor)
        elif position_spread > 10.0 * float(cfg.diversity_collapse_threshold):
            # Trajectories are well-spread — reduce noise gently
            new_scale = current_scale * float(cfg.diversity_scale_down_factor)

        new_scale = float(
            np.clip(new_scale, float(cfg.min_perturbation_scale), float(cfg.max_perturbation_scale))
        )
        return new_scale


# ---------------------------------------------------------------------------
# Helper used inside compute_ensemble_statistics (avoids repeated dtype cast)
# ---------------------------------------------------------------------------

def float_array(arr: np.ndarray) -> np.ndarray:
    return np.asarray(arr, dtype=np.float32)


# ---------------------------------------------------------------------------
# Convenience factory
# ---------------------------------------------------------------------------

def make_ensemble_config_from_dict(raw: dict) -> InferenceEnsembleConfig:
    """Build an ``InferenceEnsembleConfig`` from a plain dict (e.g. from TOML).

    Unknown keys are silently ignored so that config files can include
    additional metadata without breaking the loader.

    Expected TOML section::

        [inverse_design.ensemble_exploration]
        enabled = true
        num_trajectories = 32
        perturbation_scale = 0.01
        adaptive_scaling = true
        aggregation = "mean_variance"
        consensus_metric = "cluster"
        escape_detection = true
        trajectory_pruning = true
        parallel = true
    """
    known_fields = {f.name for f in InferenceEnsembleConfig.__dataclass_fields__.values()}
    filtered = {k: v for k, v in raw.items() if k in known_fields}
    return InferenceEnsembleConfig(**filtered)
