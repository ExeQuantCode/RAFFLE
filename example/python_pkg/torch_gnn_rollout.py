"""Replay-buffer helpers for the PyTorch inverse-design workflow."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Mapping, Sequence

import numpy as np


@dataclass(frozen=True)
class ReplaySample:
    atoms: object
    priority: float
    categories: tuple[str, ...]
    stage_index: int
    step: int
    restart_index: int
    true_target_fingerprint_mse: float
    surrogate_target_fingerprint_mse: float
    fingerprint_drift_mse: float
    position_difference: float | None = None
    repulsion_loss: float = 0.0
    cell_violation_loss: float = 0.0
    is_initial_state: bool = False


class PrioritizedReplayBuffer:
    def __init__(self, capacity: int, seed: int = 0):
        self.capacity = max(int(capacity), 0)
        self._rng = np.random.default_rng(int(seed))
        self._samples: list[ReplaySample] = []
        self._next_order = 0
        self._order: list[int] = []

    def __len__(self) -> int:
        return len(self._samples)

    @property
    def samples(self) -> tuple[ReplaySample, ...]:
        return tuple(self._samples)

    def add(self, sample: ReplaySample) -> None:
        if self.capacity <= 0:
            return
        copied_sample = ReplaySample(
            atoms=sample.atoms.copy(),
            priority=max(float(sample.priority), 1.0e-8),
            categories=tuple(sorted(set(sample.categories))),
            stage_index=int(sample.stage_index),
            step=int(sample.step),
            restart_index=int(sample.restart_index),
            true_target_fingerprint_mse=float(sample.true_target_fingerprint_mse),
            surrogate_target_fingerprint_mse=float(sample.surrogate_target_fingerprint_mse),
            fingerprint_drift_mse=float(sample.fingerprint_drift_mse),
            position_difference=(
                None if sample.position_difference is None else float(sample.position_difference)
            ),
            repulsion_loss=float(sample.repulsion_loss),
            cell_violation_loss=float(sample.cell_violation_loss),
            is_initial_state=bool(sample.is_initial_state),
        )
        self._samples.append(copied_sample)
        self._order.append(self._next_order)
        self._next_order += 1
        self._trim_to_capacity()

    def extend(self, samples: Iterable[ReplaySample]) -> None:
        for sample in samples:
            self.add(sample)

    def sample(
        self,
        sample_size: int,
        category_weights: Mapping[str, float] | None = None,
    ) -> list[ReplaySample]:
        if not self._samples or int(sample_size) <= 0:
            return []

        resolved_weights = dict(category_weights or {})
        population = len(self._samples)
        size = min(int(sample_size), population)
        sample_weights = np.asarray(
            [self._sample_weight(sample, resolved_weights) for sample in self._samples],
            dtype=np.float64,
        )
        if not np.any(sample_weights > 0.0):
            sample_weights = np.full(population, 1.0 / population, dtype=np.float64)
        else:
            sample_weights = sample_weights / np.sum(sample_weights)

        chosen = self._rng.choice(population, size=size, replace=False, p=sample_weights)
        return [self._samples[int(index)] for index in chosen]

    def stats(self) -> dict[str, float | int | dict[str, int]]:
        category_counts: dict[str, int] = {}
        for sample in self._samples:
            for category in sample.categories:
                category_counts[category] = category_counts.get(category, 0) + 1

        if not self._samples:
            return {
                "capacity": int(self.capacity),
                "size": 0,
                "category_counts": category_counts,
                "mean_priority": 0.0,
                "mean_fingerprint_drift_mse": 0.0,
                "mean_true_target_fingerprint_mse": 0.0,
            }

        priorities = np.asarray([sample.priority for sample in self._samples], dtype=np.float64)
        drifts = np.asarray(
            [sample.fingerprint_drift_mse for sample in self._samples],
            dtype=np.float64,
        )
        true_errors = np.asarray(
            [sample.true_target_fingerprint_mse for sample in self._samples],
            dtype=np.float64,
        )
        return {
            "capacity": int(self.capacity),
            "size": int(len(self._samples)),
            "category_counts": category_counts,
            "mean_priority": float(np.mean(priorities)),
            "mean_fingerprint_drift_mse": float(np.mean(drifts)),
            "mean_true_target_fingerprint_mse": float(np.mean(true_errors)),
        }

    def _sample_weight(
        self,
        sample: ReplaySample,
        category_weights: Mapping[str, float],
    ) -> float:
        category_factor = sum(float(category_weights.get(category, 1.0)) for category in sample.categories)
        if category_factor <= 0.0:
            category_factor = 1.0
        return max(float(sample.priority), 1.0e-8) * category_factor

    def _trim_to_capacity(self) -> None:
        while len(self._samples) > self.capacity:
            discard_index = self._discard_index()
            del self._samples[discard_index]
            del self._order[discard_index]

    def _discard_index(self) -> int:
        discard_index = 0
        discard_key = (self._samples[0].priority, self._order[0])
        for index, (sample, order) in enumerate(zip(self._samples, self._order)):
            candidate_key = (sample.priority, order)
            if candidate_key < discard_key:
                discard_key = candidate_key
                discard_index = index
        return discard_index


def classify_rollout_step(
    *,
    true_target_fingerprint_mse: float,
    previous_true_target_fingerprint_mse: float | None,
    fingerprint_drift_mse: float,
    repulsion_loss: float,
    cell_violation_loss: float,
    drift_threshold: float,
    high_error_threshold: float,
    instability_threshold: float,
) -> tuple[str, ...]:
    categories: set[str] = set()
    true_error = float(true_target_fingerprint_mse)
    previous_true_error = (
        None if previous_true_target_fingerprint_mse is None else float(previous_true_target_fingerprint_mse)
    )
    drift = float(fingerprint_drift_mse)

    if previous_true_error is None or true_error <= previous_true_error:
        categories.add("successful")
    else:
        categories.add("failed")

    if drift >= float(drift_threshold):
        categories.add("difficult")
    if true_error >= float(high_error_threshold):
        categories.add("high_error")
    if (
        float(repulsion_loss) >= float(instability_threshold)
        or float(cell_violation_loss) >= float(instability_threshold)
        or (previous_true_error is not None and true_error > previous_true_error)
    ):
        categories.add("unstable")

    return tuple(sorted(categories))
