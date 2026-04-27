from __future__ import annotations

import numpy as np


def wrap_atoms_to_unit_cell(atoms):
    wrapped = atoms.copy()
    wrapped.wrap(eps=1.0e-12)
    return wrapped


def _has_periodic_cell(atoms) -> bool:
    cell = atoms.cell.array
    if not np.all(np.asarray(atoms.pbc, dtype=bool)):
        return False
    return abs(float(np.linalg.det(cell))) > 1.0e-12


def _minimum_image_delta(cell: np.ndarray, delta: np.ndarray) -> np.ndarray:
    inverse_cell = np.linalg.inv(cell)
    delta_fractional = delta @ inverse_cell
    delta_fractional -= np.round(delta_fractional)
    return delta_fractional @ cell


def minimum_image_displacements(reference, candidate) -> np.ndarray:
    delta = candidate.get_positions() - reference.get_positions()
    if not _has_periodic_cell(reference):
        return delta
    return _minimum_image_delta(reference.cell.array, delta)


def _kabsch_rotation(reference_positions: np.ndarray, candidate_positions: np.ndarray) -> np.ndarray:
    covariance = candidate_positions.T @ reference_positions
    left, _, right_t = np.linalg.svd(covariance)
    rotation = left @ right_t
    if np.linalg.det(rotation) < 0.0:
        left[:, -1] *= -1.0
        rotation = left @ right_t
    return rotation


def _aligned_displacements(
    reference_positions: np.ndarray,
    candidate_positions: np.ndarray,
    cell: np.ndarray,
    use_periodic_cell: bool,
    allow_rotation: bool,
) -> np.ndarray:
    reference_centroid = reference_positions.mean(axis=0, keepdims=True)
    candidate_centroid = candidate_positions.mean(axis=0, keepdims=True)
    reference_centered = reference_positions - reference_centroid
    candidate_centered = candidate_positions - candidate_centroid

    aligned_candidate = candidate_centered
    if allow_rotation and len(reference_positions) > 1:
        rotation = _kabsch_rotation(reference_centered, candidate_centered)
        aligned_candidate = candidate_centered @ rotation
    aligned_candidate = aligned_candidate + reference_centroid

    delta = aligned_candidate - reference_positions
    if not use_periodic_cell:
        return delta
    return _minimum_image_delta(cell, delta)


def symmetry_aware_displacements(reference, candidate, allow_rotation: bool = True) -> np.ndarray:
    if len(reference) != len(candidate):
        raise ValueError("reference and candidate must contain the same number of atoms")

    reference_positions = reference.get_positions()
    candidate_positions = candidate.get_positions()
    use_periodic_cell = _has_periodic_cell(reference)
    cell = reference.cell.array

    candidates = [
        _aligned_displacements(
            reference_positions,
            candidate_positions,
            cell,
            use_periodic_cell=use_periodic_cell,
            allow_rotation=allow_rotation,
        )
    ]

    if use_periodic_cell:
        candidates.append(
            _aligned_displacements(
                reference_positions,
                reference_positions + minimum_image_displacements(reference, candidate),
                cell,
                use_periodic_cell=True,
                allow_rotation=allow_rotation,
            )
        )

    return min(candidates, key=lambda delta: float(np.mean(np.sum(delta ** 2, axis=1))))


def symmetry_aware_rmsd(reference, candidate, allow_rotation: bool = True) -> float:
    delta_cartesian = symmetry_aware_displacements(
        reference,
        candidate,
        allow_rotation=allow_rotation,
    )
    return float(np.sqrt(np.mean(np.sum(delta_cartesian ** 2, axis=1))))
