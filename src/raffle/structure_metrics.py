from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from ase import Atoms
from ase.neighborlist import neighbor_list
from scipy.optimize import linear_sum_assignment
import spglib


_ENVIRONMENT_MISMATCH_PENALTY = 5.0
_PAIR_MISMATCH_PENALTY = 10.0
_INVALID_ENVIRONMENT_COST = 1.0e6


@dataclass(frozen=True)
class _AlignmentResult:
    displacements: np.ndarray
    translation_vector: np.ndarray
    rmsd: float


def wrap_atoms_to_unit_cell(atoms):
    wrapped = atoms.copy()
    wrapped.wrap(eps=1.0e-12)
    return wrapped


def _has_periodic_cell(atoms) -> bool:
    cell = atoms.cell.array
    if not np.all(np.asarray(atoms.pbc, dtype=bool)):
        return False
    return abs(float(np.linalg.det(cell))) > 1.0e-12


def standardise_atoms(atoms, symprec: float = 1.0e-5) -> Atoms:
    if not _has_periodic_cell(atoms):
        return atoms.copy()

    wrapped = wrap_atoms_to_unit_cell(atoms)
    primitive = spglib.standardize_cell(
        (
            wrapped.cell.array,
            wrapped.get_scaled_positions(wrap=False),
            wrapped.numbers,
        ),
        to_primitive=True,
        no_idealize=False,
        symprec=symprec,
    )
    if primitive is None:
        return wrapped

    lattice, positions, numbers = primitive
    return wrap_atoms_to_unit_cell(
        Atoms(
            numbers=numbers,
            scaled_positions=positions,
            cell=lattice,
            pbc=True,
        )
    )


def minimum_image_displacement(delta_fractional: np.ndarray) -> np.ndarray:
    return delta_fractional - np.round(delta_fractional)


def _minimum_image_delta(cell: np.ndarray, delta: np.ndarray) -> np.ndarray:
    inverse_cell = np.linalg.inv(cell)
    delta_fractional = delta @ inverse_cell
    delta_fractional = minimum_image_displacement(delta_fractional)
    return delta_fractional @ cell


def minimum_image_displacements(reference, candidate) -> np.ndarray:
    if len(reference) != len(candidate):
        raise ValueError("reference and candidate must contain the same number of atoms")
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


def _translation_vector(
    reference_positions: np.ndarray,
    candidate_positions: np.ndarray,
    cell: np.ndarray,
    use_periodic_cell: bool,
) -> np.ndarray:
    delta = reference_positions - candidate_positions
    if use_periodic_cell:
        delta = _minimum_image_delta(cell, delta)
    return delta.mean(axis=0, keepdims=True)


def _direct_alignment_result(
    reference_positions: np.ndarray,
    candidate_positions: np.ndarray,
    cell: np.ndarray,
    use_periodic_cell: bool,
    allow_rotation: bool,
) -> _AlignmentResult:
    reference_centroid = reference_positions.mean(axis=0, keepdims=True)
    candidate_centroid = candidate_positions.mean(axis=0, keepdims=True)
    reference_centered = reference_positions - reference_centroid
    candidate_centered = candidate_positions - candidate_centroid

    aligned_candidate = candidate_centered
    if allow_rotation and len(reference_positions) > 1:
        rotation = _kabsch_rotation(reference_centered, candidate_centered)
        aligned_candidate = candidate_centered @ rotation
    aligned_candidate = aligned_candidate + candidate_centroid

    translation_vector = _translation_vector(
        reference_positions,
        aligned_candidate,
        cell,
        use_periodic_cell,
    )
    aligned_candidate = aligned_candidate + translation_vector

    delta = aligned_candidate - reference_positions
    if use_periodic_cell:
        delta = _minimum_image_delta(cell, delta)

    rmsd = float(np.sqrt(np.mean(np.sum(delta ** 2, axis=1)))) if len(delta) else 0.0
    return _AlignmentResult(
        displacements=delta,
        translation_vector=translation_vector.reshape(3).copy(),
        rmsd=rmsd,
    )


def _best_direct_alignment_result(reference, candidate, allow_rotation: bool) -> _AlignmentResult:
    if len(reference) != len(candidate):
        raise ValueError("reference and candidate must contain the same number of atoms")

    reference_positions = reference.get_positions()
    candidate_positions = candidate.get_positions()
    use_periodic_cell = _has_periodic_cell(reference) and _has_periodic_cell(candidate)
    cell = reference.cell.array

    candidates = [
        _direct_alignment_result(
            reference_positions,
            candidate_positions,
            cell,
            use_periodic_cell=use_periodic_cell,
            allow_rotation=allow_rotation,
        )
    ]

    if use_periodic_cell:
        candidates.append(
            _direct_alignment_result(
                reference_positions,
                reference_positions + minimum_image_displacements(reference, candidate),
                cell,
                use_periodic_cell=True,
                allow_rotation=allow_rotation,
            )
        )

    return min(candidates, key=lambda result: result.rmsd)


def _build_local_environments(atoms, rcut: float) -> list[list[tuple[int, float]]]:
    indices_i, indices_j, distances = neighbor_list("ijd", atoms, cutoff=float(rcut))
    environments: list[list[tuple[int, float]]] = [[] for _ in range(len(atoms))]
    for index_i, index_j, distance in zip(indices_i, indices_j, distances):
        environments[index_i].append((int(atoms.numbers[index_j]), round(float(distance), 4)))
    for environment in environments:
        environment.sort()
    return environments


def _environment_distance(
    environment_a: list[tuple[int, float]],
    environment_b: list[tuple[int, float]],
) -> float:
    size = max(len(environment_a), len(environment_b))
    if size == 0:
        return 0.0

    cost = np.full((size, size), _ENVIRONMENT_MISMATCH_PENALTY, dtype=float)
    for row_index, (atomic_number_a, distance_a) in enumerate(environment_a):
        for column_index, (atomic_number_b, distance_b) in enumerate(environment_b):
            if atomic_number_a != atomic_number_b:
                continue
            cost[row_index, column_index] = abs(distance_a - distance_b)

    row_ind, col_ind = linear_sum_assignment(cost)
    values = cost[row_ind, col_ind]
    return float(np.sqrt(np.mean(values ** 2)))


def _solve_square_assignment(cost_matrix: np.ndarray, pad_cost: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    row_count, column_count = cost_matrix.shape
    size = max(row_count, column_count)
    square_cost = np.full((size, size), float(pad_cost), dtype=float)
    square_cost[:row_count, :column_count] = cost_matrix
    row_ind, col_ind = linear_sum_assignment(square_cost)
    return square_cost, row_ind, col_ind


def _real_assignment_pairs(
    row_ind: np.ndarray,
    col_ind: np.ndarray,
    row_count: int,
    column_count: int,
) -> list[tuple[int, int]]:
    return [
        (int(row), int(column))
        for row, column in zip(row_ind, col_ind)
        if row < row_count and column < column_count
    ]


def _append_unique_translation(
    translations: list[np.ndarray],
    candidate: np.ndarray,
) -> None:
    for existing in translations:
        if np.allclose(existing, candidate, atol=1.0e-10, rtol=0.0):
            return
    translations.append(np.asarray(candidate, dtype=float).reshape(3).copy())


def _periodic_distance_matrix(
    fractional_positions_a: np.ndarray,
    fractional_positions_b: np.ndarray,
    cell: np.ndarray,
) -> np.ndarray:
    delta_fractional = fractional_positions_a[:, None, :] - fractional_positions_b[None, :, :]
    delta_fractional = minimum_image_displacement(delta_fractional)
    return np.linalg.norm(delta_fractional @ cell, axis=-1)


def _cartesian_distance_matrix(
    cartesian_positions_a: np.ndarray,
    cartesian_positions_b: np.ndarray,
) -> np.ndarray:
    delta = cartesian_positions_a[:, None, :] - cartesian_positions_b[None, :, :]
    return np.linalg.norm(delta, axis=-1)


def _prepare_atoms(atoms, use_primitive_representation: bool, symprec: float):
    if use_primitive_representation:
        return standardise_atoms(atoms, symprec=symprec)
    if _has_periodic_cell(atoms):
        return wrap_atoms_to_unit_cell(atoms)
    return atoms.copy()


def _rotated_candidate(reference, candidate):
    rotated = candidate.copy()
    reference_positions = reference.get_positions()
    candidate_positions = candidate.get_positions()
    reference_centroid = reference_positions.mean(axis=0, keepdims=True)
    candidate_centroid = candidate_positions.mean(axis=0, keepdims=True)
    rotation = _kabsch_rotation(
        reference_positions - reference_centroid,
        candidate_positions - candidate_centroid,
    )
    rotated.set_positions((candidate_positions - candidate_centroid) @ rotation + candidate_centroid)
    return rotated


def _translation_candidates(reference, candidate, environment_cost: np.ndarray) -> list[np.ndarray]:
    square_cost, row_ind, col_ind = _solve_square_assignment(
        environment_cost,
        _ENVIRONMENT_MISMATCH_PENALTY,
    )
    translations: list[np.ndarray] = []
    _append_unique_translation(translations, np.zeros(3, dtype=float))

    real_pairs = _real_assignment_pairs(row_ind, col_ind, len(reference), len(candidate))
    if not real_pairs:
        return translations

    if _has_periodic_cell(reference) and _has_periodic_cell(candidate):
        reference_positions = reference.get_scaled_positions(wrap=False)
        candidate_positions = candidate.get_scaled_positions(wrap=False)
        deltas = []
        for row_index, column_index in real_pairs:
            if square_cost[row_index, column_index] >= _INVALID_ENVIRONMENT_COST:
                continue
            delta = minimum_image_displacement(
                reference_positions[row_index] - candidate_positions[column_index]
            )
            deltas.append(delta)
            _append_unique_translation(translations, delta)
    else:
        reference_positions = reference.get_positions()
        candidate_positions = candidate.get_positions()
        deltas = []
        for row_index, column_index in real_pairs:
            if square_cost[row_index, column_index] >= _INVALID_ENVIRONMENT_COST:
                continue
            delta = reference_positions[row_index] - candidate_positions[column_index]
            deltas.append(delta)
            _append_unique_translation(translations, delta)

    if deltas:
        _append_unique_translation(translations, np.mean(np.asarray(deltas, dtype=float), axis=0))

    return translations


def _evaluate_alignment(reference, candidate, translation: np.ndarray) -> _AlignmentResult:
    periodic = _has_periodic_cell(reference) and _has_periodic_cell(candidate)
    reference_numbers = np.asarray(reference.numbers, dtype=int)
    candidate_numbers = np.asarray(candidate.numbers, dtype=int)

    if periodic:
        reference_fractional = reference.get_scaled_positions(wrap=False)
        shifted_fractional = (candidate.get_scaled_positions(wrap=False) + translation) % 1.0
        reference_positions = reference_fractional @ reference.cell.array
        shifted_candidate_positions = shifted_fractional @ reference.cell.array
        distances = _periodic_distance_matrix(
            reference_fractional,
            shifted_fractional,
            reference.cell.array,
        )
        translation_vector = np.asarray(translation @ reference.cell.array, dtype=float)
    else:
        reference_positions = reference.get_positions()
        shifted_candidate_positions = candidate.get_positions() + translation
        distances = _cartesian_distance_matrix(reference_positions, shifted_candidate_positions)
        translation_vector = np.asarray(translation, dtype=float)

    pair_cost = np.full((len(reference), len(candidate)), _PAIR_MISMATCH_PENALTY, dtype=float)
    matching_species = reference_numbers[:, None] == candidate_numbers[None, :]
    pair_cost[matching_species] = distances[matching_species]

    square_cost, row_ind, col_ind = _solve_square_assignment(pair_cost, _PAIR_MISMATCH_PENALTY)
    values = square_cost[row_ind, col_ind]
    displacements = np.full((len(reference), 3), np.nan, dtype=float)
    for row_index, column_index in _real_assignment_pairs(row_ind, col_ind, len(reference), len(candidate)):
        delta = shifted_candidate_positions[column_index] - reference_positions[row_index]
        if periodic:
            delta = _minimum_image_delta(reference.cell.array, delta)
        displacements[row_index] = delta

    rmsd = float(np.sqrt(np.mean(values ** 2))) if values.size else 0.0
    return _AlignmentResult(
        displacements=displacements,
        translation_vector=translation_vector.reshape(3).copy(),
        rmsd=rmsd,
    )


def _best_alignment_result(
    reference,
    candidate,
    *,
    allow_rotation: bool,
    rcut: float,
    symprec: float,
    use_primitive_representation: bool,
) -> _AlignmentResult:
    use_periodic_matching = _has_periodic_cell(reference) and _has_periodic_cell(candidate)
    prepared_reference = _prepare_atoms(
        reference,
        use_primitive_representation=use_primitive_representation and use_periodic_matching,
        symprec=symprec,
    )

    candidate_versions = [candidate.copy()]
    if allow_rotation and len(reference) == len(candidate) and len(reference) > 1:
        candidate_versions.append(_rotated_candidate(reference, candidate))

    best_result: _AlignmentResult | None = None
    for candidate_version in candidate_versions:
        prepared_candidate = _prepare_atoms(
            candidate_version,
            use_primitive_representation=use_primitive_representation and use_periodic_matching,
            symprec=symprec,
        )
        environments_reference = _build_local_environments(prepared_reference, rcut=rcut)
        environments_candidate = _build_local_environments(prepared_candidate, rcut=rcut)

        environment_cost = np.full(
            (len(prepared_reference), len(prepared_candidate)),
            _INVALID_ENVIRONMENT_COST,
            dtype=float,
        )
        for row_index in range(len(prepared_reference)):
            for column_index in range(len(prepared_candidate)):
                if prepared_reference.numbers[row_index] != prepared_candidate.numbers[column_index]:
                    continue
                environment_cost[row_index, column_index] = _environment_distance(
                    environments_reference[row_index],
                    environments_candidate[column_index],
                )

        for translation in _translation_candidates(prepared_reference, prepared_candidate, environment_cost):
            result = _evaluate_alignment(prepared_reference, prepared_candidate, translation)
            if best_result is None or result.rmsd < best_result.rmsd:
                best_result = result

    if best_result is None:
        raise RuntimeError("Unable to align structures")
    return best_result


def structure_similarity_rmsd(
    reference,
    candidate,
    rcut: float = 6.0,
    return_translation: bool = False,
    allow_rotation: bool = True,
    symprec: float = 1.0e-5,
) -> float | tuple[float, np.ndarray]:
    result = _best_alignment_result(
        reference,
        candidate,
        allow_rotation=allow_rotation,
        rcut=rcut,
        symprec=symprec,
        use_primitive_representation=True,
    )
    if return_translation:
        return result.rmsd, result.translation_vector.copy()
    return result.rmsd


def symmetry_aware_displacements(reference, candidate, allow_rotation: bool = True) -> np.ndarray:
    result = _best_direct_alignment_result(
        reference,
        candidate,
        allow_rotation=allow_rotation,
    )
    return result.displacements.copy()
