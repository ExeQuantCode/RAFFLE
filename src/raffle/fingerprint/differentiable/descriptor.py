"""PyTorch implementation of RAFFLE's 2-, 3-, and 4-body fingerprints."""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime
import math
from typing import Any, Mapping, Sequence

import numpy as np
import torch
from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes
from ase.data import atomic_numbers, covalent_radii
from ase.neighborlist import primitive_neighbor_list


@dataclass
class FingerprintResult:
    """Fingerprint distributions and their pre-normalisation contributions."""

    two_body: torch.Tensor
    three_body: torch.Tensor
    four_body: torch.Tensor
    unnormalized: dict[str, torch.Tensor]
    grids: dict[str, torch.Tensor]
    elements: tuple[str, ...]
    composition: dict[str, int] | None = None

    @property
    def df_2body(self) -> torch.Tensor:
        return self.two_body

    @property
    def df_3body(self) -> torch.Tensor:
        return self.three_body

    @property
    def df_4body(self) -> torch.Tensor:
        return self.four_body


@dataclass
class RelaxationResult:
    """Output from a coupled position/species relaxation."""

    structure: Atoms
    fingerprint: FingerprintResult
    loss: float
    converged: bool
    n_steps: int
    species_probabilities: torch.Tensor
    composition: dict[str, int]
    loss_history: list[dict[str, float]]
    trajectory: list[Atoms]
    species_probabilities_history: list[torch.Tensor] = field(default_factory=list)


class _NeighborTopology:
    """Cached periodic candidate pairs; only topology construction uses NumPy."""

    def __init__(self, cutoff: float, skin: float = 0.5):
        if cutoff <= 0.0 or skin < 0.0:
            raise ValueError("cutoff must be positive and skin must be non-negative")
        self.cutoff = float(cutoff)
        self.skin = float(skin)
        self.reference_positions: np.ndarray | None = None
        self.reference_cell: np.ndarray | None = None
        self.reference_pbc: np.ndarray | None = None
        self.i = np.empty(0, dtype=np.int64)
        self.j = np.empty(0, dtype=np.int64)
        self.shifts = np.empty((0, 3), dtype=np.int64)

    def pairs(
        self,
        positions: torch.Tensor,
        cell: torch.Tensor,
        pbc: Sequence[bool],
    ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        pos = positions.detach().cpu().numpy()
        cell_np = cell.detach().cpu().numpy()
        pbc_np = np.asarray(pbc, dtype=bool)
        rebuild = self.reference_positions is None
        if not rebuild:
            rebuild = (
                self.reference_positions.shape != pos.shape
                or not np.array_equal(self.reference_pbc, pbc_np)
                or not np.allclose(self.reference_cell, cell_np, rtol=0.0, atol=1e-12)
            )
        if not rebuild:
            displacement = np.linalg.norm(pos - self.reference_positions, axis=1).max(
                initial=0.0
            )
            rebuild = displacement > self.skin / 2.0 if self.skin else displacement > 0.0
        if rebuild:
            self.i, self.j, self.shifts = primitive_neighbor_list(
                "ijS",
                pbc_np,
                cell_np,
                pos,
                self.cutoff + self.skin,
                self_interaction=True,
            )
            order = np.lexsort(
                (
                    self.shifts[:, 2],
                    self.shifts[:, 1],
                    self.shifts[:, 0],
                    self.j,
                    np.any(self.shifts != 0, axis=1),
                    self.i,
                )
            )
            self.i = self.i[order]
            self.j = self.j[order]
            self.shifts = self.shifts[order]
            self.reference_positions = pos.copy()
            self.reference_cell = cell_np.copy()
            self.reference_pbc = pbc_np.copy()
        device = positions.device
        return (
            torch.as_tensor(self.i, dtype=torch.long, device=device),
            torch.as_tensor(self.j, dtype=torch.long, device=device),
            torch.as_tensor(self.shifts, dtype=torch.long, device=device),
        )


class DifferentiableFingerprint:
    """Direct differentiable implementation of the RAFFLE structure descriptor.

    Parameters may be a ``Raffle__Distribs_Container`` or a mapping containing
    ``width``, ``sigma``, ``cutoff_min``, ``cutoff_max``,
    ``radius_distance_tol``, and optionally ``bond_radii``.  The default
    parameters are the distribution-container defaults used by RAFFLE.
    """

    def __init__(
        self,
        elements: Sequence[str],
        parameters: Any = None,
        *,
        bond_radii: Mapping[tuple[str, str], float] | None = None,
        skin: float = 0.5,
    ):
        normalized_elements = tuple(sorted({str(element).strip() for element in elements}))
        if not normalized_elements or any(
            element not in atomic_numbers or atomic_numbers[element] == 0
            for element in normalized_elements
        ):
            raise ValueError("elements must contain valid chemical element symbols")
        self.elements = normalized_elements
        self.width = self._parameter(parameters, "width", [0.025, math.pi / 64, math.pi / 64])
        self.sigma = self._parameter(parameters, "sigma", [0.1, 0.1, 0.1])
        self.cutoff_min = self._parameter(parameters, "cutoff_min", [0.5, 0.0, 0.0])
        self.cutoff_max = self._parameter(parameters, "cutoff_max", [6.0, math.pi, math.pi])
        self.radius_distance_tol = self._parameter(
            parameters, "radius_distance_tol", [1.5, 2.5, 3.0, 6.0]
        )
        self._validate_parameters()

        configured_radii: dict[tuple[str, str], float] = {}
        if parameters is not None:
            getter = getattr(parameters, "get_bond_radii", None)
            if callable(getter):
                configured_radii.update(getter())
            elif isinstance(parameters, Mapping):
                configured_radii.update(parameters.get("bond_radii", {}))
        if bond_radii is not None:
            configured_radii.update(bond_radii)
        self.bond_radii = {
            self._pair_key(pair[0], pair[1]): float(radius)
            for pair, radius in configured_radii.items()
        }
        self._radii = torch.tensor(
            [
                [
                    self.bond_radii.get(
                        self._pair_key(a, b),
                        0.5
                        * (
                            float(covalent_radii[atomic_numbers[a]])
                            + float(covalent_radii[atomic_numbers[b]])
                        ),
                    )
                    for b in self.elements
                ]
                for a in self.elements
            ],
            dtype=torch.float64,
        )
        self._topology = _NeighborTopology(max(self.cutoff_max[0], 1e-4), skin)

    @staticmethod
    def _pair_key(a: str, b: str) -> tuple[str, str]:
        return tuple(sorted((str(a).strip(), str(b).strip())))

    @staticmethod
    def _parameter(parameters: Any, name: str, default: Sequence[float]) -> list[float]:
        if parameters is None:
            value = default
        elif isinstance(parameters, Mapping):
            value = parameters.get(name, default)
        else:
            value = getattr(parameters, name, default)
        result = [float(item) for item in value]
        if len(result) != len(default):
            raise ValueError(f"{name} must contain exactly {len(default)} values")
        return result

    def _validate_parameters(self) -> None:
        if any(width <= 0 for width in self.width):
            raise ValueError("all descriptor widths must be positive")
        if any(sigma <= 0 for sigma in self.sigma):
            raise ValueError("all Gaussian sigma values must be positive")
        if any(lo >= hi for lo, hi in zip(self.cutoff_min, self.cutoff_max)):
            raise ValueError("each cutoff_min value must be smaller than cutoff_max")
        if any(value <= 0 for value in self.radius_distance_tol):
            raise ValueError("radius_distance_tol values must be positive")

    def _grids(self, dtype: torch.dtype, device: torch.device) -> dict[str, torch.Tensor]:
        return {
            name: lo + torch.arange(
                int(math.floor((hi - lo) / width + 0.5)) + 1,
                dtype=dtype,
                device=device,
            ) * width
            for name, lo, hi, width in zip(
                ("two_body", "three_body", "four_body"),
                self.cutoff_min,
                self.cutoff_max,
                self.width,
            )
        }

    @staticmethod
    def _project(
        values: torch.Tensor,
        weights: torch.Tensor,
        denominator: torch.Tensor,
        grid: torch.Tensor,
        sigma: float,
        width: float,
        cutoff_min: float,
    ) -> torch.Tensor:
        eta = 1.0 / (2.0 * sigma * sigma)
        if values.numel() == 0:
            return grid.new_zeros(grid.shape)
        nearest = torch.floor((values - cutoff_min) / width + 0.5)
        bin_index = torch.arange(grid.numel(), device=grid.device, dtype=grid.dtype)
        support = (bin_index[None, :] - nearest[:, None]).abs() <= math.ceil(
            math.sqrt(16.0 / eta) / width
        )
        gaussian = torch.exp(-eta * (values[:, None] - grid[None, :]).square())
        projected = (gaussian * support.to(gaussian.dtype) * weights[:, None]).sum(dim=0)
        return projected * math.sqrt(eta / math.pi) / denominator.clamp_min(
            torch.finfo(grid.dtype).tiny
        )

    @staticmethod
    def _cosine_window(
        distance: torch.Tensor, lower: torch.Tensor, upper: torch.Tensor
    ) -> tuple[torch.Tensor, torch.Tensor]:
        valid = (distance[:, None] >= lower[None, :]) & (
            distance[:, None] <= upper[None, :]
        )
        span = (upper - lower).clamp_min(torch.finfo(distance.dtype).eps)
        phase = 2.0 * math.pi * (distance[:, None] - lower[None, :]) / span[None, :]
        window = 0.5 * (1.0 - torch.cos(phase))
        return (window * valid.to(window.dtype)), valid.to(window.dtype)

    @staticmethod
    def _normalise(values: torch.Tensor) -> torch.Tensor:
        totals = values.sum(dim=0, keepdim=True)
        active = values.abs().amax(dim=0, keepdim=True) > 1e-6
        return torch.where(active, values / totals.clamp_min(torch.finfo(values.dtype).tiny), values)

    def __call__(
        self,
        positions: torch.Tensor | np.ndarray,
        cell: torch.Tensor | np.ndarray,
        species_probabilities: torch.Tensor | np.ndarray,
        pbc: Sequence[bool] = (True, True, True),
    ) -> FingerprintResult:
        positions = torch.as_tensor(positions)
        if not positions.is_floating_point():
            positions = positions.to(torch.float64)
        device, dtype = positions.device, positions.dtype
        cell = torch.as_tensor(cell, dtype=dtype, device=device)
        probabilities = torch.as_tensor(species_probabilities, dtype=dtype, device=device)
        n_atoms, n_elements = probabilities.shape
        if positions.shape != (n_atoms, 3):
            raise ValueError("positions must have shape [n_atoms, 3]")
        if cell.shape != (3, 3):
            raise ValueError("cell must have shape [3, 3]")
        if n_elements != len(self.elements) or n_atoms == 0:
            raise ValueError("species_probabilities must have shape [n_atoms, len(elements)]")
        if len(pbc) != 3:
            raise ValueError("pbc must have three Boolean values")
        if not torch.isfinite(positions).all() or not torch.isfinite(probabilities).all():
            raise ValueError("positions and species probabilities must be finite")
        if (probabilities < 0).any():
            raise ValueError("species probabilities cannot be negative")
        if not torch.allclose(
            probabilities.sum(dim=1),
            torch.ones(n_atoms, dtype=dtype, device=device),
            atol=1e-6,
            rtol=1e-6,
        ):
            raise ValueError("each species-probability row must sum to one")

        radii = self._radii.to(dtype=dtype, device=device)
        grids = self._grids(dtype, device)
        pair_indices = [(a, b) for a in range(n_elements) for b in range(a, n_elements)]
        raw2 = positions.new_zeros((grids["two_body"].numel(), len(pair_indices)))
        raw3 = positions.new_zeros((grids["three_body"].numel(), n_elements))
        raw4 = positions.new_zeros((grids["four_body"].numel(), n_elements))

        center, neighbor, shifts = self._topology.pairs(positions, cell, pbc)
        vectors = positions[neighbor] - positions[center] + shifts.to(dtype) @ cell
        distances = torch.linalg.vector_norm(vectors, dim=1)
        image_mask = shifts.ne(0).any(dim=1)
        by_center = [torch.nonzero(center == atom, as_tuple=False).flatten() for atom in range(n_atoms)]
        rmin, rmax = self.cutoff_min[0], self.cutoff_max[0]
        radial = (distances >= rmin) & (distances <= rmax) & (distances > 1e-10)
        grids2 = grids["two_body"]
        eta2 = 1.0 / (2.0 * self.sigma[0] ** 2)

        for pair_channel, (alpha, beta) in enumerate(pair_indices):
            for central_element, neighbor_element in (
                ((alpha, beta),) if alpha == beta else ((alpha, beta), (beta, alpha))
            ):
                for atom in range(n_atoms):
                    ids = by_center[atom]
                    if ids.numel() == 0:
                        continue
                    valid_ids = ids[radial[ids]]
                    if valid_ids.numel() == 0:
                        continue
                    r = distances[valid_ids]
                    p_neighbor = probabilities[neighbor[valid_ids], neighbor_element]
                    denominator = p_neighbor.sum()
                    weights = probabilities[atom, central_element] * p_neighbor
                    gaussian = torch.exp(-eta2 * (r[:, None] - grids2[None, :]).square())
                    nearest = torch.floor((r - rmin) / self.width[0] + 0.5)
                    bins = torch.arange(grids2.numel(), dtype=dtype, device=device)
                    support = (bins[None, :] - nearest[:, None]).abs() <= math.ceil(
                        math.sqrt(16.0 / eta2) / self.width[0]
                    )
                    raw2[:, pair_channel] = raw2[:, pair_channel] + (
                        gaussian * support.to(dtype) * weights[:, None]
                    ).sum(dim=0) * math.sqrt(eta2 / math.pi) / denominator.clamp_min(
                        torch.finfo(dtype).tiny
                    )

        raw2 = raw2 / grids2[:, None].square()
        taper_max = self.cutoff_max[0] - 0.25
        taper_min = self.cutoff_min[0] + 0.25
        taper = torch.ones_like(grids2)
        upper = grids2 > taper_max
        lower = grids2 < taper_min
        taper = torch.where(
            upper,
            0.5 * (1.0 + torch.cos(math.pi * (grids2 - taper_max) / 0.25)),
            taper,
        )
        taper = torch.where(
            lower,
            0.5 * (1.0 + torch.cos(math.pi * (grids2 - taper_min) / 0.25)),
            taper,
        )
        raw2 = raw2 * taper[:, None]
        norm2 = self._normalise(raw2)

        grids3 = grids["three_body"]
        for atom in range(n_atoms):
            ids = by_center[atom]
            if ids.numel() < 2:
                continue
            for alpha in range(n_elements):
                d = distances[ids]
                pair_radii = radii[alpha]
                lower3 = torch.maximum(
                    d.new_full((n_elements,), self.cutoff_min[0]),
                    pair_radii * self.radius_distance_tol[0],
                )
                upper3 = torch.minimum(
                    d.new_full((n_elements,), self.cutoff_max[0]),
                    pair_radii * self.radius_distance_tol[1],
                )
                win3_by_species, valid3_by_species = self._cosine_window(
                    d, lower3, upper3
                )
                neighbor_probabilities = probabilities[neighbor[ids]]
                win3 = (win3_by_species * neighbor_probabilities).sum(dim=1)
                valid3 = (valid3_by_species * neighbor_probabilities).sum(dim=1)
                possible = valid3_by_species.any(dim=1)
                ids = ids[possible]
                d = d[possible]
                win3 = win3[possible]
                valid3 = valid3[possible]
                if ids.numel() < 2:
                    continue
                left, right = torch.triu_indices(ids.numel(), ids.numel(), offset=1, device=device)
                pair_count_weights = valid3[left] * valid3[right]
                denominator = pair_count_weights.sum()
                if denominator.detach().item() == 0.0:
                    continue
                bond_weights = (
                    win3[left] * win3[right]
                    / (d[left].square() * d[right].square()).clamp_min(1e-20)
                )
                cos_angle = (
                    vectors[ids[left]] * vectors[ids[right]]
                ).sum(dim=1) / (d[left] * d[right]).clamp_min(1e-12)
                angle_epsilon = max(torch.finfo(dtype).eps * 10.0, 1e-12)
                angles = torch.acos(
                    cos_angle.clamp(-1.0 + angle_epsilon, 1.0 - angle_epsilon)
                )
                angle_valid = torch.linalg.vector_norm(
                    torch.linalg.cross(vectors[ids[left]], vectors[ids[right]], dim=1),
                    dim=1,
                ) > 1e-12
                contribution = self._project(
                    angles,
                    bond_weights
                    * angle_valid.to(dtype)
                    * probabilities[atom, alpha],
                    denominator,
                    grids3,
                    self.sigma[1],
                    self.width[1],
                    self.cutoff_min[1],
                )
                raw3[:, alpha] = raw3[:, alpha] + contribution

        norm3 = self._normalise(raw3)

        grids4 = grids["four_body"]
        for atom in range(n_atoms):
            ids = by_center[atom]
            image_ids = ids[image_mask[ids]]
            if ids.numel() < 2 or image_ids.numel() == 0:
                continue
            for alpha in range(n_elements):
                d = distances[ids]
                pair_radii = radii[alpha]
                lower3 = torch.maximum(
                    d.new_full((n_elements,), self.cutoff_min[0]),
                    pair_radii * self.radius_distance_tol[0],
                )
                upper3 = torch.minimum(
                    d.new_full((n_elements,), self.cutoff_max[0]),
                    pair_radii * self.radius_distance_tol[1],
                )
                win3_by_species, valid3_by_species = self._cosine_window(
                    d, lower3, upper3
                )
                neighbor_probabilities = probabilities[neighbor[ids]]
                win3 = (win3_by_species * neighbor_probabilities).sum(dim=1)
                valid3 = (valid3_by_species * neighbor_probabilities).sum(dim=1)
                possible3 = valid3_by_species.any(dim=1)
                ids = ids[possible3]
                d = d[possible3]
                win3 = win3[possible3]
                valid3 = valid3[possible3]
                win3_by_species = win3_by_species[possible3]
                valid3_by_species = valid3_by_species[possible3]
                neighbor_probabilities = neighbor_probabilities[possible3]
                if ids.numel() < 2:
                    continue

                image_d = distances[image_ids]
                lower4 = torch.maximum(
                    image_d.new_full((n_elements,), self.cutoff_min[0]),
                    pair_radii * self.radius_distance_tol[2],
                )
                upper4 = torch.minimum(
                    image_d.new_full((n_elements,), self.cutoff_max[0]),
                    pair_radii * self.radius_distance_tol[3],
                )
                win4_by_species, valid4_by_species = self._cosine_window(
                    image_d, lower4, upper4
                )
                image_probabilities = probabilities[neighbor[image_ids]]
                win4 = (win4_by_species * image_probabilities).sum(dim=1)
                valid4 = (valid4_by_species * image_probabilities).sum(dim=1)
                possible4 = valid4_by_species.any(dim=1)
                image_ids = image_ids[possible4]
                image_d = image_d[possible4]
                win4 = win4[possible4]
                valid4 = valid4[possible4]
                win4_by_species = win4_by_species[possible4]
                image_probabilities = image_probabilities[possible4]
                if image_ids.numel() == 0:
                    continue

                left, right = torch.triu_indices(ids.numel(), ids.numel(), offset=1, device=device)
                pair_count_weights = valid3[left] * valid3[right]
                denominator = pair_count_weights.sum() * valid4.sum()
                if denominator.detach().item() == 0.0:
                    continue
                # Enumerate only local triplets, never a dense N^4 tensor.
                chunk = 4096
                n_combinations = left.numel() * image_ids.numel()
                angle_epsilon = max(torch.finfo(dtype).eps * 10.0, 1e-12)
                for beta in range(n_elements):
                    for gamma in range(n_elements):
                        for start in range(0, n_combinations, chunk):
                            flat_index = torch.arange(
                                start, min(start + chunk, n_combinations), device=device
                            )
                            pair_index = flat_index // image_ids.numel()
                            j, k = left[pair_index], right[pair_index]
                            l = flat_index % image_ids.numel()
                            vj = vectors[ids[j]]
                            vk = vectors[ids[k]]
                            first, second = (vk, vj) if beta > gamma else (vj, vk)
                            vl = vectors[image_ids[l]]
                            normal_a = torch.linalg.cross(first, second, dim=1)
                            normal_b = torch.linalg.cross(second, vl, dim=1)
                            norm_a = torch.linalg.vector_norm(normal_a, dim=1)
                            norm_b = torch.linalg.vector_norm(normal_b, dim=1)
                            valid_angle = (norm_a > 1e-12) & (norm_b > 1e-12)
                            cosine = (normal_a * normal_b).sum(dim=1) / (
                                norm_a * norm_b
                            ).clamp_min(1e-24)
                            angles = torch.acos(
                                cosine.clamp(
                                    -1.0 + angle_epsilon, 1.0 - angle_epsilon
                                )
                            )
                            weights = (
                                win3_by_species[j, beta]
                                * neighbor_probabilities[j, beta]
                                * win3_by_species[k, gamma]
                                * neighbor_probabilities[k, gamma]
                                * win4[l]
                                / (
                                    d[j].square()
                                    * d[k].square()
                                    * image_d[l].square()
                                ).clamp_min(1e-24)
                            )
                            raw4[:, alpha] = raw4[:, alpha] + self._project(
                                angles,
                                weights
                                * valid_angle.to(dtype)
                                * probabilities[atom, alpha],
                                denominator,
                                grids4,
                                self.sigma[2],
                                self.width[2],
                                self.cutoff_min[2],
                            )

        norm4 = self._normalise(raw4)
        composition = None
        return FingerprintResult(
            norm2,
            norm3,
            norm4,
            {"two_body": raw2, "three_body": raw3, "four_body": raw4},
            grids,
            self.elements,
            composition,
        )

    def generate(self, atoms: Atoms) -> FingerprintResult:
        """Generate a reusable target fingerprint from an ASE structure."""
        probabilities = self.probabilities_from_atoms(atoms)
        result = self(
            torch.as_tensor(atoms.positions, dtype=torch.float64),
            torch.as_tensor(np.asarray(atoms.cell), dtype=torch.float64),
            probabilities,
            atoms.pbc,
        )
        result.composition = {
            element: atoms.get_chemical_symbols().count(element) for element in self.elements
        }
        return result

    def probabilities_from_atoms(self, atoms: Atoms) -> torch.Tensor:
        missing = sorted(set(atoms.get_chemical_symbols()) - set(self.elements))
        if missing:
            raise ValueError(f"structure contains elements not configured in descriptor: {missing}")
        result = torch.zeros((len(atoms), len(self.elements)), dtype=torch.float64)
        for index, symbol in enumerate(atoms.get_chemical_symbols()):
            result[index, self.elements.index(symbol)] = 1.0
        return result


def _fingerprint_loss(
    current: FingerprintResult,
    target: FingerprintResult,
    weights: Sequence[float] = (1.0, 1.0, 1.0),
    radial_smoothing: float = 0.0,
) -> tuple[torch.Tensor, dict[str, torch.Tensor]]:
    if len(weights) != 3 or any(
        not math.isfinite(float(weight)) or weight < 0 for weight in weights
    ):
        raise ValueError("loss_weights must contain three finite non-negative values")
    radial_difference = current.two_body - target.two_body.to(current.two_body)
    if radial_smoothing:
        grid = current.grids["two_body"]
        if grid.numel() < 2:
            raise ValueError("radial smoothing requires at least two radial bins")
        width = float((grid[1] - grid[0]).detach())
        radius = math.ceil(4.0 * radial_smoothing / width)
        offsets = torch.arange(
            -radius, radius + 1, dtype=grid.dtype, device=grid.device
        ) * width
        kernel = torch.exp(-0.5 * (offsets / radial_smoothing).square())
        kernel = (kernel / kernel.sum()).reshape(1, 1, -1)
        radial_difference = torch.nn.functional.conv1d(
            radial_difference.T.unsqueeze(1), kernel, padding=radius
        ).squeeze(1).T
    losses = {
        "two_body": radial_difference.square().mean(),
        "three_body": (current.three_body - target.three_body.to(current.three_body)).square().mean(),
        "four_body": (current.four_body - target.four_body.to(current.four_body)).square().mean(),
    }
    total = sum(weight * loss for weight, loss in zip(weights, losses.values()))
    return total, losses


class FingerprintCalculator(Calculator):
    """ASE calculator for fingerprint mismatch (not a physical energy model)."""

    implemented_properties = ["energy", "forces"]

    def __init__(
        self,
        descriptor: DifferentiableFingerprint,
        target: FingerprintResult,
        *,
        loss_weights: Sequence[float] = (1.0, 1.0, 1.0),
        species_probabilities: torch.Tensor | np.ndarray | None = None,
        **kwargs: Any,
    ):
        super().__init__(**kwargs)
        self.descriptor = descriptor
        self.target = target
        self.loss_weights = tuple(float(value) for value in loss_weights)
        if len(self.loss_weights) != 3 or any(
            not math.isfinite(value) or value < 0 for value in self.loss_weights
        ):
            raise ValueError("loss_weights must contain three finite non-negative values")
        self._species_probabilities = species_probabilities

    def set_species_probabilities(
        self, probabilities: torch.Tensor | np.ndarray | None
    ) -> None:
        self._species_probabilities = probabilities
        self.reset()

    def calculate(
        self,
        atoms: Atoms | None = None,
        properties: Sequence[str] = ("energy", "forces"),
        system_changes: Sequence[str] = all_changes,
    ) -> None:
        super().calculate(atoms, properties, system_changes)
        if atoms is None:
            raise ValueError("ASE atoms are required")
        probabilities = (
            self.descriptor.probabilities_from_atoms(atoms)
            if self._species_probabilities is None
            else torch.as_tensor(self._species_probabilities, dtype=torch.float64)
        )
        positions = torch.tensor(atoms.positions, dtype=torch.float64, requires_grad=True)
        result = self.descriptor(
            positions,
            torch.as_tensor(np.asarray(atoms.cell), dtype=torch.float64),
            probabilities,
            atoms.pbc,
        )
        energy, _ = _fingerprint_loss(result, self.target, self.loss_weights)
        gradient = torch.autograd.grad(energy, positions)[0]
        self.results = {
            "energy": float(energy.detach()),
            "forces": (-gradient).detach().cpu().numpy(),
        }


class FingerprintRelaxation:
    """Jointly optimise ASE positions and continuous element assignments."""

    def __init__(
        self,
        descriptor: DifferentiableFingerprint,
        target: FingerprintResult,
        *,
        target_composition: Mapping[str, int] | None = None,
        composition_constraint: float | str = 0.0,
        fixed_position_mask: Sequence[bool] | None = None,
        fixed_species_mask: Sequence[bool] | None = None,
        loss_weights: Sequence[float] = (1.0, 1.0, 1.0),
        position_learning_rate: float = 0.02,
        species_learning_rate: float = 0.05,
        temperature: float = 1.0,
        final_temperature: float = 0.2,
        entropy_weight: float = 0.0,
        repulsion_weight: float = 0.0,
        radial_smoothing: float = 0.0,
    ):
        if isinstance(composition_constraint, str):
            if composition_constraint.lower() != "exact":
                raise ValueError('composition_constraint string alias must be "exact"')
            composition_constraint = math.inf
        if composition_constraint < 0 or math.isnan(float(composition_constraint)):
            raise ValueError("composition_constraint must be non-negative")
        self.descriptor = descriptor
        self.target = target
        self.composition_constraint = float(composition_constraint)
        self.target_composition = dict(
            target_composition or target.composition or {}
        )
        if math.isinf(self.composition_constraint) and not self.target_composition:
            raise ValueError("exact composition requires target_composition metadata")
        self.fixed_position_mask = fixed_position_mask
        self.fixed_species_mask = fixed_species_mask
        self.loss_weights = tuple(float(value) for value in loss_weights)
        if len(self.loss_weights) != 3 or any(
            not math.isfinite(value) or value < 0 for value in self.loss_weights
        ):
            raise ValueError("loss_weights must contain three finite non-negative values")
        self.position_learning_rate = float(position_learning_rate)
        self.species_learning_rate = float(species_learning_rate)
        self.temperature = float(temperature)
        self.final_temperature = float(final_temperature)
        self.entropy_weight = float(entropy_weight)
        self.repulsion_weight = float(repulsion_weight)
        self.radial_smoothing = float(radial_smoothing)
        if not math.isfinite(self.radial_smoothing) or self.radial_smoothing < 0:
            raise ValueError("radial_smoothing must be finite and non-negative")
        if min(self.position_learning_rate, self.species_learning_rate) <= 0:
            raise ValueError("learning rates must be positive")
        if min(self.temperature, self.final_temperature) <= 0:
            raise ValueError("softmax temperatures must be positive")

    def _composition_counts(self, atoms: Atoms) -> torch.Tensor:
        unknown = set(self.target_composition) - set(self.descriptor.elements)
        if unknown:
            raise ValueError(f"target composition contains elements outside the descriptor: {sorted(unknown)}")
        raw_counts = [
            self.target_composition.get(element, 0) for element in self.descriptor.elements
        ]
        if any(
            not math.isfinite(float(count))
            or float(count) < 0
            or not float(count).is_integer()
            for count in raw_counts
        ):
            raise ValueError("target composition counts must be non-negative integers")
        counts = torch.tensor(
            [int(count) for count in raw_counts],
            dtype=torch.float64,
        )
        if (counts < 0).any() or int(counts.sum()) != len(atoms):
            raise ValueError("target composition counts must be non-negative and sum to atom count")
        return counts

    @staticmethod
    def _sinkhorn(logits: torch.Tensor, counts: torch.Tensor, temperature: float) -> torch.Tensor:
        if logits.shape[0] == 0:
            return logits
        active = counts > 0
        if not active.any():
            raise ValueError("exact composition leaves no atoms for variable species")
        log_kernel = logits[:, active] / temperature
        target = counts[active]
        log_target = target.log()
        column_scale = torch.zeros_like(log_target)
        for _ in range(100):
            row_scale = -torch.logsumexp(log_kernel + column_scale[None, :], dim=1)
            column_scale = log_target - torch.logsumexp(
                log_kernel + row_scale[:, None], dim=0
            )
        probabilities = torch.zeros_like(logits)
        probabilities[:, active] = torch.exp(
            log_kernel + row_scale[:, None] + column_scale[None, :]
        )
        return probabilities

    def _probabilities(
        self,
        atoms: Atoms,
        logits: torch.Tensor,
        fixed_species: torch.Tensor,
        temperature: float,
        exact_counts: torch.Tensor | None,
    ) -> torch.Tensor:
        n_atoms = len(atoms)
        fixed = self.descriptor.probabilities_from_atoms(atoms).to(logits)
        variable = ~fixed_species
        if exact_counts is None:
            probabilities = logits.softmax(dim=-1)
        else:
            remaining = exact_counts.to(logits) - fixed[fixed_species].sum(dim=0)
            if (remaining < -1e-8).any() or abs(float(remaining.sum()) - int(variable.sum())) > 1e-8:
                raise ValueError("fixed species assignments make target composition infeasible")
            probabilities = logits.new_zeros((n_atoms, logits.shape[1]))
            probabilities[variable] = self._sinkhorn(
                logits[variable], remaining, temperature
            )
        probabilities = probabilities.clone()
        probabilities[fixed_species] = fixed[fixed_species]
        return probabilities

    def _assign_species(
        self, probabilities: torch.Tensor, fixed_species: torch.Tensor, atoms: Atoms,
        exact_counts: torch.Tensor | None,
    ) -> np.ndarray:
        from scipy.optimize import linear_sum_assignment

        assignment = probabilities.argmax(dim=1).cpu().numpy()
        if exact_counts is None:
            return assignment
        fixed = self.descriptor.probabilities_from_atoms(atoms).to(probabilities)
        remaining = exact_counts.to(probabilities) - fixed[fixed_species].sum(dim=0)
        variable_indices = torch.nonzero(~fixed_species, as_tuple=False).flatten().cpu().numpy()
        slots = np.repeat(
            np.arange(len(self.descriptor.elements)),
            remaining.round().to(torch.long).cpu().numpy(),
        )
        if len(slots) != len(variable_indices):
            raise ValueError("exact composition counts cannot be assigned to variable atoms")
        cost = -probabilities[variable_indices][:, slots].detach().cpu().numpy()
        rows, columns = linear_sum_assignment(cost)
        assignment[variable_indices[rows]] = slots[columns]
        assignment[fixed_species.cpu().numpy()] = fixed[fixed_species].argmax(dim=1).cpu().numpy()
        return assignment

    def run(
        self,
        initial_atoms: Atoms,
        *,
        optimise_positions: bool = True,
        optimise_species: bool = True,
        max_steps: int = 500,
        trajectory_interval: int = 1,
        final_refine_steps: int = 20,
        tolerance: float = 1e-8,
        callback: Any = None,
        verbose: bool = False,
    ) -> RelaxationResult:
        if len(initial_atoms) == 0:
            raise ValueError("cannot relax an empty structure")
        if max_steps < 1 or trajectory_interval < 1:
            raise ValueError("max_steps and trajectory_interval must be positive")
        if final_refine_steps < 0 or tolerance < 0 or not math.isfinite(tolerance):
            raise ValueError("final_refine_steps and tolerance must be non-negative")
        atoms = initial_atoms.copy()
        n_atoms, n_elements = len(atoms), len(self.descriptor.elements)
        fixed_positions = torch.as_tensor(
            self.fixed_position_mask if self.fixed_position_mask is not None else np.zeros(n_atoms),
            dtype=torch.bool,
        )
        fixed_species = torch.as_tensor(
            self.fixed_species_mask if self.fixed_species_mask is not None else np.zeros(n_atoms),
            dtype=torch.bool,
        )
        if fixed_positions.shape != (n_atoms,) or fixed_species.shape != (n_atoms,):
            raise ValueError("fixed masks must be Boolean arrays of length n_atoms")
        if not optimise_positions:
            fixed_positions[:] = True
        if not optimise_species:
            fixed_species[:] = True
        exact_counts = (
            self._composition_counts(atoms)
            if math.isinf(self.composition_constraint)
            else None
        )
        start_probabilities = self.descriptor.probabilities_from_atoms(atoms)
        logits = torch.nn.Parameter(
            torch.where(start_probabilities > 0, 2.0, 0.0).to(torch.float64)
        )
        position_parameter = torch.nn.Parameter(
            torch.as_tensor(atoms.positions, dtype=torch.float64).clone()
        )
        optimizer = torch.optim.Adam(
            [
                {"params": [position_parameter], "lr": self.position_learning_rate},
                {"params": [logits], "lr": self.species_learning_rate},
            ]
        )
        trajectory = [atoms.copy()]
        probability_history: list[torch.Tensor] = [start_probabilities.clone()]
        history: list[dict[str, float]] = []
        prior_positions = position_parameter.detach().clone()
        converged = False
        stable_steps = 0

        def log_step(phase: str, step: int, loss: float, gradient: torch.Tensor | None,
                     species_gradient: torch.Tensor | None, displacement: torch.Tensor) -> None:
            if not verbose:
                return
            fmax = 0.0 if gradient is None else float(
                torch.linalg.vector_norm(gradient, dim=1).max()
            )
            species_max = 0.0 if species_gradient is None else float(
                torch.linalg.vector_norm(species_gradient, dim=1).max()
            )
            max_step = float(torch.linalg.vector_norm(displacement, dim=1).max())
            print(
                f"{phase:>6} {step:5d} {datetime.now():%H:%M:%S} "
                f"{loss:12.6e} {fmax:12.6e} {species_max:12.6e} {max_step:10.4e}",
                flush=True,
            )

        if verbose:
            print(
                f"{'Phase':>6} {'Step':>5} {'Time':>8} {'Loss':>12} "
                f"{'Fmax':>12} {'SpeciesGrad':>12} {'MaxStep':>10}",
                flush=True,
            )

        for step in range(max_steps):
            fraction = step / max(1, max_steps - 1)
            temp = self.temperature * (self.final_temperature / self.temperature) ** fraction
            optimizer.zero_grad(set_to_none=True)
            position_parameter.grad = None
            probabilities = self._probabilities(
                atoms, logits, fixed_species, temp, exact_counts
            )
            evaluated_positions = torch.where(
                fixed_positions[:, None], prior_positions, position_parameter
            )
            current = self.descriptor(
                evaluated_positions,
                torch.as_tensor(np.asarray(atoms.cell), dtype=torch.float64),
                probabilities,
                atoms.pbc,
            )
            total, components = _fingerprint_loss(
                current, self.target, self.loss_weights, self.radial_smoothing
            )
            if self.composition_constraint > 0 and not math.isinf(self.composition_constraint):
                predicted = probabilities.sum(dim=0) / n_atoms
                desired = self._composition_counts(atoms).to(predicted) / n_atoms
                composition_loss = (predicted - desired).square().sum()
                total = total + self.composition_constraint * composition_loss
                components["composition"] = composition_loss
            if self.entropy_weight:
                entropy = -(probabilities * probabilities.clamp_min(1e-15).log()).sum(dim=1).mean()
                total = total + self.entropy_weight * entropy
                components["entropy"] = entropy
            if self.repulsion_weight:
                ids_i, ids_j, shifts = self.descriptor._topology.pairs(
                    evaluated_positions,
                    torch.as_tensor(np.asarray(atoms.cell), dtype=torch.float64),
                    atoms.pbc,
                )
                vectors = evaluated_positions[ids_j] - evaluated_positions[ids_i] + shifts.to(
                    evaluated_positions.dtype
                ) @ torch.as_tensor(np.asarray(atoms.cell), dtype=torch.float64)
                distances = torch.linalg.vector_norm(vectors, dim=1)
                radii = self.descriptor._radii.to(probabilities)
                minimum = 0.7 * (
                    probabilities[ids_i]
                    * (probabilities[ids_j] @ radii.T)
                ).sum(dim=1)
                repulsion = torch.relu(minimum - distances).square().mean()
                total = total + self.repulsion_weight * repulsion
                components["repulsion"] = repulsion
            total.backward()
            if not optimise_positions:
                position_parameter.grad = None
            elif position_parameter.grad is not None:
                position_parameter.grad[fixed_positions] = 0.0
            if not optimise_species:
                logits.grad = None
            previous_positions = position_parameter.detach().clone() if verbose else None
            optimizer.step()
            with torch.no_grad():
                position_parameter[fixed_positions] = prior_positions[fixed_positions]
            record = {key: float(value.detach()) for key, value in components.items()}
            record["total"] = float(total.detach())
            history.append(record)
            if verbose:
                log_step(
                    "Joint", step, record["total"], position_parameter.grad,
                    logits.grad, position_parameter.detach() - previous_positions,
                )
            save_frame = step % trajectory_interval == 0 or step == max_steps - 1
            if callback is not None or save_frame:
                frame = atoms.copy()
                frame.positions = evaluated_positions.detach().cpu().numpy()
                frame_probabilities = probabilities.detach().clone()
                assignment = self._assign_species(
                    frame_probabilities, fixed_species, atoms, exact_counts
                )
                frame.set_chemical_symbols(
                    [self.descriptor.elements[index] for index in assignment]
                )
            if callback is not None:
                callback(step, frame.copy(), frame_probabilities.clone(), record)
            if save_frame:
                trajectory.append(frame)
                probability_history.append(frame_probabilities)
            if len(history) > 1 and tolerance > 0:
                previous_loss, current_loss = history[-2]["total"], history[-1]["total"]
                scale = max(abs(previous_loss), abs(current_loss), torch.finfo(torch.float64).tiny)
                stable_steps = (
                    stable_steps + 1
                    if abs(previous_loss - current_loss) <= tolerance * scale
                    else 0
                )
                if stable_steps >= 5:
                    converged = True
                    break

        with torch.no_grad():
            probabilities = self._probabilities(
                atoms, logits, fixed_species, self.final_temperature, exact_counts
            )
            final_positions = torch.where(
                fixed_positions[:, None], prior_positions, position_parameter
            )
            assignment = self._assign_species(probabilities, fixed_species, atoms, exact_counts)
        atoms.set_positions(final_positions.detach().cpu().numpy())
        atoms.set_chemical_symbols([self.descriptor.elements[index] for index in assignment])
        final_probabilities = torch.nn.functional.one_hot(
            torch.as_tensor(assignment), num_classes=n_elements
        ).to(torch.float64)

        refine_positions = torch.nn.Parameter(torch.tensor(atoms.positions, dtype=torch.float64))
        refine_optimizer = torch.optim.Adam([refine_positions], lr=self.position_learning_rate)
        best_positions = refine_positions.detach().clone()
        best_loss = math.inf
        for refine_step in range(final_refine_steps if optimise_positions else 0):
            refine_optimizer.zero_grad(set_to_none=True)
            final_fp = self.descriptor(
                refine_positions,
                torch.as_tensor(np.asarray(atoms.cell), dtype=torch.float64),
                final_probabilities,
                atoms.pbc,
            )
            final_loss, _ = _fingerprint_loss(
                final_fp, self.target, self.loss_weights, self.radial_smoothing
            )
            loss_value = float(final_loss.detach())
            if loss_value < best_loss:
                best_loss = loss_value
                best_positions = refine_positions.detach().clone()
            final_loss.backward()
            refine_positions.grad[fixed_positions] = 0.0
            previous_positions = refine_positions.detach().clone() if verbose else None
            refine_optimizer.step()
            with torch.no_grad():
                refine_positions[fixed_positions] = prior_positions[fixed_positions]
            if verbose:
                log_step(
                    "Refine", refine_step, loss_value, refine_positions.grad, None,
                    refine_positions.detach() - previous_positions,
                )

        final_fp = self.descriptor(
            refine_positions,
            torch.as_tensor(np.asarray(atoms.cell), dtype=torch.float64),
            final_probabilities,
            atoms.pbc,
        )
        final_loss, _ = _fingerprint_loss(
            final_fp, self.target, self.loss_weights, self.radial_smoothing
        )
        if float(final_loss.detach()) > best_loss:
            with torch.no_grad():
                refine_positions.copy_(best_positions)
            final_fp = self.descriptor(
                refine_positions,
                torch.as_tensor(np.asarray(atoms.cell), dtype=torch.float64),
                final_probabilities,
                atoms.pbc,
            )
            final_loss, _ = _fingerprint_loss(
                final_fp, self.target, self.loss_weights, self.radial_smoothing
            )
        atoms.set_positions(refine_positions.detach().cpu().numpy())
        composition = {
            element: atoms.get_chemical_symbols().count(element)
            for element in self.descriptor.elements
        }
        trajectory.append(atoms.copy())
        probability_history.append(final_probabilities.detach().clone())
        if fixed_positions.any() and not np.array_equal(
            atoms.positions[fixed_positions.numpy()], initial_atoms.positions[fixed_positions.numpy()]
        ):
            raise RuntimeError("fixed-position constraint was violated")
        return RelaxationResult(
            atoms,
            final_fp,
            float(final_loss.detach()),
            converged,
            len(history),
            probabilities.detach().clone(),
            composition,
            history,
            trajectory,
            probability_history,
        )
