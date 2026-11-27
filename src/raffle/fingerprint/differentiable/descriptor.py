"""Differentiable RAFFLE fingerprints and fixed-cell inverse design.

The descriptor produces normalized, species-resolved distance and angle
histograms. The calculator exposes fingerprint mismatch to ASE, while the
relaxation class jointly optimizes positions and continuous species
probabilities. Neither objective is a physical energy or force field.
"""

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
from ase.geometry import complete_cell, find_mic
from ase.neighborlist import primitive_neighbor_list


_FOUR_BODY_CHECKPOINT_MIN_COMBINATIONS = 1024
_FOUR_BODY_UNCHECKPOINTED_BUDGET = 2_500_000
# Per-device offset vectors reused by every sparse Gaussian projection.
_SPARSE_OFFSET_CACHE: dict[tuple[int, torch.device], torch.Tensor] = {}


def _sparse_offsets(radius: int, device: torch.device) -> torch.Tensor:
    key = (radius, device)
    offsets = _SPARSE_OFFSET_CACHE.get(key)
    if offsets is None:
        offsets = torch.arange(-radius, radius + 1, device=device)
        _SPARSE_OFFSET_CACHE[key] = offsets
    return offsets


@dataclass
class FingerprintResult:
    """Species-resolved, bin-by-channel fingerprint tensors.

    Two-body channels are unordered element pairs in alphabetical order;
    three- and four-body channels correspond to the central element. ``grids``
    and ``unnormalized`` use the keys ``two_body``, ``three_body`` and
    ``four_body``. Disabled orders have zero arrays and are identified by
    ``body_orders``. ``composition`` is populated only when generated from ASE.
    """

    two_body: torch.Tensor
    three_body: torch.Tensor
    four_body: torch.Tensor
    unnormalized: dict[str, torch.Tensor]
    grids: dict[str, torch.Tensor]
    elements: tuple[str, ...]
    composition: dict[str, int] | None = None
    body_orders: tuple[int, ...] = (2, 3, 4)

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
    """Output from a coupled position/species relaxation.

    ``fingerprint`` and ``loss`` describe the final *discrete* structure;
    ``species_probabilities`` contains the final soft probabilities before
    assignment. ``loss_history`` describes joint, pre-update evaluations, not
    the final refinement. Trajectory and probability histories contain the
    initial state, saved joint states, and final discrete state.
    """

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
    """Cache cutoff-plus-skin candidate pairs, sorted by central atom.

    The cache is rebuilt when the cell, periodicity or atom count changes, or
    an atom moves more than half the skin since the last rebuild. Pair
    membership and image shifts are discrete; displacement vectors and
    distances are recomputed from PyTorch positions for every evaluation.
    Crossing a periodic boundary adjusts cached image shifts without rebuilding.
    """

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
        self.counts: tuple[int, ...] = ()
        self.reference_cutoff: float | None = None
        self.reference_skin: float | None = None
        self.wraps = np.empty((0, 3), dtype=np.int64)
        self._device_pairs: dict[
            torch.device, tuple[torch.Tensor, torch.Tensor, torch.Tensor]
        ] = {}

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
                or not np.array_equal(self.reference_cell, cell_np)
                or self.reference_cutoff != self.cutoff
                or self.reference_skin != self.skin
            )
        if not rebuild:
            delta = pos - self.reference_positions
            displacement = np.linalg.norm(delta, axis=1).max(initial=0.0)
            wraps = np.zeros_like(delta, dtype=np.int64)
            if pbc_np.any() and displacement > self.skin / 2.0:
                delta_mic, lengths = find_mic(delta, cell_np, pbc=pbc_np)
                displacement = lengths.max(initial=0.0)
                wraps = np.rint(
                    np.linalg.solve(complete_cell(cell_np).T, (delta - delta_mic).T).T
                ).astype(np.int64)
                wraps[:, ~pbc_np] = 0
            rebuild = displacement > self.skin / 2.0
            if not rebuild and not np.array_equal(wraps, self.wraps):
                wrap_delta = wraps - self.wraps
                self.shifts = self.shifts + wrap_delta[self.i] - wrap_delta[self.j]
                order = np.lexsort(
                    (self.shifts[:, 2], self.shifts[:, 1], self.shifts[:, 0],
                     self.j, np.any(self.shifts != 0, axis=1), self.i)
                )
                self.i, self.j, self.shifts = self.i[order], self.j[order], self.shifts[order]
                self.wraps = wraps
                self._device_pairs.clear()
        if rebuild:
            self.i, self.j, self.shifts = primitive_neighbor_list(
                "ijS",
                pbc_np,
                cell_np,
                pos,
                self.cutoff + self.skin,
                self_interaction=False,
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
            # Sorted centers let each atom's pairs be sliced without N scans.
            self.counts = tuple(np.bincount(self.i, minlength=len(pos)).tolist())
            self.reference_positions = pos.copy()
            self.reference_cell = cell_np.copy()
            self.reference_pbc = pbc_np.copy()
            self.reference_cutoff = self.cutoff
            self.reference_skin = self.skin
            self.wraps = np.zeros_like(pos, dtype=np.int64)
            self._device_pairs.clear()
        device = positions.device
        if device not in self._device_pairs:
            self._device_pairs[device] = (
                torch.as_tensor(self.i, dtype=torch.long, device=device),
                torch.as_tensor(self.j, dtype=torch.long, device=device),
                torch.as_tensor(self.shifts, dtype=torch.long, device=device),
            )
        return self._device_pairs[device]


class _SparseGaussianProjection(torch.autograd.Function):
    """Recompute local Gaussian support in backward instead of retaining it."""

    @staticmethod
    def forward(ctx, values, weights, grid, eta, width, minimum, radius):
        ctx.save_for_backward(values, weights, grid)
        ctx.eta, ctx.width, ctx.minimum, ctx.radius = eta, width, minimum, radius
        offsets = _sparse_offsets(radius, values.device)
        ctx.chunk = min(8192, max(1, 2_097_152 // offsets.numel()))
        result = weights.new_zeros(grid.shape)
        for start in range(0, values.numel(), ctx.chunk):
            x, weight = values[start:start + ctx.chunk], weights[start:start + ctx.chunk]
            nearest = torch.floor((x - minimum) / width + 0.5).to(torch.long)
            index = nearest[:, None] + offsets[None, :]
            valid = (index >= 0) & (index < grid.numel())
            safe = index.clamp(0, grid.numel() - 1)
            gaussian = torch.exp(-eta * (x[:, None] - grid[safe]).square()) * valid
            result.index_add_(0, safe.reshape(-1), (gaussian * weight[:, None]).reshape(-1))
        return result

    @staticmethod
    def backward(ctx, grad_output):
        values, weights, grid = ctx.saved_tensors
        offsets = _sparse_offsets(ctx.radius, values.device)
        grad_values, grad_weights = [], []
        for start in range(0, values.numel(), ctx.chunk):
            x, weight = values[start:start + ctx.chunk], weights[start:start + ctx.chunk]
            nearest = torch.floor((x - ctx.minimum) / ctx.width + 0.5).to(torch.long)
            index = nearest[:, None] + offsets[None, :]
            valid = (index >= 0) & (index < grid.numel())
            safe = index.clamp(0, grid.numel() - 1)
            difference = x[:, None] - grid[safe]
            gaussian = torch.exp(-ctx.eta * difference.square()) * valid
            weighted_gradient = grad_output[safe] * gaussian
            grad_weights.append(weighted_gradient.sum(dim=1))
            grad_values.append(
                (weighted_gradient * weight[:, None] * (-2.0 * ctx.eta * difference)).sum(dim=1)
            )
        return torch.cat(grad_values), torch.cat(grad_weights), None, None, None, None, None


class DifferentiableFingerprint:
    """Compute differentiable 2-, 3- and 4-body RAFFLE-style histograms.

    ``elements`` supplies the possible species, sorted alphabetically for
    channel ordering. ``parameters`` may be a RAFFLE distribution container
    or a mapping of ``width``, ``sigma``, ``cutoff_min``, ``cutoff_max``,
    ``radius_distance_tol`` and optionally ``bond_radii`` and
    ``gaussian_tolerance``. The latter is the per-Gaussian kernel amplitude
    bound outside the retained support and defaults to ``exp(-16)``. Explicit
    ``bond_radii`` override those in the parameters; missing pairs use the
    mean of ASE covalent radii. The four distance tolerances define the
    lower/upper radial windows for the first two neighbours and for the
    fourth-body neighbour, respectively.

    ``skin`` controls the cached neighbor-list margin in Angstrom.
    ``body_orders`` selects which histograms to compute; disabled orders
    retain correctly shaped zero tensors. ``dtype`` sets the floating-point
    precision for ASE target generation, the calculator and relaxation. Direct
    calls use the dtype and device of floating-point position inputs instead.
    ``parallel_workers`` optionally enables the compiled CPU four-body
    projection with the specified number of workers.

    Two-body channels histogram distances, with per-center species-weighted
    neighbor normalization, inverse-square radial weighting and an edge taper.
    Three-body channels histogram angles between neighbor vectors; four-body
    channels histogram angles between the normals of two neighbor planes.
    Angular contributions use distance windows and inverse-square weights.
    Channels with contributions above the normalization threshold are
    normalized independently.

    Positions and species probabilities participate in autograd, but
    neighbor membership, bin support and angular validity masks are discrete.
    Gradients are therefore piecewise differentiable, not smooth across every
    cutoff or degeneracy. Exact agreement with Fortran is not guaranteed for
    every geometry.
    """

    def __init__(
        self,
        elements: Sequence[str],
        parameters: Any = None,
        *,
        bond_radii: Mapping[tuple[str, str], float] | None = None,
        skin: float = 0.5,
        body_orders: Sequence[int] = (2, 3, 4),
        dtype: torch.dtype = torch.float64,
        parallel_workers: int | None = None,
    ):
        if dtype not in (torch.float32, torch.float64):
            raise ValueError("dtype must be torch.float32 or torch.float64")
        if parallel_workers is not None and (
            isinstance(parallel_workers, bool)
            or not isinstance(parallel_workers, int)
            or parallel_workers < 1
            or parallel_workers > 2**31 - 1
        ):
            raise ValueError("parallel_workers must be a positive 32-bit integer or None")
        self.dtype = dtype
        self.parallel_workers = parallel_workers
        normalized_elements = tuple(sorted({str(element).strip() for element in elements}))
        if not normalized_elements or any(
            element not in atomic_numbers or atomic_numbers[element] == 0
            for element in normalized_elements
        ):
            raise ValueError("elements must contain valid chemical element symbols")
        self.elements = normalized_elements
        self.body_orders = tuple(sorted(set(body_orders)))
        if not self.body_orders or any(order not in (2, 3, 4) for order in self.body_orders):
            raise ValueError("body_orders must contain one or more of 2, 3 and 4")
        self.width = self._parameter(parameters, "width", [0.025, math.pi / 64, math.pi / 64])
        self.sigma = self._parameter(parameters, "sigma", [0.1, 0.1, 0.1])
        self.cutoff_min = self._parameter(parameters, "cutoff_min", [0.5, 0.0, 0.0])
        self.cutoff_max = self._parameter(parameters, "cutoff_max", [6.0, math.pi, math.pi])
        self.radius_distance_tol = self._parameter(
            parameters, "radius_distance_tol", [1.5, 2.5, 3.0, 6.0]
        )
        if isinstance(parameters, Mapping):
            gaussian_tolerance = parameters.get("gaussian_tolerance", math.exp(-16.0))
        else:
            gaussian_tolerance = getattr(parameters, "gaussian_tolerance", math.exp(-16.0))
        self.gaussian_tolerance = float(gaussian_tolerance)
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
        # Grids, radii, taper and radial windows are configuration constants.
        # Rebuilding them on every descriptor evaluation dominates small calls,
        # so they are memoized against the parameters that define them.
        self._constants_cache: dict[tuple[Any, ...], dict[str, Any]] = {}
        self._triu_cache: dict[
            tuple[int, torch.device], tuple[torch.Tensor, torch.Tensor]
        ] = {}

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
        if not math.isfinite(self.gaussian_tolerance) or not 0 < self.gaussian_tolerance < 1:
            raise ValueError("gaussian_tolerance must be finite and between 0 and 1")

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

    def _constants(self, dtype: torch.dtype, device: torch.device) -> dict[str, Any]:
        """Return memoized grids, radii, taper and radial windows."""
        # Include the tunable parameters so that mutating ``width`` and friends
        # invalidates the cached tensors instead of returning stale grids.
        key = (
            dtype,
            device,
            tuple(self.width),
            tuple(self.sigma),
            tuple(self.cutoff_min),
            tuple(self.cutoff_max),
            tuple(self.radius_distance_tol),
            self.gaussian_tolerance,
        )
        cached = self._constants_cache.get(key)
        if cached is not None:
            return cached
        grids = self._grids(dtype, device)
        radii = self._radii.to(dtype=dtype, device=device)
        n_elements = len(self.elements)
        grids2 = grids["two_body"]

        taper_max = self.cutoff_max[0] - 0.25
        taper_min = self.cutoff_min[0] + 0.25
        taper = torch.where(
            grids2 > taper_max,
            0.5 * (1.0 + torch.cos(math.pi * (grids2 - taper_max) / 0.25)),
            torch.ones_like(grids2),
        )
        taper = torch.where(
            grids2 < taper_min,
            0.5 * (1.0 + torch.cos(math.pi * (grids2 - taper_min) / 0.25)),
            taper,
        )

        # The neighbour windows depend on the central species only, so they can
        # be built once instead of for every atom.
        lower = grids2.new_full((n_elements,), self.cutoff_min[0])
        upper = grids2.new_full((n_elements,), self.cutoff_max[0])
        near_windows = [
            (
                torch.maximum(lower, radii[alpha] * self.radius_distance_tol[0]),
                torch.minimum(upper, radii[alpha] * self.radius_distance_tol[1]),
            )
            for alpha in range(n_elements)
        ]
        fourth_windows = [
            (
                torch.maximum(lower, radii[alpha] * self.radius_distance_tol[2]),
                torch.minimum(upper, radii[alpha] * self.radius_distance_tol[3]),
            )
            for alpha in range(n_elements)
        ]

        eta2 = 1.0 / (2.0 * self.sigma[0] ** 2)
        radial_radius = math.ceil(
            math.sqrt(-math.log(self.gaussian_tolerance) / eta2) / self.width[0]
        )
        pair_indices = [
            (a, b) for a in range(n_elements) for b in range(a, n_elements)
        ]
        radial_offsets = torch.arange(
            -radial_radius, radial_radius + 1, device=device
        )
        radial_chunk = min(
            2048, max(1, 524288 // (radial_offsets.numel() * len(pair_indices)))
        )
        constants = {
            "grids": grids,
            "radii": radii,
            "pair_indices": pair_indices,
            "taper": taper,
            "bins2": torch.arange(grids2.numel(), dtype=dtype, device=device),
            "near_windows": near_windows,
            "fourth_windows": fourth_windows,
            "eta2": eta2,
            "radial_radius": radial_radius,
            "radial_offsets": radial_offsets,
            "radial_chunk": radial_chunk,
            "radial_sparse": 2 * radial_radius + 1 < grids2.numel(),
        }
        self._constants_cache[key] = constants
        return constants

    def _triu_indices(
        self, n: int, device: torch.device
    ) -> tuple[torch.Tensor, torch.Tensor]:
        key = (n, device)
        cached = self._triu_cache.get(key)
        if cached is None:
            cached = torch.triu_indices(n, n, offset=1, device=device)
            self._triu_cache[key] = cached
        return cached

    @staticmethod
    def _project(
        values: torch.Tensor,
        weights: torch.Tensor,
        denominator: torch.Tensor,
        grid: torch.Tensor,
        sigma: float,
        width: float,
        cutoff_min: float,
        tolerance: float = math.exp(-16.0),
    ) -> torch.Tensor:
        """Project weighted values onto Gaussian bins within truncated support."""
        eta = 1.0 / (2.0 * sigma * sigma)
        if values.numel() == 0:
            return grid.new_zeros(grid.shape)
        radius = math.ceil(math.sqrt(-math.log(tolerance) / eta) / width)
        if (
            values.device.type == "cpu"
            and values.numel() >= 128
            and 4 * radius + 2 <= grid.numel()
        ):
            projected = _SparseGaussianProjection.apply(
                values, weights, grid, eta, width, cutoff_min, radius
            )
        else:
            nearest = torch.floor((values - cutoff_min) / width + 0.5)
            bin_index = torch.arange(grid.numel(), device=grid.device, dtype=grid.dtype)
            support = (bin_index[None, :] - nearest[:, None]).abs() <= radius
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
        *,
        active_species: Sequence[bool] | None = None,
    ) -> FingerprintResult:
        """Evaluate the descriptor without discarding position/species gradients.

        ``positions`` has shape ``(N, 3)`` in Angstrom, ``cell`` is the ASE
        3-by-3 cell in Angstrom, and ``species_probabilities`` has shape
        ``(N, S)`` with nonnegative rows summing to one. Periodic image shifts
        follow ``pbc``. Returned histograms have one row per grid bin.

        ``active_species`` optionally marks species that can carry probability
        mass. Channels for inactive species are known to be zero and are
        skipped; their entries in the returned histograms and their gradient
        contributions are therefore omitted. Inactive species must have zero
        probability mass. Omit it (the default) to evaluate every channel,
        including channels whose probabilities happen to be zero.
        """
        positions = torch.as_tensor(positions)
        if not positions.is_floating_point():
            positions = positions.to(self.dtype)
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

        constants = self._constants(dtype, device)
        grids = constants["grids"]
        pair_indices = constants["pair_indices"]
        if active_species is None:
            channels = range(n_elements)
        else:
            mask = np.asarray(active_species, dtype=bool)
            if mask.shape != (n_elements,):
                raise ValueError("active_species must contain one Boolean per element")
            if not mask.any():
                raise ValueError("active_species must select at least one element")
            if (~mask).any():
                inactive = torch.as_tensor(~mask, dtype=torch.bool, device=device)
                if bool((probabilities[:, inactive] != 0).any()):
                    raise ValueError(
                        "inactive species must have zero probability mass"
                    )
            channels = np.flatnonzero(mask)
        raw2 = positions.new_zeros((grids["two_body"].numel(), len(pair_indices)))
        raw3 = positions.new_zeros((grids["three_body"].numel(), n_elements))
        raw4 = positions.new_zeros((grids["four_body"].numel(), n_elements))

        self._topology.cutoff = max(self.cutoff_max[0], 1e-4)
        center, neighbor, shifts = self._topology.pairs(positions, cell, pbc)
        vectors = positions[neighbor] - positions[center] + shifts.to(dtype) @ cell
        distances = torch.linalg.vector_norm(vectors, dim=1)
        by_center = (
            torch.arange(center.numel(), device=device).split(self._topology.counts)
            if 3 in self.body_orders or 4 in self.body_orders else ()
        )
        rmin, rmax = self.cutoff_min[0], self.cutoff_max[0]
        radial = (distances >= rmin) & (distances <= rmax) & (distances > 1e-10)
        grids2 = grids["two_body"]
        eta2 = constants["eta2"]

        if 2 in self.body_orders:
            # Reuse each radial Gaussian across all species-pair channels.
            # Denominators must still be computed separately for each center.
            radial_ids = torch.nonzero(radial, as_tuple=False).flatten()
            radial_centers = center[radial_ids]
            radial_neighbors = neighbor[radial_ids]
            neighbor_probabilities = probabilities[radial_neighbors]
            # Each center's denominator is the sum of its neighbors'
            # probabilities for the neighbor element in that channel.
            denominators = probabilities.new_zeros((n_atoms, n_elements)).index_add(
                0, radial_centers, neighbor_probabilities
            ).clamp_min(torch.finfo(dtype).tiny)
            bins = constants["bins2"]
            support_radius = constants["radial_radius"]
            sparse_projection = constants["radial_sparse"]
            if sparse_projection:
                offsets = constants["radial_offsets"]
                chunk = constants["radial_chunk"]
            else:
                chunk = 2048
            checkpoint_radial = (
                sparse_projection and n_atoms >= 512 and n_elements >= 4
                and torch.is_grad_enabled()
                and (distances.requires_grad or probabilities.requires_grad)
            )
            if checkpoint_radial:
                from torch.utils.checkpoint import checkpoint

            def project_radial(r, channel_weights, accumulator, width=self.width[0]):
                nearest = torch.floor((r - rmin) / width + 0.5)
                bin_indices = nearest.to(torch.long)[:, None] + offsets[None, :]
                valid_bins = (bin_indices >= 0) & (bin_indices < grids2.numel())
                safe_indices = bin_indices.clamp(0, grids2.numel() - 1)
                gaussian = torch.exp(-eta2 * (r[:, None] - grids2[safe_indices]).square())
                projected = gaussian * valid_bins.to(dtype)
                contributions = (
                    projected[:, :, None] * channel_weights[:, None, :]
                ).reshape(-1, len(pair_indices))
                return accumulator.index_add(
                    0, safe_indices.reshape(-1), contributions * math.sqrt(eta2 / math.pi)
                )

            for start in range(0, radial_ids.numel(), chunk):
                ids = radial_ids[start:start + chunk]
                centers = radial_centers[start:start + chunk]
                p_center = probabilities[centers]
                p_neighbor = neighbor_probabilities[start:start + chunk]
                denominator = denominators[centers]
                r = distances[ids]
                channel_weights = []
                for alpha, beta in pair_indices:
                    weights = p_center[:, alpha] * p_neighbor[:, beta] / denominator[:, beta]
                    if alpha != beta:
                        weights = weights + (
                            p_center[:, beta] * p_neighbor[:, alpha] / denominator[:, alpha]
                        )
                    channel_weights.append(weights)
                channel_weights = torch.stack(channel_weights, dim=1)
                if sparse_projection:
                    raw2 = (
                        checkpoint(
                            project_radial, r, channel_weights, raw2,
                            use_reentrant=False, preserve_rng_state=False,
                        )
                        if checkpoint_radial else project_radial(r, channel_weights, raw2)
                    )
                else:
                    nearest = torch.floor((r - rmin) / self.width[0] + 0.5)
                    support = (bins[None, :] - nearest[:, None]).abs() <= support_radius
                    gaussian = torch.exp(-eta2 * (r[:, None] - grids2[None, :]).square())
                    projected = gaussian * support.to(dtype)
                    raw2 = raw2 + projected.T @ channel_weights * math.sqrt(eta2 / math.pi)

        raw2 = raw2 / grids2[:, None].square()
        raw2 = raw2 * constants["taper"][:, None]
        norm2 = self._normalise(raw2)

        grids3 = grids["three_body"]
        near_cache: list[list[tuple[torch.Tensor, ...] | None]] | None = (
            [[None for _ in range(n_elements)] for _ in range(n_atoms)]
            if 3 in self.body_orders and 4 in self.body_orders else None
        )
        # Form unordered neighbor pairs around each center; distance windows
        # are averaged over the neighbors' continuous species probabilities.
        for atom in (range(n_atoms) if 3 in self.body_orders else ()):
            center_ids = by_center[atom]
            if center_ids.numel() < 2:
                continue
            for alpha in channels:
                ids = center_ids
                d = distances[ids]
                lower3, upper3 = constants["near_windows"][alpha]
                win3_by_species, valid3_by_species = self._cosine_window(
                    d, lower3, upper3
                )
                neighbor_probabilities = probabilities[neighbor[ids]]
                win3 = (win3_by_species * neighbor_probabilities).sum(dim=1)
                valid3 = (valid3_by_species * neighbor_probabilities).sum(dim=1)
                # Filter by geometric possibility, not the current species
                # probability, to retain potential species gradients.
                possible = valid3_by_species.any(dim=1)
                ids = ids[possible]
                d = d[possible]
                win3 = win3[possible]
                valid3 = valid3[possible]
                if ids.numel() < 2:
                    continue
                left, right = self._triu_indices(ids.numel(), device)
                pair_count_weights = valid3[left] * valid3[right]
                denominator = pair_count_weights.sum()
                if denominator.detach().item() == 0.0:
                    continue
                if near_cache is not None:
                    near_cache[atom][alpha] = (
                        ids, d, win3_by_species[possible], neighbor_probabilities[possible],
                        left, right, denominator,
                    )
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
                    self.gaussian_tolerance,
                )
                raw3[:, alpha] = raw3[:, alpha] + contribution

        norm3 = self._normalise(raw3)

        grids4 = grids["four_body"]
        parallel_four = (
            self.parallel_workers is not None and 4 in self.body_orders
            and device.type == "cpu" and dtype == torch.float64
        )
        if parallel_four:
            from ._parallel_four import _load_native, project

            _load_native()
            minimum_parallel_denominator = torch.finfo(dtype).tiny * max(
                1.0, 1.0 / (self.sigma[2] * math.sqrt(2.0 * math.pi))
            )
        environments = []
        checkpoint_four_body = (
            4 in self.body_orders and n_atoms >= 128 and n_elements >= 4
            and torch.is_grad_enabled()
            and (vectors.requires_grad or probabilities.requires_grad)
        )
        if checkpoint_four_body:
            from torch.utils.checkpoint import checkpoint
        uncheckpointed_combinations = 0
        sparse_four_body = (
            device.type == "cpu"
            and dtype == torch.float32
            and 4 * math.ceil(
                math.sqrt(-2.0 * math.log(self.gaussian_tolerance))
                * self.sigma[2] / self.width[2]
            ) + 2 <= grids4.numel()
        )

        # Combine each near-neighbor pair with a fourth-body neighbor and
        # project the angle between the resulting plane normals.
        for atom in (range(n_atoms) if 4 in self.body_orders else ()):
            center_ids = by_center[atom]
            if center_ids.numel() < 2:
                continue
            for alpha in channels:
                fourth_ids = center_ids
                if near_cache is not None:
                    cached = near_cache[atom][alpha]
                    if cached is None:
                        continue
                    (ids, d, win3_by_species, neighbor_probabilities,
                     left, right, near_denominator) = cached
                else:
                    ids = center_ids
                    d = distances[ids]
                    lower3, upper3 = constants["near_windows"][alpha]
                    win3_by_species, valid3_by_species = self._cosine_window(
                        d, lower3, upper3
                    )
                    neighbor_probabilities = probabilities[neighbor[ids]]
                    valid3 = (valid3_by_species * neighbor_probabilities).sum(dim=1)
                    possible3 = valid3_by_species.any(dim=1)
                    ids = ids[possible3]
                    d = d[possible3]
                    valid3 = valid3[possible3]
                    win3_by_species = win3_by_species[possible3]
                    neighbor_probabilities = neighbor_probabilities[possible3]
                    if ids.numel() < 2:
                        continue
                    left, right = self._triu_indices(ids.numel(), device)
                    near_denominator = (valid3[left] * valid3[right]).sum()

                fourth_d = distances[fourth_ids]
                lower4, upper4 = constants["fourth_windows"][alpha]
                win4_by_species, valid4_by_species = self._cosine_window(
                    fourth_d, lower4, upper4
                )
                fourth_probabilities = probabilities[neighbor[fourth_ids]]
                win4 = (win4_by_species * fourth_probabilities).sum(dim=1)
                valid4 = (valid4_by_species * fourth_probabilities).sum(dim=1)
                possible4 = valid4_by_species.any(dim=1)
                fourth_ids = fourth_ids[possible4]
                fourth_d = fourth_d[possible4]
                win4 = win4[possible4]
                valid4 = valid4[possible4]
                if fourth_ids.numel() == 0:
                    continue

                denominator = near_denominator * valid4.sum()
                denominator_value = denominator.detach().item()
                if denominator_value == 0.0:
                    continue
                # Enumerate only local triplets, never a dense N^4 tensor.
                chunk = 8192
                n_combinations = left.numel() * fourth_ids.numel()
                projected_combinations = n_combinations * (2 if n_elements > 1 else 1)
                checkpoint_environment = checkpoint_four_body and (
                    not sparse_four_body
                    or n_combinations < 128
                    or n_combinations >= _FOUR_BODY_CHECKPOINT_MIN_COMBINATIONS
                    or uncheckpointed_combinations + projected_combinations
                    > _FOUR_BODY_UNCHECKPOINTED_BUDGET
                )
                if checkpoint_four_body and not checkpoint_environment:
                    uncheckpointed_combinations += projected_combinations
                angle_epsilon = max(torch.finfo(dtype).eps * 10.0, 1e-12)
                species_weights = win3_by_species * neighbor_probabilities
                # All beta <= gamma share one plane orientation; beta > gamma
                # share the reverse. Prefix sums combine their weights before
                # projection, avoiding a separate histogram for each pair.
                if n_elements > 1:
                    upper_weights = species_weights.flip(1).cumsum(1).flip(1)
                    lower_weights = torch.cat(
                        (torch.zeros_like(species_weights[:, :1]),
                         species_weights.cumsum(1)[:, :-1]),
                        dim=1,
                    )
                else:
                    upper_weights = lower_weights = species_weights
                near_vectors = vectors[ids]
                fourth_vectors = vectors[fourth_ids]
                pair_normals = torch.linalg.cross(
                    near_vectors[left], near_vectors[right], dim=1
                )
                pair_norms = torch.linalg.vector_norm(pair_normals, dim=1)
                fourth_normals = torch.linalg.cross(
                    near_vectors[:, None, :].expand(-1, fourth_ids.numel(), -1),
                    fourth_vectors[None, :, :].expand(ids.numel(), -1, -1),
                    dim=2,
                )
                fourth_norms = torch.linalg.vector_norm(fourth_normals, dim=2)
                forward_pair_weights = (
                    species_weights[left] * upper_weights[right]
                ).sum(dim=1)
                reverse_pair_weights = (
                    (species_weights[left] * lower_weights[right]).sum(dim=1)
                    if n_elements > 1 else forward_pair_weights
                )
                if parallel_four and denominator_value >= minimum_parallel_denominator:
                    environments.append((
                        alpha, pair_normals, pair_norms, fourth_normals, fourth_norms,
                        d, fourth_d, win4, probabilities[atom, alpha],
                        forward_pair_weights, reverse_pair_weights, denominator, left, right,
                    ))
                    continue
                for start in range(0, n_combinations, chunk):
                    flat_index = torch.arange(
                        start, min(start + chunk, n_combinations), device=device
                    )
                    pair_index = flat_index // fourth_ids.numel()
                    j, k = left[pair_index], right[pair_index]
                    l = flat_index % fourth_ids.numel()
                    for reversed_order in range(2 if n_elements > 1 else 1):
                        pair_weights = (
                            reverse_pair_weights if reversed_order else forward_pair_weights
                        )

                        # Freeze loop-specific indices and settings for backward recomputation.
                        def project_chunk(
                            pair_normals, pair_norms, fourth_normals, fourth_norms,
                            d, fourth_d, win4, center_probability, pair_weights, denominator,
                            pair_index=pair_index, j=j, k=k, l=l, reversed_order=reversed_order,
                            angle_epsilon=angle_epsilon, sigma=self.sigma[2], width=self.width[2],
                            cutoff_min=self.cutoff_min[2], tolerance=self.gaussian_tolerance,
                        ):
                            normal_a = pair_normals[pair_index]
                            if reversed_order:
                                normal_a = -normal_a
                            second = j if reversed_order else k
                            normal_b = fourth_normals[second, l]
                            norm_a = pair_norms[pair_index]
                            norm_b = fourth_norms[second, l]
                            valid_angle = (norm_a > 1e-12) & (norm_b > 1e-12)
                            cosine = (normal_a * normal_b).sum(dim=1) / (
                                norm_a * norm_b
                            ).clamp_min(1e-24)
                            angles = torch.acos(
                                cosine.clamp(-1.0 + angle_epsilon, 1.0 - angle_epsilon)
                            )
                            common_weights = win4[l] * center_probability / (
                                d[j].square() * d[k].square() * fourth_d[l].square()
                            ).clamp_min(1e-24)
                            weights = (
                                pair_weights[pair_index] * common_weights * valid_angle.to(dtype)
                            )
                            return self._project(
                                angles, weights, denominator, grids4, sigma,
                                width, cutoff_min, tolerance,
                            )

                        arguments = (
                            pair_normals, pair_norms, fourth_normals, fourth_norms,
                            d, fourth_d, win4, probabilities[atom, alpha],
                            pair_weights, denominator,
                        )
                        contribution = (
                            checkpoint(
                                project_chunk, *arguments,
                                use_reentrant=False, preserve_rng_state=False,
                            )
                            if checkpoint_environment else project_chunk(*arguments)
                        )
                        raw4[:, alpha] = raw4[:, alpha] + contribution

        if parallel_four:
            raw4 = project(
                environments, grids4, self.sigma[2], self.width[2], self.cutoff_min[2],
                self.gaussian_tolerance, n_elements, self.parallel_workers, raw4,
            )
        norm4 = self._normalise(raw4)
        composition = None
        return FingerprintResult(
            norm2,
            norm3,
            norm4,
            {"two_body": raw2, "three_body": raw3, "four_body": raw4},
            # Copy the shared cached grids so callers cannot corrupt them.
            {name: value.clone() for name, value in grids.items()},
            self.elements,
            composition,
            self.body_orders,
        )

    def generate(self, atoms: Atoms) -> FingerprintResult:
        """Generate a reusable target fingerprint from an ASE structure."""
        probabilities = self.probabilities_from_atoms(atoms)
        present = set(atoms.get_chemical_symbols())
        result = self(
            torch.as_tensor(atoms.positions, dtype=self.dtype),
            torch.as_tensor(np.asarray(atoms.cell), dtype=self.dtype),
            probabilities,
            atoms.pbc,
            active_species=[element in present for element in self.elements],
        )
        result.composition = {
            element: atoms.get_chemical_symbols().count(element) for element in self.elements
        }
        return result

    def probabilities_from_atoms(self, atoms: Atoms) -> torch.Tensor:
        missing = sorted(set(atoms.get_chemical_symbols()) - set(self.elements))
        if missing:
            raise ValueError(f"structure contains elements not configured in descriptor: {missing}")
        result = torch.zeros((len(atoms), len(self.elements)), dtype=self.dtype)
        for index, symbol in enumerate(atoms.get_chemical_symbols()):
            result[index, self.elements.index(symbol)] = 1.0
        return result


def _fingerprint_loss(
    current: FingerprintResult,
    target: FingerprintResult,
    weights: Sequence[float] = (1.0, 1.0, 1.0),
    radial_smoothing: float = 0.0,
    angular_smoothing: float = 0.0,
    loss_reduction: str = "mean",
) -> tuple[torch.Tensor, dict[str, torch.Tensor]]:
    if len(weights) != 3 or any(
        not math.isfinite(float(weight)) or weight < 0 for weight in weights
    ):
        raise ValueError("loss_weights must contain three finite non-negative values")
    if loss_reduction not in ("mean", "density"):
        raise ValueError('loss_reduction must be "mean" or "density"')
    for order, weight in zip((2, 3, 4), weights):
        if weight and (order not in current.body_orders or order not in target.body_orders):
            raise ValueError(f"{order}-body loss requires that body order in both fingerprints")

    def component_loss(name: str, actual: torch.Tensor, expected: torch.Tensor,
                       smoothing: float) -> torch.Tensor:
        difference = actual - expected.to(actual)
        grid = current.grids[name]
        if smoothing or loss_reduction == "density":
            if grid.numel() < 2:
                raise ValueError("smoothing and density loss require at least two bins")
            width = float((grid[1] - grid[0]).detach())
        if smoothing:
            radius = math.ceil(4.0 * smoothing / width)
            offsets = torch.arange(
                -radius, radius + 1, dtype=grid.dtype, device=grid.device
            ) * width
            kernel = torch.exp(-0.5 * (offsets / smoothing).square())
            kernel = (kernel / kernel.sum()).reshape(1, 1, -1)
            difference = torch.nn.functional.conv1d(
                difference.T.unsqueeze(1), kernel, padding=radius
            ).squeeze(1).T
        if loss_reduction == "density":
            return difference.square().sum(dim=0).mean() / width
        return difference.square().mean()

    losses = {
        "two_body": component_loss(
            "two_body", current.two_body, target.two_body, radial_smoothing
        ),
        "three_body": component_loss(
            "three_body", current.three_body, target.three_body, angular_smoothing
        ),
        "four_body": component_loss(
            "four_body", current.four_body, target.four_body, angular_smoothing
        ),
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
        if self._species_probabilities is None:
            probabilities = self.descriptor.probabilities_from_atoms(atoms)
            present = set(atoms.get_chemical_symbols())
            active_species = [element in present for element in self.descriptor.elements]
        else:
            probabilities = torch.as_tensor(
                self._species_probabilities, dtype=self.descriptor.dtype
            )
            active_species = None
        positions = torch.tensor(atoms.positions, dtype=self.descriptor.dtype, requires_grad=True)
        result = self.descriptor(
            positions,
            torch.as_tensor(np.asarray(atoms.cell), dtype=self.descriptor.dtype),
            probabilities,
            atoms.pbc,
            active_species=active_species,
        )
        energy, _ = _fingerprint_loss(result, self.target, self.loss_weights)
        gradient = torch.autograd.grad(energy, positions)[0]
        self.results = {
            "energy": float(energy.detach()),
            "forces": (-gradient).detach().cpu().numpy(),
        }


class FingerprintRelaxation:
    """Optimize positions and species against a target fingerprint with Adam.

    The target supplies histograms, grids and optionally composition, not
    target coordinates. The candidate's cell and periodicity stay fixed.
    Generate the target with the same elements and grids as the descriptor:
    matching grid values and element identities are not validated by the loss.

    ``loss_weights`` weights the 2-, 3- and 4-body squared differences.
    ``loss_reduction="mean"`` averages over bins and channels; ``"density"``
    sums over bins, divides by bin width and averages over channels.
    ``radial_smoothing`` (Angstrom) and ``angular_smoothing`` (radians) apply
    Gaussian convolutions to fingerprint differences, not to the structures.

    ``composition_constraint`` is a composition temperature: zero or
    ``"exact"`` enforces exact counts, a finite positive value divides the
    squared fractional composition mismatch, and infinity, ``"none"`` or
    ``None`` disables the constraint. Exact mode uses Sinkhorn-constrained
    probabilities and an integer assignment. The counts come from
    ``target_composition`` or target metadata. This deterministic bias does
    not sample a Boltzmann distribution or enable thermal barrier crossing.
    Entropy and short-range repulsion are optional joint-loop penalties; both
    default to zero. The current repulsion uses 0.7 times the descriptor's pair radius
    as its onset, excludes zero-shift self-pairs, and is not included in the
    final position-only refinement. Periodic self-images are included. The
    soft penalty is not a hard minimum-distance guarantee.

    ``temperature`` and ``final_temperature`` control geometric annealing of
    Sinkhorn probabilities in exact mode; the unconstrained softmax does not
    use them. Position and species learning rates are independent; neither the
    objective nor its gradients are physical energies or forces.
    """

    def __init__(
        self,
        descriptor: DifferentiableFingerprint,
        target: FingerprintResult,
        *,
        target_composition: Mapping[str, int] | None = None,
        composition_constraint: float | str | None = None,
        loss_weights: Sequence[float] = (1.0, 1.0, 1.0),
        position_learning_rate: float = 0.02,
        species_learning_rate: float = 0.05,
        temperature: float = 1.0,
        final_temperature: float = 0.2,
        entropy_weight: float = 0.0,
        repulsion_weight: float = 0.0,
        radial_smoothing: float = 0.0,
        angular_smoothing: float = 0.0,
        loss_reduction: str = "mean",
    ):
        if composition_constraint is None:
            composition_constraint = math.inf
        elif isinstance(composition_constraint, str):
            aliases = {"exact": 0.0, "none": math.inf}
            alias = composition_constraint.strip().lower()
            if alias not in aliases:
                raise ValueError('composition_constraint string must be "exact" or "none"')
            composition_constraint = aliases[alias]
        composition_constraint = float(composition_constraint)
        if composition_constraint < 0 or math.isnan(composition_constraint):
            raise ValueError("composition_constraint must be non-negative")
        if (
            0 < composition_constraint < math.inf
            and not math.isfinite(1.0 / composition_constraint)
        ):
            raise ValueError('composition_constraint is too small; use 0 or "exact"')
        self.descriptor = descriptor
        self.target = target
        self.composition_constraint = composition_constraint
        self.target_composition = dict(
            target_composition or target.composition or {}
        )
        if not math.isinf(self.composition_constraint) and not self.target_composition:
            raise ValueError("composition constraint requires target_composition metadata")
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
        self.repulsion_alpha = 4.0
        self.radial_smoothing = float(radial_smoothing)
        if not math.isfinite(self.radial_smoothing) or self.radial_smoothing < 0:
            raise ValueError("radial_smoothing must be finite and non-negative")
        self.angular_smoothing = float(angular_smoothing)
        if not math.isfinite(self.angular_smoothing) or self.angular_smoothing < 0:
            raise ValueError("angular_smoothing must be finite and non-negative")
        if loss_reduction not in ("mean", "density"):
            raise ValueError('loss_reduction must be "mean" or "density"')
        self.loss_reduction = loss_reduction
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
        """Balance soft assignment rows and element counts in log space.

        Columns with zero requested count stay zero. Iteration stops when the
        maximum count error is below 1e-5, or fails explicitly after 1000
        iterations. The scaling iterations remain part of the autograd graph.
        """
        if logits.shape[0] == 0:
            return logits
        active = counts > 0
        if not active.any():
            raise ValueError("exact composition leaves no atoms for variable species")
        # Float32 reductions can exceed the absolute count tolerance for large
        # systems; balance in float64 before returning the requested dtype.
        log_kernel = logits[:, active].to(torch.float64) / temperature
        target = counts[active].to(torch.float64)
        log_target = target.log()
        column_scale = torch.zeros_like(log_target)
        for iteration in range(1000):
            row_scale = -torch.logsumexp(log_kernel + column_scale[None, :], dim=1)
            column_scale = log_target - torch.logsumexp(
                log_kernel + row_scale[:, None], dim=0
            )
            if iteration % 10 == 9:
                with torch.no_grad():
                    active_probabilities = (log_kernel + column_scale[None, :]).softmax(dim=1)
                    column_error = (active_probabilities.sum(dim=0) - target).abs().max()
                if float(column_error) < 1e-5:
                    break
        else:
            raise RuntimeError(
                "Sinkhorn composition constraint did not converge; "
                "try a higher temperature or lower species learning rate"
            )
        probabilities = torch.zeros_like(logits)
        probabilities[:, active] = (log_kernel + column_scale[None, :]).softmax(dim=1).to(logits.dtype)
        return probabilities

    def _probabilities(
        self,
        initial_probabilities: torch.Tensor,
        logits: torch.Tensor,
        fixed_species: torch.Tensor,
        temperature: float,
        exact_counts: torch.Tensor | None,
    ) -> torch.Tensor:
        n_atoms = len(initial_probabilities)
        fixed = initial_probabilities.to(logits)
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
        self, probabilities: torch.Tensor, fixed_species: torch.Tensor,
        initial_probabilities: torch.Tensor,
        exact_counts: torch.Tensor | None,
    ) -> np.ndarray:
        """Round probabilities, honoring fixed species and exact counts."""
        assignment = probabilities.argmax(dim=1).cpu().numpy()
        if exact_counts is None:
            return assignment
        from scipy.optimize import linear_sum_assignment

        fixed = initial_probabilities.to(probabilities)
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
        fixed_position_mask: Sequence[bool] | None = None,
        fixed_species_mask: Sequence[bool] | None = None,
        max_steps: int = 500,
        trajectory_interval: int = 1,
        final_refine_steps: int = 20,
        tolerance: float = 1e-8,
        method: str = "adam",
        exploration_temperature: float = 0.0,
        denoiser_steps: int = 0,
        denoising_position_scale: float = 0.3,
        denoising_species_scale: float = 1.0,
        random_seed: int | None = None,
        callback: Any = None,
        verbose: bool = False,
    ) -> RelaxationResult:
        """Run joint Adam updates, discretise species, and refine positions.

        ``optimise_positions`` and ``optimise_species`` can freeze all atoms
        in either variable; per-atom masks independently freeze selected
        positions and species for this run.
        The input ASE object is copied, and its cell is never optimized.

        ``method`` selects the original Adam updates, Adam followed by
        Metropolis trial moves, or overdamped Langevin updates. The separate
        ``exploration_temperature`` is in fingerprint-objective units, not
        kelvin. ``denoiser_steps`` adds a preliminary stage of annealed
        perturbations followed by Adam updates against the fingerprint loss;
        this is objective-guided denoising, not a trained diffusion model.
        Its position and logit noise scales are controlled separately.
        ``random_seed`` repeats the stochastic draws without changing the
        global PyTorch random state.

        At each joint step, logits become species probabilities by softmax or,
        in exact mode, Sinkhorn balancing at the annealed temperature. The
        descriptor and weighted loss are evaluated, optional composition,
        entropy and repulsion terms are added, and the selected method updates
        positions and logits. Fixed values are restored after each update.
        For Adam, five consecutive relative loss changes no larger than
        ``tolerance`` stop the joint loop; zero tolerance disables early
        stopping. Stochastic methods run all ``max_steps``.

        ``callback(step, frame, probabilities, record)`` receives the
        pre-update evaluated state and unweighted component losses plus the
        weighted total. Joint trajectory frames are saved at
        ``trajectory_interval`` and use the same evaluated state. The first
        frame preserves the input species, and the final frame contains the
        assigned discrete species and final positions.

        After the joint loop, soft probabilities are converted to discrete
        species. Exact mode uses ``final_temperature`` and a minimum-cost
        assignment; otherwise each atom takes its most probable species.
        Up to ``final_refine_steps`` further Adam updates optimize positions
        with those species fixed, retaining the best fingerprint loss
        encountered. This refinement uses the same smoothing and body weights,
        but not the joint composition, entropy or repulsion penalties.

        The returned ``loss`` and ``fingerprint`` describe the final discrete
        structure; ``loss_history`` describes joint pre-update states, and
        ``species_probabilities`` is the final soft state before rounding.
        ``n_steps`` includes denoising and joint steps; ``converged`` refers
        to the joint loop. Verbose gradients are objective gradients, not
        physical forces.
        """
        if len(initial_atoms) == 0:
            raise ValueError("cannot relax an empty structure")
        if max_steps < 1 or trajectory_interval < 1:
            raise ValueError("max_steps and trajectory_interval must be positive")
        if final_refine_steps < 0 or tolerance < 0 or not math.isfinite(tolerance):
            raise ValueError("final_refine_steps and tolerance must be non-negative")
        if method not in ("adam", "metropolis", "langevin"):
            raise ValueError('method must be "adam", "metropolis", or "langevin"')
        if not math.isfinite(exploration_temperature) or exploration_temperature < 0:
            raise ValueError("exploration_temperature must be finite and non-negative")
        if method == "adam" and exploration_temperature != 0:
            raise ValueError('exploration_temperature requires method="metropolis" or "langevin"')
        if isinstance(denoiser_steps, bool) or not isinstance(denoiser_steps, int) or denoiser_steps < 0:
            raise ValueError("denoiser_steps must be a non-negative integer")
        if any(not math.isfinite(value) or value < 0 for value in (
            denoising_position_scale, denoising_species_scale
        )):
            raise ValueError("denoising noise scales must be finite and non-negative")
        if random_seed is not None and (
            isinstance(random_seed, bool) or not isinstance(random_seed, int)
            or not 0 <= random_seed < 2**64
        ):
            raise ValueError("random_seed must be an integer from 0 to 2**64 - 1")
        atoms = initial_atoms.copy()
        n_atoms, n_elements = len(atoms), len(self.descriptor.elements)
        fixed_positions = torch.as_tensor(
            fixed_position_mask if fixed_position_mask is not None else np.zeros(n_atoms),
            dtype=torch.bool,
        ).clone()
        fixed_species = torch.as_tensor(
            fixed_species_mask if fixed_species_mask is not None else np.zeros(n_atoms),
            dtype=torch.bool,
        ).clone()
        if fixed_positions.shape != (n_atoms,) or fixed_species.shape != (n_atoms,):
            raise ValueError("fixed masks must be Boolean arrays of length n_atoms")
        if not optimise_positions:
            fixed_positions[:] = True
        if not optimise_species:
            fixed_species[:] = True
        exact_counts = (
            self._composition_counts(atoms)
            if self.composition_constraint == 0
            else None
        )
        # In exact mode species with a zero count never receive probability
        # mass, so their descriptor channels (and gradients) are structurally
        # zero and can be skipped.
        active_species = (
            [bool(count) for count in exact_counts.tolist()]
            if exact_counts is not None else None
        )
        start_probabilities = self.descriptor.probabilities_from_atoms(atoms)
        dtype = self.descriptor.dtype
        cell = torch.as_tensor(np.asarray(atoms.cell), dtype=dtype)
        target_fractions = (
            self._composition_counts(atoms) / n_atoms
            if 0 < self.composition_constraint < math.inf else None
        )
        # Initial labels are favored, but remain soft unless fixed or forced by
        # an exact composition constraint.
        logits = torch.nn.Parameter(
            torch.where(start_probabilities > 0, 2.0, 0.0).to(dtype)
        )
        position_parameter = torch.nn.Parameter(
            torch.as_tensor(atoms.positions, dtype=dtype).clone()
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
        generator = None
        if denoiser_steps or method in ("metropolis", "langevin"):
            generator = torch.Generator(device="cpu")
            if random_seed is None:
                generator.seed()
            else:
                generator.manual_seed(random_seed)

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

        def evaluate_state(position_values, species_logits, temperature):
            probabilities = self._probabilities(
                start_probabilities, species_logits, fixed_species, temperature, exact_counts
            )
            evaluated_positions = torch.where(
                fixed_positions[:, None], prior_positions, position_values
            )
            current = self.descriptor(
                evaluated_positions,
                cell,
                probabilities,
                atoms.pbc,
                active_species=active_species,
            )
            total, components = _fingerprint_loss(
                current, self.target, self.loss_weights, self.radial_smoothing,
                self.angular_smoothing, self.loss_reduction,
            )
            if target_fractions is not None:
                predicted = probabilities.sum(dim=0) / n_atoms
                desired = target_fractions.to(predicted)
                composition_loss = (predicted - desired).square().sum()
                penalty = composition_loss / self.composition_constraint
                if not torch.isfinite(penalty.detach()).item():
                    raise FloatingPointError(
                        'composition penalty overflowed; increase composition_constraint or use "exact"'
                    )
                total = total + penalty
                components["composition"] = composition_loss
            if self.entropy_weight:
                entropy = -(probabilities * probabilities.clamp_min(1e-15).log()).sum(dim=1).mean()
                total = total + self.entropy_weight * entropy
                components["entropy"] = entropy
            if self.repulsion_weight:
                ids_i, ids_j, shifts = self.descriptor._topology.pairs(
                    evaluated_positions,
                    cell,
                    atoms.pbc,
                )

                # Exclude zero-image self-interactions, but retain periodic self-images.
                interacting = (ids_i != ids_j) | shifts.ne(0).any(dim=1)
                ids_i = ids_i[interacting]
                ids_j = ids_j[interacting]
                shifts = shifts[interacting]

                vectors = (
                    evaluated_positions[ids_j]
                    - evaluated_positions[ids_i]
                    + shifts.to(evaluated_positions.dtype) @ cell
                )
                distances = torch.linalg.vector_norm(vectors, dim=1)

                # Pair-specific radii, including continuous species probabilities.
                radii = self.descriptor._radii.to(probabilities)

                pair_radius = (
                    probabilities[ids_i]
                    * (probabilities[ids_j] @ radii.T)
                ).sum(dim=1)

                # Convert the average atomic radii into pair separation.
                cutoff = 2.0 * pair_radius

                # Dimensionless fractional overlap: x = max(0, (R_cut - r) / R_cut).
                cutoff = cutoff.clamp_min(1e-8)
                overlap = ((cutoff - distances) / cutoff).clamp_min(0.0)

                # Smooth exponential repulsive potential.
                alpha = self.repulsion_alpha  # Suggested initial value: 4.0
                penalties = overlap.pow(3) * torch.exp(alpha * overlap)

                repulsion = penalties.mean() if penalties.numel() else penalties.sum()

                total = total + self.repulsion_weight * repulsion
                components["repulsion"] = repulsion
            return total, components, probabilities, evaluated_positions

        def add_noise(values, fixed, scale):
            if scale == 0 or bool(fixed.all()):
                return
            noise = torch.randn(values.shape, dtype=values.dtype, device=values.device, generator=generator)
            values.add_(noise * (~fixed)[:, None], alpha=scale)

        # Denoising is a separate warm-up stage; max_steps still counts joint steps.
        for step in range(denoiser_steps + max_steps):
            denoising = step < denoiser_steps
            joint_step = max(0, step - denoiser_steps)
            fraction = joint_step / max(1, max_steps - 1)
            temp = self.temperature * (self.final_temperature / self.temperature) ** fraction
            if denoising:
                noise_fraction = (denoiser_steps - step) / denoiser_steps
                with torch.no_grad():
                    add_noise(position_parameter, fixed_positions, denoising_position_scale * noise_fraction)
                    add_noise(logits, fixed_species, denoising_species_scale * noise_fraction)
            optimizer.zero_grad(set_to_none=True)
            position_parameter.grad = None
            total, components, probabilities, evaluated_positions = evaluate_state(
                position_parameter, logits, temp
            )
            total.backward()
            if (denoising or method != "adam") and any(
                gradient is not None and not torch.isfinite(gradient).all().item()
                for gradient in (position_parameter.grad, logits.grad)
            ):
                raise FloatingPointError("stochastic relaxation produced a non-finite gradient")
            if (
                target_fractions is not None
                and logits.grad is not None
                and not torch.isfinite(logits.grad.square()).all().item()
            ):
                raise FloatingPointError(
                    "species gradient is non-finite or too large for Adam; "
                    'check loss scales or increase composition_constraint'
                )
            if not optimise_positions:
                position_parameter.grad = None
            elif position_parameter.grad is not None:
                position_parameter.grad[fixed_positions] = 0.0
            if not optimise_species:
                logits.grad = None
            previous_positions = position_parameter.detach().clone() if verbose else None
            accepted = None
            if denoising or method != "langevin":
                optimizer.step()
            else:
                position_noise = math.sqrt(2.0 * exploration_temperature * self.position_learning_rate)
                species_noise = math.sqrt(2.0 * exploration_temperature * self.species_learning_rate)
                if not math.isfinite(position_noise) or not math.isfinite(species_noise):
                    raise FloatingPointError("Langevin noise scale is not finite")
                with torch.no_grad():
                    if position_parameter.grad is not None:
                        position_parameter.add_(position_parameter.grad, alpha=-self.position_learning_rate)
                    if logits.grad is not None:
                        logits.add_(logits.grad, alpha=-self.species_learning_rate)
                    add_noise(position_parameter, fixed_positions, position_noise)
                    add_noise(logits, fixed_species, species_noise)
                    if not torch.isfinite(position_parameter).all() or not torch.isfinite(logits).all():
                        raise FloatingPointError("Langevin update produced non-finite values")
            with torch.no_grad():
                position_parameter[fixed_positions] = prior_positions[fixed_positions]
                if not denoising and method == "metropolis":
                    if bool(fixed_positions.all()) and bool(fixed_species.all()):
                        accepted = True
                    else:
                        baseline, _, _, _ = evaluate_state(position_parameter, logits, temp)
                        proposed_positions = position_parameter.clone()
                        proposed_logits = logits.clone()
                        add_noise(proposed_positions, fixed_positions, self.position_learning_rate)
                        add_noise(proposed_logits, fixed_species, self.species_learning_rate)
                        proposed, _, _, _ = evaluate_state(proposed_positions, proposed_logits, temp)
                        baseline_value, proposed_value = float(baseline), float(proposed)
                        if not math.isfinite(baseline_value) or not math.isfinite(proposed_value):
                            raise FloatingPointError("Metropolis objective is not finite")
                        change = proposed_value - baseline_value
                        accepted = change <= 0 or (
                            exploration_temperature > 0
                            and float(torch.rand((), generator=generator))
                            < math.exp(-change / exploration_temperature)
                        )
                        if accepted:
                            position_parameter.copy_(proposed_positions)
                            logits.copy_(proposed_logits)
            record = {key: float(value.detach()) for key, value in components.items()}
            record["total"] = float(total.detach())
            if denoising:
                record["denoising"] = 1.0
            if accepted is not None:
                record["metropolis_accepted"] = float(accepted)
            history.append(record)
            if verbose:
                log_step(
                    "Denoise" if denoising else "Joint", step, record["total"], position_parameter.grad,
                    logits.grad, position_parameter.detach() - previous_positions,
                )
            save_frame = step % trajectory_interval == 0 or step == denoiser_steps + max_steps - 1
            if callback is not None or save_frame:
                # Keep frame positions and probabilities from the same
                # pre-update descriptor evaluation as the recorded loss.
                frame = atoms.copy()
                frame.positions = evaluated_positions.detach().cpu().numpy()
                frame.positions[fixed_positions.numpy()] = initial_atoms.positions[fixed_positions.numpy()]
                frame_probabilities = probabilities.detach().clone()
                assignment = self._assign_species(
                    frame_probabilities, fixed_species, start_probabilities, exact_counts
                )
                frame.set_chemical_symbols(
                    [self.descriptor.elements[index] for index in assignment]
                )
            if callback is not None:
                callback(step, frame.copy(), frame_probabilities.clone(), record)
            if save_frame:
                trajectory.append(frame)
                probability_history.append(frame_probabilities)
            if not denoising and method == "adam" and joint_step > 0 and tolerance > 0:
                previous_loss, current_loss = history[-2]["total"], history[-1]["total"]
                scale = max(abs(previous_loss), abs(current_loss), torch.finfo(dtype).tiny)
                stable_steps = (
                    stable_steps + 1
                    if abs(previous_loss - current_loss) <= tolerance * scale
                    else 0
                )
                if stable_steps >= 5:
                    converged = True
                    break

        # Project the final soft species onto a discrete, possibly constrained,
        # assignment before position-only refinement.
        with torch.no_grad():
            probabilities = self._probabilities(
                start_probabilities, logits, fixed_species, self.final_temperature, exact_counts
            )
            final_positions = torch.where(
                fixed_positions[:, None], prior_positions, position_parameter
            )
            assignment = self._assign_species(
                probabilities, fixed_species, start_probabilities, exact_counts
            )
        atoms.set_positions(final_positions.detach().cpu().numpy())
        atoms.set_chemical_symbols([self.descriptor.elements[index] for index in assignment])
        final_probabilities = torch.nn.functional.one_hot(
            torch.as_tensor(assignment), num_classes=n_elements
        ).to(dtype)

        refine_positions = torch.nn.Parameter(torch.tensor(atoms.positions, dtype=dtype))
        refine_optimizer = torch.optim.Adam([refine_positions], lr=self.position_learning_rate)
        best_positions = refine_positions.detach().clone()
        best_loss = math.inf
        # Species are fixed here; retain the best evaluated positions because
        # Adam's final update can increase the fingerprint loss.
        for refine_step in range(final_refine_steps if optimise_positions else 0):
            refine_optimizer.zero_grad(set_to_none=True)
            final_fp = self.descriptor(
                refine_positions,
                cell,
                final_probabilities,
                atoms.pbc,
                active_species=active_species,
            )
            final_loss, _ = _fingerprint_loss(
                final_fp, self.target, self.loss_weights, self.radial_smoothing,
                self.angular_smoothing, self.loss_reduction,
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

        with torch.no_grad():
            final_fp = self.descriptor(
                refine_positions,
                cell,
                final_probabilities,
                atoms.pbc,
                active_species=active_species,
            )
            final_loss, _ = _fingerprint_loss(
                final_fp, self.target, self.loss_weights, self.radial_smoothing,
                self.angular_smoothing, self.loss_reduction,
            )
            if float(final_loss) > best_loss:
                refine_positions.copy_(best_positions)
                final_fp = self.descriptor(
                    refine_positions,
                    cell,
                    final_probabilities,
                    atoms.pbc,
                    active_species=active_species,
                )
                final_loss, _ = _fingerprint_loss(
                    final_fp, self.target, self.loss_weights, self.radial_smoothing,
                    self.angular_smoothing, self.loss_reduction,
                )
        atoms.set_positions(refine_positions.detach().cpu().numpy())
        atoms.positions[fixed_positions.numpy()] = initial_atoms.positions[fixed_positions.numpy()]
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
