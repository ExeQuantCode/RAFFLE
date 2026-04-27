from __future__ import annotations

from dataclasses import dataclass
import math
import random
from typing import Callable, Iterable, Optional, Sequence, Tuple

import numpy as np
import torch
import torch.nn.functional as torch_functional
from torch import nn

from .gnn_fingerprint import GNNFingerprint
from .structure_metrics import wrap_atoms_to_unit_cell


ELEMENT_PROPERTIES = {
    "H": (1.0 / 100.0, 0.31),
    "He": (2.0 / 100.0, 0.28),
    "Li": (3.0 / 100.0, 1.28),
    "Be": (4.0 / 100.0, 0.96),
    "B": (5.0 / 100.0, 0.84),
    "C": (6.0 / 100.0, 0.76),
    "N": (7.0 / 100.0, 0.71),
    "O": (8.0 / 100.0, 0.66),
    "F": (9.0 / 100.0, 0.57),
    "Ne": (10.0 / 100.0, 0.58),
    "Na": (11.0 / 100.0, 1.66),
    "Mg": (12.0 / 100.0, 1.41),
    "Al": (13.0 / 100.0, 1.21),
    "Si": (14.0 / 100.0, 1.11),
    "P": (15.0 / 100.0, 1.07),
    "S": (16.0 / 100.0, 1.05),
    "Cl": (17.0 / 100.0, 1.02),
    "Ar": (18.0 / 100.0, 1.06),
    "K": (19.0 / 100.0, 2.03),
    "Ca": (20.0 / 100.0, 1.76),
    "Sc": (21.0 / 100.0, 1.70),
    "Ti": (22.0 / 100.0, 1.60),
    "V": (23.0 / 100.0, 1.53),
    "Cr": (24.0 / 100.0, 1.39),
    "Mn": (25.0 / 100.0, 1.39),
    "Fe": (26.0 / 100.0, 1.32),
    "Co": (27.0 / 100.0, 1.26),
    "Ni": (28.0 / 100.0, 1.24),
    "Cu": (29.0 / 100.0, 1.32),
    "Zn": (30.0 / 100.0, 1.22),
    "Ga": (31.0 / 100.0, 1.22),
    "Ge": (32.0 / 100.0, 1.20),
    "As": (33.0 / 100.0, 1.19),
    "Se": (34.0 / 100.0, 1.20),
    "Br": (35.0 / 100.0, 1.20),
    "Mo": (42.0 / 100.0, 1.54),
    "Ba": (56.0 / 100.0, 2.15),
    "W": (74.0 / 100.0, 1.62),
}

ARCHITECTURE_ALIASES = {
    "residual": "residual",
    "torch_gnn_residual": "residual",
    "gated": "gated",
    "torch_gnn_gated": "gated",
    "attention": "attention",
    "torch_gnn_attention": "attention",
    "attention_coupled": "attention_coupled",
    "torch_gnn_attention_coupled": "attention_coupled",
    "attention_conservative": "attention_conservative",
    "torch_gnn_attention_conservative": "attention_conservative",
}

FINGERPRINT_NEGATIVE_TAIL_BETA = 500.0


@dataclass(frozen=True)
class MultigraphTopology:
    symbols: Tuple[str, ...]
    species_index: np.ndarray
    atomic_numbers: np.ndarray
    covalent_radii: np.ndarray
    pair_index: np.ndarray
    pair_type_index: np.ndarray
    angle_index: np.ndarray
    angle_atoms: np.ndarray
    angle_species_index: np.ndarray
    triplet_index: np.ndarray
    dihedral_index: np.ndarray
    dihedral_atoms: np.ndarray
    dihedral_species_index: np.ndarray


@dataclass
class PreparedStructure:
    positions: np.ndarray
    cell: np.ndarray
    pbc: np.ndarray
    topology: MultigraphTopology
    target_2body: Optional[np.ndarray] = None
    target_3body: Optional[np.ndarray] = None
    target_4body: Optional[np.ndarray] = None


def _zero_init_linear(linear: nn.Linear) -> None:
    nn.init.zeros_(linear.weight)
    if linear.bias is not None:
        nn.init.zeros_(linear.bias)


def _float_tensor(array, device: torch.device) -> torch.Tensor:
    return torch.as_tensor(array, dtype=torch.float32, device=device)


def _long_tensor(array, device: torch.device) -> torch.Tensor:
    return torch.as_tensor(array, dtype=torch.long, device=device)


class ResidualMessageLayer(nn.Module):
    def __init__(self, hidden_dim: int, edge_dim: int, global_dim: int):
        super().__init__()
        self.message = nn.Sequential(
            nn.Linear(hidden_dim + edge_dim + global_dim, hidden_dim),
            nn.SiLU(),
            nn.Linear(hidden_dim, hidden_dim),
        )
        self.update = nn.Sequential(
            nn.Linear(2 * hidden_dim, hidden_dim),
            nn.SiLU(),
            nn.Linear(hidden_dim, hidden_dim),
        )
        _zero_init_linear(self.update[-1])

    def forward(
        self,
        hidden: torch.Tensor,
        edge_index: torch.Tensor,
        edge_attr: torch.Tensor,
        edge_weight: torch.Tensor,
        global_features: torch.Tensor,
    ) -> torch.Tensor:
        if edge_index.numel() == 0:
            return hidden

        src = edge_index[0]
        dst = edge_index[1]
        edge_global = global_features.expand(src.shape[0], -1)
        message_input = torch.cat([hidden[src], edge_attr, edge_global], dim=-1)
        message = self.message(message_input) * edge_weight.unsqueeze(-1)

        aggregated = torch.zeros_like(hidden)
        aggregated.index_add_(0, dst, message)

        normaliser = torch.zeros(hidden.shape[0], device=hidden.device, dtype=hidden.dtype)
        normaliser.index_add_(0, dst, edge_weight)
        aggregated = aggregated / normaliser.clamp_min(1.0e-6).unsqueeze(-1)

        return hidden + self.update(torch.cat([hidden, aggregated], dim=-1))


class GatedMessageLayer(nn.Module):
    def __init__(self, hidden_dim: int, edge_dim: int, global_dim: int):
        super().__init__()
        self.message = nn.Sequential(
            nn.Linear(hidden_dim + edge_dim + global_dim, hidden_dim),
            nn.SiLU(),
            nn.Linear(hidden_dim, hidden_dim),
        )
        self.gate = nn.Sequential(
            nn.Linear(hidden_dim + edge_dim + global_dim, hidden_dim),
            nn.SiLU(),
            nn.Linear(hidden_dim, hidden_dim),
        )
        self.update = nn.GRUCell(hidden_dim, hidden_dim)

    def forward(
        self,
        hidden: torch.Tensor,
        edge_index: torch.Tensor,
        edge_attr: torch.Tensor,
        edge_weight: torch.Tensor,
        global_features: torch.Tensor,
    ) -> torch.Tensor:
        if edge_index.numel() == 0:
            return hidden

        src = edge_index[0]
        dst = edge_index[1]
        edge_global = global_features.expand(src.shape[0], -1)
        message_input = torch.cat([hidden[src], edge_attr, edge_global], dim=-1)
        candidate = self.message(message_input)
        gate = torch.sigmoid(self.gate(message_input))
        message = gate * candidate * edge_weight.unsqueeze(-1)

        aggregated = torch.zeros_like(hidden)
        aggregated.index_add_(0, dst, message)

        normaliser = torch.zeros(hidden.shape[0], device=hidden.device, dtype=hidden.dtype)
        normaliser.index_add_(0, dst, edge_weight)
        aggregated = aggregated / normaliser.clamp_min(1.0e-6).unsqueeze(-1)

        return self.update(aggregated, hidden)


class AttentionMessageLayer(nn.Module):
    def __init__(self, hidden_dim: int, edge_dim: int, global_dim: int):
        super().__init__()
        input_dim = hidden_dim + edge_dim + global_dim
        self.key = nn.Linear(input_dim, hidden_dim)
        self.value = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.SiLU(),
            nn.Linear(hidden_dim, hidden_dim),
        )
        self.query = nn.Linear(hidden_dim + global_dim, hidden_dim)
        self.update = nn.Sequential(
            nn.Linear(2 * hidden_dim, hidden_dim),
            nn.SiLU(),
            nn.Linear(hidden_dim, hidden_dim),
        )
        _zero_init_linear(self.update[-1])

    def forward(
        self,
        hidden: torch.Tensor,
        edge_index: torch.Tensor,
        edge_attr: torch.Tensor,
        edge_weight: torch.Tensor,
        global_features: torch.Tensor,
    ) -> torch.Tensor:
        if edge_index.numel() == 0:
            return hidden

        src = edge_index[0]
        dst = edge_index[1]
        edge_global = global_features.expand(src.shape[0], -1)
        src_input = torch.cat([hidden[src], edge_attr, edge_global], dim=-1)
        dst_input = torch.cat([hidden[dst], edge_global], dim=-1)

        keys = self.key(src_input)
        values = self.value(src_input)
        queries = self.query(dst_input)
        scores = torch.sum(queries * keys, dim=-1) / math.sqrt(hidden.shape[-1])

        max_scores = torch.full(
            (hidden.shape[0],),
            torch.finfo(hidden.dtype).min,
            dtype=hidden.dtype,
            device=hidden.device,
        )
        max_scores.scatter_reduce_(0, dst, scores, reduce="amax", include_self=True)
        scaled_scores = torch.exp(scores - max_scores[dst]) * edge_weight.clamp_min(1.0e-6)
        score_sums = torch.zeros(hidden.shape[0], device=hidden.device, dtype=hidden.dtype)
        score_sums.index_add_(0, dst, scaled_scores)
        attention = scaled_scores / score_sums[dst].clamp_min(1.0e-6)

        aggregated = torch.zeros_like(hidden)
        aggregated.index_add_(0, dst, values * attention.unsqueeze(-1))
        return hidden + self.update(torch.cat([hidden, aggregated], dim=-1))


def _make_message_layer(layer_kind: str, hidden_dim: int, edge_dim: int, global_dim: int) -> nn.Module:
    if layer_kind == "residual":
        return ResidualMessageLayer(hidden_dim, edge_dim, global_dim)
    if layer_kind == "gated":
        return GatedMessageLayer(hidden_dim, edge_dim, global_dim)
    if layer_kind in {"attention", "attention_conservative"}:
        return AttentionMessageLayer(hidden_dim, edge_dim, global_dim)
    raise ValueError(f"Unsupported architecture '{layer_kind}'")


class GraphBranch(nn.Module):
    def __init__(
        self,
        node_dim: int,
        edge_dim: int,
        output_dim: int,
        hidden_dim: int,
        num_message_layers: int,
        global_dim: int,
        message_layer_kind: str,
    ):
        super().__init__()
        self.encoder = nn.Sequential(
            nn.Linear(node_dim + global_dim, hidden_dim),
            nn.SiLU(),
            nn.Linear(hidden_dim, hidden_dim),
            nn.SiLU(),
        )
        self.layers = nn.ModuleList(
            _make_message_layer(message_layer_kind, hidden_dim, edge_dim, global_dim)
            for _ in range(num_message_layers)
        )
        self.head = nn.Sequential(
            nn.Linear(hidden_dim + global_dim, hidden_dim),
            nn.SiLU(),
            nn.Linear(hidden_dim, output_dim),
        )
        _zero_init_linear(self.head[-1])
        self.residual_limit = 0.25 if message_layer_kind == "attention_conservative" else None
        self.base_scale = nn.Parameter(torch.ones(output_dim, dtype=torch.float32))
        self.base_bias = nn.Parameter(torch.zeros(output_dim, dtype=torch.float32))

    def forward(
        self,
        node_features: torch.Tensor,
        edge_index: torch.Tensor,
        edge_attr: torch.Tensor,
        edge_weight: torch.Tensor,
        global_features: torch.Tensor,
        base_vertex_fingerprint: torch.Tensor,
    ) -> Tuple[torch.Tensor, torch.Tensor]:
        global_on_nodes = global_features.expand(node_features.shape[0], -1)
        hidden = self.encoder(torch.cat([node_features, global_on_nodes], dim=-1))
        for layer in self.layers:
            hidden = layer(hidden, edge_index, edge_attr, edge_weight, global_features)

        residual = self.head(torch.cat([hidden, global_on_nodes], dim=-1))
        if self.residual_limit is not None:
            residual = self.residual_limit * torch.tanh(residual)
        vertex_fingerprint = (
            base_vertex_fingerprint * self.base_scale.unsqueeze(0)
            + self.base_bias.unsqueeze(0)
            + residual
        )
        graph_fingerprint = vertex_fingerprint.mean(dim=0)
        return vertex_fingerprint, graph_fingerprint


class TorchGNNFingerprint(nn.Module):
    """PyTorch multigraph surrogate for RAFFLE analytical fingerprints.

    The class preserves the fixed 2-/3-/4-body graph topology from the current
    Fortran implementation, uses the existing compiled RAFFLE interface as the
    analytical target oracle, and keeps the learnable part as a residual over a
    deterministic analytical-style multigraph base operator.
    """

    def __init__(
        self,
        species_list: Sequence[str],
        bond_cutoff: float = 6.0,
        hidden_dim: int = 128,
        num_message_layers: int = 2,
        learning_rate: float = 1.0e-3,
        lr_decay_rate: float = 1.0e-2,
        component_weight: Sequence[float] = (4.0, 1.0, 1.0),
        smooth_cutoff_width: float = 0.15,
        seed: int = 42,
        reference_hidden_sizes: Sequence[int] = (32,),
        reference_num_time_steps: int = 2,
        reference_output_dim: int = 16,
        reference_max_degree: int = 8,
        reference_layer_type: int = 1,
        reference_n_rbf: int = 12,
        reference_kernel_hidden: int = 32,
        architecture: str = "residual",
        device: Optional[str] = None,
    ):
        super().__init__()
        self.species_list = [str(symbol).strip() for symbol in species_list]
        self.num_species = len(self.species_list)
        self.bond_cutoff = float(bond_cutoff)
        self.learning_rate = float(learning_rate)
        self.lr_decay_rate = float(lr_decay_rate)
        self.smooth_cutoff_width = float(smooth_cutoff_width)
        self.seed = int(seed)
        self._rng = random.Random(seed)
        self._device = torch.device(device or ("cuda" if torch.cuda.is_available() else "cpu"))
        self.architecture = ARCHITECTURE_ALIASES.get(str(architecture).strip(), str(architecture).strip())
        self._use_component_coupling = self.architecture == "attention_coupled"
        self._message_layer_kind = "attention" if self._use_component_coupling else self.architecture

        torch.manual_seed(self.seed)
        np.random.seed(self.seed)

        self._species_to_index = {
            symbol: index for index, symbol in enumerate(self.species_list)
        }
        self._pair_to_index = {}
        pair_index = 0
        for left in range(self.num_species):
            for right in range(left, self.num_species):
                self._pair_to_index[(left, right)] = pair_index
                pair_index += 1
        self.num_pairs = pair_index

        self.reference_model = GNNFingerprint(
            species_list=self.species_list,
            bond_cutoff=self.bond_cutoff,
            gnn_hidden_sizes=list(reference_hidden_sizes),
            learning_rate=self.learning_rate,
            lr_decay_rate=self.lr_decay_rate,
            num_time_steps=int(reference_num_time_steps),
            gnn_output_dim=int(reference_output_dim),
            max_degree=int(reference_max_degree),
            layer_type=int(reference_layer_type),
            n_rbf=int(reference_n_rbf),
            kernel_hidden=int(reference_kernel_hidden),
            seed=self.seed,
        )

        self.fingerprint_dim_2body = int(self.reference_model.fingerprint_dim_2body)
        self.fingerprint_dim_3body = int(self.reference_model.fingerprint_dim_3body)
        self.fingerprint_dim_4body = int(self.reference_model.fingerprint_dim_4body)
        self.fingerprint_dim = int(self.reference_model.fingerprint_dim)

        self.nbins = (
            self.fingerprint_dim_2body // max(self.num_pairs, 1),
            self.fingerprint_dim_3body // max(self.num_species, 1),
            self.fingerprint_dim_4body // max(self.num_species, 1),
        )
        self.width = (0.025, math.pi / 64.0, math.pi / 64.0)
        self.sigma = (0.1, 0.1, 0.1)
        self.cutoff_min = (0.5, 0.0, 0.0)
        self.cutoff_max = (self.bond_cutoff, math.pi, math.pi)
        self.global_dim = 6
        self.component_weight = tuple(float(weight) for weight in component_weight)
        self.register_buffer(
            "_component_weight_tensor",
            torch.tensor(self.component_weight, dtype=torch.float32),
        )
        self.register_buffer(
            "_centers_2body",
            self.cutoff_min[0]
            + self.width[0] * torch.arange(self.nbins[0], dtype=torch.float32),
        )
        self.register_buffer(
            "_centers_3body",
            self.cutoff_min[1]
            + self.width[1] * torch.arange(self.nbins[1], dtype=torch.float32),
        )
        self.register_buffer(
            "_centers_4body",
            self.cutoff_min[2]
            + self.width[2] * torch.arange(self.nbins[2], dtype=torch.float32),
        )

        self.branch_2body = GraphBranch(
            node_dim=3 + self.num_species + 2,
            edge_dim=1,
            output_dim=self.fingerprint_dim_2body,
            hidden_dim=int(hidden_dim),
            num_message_layers=int(num_message_layers),
            global_dim=self.global_dim,
            message_layer_kind=self._message_layer_kind,
        )
        self.branch_3body = GraphBranch(
            node_dim=7 + 2 * self.num_species,
            edge_dim=1,
            output_dim=self.fingerprint_dim_3body,
            hidden_dim=int(hidden_dim),
            num_message_layers=int(num_message_layers),
            global_dim=self.global_dim,
            message_layer_kind=self._message_layer_kind,
        )
        self.branch_4body = GraphBranch(
            node_dim=6 + 3 * self.num_species,
            edge_dim=1,
            output_dim=self.fingerprint_dim_4body,
            hidden_dim=int(hidden_dim),
            num_message_layers=int(num_message_layers),
            global_dim=self.global_dim,
            message_layer_kind=self._message_layer_kind,
        )
        if self._use_component_coupling:
            context_input_dim = self.fingerprint_dim + self.global_dim
            self.component_context = nn.Sequential(
                nn.Linear(context_input_dim, int(hidden_dim)),
                nn.SiLU(),
                nn.Linear(int(hidden_dim), int(hidden_dim)),
                nn.SiLU(),
            )
            self.component_coupling_strength = 0.25
            self.component_gate_2body = nn.Linear(int(hidden_dim), self.fingerprint_dim_2body)
            self.component_gate_3body = nn.Linear(int(hidden_dim), self.fingerprint_dim_3body)
            self.component_gate_4body = nn.Linear(int(hidden_dim), self.fingerprint_dim_4body)
            self.component_bias_2body = nn.Linear(int(hidden_dim), self.fingerprint_dim_2body)
            self.component_bias_3body = nn.Linear(int(hidden_dim), self.fingerprint_dim_3body)
            self.component_bias_4body = nn.Linear(int(hidden_dim), self.fingerprint_dim_4body)
            for layer in (
                self.component_gate_2body,
                self.component_gate_3body,
                self.component_gate_4body,
                self.component_bias_2body,
                self.component_bias_3body,
                self.component_bias_4body,
            ):
                _zero_init_linear(layer)
        self._topology_cache = {}
        self._is_fitted = False
        self._optimiser = None
        self._scheduler = None
        self.to(self._device)

    @property
    def is_fitted(self) -> bool:
        return self._is_fitted

    def _project_fingerprint_tensor(self, fingerprint: torch.Tensor) -> torch.Tensor:
        # Preserve calibrated positive outputs exactly while smoothly folding any
        # negative tail back into the physical non-negative fingerprint domain.
        return torch.where(
            fingerprint >= 0.0,
            fingerprint,
            torch_functional.softplus(
                fingerprint,
                beta=FINGERPRINT_NEGATIVE_TAIL_BETA,
                threshold=20.0,
            ),
        )

    def _project_fingerprint_targets(self, fingerprint: torch.Tensor) -> torch.Tensor:
        return fingerprint.clamp_min(0.0)

    def _project_vertex_fingerprints(
        self,
        vertices: Tuple[torch.Tensor, torch.Tensor, torch.Tensor],
    ) -> Tuple[
        Tuple[torch.Tensor, torch.Tensor, torch.Tensor],
        Tuple[torch.Tensor, torch.Tensor, torch.Tensor],
    ]:
        projected_vertices = tuple(
            self._project_fingerprint_tensor(vertex) for vertex in vertices
        )
        projected_fingerprints = tuple(
            vertex.mean(dim=0) for vertex in projected_vertices
        )
        return (
            (projected_vertices[0], projected_vertices[1], projected_vertices[2]),
            (projected_fingerprints[0], projected_fingerprints[1], projected_fingerprints[2]),
        )

    def _element_properties(self, symbol: str) -> Tuple[float, float]:
        return ELEMENT_PROPERTIES.get(symbol, (0.0, 1.0))

    def _topology_key(self, symbols: Sequence[str]) -> Tuple[str, ...]:
        return tuple(str(symbol).strip() for symbol in symbols)

    def _build_topology(self, symbols: Sequence[str]) -> MultigraphTopology:
        key = self._topology_key(symbols)
        if key in self._topology_cache:
            return self._topology_cache[key]

        species_index = np.asarray(
            [self._species_to_index[str(symbol).strip()] for symbol in key],
            dtype=np.int64,
        )
        atomic_numbers = np.asarray(
            [self._element_properties(symbol)[0] for symbol in key],
            dtype=np.float32,
        )
        covalent_radii = np.asarray(
            [self._element_properties(symbol)[1] for symbol in key],
            dtype=np.float32,
        )

        num_atoms = len(key)
        pair_records = []
        pair_type_records = []
        for atom_i in range(num_atoms - 1):
            for atom_j in range(atom_i + 1, num_atoms):
                pair_records.append((atom_i, atom_j))
                pair_species = tuple(sorted((species_index[atom_i], species_index[atom_j])))
                pair_type_records.append(self._pair_to_index[pair_species])
        pair_index = np.asarray(pair_records, dtype=np.int64)
        pair_type_index = np.asarray(pair_type_records, dtype=np.int64)

        angle_index = []
        angle_atoms = []
        angle_species_index = []
        if len(pair_records) >= 2:
            for pair_i in range(len(pair_records) - 1):
                atom_a = pair_records[pair_i]
                for pair_j in range(pair_i + 1, len(pair_records)):
                    atom_b = pair_records[pair_j]
                    shared = sorted(set(atom_a).intersection(atom_b))
                    if len(shared) != 1:
                        continue
                    shared_atom = shared[0]
                    other_i = atom_a[1] if atom_a[0] == shared_atom else atom_a[0]
                    other_j = atom_b[1] if atom_b[0] == shared_atom else atom_b[0]
                    angle_index.append((pair_i, pair_j))
                    angle_atoms.append((shared_atom, other_i, other_j))
                    angle_species_index.append(species_index[shared_atom])

        triplet_records = []
        if len(pair_records) >= 2:
            for atom_j in range(num_atoms):
                incident = []
                for pair_id, (atom_i, atom_k) in enumerate(pair_records):
                    if atom_i == atom_j:
                        incident.append((pair_id, atom_k))
                    elif atom_k == atom_j:
                        incident.append((pair_id, atom_i))
                for _, atom_i in incident:
                    for _, atom_k in incident:
                        if atom_i == atom_k:
                            continue
                        triplet_records.append((atom_i, atom_j, atom_k))

        dihedral_index = []
        dihedral_atoms = []
        dihedral_species_index = []
        if len(triplet_records) >= 2:
            for triplet_i, (atom_i, atom_j, atom_k) in enumerate(triplet_records):
                for triplet_j, triplet_next in enumerate(triplet_records):
                    if triplet_i == triplet_j:
                        continue
                    if atom_j != triplet_next[0] or atom_k != triplet_next[1]:
                        continue
                    atom_l = triplet_next[2]
                    if atom_i == atom_l:
                        continue
                    dihedral_index.append((triplet_i, triplet_j))
                    dihedral_atoms.append((atom_i, atom_j, atom_k, atom_l))
                    dihedral_species_index.append(species_index[atom_j])

        topology = MultigraphTopology(
            symbols=key,
            species_index=species_index,
            atomic_numbers=atomic_numbers,
            covalent_radii=covalent_radii,
            pair_index=np.asarray(pair_index, dtype=np.int64).reshape(-1, 2),
            pair_type_index=np.asarray(pair_type_index, dtype=np.int64),
            angle_index=np.asarray(angle_index, dtype=np.int64).reshape(-1, 2),
            angle_atoms=np.asarray(angle_atoms, dtype=np.int64).reshape(-1, 3),
            angle_species_index=np.asarray(angle_species_index, dtype=np.int64),
            triplet_index=np.asarray(triplet_records, dtype=np.int64).reshape(-1, 3),
            dihedral_index=np.asarray(dihedral_index, dtype=np.int64).reshape(-1, 2),
            dihedral_atoms=np.asarray(dihedral_atoms, dtype=np.int64).reshape(-1, 4),
            dihedral_species_index=np.asarray(dihedral_species_index, dtype=np.int64),
        )
        self._topology_cache[key] = topology
        return topology

    def prepare_structure(self, atoms, include_targets: bool = True) -> PreparedStructure:
        symbols = tuple(str(symbol).strip() for symbol in atoms.get_chemical_symbols())
        topology = self._build_topology(symbols)
        targets = (None, None, None)
        if include_targets:
            targets = self.reference_model.compute_fingerprint_components(atoms)

        return PreparedStructure(
            positions=np.asarray(atoms.get_positions(), dtype=np.float32),
            cell=np.asarray(atoms.cell.array, dtype=np.float32),
            pbc=np.asarray(atoms.pbc, dtype=bool),
            topology=topology,
            target_2body=None if targets[0] is None else np.asarray(targets[0], dtype=np.float32),
            target_3body=None if targets[1] is None else np.asarray(targets[1], dtype=np.float32),
            target_4body=None if targets[2] is None else np.asarray(targets[2], dtype=np.float32),
        )

    def prepare_dataset(self, structures: Iterable, include_targets: bool = True) -> list[PreparedStructure]:
        return [self.prepare_structure(atoms, include_targets=include_targets) for atoms in structures]

    def compute_reference_components(self, atoms) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        return self.reference_model.compute_fingerprint_components(atoms)

    def compute_reference_fingerprint(self, atoms) -> np.ndarray:
        return self.reference_model.compute_fingerprint(atoms)

    def _minimum_image_delta(
        self,
        cell: torch.Tensor,
        pbc: torch.Tensor,
        delta_cart: torch.Tensor,
    ) -> torch.Tensor:
        if delta_cart.numel() == 0 or not bool(pbc.any()):
            return delta_cart
        inverse_cell = torch.linalg.inv(cell)
        delta_frac = delta_cart @ inverse_cell
        if pbc[0]:
            delta_frac[:, 0] = delta_frac[:, 0] - torch.round(delta_frac[:, 0])
        if pbc[1]:
            delta_frac[:, 1] = delta_frac[:, 1] - torch.round(delta_frac[:, 1])
        if pbc[2]:
            delta_frac[:, 2] = delta_frac[:, 2] - torch.round(delta_frac[:, 2])
        return delta_frac @ cell

    def _smooth_cutoff(self, distance: torch.Tensor) -> torch.Tensor:
        return torch.sigmoid((self.bond_cutoff - distance) / self.smooth_cutoff_width)

    def _gaussian_basis(
        self,
        value: torch.Tensor,
        centers: torch.Tensor,
        sigma: float,
    ) -> torch.Tensor:
        basis = torch.exp(-0.5 * ((value.unsqueeze(-1) - centers.unsqueeze(0)) / sigma) ** 2)
        normaliser = basis.sum(dim=-1, keepdim=True).clamp_min(1.0e-8)
        return basis / normaliser

    def _angle(self, vector_a: torch.Tensor, vector_b: torch.Tensor) -> torch.Tensor:
        numerator = torch.sum(vector_a * vector_b, dim=-1)
        denominator = vector_a.norm(dim=-1) * vector_b.norm(dim=-1)
        cosine = numerator / denominator.clamp_min(1.0e-8)
        return torch.acos(cosine.clamp(-1.0 + 1.0e-7, 1.0 - 1.0e-7))

    def _improper_dihedral(
        self,
        vector_ij: torch.Tensor,
        vector_jk: torch.Tensor,
        vector_kl: torch.Tensor,
    ) -> torch.Tensor:
        normal_1 = torch.cross(vector_ij, vector_jk, dim=-1)
        normal_2 = torch.cross(vector_jk, vector_kl, dim=-1)
        numerator = torch.sum(normal_1 * normal_2, dim=-1)
        denominator = normal_1.norm(dim=-1) * normal_2.norm(dim=-1)
        cosine = numerator / denominator.clamp_min(1.0e-8)
        return torch.acos(cosine.clamp(-1.0 + 1.0e-7, 1.0 - 1.0e-7))

    def _pair_block_vectors(
        self,
        basis: torch.Tensor,
        pair_type_index: torch.Tensor,
    ) -> torch.Tensor:
        one_hot = torch_functional.one_hot(pair_type_index, num_classes=self.num_pairs).to(basis.dtype)
        return (one_hot.unsqueeze(-1) * basis.unsqueeze(1)).reshape(basis.shape[0], -1)

    def _species_block_vectors(
        self,
        basis: torch.Tensor,
        species_index: torch.Tensor,
        output_size: int,
    ) -> torch.Tensor:
        one_hot = torch_functional.one_hot(species_index, num_classes=self.num_species).to(basis.dtype)
        return (one_hot.unsqueeze(-1) * basis.unsqueeze(1)).reshape(basis.shape[0], output_size)

    def _lattice_features(self, cell: torch.Tensor) -> torch.Tensor:
        a = torch.norm(cell[0])
        b = torch.norm(cell[1])
        c = torch.norm(cell[2])

        alpha = self._angle(cell[1].unsqueeze(0), cell[2].unsqueeze(0))[0]
        beta = self._angle(cell[0].unsqueeze(0), cell[2].unsqueeze(0))[0]
        gamma = self._angle(cell[0].unsqueeze(0), cell[1].unsqueeze(0))[0]
        return torch.stack(
            [
                a / self.bond_cutoff,
                b / self.bond_cutoff,
                c / self.bond_cutoff,
                alpha / math.pi,
                beta / math.pi,
                gamma / math.pi,
            ]
        ).unsqueeze(0)

    def _build_multigraph_tensors(
        self,
        prepared: PreparedStructure,
        positions_override: Optional[torch.Tensor] = None,
    ):
        device = self._device
        topology = prepared.topology
        positions = positions_override
        if positions is None:
            positions = _float_tensor(prepared.positions, device)
        cell = _float_tensor(prepared.cell, device)
        pbc = torch.as_tensor(prepared.pbc, dtype=torch.bool, device=device)
        species_index = _long_tensor(topology.species_index, device)
        atomic_numbers = _float_tensor(topology.atomic_numbers, device)
        covalent_radii = _float_tensor(topology.covalent_radii, device)
        species_one_hot = torch_functional.one_hot(species_index, num_classes=self.num_species).to(torch.float32)
        global_features = self._lattice_features(cell)

        atom_node_features = torch.cat(
            [
                positions / self.bond_cutoff,
                species_one_hot,
                atomic_numbers.unsqueeze(-1),
                covalent_radii.unsqueeze(-1),
            ],
            dim=-1,
        )

        num_atoms = positions.shape[0]
        pair_index = _long_tensor(topology.pair_index, device)
        if pair_index.numel() > 0:
            pair_left = pair_index[:, 0]
            pair_right = pair_index[:, 1]
            pair_delta = self._minimum_image_delta(cell, pbc, positions[pair_right] - positions[pair_left])
            pair_distance = pair_delta.norm(dim=-1)
            pair_weight = self._smooth_cutoff(pair_distance)
            pair_midpoint = 0.5 * (positions[pair_left] + positions[pair_right]) / self.bond_cutoff
            pair_unit = pair_delta / pair_distance.clamp_min(1.0e-8).unsqueeze(-1)
            pair_node_features = torch.cat(
                [
                    pair_midpoint,
                    pair_unit,
                    (pair_distance / self.bond_cutoff).unsqueeze(-1),
                    species_one_hot[pair_left],
                    species_one_hot[pair_right],
                ],
                dim=-1,
            )
            pair_weight_matrix = torch.zeros((num_atoms, num_atoms), dtype=torch.float32, device=device)
            pair_weight_matrix[pair_left, pair_right] = pair_weight
            pair_weight_matrix[pair_right, pair_left] = pair_weight

            atom_edge_index = torch.cat(
                [pair_index, pair_index[:, [1, 0]]],
                dim=0,
            ).T
            atom_edge_attr = torch.cat(
                [
                    (pair_distance / self.bond_cutoff).unsqueeze(-1),
                    (pair_distance / self.bond_cutoff).unsqueeze(-1),
                ],
                dim=0,
            )
            atom_edge_weight = torch.cat([pair_weight, pair_weight], dim=0)

            pair_basis = self._gaussian_basis(pair_distance, self._centers_2body, self.sigma[0])
            pair_vector = self._pair_block_vectors(
                pair_basis,
                _long_tensor(topology.pair_type_index, device),
            )
            pair_vector = pair_vector * pair_weight.unsqueeze(-1)
            pair_normaliser = pair_weight.sum().clamp_min(1.0e-8)
            pair_vector = pair_vector / pair_normaliser
            atom_base = torch.zeros(
                (num_atoms, self.fingerprint_dim_2body),
                dtype=torch.float32,
                device=device,
            )
            atom_scale = max(num_atoms, 1) / 2.0
            atom_base.index_add_(0, pair_left, pair_vector * atom_scale)
            atom_base.index_add_(0, pair_right, pair_vector * atom_scale)
        else:
            pair_node_features = torch.zeros(
                (1, 7 + 2 * self.num_species),
                dtype=torch.float32,
                device=device,
            )
            pair_weight_matrix = torch.zeros((num_atoms, num_atoms), dtype=torch.float32, device=device)
            atom_edge_index = torch.zeros((2, 0), dtype=torch.long, device=device)
            atom_edge_attr = torch.zeros((0, 1), dtype=torch.float32, device=device)
            atom_edge_weight = torch.zeros((0,), dtype=torch.float32, device=device)
            atom_base = torch.zeros(
                (num_atoms, self.fingerprint_dim_2body),
                dtype=torch.float32,
                device=device,
            )

        angle_index = _long_tensor(topology.angle_index, device)
        angle_atoms = _long_tensor(topology.angle_atoms, device)
        if angle_index.numel() > 0:
            shared = angle_atoms[:, 0]
            other_i = angle_atoms[:, 1]
            other_j = angle_atoms[:, 2]
            angle_vec_i = self._minimum_image_delta(cell, pbc, positions[other_i] - positions[shared])
            angle_vec_j = self._minimum_image_delta(cell, pbc, positions[other_j] - positions[shared])
            angle_value = self._angle(angle_vec_i, angle_vec_j)
            angle_weight = (
                pair_weight_matrix[shared, other_i]
                * pair_weight_matrix[shared, other_j]
            )
            pair_edge_index = torch.cat(
                [angle_index, angle_index[:, [1, 0]]],
                dim=0,
            ).T
            pair_edge_attr = torch.cat(
                [(angle_value / math.pi).unsqueeze(-1), (angle_value / math.pi).unsqueeze(-1)],
                dim=0,
            )
            pair_edge_weight = torch.cat([angle_weight, angle_weight], dim=0)

            angle_basis = self._gaussian_basis(angle_value, self._centers_3body, self.sigma[1])
            angle_vector = self._species_block_vectors(
                angle_basis,
                _long_tensor(topology.angle_species_index, device),
                self.fingerprint_dim_3body,
            )
            angle_vector = angle_vector * angle_weight.unsqueeze(-1)
            angle_normaliser = angle_weight.sum().clamp_min(1.0e-8)
            angle_vector = angle_vector / angle_normaliser
            pair_base = torch.zeros(
                (pair_node_features.shape[0], self.fingerprint_dim_3body),
                dtype=torch.float32,
                device=device,
            )
            pair_scale = max(pair_node_features.shape[0], 1) / 2.0
            pair_base.index_add_(0, angle_index[:, 0], angle_vector * pair_scale)
            pair_base.index_add_(0, angle_index[:, 1], angle_vector * pair_scale)
        else:
            pair_edge_index = torch.zeros((2, 0), dtype=torch.long, device=device)
            pair_edge_attr = torch.zeros((0, 1), dtype=torch.float32, device=device)
            pair_edge_weight = torch.zeros((0,), dtype=torch.float32, device=device)
            pair_base = torch.zeros(
                (pair_node_features.shape[0], self.fingerprint_dim_3body),
                dtype=torch.float32,
                device=device,
            )

        triplet_index = _long_tensor(topology.triplet_index, device)
        if triplet_index.numel() > 0:
            triplet_i = triplet_index[:, 0]
            triplet_j = triplet_index[:, 1]
            triplet_k = triplet_index[:, 2]
            vector_ij = self._minimum_image_delta(cell, pbc, positions[triplet_i] - positions[triplet_j])
            vector_jk = self._minimum_image_delta(cell, pbc, positions[triplet_k] - positions[triplet_j])
            distance_ij = vector_ij.norm(dim=-1)
            distance_jk = vector_jk.norm(dim=-1)
            triplet_angle = self._angle(vector_ij, vector_jk)
            triplet_centroid = (
                positions[triplet_i] + positions[triplet_j] + positions[triplet_k]
            ) / (3.0 * self.bond_cutoff)
            triplet_node_features = torch.cat(
                [
                    triplet_centroid,
                    (distance_ij / self.bond_cutoff).unsqueeze(-1),
                    (distance_jk / self.bond_cutoff).unsqueeze(-1),
                    (triplet_angle / math.pi).unsqueeze(-1),
                    species_one_hot[triplet_i],
                    species_one_hot[triplet_j],
                    species_one_hot[triplet_k],
                ],
                dim=-1,
            )
        else:
            triplet_node_features = torch.zeros(
                (1, 6 + 3 * self.num_species),
                dtype=torch.float32,
                device=device,
            )

        dihedral_index = _long_tensor(topology.dihedral_index, device)
        dihedral_atoms = _long_tensor(topology.dihedral_atoms, device)
        if dihedral_index.numel() > 0:
            atom_i = dihedral_atoms[:, 0]
            atom_j = dihedral_atoms[:, 1]
            atom_k = dihedral_atoms[:, 2]
            atom_l = dihedral_atoms[:, 3]
            vector_ij = self._minimum_image_delta(cell, pbc, positions[atom_i] - positions[atom_j])
            vector_jk = self._minimum_image_delta(cell, pbc, positions[atom_k] - positions[atom_j])
            vector_kl = self._minimum_image_delta(cell, pbc, positions[atom_l] - positions[atom_k])
            dihedral_value = self._improper_dihedral(vector_ij, vector_jk, vector_kl)
            dihedral_weight = (
                pair_weight_matrix[atom_i, atom_j]
                * pair_weight_matrix[atom_j, atom_k]
                * pair_weight_matrix[atom_k, atom_l]
            )
            triplet_edge_index = torch.cat(
                [dihedral_index, dihedral_index[:, [1, 0]]],
                dim=0,
            ).T
            triplet_edge_attr = torch.cat(
                [(dihedral_value / math.pi).unsqueeze(-1), (dihedral_value / math.pi).unsqueeze(-1)],
                dim=0,
            )
            triplet_edge_weight = torch.cat([dihedral_weight, dihedral_weight], dim=0)

            dihedral_basis = self._gaussian_basis(dihedral_value, self._centers_4body, self.sigma[2])
            dihedral_vector = self._species_block_vectors(
                dihedral_basis,
                _long_tensor(topology.dihedral_species_index, device),
                self.fingerprint_dim_4body,
            )
            dihedral_vector = dihedral_vector * dihedral_weight.unsqueeze(-1)
            dihedral_normaliser = dihedral_weight.sum().clamp_min(1.0e-8)
            dihedral_vector = dihedral_vector / dihedral_normaliser
            triplet_base = torch.zeros(
                (triplet_node_features.shape[0], self.fingerprint_dim_4body),
                dtype=torch.float32,
                device=device,
            )
            triplet_scale = max(triplet_node_features.shape[0], 1) / 2.0
            triplet_base.index_add_(0, dihedral_index[:, 0], dihedral_vector * triplet_scale)
            triplet_base.index_add_(0, dihedral_index[:, 1], dihedral_vector * triplet_scale)
        else:
            triplet_edge_index = torch.zeros((2, 0), dtype=torch.long, device=device)
            triplet_edge_attr = torch.zeros((0, 1), dtype=torch.float32, device=device)
            triplet_edge_weight = torch.zeros((0,), dtype=torch.float32, device=device)
            triplet_base = torch.zeros(
                (triplet_node_features.shape[0], self.fingerprint_dim_4body),
                dtype=torch.float32,
                device=device,
            )

        return {
            "global_features": global_features,
            "atom_node_features": atom_node_features,
            "atom_edge_index": atom_edge_index,
            "atom_edge_attr": atom_edge_attr,
            "atom_edge_weight": atom_edge_weight,
            "atom_base": atom_base,
            "pair_node_features": pair_node_features,
            "pair_edge_index": pair_edge_index,
            "pair_edge_attr": pair_edge_attr,
            "pair_edge_weight": pair_edge_weight,
            "pair_base": pair_base,
            "triplet_node_features": triplet_node_features,
            "triplet_edge_index": triplet_edge_index,
            "triplet_edge_attr": triplet_edge_attr,
            "triplet_edge_weight": triplet_edge_weight,
            "triplet_base": triplet_base,
        }

    def _forward_prepared(
        self,
        prepared: PreparedStructure,
        positions_override: Optional[torch.Tensor] = None,
        return_vertices: bool = False,
    ):
        graph = self._build_multigraph_tensors(prepared, positions_override=positions_override)
        vertex_2body, fingerprint_2body = self.branch_2body(
            graph["atom_node_features"],
            graph["atom_edge_index"],
            graph["atom_edge_attr"],
            graph["atom_edge_weight"],
            graph["global_features"],
            graph["atom_base"],
        )
        vertex_3body, fingerprint_3body = self.branch_3body(
            graph["pair_node_features"],
            graph["pair_edge_index"],
            graph["pair_edge_attr"],
            graph["pair_edge_weight"],
            graph["global_features"],
            graph["pair_base"],
        )
        vertex_4body, fingerprint_4body = self.branch_4body(
            graph["triplet_node_features"],
            graph["triplet_edge_index"],
            graph["triplet_edge_attr"],
            graph["triplet_edge_weight"],
            graph["global_features"],
            graph["triplet_base"],
        )
        if self._use_component_coupling:
            (
                (vertex_2body, vertex_3body, vertex_4body),
                (fingerprint_2body, fingerprint_3body, fingerprint_4body),
            ) = self._apply_component_coupling(
                graph["global_features"],
                (vertex_2body, vertex_3body, vertex_4body),
                (fingerprint_2body, fingerprint_3body, fingerprint_4body),
            )
        (
            (vertex_2body, vertex_3body, vertex_4body),
            (fingerprint_2body, fingerprint_3body, fingerprint_4body),
        ) = self._project_vertex_fingerprints(
            (vertex_2body, vertex_3body, vertex_4body)
        )
        if return_vertices:
            return (
                (vertex_2body, vertex_3body, vertex_4body),
                (fingerprint_2body, fingerprint_3body, fingerprint_4body),
            )
        return fingerprint_2body, fingerprint_3body, fingerprint_4body

    def _apply_component_coupling(
        self,
        global_features: torch.Tensor,
        vertices: Tuple[torch.Tensor, torch.Tensor, torch.Tensor],
        fingerprints: Tuple[torch.Tensor, torch.Tensor, torch.Tensor],
    ) -> Tuple[
        Tuple[torch.Tensor, torch.Tensor, torch.Tensor],
        Tuple[torch.Tensor, torch.Tensor, torch.Tensor],
    ]:
        context_input = torch.cat(
            [
                fingerprints[0],
                fingerprints[1],
                fingerprints[2],
                global_features.squeeze(0),
            ],
            dim=0,
        ).unsqueeze(0)
        context = self.component_context(context_input).squeeze(0)

        adjusted_vertices = []
        adjusted_fingerprints = []
        for vertex, gate_layer, bias_layer in (
            (vertices[0], self.component_gate_2body, self.component_bias_2body),
            (vertices[1], self.component_gate_3body, self.component_bias_3body),
            (vertices[2], self.component_gate_4body, self.component_bias_4body),
        ):
            gate = self.component_coupling_strength * torch.tanh(gate_layer(context))
            bias = self.component_coupling_strength * torch.tanh(bias_layer(context))
            adjusted_vertex = vertex * (1.0 + gate.unsqueeze(0)) + bias.unsqueeze(0)
            adjusted_vertices.append(adjusted_vertex)
            adjusted_fingerprints.append(adjusted_vertex.mean(dim=0))

        return (
            (adjusted_vertices[0], adjusted_vertices[1], adjusted_vertices[2]),
            (adjusted_fingerprints[0], adjusted_fingerprints[1], adjusted_fingerprints[2]),
        )

    def _component_loss(
        self,
        predicted_2body: torch.Tensor,
        predicted_3body: torch.Tensor,
        predicted_4body: torch.Tensor,
        target_2body: torch.Tensor,
        target_3body: torch.Tensor,
        target_4body: torch.Tensor,
    ) -> torch.Tensor:
        target_2body = self._project_fingerprint_targets(target_2body)
        target_3body = self._project_fingerprint_targets(target_3body)
        target_4body = self._project_fingerprint_targets(target_4body)
        loss_2body = torch.mean((predicted_2body - target_2body) ** 2)
        loss_3body = torch.mean((predicted_3body - target_3body) ** 2)
        loss_4body = torch.mean((predicted_4body - target_4body) ** 2)
        return (
            self._component_weight_tensor[0] * loss_2body
            + self._component_weight_tensor[1] * loss_3body
            + self._component_weight_tensor[2] * loss_4body
        )

    def _evaluate_entries(self, entries: Sequence[PreparedStructure]) -> float:
        self.eval()
        losses = []
        with torch.no_grad():
            for entry in entries:
                target_2body = _float_tensor(entry.target_2body, self._device)
                target_3body = _float_tensor(entry.target_3body, self._device)
                target_4body = _float_tensor(entry.target_4body, self._device)
                prediction = self._forward_prepared(entry)
                losses.append(
                    float(
                        self._component_loss(
                            prediction[0],
                            prediction[1],
                            prediction[2],
                            target_2body,
                            target_3body,
                            target_4body,
                        ).item()
                    )
                )
        return float(np.mean(losses)) if losses else 0.0

    def _calibrate_base_affine(self, entries: Sequence[PreparedStructure]) -> None:
        if not entries:
            return

        self.eval()
        with torch.no_grad():
            base_2body = []
            base_3body = []
            base_4body = []
            target_2body = []
            target_3body = []
            target_4body = []
            for entry in entries:
                graph = self._build_multigraph_tensors(entry)
                base_2body.append(graph["atom_base"].mean(dim=0))
                base_3body.append(graph["pair_base"].mean(dim=0))
                base_4body.append(graph["triplet_base"].mean(dim=0))
                target_2body.append(_float_tensor(entry.target_2body, self._device))
                target_3body.append(_float_tensor(entry.target_3body, self._device))
                target_4body.append(_float_tensor(entry.target_4body, self._device))

            mean_base_2body = torch.stack(base_2body).mean(dim=0)
            mean_base_3body = torch.stack(base_3body).mean(dim=0)
            mean_base_4body = torch.stack(base_4body).mean(dim=0)
            mean_target_2body = torch.stack(target_2body).mean(dim=0)
            mean_target_3body = torch.stack(target_3body).mean(dim=0)
            mean_target_4body = torch.stack(target_4body).mean(dim=0)

            for branch, mean_base, mean_target in (
                (self.branch_2body, mean_base_2body, mean_target_2body),
                (self.branch_3body, mean_base_3body, mean_target_3body),
                (self.branch_4body, mean_base_4body, mean_target_4body),
            ):
                scale = torch.where(
                    mean_base.abs() > 1.0e-6,
                    mean_target / mean_base,
                    torch.ones_like(mean_base),
                )
                scale = scale.clamp(0.1, 10.0)
                bias = mean_target - scale * mean_base
                branch.base_scale.copy_(scale)
                branch.base_bias.copy_(bias)

    def fit(
        self,
        structures: Sequence,
        num_epochs: int = 100,
        batch_size: int = 16,
        augment_structures: Optional[Sequence] = None,
        verbose: int = 0,
        reset_optimiser: bool = False,
        recalibrate_base: Optional[bool] = None,
    ) -> list[float]:
        combined_structures = list(structures)
        if augment_structures:
            combined_structures.extend(list(augment_structures))
        entries = self.prepare_dataset(combined_structures, include_targets=True)
        if recalibrate_base is None:
            recalibrate_base = not self._is_fitted
        if recalibrate_base:
            self._calibrate_base_affine(entries)

        initial_loss = self._evaluate_entries(entries)
        history = [initial_loss]

        if reset_optimiser or self._optimiser is None:
            self._optimiser = torch.optim.Adam(self.parameters(), lr=self.learning_rate)
            self._scheduler = torch.optim.lr_scheduler.ExponentialLR(
                self._optimiser,
                gamma=float(math.exp(-self.lr_decay_rate)),
            )
        optimiser = self._optimiser
        scheduler = self._scheduler

        for epoch in range(int(num_epochs)):
            self.train(mode=True)
            shuffled = list(entries)
            self._rng.shuffle(shuffled)
            for start in range(0, len(shuffled), max(int(batch_size), 1)):
                batch = shuffled[start:start + max(int(batch_size), 1)]
                optimiser.zero_grad()
                loss = torch.zeros((), dtype=torch.float32, device=self._device)
                for entry in batch:
                    target_2body = _float_tensor(entry.target_2body, self._device)
                    target_3body = _float_tensor(entry.target_3body, self._device)
                    target_4body = _float_tensor(entry.target_4body, self._device)
                    prediction = self._forward_prepared(entry)
                    loss = loss + self._component_loss(
                        prediction[0],
                        prediction[1],
                        prediction[2],
                        target_2body,
                        target_3body,
                        target_4body,
                    )
                loss = loss / max(len(batch), 1)
                loss.backward()
                torch.nn.utils.clip_grad_value_(self.parameters(), 1.0e-1)
                torch.nn.utils.clip_grad_norm_(self.parameters(), 1.0e-1)
                optimiser.step()
            scheduler.step()
            epoch_loss = self._evaluate_entries(entries)
            history.append(epoch_loss)
            if verbose > 0 and ((epoch + 1) % max(int(num_epochs) // 10, 1) == 0 or epoch == 0):
                print(f"epoch={epoch + 1:4d} loss={epoch_loss:.6e}")

        self._is_fitted = True
        return history

    def predict_components(self, atoms) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        self.eval()
        prepared = self.prepare_structure(atoms, include_targets=False)
        with torch.no_grad():
            fingerprint_2body, fingerprint_3body, fingerprint_4body = self._forward_prepared(prepared)
        return (
            fingerprint_2body.detach().cpu().numpy().astype(np.float32),
            fingerprint_3body.detach().cpu().numpy().astype(np.float32),
            fingerprint_4body.detach().cpu().numpy().astype(np.float32),
        )

    def predict(self, atoms) -> np.ndarray:
        fingerprint_2body, fingerprint_3body, fingerprint_4body = self.predict_components(atoms)
        return np.concatenate([fingerprint_2body, fingerprint_3body, fingerprint_4body]).astype(np.float32)

    def compute_fingerprint(self, atoms) -> np.ndarray:
        return self.predict(atoms)

    def compute_vertex_fingerprints(self, atoms) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        self.eval()
        prepared = self.prepare_structure(atoms, include_targets=False)
        with torch.no_grad():
            (vertex_2body, vertex_3body, vertex_4body), _ = self._forward_prepared(
                prepared,
                return_vertices=True,
            )
        return (
            vertex_2body.detach().cpu().numpy().astype(np.float32),
            vertex_3body.detach().cpu().numpy().astype(np.float32),
            vertex_4body.detach().cpu().numpy().astype(np.float32),
        )

    def compute_gradients(self, atoms) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        self.eval()
        prepared = self.prepare_structure(atoms, include_targets=False)
        positions = _float_tensor(prepared.positions, self._device)
        positions.requires_grad_(True)

        def branch_2body(pos: torch.Tensor) -> torch.Tensor:
            return self._forward_prepared(prepared, positions_override=pos)[0]

        def branch_3body(pos: torch.Tensor) -> torch.Tensor:
            return self._forward_prepared(prepared, positions_override=pos)[1]

        def branch_4body(pos: torch.Tensor) -> torch.Tensor:
            return self._forward_prepared(prepared, positions_override=pos)[2]

        jacobian_2body = torch.autograd.functional.jacobian(branch_2body, positions)
        jacobian_3body = torch.autograd.functional.jacobian(branch_3body, positions)
        jacobian_4body = torch.autograd.functional.jacobian(branch_4body, positions)
        return (
            jacobian_2body.detach().cpu().numpy().transpose(1, 2, 0).astype(np.float32),
            jacobian_3body.detach().cpu().numpy().transpose(1, 2, 0).astype(np.float32),
            jacobian_4body.detach().cpu().numpy().transpose(1, 2, 0).astype(np.float32),
        )

    def _positions_to_loss(
        self,
        prepared: PreparedStructure,
        positions: torch.Tensor,
        target_2body: torch.Tensor,
        target_3body: torch.Tensor,
        target_4body: torch.Tensor,
        reference_positions: Optional[torch.Tensor],
        target_vertex_fingerprints: Optional[Tuple[torch.Tensor, torch.Tensor, torch.Tensor]] = None,
        target_positions: Optional[torch.Tensor] = None,
        fixed_mask: Optional[torch.Tensor] = None,
        fingerprint_loss_weight: float = 1.0,
        target_vertex_weight: float = 0.0,
        target_position_weight: float = 0.0,
    ) -> Tuple[torch.Tensor, torch.Tensor]:
        (
            (vertex_2body, vertex_3body, vertex_4body),
            (prediction_2body, prediction_3body, prediction_4body),
        ) = self._forward_prepared(
            prepared,
            positions_override=positions,
            return_vertices=True,
        )
        fingerprint_loss = self._component_loss(
            prediction_2body,
            prediction_3body,
            prediction_4body,
            target_2body,
            target_3body,
            target_4body,
        )
        total_loss = float(fingerprint_loss_weight) * fingerprint_loss

        topology = prepared.topology
        pair_index = _long_tensor(topology.pair_index, self._device)
        cell = _float_tensor(prepared.cell, self._device)
        pbc = torch.as_tensor(prepared.pbc, dtype=torch.bool, device=self._device)
        covalent_radii = _float_tensor(topology.covalent_radii, self._device)

        regularisation = torch.zeros((), dtype=torch.float32, device=self._device)
        if pair_index.numel() > 0:
            pair_left = pair_index[:, 0]
            pair_right = pair_index[:, 1]
            pair_delta = self._minimum_image_delta(cell, pbc, positions[pair_right] - positions[pair_left])
            pair_distance = pair_delta.norm(dim=-1)

            min_distance = 0.75 * (covalent_radii[pair_left] + covalent_radii[pair_right])
            overlap = torch.relu(min_distance - pair_distance) / min_distance.clamp_min(1.0e-6)
            regularisation = regularisation + 10.0 * torch.mean(overlap ** 2)

        total_loss = total_loss + regularisation

        if target_vertex_fingerprints is not None and float(target_vertex_weight) > 0.0:
            target_vertex_2body, target_vertex_3body, target_vertex_4body = target_vertex_fingerprints
            vertex_loss = (
                torch.mean((vertex_2body - target_vertex_2body) ** 2)
                + torch.mean((vertex_3body - target_vertex_3body) ** 2)
                + torch.mean((vertex_4body - target_vertex_4body) ** 2)
            )
            total_loss = total_loss + float(target_vertex_weight) * vertex_loss

        if target_positions is not None and float(target_position_weight) > 0.0:
            if fixed_mask is None:
                position_loss = torch.mean((positions - target_positions) ** 2)
            else:
                movable_mask = ~fixed_mask
                if bool(movable_mask.any()):
                    position_loss = torch.mean((positions[movable_mask] - target_positions[movable_mask]) ** 2)
                else:
                    position_loss = torch.zeros((), dtype=torch.float32, device=self._device)
            total_loss = total_loss + float(target_position_weight) * position_loss

        return total_loss, fingerprint_loss

    def inverse_design(
        self,
        target_fingerprint: np.ndarray,
        atoms,
        fixed_atoms: np.ndarray,
        num_steps: int = 200,
        step_size: float = 1.0e-2,
        verbose: int = 0,
        target_atoms=None,
        fingerprint_loss_weight: float = 1.0,
        target_vertex_weight: float = 0.0,
        target_position_weight: float = 0.0,
        inverse_lr_decay_rate: Optional[float] = None,
        num_restarts: int = 1,
        restart_noise_scale: float = 0.0,
        step_observer: Optional[Callable[[dict[str, object]], None]] = None,
    ):
        self.eval()
        prepared = self.prepare_structure(atoms, include_targets=False)
        positions_initial = _float_tensor(prepared.positions, self._device)
        fixed_mask = torch.as_tensor(np.asarray(fixed_atoms, dtype=bool), dtype=torch.bool, device=self._device)
        movable_mask = ~fixed_mask
        target = self._project_fingerprint_targets(
            _float_tensor(target_fingerprint, self._device)
        )
        target_2body = target[:self.fingerprint_dim_2body]
        offset = self.fingerprint_dim_2body
        target_3body = target[offset:offset + self.fingerprint_dim_3body]
        offset += self.fingerprint_dim_3body
        target_4body = target[offset:offset + self.fingerprint_dim_4body]

        target_vertex_fingerprints = None
        target_positions = None
        if target_atoms is not None:
            prepared_target = self.prepare_structure(target_atoms, include_targets=False)
            if prepared_target.positions.shape != prepared.positions.shape:
                raise ValueError("target_atoms must have the same number of atoms as atoms")
            with torch.no_grad():
                (target_vertex_2body, target_vertex_3body, target_vertex_4body), _ = self._forward_prepared(
                    prepared_target,
                    return_vertices=True,
                )
            target_vertex_fingerprints = (
                target_vertex_2body.detach(),
                target_vertex_3body.detach(),
                target_vertex_4body.detach(),
            )
            target_positions = _float_tensor(prepared_target.positions, self._device)

        if inverse_lr_decay_rate is None:
            inverse_lr_decay_rate = self.lr_decay_rate

        reference_positions = positions_initial
        if target_atoms is not None and (
            float(target_vertex_weight) > 0.0 or float(target_position_weight) > 0.0
        ):
            reference_positions = None

        num_restarts = max(int(num_restarts), 1)
        restart_noise_scale = max(float(restart_noise_scale), 0.0)
        best_positions = positions_initial.clone()
        best_loss = float("inf")
        for restart_index in range(num_restarts):
            restart_positions = positions_initial.clone()
            if restart_index > 0 and restart_noise_scale > 0.0 and bool(movable_mask.any()):
                restart_generator = torch.Generator(device="cpu")
                restart_generator.manual_seed(self.seed + restart_index)
                restart_noise = torch.randn(
                    restart_positions.shape,
                    generator=restart_generator,
                    dtype=restart_positions.dtype,
                ).to(self._device)
                restart_positions = restart_positions + restart_noise_scale * restart_noise
                restart_positions[fixed_mask] = positions_initial[fixed_mask]

            if step_observer is not None:
                initial_atoms = atoms.copy()
                initial_atoms.set_positions(restart_positions.detach().cpu().numpy())
                step_observer(
                    {
                        "restart_index": int(restart_index),
                        "num_restarts": int(num_restarts),
                        "step": 0,
                        "num_steps": int(num_steps),
                        "is_initial_state": True,
                        "atoms": initial_atoms,
                    }
                )

            positions_parameter = nn.Parameter(restart_positions)
            optimiser = torch.optim.Adam([positions_parameter], lr=float(step_size))
            scheduler = None
            if float(inverse_lr_decay_rate) > 0.0:
                scheduler = torch.optim.lr_scheduler.ExponentialLR(
                    optimiser,
                    gamma=float(math.exp(-float(inverse_lr_decay_rate))),
                )

            for step in range(int(num_steps)):
                optimiser.zero_grad()
                candidate_positions = torch.where(
                    fixed_mask.unsqueeze(-1),
                    positions_initial,
                    positions_parameter,
                )
                total_loss, fingerprint_loss = self._positions_to_loss(
                    prepared,
                    candidate_positions,
                    target_2body,
                    target_3body,
                    target_4body,
                    reference_positions,
                    target_vertex_fingerprints=target_vertex_fingerprints,
                    target_positions=target_positions,
                    fixed_mask=fixed_mask,
                    fingerprint_loss_weight=fingerprint_loss_weight,
                    target_vertex_weight=target_vertex_weight,
                    target_position_weight=target_position_weight,
                )
                total_loss.backward()
                if positions_parameter.grad is not None:
                    positions_parameter.grad[fixed_mask] = 0.0
                torch.nn.utils.clip_grad_value_([positions_parameter], 1.0e-1)
                optimiser.step()
                if scheduler is not None:
                    scheduler.step()
                with torch.no_grad():
                    positions_parameter.data[fixed_mask] = positions_initial[fixed_mask]
                if step_observer is not None:
                    observed_atoms = atoms.copy()
                    observed_atoms.set_positions(
                        torch.where(
                            fixed_mask.unsqueeze(-1),
                            positions_initial,
                            positions_parameter,
                        )
                        .detach()
                        .cpu()
                        .numpy()
                    )
                    step_observer(
                        {
                            "restart_index": int(restart_index),
                            "num_restarts": int(num_restarts),
                            "step": int(step + 1),
                            "num_steps": int(num_steps),
                            "is_initial_state": False,
                            "atoms": observed_atoms,
                        }
                    )
                current_loss = float(total_loss.item())
                if current_loss < best_loss:
                    best_loss = current_loss
                    best_positions = candidate_positions.detach().clone()
                if verbose > 0 and ((step + 1) % max(int(num_steps) // 10, 1) == 0 or step == 0):
                    print(
                        f"restart={restart_index + 1:2d}/{num_restarts:2d} "
                        f"step={step + 1:4d} total_loss={current_loss:.6e} "
                        f"fingerprint_loss={float(fingerprint_loss.item()):.6e}"
                    )

            with torch.no_grad():
                final_positions = torch.where(
                    fixed_mask.unsqueeze(-1),
                    positions_initial,
                    positions_parameter,
                )
                final_total_loss, _ = self._positions_to_loss(
                    prepared,
                    final_positions,
                    target_2body,
                    target_3body,
                    target_4body,
                    reference_positions,
                    target_vertex_fingerprints=target_vertex_fingerprints,
                    target_positions=target_positions,
                    fixed_mask=fixed_mask,
                    fingerprint_loss_weight=fingerprint_loss_weight,
                    target_vertex_weight=target_vertex_weight,
                    target_position_weight=target_position_weight,
                )
                final_loss = float(final_total_loss.item())
            if final_loss < best_loss:
                best_loss = final_loss
                best_positions = final_positions.detach().clone()

        optimised = atoms.copy()
        optimised.set_positions(best_positions.detach().cpu().numpy())
        return wrap_atoms_to_unit_cell(optimised)
