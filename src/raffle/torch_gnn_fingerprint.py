from __future__ import annotations

from dataclasses import dataclass, field
import itertools
import math
import random
from typing import Callable, Iterable, Optional, Sequence, Tuple

import numpy as np
import torch
import torch.nn.functional as torch_functional
from torch import nn

from ase.constraints import FixAtoms  # added import

from .gnn_fingerprint import GNNFingerprint
from .structure_metrics import wrap_atoms_to_unit_cell

from .raffle import generator as _generator_class


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
    "transformer": "graph_transformer",
    "graph_transformer": "graph_transformer",
    "torch_gnn_transformer": "graph_transformer",
    "torch_gnn_graph_transformer": "graph_transformer",
    "transformer_coupled": "graph_transformer_coupled",
    "graph_transformer_coupled": "graph_transformer_coupled",
    "torch_gnn_transformer_coupled": "graph_transformer_coupled",
    "torch_gnn_graph_transformer_coupled": "graph_transformer_coupled",
    "graph_operator": "graph_operator",
    "graph_neural_operator": "graph_operator",
    "torch_gnn_graph_operator": "graph_operator",
    "torch_gnn_graph_neural_operator": "graph_operator",
    "graph_operator_coupled": "graph_operator_coupled",
    "graph_neural_operator_coupled": "graph_operator_coupled",
    "torch_gnn_graph_operator_coupled": "graph_operator_coupled",
    "torch_gnn_graph_neural_operator_coupled": "graph_operator_coupled",
    "kan": "multkan",
    "multkan": "multkan",
    "torch_gnn_kan": "multkan",
    "torch_gnn_multkan": "multkan",
    "kan_coupled": "multkan_coupled",
    "multkan_coupled": "multkan_coupled",
    "torch_gnn_kan_coupled": "multkan_coupled",
    "torch_gnn_multkan_coupled": "multkan_coupled",
}

COUPLED_ARCHITECTURES = {
    "attention_coupled": "attention",
    "graph_transformer_coupled": "graph_transformer",
    "graph_operator_coupled": "graph_operator",
    "multkan_coupled": "multkan",
}

FINGERPRINT_NEGATIVE_TAIL_BETA = 500.0


@dataclass(frozen=True)
class MultigraphTopology:
    symbols: Tuple[str, ...]
    species_index: np.ndarray
    atomic_numbers: np.ndarray
    covalent_radii: np.ndarray
    pair_image_shift: np.ndarray
    pair_target_species_index: np.ndarray
    pair_index: np.ndarray
    pair_type_index: np.ndarray
    pair_cutoff_weight_3body: np.ndarray
    pair_cutoff_weight_4body: np.ndarray
    angle_index: np.ndarray
    angle_species_index: np.ndarray
    triplet_index: np.ndarray
    triplet_pair_ids: np.ndarray
    triplet_center_index: np.ndarray
    dihedral_index: np.ndarray
    dihedral_pair_id: np.ndarray
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
    graph_stats: dict[str, float] = field(default_factory=dict)


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


def _edge_softmax(
    scores: torch.Tensor,
    dst: torch.Tensor,
    num_nodes: int,
    edge_weight: torch.Tensor,
) -> torch.Tensor:
    if scores.ndim == 1:
        scores = scores.unsqueeze(-1)
        squeeze_last = True
    else:
        squeeze_last = False

    expanded_dst = dst.unsqueeze(-1).expand(-1, scores.shape[-1])
    max_scores = torch.full(
        (int(num_nodes), scores.shape[-1]),
        torch.finfo(scores.dtype).min,
        dtype=scores.dtype,
        device=scores.device,
    )
    max_scores.scatter_reduce_(0, expanded_dst, scores, reduce="amax", include_self=True)

    scaled_scores = (
        torch.exp(scores - max_scores[dst])
        * edge_weight.clamp_min(1.0e-6).unsqueeze(-1)
    )
    score_sums = torch.zeros_like(max_scores)
    score_sums.index_add_(0, dst, scaled_scores)
    attention = scaled_scores / score_sums[dst].clamp_min(1.0e-6)
    if squeeze_last:
        return attention.squeeze(-1)
    return attention


class RadialKANLinear(nn.Module):
    def __init__(
        self,
        input_dim: int,
        output_dim: int,
        basis_size: int = 8,
        zero_init: bool = False,
    ):
        super().__init__()
        self.input_dim = int(input_dim)
        self.output_dim = int(output_dim)
        self.base = nn.Linear(self.input_dim, self.output_dim)
        self.coefficients = nn.Parameter(
            torch.zeros(self.input_dim, int(basis_size), self.output_dim, dtype=torch.float32)
        )
        self.register_buffer(
            "centers",
            torch.linspace(-1.5, 1.5, int(basis_size), dtype=torch.float32),
        )
        self.register_buffer(
            "widths",
            torch.full(
                (self.input_dim, int(basis_size)),
                2.5 / max(int(basis_size) - 1, 1),
                dtype=torch.float32,
            ),
        )
        if zero_init:
            _zero_init_linear(self.base)

    def forward(self, inputs: torch.Tensor) -> torch.Tensor:
        flattened = inputs.reshape(-1, self.input_dim)
        bounded = 1.5 * torch.tanh(flattened).unsqueeze(-1)
        basis = torch.exp(
            -(
                (bounded - self.centers.view(1, 1, -1))
                / self.widths.unsqueeze(0).clamp_min(1.0e-3)
            ) ** 2
        )
        spline = torch.einsum("bik,iko->bo", basis, self.coefficients)
        outputs = self.base(flattened) + spline
        return outputs.view(*inputs.shape[:-1], self.output_dim)


class GraphTransformerMessageLayer(nn.Module):
    def __init__(self, hidden_dim: int, edge_dim: int, global_dim: int):
        super().__init__()
        if hidden_dim % 4 == 0:
            self.num_heads = 4
        elif hidden_dim % 2 == 0:
            self.num_heads = 2
        else:
            self.num_heads = 1
        self.head_dim = hidden_dim // self.num_heads
        input_dim = hidden_dim + edge_dim + global_dim
        self.key = nn.Linear(input_dim, hidden_dim)
        self.value = nn.Linear(input_dim, hidden_dim)
        self.query = nn.Linear(hidden_dim + global_dim, hidden_dim)
        self.output = nn.Linear(hidden_dim, hidden_dim)
        self.feed_forward = nn.Sequential(
            nn.Linear(hidden_dim, 2 * hidden_dim),
            nn.SiLU(),
            nn.Linear(2 * hidden_dim, hidden_dim),
        )
        self.norm1 = nn.LayerNorm(hidden_dim)
        self.norm2 = nn.LayerNorm(hidden_dim)
        _zero_init_linear(self.output)
        _zero_init_linear(self.feed_forward[-1])

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

        keys = self.key(src_input).view(-1, self.num_heads, self.head_dim)
        values = self.value(src_input).view(-1, self.num_heads, self.head_dim)
        queries = self.query(dst_input).view(-1, self.num_heads, self.head_dim)
        scores = torch.sum(queries * keys, dim=-1) / math.sqrt(self.head_dim)
        attention = _edge_softmax(scores, dst, hidden.shape[0], edge_weight)

        aggregated = torch.zeros(
            hidden.shape[0],
            self.num_heads,
            self.head_dim,
            dtype=hidden.dtype,
            device=hidden.device,
        )
        aggregated.index_add_(0, dst, values * attention.unsqueeze(-1))
        attended = self.output(aggregated.reshape(hidden.shape[0], hidden.shape[1]))
        hidden = self.norm1(hidden + attended)
        return self.norm2(hidden + self.feed_forward(hidden))


class OperatorMessageLayer(nn.Module):
    def __init__(self, hidden_dim: int, edge_dim: int, global_dim: int):
        super().__init__()
        input_dim = 2 * hidden_dim + edge_dim + global_dim
        self.kernel = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.SiLU(),
            nn.Linear(hidden_dim, hidden_dim),
        )
        self.gate = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.SiLU(),
            nn.Linear(hidden_dim, 1),
        )
        self.update = nn.Sequential(
            nn.Linear(2 * hidden_dim + global_dim, hidden_dim),
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
        global_on_nodes = global_features.expand(hidden.shape[0], -1)
        if edge_index.numel() == 0:
            global_context = hidden.mean(dim=0, keepdim=True).expand_as(hidden)
            return hidden + self.update(torch.cat([hidden, global_context, global_on_nodes], dim=-1))

        src = edge_index[0]
        dst = edge_index[1]
        edge_global = global_features.expand(src.shape[0], -1)
        operator_input = torch.cat([hidden[src], hidden[dst], edge_attr, edge_global], dim=-1)
        kernel_values = self.kernel(operator_input)
        kernel_scale = torch_functional.softplus(self.gate(operator_input)).squeeze(-1)
        weighted_scale = kernel_scale * edge_weight.clamp_min(1.0e-6)

        aggregated = torch.zeros_like(hidden)
        aggregated.index_add_(0, dst, kernel_values * weighted_scale.unsqueeze(-1))
        normaliser = torch.zeros(hidden.shape[0], device=hidden.device, dtype=hidden.dtype)
        normaliser.index_add_(0, dst, weighted_scale)
        aggregated = aggregated / normaliser.clamp_min(1.0e-6).unsqueeze(-1)
        global_context = aggregated.mean(dim=0, keepdim=True).expand_as(hidden)
        return hidden + self.update(
            torch.cat([hidden, aggregated + global_context, global_on_nodes], dim=-1)
        )


class MultKANMessageLayer(nn.Module):
    def __init__(self, hidden_dim: int, edge_dim: int, global_dim: int):
        super().__init__()
        input_dim = hidden_dim + edge_dim + global_dim
        self.message = RadialKANLinear(input_dim, hidden_dim)
        self.gate = RadialKANLinear(input_dim, hidden_dim)
        self.update = RadialKANLinear(2 * hidden_dim, hidden_dim, zero_init=True)

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
        message = candidate * gate * edge_weight.unsqueeze(-1)

        aggregated = torch.zeros_like(hidden)
        aggregated.index_add_(0, dst, message)
        normaliser = torch.zeros(hidden.shape[0], device=hidden.device, dtype=hidden.dtype)
        normaliser.index_add_(0, dst, edge_weight.clamp_min(1.0e-6))
        aggregated = aggregated / normaliser.clamp_min(1.0e-6).unsqueeze(-1)
        return hidden + self.update(torch.cat([hidden, aggregated], dim=-1))


def _make_message_layer(layer_kind: str, hidden_dim: int, edge_dim: int, global_dim: int) -> nn.Module:
    if layer_kind == "residual":
        return ResidualMessageLayer(hidden_dim, edge_dim, global_dim)
    if layer_kind == "gated":
        return GatedMessageLayer(hidden_dim, edge_dim, global_dim)
    if layer_kind in {"attention", "attention_conservative"}:
        return AttentionMessageLayer(hidden_dim, edge_dim, global_dim)
    if layer_kind == "graph_transformer":
        return GraphTransformerMessageLayer(hidden_dim, edge_dim, global_dim)
    if layer_kind == "graph_operator":
        return OperatorMessageLayer(hidden_dim, edge_dim, global_dim)
    if layer_kind == "multkan":
        return MultKANMessageLayer(hidden_dim, edge_dim, global_dim)
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
        hidden_dim: int = 128,  # Now used as default/fallback
        hidden_dim_2body: Optional[int] = None,
        hidden_dim_3body: Optional[int] = None,
        hidden_dim_4body: Optional[int] = None,
        num_message_layers: int = 2,  # Now used as default/fallback
        num_message_layers_2body: Optional[int] = None,
        num_message_layers_3body: Optional[int] = None,
        num_message_layers_4body: Optional[int] = None,
        component_weight: Sequence[float] = (4.0, 1.0, 1.0),
        smooth_cutoff_width: float = 0.15,
        seed: int = 42,
        architecture: str = "residual",
        device: Optional[str] = None,
    ):
        super().__init__()
        self.species_list = [str(symbol).strip() for symbol in species_list]
        self.num_species = len(self.species_list)
        self.bond_cutoff = float(bond_cutoff)
        self.smooth_cutoff_width = float(smooth_cutoff_width)
        self.seed = int(seed)
        self._rng = random.Random(seed)
        self._device = torch.device(device or ("cuda" if torch.cuda.is_available() else "cpu"))
        self.architecture = ARCHITECTURE_ALIASES.get(str(architecture).strip(), str(architecture).strip())
        self._coupled_message_layer_kind = COUPLED_ARCHITECTURES.get(self.architecture)
        self._use_component_coupling = self._coupled_message_layer_kind is not None
        self._message_layer_kind = self._coupled_message_layer_kind or self.architecture

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

        self.reference_model = _generator_class.raffle_generator(
            seed=self.seed,
        )
        self.nbins = self.reference_model.distributions.get_nbins()

        self.fingerprint_dim_2body = int(self.nbins[0] * self.num_pairs)
        self.fingerprint_dim_3body = int(self.nbins[1] * self.num_species)
        self.fingerprint_dim_4body = int(self.nbins[2] * self.num_species)

        self.fingerprint_dim = self.fingerprint_dim_2body + self.fingerprint_dim_3body + self.fingerprint_dim_4body

        self.nbins = (
            self.fingerprint_dim_2body // max(self.num_pairs, 1),
            self.fingerprint_dim_3body // max(self.num_species, 1),
            self.fingerprint_dim_4body // max(self.num_species, 1),
        )
        self.width = self.reference_model.distributions.width
        self.sigma = self.reference_model.distributions.sigma
        self.cutoff_min = self.reference_model.distributions.cutoff_min
        self.cutoff_max = self.reference_model.distributions.cutoff_max
        self.radius_distance_tol = self.radius_distance_tol = self.reference_model.distributions.radius_distance_tol
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
        self._two_body_eta = 1.0 / (2.0 * (self.sigma[0] ** 2))
        pair_type_lookup = torch.full(
            (self.num_species, self.num_species),
            -1,
            dtype=torch.long,
        )
        for (left, right), pair_index in self._pair_to_index.items():
            pair_type_lookup[left, right] = pair_index
            pair_type_lookup[right, left] = pair_index
        self.register_buffer("_pair_type_lookup", pair_type_lookup)
        self.register_buffer(
            "_reference_2body_bin_weights",
            self._build_reference_2body_bin_weights(),
        )

        # Set branch-specific parameters with fallback to defaults
        hidden_dim_2body = hidden_dim_2body if hidden_dim_2body is not None else hidden_dim
        hidden_dim_3body = hidden_dim_3body if hidden_dim_3body is not None else hidden_dim
        hidden_dim_4body = hidden_dim_4body if hidden_dim_4body is not None else hidden_dim

        num_message_layers_2body = num_message_layers_2body if num_message_layers_2body is not None else num_message_layers
        num_message_layers_3body = num_message_layers_3body if num_message_layers_3body is not None else num_message_layers
        num_message_layers_4body = num_message_layers_4body if num_message_layers_4body is not None else num_message_layers

        # Create branches with independent parameters
        self.branch_2body = GraphBranch(
            node_dim=3 + self.num_species + 2,
            edge_dim=1,
            output_dim=self.fingerprint_dim_2body,
            hidden_dim=hidden_dim_2body,
            num_message_layers=num_message_layers_2body,
            global_dim=self.global_dim,
            message_layer_kind=self._message_layer_kind,
        )
        self.branch_3body = GraphBranch(
            node_dim=7 + 2 * self.num_species,
            edge_dim=1,
            output_dim=self.fingerprint_dim_3body,
            hidden_dim=hidden_dim_3body,
            num_message_layers=num_message_layers_3body,
            global_dim=self.global_dim,
            message_layer_kind=self._message_layer_kind,
        )
        self.branch_4body = GraphBranch(
            node_dim=6 + 3 * self.num_species,
            edge_dim=1,
            output_dim=self.fingerprint_dim_4body,
            hidden_dim=hidden_dim_4body,
            num_message_layers=num_message_layers_4body,
            global_dim=self.global_dim,
            message_layer_kind=self._message_layer_kind,
        )

        # Store configurations for reference
        self.hidden_dims = {
            '2body': hidden_dim_2body,
            '3body': hidden_dim_3body,
            '4body': hidden_dim_4body,
        }
        self.num_message_layers = {
            '2body': num_message_layers_2body,
            '3body': num_message_layers_3body,
            '4body': num_message_layers_4body,
        }

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

    @is_fitted.setter
    def is_fitted(self, value: bool) -> None:
        self._is_fitted = value

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
        projected_vertices = (
            self._pair_block_normalise(self._project_fingerprint_tensor(vertices[0])),
            self._project_fingerprint_tensor(vertices[1]),
            self._project_fingerprint_tensor(vertices[2]),
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

    def _fortran_shell_weight(self, distance: float, lower: float, upper: float) -> float:
        if distance < lower or distance > upper:
            return 0.0
        smooth_upper = min(self.cutoff_max[0], upper)
        span = max(smooth_upper - lower, 1.0e-8)
        phase = (2.0 * math.pi) * (distance - lower) / span
        return -0.5 * (math.cos(phase) - 1.0)

    def _build_topology(
        self,
        symbols: Sequence[str],
        positions: np.ndarray,
        cell: np.ndarray,
        pbc: np.ndarray,
    ) -> MultigraphTopology:
        key = self._topology_key(symbols)
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
        image_shifts = self._periodic_image_shifts(
            torch.as_tensor(cell, dtype=torch.float32),
            torch.as_tensor(pbc, dtype=torch.bool),
        ).detach().cpu().numpy().astype(np.int64)

        pair_image_shift = []
        pair_target_species_index = []
        pair_index = []
        pair_type_index = []
        pair_cutoff_weight_3body = []
        pair_cutoff_weight_4body = []
        centre_species_neighbours: dict[tuple[int, int], list[int]] = {}

        for center_atom in range(num_atoms):
            center_species = int(species_index[center_atom])
            for target_atom in range(num_atoms):
                target_species = int(species_index[target_atom])
                pair_species = tuple(sorted((center_species, target_species)))
                pair_type = int(self._pair_to_index[pair_species])
                pair_radius = 0.5 * (covalent_radii[center_atom] + covalent_radii[target_atom])
                tol_3_low = max(self.cutoff_min[0], self.radius_distance_tol[0] * pair_radius)
                tol_3_high = min(self.cutoff_max[0], self.radius_distance_tol[1] * pair_radius)
                tol_4_low = max(self.cutoff_min[0], self.radius_distance_tol[2] * pair_radius)
                tol_4_high = min(self.cutoff_max[0], self.radius_distance_tol[3] * pair_radius)

                for shift in image_shifts:
                    if center_atom == target_atom and int(shift[0]) == 0 and int(shift[1]) == 0 and int(shift[2]) == 0:
                        continue
                    shifted_target = positions[target_atom] + np.matmul(shift.astype(np.float32), cell)
                    delta = shifted_target - positions[center_atom]
                    distance = float(np.linalg.norm(delta))
                    if distance < self.cutoff_min[0] or distance > self.cutoff_max[0]:
                        continue

                    record_id = len(pair_index)
                    pair_image_shift.append((int(shift[0]), int(shift[1]), int(shift[2])))
                    pair_target_species_index.append(target_species)
                    pair_index.append((center_atom, target_atom))
                    pair_type_index.append(pair_type)
                    pair_cutoff_weight_3body.append(
                        self._fortran_shell_weight(distance, tol_3_low, tol_3_high)
                    )
                    pair_cutoff_weight_4body.append(
                        self._fortran_shell_weight(distance, tol_4_low, tol_4_high)
                    )
                    centre_species_neighbours.setdefault((center_atom, target_species), []).append(record_id)

        angle_index = []
        angle_species_index = []
        triplet_index = []
        triplet_pair_ids = []
        triplet_center_index = []
        dihedral_index = []
        dihedral_pair_id = []
        dihedral_species_index = []
        triplet_lookup: dict[tuple[int, int, int, int], int] = {}

        # For each triplet (defined by center_atom and two neighbors)
        # Find all other triplets with the same center atom and form dihedrals
        for center_atom in range(num_atoms):
            center_species = int(species_index[center_atom])
            for neighbour_species in range(self.num_species):
                all_records = centre_species_neighbours.get((center_atom, neighbour_species), [])
                if not all_records:
                    continue

                records_3body = [r for r in all_records if pair_cutoff_weight_3body[r] > 0.0]
                records_4body = [r for r in all_records if pair_cutoff_weight_4body[r] > 0.0]

                if len(records_3body) < 2:
                    continue

                # First, collect all triplets for this center atom
                triplets_for_center = []
                for left_idx in range(len(records_3body) - 1):
                    for right_idx in range(left_idx + 1, len(records_3body)):
                        pair_left = records_3body[left_idx]
                        pair_right = records_3body[right_idx]

                        triplet_key = (center_atom, neighbour_species, pair_left, pair_right)
                        triplet_id = triplet_lookup.get(triplet_key)
                        if triplet_id is None:
                            triplet_id = len(triplet_index)
                            triplet_lookup[triplet_key] = triplet_id
                            triplet_index.append((center_atom, neighbour_species, triplet_id))
                            triplet_pair_ids.append((pair_left, pair_right))
                            triplet_center_index.append(center_atom)

                        triplets_for_center.append(triplet_id)

                        # For angle indexing (3-body)
                        angle_index.append((pair_left, pair_right))
                        angle_species_index.append(center_species)

                # Now generate 4-body (dihedral) interactions
                # For each pair of triplets with the same center atom
                if len(triplets_for_center) >= 2:
                    for i in range(len(triplets_for_center) - 1):
                        for j in range(i + 1, len(triplets_for_center)):
                            triplet_1 = triplets_for_center[i]
                            triplet_2 = triplets_for_center[j]

                            # Each dihedral pairs two triplets
                            # The 4-body vertex will be computed from these two triplets
                            dihedral_index.append((triplet_1, triplet_2))

                            # The additional pair for 4-body (could be any record from 4body list)
                            # You might want to use the pair that defines the dihedral
                            for pair_l in records_4body:
                                dihedral_pair_id.append(pair_l)
                                dihedral_species_index.append(center_species)
                                break  # Only add one pair per dihedral? Or multiple?

        return MultigraphTopology(
            symbols=key,
            species_index=species_index,
            atomic_numbers=atomic_numbers,
            covalent_radii=covalent_radii,
            pair_image_shift=np.asarray(pair_image_shift, dtype=np.int64).reshape(-1, 3),
            pair_target_species_index=np.asarray(pair_target_species_index, dtype=np.int64),
            pair_index=np.asarray(pair_index, dtype=np.int64).reshape(-1, 2),
            pair_type_index=np.asarray(pair_type_index, dtype=np.int64),
            pair_cutoff_weight_3body=np.asarray(pair_cutoff_weight_3body, dtype=np.float32),
            pair_cutoff_weight_4body=np.asarray(pair_cutoff_weight_4body, dtype=np.float32),
            angle_index=np.asarray(angle_index, dtype=np.int64).reshape(-1, 2),
            angle_species_index=np.asarray(angle_species_index, dtype=np.int64),
            triplet_index=np.asarray(triplet_index, dtype=np.int64).reshape(-1, 3),
            triplet_pair_ids=np.asarray(triplet_pair_ids, dtype=np.int64).reshape(-1, 2),
            triplet_center_index=np.asarray(triplet_center_index, dtype=np.int64),
            dihedral_index=np.asarray(dihedral_index, dtype=np.int64).reshape(-1, 2),
            dihedral_pair_id=np.asarray(dihedral_pair_id, dtype=np.int64),
            dihedral_species_index=np.asarray(dihedral_species_index, dtype=np.int64),
        )

    def prepare_structure(self, atoms, include_targets: bool = True) -> PreparedStructure:
        symbols = tuple(str(symbol).strip() for symbol in atoms.get_chemical_symbols())
        positions = np.asarray(atoms.get_positions(), dtype=np.float32)
        cell = np.asarray(atoms.cell.array, dtype=np.float32)
        pbc = np.asarray(atoms.pbc, dtype=bool)
        topology = self._build_topology(symbols, positions, cell, pbc)
        targets = (None, None, None)
        if include_targets:
            targets = self._compute_reference_fingerprint_components(atoms)

        graph_stats = {
            "num_atoms": float(len(symbols)),
            "num_2body_edges": float(topology.pair_index.shape[0]),
            "num_3body_edges": float(topology.angle_index.shape[0]),
            "num_4body_edges": float(topology.dihedral_index.shape[0]),
        }

        return PreparedStructure(
            positions=positions,
            cell=cell,
            pbc=pbc,
            topology=topology,
            target_2body=None if targets[0] is None else np.asarray(targets[0], dtype=np.float32),
            target_3body=None if targets[1] is None else np.asarray(targets[1], dtype=np.float32),
            target_4body=None if targets[2] is None else np.asarray(targets[2], dtype=np.float32),
            graph_stats=graph_stats,
        )

    def _update_prepared_structure(
        self,
        prepared: PreparedStructure,
        positions: torch.Tensor,
    ) -> PreparedStructure:
        """
        Update a PreparedStructure with new positions, rebuilding the topology.
        This should be called after each optimization step.
        """
        # Convert positions back to numpy for topology building
        positions_np = positions.detach().cpu().numpy().astype(np.float32)

        # Rebuild topology with new positions
        symbols = prepared.topology.symbols
        cell = prepared.cell
        pbc = prepared.pbc
        new_topology = self._build_topology(symbols, positions_np, cell, pbc)

        # Create new PreparedStructure with updated topology
        return PreparedStructure(
            positions=positions_np,
            cell=cell,
            pbc=pbc,
            topology=new_topology,
            target_2body=prepared.target_2body,
            target_3body=prepared.target_3body,
            target_4body=prepared.target_4body,
            graph_stats=prepared.graph_stats,
        )

    def prepare_dataset(self, structures: Iterable, include_targets: bool = True) -> list[PreparedStructure]:
        return [self.prepare_structure(atoms, include_targets=include_targets) for atoms in structures]

    def _ensure_finite(self, label: str, values: np.ndarray) -> np.ndarray:
        array = np.asarray(values, dtype=np.float32)
        if not np.all(np.isfinite(array)):
            raise RuntimeError(f"Non-finite values detected in {label}.")
        return array

    def _compute_reference_fingerprint(self, atoms) -> np.ndarray:
        fp2, fp3, fp4 = self._compute_reference_fingerprint_components(atoms)
        return np.concatenate([fp2, fp3, fp4])

    def _compute_reference_fingerprint_components(self, atoms) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        # Get reference components
        fp2, fp3, fp4 = self.reference_model.distributions._compute_fingerprint_components(atoms)

        # Get species information
        symbols = tuple(str(symbol).strip() for symbol in atoms.get_chemical_symbols())
        present_species = sorted(set(symbols))
        num_present = len(present_species)
        num_total = self.num_species

        # Only reshape if needed (same logic for all fingerprint types)
        if num_present != num_total:
            # Helper to expand species-based fingerprints
            def expand_species_fingerprint(fp, nbins):
                fp_reshaped = fp.reshape(num_present, nbins)
                fp_full = np.zeros((num_total, nbins), dtype=np.float32)
                for i, species in enumerate(present_species):
                    fp_full[self._species_to_index[species]] = fp_reshaped[i]
                return fp_full.flatten()

            # Expand 3-body and 4-body fingerprints
            fp3 = expand_species_fingerprint(fp3, int(self.nbins[1]))
            fp4 = expand_species_fingerprint(fp4, int(self.nbins[2]))

            # Expand 2-body fingerprint (pairs, different mapping)
            nbins_2 = int(self.nbins[0])
            num_pairs_present = num_present * (num_present + 1) // 2
            num_pairs_total = self.num_pairs

            if num_pairs_present != num_pairs_total:
                fp2_reshaped = fp2.reshape(num_pairs_present, nbins_2)
                fp2_full = np.zeros((num_pairs_total, nbins_2), dtype=np.float32)

                # Map present species pairs to global pair indices
                for i, s1 in enumerate(present_species):
                    for j, s2 in enumerate(present_species[i:], start=i):
                        pair_key = tuple(sorted((self._species_to_index[s1], self._species_to_index[s2])))
                        global_pair_idx = self._pair_to_index[pair_key]
                        present_pair_idx = i * num_present - i*(i-1)//2 + (j - i)
                        fp2_full[global_pair_idx] = fp2_reshaped[present_pair_idx]

                fp2 = fp2_full.flatten()

        return (
            self._ensure_finite("2-body fingerprint", fp2),
            self._ensure_finite("3-body fingerprint", fp3),
            self._ensure_finite("4-body fingerprint", fp4),
        )

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

    def _wrap_positions_into_cell(
        self,
        positions: torch.Tensor,
        cell: torch.Tensor,
        pbc: torch.Tensor,
    ) -> torch.Tensor:
        if positions.numel() == 0 or not bool(pbc.any()):
            return positions
        inverse_cell = torch.linalg.inv(cell)
        positions_frac = positions @ inverse_cell
        wrapped_frac = positions_frac.clone()
        if pbc[0]:
            wrapped_frac[:, 0] = torch.remainder(wrapped_frac[:, 0], 1.0)
        if pbc[1]:
            wrapped_frac[:, 1] = torch.remainder(wrapped_frac[:, 1], 1.0)
        if pbc[2]:
            wrapped_frac[:, 2] = torch.remainder(wrapped_frac[:, 2], 1.0)
        return wrapped_frac @ cell

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

    def _build_reference_2body_bin_weights(self) -> torch.Tensor:
        weights = torch.ones_like(self._centers_2body)
        upper_start = self.cutoff_max[0] - 0.25
        lower_end = self.cutoff_min[0] + 0.25
        if upper_start < self.cutoff_max[0]:
            upper_weight = 0.5 * (
                1.0
                + torch.cos(
                    math.pi
                    * (self._centers_2body - upper_start)
                    / (self.cutoff_max[0] - upper_start)
                )
            )
            weights = torch.where(self._centers_2body > upper_start, upper_weight, weights)
        if lower_end > self.cutoff_min[0]:
            lower_weight = 0.5 * (
                1.0
                + torch.cos(
                    math.pi
                    * (self._centers_2body - lower_end)
                    / (lower_end - self.cutoff_min[0])
                )
            )
            weights = torch.where(self._centers_2body < lower_end, lower_weight, weights)
        return weights

    def _reference_distribution(
        self,
        values: torch.Tensor,
        centers: torch.Tensor,
        eta: float,
    ) -> torch.Tensor:
        if values.numel() == 0:
            return torch.zeros_like(centers)
        basis = torch.exp(-eta * (values.unsqueeze(-1) - centers.unsqueeze(0)) ** 2)
        histogram = basis.sum(dim=0)
        histogram = histogram * math.sqrt(eta / math.pi) / float(values.numel())
        return histogram

    def _pair_block_normalise(self, fingerprint: torch.Tensor) -> torch.Tensor:
        if self.num_pairs <= 0 or self.nbins[0] <= 0:
            return fingerprint
        if fingerprint.ndim == 1:
            blocks = fingerprint.reshape(self.num_pairs, self.nbins[0])
            block_sum = blocks.sum(dim=-1, keepdim=True)
            normalised = torch.where(
                block_sum > 1.0e-8,
                blocks / block_sum.clamp_min(1.0e-8),
                blocks,
            )
            return normalised.reshape(-1)
        if fingerprint.ndim == 2:
            blocks = fingerprint.reshape(fingerprint.shape[0], self.num_pairs, self.nbins[0])
            block_sum = blocks.sum(dim=-1, keepdim=True)
            normalised = torch.where(
                block_sum > 1.0e-8,
                blocks / block_sum.clamp_min(1.0e-8),
                blocks,
            )
            return normalised.reshape(fingerprint.shape[0], -1)
        raise ValueError(f"Unsupported 2-body fingerprint rank: {fingerprint.ndim}")

    def _periodic_image_shifts(self, cell: torch.Tensor, pbc: torch.Tensor) -> torch.Tensor:
        cell_lengths = torch.linalg.norm(cell, dim=1)
        shift_ranges = []
        for axis in range(3):
            if bool(pbc[axis].item()):
                axis_length = max(float(cell_lengths[axis].item()), 1.0e-8)
                max_shift = int(math.ceil(self.cutoff_max[0] / axis_length)) + 1
                shift_ranges.append(range(-max_shift, max_shift + 1))
            else:
                shift_ranges.append(range(0, 1))
        shifts = list(itertools.product(*shift_ranges))
        return torch.tensor(shifts, dtype=cell.dtype, device=cell.device)

    def _reference_style_2body_base(
        self,
        positions: torch.Tensor,
        cell: torch.Tensor,
        pbc: torch.Tensor,
        species_index: torch.Tensor,
    ) -> torch.Tensor:
        num_atoms = positions.shape[0]
        if num_atoms == 0:
            return torch.zeros(
                (0, self.fingerprint_dim_2body),
                dtype=positions.dtype,
                device=positions.device,
            )

        shifts = self._periodic_image_shifts(cell, pbc)
        shift_cart = shifts @ cell
        shifted_positions = positions.unsqueeze(0) + shift_cart.unsqueeze(1)
        zero_shift = torch.all(shifts == 0, dim=1)

        histogram = torch.zeros(
            (self.nbins[0], self.num_pairs),
            dtype=positions.dtype,
            device=positions.device,
        )
        centers = self._centers_2body.to(dtype=positions.dtype)
        bin_weights = self._reference_2body_bin_weights.to(dtype=positions.dtype)

        for center_index in range(num_atoms):
            center_species = species_index[center_index]
            pair_types = self._pair_type_lookup[center_species, species_index]
            distances = torch.linalg.norm(
                shifted_positions - positions[center_index].view(1, 1, 3),
                dim=-1,
            )
            valid = (distances >= self.cutoff_min[0]) & (distances <= self.cutoff_max[0])
            valid[zero_shift, center_index] = False

            for pair_type in range(self.num_pairs):
                neighbour_mask = pair_types == pair_type
                if not torch.any(neighbour_mask):
                    continue
                pair_values = distances[:, neighbour_mask]
                pair_valid = valid[:, neighbour_mask]
                values = pair_values[pair_valid]
                if values.numel() == 0:
                    continue
                histogram[:, pair_type] = histogram[:, pair_type] + self._reference_distribution(
                    values,
                    centers,
                    self._two_body_eta,
                )

        histogram = histogram / centers.pow(2).clamp_min(1.0e-8).unsqueeze(-1)
        histogram = histogram * bin_weights.unsqueeze(-1)
        block_sum = histogram.sum(dim=0, keepdim=True)
        histogram = torch.where(
            block_sum > 1.0e-8,
            histogram / block_sum.clamp_min(1.0e-8),
            histogram,
        )
        graph_base = histogram.transpose(0, 1).reshape(-1)
        return graph_base.unsqueeze(0).expand(num_atoms, -1)

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
        species_probabilities: Optional[torch.Tensor] = None,
    ):
        """Build graph tensors; if species_probabilities is given, use it for one‑hot features."""
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

        # ---- CENTERING ----
        # Compute center of mass (average of all atoms) and subtract from positions
        # used in feature construction. This removes the global translation component.
        com = positions.mean(dim=0, keepdim=True)
        positions_centered = positions - com  # used for features only
        # ------------------

        # Species one-hot: use probabilities if provided, else fallback to topology species
        if species_probabilities is not None:
            species_one_hot = species_probabilities.to(torch.float32)  # shape (num_atoms, num_species)
        else:
            species_one_hot = torch_functional.one_hot(species_index, num_classes=self.num_species).to(torch.float32)

        global_features = self._lattice_features(cell)

        # ---- 2-body node features (atom nodes) ----
        # Use centered positions, species one-hot, atomic number, and covalent radius
        atom_node_features = torch.cat(
            [
                positions_centered / self.bond_cutoff,      # centered position
                species_one_hot,
                atomic_numbers.unsqueeze(-1),
                covalent_radii.unsqueeze(-1),
            ],
            dim=-1,
        )

        num_atoms = positions.shape[0]
        pair_index = _long_tensor(topology.pair_index, device)
        atom_base = self._reference_style_2body_base(positions, cell, pbc, species_index)

        # ---- 2-body graph (pair nodes) ----
        if pair_index.numel() > 0:
            pair_left = _long_tensor(topology.pair_index[:, 0], device)
            pair_right = _long_tensor(topology.pair_index[:, 1], device)
            pair_shift = _float_tensor(topology.pair_image_shift, device)

            # Pair displacement and distance (use original positions for PBC correctness)
            pair_delta = positions[pair_right] + pair_shift @ cell - positions[pair_left]
            pair_distance = pair_delta.norm(dim=-1)
            pair_weight = self._smooth_cutoff(pair_distance)

            # Pair midpoint: use centered positions (with image shift)
            pair_midpoint = 0.5 * (positions_centered[pair_left] + positions_centered[pair_right] + pair_shift @ cell) / self.bond_cutoff
            pair_unit = pair_delta / pair_distance.clamp_min(1.0e-8).unsqueeze(-1)

            # Pair node features: midpoint (centered), unit vector, distance, species of both atoms
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

            atom_edge_index = pair_index.T
            atom_edge_attr = (pair_distance / self.bond_cutoff).unsqueeze(-1)
            atom_edge_weight = pair_weight
        else:
            pair_node_features = torch.zeros(
                (1, 7 + 2 * self.num_species),   # will be adjusted by node_dim
                dtype=torch.float32,
                device=device,
            )
            atom_edge_index = torch.zeros((2, 0), dtype=torch.long, device=device)
            atom_edge_attr = torch.zeros((0, 1), dtype=torch.float32, device=device)
            atom_edge_weight = torch.zeros((0,), dtype=torch.float32, device=device)

        # ---- 3-body graph (angle nodes) ----
        angle_index = _long_tensor(topology.angle_index, device)
        if angle_index.numel() > 0:
            # For angle features, we need pair_delta for all pair records
            pair_delta_all = positions[_long_tensor(topology.pair_index[:,1], device)] + _float_tensor(topology.pair_image_shift, device) @ cell - positions[_long_tensor(topology.pair_index[:,0], device)]
            angle_pairs = _long_tensor(topology.angle_index, device)
            left_pair_id = angle_pairs[:, 0]
            right_pair_id = angle_pairs[:, 1]
            angle_vec_i = pair_delta_all[left_pair_id]
            angle_vec_j = pair_delta_all[right_pair_id]
            angle_value = self._angle(angle_vec_i, angle_vec_j)

            weight_3body = _float_tensor(topology.pair_cutoff_weight_3body, device)
            angle_weight = (
                weight_3body[left_pair_id]
                * weight_3body[right_pair_id]
                / (
                    angle_vec_i.norm(dim=-1).pow(2).clamp_min(1.0e-8)
                    * angle_vec_j.norm(dim=-1).pow(2).clamp_min(1.0e-8)
                )
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

            # Build angle basis and embed species (center species)
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
            pair_base.index_add_(0, left_pair_id, angle_vector)
            pair_base.index_add_(0, right_pair_id, angle_vector)
        else:
            pair_edge_index = torch.zeros((2, 0), dtype=torch.long, device=device)
            pair_edge_attr = torch.zeros((0, 1), dtype=torch.float32, device=device)
            pair_edge_weight = torch.zeros((0,), dtype=torch.float32, device=device)
            pair_base = torch.zeros(
                (pair_node_features.shape[0], self.fingerprint_dim_3body),
                dtype=torch.float32,
                device=device,
            )

        # ---- 4-body graph (triplet nodes) ----
        triplet_index = _long_tensor(topology.triplet_index, device)
        if triplet_index.numel() > 0:
            # Recompute pair_delta for all pairs (use cached if available)
            pair_delta_all = positions[_long_tensor(topology.pair_index[:,1], device)] + _float_tensor(topology.pair_image_shift, device) @ cell - positions[_long_tensor(topology.pair_index[:,0], device)]
            triplet_pair_ids = _long_tensor(topology.triplet_pair_ids, device)
            triplet_centers = _long_tensor(topology.triplet_center_index, device)
            left_pair_id = triplet_pair_ids[:, 0]
            right_pair_id = triplet_pair_ids[:, 1]
            vector_ij = pair_delta_all[left_pair_id]
            vector_jk = pair_delta_all[right_pair_id]
            distance_ij = vector_ij.norm(dim=-1)
            distance_jk = vector_jk.norm(dim=-1)
            triplet_angle = self._angle(vector_ij, vector_jk)

            # Triplet centroid: use centered positions
            triplet_centroid = (
                positions_centered[triplet_centers]
                + positions_centered[_long_tensor(topology.pair_index[:,1], device)][left_pair_id]
                + positions_centered[_long_tensor(topology.pair_index[:,1], device)][right_pair_id]
            ) / (3.0 * self.bond_cutoff)

            # Get species for the three atoms in the triplet
            left_species = species_one_hot[_long_tensor(topology.pair_target_species_index, device)[left_pair_id]]
            right_species = species_one_hot[_long_tensor(topology.pair_target_species_index, device)[right_pair_id]]

            triplet_node_features = torch.cat(
                [
                    triplet_centroid,
                    (distance_ij / self.bond_cutoff).unsqueeze(-1),
                    (distance_jk / self.bond_cutoff).unsqueeze(-1),
                    (triplet_angle / math.pi).unsqueeze(-1),
                    species_one_hot[triplet_centers],
                    left_species,
                    right_species,
                ],
                dim=-1,
            )
        else:
            triplet_node_features = torch.zeros(
                (1, 6 + 3 * self.num_species),
                dtype=torch.float32,
                device=device,
            )

        # ---- 4-body dihedral part (edges between triplets) ----
        dihedral_index = _long_tensor(topology.dihedral_index, device)
        if dihedral_index.numel() > 0:
            pair_delta_all = positions[_long_tensor(topology.pair_index[:,1], device)] + _float_tensor(topology.pair_image_shift, device) @ cell - positions[_long_tensor(topology.pair_index[:,0], device)]
            triplet_pair_ids = _long_tensor(topology.triplet_pair_ids, device)
            pair_l = _long_tensor(topology.dihedral_pair_id, device)
            triplet_ids = dihedral_index[:, 0]
            pair_i = triplet_pair_ids[triplet_ids, 0]
            pair_k = triplet_pair_ids[triplet_ids, 1]
            vector_ij = pair_delta_all[pair_i]
            vector_jk = pair_delta_all[pair_k]
            vector_kl = pair_delta_all[pair_l]
            dihedral_value = self._improper_dihedral(vector_ij, vector_jk, vector_kl)

            weight_3body = _float_tensor(topology.pair_cutoff_weight_3body, device)
            weight_4body = _float_tensor(topology.pair_cutoff_weight_4body, device)
            dihedral_weight = (
                weight_3body[pair_i]
                * weight_3body[pair_k]
                * weight_4body[pair_l]
                / (
                    vector_ij.norm(dim=-1).pow(2).clamp_min(1.0e-8)
                    * vector_jk.norm(dim=-1).pow(2).clamp_min(1.0e-8)
                    * vector_kl.norm(dim=-1).pow(2).clamp_min(1.0e-8)
                )
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
            triplet_base.index_add_(0, dihedral_index[:, 0], dihedral_vector)
            triplet_base.index_add_(0, dihedral_index[:, 1], dihedral_vector)
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
        species_probabilities: Optional[torch.Tensor] = None,   # NEW
        return_vertices: bool = False,
    ):
        graph = self._build_multigraph_tensors(
            prepared,
            positions_override=positions_override,
            species_probabilities=species_probabilities,
        )
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

        # Normalize by target energy so sparse, low-amplitude spectra do not
        # make the trivial all-zero predictor artificially cheap.
        def relative_component_mse(predicted: torch.Tensor, target: torch.Tensor) -> torch.Tensor:
            target_energy = torch.mean(target ** 2).clamp_min(1.0e-8)
            return torch.mean((predicted - target) ** 2) / target_energy

        loss_2body = relative_component_mse(predicted_2body, target_2body)
        loss_3body = relative_component_mse(predicted_3body, target_3body)
        loss_4body = relative_component_mse(predicted_4body, target_4body)
        return (
            self._component_weight_tensor[0] * loss_2body
            + self._component_weight_tensor[1] * loss_3body
            + self._component_weight_tensor[2] * loss_4body
        )

    def recommended_parameter_count_target(self) -> int:
        # Keep model size proportional to descriptor dimensionality while
        # avoiding unnecessary growth for small single-species systems.
        descriptor_dim = max(int(self.fingerprint_dim), 1)
        return int(max(25_000, min(400_000, round(24.0 * descriptor_dim))))

    def parameter_counts(self) -> dict[str, int]:
        total = int(sum(parameter.numel() for parameter in self.parameters()))
        trainable = int(
            sum(parameter.numel() for parameter in self.parameters() if parameter.requires_grad)
        )
        return {
            "recommended_target": int(self.recommended_parameter_count_target()),
            "total": total,
            "trainable": trainable,
        }

    def graph_statistics(self, atoms) -> dict[str, float]:
        prepared = self.prepare_structure(atoms, include_targets=False)
        topology = prepared.topology
        atom_count = max(float(len(prepared.positions)), 1.0)
        return {
            **prepared.graph_stats,
            "mean_2body_edges_per_atom": float(topology.pair_index.shape[0]) / atom_count,
            "mean_3body_edges_per_atom": float(topology.angle_index.shape[0]) / atom_count,
            "mean_4body_edges_per_atom": float(topology.dihedral_index.shape[0]) / atom_count,
        }

    def descriptor_surrogate_agreement(self, atoms) -> dict[str, float]:
        prepared = self.prepare_structure(atoms, include_targets=True)
        graph = self._build_multigraph_tensors(prepared)
        target_2body = _float_tensor(prepared.target_2body, self._device)
        target_3body = _float_tensor(prepared.target_3body, self._device)
        target_4body = _float_tensor(prepared.target_4body, self._device)

        base_2body = graph["atom_base"].mean(dim=0)
        base_3body = graph["pair_base"].mean(dim=0)
        base_4body = graph["triplet_base"].mean(dim=0)

        def _mae(prediction: torch.Tensor, target: torch.Tensor) -> float:
            return float(torch.mean(torch.abs(prediction - target)).item())

        def _rmse(prediction: torch.Tensor, target: torch.Tensor) -> float:
            return float(torch.sqrt(torch.mean((prediction - target) ** 2)).item())

        return {
            "mae_2body": _mae(base_2body, target_2body),
            "mae_3body": _mae(base_3body, target_3body),
            "mae_4body": _mae(base_4body, target_4body),
            "rmse_2body": _rmse(base_2body, target_2body),
            "rmse_3body": _rmse(base_3body, target_3body),
            "rmse_4body": _rmse(base_4body, target_4body),
        }

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
        learning_rate: float = 1.0e-3,
        lr_decay_rate: float = 1.0e-2,
        augment_structures: Optional[Sequence] = None,
        verbose: int = 0,
        reset_optimiser: bool = False,
        recalibrate_base: Optional[bool] = None,
    ) -> list[float]:
        """Train the surrogate model on reference RAFFLE fingerprints.

        Procedure:
        1. Build a training set from ``structures`` plus optional
           ``augment_structures`` and precompute prepared graph entries.
        2. Report parameter-count diagnostics.
        3. Optionally recalibrate per-branch affine base corrections
           (enabled by default on first fit).
        4. Record an initial full-dataset loss before any gradient updates.
        5. Create (or reuse) Adam + exponential LR scheduler.
        6. Run shuffled mini-batch training for ``num_epochs`` with gradient
           clipping and append a full-dataset evaluation loss after each epoch.

        Args:
            structures: Base structures used for supervised fitting.
            num_epochs: Number of training epochs.
            batch_size: Mini-batch size used in each epoch.
            augment_structures: Optional additional structures concatenated to
                ``structures`` before training.
            verbose: Logging level. ``0`` prints parameter counts only,
                values ``> 0`` also print periodic epoch losses.
            reset_optimiser: If ``True``, recreate optimizer and scheduler even
                if previous state exists.
            recalibrate_base: Controls base-affine recalibration. If ``None``,
                recalibration runs only when the model has not been fitted yet.

        Returns:
            Loss history where index ``0`` is the pre-training baseline and
            each subsequent value is the post-epoch dataset loss.
        """
        combined_structures = list(structures)
        if augment_structures:
            combined_structures.extend(list(augment_structures))
        entries = self.prepare_dataset(combined_structures, include_targets=True)
        counts = self.parameter_counts()
        if verbose >= 0:
            print(
                "[torch_gnn] parameter_count_target="
                f"{counts['recommended_target']:,} total={counts['total']:,} "
                f"trainable={counts['trainable']:,}"
            )
        if recalibrate_base is None:
            recalibrate_base = not self._is_fitted
        if recalibrate_base:
            self._calibrate_base_affine(entries)

        initial_loss = self._evaluate_entries(entries)
        history = [initial_loss]

        if reset_optimiser or self._optimiser is None:
            self._optimiser = torch.optim.Adam(self.parameters(), lr=learning_rate)
            self._scheduler = torch.optim.lr_scheduler.ExponentialLR(
                self._optimiser,
                gamma=float(math.exp(-lr_decay_rate)),
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
        fixed_mask: Optional[torch.Tensor] = None,
        species_probabilities: Optional[torch.Tensor] = None,   # NEW
        minimum_distance_scale: float = 0.75,
        repulsion_max: float = 100.0,
        repulsion_cutoff_scale: float = 1.0,
    ) -> dict[str, torch.Tensor]:
        (
            (vertex_2body, vertex_3body, vertex_4body),
            (prediction_2body, prediction_3body, prediction_4body),
        ) = self._forward_prepared(
            prepared,
            positions_override=positions,
            species_probabilities=species_probabilities,   # pass through
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
        # # # Use log1p to compress large values
        # # fingerprint_loss_compressed = torch.log1p(fingerprint_loss)
        # # # Or use a smooth cap:
        # fingerprint_loss_compressed = fingerprint_loss / (1.0 + fingerprint_loss / 1e3)  # saturates at 1e3

        topology = prepared.topology
        pair_index = _long_tensor(topology.pair_index, self._device)
        pair_image_shift = _float_tensor(topology.pair_image_shift, self._device)   # NEW
        cell = _float_tensor(prepared.cell, self._device)
        pbc = torch.as_tensor(prepared.pbc, dtype=torch.bool, device=self._device)
        covalent_radii = _float_tensor(topology.covalent_radii, self._device)

        repulsion_loss = torch.zeros((), dtype=torch.float32, device=self._device)
        if pair_index.numel() > 0:
            pair_left = pair_index[:, 0]
            pair_right = pair_index[:, 1]

            # Exclude self-pairs (atoms with themselves)
            is_self = (pair_left == pair_right)
            if torch.all(is_self):
                repulsion_loss = torch.zeros((), dtype=torch.float32, device=self._device)
            else:
                non_self_mask = ~is_self
                pair_left = pair_left[non_self_mask]
                pair_right = pair_right[non_self_mask]
                pair_image_shift = pair_image_shift[non_self_mask]

                # Use the stored image shift to compute displacement
                pair_delta = positions[pair_right] + pair_image_shift @ cell - positions[pair_left]
                pair_distance = pair_delta.norm(dim=-1)

                # print the minimum pair distance for debugging
                min_distance = pair_distance.min().item()
                if min_distance < 1e-3:
                    print(f"[DEBUG] Minimum pair distance: {min_distance:.6e}")

                covalent_sum = covalent_radii[pair_left] + covalent_radii[pair_right]
                r_min = float(minimum_distance_scale) * covalent_sum
                r_cutoff = float(repulsion_cutoff_scale) * covalent_sum

                is_active = pair_distance < r_cutoff
                r_ratio = r_min / pair_distance.clamp_min(1e-6)
                repulsion_value = (r_ratio ** 2) * (1 - pair_distance / r_cutoff.clamp_min(1e-6)) ** 2
                repulsion_value = torch.clamp(repulsion_value, max=float(repulsion_max))
                repulsion_value = repulsion_value * is_active.float()
                repulsion_loss = repulsion_value.sum()   # sum, not mean

        return {
            "fingerprint_loss": fingerprint_loss,
            "repulsion_loss": repulsion_loss,
        }

    def _discretize_species(
        self,
        probabilities: torch.Tensor,
        mode: str = 'argmax'
    ) -> torch.Tensor:
        """Convert continuous probabilities to integer species indices."""
        if mode == 'argmax':
            return torch.argmax(probabilities, dim=-1)
        elif mode == 'sample':
            # sample from multinomial distribution
            # probabilities shape (num_atoms, num_species)
            dist = torch.distributions.Categorical(probs=probabilities)
            return dist.sample()
        else:
            raise ValueError(f"Unknown species discretization mode: {mode}")

    def _rebuild_topology_with_species(
        self,
        prepared: PreparedStructure,
        positions: torch.Tensor,
        species_indices: torch.Tensor,
    ) -> PreparedStructure:
        """Rebuild topology using new positions and species indices."""
        positions_np = positions.detach().cpu().numpy().astype(np.float32)
        symbols = [self.species_list[int(idx)] for idx in species_indices.cpu().numpy()]
        cell = prepared.cell
        pbc = prepared.pbc
        new_topology = self._build_topology(symbols, positions_np, cell, pbc)
        return PreparedStructure(
            positions=positions_np,
            cell=cell,
            pbc=pbc,
            topology=new_topology,
            target_2body=prepared.target_2body,
            target_3body=prepared.target_3body,
            target_4body=prepared.target_4body,
            graph_stats=prepared.graph_stats,
        )

    def inverse_design(
        self,
        target_fingerprint: np.ndarray,
        atoms,
        fixed_atoms: Optional[np.ndarray] = None,           # deprecated
        num_steps: int = 200,
        step_size: float = 1.0e-2,
        verbose: int = 0,
        fingerprint_loss_weight: float = 1.0,
        inverse_lr_decay_rate: Optional[float] = None,
        num_restarts: int = 1,
        restart_noise_scale: float = 0.0,
        repulsion_weight: float = 10.0,
        minimum_distance_scale: float = 0.75,
        coordinate_clip_value: Optional[float] = None,
        wrap_positions_to_cell: bool = True,
        step_observer: Optional[Callable[[dict[str, object]], None]] = None,
        use_augmented_lagrangian: bool = True,
        initial_multiplier: float = 0.0,
        penalty_parameter: float = 1.0,
        multiplier_increase_factor: float = 2.0,
        max_multiplier: float = 1e6,
        constraint_tolerance: float = 1e-6,
        update_topology_every_n_steps: int = 10,
        # NEW parameters for species optimisation
        optimize_species: bool = False,
        species_optimization_mode: str = 'argmax',          # 'argmax' or 'sample'
        species_learning_rate: float = 1.0e-2,              # LR for species logits
        fixed_species: Optional[np.ndarray] = None,        # boolean mask for fixed species
        species_initial: Optional[np.ndarray] = None,      # initial species indices (optional)
        return_trajectory: bool = False,
    ):
        # --- Prepare initial structure ---
        self.eval()
        working_atoms = atoms.copy()

        # --- Handle constraints and fixed_atoms ---
        if fixed_atoms is not None:
            import warnings
            warnings.warn("`fixed_atoms` is deprecated; use `constraints` with ASE FixAtoms.", DeprecationWarning)
            # convert fixed_atoms to FixAtoms constraint
            indices = np.where(fixed_atoms)[0].tolist()
            fix_constraint = FixAtoms(indices=indices)
            working_atoms.set_constraint(fix_constraint)

        # Extract fixed position mask from constraints on the atoms object
        fixed_pos_mask = torch.zeros(len(working_atoms), dtype=torch.bool, device=self._device)
        if working_atoms.constraints is not None:
            for con in working_atoms.constraints:
                if isinstance(con, FixAtoms):
                    indices = con.get_indices()
                    fixed_pos_mask[indices] = True

        # --- Prepare initial structure ---
        self.eval()
        # If species_initial is provided, set them
        if species_initial is not None:
            symbols = [self.species_list[int(idx)] for idx in species_initial]
            working_atoms.set_chemical_symbols(symbols)

        prepared = self.prepare_structure(working_atoms, include_targets=False)
        positions_initial = _float_tensor(prepared.positions, self._device)
        cell = _float_tensor(prepared.cell, self._device)
        pbc = torch.as_tensor(prepared.pbc, dtype=torch.bool, device=self._device)
        fixed_mask = fixed_pos_mask
        movable_mask = ~fixed_mask

        # --- Setup target fingerprints ---
        target = self._project_fingerprint_targets(
            _float_tensor(target_fingerprint, self._device)
        )
        target_2body = target[:self.fingerprint_dim_2body]
        offset = self.fingerprint_dim_2body
        target_3body = target[offset:offset + self.fingerprint_dim_3body]
        offset += self.fingerprint_dim_3body
        target_4body = target[offset:offset + self.fingerprint_dim_4body]

        num_restarts = max(int(num_restarts), 1)
        restart_noise_scale = max(float(restart_noise_scale), 0.0)

        # Store best restart information (by final loss)
        best_restart_loss = float("inf")
        best_restart_trajectory = None

        # --- Initialise augmented Lagrangian ---
        if use_augmented_lagrangian:
            lambda_mult = float(initial_multiplier)
            mu = float(penalty_parameter)
            effective_repulsion_weight = 0.0
        else:
            lambda_mult = 0.0
            mu = 0.0
            effective_repulsion_weight = repulsion_weight

        update_topology_every_n_steps = max(int(update_topology_every_n_steps), 1)

        # --- Prepare species optimisation ---
        if optimize_species:
            # Create species logits parameter
            # initial species from working_atoms
            initial_symbols = working_atoms.get_chemical_symbols()
            initial_indices = torch.tensor(
                [self._species_to_index[sym] for sym in initial_symbols],
                dtype=torch.long,
                device=self._device
            )
            # one-hot as initial probabilities
            init_prob = torch_functional.one_hot(initial_indices, num_classes=self.num_species).float()
            # add small noise to break symmetry if desired
            init_logits = torch.log(init_prob + 1e-8)
            species_logits = nn.Parameter(init_logits, requires_grad=True)
            # fixed species mask
            if fixed_species is not None:
                fixed_species_mask = torch.as_tensor(fixed_species, dtype=torch.bool, device=self._device)
                # We'll zero out gradients for fixed species
            else:
                fixed_species_mask = torch.zeros(len(initial_symbols), dtype=torch.bool, device=self._device)
        else:
            species_logits = None
            fixed_species_mask = None

        # --- Restart loop ---
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

            if optimize_species and restart_index > 0:
                # re‑initialize species logits? not necessary, we keep global parameter
                pass

            positions_parameter = nn.Parameter(restart_positions)

            # --- Prepare parameters ---
            params = [positions_parameter] if not optimize_species else [positions_parameter, species_logits]
            # if optimizing species, we could use a different LR for species, but we use same for simplicity
            # Actually we can set different LRs by passing param groups
            if optimize_species:
                optimiser = torch.optim.Adam([
                    {'params': [positions_parameter], 'lr': float(step_size)},
                    {'params': [species_logits], 'lr': float(species_learning_rate)}
                ])
            else:
                optimiser = torch.optim.Adam(params, lr=float(step_size))
            scheduler = None
            if inverse_lr_decay_rate is not None and inverse_lr_decay_rate > 0.0:
                scheduler = torch.optim.lr_scheduler.ExponentialLR(
                    optimiser,
                    gamma=float(math.exp(-float(inverse_lr_decay_rate))),
                )

            restart_trajectory = [working_atoms.copy()] if return_trajectory else None

            # --- Step loop ---
            for step in range(int(num_steps)):
                # Update topology periodically if positions or species have changed
                if step % update_topology_every_n_steps == 0 and step > 0:
                    # Get current positions and species (discretized)
                    current_positions = torch.where(
                        fixed_mask.unsqueeze(-1),
                        positions_initial,
                        positions_parameter,
                    )
                    if optimize_species:
                        probs = torch_functional.softmax(species_logits, dim=-1)
                        species_idx = self._discretize_species(probs, mode='argmax')  # use argmax for topology
                    else:
                        # use topology's species
                        species_idx = torch.tensor(prepared.topology.species_index, device=self._device)
                    prepared = self._rebuild_topology_with_species(
                        prepared,
                        current_positions,
                        species_idx
                    )
                    # update cell in case it changed (shouldn't)
                    cell = _float_tensor(prepared.cell, self._device)
                    pbc = torch.as_tensor(prepared.pbc, dtype=torch.bool, device=self._device)

                optimiser.zero_grad()
                candidate_positions = torch.where(
                    fixed_mask.unsqueeze(-1),
                    positions_initial,
                    positions_parameter,
                )

                # Compute species probabilities if optimizing
                if optimize_species:
                    probs = torch_functional.softmax(species_logits, dim=-1)
                    # apply fixed species mask: set gradients to zero later
                else:
                    probs = None

                # Get loss components
                loss_components = self._positions_to_loss(
                    prepared,
                    candidate_positions,
                    target_2body,
                    target_3body,
                    target_4body,
                    fixed_mask=fixed_mask,
                    species_probabilities=probs,   # pass probabilities
                    minimum_distance_scale=minimum_distance_scale,
                )

                fingerprint_loss = loss_components["fingerprint_loss"]
                repulsion_loss = loss_components["repulsion_loss"]

                total_loss = self._compute_inverse_design_loss(
                    fingerprint_loss,
                    repulsion_loss,
                    fingerprint_loss_weight,
                    effective_repulsion_weight,
                    use_augmented_lagrangian,
                    lambda_mult,
                    mu,
                    num_atoms=len(working_atoms)
                )

                total_loss.backward()
                # Zero gradients for fixed atoms positions
                if positions_parameter.grad is not None:
                    positions_parameter.grad[fixed_mask] = 0.0
                # Zero gradients for fixed species
                if optimize_species and species_logits.grad is not None:
                    # if fixed_species_mask is provided, set grad to 0 for those atoms
                    if fixed_species_mask is not None:
                        species_logits.grad[fixed_species_mask] = 0.0

                torch.nn.utils.clip_grad_value_([positions_parameter], 1.0e-1)
                if optimize_species:
                    torch.nn.utils.clip_grad_value_([species_logits], 1.0e-1)
                optimiser.step()
                if scheduler is not None:
                    scheduler.step()

                # Update augmented Lagrangian multipliers after step
                if use_augmented_lagrangian:
                    with torch.no_grad():
                        # recompute repulsion with updated positions and species
                        updated_components = self._positions_to_loss(
                            prepared,
                            positions_parameter,
                            target_2body,
                            target_3body,
                            target_4body,
                            fixed_mask=fixed_mask,
                            species_probabilities=probs,
                            minimum_distance_scale=minimum_distance_scale,
                        )
                        rep_val = updated_components["repulsion_loss"].detach().item()
                        if rep_val > 0.0:
                            lambda_mult = max(0.0, lambda_mult + mu * rep_val)
                            lambda_mult = min(lambda_mult, max_multiplier)
                            if rep_val > constraint_tolerance:
                                mu = min(mu * multiplier_increase_factor, 1e6)

                # --- Wrapping and clipping (unchanged) ---
                with torch.no_grad():
                    positions_parameter.data[fixed_mask] = positions_initial[fixed_mask]
                    if bool(wrap_positions_to_cell) and bool(pbc.any()) and bool(movable_mask.any()):
                        wrapped_positions = self._wrap_positions_into_cell(
                            positions_parameter.data,
                            cell,
                            pbc,
                        )
                        positions_parameter.data[movable_mask] = wrapped_positions[movable_mask]
                        positions_parameter.data[fixed_mask] = positions_initial[fixed_mask]
                    if coordinate_clip_value is not None and bool(movable_mask.any()):
                        max_delta = float(coordinate_clip_value)
                        movable_delta = positions_parameter.data[movable_mask] - positions_initial[movable_mask]
                        if bool(wrap_positions_to_cell) and bool(pbc.any()):
                            movable_delta = self._minimum_image_delta(
                                cell,
                                pbc,
                                movable_delta,
                            )
                        positions_parameter.data[movable_mask] = positions_initial[movable_mask] + movable_delta.clamp(
                            min=-max_delta,
                            max=max_delta,
                        )
                    if bool(wrap_positions_to_cell) and bool(pbc.any()) and bool(movable_mask.any()):
                        wrapped_positions = self._wrap_positions_into_cell(
                            positions_parameter.data,
                            cell,
                            pbc,
                        )
                        positions_parameter.data[movable_mask] = wrapped_positions[movable_mask]
                        positions_parameter.data[fixed_mask] = positions_initial[fixed_mask]
                    current_positions = torch.where(
                        fixed_mask.unsqueeze(-1),
                        positions_initial,
                        positions_parameter,
                    )
                    # recompute loss for observer and best tracking
                    current_loss_components = self._positions_to_loss(
                        prepared,
                        current_positions,
                        target_2body,
                        target_3body,
                        target_4body,
                        fixed_mask=fixed_mask,
                        species_probabilities=probs,
                        minimum_distance_scale=minimum_distance_scale,
                    )
                    current_total_loss = self._compute_inverse_design_loss(
                        current_loss_components["fingerprint_loss"],
                        current_loss_components["repulsion_loss"],
                        fingerprint_loss_weight=fingerprint_loss_weight,
                        repulsion_weight=effective_repulsion_weight,
                        use_augmented_lagrangian=use_augmented_lagrangian,
                        lambda_mult=lambda_mult,
                        mu=mu,
                        num_atoms=len(working_atoms)
                    )
                if step_observer is not None:
                    observer_positions = current_positions
                    observed_atoms = working_atoms.copy()
                    observed_atoms.set_positions(observer_positions.detach().cpu().numpy())
                    current_learning_rate = float(optimiser.param_groups[0]["lr"])
                    step_observer(
                        {
                            "restart_index": int(restart_index),
                            "num_restarts": int(num_restarts),
                            "step": int(step + 1),
                            "num_steps": int(num_steps),
                            "is_initial_state": False,
                            "atoms": observed_atoms,
                            "learning_rate": current_learning_rate,
                            "total_loss": float(current_total_loss.item()),
                            "fingerprint_loss": float(
                                current_loss_components["fingerprint_loss"].item()
                            ),
                            "repulsion_loss": float(
                                current_loss_components["repulsion_loss"].item()
                            ),
                            "minimum_distance_scale": float(minimum_distance_scale),
                            "lambda_multiplier": float(lambda_mult) if use_augmented_lagrangian else 0.0,
                            "penalty_parameter": float(mu) if use_augmented_lagrangian else 0.0,
                        }
                    )
                if verbose > 0 and ((step + 1) % max(int(num_steps) // 10, 1) == 0 or step == 0):
                    extra = ""
                    if use_augmented_lagrangian:
                        extra = f" λ={lambda_mult:.2e} μ={mu:.2e}"
                    print(
                        f"restart={restart_index + 1:2d}/{num_restarts:2d} "
                        f"step={step + 1:4d} total_loss={current_total_loss:.6e} "
                        f"fingerprint_loss={float(current_loss_components['fingerprint_loss'].item()):.6e} "
                        f"repulsion_loss={float(current_loss_components['repulsion_loss'].item()):.6e}{extra}"
                    )

                if return_trajectory:
                    traj_atoms = working_atoms.copy()
                    traj_atoms.set_positions(current_positions.detach().cpu().numpy())
                    if optimize_species:
                        # Get current species from probabilities
                        current_probs = torch_functional.softmax(species_logits, dim=-1)
                        current_species_idx = torch.argmax(current_probs, dim=-1)
                        current_symbols = [self.species_list[int(idx)] for idx in current_species_idx.cpu().numpy()]
                        traj_atoms.set_chemical_symbols(current_symbols)
                    restart_trajectory.append(traj_atoms)

            # --- End of step loop for this restart ---

            # Get final structure for this restart
            with torch.no_grad():
                final_positions = torch.where(
                    fixed_mask.unsqueeze(-1),
                    positions_initial,
                    positions_parameter,
                )
                final_loss_components = self._positions_to_loss(
                    prepared,
                    final_positions,
                    target_2body,
                    target_3body,
                    target_4body,
                    fixed_mask=fixed_mask,
                    species_probabilities=probs,
                    minimum_distance_scale=minimum_distance_scale,
                )
                final_total_loss = self._compute_inverse_design_loss(
                    final_loss_components["fingerprint_loss"],
                    final_loss_components["repulsion_loss"],
                    fingerprint_loss_weight=fingerprint_loss_weight,
                    repulsion_weight=effective_repulsion_weight,
                    use_augmented_lagrangian=use_augmented_lagrangian,
                    lambda_mult=lambda_mult,
                    mu=mu,
                    num_atoms=len(working_atoms)
                )
                final_loss = float(final_total_loss.item())

                # Store final species probabilities for this restart
                if optimize_species:
                    restart_final_probs = probs.detach().clone()

                # Create final atoms object for this restart
                restart_final_atoms = working_atoms.copy()

                # Apply final species if optimizing
                if optimize_species:
                    if species_optimization_mode == 'argmax':
                        species_idx = torch.argmax(restart_final_probs, dim=-1)
                    elif species_optimization_mode == 'sample':
                        dist = torch.distributions.Categorical(probs=restart_final_probs)
                        species_idx = dist.sample()
                    else:
                        species_idx = torch.argmax(restart_final_probs, dim=-1)
                    final_symbols = [self.species_list[int(idx)] for idx in species_idx.cpu().numpy()]
                    restart_final_atoms.set_chemical_symbols(final_symbols)

                # Set final positions (with wrapping if specified)
                final_positions_np = final_positions.detach().cpu().numpy()
                if wrap_positions_to_cell:
                    restart_final_atoms.set_positions(final_positions_np)
                    restart_final_atoms = wrap_atoms_to_unit_cell(restart_final_atoms)
                else:
                    restart_final_atoms.set_positions(final_positions_np)

            # Track best restart by final loss
            if final_loss < best_restart_loss:
                best_restart_loss = final_loss
                best_restart_atoms = restart_final_atoms
                best_restart_trajectory = restart_trajectory

        # --- End of restarts ---

        # Return results
        if return_trajectory:
            # Return trajectory from the best restart
            if best_restart_trajectory is not None:
                return best_restart_trajectory
            else:
                return []
        else:
            # Return final atoms from the best restart
            return best_restart_atoms

    def compute_per_atom_loss(
        self,
        atoms: Atoms,
        target_fingerprint: np.ndarray,
        fingerprint_loss_weight: float = 1.0,
    ) -> Dict[str, Union[float, np.ndarray]]:
        """
        Compute per-atom fingerprint loss for a given structure.
        """
        self.eval()

        # Prepare structure
        prepared = self.prepare_structure(atoms, include_targets=False)
        positions = _float_tensor(prepared.positions, self._device)

        # Project target
        target = self._project_fingerprint_targets(
            _float_tensor(target_fingerprint, self._device)
        )

        # Get fingerprints and vertex contributions
        with torch.no_grad():
            (
                (vertex_2body, vertex_3body, vertex_4body),
                (prediction_2body, prediction_3body, prediction_4body),
            ) = self._forward_prepared(
                prepared,
                positions_override=positions,
                return_vertices=True,
            )

        # Split target
        target_2body = target[:self.fingerprint_dim_2body]
        target_3body = target[self.fingerprint_dim_2body:self.fingerprint_dim_2body + self.fingerprint_dim_3body]
        target_4body = target[self.fingerprint_dim_2body + self.fingerprint_dim_3body:]

        # Per-vertex squared errors (mean over fingerprint dimensions)
        per_vertex_2body = ((vertex_2body - target_2body) ** 2).mean(dim=1)  # (n_atoms,)
        per_vertex_3body = ((vertex_3body - target_3body) ** 2).mean(dim=1)  # (n_angles,)
        per_vertex_4body = ((vertex_4body - target_4body) ** 2).mean(dim=1)  # (n_triplets,)

        n_atoms = len(atoms)

        # --- 2-body: already per-atom ---
        per_atom_2body = per_vertex_2body

        # --- 3-body: Find center atom for each angle ---
        # angle_index is (n_angles, 2) where each entry is a pair index
        # The center atom is the common atom between the two pairs
        angle_index = torch.as_tensor(prepared.topology.angle_index, device=self._device)
        pair_index = torch.as_tensor(prepared.topology.pair_index, device=self._device)  # (n_pairs, 2)

        # For each angle, find the center atom (the one shared by both pairs)
        angle_center = torch.zeros(len(angle_index), dtype=torch.long, device=self._device)
        for i, (pair1, pair2) in enumerate(angle_index):
            # Get the two atoms in each pair
            atoms1 = pair_index[pair1]  # [atom_a, atom_b]
            atoms2 = pair_index[pair2]  # [atom_c, atom_d]
            # Find the common atom
            common = torch.where(atoms1.unsqueeze(1) == atoms2.unsqueeze(0))[0]
            angle_center[i] = atoms1[common[0]] if len(common) > 0 else atoms1[0]

        # --- 4-body: triplet_center_index directly gives center atom ---
        triplet_center = torch.as_tensor(prepared.topology.triplet_center_index, device=self._device)

        # Debug info
        print(f"[DEBUG] n_atoms: {n_atoms}")
        print(f"[DEBUG] per_vertex_2body shape: {per_vertex_2body.shape}")
        print(f"[DEBUG] per_vertex_3body shape: {per_vertex_3body.shape}")
        print(f"[DEBUG] per_vertex_4body shape: {per_vertex_4body.shape}")
        print(f"[DEBUG] angle_center shape: {angle_center.shape}, min/max: {angle_center.min()}/{angle_center.max()}")
        print(f"[DEBUG] triplet_center shape: {triplet_center.shape}, min/max: {triplet_center.min()}/{triplet_center.max()}")

        # Ensure indices are valid
        assert angle_center.max() < n_atoms, f"angle_center has index {angle_center.max()} >= {n_atoms}"
        assert triplet_center.max() < n_atoms, f"triplet_center has index {triplet_center.max()} >= {n_atoms}"

        # Accumulate errors per atom
        per_atom_3body = torch.zeros(n_atoms, device=self._device)
        per_atom_4body = torch.zeros(n_atoms, device=self._device)
        count_3body = torch.zeros(n_atoms, device=self._device)
        count_4body = torch.zeros(n_atoms, device=self._device)

        # Accumulate 3-body errors
        if len(angle_center) > 0:
            per_atom_3body.scatter_add_(0, angle_center, per_vertex_3body)
            count_3body.scatter_add_(0, angle_center, torch.ones_like(per_vertex_3body))

        # Accumulate 4-body errors
        if len(triplet_center) > 0:
            per_atom_4body.scatter_add_(0, triplet_center, per_vertex_4body)
            count_4body.scatter_add_(0, triplet_center, torch.ones_like(per_vertex_4body))

        # Average (avoid division by zero)
        per_atom_3body = per_atom_3body / count_3body.clamp_min(1.0)
        per_atom_4body = per_atom_4body / count_4body.clamp_min(1.0)

        # Total per-atom loss
        per_atom_loss = fingerprint_loss_weight * (per_atom_2body + per_atom_3body + per_atom_4body)

        return {
            'total_loss': float(per_atom_loss.sum().item()),
            'per_atom_loss': per_atom_loss.cpu().numpy(),
            'per_atom_2body': per_atom_2body.cpu().numpy(),
            'per_atom_3body': per_atom_3body.cpu().numpy(),
            'per_atom_4body': per_atom_4body.cpu().numpy(),
            'n_angles': len(per_vertex_3body),
            'n_triplets': len(per_vertex_4body),
        }

    def _compute_inverse_design_loss(
        self,
        fingerprint_loss: torch.Tensor,
        repulsion_loss: torch.Tensor,
        fingerprint_loss_weight: float,
        repulsion_weight: float,
        use_augmented_lagrangian: bool,
        lambda_mult: float,
        mu: float,
        num_atoms: int,  # NEW: number of atoms in the system
        # New parameters for dynamic repulsion scaling
        repulsion_effect_threshold: float = 0.1,  # Repulsion starts having effect here
        repulsion_only_threshold: float = 10.0,  # Only repulsion matters here
        fingerprint_suppression_strength: float = 1.0,  # How much to suppress fingerprint (0-1)
    ) -> torch.Tensor:
        """
        Combine loss components with dynamic scaling based on repulsion magnitude.

        The repulsion loss is normalized by the number of atoms since it's a sum
        rather than a mean. The scaling thresholds are adjusted accordingly.
        """

        # Normalize repulsion by number of atoms to make thresholds atom-count independent
        # This converts sum-based repulsion to a per-atom average

        with torch.no_grad():
            norm_repulsion = repulsion_loss / num_atoms
            rep_val = norm_repulsion.detach().item()
            # print(f"[DEBUG] Normalized repulsion value: {rep_val:.6e}")

            # --- Dynamic fingerprint suppression ---
            # Suppress fingerprint loss when repulsion is high
            if rep_val <= repulsion_effect_threshold:
                # No suppression when repulsion is negligible
                fingerprint_scale = 1.0
            elif rep_val >= repulsion_only_threshold:
                # Only repulsion matters (completely suppress fingerprint)
                fingerprint_scale = 0.0
            else:
                # Sigmoid-like decay between effect_threshold and only_threshold
                # Maps rep_val from [effect_threshold, only_threshold] to [1.0, 0.0]
                t = (rep_val - repulsion_effect_threshold) / (repulsion_only_threshold - repulsion_effect_threshold)
                # Use smooth sigmoid: scale = 1.0 / (1.0 + exp(5.0 * (t - 0.5)))
                # This gives S-curve: stays near 1.0 until ~0.3, drops sharply, stays near 0.0 after ~0.7
                sigmoid_t = 1.0 / (1.0 + torch.exp(torch.tensor(5.0 * (t - 0.5))))
                fingerprint_scale = 1.0 - fingerprint_suppression_strength * (1.0 - sigmoid_t)
                fingerprint_scale = float(fingerprint_scale)

        # Apply scaling to fingerprint loss
        scaled_fingerprint_weight = fingerprint_loss_weight * fingerprint_scale

        # --- Augmented Lagrangian with normalized repulsion ---
        if use_augmented_lagrangian:
            # Use normalized repulsion for the augmented Lagrangian terms
            # This ensures consistent behavior regardless of system size

            return (
                scaled_fingerprint_weight * fingerprint_loss
                + lambda_mult * repulsion_loss
                + 0.5 * mu * (repulsion_loss ** 2)
            )
        else:
            # Fixed weighting (also normalized)
            return (
                scaled_fingerprint_weight * fingerprint_loss
                + repulsion_weight * repulsion_loss
            )
