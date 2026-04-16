"""
Graph neural network fingerprint module for RAFFLE.

This module provides a GNN-based framework that:
  1. Converts atomic structures to molecular graphs (atoms→vertices, bonds→edges)
  2. Learns RAFFLE descriptor fingerprints via message-passing
  3. Supports forward inference for predicting descriptors
  4. Enables inverse design (generating structures from target descriptors)
  5. Allows partial atomic optimisation via boolean atom masks

The GNN is implemented in pure numpy to remain within the RAFFLE ecosystem
without external ML frameworks. The graph representation naturally respects
the topology of atomic structures.
"""
import numpy as np
from typing import Optional, List, Tuple

from . import _raffle
from .raffle import geom_rw


class SimpleGraphConv:
    """A single graph convolution layer using numpy.

    Implements a simplified message-passing update:
        h' = sigma(A_hat @ h @ W + b)
    where A_hat = D^{-1/2} A D^{-1/2} is the symmetrically normalised
    adjacency matrix (with self-loops).
    """

    def __init__(
        self,
        in_features: int,
        out_features: int,
        activation: str = "relu",
        rng: Optional[np.random.Generator] = None,
    ):
        if rng is None:
            rng = np.random.default_rng(42)

        # Glorot initialisation
        scale = np.sqrt(2.0 / (in_features + out_features))
        self.weight = rng.normal(0, scale, (in_features, out_features)).astype(
            np.float32
        )
        self.bias = np.zeros((1, out_features), dtype=np.float32)
        self.activation = activation

        # Adam state
        self._m_w = np.zeros_like(self.weight)
        self._v_w = np.zeros_like(self.weight)
        self._m_b = np.zeros_like(self.bias)
        self._v_b = np.zeros_like(self.bias)

    @staticmethod
    def _relu(x):
        return np.maximum(0, x)

    @staticmethod
    def _relu_grad(x):
        return (x > 0).astype(np.float32)

    def forward(self, h: np.ndarray, a_hat: np.ndarray) -> np.ndarray:
        """Forward pass.

        Parameters
        ----------
        h : np.ndarray, shape (num_vertices, in_features)
            Node feature matrix.
        a_hat : np.ndarray, shape (num_vertices, num_vertices)
            Normalised adjacency matrix.

        Returns
        -------
        np.ndarray, shape (num_vertices, out_features)
        """
        self._h_in = h
        self._a_hat = a_hat
        self._z = a_hat @ h @ self.weight + self.bias
        if self.activation == "relu":
            self._out = self._relu(self._z)
        else:
            self._out = self._z
        return self._out

    def backward(self, grad_out: np.ndarray) -> np.ndarray:
        """Backward pass. Returns gradient w.r.t. input h.

        Parameters
        ----------
        grad_out : np.ndarray, shape (num_vertices, out_features)
            Gradient from upstream.

        Returns
        -------
        np.ndarray, shape (num_vertices, in_features)
            Gradient w.r.t. input h.
        """
        if self.activation == "relu":
            grad_z = grad_out * self._relu_grad(self._z)
        else:
            grad_z = grad_out

        ah = self._a_hat @ self._h_in  # (N, in_features)
        self._grad_w = ah.T @ grad_z  # (in_features, out_features)
        self._grad_b = np.sum(grad_z, axis=0, keepdims=True)

        # Gradient w.r.t. h
        grad_h = self._a_hat.T @ (grad_z @ self.weight.T)
        return grad_h


class SimpleGNN:
    """A minimal graph neural network using numpy.

    Architecture:
      - Graph convolution layers for message passing
      - Global mean readout pooling
      - Dense output layer (linear activation)

    Supports forward pass, backpropagation, and Adam optimisation.
    """

    def __init__(
        self,
        num_vertex_features: int,
        gnn_hidden_sizes: List[int],
        output_dim: int,
        learning_rate: float = 0.001,
    ):
        """Initialise the GNN.

        Parameters
        ----------
        num_vertex_features : int
            Number of input features per vertex.
        gnn_hidden_sizes : list of int
            Sizes of graph convolution hidden layers.
        output_dim : int
            Dimension of the output vector (fingerprint).
        learning_rate : float
            Learning rate for Adam optimiser.
        """
        self.learning_rate = learning_rate
        self.output_dim = output_dim

        rng = np.random.default_rng(42)

        # Graph convolution layers
        self.conv_layers: List[SimpleGraphConv] = []
        in_dim = num_vertex_features
        for h_dim in gnn_hidden_sizes:
            self.conv_layers.append(
                SimpleGraphConv(in_dim, h_dim, activation="relu", rng=rng)
            )
            in_dim = h_dim

        # Dense output layer (from pooled graph-level vector to fingerprint)
        readout_dim = gnn_hidden_sizes[-1] if gnn_hidden_sizes else num_vertex_features
        scale = np.sqrt(2.0 / (readout_dim + output_dim))
        self.out_weight = rng.normal(0, scale, (readout_dim, output_dim)).astype(
            np.float32
        )
        self.out_bias = np.zeros((1, output_dim), dtype=np.float32)

        # Adam state for output layer
        self._m_ow = np.zeros_like(self.out_weight)
        self._v_ow = np.zeros_like(self.out_weight)
        self._m_ob = np.zeros_like(self.out_bias)
        self._v_ob = np.zeros_like(self.out_bias)

        self._t = 0

    def forward(
        self, vertex_features: np.ndarray, adjacency: np.ndarray
    ) -> np.ndarray:
        """Forward pass through the GNN.

        Parameters
        ----------
        vertex_features : np.ndarray, shape (num_vertices, num_features)
            Node feature matrix.
        adjacency : np.ndarray, shape (num_vertices, num_vertices)
            Adjacency matrix (with or without self-loops).

        Returns
        -------
        np.ndarray, shape (output_dim,)
            Graph-level output vector.
        """
        # Compute normalised adjacency: A_hat = D^{-1/2} A D^{-1/2}
        a = adjacency.copy()
        d = np.sum(a, axis=1)
        d_inv_sqrt = np.where(d > 0, 1.0 / np.sqrt(d), 0.0)
        d_mat = np.diag(d_inv_sqrt)
        a_hat = d_mat @ a @ d_mat
        self._a_hat = a_hat.astype(np.float32)

        # Message passing
        h = vertex_features
        for layer in self.conv_layers:
            h = layer.forward(h, self._a_hat)

        # Global mean pooling
        self._h_pool_input = h
        self._num_vertices = h.shape[0]
        graph_vec = np.mean(h, axis=0, keepdims=True)  # (1, readout_dim)
        self._graph_vec = graph_vec

        # Dense output (linear activation)
        output = graph_vec @ self.out_weight + self.out_bias  # (1, output_dim)
        return output[0]  # (output_dim,)

    def backward(self, y_true: np.ndarray) -> float:
        """Backward pass and Adam weight update.

        Parameters
        ----------
        y_true : np.ndarray, shape (output_dim,)
            Target output vector.

        Returns
        -------
        float
            Mean squared error loss.
        """
        y_pred = self._graph_vec @ self.out_weight + self.out_bias  # (1, out)
        y_pred = y_pred[0]

        # MSE loss
        diff = y_pred - y_true
        loss = float(np.mean(diff ** 2))

        # Gradient of MSE w.r.t. output
        grad_out = 2.0 * diff / len(y_true)  # (output_dim,)

        # Gradient w.r.t. output layer
        grad_ow = self._graph_vec.T @ grad_out.reshape(1, -1)
        grad_ob = grad_out.reshape(1, -1)

        # Gradient w.r.t. graph_vec
        grad_graph_vec = (grad_out.reshape(1, -1) @ self.out_weight.T)[0]

        # Gradient through mean pooling
        grad_h = (
            np.ones((self._num_vertices, 1), dtype=np.float32)
            * grad_graph_vec.reshape(1, -1)
            / self._num_vertices
        )

        # Backprop through conv layers (reverse order)
        for layer in reversed(self.conv_layers):
            grad_h = layer.backward(grad_h)

        # Adam update
        self._t += 1
        beta1, beta2, eps = 0.9, 0.999, 1e-8

        # Update output layer
        self._m_ow = beta1 * self._m_ow + (1 - beta1) * grad_ow
        self._v_ow = beta2 * self._v_ow + (1 - beta2) * grad_ow ** 2
        m_hat = self._m_ow / (1 - beta1 ** self._t)
        v_hat = self._v_ow / (1 - beta2 ** self._t)
        self.out_weight -= self.learning_rate * m_hat / (np.sqrt(v_hat) + eps)

        self._m_ob = beta1 * self._m_ob + (1 - beta1) * grad_ob
        self._v_ob = beta2 * self._v_ob + (1 - beta2) * grad_ob ** 2
        m_hat = self._m_ob / (1 - beta1 ** self._t)
        v_hat = self._v_ob / (1 - beta2 ** self._t)
        self.out_bias -= self.learning_rate * m_hat / (np.sqrt(v_hat) + eps)

        # Update conv layers
        for layer in self.conv_layers:
            layer._m_w = beta1 * layer._m_w + (1 - beta1) * layer._grad_w
            layer._v_w = beta2 * layer._v_w + (1 - beta2) * layer._grad_w ** 2
            m_hat = layer._m_w / (1 - beta1 ** self._t)
            v_hat = layer._v_w / (1 - beta2 ** self._t)
            layer.weight -= self.learning_rate * m_hat / (
                np.sqrt(v_hat) + eps
            )

            layer._m_b = beta1 * layer._m_b + (1 - beta1) * layer._grad_b
            layer._v_b = beta2 * layer._v_b + (1 - beta2) * layer._grad_b ** 2
            m_hat = layer._m_b / (1 - beta1 ** self._t)
            v_hat = layer._v_b / (1 - beta2 ** self._t)
            layer.bias -= self.learning_rate * m_hat / (np.sqrt(v_hat) + eps)

        return loss

    def train_batch(
        self,
        graphs: List[Tuple[np.ndarray, np.ndarray]],
        targets: List[np.ndarray],
    ) -> float:
        """Train on a batch of graphs.

        Parameters
        ----------
        graphs : list of (vertex_features, adjacency) tuples
        targets : list of np.ndarray target fingerprints

        Returns
        -------
        float
            Average loss over the batch.
        """
        total_loss = 0.0
        for (vf, adj), target in zip(graphs, targets):
            self.forward(vf, adj)
            total_loss += self.backward(target)
        return total_loss / len(graphs)



class GNNFingerprint:
    """Fortran-backed GNN fingerprint interface.

    This adapter exposes the ATHENA-based Fortran implementation in
    ``raffle__gnn_fingerprint`` through the Python package API used by the
    examples. It keeps the public constructor and core methods aligned with the
    previous Python implementation, but all training, prediction, and inverse
    design work is delegated to the compiled Fortran model.
    """

    def __init__(
        self,
        species_list: List[str],
        bond_cutoff: float = 6.0,
        gnn_hidden_sizes: Optional[List[int]] = None,
        learning_rate: float = 0.001,
        lr_decay_rate: float = 1.E-2,
        num_time_steps: int = 3,
        gnn_output_dim: int = 32,
        max_degree: int = 12,
        use_mlip_layer: bool = False,
        n_rbf: int = 20,
        kernel_hidden: int = 64,
        layer_type: int = -1,
        seed: int = 42,
    ):
        if gnn_hidden_sizes is None:
            gnn_hidden_sizes = [64]

        self.species_list = [str(species).strip()[:3].ljust(3) for species in species_list]
        self.num_species = len(self.species_list)
        self.bond_cutoff = float(bond_cutoff)
        self._gnn_hidden_sizes = [int(size) for size in gnn_hidden_sizes]
        self._learning_rate = float(learning_rate)
        self._lr_decay_rate = float(lr_decay_rate)
        self._num_time_steps = int(num_time_steps)
        self._gnn_output_dim = int(gnn_output_dim)
        self._max_degree = int(max_degree)
        self._use_mlip_layer = bool(use_mlip_layer)
        self._n_rbf = int(n_rbf)
        self._kernel_hidden = int(kernel_hidden)
        self._last_train_summary = None
        self._seed = seed

        # Resolve layer_type: -1 means infer from use_mlip_layer
        if layer_type >= 0:
            self._layer_type = int(layer_type)
        elif use_mlip_layer:
            self._layer_type = 1
        else:
            self._layer_type = 0

        self._handle = _raffle.f90wrap_gnn_fingerprint_type_initialise()
        hidden_sizes = np.asarray(self._gnn_hidden_sizes, dtype=np.int32)
        species_array = np.asarray(self.species_list, dtype="U3")
        _raffle.f90wrap_gnn_fingerprint_type__initialise(
            this=self._handle,
            species_list=species_array,
            n_species=len(self.species_list),
            num_time_steps=self._num_time_steps,
            gnn_output_dim=self._gnn_output_dim,
            max_degree=self._max_degree,
            hidden_sizes=hidden_sizes,
            n_hidden=hidden_sizes.size,
            learning_rate=self._learning_rate,
            lr_decay_rate=self._lr_decay_rate,
            bond_cutoff=self.bond_cutoff,
            use_mlip_layer=self._use_mlip_layer,
            n_rbf=self._n_rbf,
            kernel_hidden=self._kernel_hidden,
            layer_type_in=self._layer_type,
            seed=self._seed,
        )

    def __del__(self):
        handle = getattr(self, "_handle", None)
        if handle is None:
            return
        try:
            _raffle.f90wrap_gnn_fingerprint_type_finalise(handle)
        except Exception:
            pass
        self._handle = None

    @property
    def fingerprint_dim(self) -> int:
        return int(_raffle.f90wrap_gnn_fingerprint_type__get__fingerprint_dim(self._handle))

    @property
    def fingerprint_dim_2body(self) -> int:
        return int(
            _raffle.f90wrap_gnn_fingerprint_type__get__fingerprint_dim_2body(
                self._handle
            )
        )

    @property
    def fingerprint_dim_3body(self) -> int:
        return int(
            _raffle.f90wrap_gnn_fingerprint_type__get__fingerprint_dim_3body(
                self._handle
            )
        )

    @property
    def fingerprint_dim_4body(self) -> int:
        return int(
            _raffle.f90wrap_gnn_fingerprint_type__get__fingerprint_dim_4body(
                self._handle
            )
        )

    @property
    def component_dims(self) -> Tuple[int, int, int]:
        return (
            self.fingerprint_dim_2body,
            self.fingerprint_dim_3body,
            self.fingerprint_dim_4body,
        )

    @property
    def is_trained(self) -> bool:
        return bool(_raffle.f90wrap_gnn_fingerprint_type__get__is_trained(self._handle))

    @property
    def use_mlip_layer(self) -> bool:
        return bool(_raffle.f90wrap_gnn_fingerprint_type__get__use_mlip_layer(self._handle))

    @property
    def layer_type(self) -> int:
        return int(_raffle.f90wrap_gnn_fingerprint_type__get__layer_type(self._handle))

    LAYER_NAMES = {
        0: "duvenaud",
        1: "raffle_mlip",
        2: "schnet",
        3: "dimenet",
        4: "hybrid",
    }

    @property
    def layer_name(self) -> str:
        return self.LAYER_NAMES.get(self.layer_type, f"unknown({self.layer_type})")

    def _atoms_to_basis(self, atoms):
        basis = geom_rw.basis()
        basis.fromase(atoms)
        return basis

    def _structures_to_basis_handles(self, structures):
        basis_objects = [self._atoms_to_basis(atoms) for atoms in structures]
        basis_handles = np.asarray([basis._handle for basis in basis_objects], dtype=np.int32).T
        return basis_objects, basis_handles

    def _ensure_finite(self, label: str, values: np.ndarray) -> np.ndarray:
        array = np.asarray(values, dtype=np.float32)
        if not np.all(np.isfinite(array)):
            raise RuntimeError(f"Non-finite values detected in {label}.")
        return array

    def _dataset_mse(self, structures) -> float:
        losses = []
        for atoms in structures:
            target = self.compute_fingerprint(atoms)
            predicted = self.predict(atoms)
            losses.append(float(np.mean((predicted - target) ** 2)))
        return float(np.mean(losses)) if losses else 0.0

    def compute_fingerprint(self, atoms) -> np.ndarray:
        basis = self._atoms_to_basis(atoms)
        fingerprint = _raffle.f90wrap_gnn_fingerprint_type__compute_fingerprint(
            this=self._handle,
            basis=basis._handle,
            fp_dim=self.fingerprint_dim,
        )
        return self._ensure_finite("fingerprint", fingerprint)

    def compute_fingerprint_components(self, atoms) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        basis = self._atoms_to_basis(atoms)
        fp2, fp3, fp4 = _raffle.f90wrap_gnn_fingerprint_type__compute_fingerprint_components(
            this=self._handle,
            basis=basis._handle,
            fp_dim_2body=self.fingerprint_dim_2body,
            fp_dim_3body=self.fingerprint_dim_3body,
            fp_dim_4body=self.fingerprint_dim_4body,
        )
        return (
            self._ensure_finite("2-body fingerprint", fp2),
            self._ensure_finite("3-body fingerprint", fp3),
            self._ensure_finite("4-body fingerprint", fp4),
        )

    def compute_gradients(self, atoms) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        basis = self._atoms_to_basis(atoms)
        grad2, grad3, grad4 = _raffle.f90wrap_gnn_fingerprint_type__compute_gradients(
            this=self._handle,
            basis=basis._handle,
            fp_dim_2body=self.fingerprint_dim_2body,
            fp_dim_3body=self.fingerprint_dim_3body,
            fp_dim_4body=self.fingerprint_dim_4body,
            n_atoms=len(atoms),
        )
        grad2 = self._ensure_finite("2-body gradients", grad2).transpose(1, 2, 0)
        grad3 = self._ensure_finite("3-body gradients", grad3).transpose(1, 2, 0)
        grad4 = self._ensure_finite("4-body gradients", grad4).transpose(1, 2, 0)
        return grad2, grad3, grad4

    def compute_fingerprint_direct(self, atoms) -> np.ndarray:
        from ase.geometry import get_distances

        r_min, r_max = 0.5, 6.0
        n_bins = 220
        sigma = 0.1

        bin_edges = np.linspace(r_min, r_max, n_bins + 1)
        bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
        positions = atoms.get_positions()
        cell = atoms.get_cell() if any(atoms.pbc) else None
        _, distances = get_distances(positions, cell=cell, pbc=atoms.pbc)

        n_atoms = len(atoms)
        rdf = np.zeros(n_bins, dtype=np.float32)
        eta = 1.0 / (2.0 * sigma ** 2)

        for i in range(n_atoms):
            for j in range(n_atoms):
                if i == j:
                    continue
                d = distances[i, j]
                if r_min <= d <= r_max:
                    rdf += np.exp(-eta * (d - bin_centers) ** 2)

        if n_atoms > 1:
            rdf *= np.sqrt(eta / np.pi) / n_atoms

        return rdf

    def train(
        self,
        structures: list,
        num_epochs: int = 100,
        batch_size: int = None,
        verbose: int = 0,
        use_simple_fingerprint: bool = False,
    ) -> List[float]:
        if use_simple_fingerprint:
            raise NotImplementedError(
                "The Fortran-backed GNN trains against the RAFFLE descriptor fingerprint only."
            )

        if not structures:
            raise ValueError("At least one structure is required for training.")

        initial_loss = self._dataset_mse(structures)
        basis_objects, basis_handles = self._structures_to_basis_handles(structures)
        effective_batch_size = min(len(basis_objects), 32) if batch_size is None else max(1, min(int(batch_size), len(basis_objects)))
        _raffle.f90wrap_gnn_fingerprint_type__train(
            this=self._handle,
            basis_handles=basis_handles,
            num_epochs=int(num_epochs),
            batch_size=effective_batch_size,
            verbose=int(verbose),
            n_structures=len(basis_objects),
        )
        final_loss = self._dataset_mse(structures)
        if not np.isfinite(final_loss):
            raise RuntimeError("Training produced a non-finite loss.")
        self._last_train_summary = [initial_loss, final_loss]
        return [initial_loss, final_loss]

    def predict(self, atoms, use_simple_fingerprint: bool = False) -> np.ndarray:
        del use_simple_fingerprint
        basis = self._atoms_to_basis(atoms)
        prediction = _raffle.f90wrap_gnn_fingerprint_type__predict(
            this=self._handle,
            basis=basis._handle,
            fp_dim=self.fingerprint_dim,
        )
        return self._ensure_finite("prediction", prediction)

    def inverse_design(
        self,
        target_fingerprint: np.ndarray,
        atoms,
        fixed_atoms: np.ndarray,
        num_steps: int = 500,
        step_size: float = 1.0,
        verbose: int = 0,
        use_simple_fingerprint: bool = False,
        use_predict: bool = True,
    ):
        if use_simple_fingerprint:
            raise NotImplementedError(
                "The Fortran-backed GNN inverse design path uses the RAFFLE descriptor fingerprint only."
            )

        basis = self._atoms_to_basis(atoms)
        target = self._ensure_finite("inverse-design target fingerprint", target_fingerprint)
        fixed = np.asarray(fixed_atoms, dtype=bool)
        _raffle.f90wrap_gnn_fingerprint_type__inverse_design(
            this=self._handle,
            target_fp=target,
            basis=basis._handle,
            fixed_atoms=fixed,
            num_steps=int(num_steps),
            step_size=float(step_size),
            verbose=int(verbose),
            use_predict=bool(use_predict),
            fp_dim=target.size,
            n_atoms=fixed.size,
        )
        optimised = basis.toase()
        optimised.info.update(getattr(atoms, "info", {}))
        return optimised
