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
    """Graph neural network for learning and predicting RAFFLE descriptor
    fingerprints from atomic structures represented as molecular graphs.

    This class wraps a simple graph neural network that maps atomic structure
    graphs (atoms = vertices, bonds = edges) to RAFFLE distribution function
    descriptors via message passing.

    Parameters
    ----------
    species_list : list of str
        List of chemical species (e.g. ['C', 'Si']).
    bond_cutoff : float, optional
        Bond cutoff distance in Angstroms. Default: 6.0.
    gnn_hidden_sizes : list of int, optional
        GNN hidden layer sizes. Default: [32, 32].
    learning_rate : float, optional
        Learning rate for Adam optimiser. Default: 0.001.
    """

    def __init__(
        self,
        species_list: List[str],
        bond_cutoff: float = 6.0,
        gnn_hidden_sizes: Optional[List[int]] = None,
        learning_rate: float = 0.001,
    ):
        self.species_list = [s.strip() for s in species_list]
        self.num_species = len(self.species_list)
        self.bond_cutoff = bond_cutoff

        if gnn_hidden_sizes is None:
            gnn_hidden_sizes = [32, 32]

        # Features per vertex: 3 coords + one-hot species
        self.num_vertex_features = 3 + self.num_species

        self._gnn_hidden_sizes = gnn_hidden_sizes
        self._learning_rate = learning_rate
        self._fingerprint_dim = None
        self._network = None
        self._is_trained = False

    @property
    def fingerprint_dim(self) -> Optional[int]:
        """Dimension of the fingerprint vector."""
        return self._fingerprint_dim

    @property
    def is_trained(self) -> bool:
        """Whether the network has been trained."""
        return self._is_trained

    def atoms_to_graph(
        self, atoms
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Convert an ASE Atoms object to a molecular graph.

        Parameters
        ----------
        atoms : ase.Atoms
            Atomic structure.

        Returns
        -------
        vertex_features : np.ndarray, shape (num_atoms, num_vertex_features)
            Node features [x, y, z, one_hot_species].
        adjacency : np.ndarray, shape (num_atoms, num_atoms)
            Adjacency matrix with self-loops (1 for connected, 0 otherwise).
        """
        n_atoms = len(atoms)
        vertex_features = np.zeros(
            (n_atoms, self.num_vertex_features), dtype=np.float32
        )

        # Get positions
        if any(atoms.pbc):
            positions = atoms.get_scaled_positions()
        else:
            positions = atoms.get_positions()

        symbols = atoms.get_chemical_symbols()

        for i, (sym, pos) in enumerate(zip(symbols, positions)):
            vertex_features[i, :3] = pos[:3]
            sym_stripped = sym.strip()
            if sym_stripped in self.species_list:
                sp_idx = self.species_list.index(sym_stripped)
                vertex_features[i, 3 + sp_idx] = 1.0

        # Build adjacency matrix
        cart_pos = atoms.get_positions()
        adjacency = np.eye(n_atoms, dtype=np.float32)  # self-loops

        for i in range(n_atoms):
            for j in range(i + 1, n_atoms):
                # Minimum image convention for periodic systems
                diff = cart_pos[j] - cart_pos[i]
                if any(atoms.pbc):
                    cell = atoms.get_cell()
                    scaled = np.linalg.solve(cell.T, diff)
                    scaled -= np.round(scaled)
                    diff = cell.T @ scaled

                dist = np.linalg.norm(diff)
                if 0 < dist <= self.bond_cutoff:
                    adjacency[i, j] = 1.0
                    adjacency[j, i] = 1.0

        return vertex_features, adjacency

    def compute_fingerprint(self, atoms) -> np.ndarray:
        """Compute the RAFFLE descriptor fingerprint for a structure.

        Uses the RAFFLE Fortran library to compute 2-body, 3-body, and
        4-body distribution functions and flattens them into a vector.

        Parameters
        ----------
        atoms : ase.Atoms
            Atomic structure.

        Returns
        -------
        np.ndarray
            Flattened fingerprint vector.
        """
        from raffle.generator import raffle_generator

        gen = raffle_generator()
        gen.distributions.set_element_energies(
            {s: 0.0 for s in self.species_list}
        )

        container = gen.distributions
        container.create(
            atoms,
            energy=(
                atoms.info.get("energy", 0.0) if hasattr(atoms, "info") else 0.0
            ),
        )

        df_2body = container.get_2body()
        df_3body = container.get_3body()
        df_4body = container.get_4body()

        fp_parts = []
        if df_2body is not None:
            fp_parts.append(df_2body.flatten())
        if df_3body is not None:
            fp_parts.append(df_3body.flatten())
        if df_4body is not None:
            fp_parts.append(df_4body.flatten())

        if not fp_parts:
            raise RuntimeError("Failed to compute descriptor fingerprint")

        fingerprint = np.concatenate(fp_parts).astype(np.float32)

        if self._fingerprint_dim is None:
            self._fingerprint_dim = len(fingerprint)
            self._network = SimpleGNN(
                num_vertex_features=self.num_vertex_features,
                gnn_hidden_sizes=self._gnn_hidden_sizes,
                output_dim=self._fingerprint_dim,
                learning_rate=self._learning_rate,
            )

        return fingerprint

    def compute_fingerprint_direct(self, atoms) -> np.ndarray:
        """Compute a simplified radial distribution fingerprint using numpy.

        Parameters
        ----------
        atoms : ase.Atoms
            Atomic structure.

        Returns
        -------
        np.ndarray
            Simplified radial distribution fingerprint.
        """
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

        if self._fingerprint_dim is None:
            self._fingerprint_dim = len(rdf)
            self._network = SimpleGNN(
                num_vertex_features=self.num_vertex_features,
                gnn_hidden_sizes=self._gnn_hidden_sizes,
                output_dim=self._fingerprint_dim,
                learning_rate=self._learning_rate,
            )

        return rdf

    def train(
        self,
        structures: list,
        num_epochs: int = 100,
        verbose: int = 0,
        use_simple_fingerprint: bool = False,
    ) -> List[float]:
        """Train the GNN on a set of atomic structures.

        Parameters
        ----------
        structures : list of ase.Atoms
            Training structures.
        num_epochs : int
            Number of training epochs.
        verbose : int
            Verbosity level.
        use_simple_fingerprint : bool
            If True, use simplified direct fingerprint computation.

        Returns
        -------
        list of float
            Training loss history.
        """
        compute_fp = (
            self.compute_fingerprint_direct if use_simple_fingerprint
            else self.compute_fingerprint
        )

        graphs = []
        targets = []
        for atoms in structures:
            vf, adj = self.atoms_to_graph(atoms)
            fp = compute_fp(atoms)
            graphs.append((vf, adj))
            targets.append(fp)

        if self._network is None:
            raise RuntimeError(
                "Network not initialised. "
                "Call compute_fingerprint first."
            )

        loss_history = []
        for epoch in range(num_epochs):
            epoch_loss = self._network.train_batch(graphs, targets)
            loss_history.append(epoch_loss)

            if verbose > 0 and (epoch + 1) % max(1, num_epochs // 10) == 0:
                print(f"  epoch={epoch+1}, loss={epoch_loss:.6f}")

        self._is_trained = True
        return loss_history

    def predict(self, atoms) -> np.ndarray:
        """Forward inference: predict fingerprint from atomic structure.

        Parameters
        ----------
        atoms : ase.Atoms
            Atomic structure.

        Returns
        -------
        np.ndarray
            Predicted fingerprint vector.
        """
        if self._network is None:
            raise RuntimeError("Network not initialised. Train first.")

        vf, adj = self.atoms_to_graph(atoms)
        return self._network.forward(vf, adj)

    def inverse_design(
        self,
        target_fingerprint: np.ndarray,
        atoms,
        fixed_atoms: np.ndarray,
        num_steps: int = 200,
        step_size: float = 0.01,
        verbose: int = 0,
        use_simple_fingerprint: bool = False,
    ):
        """Inverse design: optimise atomic positions to match target descriptor.

        Parameters
        ----------
        target_fingerprint : np.ndarray
            Target descriptor fingerprint.
        atoms : ase.Atoms
            Initial atomic structure (modified in place).
        fixed_atoms : np.ndarray of bool
            Mask array. True = atom is fixed, False = atom is optimisable.
        num_steps : int
            Number of optimisation steps.
        step_size : float
            Step size for coordinate perturbation.
        verbose : int
            Verbosity level.
        use_simple_fingerprint : bool
            If True, use simplified fingerprint computation.

        Returns
        -------
        ase.Atoms
            Optimised structure.
        """
        compute_fp = (
            self.compute_fingerprint_direct if use_simple_fingerprint
            else self.compute_fingerprint
        )

        delta = 1e-4
        best_loss = float("inf")
        best_positions = atoms.get_positions().copy()

        current_fp = compute_fp(atoms)
        loss = float(np.mean((current_fp - target_fingerprint) ** 2))
        best_loss = loss
        best_positions = atoms.get_positions().copy()

        if verbose > 0:
            print(f"  Inverse design step 0, loss = {loss:.6e}")

        for step in range(1, num_steps + 1):
            positions = atoms.get_positions().copy()
            grad_positions = np.zeros_like(positions)

            for i in range(len(atoms)):
                if i < len(fixed_atoms) and fixed_atoms[i]:
                    continue

                for coord in range(3):
                    perturbed_pos = positions.copy()
                    perturbed_pos[i, coord] += delta
                    atoms.set_positions(perturbed_pos)

                    perturbed_fp = compute_fp(atoms)
                    loss_perturbed = float(
                        np.mean((perturbed_fp - target_fingerprint) ** 2)
                    )
                    grad_positions[i, coord] = (loss_perturbed - loss) / delta

            new_positions = positions - step_size * grad_positions
            atoms.set_positions(new_positions)

            current_fp = compute_fp(atoms)
            loss = float(np.mean((current_fp - target_fingerprint) ** 2))

            if loss < best_loss:
                best_loss = loss
                best_positions = atoms.get_positions().copy()

            if verbose > 0 and step % 10 == 0:
                print(f"  Inverse design step {step}, loss = {loss:.6e}")

            if loss < 1e-8:
                if verbose > 0:
                    print(f"  Converged at step {step}")
                break

        atoms.set_positions(best_positions)

        if verbose > 0:
            print(f"  Inverse design final loss = {best_loss:.6e}")

        return atoms
