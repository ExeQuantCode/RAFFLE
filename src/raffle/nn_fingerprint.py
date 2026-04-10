"""
Neural network fingerprint module for RAFFLE.

This module provides a neural network framework that:
  1. Learns RAFFLE descriptor fingerprints from atomic structures
  2. Supports forward inference for predicting descriptors
  3. Enables inverse design (generating structures from target descriptors)
  4. Allows partial atomic optimisation via boolean atom masks

The NN is implemented in pure numpy to remain within the RAFFLE ecosystem
without external ML frameworks. The descriptor computation uses the existing
RAFFLE Fortran library via f90wrap bindings.
"""
import numpy as np
from typing import Optional, List, Tuple


class SimpleNN:
    """A minimal feedforward neural network using numpy.

    Supports forward pass, backpropagation, and Adam optimisation.
    """

    def __init__(self, layer_sizes: List[int], learning_rate: float = 0.001):
        """Initialise the network.

        Parameters
        ----------
        layer_sizes : list of int
            Sizes of each layer including input and output.
            e.g. [input_dim, hidden1, hidden2, output_dim]
        learning_rate : float
            Learning rate for Adam optimiser.
        """
        self.layer_sizes = layer_sizes
        self.learning_rate = learning_rate
        self.weights = []
        self.biases = []

        # Xavier initialisation
        rng = np.random.default_rng(42)
        for i in range(len(layer_sizes) - 1):
            scale = np.sqrt(2.0 / (layer_sizes[i] + layer_sizes[i + 1]))
            w = rng.normal(0, scale, (layer_sizes[i], layer_sizes[i + 1]))
            b = np.zeros((1, layer_sizes[i + 1]))
            self.weights.append(w.astype(np.float32))
            self.biases.append(b.astype(np.float32))

        # Adam state
        self._m_w = [np.zeros_like(w) for w in self.weights]
        self._v_w = [np.zeros_like(w) for w in self.weights]
        self._m_b = [np.zeros_like(b) for b in self.biases]
        self._v_b = [np.zeros_like(b) for b in self.biases]
        self._t = 0

    @staticmethod
    def _relu(x):
        return np.maximum(0, x)

    @staticmethod
    def _relu_grad(x):
        return (x > 0).astype(np.float32)

    def forward(self, x: np.ndarray) -> np.ndarray:
        """Forward pass through the network.

        Parameters
        ----------
        x : np.ndarray
            Input array of shape (batch_size, input_dim).

        Returns
        -------
        np.ndarray
            Output array of shape (batch_size, output_dim).
        """
        self._activations = [x]
        self._pre_activations = []

        for i in range(len(self.weights)):
            z = self._activations[-1] @ self.weights[i] + self.biases[i]
            self._pre_activations.append(z)
            if i < len(self.weights) - 1:
                # ReLU for hidden layers
                a = self._relu(z)
            else:
                # Linear for output layer
                a = z
            self._activations.append(a)

        return self._activations[-1]

    def backward(self, y_true: np.ndarray) -> float:
        """Backward pass and weight update using Adam.

        Parameters
        ----------
        y_true : np.ndarray
            Target output of shape (batch_size, output_dim).

        Returns
        -------
        float
            Mean squared error loss.
        """
        batch_size = y_true.shape[0]
        y_pred = self._activations[-1]

        # MSE loss
        loss = np.mean((y_pred - y_true) ** 2)

        # Output layer gradient (linear activation, MSE loss)
        delta = 2.0 * (y_pred - y_true) / (batch_size * y_true.shape[1])

        grad_w = []
        grad_b = []

        for i in range(len(self.weights) - 1, -1, -1):
            gw = self._activations[i].T @ delta
            gb = np.sum(delta, axis=0, keepdims=True)
            grad_w.insert(0, gw)
            grad_b.insert(0, gb)

            if i > 0:
                delta = (delta @ self.weights[i].T) * \
                    self._relu_grad(self._pre_activations[i - 1])

        # Adam update
        self._t += 1
        beta1, beta2, eps = 0.9, 0.999, 1e-8

        for i in range(len(self.weights)):
            self._m_w[i] = beta1 * self._m_w[i] + (1 - beta1) * grad_w[i]
            self._v_w[i] = beta2 * self._v_w[i] + (1 - beta2) * grad_w[i]**2
            m_hat_w = self._m_w[i] / (1 - beta1**self._t)
            v_hat_w = self._v_w[i] / (1 - beta2**self._t)
            self.weights[i] -= self.learning_rate * m_hat_w / \
                (np.sqrt(v_hat_w) + eps)

            self._m_b[i] = beta1 * self._m_b[i] + (1 - beta1) * grad_b[i]
            self._v_b[i] = beta2 * self._v_b[i] + (1 - beta2) * grad_b[i]**2
            m_hat_b = self._m_b[i] / (1 - beta1**self._t)
            v_hat_b = self._v_b[i] / (1 - beta2**self._t)
            self.biases[i] -= self.learning_rate * m_hat_b / \
                (np.sqrt(v_hat_b) + eps)

        return float(loss)


class NNFingerprint:
    """Neural network for learning and predicting RAFFLE descriptor
    fingerprints from atomic structures, and for inverse design.

    This class wraps a simple feedforward neural network that maps
    atomic structure features to RAFFLE distribution function descriptors.

    Parameters
    ----------
    species_list : list of str
        List of chemical species (e.g. ['C', 'Si']).
    max_atoms : int
        Maximum number of atoms per structure (for input padding).
    hidden_layer_sizes : list of int, optional
        Sizes of hidden layers. Default: [128, 64].
    learning_rate : float, optional
        Learning rate for Adam optimiser. Default: 0.001.
    """

    def __init__(
        self,
        species_list: List[str],
        max_atoms: int,
        hidden_layer_sizes: Optional[List[int]] = None,
        learning_rate: float = 0.001,
    ):
        self.species_list = [s.strip() for s in species_list]
        self.num_species = len(self.species_list)
        self.max_atoms = max_atoms

        if hidden_layer_sizes is None:
            hidden_layer_sizes = [128, 64]

        # Features per atom: 3 coords + one-hot species
        self.features_per_atom = 3 + self.num_species
        self.input_dim = max_atoms * self.features_per_atom

        # Fingerprint dimension will be set on first compute
        self._fingerprint_dim = None
        self._network = None
        self._hidden_sizes = hidden_layer_sizes
        self._learning_rate = learning_rate
        self._is_trained = False

    @property
    def fingerprint_dim(self) -> Optional[int]:
        """Dimension of the fingerprint vector."""
        return self._fingerprint_dim

    @property
    def is_trained(self) -> bool:
        """Whether the network has been trained."""
        return self._is_trained

    def _atoms_to_input(self, atoms) -> np.ndarray:
        """Convert an ASE Atoms object to a flat input vector.

        Parameters
        ----------
        atoms : ase.Atoms
            Atomic structure.

        Returns
        -------
        np.ndarray
            Flat input vector of shape (input_dim,).
        """
        input_vec = np.zeros(self.input_dim, dtype=np.float32)

        # Get scaled (fractional) positions if periodic, else positions
        if any(atoms.pbc):
            positions = atoms.get_scaled_positions()
        else:
            positions = atoms.get_positions()

        symbols = atoms.get_chemical_symbols()

        for idx, (sym, pos) in enumerate(zip(symbols, positions)):
            if idx >= self.max_atoms:
                break

            offset = idx * self.features_per_atom
            input_vec[offset:offset + 3] = pos[:3]

            # One-hot species encoding
            sym_stripped = sym.strip()
            if sym_stripped in self.species_list:
                sp_idx = self.species_list.index(sym_stripped)
                input_vec[offset + 3 + sp_idx] = 1.0

        return input_vec

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

        df_2body, df_3body, df_4body = gen.distributions.generate_fingerprint(
            atoms
        )

        fingerprint = np.concatenate(
            [
                np.asarray(df_2body, dtype=np.float32).flatten(order="F"),
                np.asarray(df_3body, dtype=np.float32).flatten(order="F"),
                np.asarray(df_4body, dtype=np.float32).flatten(order="F"),
            ]
        )

        # Set fingerprint dim on first call
        if self._fingerprint_dim is None:
            self._fingerprint_dim = len(fingerprint)
            # Now initialise the network
            layer_sizes = (
                [self.input_dim]
                + self._hidden_sizes
                + [self._fingerprint_dim]
            )
            self._network = SimpleNN(layer_sizes, self._learning_rate)
        elif self._fingerprint_dim != len(fingerprint):
            raise RuntimeError(
                "Fingerprint dimension changed. Use a consistent fingerprint "
                "mode for a single NNFingerprint instance."
            )

        return fingerprint

    def compute_fingerprint_direct(self, atoms) -> np.ndarray:
        """Compute fingerprint directly from atomic structure using numpy.

        This is a simplified version that computes 2-body (radial)
        distribution functions without the full RAFFLE Fortran library,
        for cases where the Fortran bindings are not available.

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

        # Parameters
        r_min, r_max = 0.5, 6.0
        n_bins = 220
        sigma = 0.1

        bin_edges = np.linspace(r_min, r_max, n_bins + 1)
        bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

        # Get all pairwise distances
        positions = atoms.get_positions()
        cell = atoms.get_cell() if any(atoms.pbc) else None
        _, distances = get_distances(
            positions, cell=cell, pbc=atoms.pbc
        )

        # Build radial distribution
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
            layer_sizes = (
                [self.input_dim]
                + self._hidden_sizes
                + [self._fingerprint_dim]
            )
            self._network = SimpleNN(layer_sizes, self._learning_rate)

        return rdf

    def train(
        self,
        structures: list,
        num_epochs: int = 100,
        batch_size: int = 1,
        verbose: int = 0,
        use_simple_fingerprint: bool = False,
    ) -> List[float]:
        """Train the neural network on a set of atomic structures.

        Parameters
        ----------
        structures : list of ase.Atoms
            Training structures.
        num_epochs : int
            Number of training epochs.
        batch_size : int
            Batch size (currently trains on full batch).
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

        # Build training data
        inputs = []
        targets = []
        for atoms in structures:
            inp = self._atoms_to_input(atoms)
            fp = compute_fp(atoms)
            inputs.append(inp)
            targets.append(fp)

        input_data = np.array(inputs, dtype=np.float32)
        target_data = np.array(targets, dtype=np.float32)

        if self._network is None:
            raise RuntimeError(
                "Network not initialised. "
                "Call compute_fingerprint first or ensure fingerprint_dim "
                "is set."
            )

        # Training loop
        loss_history = []
        for epoch in range(num_epochs):
            output = self._network.forward(input_data)
            loss = self._network.backward(target_data)
            loss_history.append(loss)

            if verbose > 0 and (epoch + 1) % max(1, num_epochs // 10) == 0:
                print(f"  epoch={epoch+1}, loss={loss:.6f}")

        self._is_trained = True
        return loss_history

    def predict(self, atoms, use_simple_fingerprint: bool = False) -> np.ndarray:
        """Forward inference: predict fingerprint from atomic structure.

        Parameters
        ----------
        atoms : ase.Atoms
            Atomic structure.
        use_simple_fingerprint : bool
            If True, use simplified fingerprint (for initialisation check).

        Returns
        -------
        np.ndarray
            Predicted fingerprint vector.
        """
        if self._network is None:
            raise RuntimeError("Network not initialised. Train first.")

        input_vec = self._atoms_to_input(atoms)
        output = self._network.forward(input_vec.reshape(1, -1))
        return output[0]

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
        best_loss = float('inf')
        best_positions = atoms.get_positions().copy()

        # Compute initial loss
        current_fp = compute_fp(atoms)
        loss = np.mean((current_fp - target_fingerprint) ** 2)
        best_loss = loss
        best_positions = atoms.get_positions().copy()

        if verbose > 0:
            print(f"  Inverse design step 0, loss = {loss:.6e}")

        for step in range(1, num_steps + 1):
            positions = atoms.get_positions().copy()

            # Compute all gradients first
            grad_positions = np.zeros_like(positions)
            for i in range(len(atoms)):
                if i < len(fixed_atoms) and fixed_atoms[i]:
                    continue

                for coord in range(3):
                    # Perturb
                    perturbed_pos = positions.copy()
                    perturbed_pos[i, coord] += delta
                    atoms.set_positions(perturbed_pos)

                    # Compute perturbed fingerprint
                    perturbed_fp = compute_fp(atoms)
                    loss_perturbed = np.mean(
                        (perturbed_fp - target_fingerprint) ** 2
                    )

                    # Finite-difference gradient
                    grad_positions[i, coord] = (loss_perturbed - loss) / delta

            # Apply all gradient updates at once
            new_positions = positions - step_size * grad_positions
            atoms.set_positions(new_positions)

            # Recompute loss
            current_fp = compute_fp(atoms)
            loss = np.mean((current_fp - target_fingerprint) ** 2)

            if loss < best_loss:
                best_loss = loss
                best_positions = atoms.get_positions().copy()

            if verbose > 0 and step % 10 == 0:
                print(f"  Inverse design step {step}, loss = {loss:.6e}")

            if loss < 1e-8:
                if verbose > 0:
                    print(f"  Converged at step {step}")
                break

        # Restore best positions
        atoms.set_positions(best_positions)

        if verbose > 0:
            print(f"  Inverse design final loss = {best_loss:.6e}")

        return atoms
