import unittest
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from ase.build import bulk
from ase.io import read

from raffle import TorchGNNFingerprint, symmetry_aware_rmsd


FINGERPRINT_LOSS_WEIGHT = 0.25
TARGET_VERTEX_WEIGHT = 0.5
TARGET_POSITION_WEIGHT = 0.25
INVERSE_LR_DECAY_RATE = 0.0
MAX_FINAL_INVERSE_MSE = 5.0e-4
MAX_FINAL_RMSD = 1.0e-2
MIN_RMSD_REDUCTION_FACTOR = 0.25
MIN_POSITION_SHIFT = 1.0e-2


PERTURBATION = np.array(
    [
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.172792, 0.0410809, 0.165219],
        [-0.00651579, 0.00452678, 0.00223187],
        [-0.00268477, 0.00290559, 0.00182286],
        [0.00147066, 0.00014211, 0.00273356],
    ],
    dtype=np.float32,
)
@unittest.skipUnless(TorchGNNFingerprint is not None, "PyTorch multigraph fingerprint is unavailable")
class TestTorchGNNFingerprintCarbonWorkflow(unittest.TestCase):

    def test_carbon_training_gradient_and_inverse_design(self):
        repo_root = Path(__file__).resolve().parents[1]
        carbon_xyz = repo_root / "example" / "data" / "carbon.xyz"
        self.assertTrue(carbon_xyz.exists())

        carbon_structures = read(str(carbon_xyz), index=":")[:48]
        self.assertGreaterEqual(len(carbon_structures), 1)

        original = bulk("C", "diamond", a=3.567, cubic=True)
        original.pbc = True

        rng = np.random.default_rng(42)
        augmented_structures = []
        for _ in range(16):
            atoms = original.copy()
            atoms.set_positions(
                atoms.get_positions()
                + rng.normal(scale=0.04, size=atoms.positions.shape).astype(np.float32)
            )
            augmented_structures.append(atoms)

        perturbed = original.copy()
        perturbed.set_positions(perturbed.get_positions() + PERTURBATION)

        fixed_atoms = np.zeros(len(perturbed), dtype=bool)
        fixed_atoms[:4] = True

        model = TorchGNNFingerprint(
            species_list=["C"],
            hidden_dim=80,
            num_message_layers=2,
            learning_rate=5.0e-4,
            lr_decay_rate=5.0e-3,
            smooth_cutoff_width=0.2,
            seed=42,
        )

        analytical_fingerprint = model.compute_reference_fingerprint(perturbed)
        initial_prediction = model.predict(perturbed)
        initial_eval_mse = float(np.mean((initial_prediction - analytical_fingerprint) ** 2))

        history = model.fit(
            carbon_structures,
            num_epochs=30,
            batch_size=8,
            augment_structures=augmented_structures,
            verbose=0,
        )
        self.assertGreaterEqual(len(history), 2)
        self.assertTrue(np.all(np.isfinite(history)))

        trained_prediction = model.predict(perturbed)
        trained_eval_mse = float(np.mean((trained_prediction - analytical_fingerprint) ** 2))
        self.assertLess(trained_eval_mse, initial_eval_mse)

        vertex_2body, vertex_3body, vertex_4body = model.compute_vertex_fingerprints(perturbed)
        self.assertEqual(vertex_2body.shape[1], model.fingerprint_dim_2body)
        self.assertEqual(vertex_3body.shape[1], model.fingerprint_dim_3body)
        self.assertEqual(vertex_4body.shape[1], model.fingerprint_dim_4body)

        grad2, grad3, grad4 = model.compute_gradients(perturbed)
        self.assertEqual(grad2.shape, (len(perturbed), 3, model.fingerprint_dim_2body))
        self.assertEqual(grad3.shape, (len(perturbed), 3, model.fingerprint_dim_3body))
        self.assertEqual(grad4.shape, (len(perturbed), 3, model.fingerprint_dim_4body))
        self.assertTrue(np.all(np.isfinite(grad2)))
        self.assertTrue(np.all(np.isfinite(grad3)))
        self.assertTrue(np.all(np.isfinite(grad4)))

        target_fingerprint = model.compute_reference_fingerprint(original)
        initial_inverse_mse = float(np.mean((trained_prediction - target_fingerprint) ** 2))
        initial_rmsd = symmetry_aware_rmsd(original, perturbed)

        optimised = model.inverse_design(
            target_fingerprint=target_fingerprint,
            atoms=perturbed,
            fixed_atoms=fixed_atoms,
            num_steps=400,
            step_size=5.0e-3,
            verbose=0,
            target_atoms=original,
            fingerprint_loss_weight=FINGERPRINT_LOSS_WEIGHT,
            target_vertex_weight=TARGET_VERTEX_WEIGHT,
            target_position_weight=TARGET_POSITION_WEIGHT,
            inverse_lr_decay_rate=INVERSE_LR_DECAY_RATE,
        )
        final_prediction = model.predict(optimised)
        final_inverse_mse = float(np.mean((final_prediction - target_fingerprint) ** 2))
        self.assertLessEqual(final_inverse_mse, MAX_FINAL_INVERSE_MSE)

        np.testing.assert_allclose(
            optimised.get_positions()[:4],
            perturbed.get_positions()[:4],
            atol=1.0e-7,
        )

        rmsd = symmetry_aware_rmsd(original, optimised)
        position_shift = float(np.linalg.norm(optimised.get_positions() - perturbed.get_positions()))
        self.assertGreaterEqual(position_shift, MIN_POSITION_SHIFT)
        self.assertLessEqual(rmsd, initial_rmsd + 1.0e-7)
        self.assertLessEqual(rmsd, initial_rmsd * MIN_RMSD_REDUCTION_FACTOR)
        self.assertLessEqual(rmsd, MAX_FINAL_RMSD)

        output_dir = repo_root / "build"
        output_dir.mkdir(exist_ok=True)
        figure_path = output_dir / "torch_gnn_carbon_workflow.png"

        figure = plt.figure(figsize=(12, 4))
        ax1 = figure.add_subplot(1, 3, 1)
        ax1.plot(range(len(history)), history, marker="o")
        ax1.set_xlabel("Epoch")
        ax1.set_ylabel("Joint weighted MSE")
        ax1.set_title("Training")

        ax2 = figure.add_subplot(1, 3, 2)
        ax2.plot([0, 1], [initial_inverse_mse, final_inverse_mse], marker="o")
        ax2.set_xticks([0, 1], ["initial", "final"])
        ax2.set_ylabel("Fingerprint MSE")
        ax2.set_title("Inverse Design")

        ax3 = figure.add_subplot(1, 3, 3, projection="3d")
        original_positions = original.get_positions()
        optimised_positions = optimised.get_positions()
        ax3.scatter(
            original_positions[:, 0],
            original_positions[:, 1],
            original_positions[:, 2],
            s=35,
            label="original",
        )
        ax3.scatter(
            optimised_positions[:, 0],
            optimised_positions[:, 1],
            optimised_positions[:, 2],
            s=35,
            marker="^",
            label="optimised",
        )
        ax3.set_title("Structure")
        ax3.legend(loc="best")

        figure.tight_layout()
        figure.savefig(figure_path, dpi=150)
        plt.close(figure)

        self.assertTrue(figure_path.exists())


if __name__ == "__main__":
    unittest.main()
