import unittest
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from ase.build import bulk
from ase.io import read

from raffle.gnn_fingerprint import GNNFingerprint


class TestGNNFingerprintCarbonWorkflow(unittest.TestCase):

    def test_carbon_training_gradient_and_inverse_design(self):
        repo_root = Path(__file__).resolve().parents[1]
        carbon_xyz = repo_root / "example" / "data" / "carbon.xyz"
        self.assertTrue(carbon_xyz.exists())

        structures = read(str(carbon_xyz), index=":")
        self.assertGreaterEqual(len(structures), 1)
        structures = structures[: min(len(structures), 4)]

        model = GNNFingerprint(
            species_list=["C"],
            bond_cutoff=6.0,
            gnn_hidden_sizes=[32],
            learning_rate=5.0e-4,
            lr_decay_rate=1.0e-2,
            num_time_steps=2,
            gnn_output_dim=16,
            max_degree=8,
            layer_type=1,
            n_rbf=12,
            kernel_hidden=32,
            seed=42,
        )

        analytical_fingerprints = [model.compute_fingerprint(atoms) for atoms in structures]
        self.assertEqual(len(analytical_fingerprints), len(structures))
        self.assertTrue(all(np.all(np.isfinite(fp)) for fp in analytical_fingerprints))

        original = bulk("C", "diamond", a=3.567, cubic=True)
        original.pbc = True
        perturbed = original.copy()
        perturbation = np.array(
            [
                [0.020, -0.015, 0.000],
                [-0.010, 0.010, 0.015],
                [0.012, 0.000, -0.018],
                [-0.014, -0.010, 0.010],
                [0.018, 0.015, -0.010],
                [-0.020, 0.012, 0.000],
                [0.010, -0.018, 0.016],
                [0.000, 0.015, -0.012],
            ],
            dtype=np.float32,
        )
        perturbed.set_positions(perturbed.get_positions() + perturbation)

        analytical_eval = model.compute_fingerprint(perturbed)
        initial_prediction = model.predict(perturbed)
        initial_eval_mse = float(np.mean((initial_prediction - analytical_eval) ** 2))

        train_summary = model.train(structures, num_epochs=8, batch_size=2, verbose=0)
        self.assertEqual(len(train_summary), 2)
        self.assertTrue(np.all(np.isfinite(train_summary)))

        trained_prediction = model.predict(perturbed)
        trained_eval_mse = float(np.mean((trained_prediction - analytical_eval) ** 2))
        self.assertTrue(np.isfinite(trained_eval_mse))
        self.assertLess(trained_eval_mse, initial_eval_mse)

        fp2, fp3, fp4 = model.compute_fingerprint_components(perturbed)
        self.assertEqual(fp2.shape[0] + fp3.shape[0] + fp4.shape[0], model.fingerprint_dim)

        grad2, grad3, grad4 = model.compute_gradients(perturbed)
        self.assertEqual(grad2.shape, (len(perturbed), 3, model.fingerprint_dim_2body))
        self.assertEqual(grad3.shape, (len(perturbed), 3, model.fingerprint_dim_3body))
        self.assertEqual(grad4.shape, (len(perturbed), 3, model.fingerprint_dim_4body))
        self.assertTrue(np.all(np.isfinite(grad2)))
        self.assertTrue(np.all(np.isfinite(grad3)))
        self.assertTrue(np.all(np.isfinite(grad4)))

        fixed_atoms = np.zeros(len(perturbed), dtype=bool)
        fixed_atoms[:4] = True
        target_fingerprint = model.compute_fingerprint(original)
        initial_inverse_mse = float(np.mean((trained_prediction - target_fingerprint) ** 2))

        optimised = model.inverse_design(
            target_fingerprint=target_fingerprint,
            atoms=perturbed,
            fixed_atoms=fixed_atoms,
            num_steps=4,
            step_size=2.0e-3,
            verbose=0,
            use_predict=True,
        )
        final_prediction = model.predict(optimised)
        final_inverse_mse = float(np.mean((final_prediction - target_fingerprint) ** 2))

        self.assertTrue(np.isfinite(final_inverse_mse))
        self.assertLessEqual(final_inverse_mse, initial_inverse_mse)
        wrapped_reference = perturbed.copy()
        wrapped_reference.wrap()
        wrapped_optimised = optimised.copy()
        wrapped_optimised.wrap()
        np.testing.assert_allclose(
            wrapped_optimised.get_positions()[:4],
            wrapped_reference.get_positions()[:4],
            atol=1.0e-6,
        )

        output_dir = repo_root / "build"
        output_dir.mkdir(exist_ok=True)
        figure_path = output_dir / "gnn_carbon_workflow.png"

        fig = plt.figure(figsize=(12, 4))
        ax1 = fig.add_subplot(1, 3, 1)
        ax1.plot([0, 1], train_summary, marker="o")
        ax1.set_xticks([0, 1], ["initial", "final"])
        ax1.set_ylabel("Dataset MSE")
        ax1.set_title("Training")

        ax2 = fig.add_subplot(1, 3, 2)
        ax2.plot([0, 1], [initial_inverse_mse, final_inverse_mse], marker="o")
        ax2.set_xticks([0, 1], ["initial", "final"])
        ax2.set_ylabel("Fingerprint MSE")
        ax2.set_title("Inverse Design")

        ax3 = fig.add_subplot(1, 3, 3, projection="3d")
        original_positions = original.get_positions()
        optimised_positions = optimised.get_positions()
        ax3.scatter(
            original_positions[:, 0],
            original_positions[:, 1],
            original_positions[:, 2],
            label="original",
            s=35,
        )
        ax3.scatter(
            optimised_positions[:, 0],
            optimised_positions[:, 1],
            optimised_positions[:, 2],
            label="optimised",
            marker="^",
            s=35,
        )
        ax3.set_title("Structure")
        ax3.legend(loc="best")

        fig.tight_layout()
        fig.savefig(figure_path, dpi=150)
        plt.close(fig)

        self.assertTrue(figure_path.exists())


if __name__ == "__main__":
    unittest.main()
