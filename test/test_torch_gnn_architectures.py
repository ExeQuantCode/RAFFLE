import unittest

import numpy as np
import torch
from ase.build import bulk

from raffle import TorchGNNFingerprint


@unittest.skipUnless(TorchGNNFingerprint is not None, "PyTorch multigraph fingerprint is unavailable")
class TestTorchGNNArchitectures(unittest.TestCase):

    def test_supported_architectures_produce_nonnegative_predictions(self):
        original = bulk("C", "diamond", a=3.567, cubic=True)
        original.pbc = True

        for architecture in (
            "residual",
            "gated",
            "attention",
            "attention_coupled",
            "attention_conservative",
            "graph_transformer",
            "graph_transformer_coupled",
            "graph_operator",
            "graph_operator_coupled",
            "multkan",
            "multkan_coupled",
        ):
            with self.subTest(architecture=architecture):
                model = TorchGNNFingerprint(
                    species_list=["C"],
                    architecture=architecture,
                    hidden_dim=32,
                    num_message_layers=2,
                    seed=42,
                )
                prediction = model.predict(original)
                component_predictions = model.predict_components(original)

                self.assertEqual(prediction.shape, (model.fingerprint_dim,))
                self.assertTrue(np.all(prediction >= 0.0))
                self.assertEqual(component_predictions[0].shape, (model.fingerprint_dim_2body,))
                self.assertEqual(component_predictions[1].shape, (model.fingerprint_dim_3body,))
                self.assertEqual(component_predictions[2].shape, (model.fingerprint_dim_4body,))
                self.assertTrue(np.all(component_predictions[0] >= 0.0))
                self.assertTrue(np.all(component_predictions[1] >= 0.0))
                self.assertTrue(np.all(component_predictions[2] >= 0.0))

    def test_component_loss_penalises_zero_sparse_predictions(self):
        original = bulk("C", "diamond", a=3.567, cubic=True)
        original.pbc = True
        model = TorchGNNFingerprint(
            species_list=["C"],
            architecture="residual",
            hidden_dim=32,
            num_message_layers=2,
            component_weight=(2.0, 0.0, 0.0),
            seed=42,
        )

        target_2body, target_3body, target_4body = model.compute_reference_components(original)
        target_2body_tensor = torch.as_tensor(target_2body, dtype=torch.float32, device=model._device)
        target_3body_tensor = torch.as_tensor(target_3body, dtype=torch.float32, device=model._device)
        target_4body_tensor = torch.as_tensor(target_4body, dtype=torch.float32, device=model._device)
        zero_2body = torch.zeros_like(target_2body_tensor)
        zero_3body = torch.zeros_like(target_3body_tensor)
        zero_4body = torch.zeros_like(target_4body_tensor)

        zero_loss = model._component_loss(
            zero_2body,
            zero_3body,
            zero_4body,
            target_2body_tensor,
            target_3body_tensor,
            target_4body_tensor,
        )
        matched_loss = model._component_loss(
            target_2body_tensor,
            zero_3body,
            zero_4body,
            target_2body_tensor,
            target_3body_tensor,
            target_4body_tensor,
        )

        self.assertAlmostEqual(float(zero_loss.item()), 2.0, places=5)
        self.assertAlmostEqual(float(matched_loss.item()), 0.0, places=6)

    def test_untrained_2body_component_matches_reference_over_lattice_sweep(self):
        model = TorchGNNFingerprint(
            species_list=["C"],
            architecture="residual",
            hidden_dim=32,
            num_message_layers=2,
            seed=42,
        )

        for lattice_constant in np.linspace(3.3, 3.9, 5):
            with self.subTest(lattice_constant=float(lattice_constant)):
                atoms = bulk("C", "diamond", a=float(lattice_constant), cubic=True)
                atoms.pbc = True
                predicted_2body, _, _ = model.predict_components(atoms)
                reference_2body, _, _ = model.compute_reference_components(atoms)

                self.assertTrue(
                    np.allclose(predicted_2body, reference_2body, atol=1.0e-6, rtol=1.0e-6)
                )


if __name__ == "__main__":
    unittest.main()
