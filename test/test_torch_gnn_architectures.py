import unittest

import numpy as np
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


if __name__ == "__main__":
    unittest.main()
