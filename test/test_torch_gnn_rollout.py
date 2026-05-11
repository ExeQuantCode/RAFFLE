import sys
import unittest
from pathlib import Path

from ase.build import bulk


REPO_ROOT = Path(__file__).resolve().parents[1]
PYTHON_PKG_DIR = REPO_ROOT / "example" / "python_pkg"
if str(PYTHON_PKG_DIR) not in sys.path:
    sys.path.insert(0, str(PYTHON_PKG_DIR))

from torch_gnn_rollout import (  # noqa: E402
    PrioritizedReplayBuffer,
    ReplaySample,
    classify_rollout_step,
)


class TestTorchGNNRolloutUtilities(unittest.TestCase):

    def test_prioritized_replay_buffer_enforces_capacity_and_reports_stats(self):
        atoms = bulk("C", "diamond", a=3.567, cubic=True)
        atoms.pbc = True

        buffer = PrioritizedReplayBuffer(capacity=2, seed=42)
        buffer.add(
            ReplaySample(
                atoms=atoms,
                priority=0.1,
                categories=("successful",),
                stage_index=1,
                step=0,
                restart_index=0,
                true_target_fingerprint_mse=0.5,
                surrogate_target_fingerprint_mse=0.4,
                fingerprint_drift_mse=0.1,
                position_difference=0.2,
            )
        )
        buffer.add(
            ReplaySample(
                atoms=atoms,
                priority=1.0,
                categories=("failed", "high_error"),
                stage_index=1,
                step=5,
                restart_index=0,
                true_target_fingerprint_mse=0.3,
                surrogate_target_fingerprint_mse=0.2,
                fingerprint_drift_mse=0.1,
                position_difference=0.1,
            )
        )
        buffer.add(
            ReplaySample(
                atoms=atoms,
                priority=2.0,
                categories=("unstable",),
                stage_index=2,
                step=10,
                restart_index=0,
                true_target_fingerprint_mse=0.2,
                surrogate_target_fingerprint_mse=0.1,
                fingerprint_drift_mse=0.05,
                position_difference=0.05,
            )
        )

        self.assertEqual(len(buffer), 2)
        stats = buffer.stats()
        self.assertEqual(stats["capacity"], 2)
        self.assertEqual(stats["size"], 2)
        self.assertIn("failed", stats["category_counts"])
        self.assertIn("unstable", stats["category_counts"])
        self.assertGreater(stats["mean_priority"], 0.0)

        sampled = buffer.sample(
            sample_size=2,
            category_weights={"failed": 3.0, "unstable": 2.0},
        )
        self.assertEqual(len(sampled), 2)
        self.assertTrue(all(isinstance(entry, ReplaySample) for entry in sampled))

    def test_classify_rollout_step_labels_failures_and_difficult_examples(self):
        categories = classify_rollout_step(
            true_target_fingerprint_mse=1.0e-2,
            previous_true_target_fingerprint_mse=5.0e-3,
            fingerprint_drift_mse=2.0e-3,
            repulsion_loss=5.0e-4,
            cell_violation_loss=0.0,
            drift_threshold=1.0e-3,
            high_error_threshold=5.0e-3,
            instability_threshold=1.0e-4,
        )

        self.assertIn("failed", categories)
        self.assertIn("difficult", categories)
        self.assertIn("high_error", categories)
        self.assertIn("unstable", categories)


if __name__ == "__main__":
    unittest.main()
