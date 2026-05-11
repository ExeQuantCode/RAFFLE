import math
import json
import sys
import unittest
from pathlib import Path

import numpy as np
from ase.build import bulk
from ase.io import read

from raffle import TorchGNNFingerprint


REPO_ROOT = Path(__file__).resolve().parents[1]
PYTHON_PKG_DIR = REPO_ROOT / "example" / "python_pkg"
if str(PYTHON_PKG_DIR) not in sys.path:
    sys.path.insert(0, str(PYTHON_PKG_DIR))

from torch_gnn_carbon_workflow_example import (  # noqa: E402
    DEFAULT_MODEL_CONFIG,
    PERTURBATION,
    DEFAULT_REPLAY_CATEGORY_WEIGHTS,
    run_example,
)


@unittest.skipUnless(TorchGNNFingerprint is not None, "PyTorch multigraph fingerprint is unavailable")
class TestTorchGNNFingerprintCarbonWorkflow(unittest.TestCase):

    def test_inverse_design_rejects_target_supervision_losses(self):
        original = bulk("C", "diamond", a=3.567, cubic=True)
        original.pbc = True
        perturbed = original.copy()
        perturbed.set_positions(perturbed.get_positions() + PERTURBATION)
        fixed_atoms = np.zeros(len(perturbed), dtype=bool)

        model = TorchGNNFingerprint(
            species_list=["C"],
            hidden_dim=32,
            num_message_layers=1,
            seed=42,
        )
        target_fingerprint = model.compute_reference_fingerprint(original)

        with self.assertRaisesRegex(ValueError, "target_vertex_weight"):
            model.inverse_design(
                target_fingerprint=target_fingerprint,
                atoms=perturbed,
                fixed_atoms=fixed_atoms,
                num_steps=1,
                step_size=1.0e-2,
                target_atoms=original,
                target_vertex_weight=0.1,
            )

        with self.assertRaisesRegex(ValueError, "target_position_weight"):
            model.inverse_design(
                target_fingerprint=target_fingerprint,
                atoms=perturbed,
                fixed_atoms=fixed_atoms,
                num_steps=1,
                step_size=1.0e-2,
                target_atoms=original,
                target_position_weight=0.1,
            )

    def test_run_example_emits_rollout_and_convergence_metrics(self):
        carbon_xyz = REPO_ROOT / "example" / "data" / "carbon.xyz"
        all_carbon_structures = read(str(carbon_xyz), index=":")
        minimum_carbon_count = max(1, int(math.ceil(0.3 * len(all_carbon_structures))))

        output_dir = REPO_ROOT / "build" / "test_torch_gnn_carbon_workflow"
        metrics = run_example(
            repo_root=REPO_ROOT,
            output_dir=output_dir,
            carbon_count=minimum_carbon_count,
            augmented_count=0,
            num_epochs=2,
            batch_size=8,
            inverse_steps=10,
            inverse_step_size=0.02,
            epoch_values=[0, 1, 2],
            inverse_step_values=[0, 5, 10],
            step_size_values=[0.01, 0.02],
            fixed_leading_atoms=0,
            fingerprint_loss_weight=1.0,
            target_vertex_weight=0.0,
            target_position_weight=0.0,
            inverse_lr_decay_rate=0.0,
            inverse_restarts=1,
            inverse_restart_noise_scale=0.0,
            repulsion_weight=10.0,
            minimum_distance_scale=0.75,
            cell_violation_weight=0.01,
            coordinate_clip_value=0.2,
            rollout_stages=1,
            rollout_epochs_per_stage=1,
            rollout_step_stride=5,
            replay_buffer_capacity=8,
            replay_sample_size=4,
            replay_category_weights=dict(DEFAULT_REPLAY_CATEGORY_WEIGHTS),
            rollout_drift_threshold=1.0e-4,
            rollout_high_error_threshold=1.0e-4,
            rollout_instability_threshold=1.0e-5,
            seed=42,
            model_config={
                **DEFAULT_MODEL_CONFIG,
                "architecture": "residual",
                "hidden_dim": 32,
                "num_message_layers": 1,
                "learning_rate": 5.0e-4,
                "lr_decay_rate": 1.0e-3,
                "smooth_cutoff_width": 0.2,
                "reference_layer_type": 1,
                "component_weight": (4.0, 1.0, 1.0),
            },
            architecture_name="torch_gnn_residual",
            enable_checkpoint_step_size_sweep=False,
            enable_checkpoint_step_schedule_sweep=False,
        )

        self.assertGreaterEqual(
            metrics["num_training_carbon_structures"],
            minimum_carbon_count,
        )
        self.assertEqual(metrics["inverse_design_config"]["target_vertex_weight"], 0.0)
        self.assertEqual(metrics["inverse_design_config"]["target_position_weight"], 0.0)
        self.assertIn("convergence_summary", metrics)
        self.assertIn("rollout", metrics)
        self.assertTrue(metrics["rollout"]["enabled"])
        self.assertEqual(len(metrics["rollout"]["stages"]), 1)
        self.assertGreater(metrics["rollout"]["replay_buffer"]["size"], 0)
        self.assertGreater(metrics["rollout"]["stages"][0]["sampled_replay_count"], 0)
        self.assertIn(
            metrics["convergence_summary"]["best_step"],
            [entry["step"] for entry in metrics["inverse_step_sweep"]],
        )
        self.assertIn("fingerprint_mse", metrics["inverse_step_sweep"][0])
        self.assertIn("fingerprint_l2", metrics["inverse_step_sweep"][0])
        self.assertLessEqual(metrics["final_rmsd"], metrics["initial_rmsd"] + 1.0e-8)
        self.assertIn("descriptor_comparisons", metrics)
        descriptor_report = metrics["descriptor_comparisons"]["final"]
        self.assertEqual(descriptor_report["structure_file"], metrics["output_files"]["final"])
        self.assertIn("global_metrics", descriptor_report)
        self.assertIn("components", descriptor_report)
        self.assertIn("2body", descriptor_report["components"])
        self.assertEqual(
            len(descriptor_report["global_fingerprints"]["target_raffle"]),
            len(descriptor_report["global_fingerprints"]["ml_predicted"]),
        )
        self.assertTrue(Path(metrics["output_files"]["plot"]).exists())
        self.assertTrue(Path(metrics["output_files"]["metrics"]).exists())
        self.assertTrue(Path(metrics["output_files"]["configured_inverse_design_traj"]).exists())
        descriptor_path = Path(metrics["output_files"]["final_descriptor_comparison"])
        descriptor_plot_path = Path(metrics["output_files"]["final_descriptor_comparison_plot"])
        self.assertTrue(descriptor_path.exists())
        self.assertTrue(descriptor_plot_path.exists())
        saved_descriptor_report = json.loads(descriptor_path.read_text())
        self.assertEqual(saved_descriptor_report["structure_file"], metrics["output_files"]["final"])
        self.assertEqual(saved_descriptor_report["plot_file"], metrics["output_files"]["final_descriptor_comparison_plot"])
        configured_path = metrics["configured_inverse_design_path"]
        self.assertTrue(Path(configured_path["step_structure_dir"]).exists())
        self.assertGreater(len(configured_path["step_structure_files"]), 0)
        self.assertTrue(all(Path(path).exists() for path in configured_path["step_structure_files"]))
        saved_path_frames = read(metrics["output_files"]["configured_inverse_design_traj"], index=":")
        self.assertEqual(len(saved_path_frames), len(configured_path["step_structure_files"]))


if __name__ == "__main__":
    unittest.main()
