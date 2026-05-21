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
    DEFAULT_REPLAY_CATEGORY_WEIGHTS,
    build_perturbed_structure,
    run_example,
)
from torch_gnn_workflow_common import (  # noqa: E402
    compute_inverse_design_metrics,
    load_model_from_checkpoint,
    load_target_fingerprint,
)


@unittest.skipUnless(TorchGNNFingerprint is not None, "PyTorch multigraph fingerprint is unavailable")
class TestTorchGNNFingerprintCarbonWorkflow(unittest.TestCase):

    def test_inverse_design_rejects_target_supervision_losses(self):
        original = bulk("C", "diamond", a=3.567, cubic=True)
        original.pbc = True
        perturbed, fixed_atoms = build_perturbed_structure(
            original,
            fixed_leading_atoms=0,
            seed=42,
        )

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

    def test_inverse_design_allows_different_target_atom_count(self):
        original = bulk("C", "diamond", a=3.567, cubic=True)
        original.pbc = True
        perturbed, fixed_atoms = build_perturbed_structure(
            original,
            fixed_leading_atoms=0,
            seed=7,
        )
        different_target = bulk("C", "diamond", a=3.567, cubic=False)
        different_target = different_target[[0, 1]]
        different_target.pbc = True

        model = TorchGNNFingerprint(
            species_list=["C"],
            hidden_dim=32,
            num_message_layers=1,
            seed=42,
        )
        target_fingerprint = model.compute_reference_fingerprint(original)

        optimised = model.inverse_design(
            target_fingerprint=target_fingerprint,
            atoms=perturbed,
            fixed_atoms=fixed_atoms,
            num_steps=1,
            step_size=1.0e-2,
            target_atoms=different_target,
        )

        self.assertEqual(len(optimised), len(perturbed))

    def test_run_example_emits_rollout_and_convergence_metrics(self):
        carbon_xyz = REPO_ROOT / "example" / "data" / "carbon.xyz"
        all_carbon_structures = read(str(carbon_xyz), index=":")
        minimum_carbon_count = max(1, int(math.ceil(0.3 * len(all_carbon_structures))))

        output_dir = REPO_ROOT / "build" / "test_torch_gnn_carbon_workflow"
        metrics = run_example(
            repo_root=REPO_ROOT,
            output_dir=output_dir,
            carbon_count=minimum_carbon_count,
            augmented_count=2,
            num_epochs=2,
            batch_size=8,
            inverse_steps=10,
            inverse_step_size=0.02,
            inverse_step_values=[0, 5, 10],
            step_size_values=[0.01, 0.02],
            fixed_leading_atoms=0,
            fingerprint_loss_weight=1.0,
            target_vertex_weight=0.0,
            target_position_weight=0.0,
            inverse_lr_decay_rate=0.0,
            inverse_restarts=None,
            inverse_restart_noise_scale=None,
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
        self.assertEqual(
            metrics["num_perturbed_training_structures_per_epoch"],
            metrics["num_training_carbon_structures"] * 2,
        )
        self.assertEqual(metrics["inverse_design_config"]["target_vertex_weight"], 0.0)
        self.assertEqual(metrics["inverse_design_config"]["target_position_weight"], 0.0)
        self.assertNotIn("reference_structure_weight", metrics["inverse_design_config"])
        self.assertNotIn("epoch_values", metrics["inverse_design_config"])
        self.assertNotIn("inverse_restarts", metrics["inverse_design_config"])
        self.assertNotIn("inverse_restart_noise_scale", metrics["inverse_design_config"])
        self.assertEqual(metrics["inverse_design_config"]["resolved_inverse_steps"], 10)
        self.assertEqual(metrics["inverse_design_config"]["inverse_step_log_interval"], 10)
        self.assertIn("convergence_summary", metrics)
        self.assertIn("best_total_loss_step", metrics["convergence_summary"])
        self.assertIn(
            "position_difference_at_best_total_loss",
            metrics["convergence_summary"],
        )
        self.assertIn("rollout", metrics)
        self.assertTrue(metrics["rollout"]["enabled"])
        self.assertEqual(len(metrics["rollout"]["stages"]), 1)
        self.assertGreater(metrics["rollout"]["replay_buffer"]["size"], 0)
        self.assertGreater(metrics["rollout"]["stages"][0]["sampled_replay_count"], 0)
        self.assertIn("evaluation_pair", metrics)
        self.assertIn("structure_index", metrics["evaluation_pair"])
        self.assertIn("perturbation_settings", metrics)
        self.assertEqual(
            [entry["step"] for entry in metrics["inverse_step_sweep"]],
            [0, 10],
        )
        self.assertIn(
            metrics["convergence_summary"]["best_step"],
            [entry["step"] for entry in metrics["inverse_step_sweep"]],
        )
        self.assertIn("fingerprint_mse", metrics["inverse_step_sweep"][0])
        self.assertIn("fingerprint_l2", metrics["inverse_step_sweep"][0])
        self.assertLessEqual(
            metrics["final_fingerprint_mse"],
            metrics["initial_fingerprint_mse"] + 1.0e-12,
        )
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
        self.assertTrue(Path(metrics["output_files"]["checkpoint"]).exists())
        self.assertTrue(Path(metrics["output_files"]["target_fingerprint"]).exists())
        self.assertTrue(Path(metrics["output_files"]["reference_structure"]).exists())
        self.assertTrue(Path(metrics["output_files"]["metrics"]).exists())
        self.assertTrue(Path(metrics["output_files"]["configured_inverse_design_traj"]).exists())
        descriptor_path = Path(metrics["output_files"]["final_descriptor_comparison"])
        descriptor_plot_path = Path(metrics["output_files"]["final_descriptor_comparison_plot"])
        self.assertTrue(descriptor_path.exists())
        self.assertTrue(descriptor_plot_path.exists())
        reloaded_model, checkpoint_payload = load_model_from_checkpoint(
            Path(metrics["output_files"]["checkpoint"])
        )
        self.assertEqual(checkpoint_payload["training_history"], metrics["training_losses"])
        self.assertEqual(checkpoint_payload["training_config"]["seed"], 42)
        self.assertEqual(checkpoint_payload["model_config"]["seed"], 42)
        self.assertEqual(
            checkpoint_payload["training_config"]["workflow_config"][
                "ignored_deprecated_config_fields"
            ],
            {},
        )
        np.testing.assert_allclose(
            load_target_fingerprint(Path(metrics["output_files"]["target_fingerprint"])),
            reloaded_model.compute_reference_fingerprint(
                read(metrics["output_files"]["reference_structure"])
            ),
        )
        saved_descriptor_report = json.loads(descriptor_path.read_text())
        self.assertEqual(saved_descriptor_report["structure_file"], metrics["output_files"]["final"])
        self.assertEqual(saved_descriptor_report["plot_file"], metrics["output_files"]["final_descriptor_comparison_plot"])
        configured_path = metrics["configured_inverse_design_path"]
        self.assertTrue(Path(configured_path["step_structure_dir"]).exists())
        self.assertGreater(len(configured_path["step_structure_files"]), 0)
        self.assertTrue(all(Path(path).exists() for path in configured_path["step_structure_files"]))
        saved_path_frames = read(metrics["output_files"]["configured_inverse_design_traj"], index=":")
        self.assertEqual(len(saved_path_frames), len(configured_path["step_structure_files"]))

    def test_compute_inverse_design_metrics_skips_rmsd_for_atom_count_mismatch(self):
        initial_atoms = bulk("C", "diamond", a=3.567, cubic=True)
        initial_atoms.pbc = True
        optimised_atoms = initial_atoms.copy()
        optimised_atoms.set_positions(initial_atoms.get_positions() + 0.01)
        target_atoms = initial_atoms[[0, 1]].copy()
        target_atoms.pbc = True

        class StubModel:

            def predict(self, atoms):
                return np.full(3, float(len(atoms)), dtype=np.float32)

        metrics = compute_inverse_design_metrics(
            model=StubModel(),
            initial_atoms=initial_atoms,
            optimised_atoms=optimised_atoms,
            target_fingerprint=np.asarray([1.0, 1.0, 1.0], dtype=np.float32),
            target_atoms=target_atoms,
        )

        self.assertEqual(metrics["structure_comparison_skipped_reason"], "atom_count_mismatch")
        self.assertEqual(metrics["target_atom_count"], 2)
        self.assertNotIn("initial_rmsd", metrics)
        self.assertNotIn("final_rmsd", metrics)


if __name__ == "__main__":
    unittest.main()
