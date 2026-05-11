import io
import importlib
import sys
import types
import unittest
from pathlib import Path
from unittest import mock


REPO_ROOT = Path(__file__).resolve().parents[1]
PYTHON_PKG_DIR = REPO_ROOT / "example" / "python_pkg"
if str(PYTHON_PKG_DIR) not in sys.path:
    sys.path.insert(0, str(PYTHON_PKG_DIR))

wandb_stub = types.ModuleType("wandb")
wandb_stub.Api = mock.Mock()

carbon_wandb_stub = types.ModuleType("torch_gnn_carbon_wandb")
carbon_wandb_stub.DEFAULT_ARCHITECTURE = "torch_gnn_residual"
carbon_wandb_stub.WANDB_PROJECT = "raffle-inverse-design-new"
carbon_wandb_stub.build_model_config = lambda config: {
    "architecture": str(config["architecture"]),
    "hidden_dim": int(config["hidden_dim"]),
    "num_message_layers": int(config["num_message_layers"]),
    "learning_rate": float(config["learning_rate"]),
    "lr_decay_rate": float(config["model_lr_decay_rate"]),
    "smooth_cutoff_width": float(config["smooth_cutoff_width"]),
    "reference_layer_type": int(config["reference_layer_type"]),
    "component_weight": (
        float(config["component_weight_2body"]),
        float(config["component_weight_3body"]),
        float(config["component_weight_4body"]),
    ),
}
carbon_wandb_stub.validate_plan_constraints = mock.Mock()

workflow_stub = types.ModuleType("torch_gnn_carbon_workflow_example")
workflow_stub.CELL_VIOLATION_WEIGHT = 0.0
workflow_stub.COORDINATE_CLIP_VALUE = None
workflow_stub.DEFAULT_MODEL_CONFIG = {
    "hidden_dim": 80,
    "num_message_layers": 2,
    "learning_rate": 5.0e-4,
    "lr_decay_rate": 5.0e-3,
    "smooth_cutoff_width": 0.2,
    "reference_layer_type": 1,
    "component_weight": (4.0, 1.0, 1.0),
}
workflow_stub.DEFAULT_REPLAY_CATEGORY_WEIGHTS = {
    "successful": 1.0,
    "failed": 2.0,
    "unstable": 2.0,
    "high_error": 2.0,
    "difficult": 2.0,
}
workflow_stub.FINGERPRINT_LOSS_WEIGHT = 0.5
workflow_stub.FIXED_LEADING_ATOMS = 1
workflow_stub.INVERSE_LR_DECAY_RATE = 0.0
workflow_stub.INVERSE_RESTART_NOISE_SCALE = 0.0
workflow_stub.INVERSE_RESTARTS = 1
workflow_stub.MINIMUM_DISTANCE_SCALE = 0.75
workflow_stub.REPLAY_BUFFER_CAPACITY = 64
workflow_stub.REPLAY_SAMPLE_SIZE = 8
workflow_stub.REPULSION_WEIGHT = 10.0
workflow_stub.ROLLOUT_DRIFT_THRESHOLD = 5.0e-4
workflow_stub.ROLLOUT_EPOCHS_PER_STAGE = 1
workflow_stub.ROLLOUT_HIGH_ERROR_THRESHOLD = 1.0e-3
workflow_stub.ROLLOUT_INSTABILITY_THRESHOLD = 1.0e-4
workflow_stub.ROLLOUT_STAGES = 1
workflow_stub.ROLLOUT_STEP_STRIDE = 25
workflow_stub.TARGET_POSITION_WEIGHT = 0.0
workflow_stub.TARGET_VERTEX_WEIGHT = 0.0
workflow_stub.build_inverse_design_options = mock.Mock(return_value={"inverse": "options"})
workflow_stub.build_perturbed_structure = mock.Mock(return_value=("perturbed", "fixed"))
workflow_stub.minimum_training_carbon_count = mock.Mock(return_value=2)
workflow_stub.parse_category_weights = lambda value: {
    key.strip(): float(raw_value.strip())
    for key, raw_value in (
        item.split("=", 1) for item in value.split(",") if item.strip()
    )
}
workflow_stub.run_rollout_retraining = mock.Mock()
workflow_stub.select_carbon_structures = mock.Mock(side_effect=lambda structures, _: list(structures))
workflow_stub.sweep_epochs = mock.Mock()

workflow_common_stub = types.ModuleType("torch_gnn_workflow_common")
workflow_common_stub.default_reference_structure = mock.Mock(return_value="reference-structure")
workflow_common_stub.load_structures = mock.Mock(return_value=["s1", "s2", "s3"])
workflow_common_stub.save_model_checkpoint = mock.Mock()
workflow_common_stub.save_target_fingerprint = mock.Mock()
workflow_common_stub.write_json = mock.Mock()
workflow_common_stub.write_structure = mock.Mock()

with mock.patch.dict(
    sys.modules,
    {
        "wandb": wandb_stub,
        "torch_gnn_carbon_wandb": carbon_wandb_stub,
        "torch_gnn_carbon_workflow_example": workflow_stub,
        "torch_gnn_workflow_common": workflow_common_stub,
    },
):
    checkpoint_cli = importlib.import_module("torch_gnn_carbon_wandb_checkpoint")


class TestTorchGNNWandbCheckpointCLI(unittest.TestCase):

    def setUp(self) -> None:
        self.model = mock.Mock()
        self.model.compute_reference_fingerprint.return_value = [0.1, 0.2, 0.3]
        workflow_stub.sweep_epochs.reset_mock()
        workflow_stub.run_rollout_retraining.reset_mock()
        workflow_common_stub.save_model_checkpoint.reset_mock()
        workflow_common_stub.save_target_fingerprint.reset_mock()
        workflow_common_stub.write_json.reset_mock()
        workflow_common_stub.write_structure.reset_mock()
        carbon_wandb_stub.validate_plan_constraints.reset_mock()

        def record_training_observer(**kwargs):
            observer = kwargs["training_observer"]
            observer(1, 0.25)
            observer(2, 0.15)
            return self.model, [], [], None, []

        workflow_stub.sweep_epochs.side_effect = record_training_observer

        def record_rollout_observer(**kwargs):
            observer = kwargs["training_observer"]
            observer(3, 0.05)
            return self.model, {"enabled": True, "num_rollout_epochs": 1, "stages": []}

        workflow_stub.run_rollout_retraining.side_effect = record_rollout_observer

    def test_resolve_run_record_uses_entity_project_path(self):
        run = mock.Mock(
            entity="test-entity",
            project=checkpoint_cli.WANDB_PROJECT,
            id="abc123",
            name="test-run",
            url="https://wandb.ai/test-entity/raffle-inverse-design-new/runs/abc123",
            config={"epochs": 10},
        )
        api = mock.Mock()
        api.run.return_value = run

        with mock.patch.object(checkpoint_cli.wandb, "Api", return_value=api) as api_factory:
            record = checkpoint_cli.resolve_run_record(
                "abc123",
                project=checkpoint_cli.WANDB_PROJECT,
                entity="test-entity",
            )

        api_factory.assert_called_once_with(overrides={"project": checkpoint_cli.WANDB_PROJECT})
        api.run.assert_called_once_with("test-entity/raffle-inverse-design-new/abc123")
        self.assertEqual(record["run_id"], "abc123")
        self.assertEqual(record["source"], "wandb-api")
        self.assertEqual(record["config"], {"epochs": 10})

    def test_parse_args_defaults_output_dir_to_current_directory(self):
        default_output_dir = REPO_ROOT / "test_current_dir"

        with mock.patch.object(checkpoint_cli.Path, "cwd", return_value=default_output_dir):
            args = checkpoint_cli.parse_args(["abc123"])

        self.assertEqual(args.output_dir, default_output_dir)

    def test_main_replays_training_and_exports_checkpoint(self):
        run = mock.Mock(
            entity="test-entity",
            project=checkpoint_cli.WANDB_PROJECT,
            id="abc123",
            name="test-run",
            url="https://wandb.ai/test-entity/raffle-inverse-design-new/runs/abc123",
            config={
                "architecture": "torch_gnn_residual",
                "hidden_dim": 96,
                "num_message_layers": 3,
                "learning_rate": 2.5e-4,
                "model_lr_decay_rate": 1.0e-3,
                "smooth_cutoff_width": 0.15,
                "reference_layer_type": 2,
                "component_weight_2body": 5.0,
                "component_weight_3body": 1.5,
                "component_weight_4body": 0.5,
                "carbon_count": -1,
                "augmented_count": 0,
                "epochs": 2,
                "batch_size": 4,
                "inverse_steps": 20,
                "inverse_step_size": 0.01,
                "fixed_leading_atoms": 1,
                "fingerprint_loss_weight": 0.75,
                "target_vertex_weight": 0.0,
                "target_position_weight": 0.0,
                "inverse_lr_decay_rate": 0.0,
                "inverse_restarts": 1,
                "inverse_restart_noise_scale": 0.0,
                "repulsion_weight": 5.0,
                "minimum_distance_scale": 0.8,
                "cell_violation_weight": 0.01,
                "coordinate_clip_value": 0.25,
                "rollout_stages": 1,
                "rollout_epochs_per_stage": 1,
                "rollout_step_stride": 10,
                "replay_buffer_capacity": 16,
                "replay_sample_size": 4,
                "replay_category_weights": "successful=1.0,failed=3.0",
                "rollout_drift_threshold": 1.0e-4,
                "rollout_high_error_threshold": 1.0e-3,
                "rollout_instability_threshold": 1.0e-5,
                "seed": 7,
            },
        )
        api = mock.Mock()
        api.run.return_value = run

        output_dir = REPO_ROOT / "build" / "test_wandb_checkpoint_cli"
        argv = [
            "abc123",
            "--entity",
            "test-entity",
            "--epochs",
            "5",
            "--output-dir",
            str(output_dir),
        ]

        with (
            mock.patch.object(checkpoint_cli.wandb, "Api", return_value=api),
            mock.patch("sys.stdout", new_callable=io.StringIO) as stdout,
        ):
            checkpoint_cli.main(argv)

        carbon_wandb_stub.validate_plan_constraints.assert_called_once()
        workflow_stub.sweep_epochs.assert_called_once()
        workflow_stub.run_rollout_retraining.assert_called_once()
        self.assertEqual(workflow_stub.sweep_epochs.call_args.kwargs["epoch_values"], [5])

        checkpoint_call = workflow_common_stub.save_model_checkpoint.call_args.kwargs
        self.assertEqual(checkpoint_call["model"], self.model)
        self.assertEqual(checkpoint_call["species_list"], ["C"])
        self.assertEqual(checkpoint_call["training_history"], [0.25, 0.15, 0.05])
        self.assertEqual(checkpoint_call["model_config"]["architecture"], "torch_gnn_residual")
        self.assertEqual(checkpoint_call["model_config"]["component_weight"], (5.0, 1.5, 0.5))

        training_config = checkpoint_call["training_config"]
        self.assertEqual(training_config["source_wandb_run"]["run_id"], "abc123")
        self.assertEqual(training_config["workflow_config"]["batch_size"], 4)
        self.assertEqual(training_config["workflow_config"]["epochs"], 5)
        self.assertEqual(
            training_config["workflow_config"]["replay_category_weights"],
            {"successful": 1.0, "failed": 3.0},
        )

        metrics_payload = workflow_common_stub.write_json.call_args_list[-1].args[1]
        self.assertEqual(metrics_payload["final_training_loss"], 0.05)
        self.assertEqual(
            metrics_payload["output_files"]["checkpoint"],
            str(output_dir / "torch_gnn_model_checkpoint.pt"),
        )
        self.assertIn("torch_gnn_model_checkpoint.pt", stdout.getvalue())


if __name__ == "__main__":
    unittest.main()
