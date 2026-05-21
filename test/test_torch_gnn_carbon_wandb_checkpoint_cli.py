import io
import importlib
import json
import sys
import tempfile
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
carbon_wandb_stub.REPRODUCIBILITY_METADATA_KEY = "reproducibility_metadata_json"
carbon_wandb_stub.RESOLVED_WORKFLOW_CONFIG_KEY = "resolved_workflow_config_json"
carbon_wandb_stub.WANDB_PROJECT = "raffle-inverse-design-ensemble-new2"
carbon_wandb_stub.parse_serialized_run_payload = lambda value, *, key: (
    json.loads(value) if isinstance(value, str) else dict(value)
)
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
workflow_stub.select_carbon_structures = mock.Mock(side_effect=lambda structures, _: list(structures))
workflow_stub.sweep_epochs = mock.Mock()

workflow_common_stub = types.ModuleType("torch_gnn_workflow_common")
workflow_common_stub.default_reference_structure = mock.Mock(return_value="reference-structure")
workflow_common_stub.load_structures = mock.Mock(return_value=["s1", "s2", "s3"])
workflow_common_stub.load_model_from_checkpoint = mock.Mock(
    return_value=(
        mock.Mock(),
        {
            "training_history": [1.0, 0.5],
            "species_list": ["C"],
            "model_config": {"architecture": "torch_gnn_residual"},
            "training_config": {"seed": 42},
        },
    )
)
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
        workflow_common_stub.save_model_checkpoint.reset_mock()
        workflow_common_stub.save_target_fingerprint.reset_mock()
        workflow_common_stub.load_model_from_checkpoint.reset_mock()
        workflow_common_stub.write_json.reset_mock()
        workflow_common_stub.write_structure.reset_mock()
        carbon_wandb_stub.validate_plan_constraints.reset_mock()
        self.restore_reference_structure = mock.patch.object(
            checkpoint_cli,
            "_restore_reference_structure",
            return_value="reference-structure",
        )
        self.restore_perturbed_structure = mock.patch.object(
            checkpoint_cli,
            "_restore_perturbed_structure",
            return_value=("perturbed", "fixed"),
        )
        self.restore_reference_structure.start()
        self.restore_perturbed_structure.start()
        self.addCleanup(self.restore_reference_structure.stop)
        self.addCleanup(self.restore_perturbed_structure.stop)

        def record_training_observer(**kwargs):
            observer = kwargs["training_observer"]
            observer(1, 0.25)
            observer(2, 0.15)
            observer(3, 0.05)
            return self.model, [], [], None, [], {
                "enabled": True,
                "num_rollout_epochs": 1,
                "stages": [],
            }

        workflow_stub.sweep_epochs.side_effect = record_training_observer

    def test_resolve_run_record_uses_entity_project_path(self):
        project = "custom-project"
        run = mock.Mock(
            entity="test-entity",
            project=project,
            id="abc123",
            name="test-run",
            url=f"https://wandb.ai/test-entity/{project}/runs/abc123",
            config={"epochs": 10},
        )
        api = mock.Mock()
        api.run.return_value = run

        with mock.patch.object(checkpoint_cli.wandb, "Api", return_value=api) as api_factory:
            record = checkpoint_cli.resolve_run_record(
                "abc123",
                project=project,
                entity="test-entity",
            )

        api_factory.assert_called_once_with(overrides={"project": project})
        api.run.assert_called_once_with(f"test-entity/{project}/abc123")
        self.assertEqual(record["run_id"], "abc123")
        self.assertEqual(record["source"], "wandb-api")
        self.assertEqual(record["config"], {"epochs": 10})

    def test_parse_args_defaults_output_dir_to_current_directory(self):
        default_output_dir = REPO_ROOT / "test_current_dir"

        with mock.patch.object(checkpoint_cli.Path, "cwd", return_value=default_output_dir):
            args = checkpoint_cli.parse_args(["abc123"])

        self.assertEqual(args.output_dir, default_output_dir)

    def test_resolve_run_record_prefers_persisted_local_wandb_config(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            script_dir = Path(temp_dir)
            run_dir = script_dir / "wandb" / "run-20260518_084407-abc123"
            files_dir = run_dir / "files"
            logs_dir = run_dir / "logs"
            files_dir.mkdir(parents=True)
            logs_dir.mkdir(parents=True)

            (files_dir / "config.yaml").write_text(
                "architecture:\n"
                "    value: torch_gnn_residual\n"
                "hidden_dim:\n"
                "    value: 80\n"
                "num_message_layers:\n"
                "    value: 4\n"
                "learning_rate:\n"
                "    value: 0.0005\n"
                "model_lr_decay_rate:\n"
                "    value: 0.01\n"
                "smooth_cutoff_width:\n"
                "    value: 0.3\n"
                "reference_layer_type:\n"
                "    value: 1\n"
                "component_weight_2body:\n"
                "    value: 2\n"
                "component_weight_3body:\n"
                "    value: 2\n"
                "component_weight_4body:\n"
                "    value: 3\n"
                "epochs:\n"
                "    value: 50\n"
                "ensemble_enabled:\n"
                "    value: false\n"
                "coordinate_clip_value:\n"
                "    value: null\n"
            )
            (logs_dir / "debug.log").write_text(
                "config: {'hidden_dim': 80, 'num_message_layers': 2, 'epochs': 25}\n"
            )

            with mock.patch.object(checkpoint_cli, "SCRIPT_DIR", script_dir):
                record = checkpoint_cli.resolve_run_record(
                    "abc123",
                    project=checkpoint_cli.WANDB_PROJECT,
                )

        self.assertEqual(record["source"], "local-cache")
        self.assertEqual(record["config"]["num_message_layers"], 4)
        self.assertEqual(record["config"]["component_weight_2body"], 2)
        self.assertFalse(record["config"]["ensemble_enabled"])
        self.assertIsNone(record["config"]["coordinate_clip_value"])

    def test_resolve_run_record_uses_remote_config_when_entity_is_provided(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            project = "custom-project"
            script_dir = Path(temp_dir)
            run_dir = script_dir / "wandb" / "run-20260518_084407-abc123"
            files_dir = run_dir / "files"
            files_dir.mkdir(parents=True)
            (files_dir / "config.yaml").write_text(
                "epochs:\n"
                "    value: 50\n"
            )

            run = mock.Mock(
                entity="test-entity",
                project=project,
                id="abc123",
                name="remote-run",
                url=f"https://wandb.ai/test-entity/{project}/runs/abc123",
                config={"epochs": 10},
            )
            api = mock.Mock()
            api.run.return_value = run

            with (
                mock.patch.object(checkpoint_cli, "SCRIPT_DIR", script_dir),
                mock.patch.object(checkpoint_cli.wandb, "Api", return_value=api),
            ):
                record = checkpoint_cli.resolve_run_record(
                    "abc123",
                    project=project,
                    entity="test-entity",
                )

        api.run.assert_called_once_with(f"test-entity/{project}/abc123")
        self.assertEqual(record["source"], "wandb-api")
        self.assertEqual(record["config"], {"epochs": 10})

    def test_build_replay_config_requires_serialized_resolved_workflow_config(self):
        with self.assertRaisesRegex(ValueError, "resolved_workflow_config_json"):
            checkpoint_cli.build_replay_config({"epochs": 10})

    def test_main_exports_exact_local_checkpoint_by_default(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            script_dir = Path(temp_dir)
            run_dir = script_dir / "wandb" / "run-20260518_084407-abc123"
            files_dir = run_dir / "files"
            logs_dir = run_dir / "logs"
            source_output_dir = (
                script_dir / "build" / "wandb_inverse_design" / "torch_gnn_residual" / "abc123"
            )
            export_dir = script_dir / "exported"
            files_dir.mkdir(parents=True)
            logs_dir.mkdir(parents=True)
            source_output_dir.mkdir(parents=True)

            (files_dir / "config.yaml").write_text(
                "architecture:\n"
                "    value: torch_gnn_residual\n"
                "hidden_dim:\n"
                "    value: 80\n"
                "num_message_layers:\n"
                "    value: 4\n"
            )
            (files_dir / "wandb-summary.json").write_text(
                '{"output_dir": "build/wandb_inverse_design/torch_gnn_residual/abc123"}'
            )
            (logs_dir / "debug.log").write_text(
                "2026-05-18 08:45:53,798 INFO finishing run test-entity/"
                "raffle-inverse-design-ensemble-new/abc123\n"
            )
            (source_output_dir / "torch_gnn_model_checkpoint.pt").write_bytes(
                b"exact-checkpoint"
            )
            (source_output_dir / "torch_gnn_target_fingerprint.npy").write_bytes(
                b"fingerprint"
            )
            (source_output_dir / "torch_gnn_reference_structure.xyz").write_text(
                "reference"
            )

            workflow_common_stub.load_model_from_checkpoint.return_value = (
                mock.Mock(),
                {
                    "training_history": [0.25, 0.05],
                    "species_list": ["C"],
                    "model_config": {"architecture": "torch_gnn_residual"},
                    "training_config": {"seed": 42},
                },
            )

            with (
                mock.patch.object(checkpoint_cli, "SCRIPT_DIR", script_dir),
                mock.patch("sys.stdout", new_callable=io.StringIO) as stdout,
            ):
                checkpoint_cli.main(["abc123", "--output-dir", str(export_dir)])

            workflow_stub.sweep_epochs.assert_not_called()
            workflow_common_stub.save_model_checkpoint.assert_not_called()
            self.assertEqual(
                workflow_common_stub.load_model_from_checkpoint.call_args.args[0].resolve(),
                (export_dir / "torch_gnn_model_checkpoint.pt").resolve(),
            )
            self.assertEqual(
                (export_dir / "torch_gnn_model_checkpoint.pt").read_bytes(),
                b"exact-checkpoint",
            )
            self.assertEqual(
                (export_dir / "torch_gnn_target_fingerprint.npy").read_bytes(),
                b"fingerprint",
            )
            self.assertEqual(
                (export_dir / "torch_gnn_reference_structure.xyz").read_text(),
                "reference",
            )
            exact_metrics_payload = workflow_common_stub.write_json.call_args_list[-1].args[1]
            self.assertEqual(exact_metrics_payload["recovery_mode"], "exact-checkpoint")
            self.assertEqual(
                exact_metrics_payload["source_wandb_run"]["project"],
                "raffle-inverse-design-ensemble-new",
            )
            self.assertIn("torch_gnn_model_checkpoint.pt", stdout.getvalue())

    def test_main_downloads_exact_remote_checkpoint_artifact(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            artifact_dir = Path(temp_dir) / "downloaded-artifact"
            export_dir = Path(temp_dir) / "exported"
            artifact_dir.mkdir(parents=True)
            (artifact_dir / "torch_gnn_model_checkpoint.pt").write_bytes(b"remote-checkpoint")
            (artifact_dir / "torch_gnn_target_fingerprint.npy").write_bytes(b"fingerprint")
            (artifact_dir / "torch_gnn_reference_structure.xyz").write_text("reference")

            run = mock.Mock(
                entity="test-entity",
                project=checkpoint_cli.WANDB_PROJECT,
                id="abc123",
                name="remote-run",
                url=(
                    "https://wandb.ai/test-entity/"
                    f"{checkpoint_cli.WANDB_PROJECT}/runs/abc123"
                ),
                config={"epochs": 10},
            )
            artifact = mock.Mock(type="inverse-design-results")
            artifact.name = "inverse-design-results:latest"
            artifact.download.return_value = str(artifact_dir)
            run.logged_artifacts.return_value = [artifact]

            workflow_common_stub.load_model_from_checkpoint.return_value = (
                mock.Mock(),
                {
                    "training_history": [0.25, 0.05],
                    "species_list": ["C"],
                    "model_config": {"architecture": "torch_gnn_residual"},
                    "training_config": {"seed": 42},
                },
            )

            api = mock.Mock()
            api.run.return_value = run

            with (
                mock.patch.object(checkpoint_cli.wandb, "Api", return_value=api),
                mock.patch("sys.stdout", new_callable=io.StringIO) as stdout,
            ):
                checkpoint_cli.main(
                    [
                        "abc123",
                        "--entity",
                        "test-entity",
                        "--output-dir",
                        str(export_dir),
                    ]
                )

            workflow_stub.sweep_epochs.assert_not_called()
            artifact.download.assert_called_once()
            self.assertEqual(
                workflow_common_stub.load_model_from_checkpoint.call_args.args[0].resolve(),
                (export_dir / "torch_gnn_model_checkpoint.pt").resolve(),
            )
            self.assertEqual(
                (export_dir / "torch_gnn_model_checkpoint.pt").read_bytes(),
                b"remote-checkpoint",
            )
            exact_metrics_payload = workflow_common_stub.write_json.call_args_list[-1].args[1]
            self.assertEqual(exact_metrics_payload["checkpoint_source"]["mode"], "wandb-artifact")
            self.assertEqual(
                exact_metrics_payload["checkpoint_source"]["artifact_name"],
                "inverse-design-results:latest",
            )
            self.assertIn("torch_gnn_model_checkpoint.pt", stdout.getvalue())

    def test_main_replays_training_and_exports_checkpoint(self):
        resolved_workflow_config = {
            "spec_version": 1,
            "architecture": "torch_gnn_residual",
            "architecture_name": "torch_gnn_residual",
            "carbon_dataset_path": "example/data/carbon.xyz",
            "carbon_count": -1,
            "augmented_count": 0,
            "epochs": 2,
            "batch_size": 4,
            "inverse_steps": 20,
            "inverse_step_size": 0.01,
            "inverse_step_values": [0, 20],
            "step_size_values": [0.01],
            "fixed_leading_atoms": 1,
            "fingerprint_loss_weight": 0.75,
            "target_vertex_weight": 0.0,
            "target_position_weight": 0.0,
            "reference_structure_weight": 5.0,
            "inverse_lr_decay_rate": 0.0,
            "repulsion_weight": 5.0,
            "minimum_distance_scale": 0.8,
            "cell_violation_weight": 0.01,
            "coordinate_clip_value": 0.25,
            "rollout_stages": 1,
            "rollout_epochs_per_stage": 1,
            "rollout_step_stride": 10,
            "replay_buffer_capacity": 16,
            "replay_sample_size": 4,
            "replay_category_weights": {"successful": 1.0, "failed": 3.0},
            "rollout_drift_threshold": 1.0e-4,
            "rollout_high_error_threshold": 1.0e-3,
            "rollout_instability_threshold": 1.0e-5,
            "ignored_deprecated_config_fields": {
                "inverse_restarts": 1,
                "inverse_restart_noise_scale": 0.0,
            },
            "seed": 7,
            "model_config": {
                "architecture": "torch_gnn_residual",
                "hidden_dim": 96,
                "num_message_layers": 3,
                "learning_rate": 2.5e-4,
                "lr_decay_rate": 1.0e-3,
                "smooth_cutoff_width": 0.15,
                "reference_layer_type": 2,
                "component_weight": [5.0, 1.5, 0.5],
            },
            "reference_structure": {
                "symbols": ["C"],
                "positions": [[0.0, 0.0, 0.0]],
                "cell": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
                "pbc": [True, True, True],
            },
            "perturbation_matrix": [[0.0, 0.0, 0.0]],
            "enable_checkpoint_step_size_sweep": True,
            "enable_checkpoint_step_schedule_sweep": True,
        }
        run = mock.Mock(
            entity="test-entity",
            project=checkpoint_cli.WANDB_PROJECT,
            id="abc123",
            name="test-run",
            url=(
                "https://wandb.ai/test-entity/"
                f"{checkpoint_cli.WANDB_PROJECT}/runs/abc123"
            ),
            config={
                checkpoint_cli.RESOLVED_WORKFLOW_CONFIG_KEY: json.dumps(
                    resolved_workflow_config
                ),
                checkpoint_cli.REPRODUCIBILITY_METADATA_KEY: json.dumps(
                    {"git_commit": "abc123"}
                ),
            },
        )
        api = mock.Mock()
        api.run.return_value = run

        output_dir = REPO_ROOT / "build" / "test_wandb_checkpoint_cli"
        argv = [
            "abc123",
            "--entity",
            "test-entity",
            "--replay-training",
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
        self.assertEqual(workflow_stub.sweep_epochs.call_args.kwargs["num_epochs"], 5)
        self.assertNotIn(
            "reference_structure_weight",
            workflow_stub.build_inverse_design_options.call_args.kwargs,
        )
        self.assertNotIn(
            "original",
            workflow_stub.build_inverse_design_options.call_args.kwargs,
        )

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
        self.assertNotIn("reference_structure_weight", training_config["workflow_config"])
        self.assertEqual(
            training_config["workflow_config"]["ignored_deprecated_config_fields"],
            {"inverse_restarts": 1, "inverse_restart_noise_scale": 0.0},
        )
        self.assertEqual(
            training_config["workflow_config"]["replay_category_weights"],
            {"successful": 1.0, "failed": 3.0},
        )
        self.assertEqual(training_config["source_reproducibility_metadata"], {"git_commit": "abc123"})

        metrics_payload = workflow_common_stub.write_json.call_args_list[-1].args[1]
        self.assertEqual(metrics_payload["final_training_loss"], 0.05)
        self.assertEqual(metrics_payload["recovery_mode"], "config-replay")
        self.assertEqual(
            metrics_payload["output_files"]["checkpoint"],
            str(output_dir / "torch_gnn_model_checkpoint.pt"),
        )
        self.assertIn("torch_gnn_model_checkpoint.pt", stdout.getvalue())


if __name__ == "__main__":
    unittest.main()
