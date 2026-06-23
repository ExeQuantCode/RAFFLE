import importlib
import io
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
wandb_stub.agent = mock.Mock()
wandb_stub.sweep = mock.Mock(return_value="sweep-123")
wandb_stub.init = mock.Mock()
wandb_stub.log = mock.Mock()
wandb_stub.Table = mock.Mock()
wandb_stub.Image = mock.Mock()
wandb_stub.Artifact = mock.Mock()
wandb_stub.__version__ = "0.0.0-test"

shared_wandb_stub = types.ModuleType("torch_gnn_carbon_wandb")
shared_wandb_stub.build_final_validation_logging_payload = mock.Mock(return_value={})
shared_wandb_stub.build_reproducibility_metadata = mock.Mock(return_value={})
shared_wandb_stub.build_resolved_workflow_config = mock.Mock(return_value={})
shared_wandb_stub.iter_result_artifact_paths = mock.Mock(return_value=[])
shared_wandb_stub.persist_reproducibility_payloads = mock.Mock()
shared_wandb_stub.summarise_final_validation_cases = mock.Mock(return_value=[])

workflow_stub = types.ModuleType("torch_gnn_carbon_workflow_example")
workflow_stub.CELL_VIOLATION_WEIGHT = 0.0
workflow_stub.COORDINATE_CLIP_VALUE = None
workflow_stub.WRAP_POSITIONS_TO_CELL = True
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
workflow_stub.INVERSE_RESTARTS = 0
workflow_stub.MINIMUM_DISTANCE_SCALE = 0.75
workflow_stub.REPLAY_BUFFER_CAPACITY = 8
workflow_stub.REPLAY_SAMPLE_SIZE = 4
workflow_stub.REPULSION_WEIGHT = 10.0
workflow_stub.ROLLOUT_DRIFT_THRESHOLD = 1.0e-4
workflow_stub.ROLLOUT_EPOCHS_PER_STAGE = 1
workflow_stub.ROLLOUT_HIGH_ERROR_THRESHOLD = 1.0e-4
workflow_stub.ROLLOUT_INSTABILITY_THRESHOLD = 1.0e-5
workflow_stub.ROLLOUT_STAGES = 1
workflow_stub.ROLLOUT_STEP_STRIDE = 5
workflow_stub.TARGET_VERTEX_WEIGHT = 0.0
workflow_stub.default_inverse_step_values = lambda steps: [steps]
workflow_stub.default_step_sizes = lambda step_size: [step_size]
workflow_stub.parse_category_weights = lambda value: {}
workflow_stub.parse_float_list = lambda value: []
workflow_stub.parse_int_list = lambda value: []
workflow_stub.run_example = mock.Mock()

with mock.patch.dict(
    sys.modules,
    {
        "wandb": wandb_stub,
        "torch_gnn_carbon_wandb": shared_wandb_stub,
        "torch_gnn_carbon_workflow_example": workflow_stub,
    },
):
    ensemble_wandb = importlib.import_module("torch_gnn_carbon_ensemble_wandb")


class TestTorchGNNEnsembleSweepCLI(unittest.TestCase):

    def setUp(self) -> None:
        ensemble_wandb.wandb.sweep.reset_mock(return_value=True)
        ensemble_wandb.wandb.sweep.return_value = "sweep-123"

    def test_resolve_project_name_prefers_input_argument(self):
        self.assertEqual(
            ensemble_wandb.resolve_project_name("cli-project"),
            "cli-project",
        )

    def test_resolve_project_name_falls_back_to_script_constant(self):
        self.assertEqual(
            ensemble_wandb.resolve_project_name(None),
            ensemble_wandb.WANDB_PROJECT,
        )

    def test_resolve_project_name_errors_without_input_or_script_value(self):
        with mock.patch.object(ensemble_wandb, "WANDB_PROJECT", ""):
            with self.assertRaisesRegex(ValueError, "W&B project"):
                ensemble_wandb.resolve_project_name(None)

    def test_main_launch_sweep_uses_cli_project(self):
        with mock.patch.object(ensemble_wandb, "launch_sweep_agent") as launch_sweep_agent:
            stdout = io.StringIO()
            with mock.patch.object(
                sys,
                "argv",
                [
                    "torch_gnn_carbon_ensemble_wandb.py",
                    "--launch-sweep",
                    "--project",
                    "cli-project",
                ],
            ), mock.patch("sys.stdout", stdout):
                ensemble_wandb.main()

        ensemble_wandb.wandb.sweep.assert_called_once()
        self.assertEqual(
            ensemble_wandb.wandb.sweep.call_args.kwargs["project"],
            "cli-project",
        )
        self.assertEqual(
            launch_sweep_agent.call_args.kwargs["project"],
            "cli-project",
        )

    def test_main_launch_sweep_falls_back_to_script_project(self):
        with mock.patch.object(ensemble_wandb, "launch_sweep_agent") as launch_sweep_agent:
            with mock.patch.object(ensemble_wandb, "WANDB_PROJECT", "script-project"):
                stdout = io.StringIO()
                with mock.patch.object(
                    sys,
                    "argv",
                    [
                        "torch_gnn_carbon_ensemble_wandb.py",
                        "--launch-sweep",
                    ],
                ), mock.patch("sys.stdout", stdout):
                    ensemble_wandb.main()

        self.assertEqual(
            ensemble_wandb.wandb.sweep.call_args.kwargs["project"],
            "script-project",
        )
        self.assertEqual(
            launch_sweep_agent.call_args.kwargs["project"],
            "script-project",
        )

    def test_main_launch_sweep_reads_project_from_yaml_config(self):
        with mock.patch.object(ensemble_wandb, "launch_sweep_agent"):
            with mock.patch.object(
                ensemble_wandb,
                "load_external_wandb_config",
                return_value=({"project": "yaml-project"}, None),
            ):
                stdout = io.StringIO()
                with mock.patch.object(
                    sys,
                    "argv",
                    [
                        "torch_gnn_carbon_ensemble_wandb.py",
                        "--launch-sweep",
                        "--config-yaml",
                        "config.yaml",
                    ],
                ), mock.patch("sys.stdout", stdout):
                    ensemble_wandb.main()

        self.assertEqual(
            ensemble_wandb.wandb.sweep.call_args.kwargs["project"],
            "yaml-project",
        )

    def test_main_launch_sweep_cli_project_overrides_yaml(self):
        with mock.patch.object(ensemble_wandb, "launch_sweep_agent"):
            with mock.patch.object(
                ensemble_wandb,
                "load_external_wandb_config",
                return_value=({"project": "yaml-project"}, None),
            ):
                stdout = io.StringIO()
                with mock.patch.object(
                    sys,
                    "argv",
                    [
                        "torch_gnn_carbon_ensemble_wandb.py",
                        "--launch-sweep",
                        "--config-yaml",
                        "config.yaml",
                        "--project",
                        "cli-project",
                    ],
                ), mock.patch("sys.stdout", stdout):
                    ensemble_wandb.main()

        self.assertEqual(
            ensemble_wandb.wandb.sweep.call_args.kwargs["project"],
            "cli-project",
        )

    def test_build_sweep_config_broad_uses_reduced_model_sizes(self):
        args = ensemble_wandb.parse_args([])
        args.sweep_profile = "broad"
        sweep = ensemble_wandb.build_sweep_config(args)
        self.assertEqual(
            sweep["parameters"]["hidden_dim"]["values"],
            [32, 48, 64, 80],
        )
        self.assertEqual(
            sweep["parameters"]["num_message_layers"]["values"],
            [1, 2, 3],
        )

    def test_build_sweep_config_architecture_frontier_uses_compact_ranges(self):
        args = ensemble_wandb.parse_args([])
        args.sweep_profile = "architecture-frontier"
        sweep = ensemble_wandb.build_sweep_config(args)
        self.assertEqual(
            sweep["parameters"]["hidden_dim"]["values"],
            [16, 24, 32],
        )
        self.assertEqual(
            sweep["parameters"]["num_message_layers"]["values"],
            [1, 2, 3],
        )

    def test_main_launch_sweep_applies_yaml_sweep_override(self):
        override = {
            "method": "grid",
            "parameters": {
                "hidden_dim": {"values": [24]},
            },
        }
        with mock.patch.object(ensemble_wandb, "launch_sweep_agent"):
            with mock.patch.object(
                ensemble_wandb,
                "load_external_wandb_config",
                return_value=({}, override),
            ):
                stdout = io.StringIO()
                with mock.patch.object(
                    sys,
                    "argv",
                    [
                        "torch_gnn_carbon_ensemble_wandb.py",
                        "--launch-sweep",
                        "--config-yaml",
                        "config.yaml",
                    ],
                ), mock.patch("sys.stdout", stdout):
                    ensemble_wandb.main()

        sweep_payload = ensemble_wandb.wandb.sweep.call_args.kwargs["sweep"]
        self.assertEqual(sweep_payload["method"], "grid")
        self.assertEqual(
            sweep_payload["parameters"]["hidden_dim"]["values"],
            [24],
        )


if __name__ == "__main__":
    unittest.main()
