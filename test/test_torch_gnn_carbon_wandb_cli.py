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
wandb_stub.agent = mock.Mock()
wandb_stub.sweep = mock.Mock()
wandb_stub.init = mock.Mock()
wandb_stub.log = mock.Mock()
wandb_stub.Table = mock.Mock()
wandb_stub.Image = mock.Mock()
wandb_stub.Artifact = mock.Mock()

workflow_stub = types.ModuleType("torch_gnn_carbon_workflow_example")
workflow_stub.CELL_VIOLATION_WEIGHT = 0.0
workflow_stub.COORDINATE_CLIP_VALUE = None
workflow_stub.DEFAULT_MODEL_CONFIG = {
    "hidden_dim": 80,
    "num_message_layers": 2,
    "learning_rate": 5.0e-4,
    "lr_decay_rate": 5.0e-3,
    "smooth_cutoff_width": 0.2,
    "reference_layer_type": 0,
    "component_weight": (1.0, 1.0, 1.0),
}
workflow_stub.DEFAULT_REPLAY_CATEGORY_WEIGHTS = {
    "successful": 1.0,
    "failed": 2.0,
    "unstable": 2.0,
    "high_error": 2.0,
    "difficult": 2.0,
}
workflow_stub.FINGERPRINT_LOSS_WEIGHT = 1.0
workflow_stub.FIXED_LEADING_ATOMS = 0
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
workflow_stub.default_epoch_values = lambda epochs: [epochs]
workflow_stub.default_inverse_step_values = lambda steps: [steps]
workflow_stub.default_step_sizes = lambda step_size: [step_size]
workflow_stub.parse_category_weights = lambda value: {
    key.strip(): float(raw_value.strip())
    for key, raw_value in (
        item.split("=", 1) for item in value.split(",") if item.strip()
    )
}
workflow_stub.parse_float_list = lambda value: [float(item) for item in value.split(",") if item]
workflow_stub.parse_int_list = lambda value: [int(item) for item in value.split(",") if item]
workflow_stub.run_example = mock.Mock()

with mock.patch.dict(
    sys.modules,
    {
        "wandb": wandb_stub,
        "torch_gnn_carbon_workflow_example": workflow_stub,
    },
):
    carbon_wandb = importlib.import_module("torch_gnn_carbon_wandb")


class TestTorchGNNSweepCLI(unittest.TestCase):

    def test_log_series_tables_includes_descriptor_tables(self):
        carbon_wandb.wandb.log.reset_mock()
        carbon_wandb.wandb.Table.reset_mock()
        metrics = {
            "epoch_sweep": [{"epochs": 1, "position_difference": 0.5}],
            "inverse_step_sweep": [
                {"step": 0, "position_difference": 0.5, "fingerprint_mse": 0.1, "fingerprint_l2": 0.2}
            ],
            "step_size_sweep": [{"step_size": 0.02, "position_difference": 0.4}],
            "checkpoint_step_size_sweep": [],
            "checkpoint_step_schedule_sweep": [],
            "rollout": {},
            "descriptor_comparisons": {
                "final": {
                    "component_dimensions": {"2body": 1, "3body": 1, "4body": 1},
                    "component_summaries": {
                        "2body": {
                            "ml_predicted_vs_target": {
                                "mae": 0.1,
                                "rmse": 0.1,
                                "mse": 0.01,
                                "l2": 0.1,
                                "max_abs_error": 0.1,
                            },
                            "true_raffle_vs_target": {
                                "mae": 0.0,
                                "rmse": 0.0,
                                "mse": 0.0,
                                "l2": 0.0,
                                "max_abs_error": 0.0,
                            },
                            "ml_predicted_vs_true_raffle": {
                                "mae": 0.1,
                                "rmse": 0.1,
                                "mse": 0.01,
                                "l2": 0.1,
                                "max_abs_error": 0.1,
                            },
                        },
                        "3body": {
                            "ml_predicted_vs_target": {
                                "mae": 0.2,
                                "rmse": 0.2,
                                "mse": 0.04,
                                "l2": 0.2,
                                "max_abs_error": 0.2,
                            },
                            "true_raffle_vs_target": {
                                "mae": 0.0,
                                "rmse": 0.0,
                                "mse": 0.0,
                                "l2": 0.0,
                                "max_abs_error": 0.0,
                            },
                            "ml_predicted_vs_true_raffle": {
                                "mae": 0.2,
                                "rmse": 0.2,
                                "mse": 0.04,
                                "l2": 0.2,
                                "max_abs_error": 0.2,
                            },
                        },
                        "4body": {
                            "ml_predicted_vs_target": {
                                "mae": 0.3,
                                "rmse": 0.3,
                                "mse": 0.09,
                                "l2": 0.3,
                                "max_abs_error": 0.3,
                            },
                            "true_raffle_vs_target": {
                                "mae": 0.0,
                                "rmse": 0.0,
                                "mse": 0.0,
                                "l2": 0.0,
                                "max_abs_error": 0.0,
                            },
                            "ml_predicted_vs_true_raffle": {
                                "mae": 0.3,
                                "rmse": 0.3,
                                "mse": 0.09,
                                "l2": 0.3,
                                "max_abs_error": 0.3,
                            },
                        },
                    },
                    "components": {
                        "2body": {
                            "target_raffle": [1.0],
                            "ml_predicted": [1.1],
                            "true_raffle": [1.0],
                        },
                        "3body": {
                            "target_raffle": [2.0],
                            "ml_predicted": [2.2],
                            "true_raffle": [2.0],
                        },
                        "4body": {
                            "target_raffle": [3.0],
                            "ml_predicted": [3.3],
                            "true_raffle": [3.0],
                        },
                    },
                }
            },
        }

        carbon_wandb.log_series_tables(metrics)

        carbon_wandb.wandb.log.assert_called_once()
        payload = carbon_wandb.wandb.log.call_args.args[0]
        self.assertIn("tables/final_descriptor_components", payload)
        self.assertIn("tables/final_descriptor_values", payload)

    def test_resolve_existing_sweep_rejects_other_project(self):
        sweep = mock.Mock(
            project="different-project",
            entity="test-entity",
            id="abc123",
            url="https://wandb.ai/test-entity/different-project/sweeps/abc123",
        )
        api = mock.Mock()
        api.sweep.return_value = sweep

        with mock.patch.object(carbon_wandb.wandb, "Api", return_value=api):
            with self.assertRaisesRegex(ValueError, "belongs to project"):
                carbon_wandb.resolve_existing_sweep("abc123")

    def test_main_continues_existing_sweep(self):
        args = carbon_wandb.parse_args(["--sweep-id", "abc123", "--sweep-count", "3"])
        sweep = mock.Mock(
            project=carbon_wandb.WANDB_PROJECT,
            entity="test-entity",
            id="abc123",
            url=(
                "https://wandb.ai/test-entity/"
                f"{carbon_wandb.WANDB_PROJECT}/sweeps/abc123"
            ),
        )
        api = mock.Mock()
        api.sweep.return_value = sweep

        with (
            mock.patch.object(carbon_wandb, "parse_args", return_value=args),
            mock.patch.object(carbon_wandb.wandb, "Api", return_value=api) as api_factory,
            mock.patch.object(carbon_wandb.wandb, "agent") as agent,
            mock.patch.object(carbon_wandb, "execute_run") as execute_run,
            mock.patch("sys.stdout", new_callable=io.StringIO) as stdout,
        ):
            carbon_wandb.main()

            api_factory.assert_called_once_with(overrides={"project": carbon_wandb.WANDB_PROJECT})
            api.sweep.assert_called_once_with("abc123")

            agent.assert_called_once()
            self.assertEqual(agent.call_args.args[0], "abc123")
            self.assertEqual(agent.call_args.kwargs["entity"], "test-entity")
            self.assertEqual(agent.call_args.kwargs["project"], carbon_wandb.WANDB_PROJECT)
            self.assertEqual(agent.call_args.kwargs["count"], 3)

            execute_run.assert_not_called()
            agent.call_args.kwargs["function"]()
            execute_run.assert_called_once()
            self.assertTrue(execute_run.call_args.kwargs["sweep_run"])
            self.assertEqual(
                execute_run.call_args.kwargs["config"]["architecture"],
                carbon_wandb.DEFAULT_ARCHITECTURE,
            )
            self.assertNotIn("sweep_id", execute_run.call_args.kwargs["config"])

            output = stdout.getvalue()
            self.assertIn('"action": "continue"', output)
            self.assertIn('"sweep_id": "abc123"', output)


if __name__ == "__main__":
    unittest.main()
