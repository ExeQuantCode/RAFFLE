import io
import importlib
import json
import sys
import types
import unittest
from pathlib import Path
from unittest import mock

import numpy as np


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
wandb_stub.__version__ = "0.0.0-test"


class DummyCell:

    def __init__(self, values):
        self.array = np.asarray(values, dtype=float)


class DummyAtoms:

    def __init__(self, symbols=None, positions=None, cell=None, pbc=None):
        self._symbols = list(symbols or ["C", "C"])
        self._positions = np.asarray(
            [[0.0, 0.0, 0.0], [0.25, 0.25, 0.25]] if positions is None else positions,
            dtype=float,
        )
        self.cell = DummyCell(np.eye(3) if cell is None else cell)
        self.pbc = list([True, True, True] if pbc is None else pbc)

    def get_chemical_symbols(self):
        return list(self._symbols)

    def get_positions(self):
        return np.asarray(self._positions, dtype=float)

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
workflow_stub.select_carbon_structures = mock.Mock(
    side_effect=lambda structures, requested_count: (
        list(structures)
        if int(requested_count) <= 0 or int(requested_count) >= len(structures)
        else list(structures[: int(requested_count)])
    )
)
workflow_stub.run_example = mock.Mock()

workflow_common_stub = types.ModuleType("torch_gnn_workflow_common")
workflow_common_stub.build_inverse_design_pair = mock.Mock(
    return_value=(
        DummyAtoms(),
        DummyAtoms(positions=[[0.1, 0.0, 0.0], [0.25, 0.35, 0.25]]),
        np.asarray([True, False], dtype=bool),
        {
            "structure_index": 1,
            "perturbation_settings": {
                "min_displacement": 0.0,
                "max_displacement": 1.0,
                "minimum_interatomic_distance": 0.8,
                "max_resamples": 64,
            },
        },
    )
)
workflow_common_stub.load_structures = mock.Mock(return_value=[DummyAtoms(), DummyAtoms()])

with mock.patch.dict(
    sys.modules,
    {
        "wandb": wandb_stub,
        "torch_gnn_carbon_workflow_example": workflow_stub,
        "torch_gnn_workflow_common": workflow_common_stub,
    },
):
    carbon_wandb = importlib.import_module("torch_gnn_carbon_wandb")


class TestTorchGNNSweepCLI(unittest.TestCase):

    def setUp(self) -> None:
        workflow_stub.run_example.reset_mock()
        workflow_stub.select_carbon_structures.reset_mock(side_effect=True)
        workflow_common_stub.build_inverse_design_pair.reset_mock(return_value=True)
        workflow_common_stub.load_structures.reset_mock(return_value=True)
        workflow_common_stub.build_inverse_design_pair.return_value = (
            DummyAtoms(),
            DummyAtoms(positions=[[0.1, 0.0, 0.0], [0.25, 0.35, 0.25]]),
            np.asarray([True, False], dtype=bool),
            {
                "structure_index": 1,
                "perturbation_settings": {
                    "min_displacement": 0.0,
                    "max_displacement": 1.0,
                    "minimum_interatomic_distance": 0.8,
                    "max_resamples": 64,
                },
            },
        )
        workflow_common_stub.load_structures.return_value = [DummyAtoms(), DummyAtoms()]
        carbon_wandb.wandb.log.reset_mock()
        carbon_wandb.wandb.Artifact.reset_mock()

    def test_execute_run_logs_resolved_reproducibility_config(self):
        class DummyConfig(dict):

            def update(self, payload=None, allow_val_change=False, **kwargs):
                del allow_val_change
                if payload:
                    super().update(payload)
                if kwargs:
                    super().update(kwargs)

        class DummyRun:

            def __init__(self, config):
                self.id = "run-123"
                self.name = "dummy-run"
                self.config = DummyConfig(config)
                self.summary = {}
                self.log_artifact = mock.Mock()

        run = DummyRun(
            {
                "project": "custom-project",
                "architecture": carbon_wandb.DEFAULT_ARCHITECTURE,
                "variant_tags": "baseline",
                "run_name": "",
                "carbon_count": -1,
                "augmented_count": 2,
                "epochs": 12,
                "batch_size": 4,
                "inverse_steps": 40,
                "inverse_step_size": 0.02,
                "inverse_step_values": "",
                "step_size_values": "",
                "fixed_leading_atoms": 1,
                "fingerprint_loss_weight": 0.75,
                "target_vertex_weight": 0.0,
                "target_position_weight": 0.0,
                "inverse_lr_decay_rate": 0.0,
                "repulsion_weight": 5.0,
                "minimum_distance_scale": 0.8,
                "cell_violation_weight": 0.01,
                "coordinate_clip_value": 0.5,
                "rollout_stages": 2,
                "rollout_epochs_per_stage": 1,
                "rollout_step_stride": 5,
                "replay_buffer_capacity": 8,
                "replay_sample_size": 4,
                "replay_category_weights": "successful=1.0,failed=3.0",
                "rollout_drift_threshold": 1.0e-4,
                "rollout_high_error_threshold": 1.0e-3,
                "rollout_instability_threshold": 1.0e-5,
                "hidden_dim": 96,
                "num_message_layers": 3,
                "learning_rate": 2.5e-4,
                "model_lr_decay_rate": 1.0e-3,
                "smooth_cutoff_width": 0.15,
                "reference_layer_type": 0,
                "component_weight_2body": 5.0,
                "component_weight_3body": 1.5,
                "component_weight_4body": 0.5,
                "output_dir": "build/wandb_inverse_design",
                "skip_checkpoint_step_size_sweep": False,
                "skip_checkpoint_step_schedule_sweep": True,
                "seed": 7,
            }
        )
        context_manager = mock.MagicMock()
        context_manager.__enter__.return_value = run
        context_manager.__exit__.return_value = False
        carbon_wandb.wandb.init.return_value = context_manager
        carbon_wandb.wandb.Artifact.return_value = mock.Mock()

        workflow_stub.run_example.return_value = {
            "epoch_sweep": [],
            "inverse_step_sweep": [],
            "step_size_sweep": [],
            "checkpoint_step_size_sweep": [],
            "checkpoint_step_schedule_sweep": [],
            "rollout": {
                "enabled": False,
                "replay_buffer": {"size": 0, "mean_fingerprint_drift_mse": 0.0},
                "stages": [],
            },
            "initial_fingerprint_mse": 1.0,
            "final_fingerprint_mse": 0.2,
            "initial_rmsd": 1.5,
            "final_rmsd": 0.5,
            "rmsd_reduction_fraction": 2.0 / 3.0,
            "convergence_summary": {
                "best_step": 10,
                "tail_position_range": 0.1,
                "tail_fingerprint_range": 0.01,
                "fingerprint_nonincreasing_fraction": 1.0,
                "position_nonincreasing_fraction": 1.0,
                "converges_to_best_within_5pct": True,
            },
            "descriptor_comparisons": {},
            "architecture_name": carbon_wandb.DEFAULT_ARCHITECTURE,
            "output_files": {
                "plot": "/tmp/plot.png",
                "checkpoint": "/tmp/checkpoint.pt",
                "target_fingerprint": "/tmp/target.npy",
                "reference_structure": "/tmp/reference.xyz",
            },
        }

        carbon_wandb.execute_run(config=dict(run.config), sweep_run=False)

        self.assertEqual(
            carbon_wandb.wandb.init.call_args.kwargs["project"],
            "custom-project",
        )
        resolved_payload = json.loads(run.config[carbon_wandb.RESOLVED_WORKFLOW_CONFIG_KEY])
        metadata_payload = json.loads(run.config[carbon_wandb.REPRODUCIBILITY_METADATA_KEY])
        self.assertNotIn("reference_structure_weight", resolved_payload)
        self.assertEqual(resolved_payload["inverse_step_values"], [40])
        self.assertEqual(resolved_payload["step_size_values"], [0.02])
        self.assertEqual(
            resolved_payload["ignored_deprecated_config_fields"],
            {},
        )
        self.assertEqual(
            resolved_payload["reference_structure"]["symbols"],
            ["C", "C"],
        )
        self.assertEqual(
            resolved_payload["input_structure"]["symbols"],
            ["C", "C"],
        )
        self.assertEqual(resolved_payload["evaluation_pair"]["structure_index"], 1)
        self.assertEqual(resolved_payload["perturbation_settings"]["max_displacement"], 1.0)
        self.assertEqual(metadata_payload["carbon_dataset"]["relative_path"], "example/data/carbon.xyz")
        self.assertIn("resolved_workflow_config", run.summary)
        self.assertIn("reproducibility_metadata", run.summary)

        run_kwargs = workflow_stub.run_example.call_args.kwargs
        self.assertNotIn("reference_structure_weight", run_kwargs)
        self.assertEqual(run_kwargs["augmented_count"], 2)
        self.assertEqual(run_kwargs["inverse_step_values"], [40])
        self.assertEqual(run_kwargs["step_size_values"], [0.02])
        self.assertFalse(run_kwargs["enable_checkpoint_step_schedule_sweep"])
        self.assertEqual(run_kwargs["model_config"]["component_weight"], (5.0, 1.5, 0.5))

    def test_build_sweep_config_uses_epochs_hyperparameter(self):
        args = carbon_wandb.parse_args(["--sweep-profile", "plan-attention-refine"])

        sweep_config = carbon_wandb.build_sweep_config(args)

        self.assertEqual(
            sweep_config["parameters"]["epochs"]["values"],
            carbon_wandb.SWEEP_EPOCH_COUNTS,
        )
        self.assertNotIn("epoch_values", sweep_config["parameters"])
        self.assertNotIn("inverse_restarts", sweep_config["parameters"])
        self.assertNotIn("inverse_restart_noise_scale", sweep_config["parameters"])

    def test_parse_args_does_not_accept_epoch_values(self):
        with self.assertRaises(SystemExit):
            carbon_wandb.parse_args(["--epoch-values", "0,10,75"])

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
            "final_validation_cases": [
                {
                    "name": "diamond",
                    "atom_count": 8,
                    "fixed_atom_indices": [0],
                    "initial_rmsd": 0.2,
                    "final_rmsd": 0.05,
                    "rmsd_reduction_fraction": 0.75,
                    "inverse_steps": 30,
                    "step_size": 0.02,
                },
                {
                    "name": "graphite",
                    "atom_count": 4,
                    "fixed_atom_indices": [0],
                    "initial_rmsd": 0.1,
                    "final_rmsd": 0.04,
                    "rmsd_reduction_fraction": 0.6,
                    "inverse_steps": 30,
                    "step_size": 0.02,
                },
            ],
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
        self.assertIn("tables/final_model_validation", payload)

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
        args = carbon_wandb.parse_args(
            [
                "--project",
                "custom-project",
                "--sweep-id",
                "abc123",
                "--sweep-count",
                "3",
            ]
        )
        sweep = mock.Mock(
            project="custom-project",
            entity="test-entity",
            id="abc123",
            url=(
                "https://wandb.ai/test-entity/"
                "custom-project/sweeps/abc123"
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

            api_factory.assert_called_once_with(overrides={"project": "custom-project"})
            api.sweep.assert_called_once_with("abc123")

            agent.assert_called_once()
            self.assertEqual(agent.call_args.args[0], "abc123")
            self.assertEqual(agent.call_args.kwargs["entity"], "test-entity")
            self.assertEqual(agent.call_args.kwargs["project"], "custom-project")
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
