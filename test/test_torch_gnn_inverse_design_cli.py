import importlib
import io
import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest import mock

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[1]
PYTHON_PKG_DIR = REPO_ROOT / "example" / "python_pkg"
if str(PYTHON_PKG_DIR) not in sys.path:
    sys.path.insert(0, str(PYTHON_PKG_DIR))

ase_stub = types.ModuleType("ase")
ase_io_stub = types.ModuleType("ase.io")
ase_io_stub.write = mock.Mock()
ase_stub.io = ase_io_stub

workflow_common_stub = types.ModuleType("torch_gnn_workflow_common")
workflow_common_stub.build_inverse_design_options = mock.Mock(return_value={"option": True})
workflow_common_stub.compute_inverse_design_metrics = mock.Mock(
    return_value={"initial_fingerprint_mse": 1.0, "final_fingerprint_mse": 0.5}
)
workflow_common_stub.load_model_from_checkpoint = mock.Mock()
workflow_common_stub.load_target_fingerprint = mock.Mock()
workflow_common_stub.print_position_differences = mock.Mock()
workflow_common_stub.read_single_structure = mock.Mock()
workflow_common_stub.save_descriptor_comparison_report = mock.Mock(
    return_value={"report_file": "descriptor_report.json", "plot_file": "descriptor_plot.png"}
)
workflow_common_stub.structures_have_matching_atom_count = mock.Mock(return_value=True)
workflow_common_stub.write_inverse_design_log = mock.Mock()
workflow_common_stub.write_json = mock.Mock()
workflow_common_stub.write_structure = mock.Mock()

with mock.patch.dict(
    sys.modules,
    {
        "ase": ase_stub,
        "ase.io": ase_io_stub,
        "torch_gnn_workflow_common": workflow_common_stub,
    },
):
    inverse_design_cli = importlib.import_module("torch_gnn_inverse_design")


class TestTorchGNNInverseDesignCLI(unittest.TestCase):

    class _FakeAtoms:
        def __init__(self, label: str):
            self.label = label
            self.info: dict[str, object] = {}

        def copy(self):
            copied = TestTorchGNNInverseDesignCLI._FakeAtoms(self.label)
            copied.info = dict(self.info)
            return copied

    def setUp(self) -> None:
        self.model = mock.Mock()
        self.model.fingerprint_dim = 3
        self.model.fingerprint_dim_2body = 2
        self.model.compute_reference_fingerprint.return_value = [0.2, 0.4, 0.6]
        self.model.predict_components.return_value = (
            np.asarray([0.1, 0.3], dtype=np.float32),
            None,
            None,
        )

        def inverse_design_side_effect(**kwargs):
            observer = kwargs.get("step_observer")
            if observer is not None:
                observer(
                    {
                        "atoms": self._FakeAtoms("initial"),
                        "restart_index": 0,
                        "num_restarts": 1,
                        "step": 0,
                        "num_steps": 4,
                        "is_initial_state": True,
                    }
                )
                observer(
                    {
                        "atoms": self._FakeAtoms("final"),
                        "restart_index": 0,
                        "num_restarts": 1,
                        "step": 4,
                        "num_steps": 4,
                        "is_initial_state": False,
                    }
                )
            return "optimised-atoms"

        self.model.inverse_design.side_effect = inverse_design_side_effect

        workflow_common_stub.build_inverse_design_options.reset_mock(return_value=True)
        workflow_common_stub.build_inverse_design_options.return_value = {"option": True}
        workflow_common_stub.compute_inverse_design_metrics.reset_mock(return_value=True)
        workflow_common_stub.compute_inverse_design_metrics.return_value = {
            "initial_fingerprint_mse": 1.0,
            "final_fingerprint_mse": 0.5,
        }
        workflow_common_stub.load_model_from_checkpoint.reset_mock(return_value=True)
        workflow_common_stub.load_model_from_checkpoint.return_value = (
            self.model,
            {"training_config": {"seed": 7}},
        )
        workflow_common_stub.load_target_fingerprint.reset_mock()
        workflow_common_stub.print_position_differences.reset_mock()
        workflow_common_stub.read_single_structure.reset_mock()
        self.input_atoms = [object(), object(), object(), object(), object()]
        workflow_common_stub.read_single_structure.side_effect = [
            self.input_atoms,
            "target-atoms",
        ]
        workflow_common_stub.save_descriptor_comparison_report.reset_mock(return_value=True)
        workflow_common_stub.save_descriptor_comparison_report.return_value = {
            "report_file": "descriptor_report.json",
            "plot_file": "descriptor_plot.png",
        }
        workflow_common_stub.write_inverse_design_log.reset_mock()
        workflow_common_stub.write_json.reset_mock()
        workflow_common_stub.write_structure.reset_mock()
        ase_io_stub.write.reset_mock()

    def test_parse_args_accepts_reference_structure_alias(self):
        args = inverse_design_cli.parse_args(
            [
                "--reference-structure",
                "reference.xyz",
                "--target-structure",
                "target.xyz",
                "--model-checkpoint",
                "checkpoint.pt",
            ]
        )

        self.assertEqual(args.input_structure, Path("reference.xyz"))
        self.assertEqual(args.target_structure, Path("target.xyz"))
        self.assertIsNone(args.target_fingerprint)

    def test_parse_args_accepts_fixed_atom_indices(self):
        args = inverse_design_cli.parse_args(
            [
                "--reference-structure",
                "reference.xyz",
                "--target-structure",
                "target.xyz",
                "--model-checkpoint",
                "checkpoint.pt",
                "--fixed-atoms",
                "1,3,4",
            ]
        )

        self.assertEqual(args.fixed_atoms, [1, 3, 4])

    def test_main_uses_target_structure_when_target_fingerprint_is_omitted(self):
        with tempfile.TemporaryDirectory() as tmp_dir_name:
            tmp_dir = Path(tmp_dir_name)
            stdout = io.StringIO()

            with mock.patch("sys.stdout", stdout):
                inverse_design_cli.main(
                    [
                        "--reference-structure",
                        str(tmp_dir / "reference.xyz"),
                        "--target-structure",
                        str(tmp_dir / "target.xyz"),
                        "--model-checkpoint",
                        str(tmp_dir / "checkpoint.pt"),
                        "--fixed-leading-atoms",
                        "1",
                        "--fixed-atoms",
                        "3,4",
                        "--output-dir",
                        str(tmp_dir / "outputs"),
                    ]
                )

        workflow_common_stub.load_target_fingerprint.assert_not_called()
        self.model.compute_reference_fingerprint.assert_called_once_with("target-atoms")
        inverse_call = self.model.inverse_design.call_args.kwargs
        np.testing.assert_allclose(
            inverse_call["target_fingerprint"],
            np.asarray([0.2, 0.4, 0.6], dtype=np.float32),
        )
        np.testing.assert_array_equal(
            inverse_call["fixed_atoms"],
            np.asarray([True, False, False, True, True]),
        )

        metrics_payload = workflow_common_stub.write_json.call_args.args[1]
        self.assertEqual(metrics_payload["target_fingerprint_file"], None)
        self.assertEqual(metrics_payload["target_fingerprint_source"]["type"], "target_structure")
        self.assertTrue(metrics_payload["target_fingerprint_source"]["path"].endswith("target.xyz"))
        self.assertEqual(metrics_payload["inverse_design_config"]["fixed_atoms"], [0, 3, 4])
        self.assertTrue(metrics_payload["output_files"]["optimisation_traj"].endswith("torch_gnn_inverse_design_path.traj"))
        self.assertEqual(len(metrics_payload["optimisation_history"]), 2)
        ase_io_stub.write.assert_called_once()
        write_path, write_frames = ase_io_stub.write.call_args.args
        self.assertTrue(str(write_path).endswith("torch_gnn_inverse_design_path.traj"))
        self.assertEqual(len(write_frames), 2)
        self.assertIn("torch_gnn_inverse_design_final.xyz", stdout.getvalue())

    def test_main_allows_disabling_optimisation_traj_export(self):
        with tempfile.TemporaryDirectory() as tmp_dir_name:
            tmp_dir = Path(tmp_dir_name)

            inverse_design_cli.main(
                [
                    "--reference-structure",
                    str(tmp_dir / "reference.xyz"),
                    "--target-structure",
                    str(tmp_dir / "target.xyz"),
                    "--model-checkpoint",
                    str(tmp_dir / "checkpoint.pt"),
                    "--output-dir",
                    str(tmp_dir / "outputs"),
                    "--no-save-optimisation-traj",
                ]
            )

        metrics_payload = workflow_common_stub.write_json.call_args.args[1]
        self.assertNotIn("optimisation_traj", metrics_payload["output_files"])
        self.assertNotIn("optimisation_history", metrics_payload)
        ase_io_stub.write.assert_not_called()


if __name__ == "__main__":
    unittest.main()
