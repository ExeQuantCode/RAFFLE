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

raffle_stub = types.ModuleType("raffle")
raffle_stub.__path__ = [str(REPO_ROOT / "src" / "raffle")]
raffle_stub.TorchGNNFingerprint = None
raffle_stub.minimum_image_displacements = mock.Mock()
raffle_stub.structure_similarity_rmsd = mock.Mock()
raffle_stub.symmetry_aware_displacements = mock.Mock()

torch_gnn_rollout_stub = types.ModuleType("torch_gnn_rollout")
torch_gnn_rollout_stub.PrioritizedReplayBuffer = object
torch_gnn_rollout_stub.ReplaySample = object
torch_gnn_rollout_stub.classify_rollout_step = mock.Mock(return_value="successful")

sys.modules.pop("torch_gnn_workflow_common", None)
sys.modules.pop("torch_gnn_carbon_workflow_example", None)

with mock.patch.dict(
    sys.modules,
    {
        "raffle": raffle_stub,
        "torch_gnn_rollout": torch_gnn_rollout_stub,
    },
):
    workflow_example = importlib.import_module("torch_gnn_carbon_workflow_example")


class TestTorchGNNCarbonValidationCases(unittest.TestCase):

    def test_build_final_validation_cases_load_file_backed_inputs(self):
        cases_a = workflow_example.build_final_validation_cases(
            perturbation_settings={
                "min_displacement": 0.0,
                "max_displacement": 0.2,
                "minimum_interatomic_distance": 0.8,
                "max_resamples": 64,
            },
            seed=17,
        )
        cases_b = workflow_example.build_final_validation_cases(
            perturbation_settings={
                "min_displacement": 0.0,
                "max_displacement": 10.0,
                "minimum_interatomic_distance": 0.1,
                "max_resamples": 1,
            },
            seed=101,
        )

        self.assertEqual([case["name"] for case in cases_a], ["diamond", "graphite"])
        self.assertEqual([len(case["reference_atoms"]) for case in cases_a], [8, 4])
        for case_a, case_b in zip(cases_a, cases_b, strict=True):
            reference_atoms = case_a["reference_atoms"]
            perturbed_atoms = case_a["perturbed_atoms"]
            self.assertEqual(case_a["fixed_atom_indices"], [0])
            self.assertTrue(bool(case_a["fixed_atoms"][0]))
            np.testing.assert_allclose(
                case_a["reference_atoms"].get_positions(),
                case_b["reference_atoms"].get_positions(),
            )
            np.testing.assert_allclose(
                case_a["perturbed_atoms"].get_positions(),
                case_b["perturbed_atoms"].get_positions(),
            )
            self.assertIn(f"{case_a['name']}.xyz", case_a["input_files"]["reference"])
            self.assertIn(
                f"perturbed_{case_a['name']}.xyz",
                case_a["input_files"]["perturbed"],
            )

        diamond_case = cases_a[0]
        graphite_case = cases_a[1]
        np.testing.assert_allclose(
            diamond_case["reference_atoms"].get_positions()[0],
            [0.0, 0.0, 0.0],
        )
        np.testing.assert_allclose(
            diamond_case["perturbed_atoms"].get_positions()[0],
            [0.5, 0.0, 0.0],
        )
        np.testing.assert_allclose(
            graphite_case["reference_atoms"].get_positions()[0],
            [0.0, 0.0, 1.95076825],
        )
        np.testing.assert_allclose(
            graphite_case["perturbed_atoms"].get_positions()[0],
            [0.5, 0.0, 2.45076825],
        )

        for case in cases_a:
            reference_atoms = case["reference_atoms"]
            perturbed_atoms = case["perturbed_atoms"]
            np.testing.assert_allclose(
                np.asarray(reference_atoms.cell.array, dtype=float),
                np.asarray(perturbed_atoms.cell.array, dtype=float),
            )
            self.assertFalse(
                np.allclose(reference_atoms.get_positions(), perturbed_atoms.get_positions())
            )

    def test_fixed_atom_aligned_rmsd_reports_diamond_improvement(self):
        diamond_case = workflow_example.build_final_validation_cases()[0]
        reference_atoms = diamond_case["reference_atoms"]
        perturbed_atoms = diamond_case["perturbed_atoms"]
        fixed_atoms = diamond_case["fixed_atoms"]

        improved_atoms = workflow_example.build_fixed_atom_aligned_reference(
            reference_atoms,
            perturbed_atoms,
            fixed_atoms,
        )
        improved_atoms.set_positions(
            improved_atoms.get_positions()
            + np.asarray(
                [
                    [0.0, 0.0, 0.0],
                    [0.01, -0.01, 0.0],
                    [-0.01, 0.01, 0.0],
                    [0.01, 0.0, -0.01],
                    [0.0, 0.01, 0.01],
                    [-0.01, 0.0, 0.01],
                    [0.0, -0.01, 0.0],
                    [0.01, 0.01, -0.01],
                ],
                dtype=float,
            )
        )

        with mock.patch.object(
            workflow_example,
            "symmetry_aware_displacements",
            side_effect=lambda reference, candidate, allow_rotation=False: (
                candidate.get_positions() - reference.get_positions()
            ),
        ):
            initial_rmsd = workflow_example.compute_fixed_atom_aligned_rmsd(
                reference_atoms,
                perturbed_atoms,
                fixed_atoms,
            )
            final_rmsd = workflow_example.compute_fixed_atom_aligned_rmsd(
                reference_atoms,
                improved_atoms,
                fixed_atoms,
            )

        self.assertGreater(initial_rmsd, final_rmsd)

    def test_print_final_validation_result_matches_requested_format(self):
        stdout = io.StringIO()
        comparison_reference_path = Path("/tmp/diamond_comparison_reference.xyz")
        initial_structure_path = Path("/tmp/diamond_initial.xyz")
        final_structure_path = Path("/tmp/diamond_final.xyz")
        with mock.patch("sys.stdout", stdout):
            workflow_example.print_final_validation_result(
                "diamond",
                initial_rmsd=0.195154,
                final_rmsd=0.050326,
                comparison_reference_path=comparison_reference_path,
                initial_structure_path=initial_structure_path,
                final_structure_path=final_structure_path,
            )

        self.assertEqual(
            stdout.getvalue().splitlines(),
            [
                "Final model validation: diamond",
                f"#sym:structure_similarity_rmsd {comparison_reference_path} {initial_structure_path}",
                "Initial symmetry-aware RMSD: 0.195154 A",
                f"#sym:structure_similarity_rmsd {comparison_reference_path} {final_structure_path}",
                "Final symmetry-aware RMSD:   0.050326 A",
                "RMSD reduction:              74.21%",
            ],
        )

    def test_evaluate_final_validation_cases_serialises_file_metadata(self):
        class DummyModel:

            def compute_reference_fingerprint(self, atoms):
                del atoms
                return np.zeros(1, dtype=float)

        fake_inverse_path = {
            "traj_file": "/tmp/final-validation.traj",
            "step_structure_dir": "/tmp/final-validation-steps",
            "step_structure_files": [],
            "steps": [],
        }

        with tempfile.TemporaryDirectory() as tmp_dir:
            with (
                mock.patch.object(
                    workflow_example,
                    "inverse_design_trace",
                    side_effect=lambda **kwargs: (
                        kwargs["perturbed"].copy(),
                        [],
                        [],
                    ),
                ),
                mock.patch.object(
                    workflow_example,
                    "compute_fixed_atom_aligned_rmsd",
                    side_effect=[0.2, 0.05, 0.1, 0.04],
                ),
                mock.patch.object(
                    workflow_example,
                    "save_inverse_design_path",
                    return_value=fake_inverse_path,
                ),
                mock.patch.object(workflow_example, "write"),
            ):
                results, summary = workflow_example.evaluate_final_validation_cases(
                    model=DummyModel(),
                    output_dir=Path(tmp_dir),
                    seed=42,
                    inverse_steps=30,
                    inverse_step_size=0.02,
                    inverse_design_options={},
                )

        self.assertEqual([entry["name"] for entry in results], ["diamond", "graphite"])
        self.assertTrue(results[0]["input_files"]["reference"].endswith("diamond.xyz"))
        self.assertTrue(
            results[0]["input_files"]["perturbed"].endswith("perturbed_diamond.xyz")
        )
        self.assertEqual(results[0]["comparison_metric"], "#sym:structure_similarity_rmsd")
        self.assertIn("comparison_reference", results[0]["output_files"])
        self.assertEqual(results[0]["inverse_design_path"], fake_inverse_path)
        self.assertEqual(summary["case_count"], 2)


if __name__ == "__main__":
    unittest.main()
