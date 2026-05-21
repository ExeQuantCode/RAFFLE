import importlib
import sys
import types
import unittest
from pathlib import Path
from unittest import mock

import numpy as np
from ase import Atoms


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

with mock.patch.dict(sys.modules, {"raffle": raffle_stub}):
    workflow_common = importlib.import_module("torch_gnn_workflow_common")

build_inverse_design_pair = workflow_common.build_inverse_design_pair
measure_minimum_interatomic_distance = workflow_common.measure_minimum_interatomic_distance
perturb_structure = workflow_common.perturb_structure


class TestTorchGNNWorkflowCommon(unittest.TestCase):

    def test_perturb_structure_is_seeded_and_preserves_fixed_atoms(self):
        atoms = Atoms(
            "C2",
            positions=[[0.1, 0.2, 0.3], [1.2, 1.1, 1.0]],
            cell=np.eye(3) * 3.0,
            pbc=True,
        )
        fixed_atoms = np.asarray([True, False], dtype=bool)
        settings = {
            "min_displacement": 0.0,
            "max_displacement": 0.25,
            "minimum_interatomic_distance": 0.5,
            "max_resamples": 32,
        }

        perturbed_a = perturb_structure(
            atoms,
            seed=17,
            fixed_atoms=fixed_atoms,
            perturbation_settings=settings,
        )
        perturbed_b = perturb_structure(
            atoms,
            seed=17,
            fixed_atoms=fixed_atoms,
            perturbation_settings=settings,
        )

        np.testing.assert_allclose(
            perturbed_a.get_positions(),
            perturbed_b.get_positions(),
        )
        np.testing.assert_allclose(
            perturbed_a.get_positions()[0],
            atoms.get_positions()[0],
        )
        self.assertGreaterEqual(
            measure_minimum_interatomic_distance(perturbed_a),
            settings["minimum_interatomic_distance"],
        )

    def test_build_inverse_design_pair_wraps_pbc_and_keeps_fixed_atoms_aligned(self):
        atoms = Atoms(
            "C3",
            positions=[[1.9, 0.1, 0.1], [0.6, 0.6, 0.6], [1.2, 1.4, 1.6]],
            cell=np.eye(3) * 2.0,
            pbc=True,
        )
        settings = {
            "min_displacement": 0.0,
            "max_displacement": 0.3,
            "minimum_interatomic_distance": 0.5,
            "max_resamples": 64,
        }

        target_atoms, input_atoms, fixed_atoms, metadata = build_inverse_design_pair(
            [atoms],
            fixed_leading_atoms=1,
            seed=23,
            perturbation_settings=settings,
        )

        self.assertEqual(metadata["structure_index"], 0)
        self.assertEqual(metadata["fixed_atom_indices"], [0])
        self.assertTrue(bool(fixed_atoms[0]))
        self.assertTrue(np.all(target_atoms.get_scaled_positions(wrap=True) >= 0.0))
        self.assertTrue(np.all(target_atoms.get_scaled_positions(wrap=True) < 1.0))
        self.assertTrue(np.all(input_atoms.get_scaled_positions(wrap=True) >= 0.0))
        self.assertTrue(np.all(input_atoms.get_scaled_positions(wrap=True) < 1.0))
        np.testing.assert_allclose(
            target_atoms.get_positions()[fixed_atoms],
            input_atoms.get_positions()[fixed_atoms],
        )
        self.assertGreaterEqual(
            measure_minimum_interatomic_distance(target_atoms),
            settings["minimum_interatomic_distance"],
        )
        self.assertGreaterEqual(
            measure_minimum_interatomic_distance(input_atoms),
            settings["minimum_interatomic_distance"],
        )


if __name__ == "__main__":
    unittest.main()
