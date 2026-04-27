import unittest

import numpy as np
from ase.build import bulk

from raffle import minimum_image_displacements, symmetry_aware_rmsd, wrap_atoms_to_unit_cell


class TestStructureMetrics(unittest.TestCase):

    def setUp(self):
        self.original = bulk("C", "diamond", a=3.567, cubic=True)
        self.original.pbc = True

    def test_symmetry_aware_rmsd_ignores_global_translation(self):
        translated = self.original.copy()
        translated.set_positions(translated.get_positions() + np.array([0.1, 0.2, 0.3]))

        self.assertGreater(np.linalg.norm(minimum_image_displacements(self.original, translated)), 0.0)
        self.assertLess(symmetry_aware_rmsd(self.original, translated), 1.0e-10)

    def test_symmetry_aware_rmsd_ignores_global_rotation(self):
        rotated = self.original.copy()
        centroid = self.original.get_positions().mean(axis=0, keepdims=True)
        rotation = np.array(
            [
                [0.0, -1.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0],
            ]
        )
        rotated.set_positions((self.original.get_positions() - centroid) @ rotation + centroid)

        self.assertLess(symmetry_aware_rmsd(self.original, rotated), 1.0e-10)

    def test_wrap_atoms_to_unit_cell_keeps_fractional_positions_bounded(self):
        displaced = self.original.copy()
        displaced.set_positions(displaced.get_positions() + np.array([5.0, -3.5, 7.25]))

        wrapped = wrap_atoms_to_unit_cell(displaced)
        scaled_positions = wrapped.get_scaled_positions(wrap=False)
        self.assertTrue(np.all(scaled_positions >= -1.0e-10))
        self.assertTrue(np.all(scaled_positions < 1.0 + 1.0e-10))


if __name__ == "__main__":
    unittest.main()
