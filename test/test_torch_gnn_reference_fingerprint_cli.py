import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
from ase.build import bulk
from ase.io import write

from raffle import TorchGNNFingerprint


REPO_ROOT = Path(__file__).resolve().parents[1]
PYTHON_PKG_DIR = REPO_ROOT / "example" / "python_pkg"
if str(PYTHON_PKG_DIR) not in sys.path:
    sys.path.insert(0, str(PYTHON_PKG_DIR))

from torch_gnn_workflow_common import (  # noqa: E402
    DEFAULT_MODEL_CONFIG,
    create_model,
    load_target_fingerprint,
    save_model_checkpoint,
)


@unittest.skipUnless(TorchGNNFingerprint is not None, "PyTorch multigraph fingerprint is unavailable")
class TestTorchGNNReferenceFingerprintCLI(unittest.TestCase):

    def test_cli_writes_npy_loadable_by_inverse_design(self):
        with tempfile.TemporaryDirectory() as tmp_dir_name:
            tmp_dir = Path(tmp_dir_name)
            structure_path = tmp_dir / "reference_structure.xyz"
            checkpoint_path = tmp_dir / "torch_gnn_model_checkpoint.pt"
            output_path = tmp_dir / "target_fingerprint.npy"

            structure = bulk("C", "diamond", a=3.567, cubic=True)
            structure.pbc = True
            write(structure_path, structure, format="extxyz")

            model = create_model(
                seed=42,
                species_list=["C"],
                model_config=DEFAULT_MODEL_CONFIG,
            )
            expected_fingerprint = np.asarray(
                model.compute_reference_fingerprint(structure),
                dtype=np.float32,
            ).reshape(-1)
            save_model_checkpoint(
                model=model,
                checkpoint_path=checkpoint_path,
                species_list=["C"],
                model_config=DEFAULT_MODEL_CONFIG,
                training_config={"seed": 42},
                training_history=[0.0],
            )

            script_path = PYTHON_PKG_DIR / "torch_gnn_compute_reference_fingerprint.py"
            result = subprocess.run(
                [
                    sys.executable,
                    str(script_path),
                    "--input-structure",
                    str(structure_path),
                    "--model-checkpoint",
                    str(checkpoint_path),
                    "--output-path",
                    str(output_path),
                ],
                cwd=PYTHON_PKG_DIR,
                capture_output=True,
                text=True,
            )

            if result.returncode != 0:
                self.fail(
                    "CLI failed with non-zero exit code\n"
                    f"stdout:\n{result.stdout}\n"
                    f"stderr:\n{result.stderr}"
                )

            self.assertTrue(output_path.exists())
            saved_fingerprint = load_target_fingerprint(output_path)
            self.assertEqual(saved_fingerprint.dtype, np.float32)
            self.assertEqual(saved_fingerprint.size, int(model.fingerprint_dim))
            np.testing.assert_allclose(saved_fingerprint, expected_fingerprint)


if __name__ == "__main__":
    unittest.main()
