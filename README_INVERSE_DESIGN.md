
There is an inverse-design example in `example/python_pkg/` that demonstrates how a neural network can be trained to learn RAFFLE descriptors and then used to perform inverse design as a surrogate atomic relaxation that searches for structures matching a target RAFFLE fingerprint.

The current workflow uses PyTorch, but the code remains modular enough to swap in other backends such as ATHENA once an architecture is identified that performs well on this problem.

The original all-in-one script is `example/python_pkg/torch_gnn_carbon_workflow_example.py`. The same workflow is now also split into two standalone scripts:

- `example/python_pkg/torch_gnn_train_model.py`: build, train, and save a single model checkpoint.
- `example/python_pkg/torch_gnn_inverse_design.py`: load a saved checkpoint and run inverse design against a supplied structure and target fingerprint.

## Environment setup

Use the repository version of `raffle`, not an older published wheel, before running the new scripts:

```bash
python -m pip install -e .
```

This was the required setup during validation because the older installed package in the default Conda environment did not expose `TorchGNNFingerprint`.

## Training script

The training entrypoint accepts either CLI flags or a JSON config file. In both cases it writes the effective run configuration back out to the output directory so the exact training settings are preserved with the checkpoint.

Example flag-based run:

```bash
python example/python_pkg/torch_gnn_train_model.py \
	--training-structure-limit 100 \
	--hidden-dim 128 \
	--num-message-layers 4 \
	--batch-size 8 \
	--reference-layer-type 1 \
	--component-weight-2body 3 \
	--component-weight-3body 1 \
	--component-weight-4body 0 \
	--epochs 100 \
	--learning-rate 0.0005 \
	--model-lr-decay-rate 0.005 \
	--augmented-count 16 \
	--seed 11 \
	--output-dir build/torch_gnn_carbon_train
```

Example JSON config:

```json
{
	"training_structures": "example/data/carbon.xyz",
	"training_structure_limit": 0,
	"augmented_count": 16,
	"augmentation_noise_scale": 0.04,
	"epochs": 100,
	"batch_size": 8,
	"seed": 11,
	"output_dir": "build/torch_gnn_carbon_train",
	"model_config": {
		"architecture": "residual",
		"hidden_dim": 128,
		"num_message_layers": 4,
		"learning_rate": 0.0005,
		"lr_decay_rate": 0.005,
		"smooth_cutoff_width": 0.2,
		"reference_layer_type": 1,
		"component_weight": [3.0, 1.0, 0.0]
	}
}
```

Run from config:

```bash
python example/python_pkg/torch_gnn_train_model.py --config path/to/train_config.json
```

The training script writes these outputs into `--output-dir`:

- `torch_gnn_model_checkpoint.pt`: model weights plus the constructor and training metadata needed to reload the model.
- `torch_gnn_target_fingerprint.npy`: analytical target fingerprint for the reference structure used during training.
- `torch_gnn_reference_structure.xyz`: the reference structure as ExtXYZ.
- `torch_gnn_train_config.json`: the effective config for the run.
- `torch_gnn_training_metrics.json`: training history and output paths.

If `--reference-structure` is omitted, the training script defaults to the periodic diamond carbon structure used by the original workflow. For any non-carbon workflow, pass an explicit reference structure file.

## Inverse-design script

The inverse-design entrypoint consumes:

- an input structure file,
- a saved target fingerprint file,
- a saved model checkpoint,
- one set of inverse-design hyperparameters.

Example run:

```bash
python example/python_pkg/torch_gnn_inverse_design.py \
	--input-structure build/my_inputs/perturbed_diamond.xyz \
	--target-fingerprint build/torch_gnn_carbon_train/torch_gnn_target_fingerprint.npy \
	--model-checkpoint build/torch_gnn_carbon_train/torch_gnn_model_checkpoint.pt \
	--target-structure build/torch_gnn_carbon_train/torch_gnn_reference_structure.xyz \
	--inverse-steps 100 \
	--inverse-step-size 0.1 \
	--fingerprint-loss-weight 0.275 \
	--target-vertex-weight 0.55 \
	--save-optimisation-traj \
	--plot-2body-fingerprint-comparison \
	--fixed-leading-atoms 0 \
	--output-dir build/torch_gnn_carbon_inverse
```

The inverse-design script writes these outputs into `--output-dir`:

- `torch_gnn_inverse_design_final.xyz`: the optimised structure as ExtXYZ.
- `torch_gnn_inverse_design_metrics.json`: machine-readable success metrics and run settings.
- `torch_gnn_inverse_design_metrics.log`: a short human-readable metric summary.

Optional outputs controlled by flags:

- `--save-optimisation-traj`: writes `torch_gnn_inverse_design_path.traj`, containing the initial structure for each restart and the structure after every inverse-design optimisation step.
- `--plot-2body-fingerprint-comparison`: writes `torch_gnn_inverse_design_2body_fingerprint.png`, comparing the target 2-body fingerprint against the model-inferred 2-body fingerprint for the final optimised structure.

If `--target-structure` is provided, the script also records symmetry-aware RMSD and per-atom displacement metrics. Non-zero `--target-vertex-weight` or `--target-position-weight` requires `--target-structure`, because those losses depend on the target structure itself rather than only the fingerprint.

For periodic systems, prefer ExtXYZ-compatible files so the lattice and PBC are preserved across the training and inverse-design steps.

## Validated split workflow

The split workflow was validated end to end after `python -m pip install -e .` using:

- `torch_gnn_train_model.py` with a 1-epoch smoke test, 8 training structures, and 2 augmented structures.
- `torch_gnn_inverse_design.py` with a saved checkpoint, a saved target fingerprint, and a perturbed diamond input structure.

That smoke test reduced fingerprint MSE from `2.235e-3` to `1.282e-3` and wrote the expected checkpoint, fingerprint, config, structure, metrics JSON, and metrics log artifacts. The short smoke test was intended to validate the split workflow and file formats rather than to match the best structural RMSD.

## Best carbon result from the original workflow

The strongest result still comes from the original carbon example with the following settings:

```bash
python torch_gnn_carbon_workflow_example.py --carbon-count 0 --hidden-dim 128 --num-message-layers 4 --batch-size 8 --reference-layer-type 1 --component-weight-2body 3 --component-weight-3body 1 --component-weight-4body 0 --epochs 100 --epoch-values 10,100 --inverse-steps 100 --inverse-step-size 0.1 --inverse-step-values 0,25,100 --step-size-values 0.021,0.1 --fingerprint-loss-weight 0.275 --target-vertex-weight 0.55 --learning-rate 0.0005 --model-lr-decay-rate 0.005 --augmented-count 16 --seed 11 --fixed-leading-atoms 0
```

That run reaches a final RMSD of about `0.02 A` from an initial RMSD of about `0.5 A`.
