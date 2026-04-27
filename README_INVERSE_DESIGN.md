
There is an inverse design example in `example/python_pkg/` that demonstrates a neural network can be trained to learn RAFFLE descriptors and then used to perform inverse design as a form of surrogate atomic relaxation to find structurally similar structures to the target RAFFLE fingerprint.

The current approach uses PyTorch to implement and train these neural networks, but the code is modular and can be adapted to other frameworks such as ATHENA (with the framework already implemented, just need to find a network architecture that works for this problem).

The example script is `example/python_pkg/torch_gnn_carbon_workflow_example.py`.
It takes a dataset of carbon structures, generates their RAFFLE descriptors and trains a graph neural network as a surrogate to RAFFLE (learning the mapping from structure to RAFFLE descriptor).
Then, the trained model is used to perform inverse design by taking a perturbed diamond structure, calculating its RAFFLE descriptor, then using the trained model to optimise the atomic positions to best match the perfect diamond structure (i.e. relax a perturbed diamond structure to a perfect diamond structure).

The best achieving model is currently run using:

```python
python torch_gnn_carbon_workflow_example.py --carbon-count 0 --hidden-dim 128 --num-message-layers 4 --batch-size 8 --reference-layer-type 1 --component-weight-2body 3 --component-weight-3body 1 --component-weight-4body 0 --epochs 100 --epoch-values 10,100 --inverse-steps 100 --inverse-step-size 0.1 --inverse-step-values 0,25,100 --step-size-values 0.021,0.1 --fingerprint-loss-weight 0.275 --target-vertex-weight 0.55 --learning-rate 0.0005 --model-lr-decay-rate 0.005 --augmented-count 16 --seed 11 --fixed-leading-atoms 0
```

This returns a final RMSD for the atomic positions of 0.02 Angstrom, where the initial RMSD is 0.5 Angstrom.
