Differentiable fingerprints and inverse design
==============================================

The optional PyTorch implementation can calculate RAFFLE fingerprints directly
from positions and continuous element probabilities. Install the extra with:

.. code-block:: bash

   pip install "raffle[differentiable]"

The differentiable modules are exposed from
``raffle.fingerprint.differentiable``. The descriptor uses the standard RAFFLE
distribution parameters by default, or accepts an existing
``Raffle__Distribs_Container`` or a mapping of descriptor parameters.

Generating fingerprints
-----------------------

.. code-block:: python

   from raffle.fingerprint.differentiable import DifferentiableFingerprint

   descriptor = DifferentiableFingerprint(elements=["C", "Si"])
   target = descriptor.generate(target_atoms)
   fingerprint = descriptor(
       positions=positions,
       cell=target_atoms.cell.array,
       species_probabilities=probabilities,
       pbc=target_atoms.pbc,
   )

Element channels are ordered alphabetically, matching RAFFLE's Fortran
structure conversion. Each body order is returned as a bin-by-channel tensor.
``FingerprintResult.unnormalized`` contains the pre-channel-normalisation
contributions; ``grids`` contains each bin grid. A target generated from ASE
also stores its element counts for composition-constrained relaxations.

ASE position-only calculator
----------------------------

``FingerprintCalculator`` supplies the fingerprint mismatch as an ASE energy
and its negative position gradient as forces. These values are descriptor
objective units, **not** physical eV or eV/Angstrom. Consequently, ASE's
``fmax`` threshold is a gradient tolerance in descriptor units, not a physical
force criterion. No stress is provided.

.. code-block:: python

   from ase.optimize import FIRE
   from raffle.fingerprint.differentiable import FingerprintCalculator

   atoms.calc = FingerprintCalculator(descriptor, target)
   FIRE(atoms).run(fmax=1e-4)

The calculator keeps chemical identities fixed unless soft probabilities are
explicitly supplied with ``set_species_probabilities``. Calling that method
invalidates ASE's cached results.

Joint position and species relaxation
--------------------------------------

ASE's normal optimizers only change positions. ``FingerprintRelaxation`` uses
an Adam loop for positions and species logits, and performs a global
minimum-cost species assignment before optional final positional refinement.
``fixed_position_mask`` and ``fixed_species_mask`` are independent Boolean
arrays of length ``n_atoms``.

The initial trajectory frame retains the input species. Subsequent joint frames
and callbacks show a discrete assignment projected from the probabilities used
to evaluate that frame; exact mode enforces the requested composition in that
projection. The corresponding continuous probabilities are available in
``species_probabilities_history``. For example, an exact all-carbon target
forces every unfixed atom to carbon from the first evaluation, so its species
gradient is zero even if the input structure contained silicon.

Set ``verbose=True`` in ``run`` to print the joint and final-refinement steps,
including the objective, maximum position-gradient norm (``Fmax``), maximum
species-logit gradient norm, and maximum atomic step in Angstrom. Gradients
are in objective units, not physical forces. The loss is measured before the
displayed step. For example:

.. code-block:: python

   result = relaxation.run(initial_atoms, max_steps=500, verbose=True)

The joint loop stops after five consecutive relative loss changes at or below
``tolerance`` (default ``1e-8``); set ``tolerance=0`` to run all ``max_steps``.
This avoids premature stopping when the absolute fingerprint loss is small.
Position and species step sizes can be adjusted independently using
``position_learning_rate`` and ``species_learning_rate``. The final
``final_refine_steps`` position-only steps use Adam after species assignment,
and return the lowest-loss positions encountered. ``converged`` and ``n_steps``
refer to the joint loop.

The fingerprint objective is non-convex: small gradients and a low loss do not
guarantee the target structure or fingerprint has been recovered. For larger
position perturbations, ``radial_smoothing`` optionally smooths the 2-body
fingerprint difference with a Gaussian of the given standard deviation in
Angstrom. Its default is zero, leaving the original objective unchanged. A
coarse radial-only run can use this smoother objective, followed by a run with
``radial_smoothing=0`` against the original fingerprint. Because smoothing
changes the loss scale, coarse and fine losses are not directly comparable.
This can help escape narrow radial peaks but does not guarantee a global
minimum. ``inverse_design.py`` demonstrates this two-stage approach.

``composition_constraint`` has three modes:

* ``0.0``: unrestricted composition, with no composition penalty.
* A finite positive value: a penalty strength applied to the squared
  difference between predicted and requested element fractions.
* ``float("inf")`` or ``"exact"``: Sinkhorn-constrained probabilities and an
  exact final integer assignment.

Specify ``target_composition`` as element counts when it is unavailable in the
target fingerprint metadata. Fixed-species counts are subtracted before the
exact constraint is applied, and infeasible requests raise ``ValueError``.

Runnable examples are in
``example/python_pkg/differentiable_fingerprint/``. They compare against the
Fortran wrapper and demonstrate joint diamond inverse design. The
``inverse_design.py`` example writes ``diamond_fingerprint_comparison.png``,
which overlays the final and target fingerprints and plots their differences
for every element channel. The plotting examples require Matplotlib.

To plot the Fortran and PyTorch fingerprints and their difference, install
Matplotlib and run:

.. code-block:: bash

   python example/python_pkg/differentiable_fingerprint/plot_comparison.py
   python example/python_pkg/differentiable_fingerprint/plot_comparison.py structure.cif --output comparison.png --show

Without a structure argument the script uses cubic diamond. It writes
``fingerprint_comparison.png`` by default; ``--show`` also opens the plot.
The structure can be any ASE-readable file.

The current comparison on conventional cubic diamond agrees with the Fortran
2- and 3-body fingerprints to below ``2e-6`` maximum absolute error. The
4-body fingerprint differs by about ``1.9e-3`` maximum absolute error; exact
4-body parity has not yet been established.
