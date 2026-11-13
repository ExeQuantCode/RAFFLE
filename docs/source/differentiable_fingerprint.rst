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
For expensive systems, pass ``body_orders=(2, 3)`` to the descriptor to skip
4-body computation. Disabled arrays contain zeros and are marked in
``FingerprintResult.body_orders``; requesting a nonzero loss weight for a
disabled order raises an error.

At fixed element and bin counts, the 2-, 3- and 4-body histogram calculations
scale approximately as ``N*k``, ``N*k**2`` and ``N*k**3``, respectively, for
``N`` atoms and roughly ``k`` neighbors per atom. Sorted neighbor groups are
cached until the topology is rebuilt. The worst case can still be fourth-order
in ``N`` when all atoms are neighbors. Autograd retains intermediate
calculations for the backward pass, and storing trajectory and probability
history consumes memory proportional to the number of saved frames. Exact
species assignment can also be costly for large systems.

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
Generate target and candidate fingerprints with matching elements and grids;
the loss does not verify their element identities or grid values.

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

The optional ``repulsion_weight`` defaults to zero. Its current joint-loop
penalty uses 0.7 times the descriptor's pair radius as its onset. It excludes
zero-shift self-pairs but includes periodic images of the same atom. For default
carbon radii, its onset is about 0.532 Angstrom.
The final position-only refinement does not include this penalty. It is not
a physical potential or a reliable hard minimum-distance constraint.
``result.loss`` is the final discrete fingerprint loss and excludes the
joint-loop penalties; compare it with joint ``loss_history`` accordingly.

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

``angular_smoothing`` applies the same kind of Gaussian smoothing to the
3- and 4-body differences, with its standard deviation in radians. The
convolution uses zero padding at the grid boundaries. Angular terms provide
coarse coordination information, but their gradients can compete with the
radial term and their contribution should be checked for each system.

The default ``loss_reduction="mean"`` averages squared bin differences.
``loss_reduction="density"`` instead divides the sum of squared probability
mass differences by the bin width and averages over channels, corresponding
to the integral of a squared density difference. This compensates for the
small per-bin values without a large arbitrary multiplier; the radial and
angular terms have different units, so weights are still meaningful.

Fine relaxation optimizes the original, unsmoothed fingerprint. It is useful
when the unsmoothed per-body mismatch remains significant after smoothing;
there is no universal distance threshold or guarantee that fine relaxation
will escape a local minimum. Compare the unsmoothed fingerprints, not coarse
and fine loss values, which use different objectives.

``composition_constraint`` has three modes:

* ``0.0``: unrestricted composition, with no composition penalty.
* A finite positive value: a penalty strength applied to the squared
  difference between predicted and requested element fractions.
* ``float("inf")`` or ``"exact"``: Sinkhorn-constrained probabilities and an
  exact final integer assignment.

Specify ``target_composition`` as element counts when it is unavailable in the
target fingerprint metadata. Fixed-species counts are subtracted before the
exact constraint is applied, and infeasible requests raise ``ValueError``.
The soft probabilities satisfy the requested counts to an absolute tolerance
of ``1e-5`` per element; the final discrete assignment uses exact counts.

Moving one atom
---------------

Use ``fixed_position_mask`` to keep every other atom fixed, and set
``optimise_species=False`` to preserve species:

.. code-block:: python

   fixed = [i != mobile_index for i in range(len(initial_atoms))]
   relaxation = FingerprintRelaxation(descriptor, target, fixed_position_mask=fixed)
   result = relaxation.run(initial_atoms, optimise_species=False)

This minimizes the **global** fingerprint mismatch with respect to one atom;
it is not a fingerprint restricted to that atom's local neighborhood. The
cell is fixed. ``local_atom.py`` gives a runnable example.

Runnable examples are in
``example/python_pkg/differentiable_fingerprint/``. They compare against the
Fortran wrapper and demonstrate joint diamond inverse design. The
``inverse_design.py`` example writes ``diamond_fingerprint_comparison.png``,
which overlays the final and target fingerprints and plots their differences
for every element channel. The plotting examples require Matplotlib.

``diamond_to_carbon_phases.py`` starts each run from the same perturbed
diamond structure and attempts to match separate graphite, FCC and BCC
fingerprints, using all three smoothed body orders. The FCC and BCC targets
are constructed at diamond's volume per atom; the graphite target uses
``a = 2.46`` and ``c = 6.70`` Angstrom. The 3-body upper distance tolerance is
expanded from 2.5 to 2.8 times the carbon covalent radius so the more distant
FCC and BCC neighbors contribute to the angular fingerprints. The diamond
cell stays fixed, so exact phase recovery is not guaranteed. The script prints unsmoothed errors and
writes structures, trajectories and comparison plots. ``--steps`` and
``--fine-steps`` control the two stages.

``complex_species.py`` reads the supplied CaCu3Ti4O12 POSCAR, replaces all
species with Ca, and optimizes species using 2- and 3-body fingerprints and
the known exact composition. By default it keeps the POSCAR's positions and
cell fixed; ``--relax-positions`` also moves atoms. Thus it is a site-label
reconstruction with target geometry already provided, not an independent
structure prediction. It skips 4-body for cost and reports both the soft loss
and the final discrete assignment. Very sharp species probabilities can
make the Sinkhorn constraint difficult to converge; increasing the temperature
or lowering the species learning rate can help.

.. code-block:: bash

   python example/python_pkg/differentiable_fingerprint/diamond_to_carbon_phases.py
   python example/python_pkg/differentiable_fingerprint/complex_species.py
   python example/python_pkg/differentiable_fingerprint/local_atom.py --atom 0

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
4-body fingerprint differs by about ``4.5e-5`` maximum absolute error; exact
4-body parity has not yet been established. A non-periodic four-carbon test
also showed a 2-body maximum difference of about ``3.5e-3``. These comparisons
are not a general parity guarantee.
