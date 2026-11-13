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
Fortran wrapper and demonstrate joint diamond inverse design.

The current comparison on conventional cubic diamond agrees with the Fortran
2- and 3-body fingerprints to below ``2e-6`` maximum absolute error. The
4-body fingerprint differs by about ``1.9e-3`` maximum absolute error; exact
4-body parity has not yet been established.
