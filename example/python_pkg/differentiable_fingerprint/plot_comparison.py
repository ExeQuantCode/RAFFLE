"""Plot Fortran and PyTorch fingerprints for a structure."""

import argparse
from itertools import combinations_with_replacement

import matplotlib.pyplot as plt
import numpy as np
from ase.build import bulk
from ase.io import read

from raffle.fingerprint.differentiable import DifferentiableFingerprint
from raffle.raffle import Raffle__Distribs_Container


def plot_fingerprints(reference, actual, grids, elements, reference_label, actual_label,
                      output, show=False, body_orders=(2, 3, 4)):
    series = (
        ("2-body", "two_body", reference[0], actual[0],
         ["-".join(pair) for pair in combinations_with_replacement(elements, 2)], r"Distance ($\AA$)"),
        ("3-body", "three_body", reference[1], actual[1], elements, "Angle (rad)"),
        ("4-body", "four_body", reference[2], actual[2], elements, "Angle (rad)"),
    )
    series = tuple(item for order, item in zip((2, 3, 4), series) if order in body_orders)
    nrows = sum(len(labels) for _, _, _, _, labels, _ in series)
    figure, axes = plt.subplots(nrows, 2, figsize=(12, 2.8 * nrows), squeeze=False)
    row = 0
    for body, key, expected, values, labels, xlabel in series:
        grid = grids[key]
        if expected.shape != values.shape or values.shape[1] != len(labels):
            raise ValueError(f"{body} fingerprint channel shapes do not match")
        for channel, label in enumerate(labels):
            difference = values[:, channel] - expected[:, channel]
            axes[row, 0].plot(grid, expected[:, channel], label=reference_label)
            axes[row, 0].plot(grid, values[:, channel], "--", label=actual_label)
            axes[row, 0].set_title(f"{body}: {label}")
            axes[row, 0].set_ylabel("Normalized fingerprint")
            axes[row, 0].legend()
            axes[row, 1].plot(grid, difference)
            axes[row, 1].axhline(0, color="black", linewidth=0.6)
            axes[row, 1].set_title(
                f"{actual_label} - {reference_label} (max |difference| = {np.max(np.abs(difference)):.3g})"
            )
            axes[row, 1].set_ylabel("Difference")
            for axis in axes[row]:
                axis.set_xlabel(xlabel)
            row += 1
    figure.tight_layout()
    figure.savefig(output, dpi=150)
    print(f"Saved fingerprint comparison to {output}")
    if show:
        plt.show()
    plt.close(figure)


def plot_fingerprint_results(target, actual, output):
    if target.elements != actual.elements or target.body_orders != actual.body_orders:
        raise ValueError("fingerprints must have matching elements and body orders")
    for name in target.grids:
        if not np.array_equal(
            target.grids[name].detach().cpu().numpy(),
            actual.grids[name].detach().cpu().numpy(),
        ):
            raise ValueError(f"{name} fingerprint grids do not match")
    plot_fingerprints(
        [value.detach().cpu().numpy() for value in (target.two_body, target.three_body, target.four_body)],
        [value.detach().cpu().numpy() for value in (actual.two_body, actual.three_body, actual.four_body)],
        {name: grid.detach().cpu().numpy() for name, grid in target.grids.items()},
        target.elements,
        "Target",
        "Final",
        output,
        body_orders=target.body_orders,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "structure", nargs="?", help="ASE-readable structure (default: cubic diamond)"
    )
    parser.add_argument(
        "--output", default="fingerprint_comparison.png", help="Output image path"
    )
    parser.add_argument("--show", action="store_true", help="Also display the figure")
    args = parser.parse_args()

    atoms = read(args.structure) if args.structure else bulk("C", "diamond", a=3.567, cubic=True)
    elements = sorted(set(atoms.get_chemical_symbols()))
    parameters = Raffle__Distribs_Container.distribs_container_type()
    descriptor = DifferentiableFingerprint(elements, parameters=parameters)
    fortran = parameters.generate_fingerprint(atoms)
    python = descriptor.generate(atoms)
    plot_fingerprints(
        fortran,
        [value.detach().cpu().numpy() for value in (python.two_body, python.three_body, python.four_body)],
        {name: grid.detach().cpu().numpy() for name, grid in python.grids.items()},
        elements,
        "Fortran",
        "PyTorch",
        args.output,
        args.show,
    )


if __name__ == "__main__":
    main()
