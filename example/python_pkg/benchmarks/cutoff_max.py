import os
import pytest
import time
import numpy as np
from ase import Atoms
from ase.io import read
from chgnet.model import CHGNetCalculator
from raffle.generator import raffle_generator


@pytest.mark.parametrize("kBT", [0.4])
@pytest.mark.parametrize("width", [
    [0.025, np.pi/200.0, np.pi/200.0]
])
@pytest.mark.parametrize("grid_spacing", [0.1])
@pytest.mark.parametrize("method_ratio", [
  {'void': 0.0, 'rand': 0.0, 'walk': 0.0, 'grow': 0.0, 'min': 1.0}
])
@pytest.mark.parametrize("cutoff_max", [4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0])
@pytest.mark.benchmark(
    group="cutoff_in_evaluator",
    min_rounds=5,
    timer=time.time,
    disable_gc=True,
    warmup=False
)
def test_cutoff_in_evaluator(benchmark, grid_spacing, kBT, width, method_ratio, cutoff_max):
    
    calc = CHGNetCalculator()
    current_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(current_dir, '..', '..', '..', 'example', 'data')
    filename = os.path.join(data_dir, 'diamond.xyz')
    atoms = read(filename, index=":")


    # Initialise generator
    generator = raffle_generator()

    # set energy scale
    generator.distributions.set_kBT(kBT)

    # set the distribution function widths (2-body, 3-body, 4-body)
    generator.distributions.set_width(width)

    # get list of elements in list of atoms
    elements = []
    for atom in atoms:
        elements += atom.get_chemical_symbols()
    elements = list(set(elements))
    energies = {element: 0.0 for element in elements}
    generator.distributions.set_element_energies(energies)

    # set the cutoff_max for the distributions
    generator.distributions.set_cutoff_max([cutoff_max, np.pi, np.pi])

    generator.set_grid(grid_spacing=grid_spacing)

    generator.distributions.create(atoms)
    
    benchmark.extra_info = {
        "kBT": kBT,
        "width": width,
        "grid_spacing": grid_spacing,
        "method_ratio": method_ratio,
        "cutoff_max": cutoff_max
    }

    # Measure time for distributions.evaluate()
    result = benchmark(generator.evaluate, atoms[0])

    return result


@pytest.mark.parametrize("kBT", [0.4])
@pytest.mark.parametrize("width", [
    [0.025, np.pi/200.0, np.pi/200.0]
])
@pytest.mark.parametrize("grid_spacing", [0.1])
@pytest.mark.parametrize("method_ratio", [
  {'void': 0.0, 'rand': 0.0, 'walk': 0.0, 'grow': 0.0, 'min': 1.0}
])
@pytest.mark.parametrize("cutoff_max", [4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0])
@pytest.mark.benchmark(
    group="cutoff_in_database",
    min_rounds=5,
    timer=time.time,
    disable_gc=True,
    warmup=False
)
def test_cutoff_in_database(benchmark, grid_spacing, kBT, width, method_ratio, cutoff_max):
    
    calc = CHGNetCalculator()
    current_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(current_dir, '..', '..', '..', 'example', 'data')
    filename = os.path.join(data_dir, 'diamond.xyz')
    atoms = read(filename, index=":")


    # Initialise generator
    generator = raffle_generator()

    # set energy scale
    generator.distributions.set_kBT(kBT)

    # set the distribution function widths (2-body, 3-body, 4-body)
    generator.distributions.set_width(width)

    # get list of elements in list of atoms
    elements = []
    for atom in atoms:
        elements += atom.get_chemical_symbols()
    elements = list(set(elements))
    energies = {element: 0.0 for element in elements}
    generator.distributions.set_element_energies(energies)

    # set the cutoff_max for the distributions
    generator.distributions.set_cutoff_max([cutoff_max, np.pi, np.pi])

    generator.set_grid(grid_spacing=grid_spacing)

    benchmark.extra_info = {
        "kBT": kBT,
        "width": width,
        "grid_spacing": grid_spacing,
        "method_ratio": method_ratio,
        "cutoff_max": cutoff_max
    }

    # Measure time for distributions.create()
    result = benchmark(generator.distributions.create, atoms)

    return result

if __name__ == '__main__':
    cutoff_max = [6.0]
    results = []

    for cutoff in cutoff_max:
        for _ in range(5):  # Run the benchmark 5 times for each cutoff
            result = test_cutoff_in_database(lambda x: x, 0.1, 0.4, [0.025, np.pi/200.0, np.pi/200.0], {'void': 0.0, 'rand': 0.0, 'walk': 0.0, 'grow': 0.0, 'min': 1.0}, cutoff)
    
    for cutoff in cutoff_max:
        for _ in range(5):  # Run the benchmark 5 times for each cutoff
            result = test_cutoff_in_evaluator(lambda x: x, 0.1, 0.4, [0.025, np.pi/200.0, np.pi/200.0], {'void': 0.0, 'rand': 0.0, 'walk': 0.0, 'grow': 0.0, 'min': 1.0}, cutoff)

    # Print results
    print("\nBenchmark Results:")
    print("-" * 40)
    print(f"{'Grid Size':<10} {'Time (s)':<10}")
    print("-" * 40)
    for result in results:
        print(f"{result['grid_spacing']:<10} {result['time']:<10.4f}")
