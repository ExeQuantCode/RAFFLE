import pytest
import time
import numpy as np
from ase import Atoms
from chgnet.model import CHGNetCalculator
from raffle.generator import raffle_generator


@pytest.mark.parametrize("kBT", [0.4])
@pytest.mark.parametrize("width", [
    [0.025, np.pi/200.0, np.pi/200.0]
])
@pytest.mark.parametrize("grid_spacing", [0.5, 0.25, 0.1, 0.075, 0.05])
@pytest.mark.parametrize("method_ratio", [
  {'void': 0.0, 'rand': 0.0, 'walk': 0.0, 'grow': 0.0, 'min': 1.0},
  {'void': 0.0, 'rand': 0.0, 'walk': 0.0, 'grow': 1.0, 'min': 0.0},
  {'void': 0.0, 'rand': 0.0, 'walk': 1.0, 'grow': 0.0, 'min': 0.0},
  {'void': 0.0, 'rand': 1.0, 'walk': 0.0, 'grow': 0.0, 'min': 0.0},
  {'void': 1.0, 'rand': 0.0, 'walk': 0.0, 'grow': 0.0, 'min': 0.0}
])
@pytest.mark.benchmark(
    group="placement_methods",
    min_rounds=5,
    timer=time.time,
    disable_gc=True,
    warmup=False
)
def test_placement_methods(benchmark, grid_spacing, kBT, width, method_ratio):
    
    calc = CHGNetCalculator()
    host = Atoms('C', positions=[(0, 0, 0)], cell=[3.567, 3.567, 3.567], pbc=True, calculator=calc)

    # Initialize generator
    generator = raffle_generator()

    # set energy scale
    generator.distributions.set_kBT(kBT)

    # set the distribution function widths (2-body, 3-body, 4-body)
    generator.distributions.set_width(width)

    energies = {element: 0.0 for element in ['C']}
    generator.distributions.set_element_energies(energies)

    generator.set_grid(grid_spacing=grid_spacing)
    generator.set_host(host)

    initial_database = [Atoms('C', positions=[(0, 0, 0)], cell=[8, 8, 8], pbc=True, calculator=calc)]

    generator.distributions.create(initial_database)
    
    benchmark.extra_info = {
        "kBT": kBT,
        "width": width,
        "grid_spacing": grid_spacing,
        "method_ratio": method_ratio
    }

    # Measure time for distributions.create()
    result = benchmark(generator.generate, 
        num_structures = 1,
        stoichiometry = {'C': 1},
        seed = 0,
        method_ratio = method_ratio,
        verbose = 0
    )

    return result

if __name__ == '__main__':
    grid_spacings = [0.5, 0.25, 0.1, 0.075, 0.05]
    results = []

    for grid_spacing in grid_spacings:
        for _ in range(5):  # Run the benchmark 5 times for each grid size
            result = test_placement_methods(lambda x: x, grid_spacing, 0.4, [0.025, np.pi/200.0, np.pi/200.0])
            results.append({
                'grid_spacing': grid_spacing,
                'time': result
            })
    
    # Print results
    print("\nBenchmark Results:")
    print("-" * 40)
    print(f"{'Grid Size':<10} {'Time (s)':<10}")
    print("-" * 40)
    for result in results:
        print(f"{result['grid_spacing']:<10} {result['time']:<10.4f}")
