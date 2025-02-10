import os
import glob
from raffle.generator import raffle_generator
import pytest
from ase.io import read
import time
import numpy as np


def get_xyz_files():
    """Get all .xyz files from the ../../data directory"""
    current_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(current_dir, '..', '..', '..', 'example', 'data')
    print("Data directory:", data_dir)
    return glob.glob(os.path.join(data_dir, '*.xyz'))


@pytest.mark.parametrize("xyz_file", get_xyz_files())
@pytest.mark.parametrize("kBT", [0.2, 0.4, 0.6])
@pytest.mark.parametrize("width", [
    [0.025, np.pi/200.0, np.pi/200.0],
    [0.05, np.pi/100.0, np.pi/100.0],
    [0.1, np.pi/50.0, np.pi/50.0]
])
@pytest.mark.benchmark(
    group="distributions",
    min_rounds=5,
    timer=time.time,
    disable_gc=True,
    warmup=False
)
def test_distributions_create(benchmark, xyz_file, kBT, width):
    atoms = read(xyz_file, index=":")
    
    # Initialize generator
    generator = raffle_generator()

    # get list of elements in list of atoms
    elements = []
    for atom in atoms:
        elements += atom.get_chemical_symbols()
    elements = list(set(elements))
    energies = {element: 0.0 for element in elements}
    generator.distributions.set_element_energies(energies)

    avg_num_atoms = 0
    for atom in atoms:
        avg_num_atoms += len(atom)
    avg_num_atoms = round(avg_num_atoms / len(atoms))

    benchmark.extra_info = {
        'file': os.path.basename(xyz_file),
        'num_structures': len(atoms),
        'num_species': len(elements),
        "avg_num_atoms": avg_num_atoms,
        "kBT": kBT,
        "width": width
    }

    # set energy scale
    generator.distributions.set_kBT(kBT)

    # set the distribution function widths (2-body, 3-body, 4-body)
    generator.distributions.set_width(width)
    
    # Measure time for distributions.create()
    result = benchmark(generator.distributions.create, atoms)

    
    return {
        'file': os.path.basename(xyz_file),
        'time': result,
        'num_structures': len(atoms),
        'num_species': len(elements),
        "avg_num_atoms": avg_num_atoms
    }

@pytest.fixture
def make_customer_record():
    def _make_customer_record(name):
        return {"name": name, "orders": []}

    return _make_customer_record

if __name__ == '__main__':
    xyz_files = get_xyz_files()
    results = []
    
    result = test_distributions_create(lambda x: x)
    # for xyz_file in xyz_files:
    #     # for _ in range(5):  # Run the benchmark 5 times for each file
    #     result = test_distributions_create(lambda x: x, xyz_file)
    #     results.append(result)
    
    # Print results
    print("\nBenchmark Results:")
    print("-" * 60)
    print(f"{'File':<20} {'Time (s)':<10} {'Num Structures':<10} {'Num Species':<10} {'Avg Num Atoms':<10}")
    print("-" * 60)
    for result in results:
        print(f"{result['file']:<20} {result['time']:<10.4f} {result['num_structures']:<10} {result['num_species']:<10} {result['avg_num_atoms']:<10}")