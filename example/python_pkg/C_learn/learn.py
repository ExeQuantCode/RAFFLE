from mace.calculators import mace_mp
from raffle.generator import raffle_generator
from ase import build, Atoms
from ase.optimize import BFGS, FIRE
from ase.io import write
import numpy as np
import os
from concurrent.futures import ProcessPoolExecutor, wait, as_completed
from multiprocessing import Process
from copy import deepcopy
from multiprocessing import Queue
from joblib import Parallel, delayed

import logging
logging.basicConfig(level=logging.DEBUG)

def runInParallel(*fns):
    proc = []
    results = []
    for fn in fns:
        p = Process(target=fn)
        p.start()
        proc.append(p)
    for p in proc:
        results.append(p.join())

    print("All processes finished")
    print(results)


def process_structure_with_queue(i, structure, num_old, calc_params, optimise_structure, queue):
    # Perform the computation
    result = process_structure(i, structure, num_old, calc_params, optimise_structure)
    queue.put(result)  # Report completion

def process_structure(i, atoms, num_structures_old, calc_params, optimise_structure):
    if i < num_structures_old:
        return
    
    inew = i - num_structures_old
    atoms.calc = calc
    # positions_initial = atoms.get_positions()

    # Calculate and save the initial energy per atom
    energy_unrlxd = atoms.get_potential_energy() / len(atoms)
    print(f"Initial energy per atom: {energy_unrlxd}")

    # Optimise the structure if requested
    if optimise_structure:
        optimizer = FIRE(atoms, trajectory = f"traje{inew}.traj", logfile=f"optimisation{inew}.log")
        try:
            optimizer.run(fmax=0.05, steps=500)
        except Exception as e:
            print(f"Optimisation failed: {e}")
            return None, None, None
    
    # Save the optimised structure and its energy per atom
    energy_rlxd = atoms.get_potential_energy() / len(atoms)

    # Get the distance matrix
    distances = atoms.get_all_distances(mic=True)

    # Set self-distances (diagonal entries) to infinity to ignore them
    np.fill_diagonal(distances, np.inf)

    # Check if the minimum non-self distance is below 1.5
    if distances.min() < 1.0:
        print(f"Distance too small: {atoms.get_all_distances(mic=True).min()}")
        return None, None, None
    
    if abs(energy_rlxd - energy_unrlxd) > 5.0:
        print(f"Energy difference too large: {energy_rlxd} vs {energy_unrlxd}")
        return None, None, None
    
    return atoms, energy_unrlxd, energy_rlxd

# crystal_structures = [
#     'sc', 'fcc', 'bcc', 'hcp',
#     'diamond', 'zincblende', 'rocksalt', 'cesiumchloride',
#     'fluorite', 'wurtzite', 'tetragonal', 'orthorhombic',
#     'bct', 'rhombohedral', 'mcl'
# ]



if __name__ == "__main__":

    calc_params = {
        "model": "medium",
        "dispersion": False,
        "default_dtype": "float32",
        "device": 'cpu'
    }
    calc = mace_mp(model="medium", dispersion=False, default_dtype="float32", device='cpu')

    crystal_structures = [
        'orthorhombic', 'diamond',
        # 'bct', 'sc',
        # 'fcc', 'bcc', 'hcp',
    ]

    hosts = []
    for crystal_structure in crystal_structures:
        print(f'Crystal structure: {crystal_structure}')
        for a in [5.0, 6.0, 7.0, 8.0, 9.0, 10.0]:
            b = a
            c = a
            atom = build.bulk(
                    name = 'C',
                    crystalstructure = crystal_structure,
                    a = a,
                    b = b,
                    c = c,
            )
            hosts.append(Atoms('C', positions=[(0, 0, 0)]))
            hosts[-1].set_cell(atom.get_cell())
            hosts[-1].calc = calc
            print(hosts[-1])


    generator = raffle_generator()
    generator.distributions.set_element_energies(
        {
            'C': 0.0
        }
    )
    # set energy scale
    generator.distributions.set_kBT(0.2)
    # set the distribution function widths (2-body, 3-body, 4-body)
    generator.distributions.set_width([0.025, np.pi/200.0, np.pi/200.0])

    initial_database = [Atoms('C', positions=[(0, 0, 0)], cell=[5, 5, 5])]
    initial_database[0].calc = calc

    # generator.distributions.create(initial_database)
    generator.distributions.read_gdfs("DATTEMPT4/gdfs.txt")

    if os.path.exists("energies_rlxd.txt"):
        with open("energies_rlxd.txt", "w") as energy_file:
            pass
    else:
        open("energies_rlxd.txt", "w").close()

    if os.path.exists("energies_unrlxd.txt"):
        with open("energies_unrlxd.txt", "w") as energy_file:
            pass
    else:
        open("energies_unrlxd.txt", "w").close()
    seed = -1
    num_structures_old = 0
    optimise_structure = True
    mass = 12.011
    for host in hosts:
        print("setting host")
        generator.set_host(host)
        volume = host.get_volume()
        for num_atoms in range(10, 80):

            density = ( num_atoms + 1 ) * mass / volume
            if density > 2.4: # aluminium density is 2.7 g/cm^3 ~= 1.61 u/A^3
                # print("Density too high:", density, "u/A^3")
                continue
            elif density < 2.0:
                # print("Density too low", density, "u/A^3")
                continue

            for i in range(2):
                seed += 1
                print(f"Seed: {seed}")
                generator.generate(
                    num_structures = 5,
                    stoichiometry = { 'C': num_atoms },
                    seed = seed,
                    method_probab = {"void": 0.5, "rand": 0.0, "walk": 0.5, "grow": 0.0, "min": 1.0},
                    verbose = 0,
                )

                # print the number of structures generated
                print("Total number of structures generated: ", generator.num_structures)
                generated_structures = generator.get_structures(calc)
                num_structures_new = len(generated_structures)

                # check if directory iteration[iter] exists, if not create it
                iterdir = f"iteration{seed}/"
                if not os.path.exists(iterdir):
                    os.makedirs(iterdir)
                generator.print_settings(iterdir+"generator_settings.txt")

                # set up list of energies
                energy_unrlxd = np.zeros(num_structures_new - num_structures_old)
                energy_rlxd = np.zeros(num_structures_new - num_structures_old)
                for i in range(num_structures_new - num_structures_old):
                    write(iterdir+f"POSCAR_unrlxd_{i}", generated_structures[num_structures_old + i])
                    print(f"Structure {i} energy per atom: {generated_structures[num_structures_old + i].get_potential_energy() / len(generated_structures[num_structures_old + i])}")

                # Start parallel execution
                print("Starting parallel execution")
                results = Parallel(n_jobs=5)(
                    delayed(process_structure)(i, deepcopy(generated_structures[i]), num_structures_old, calc_params, optimise_structure)
                    for i in range(num_structures_old, num_structures_new)
                )

                # Wait for all futures to complete
                for j, result in enumerate(results):
                    generated_structures[j+num_structures_old], energy_unrlxd[j], energy_rlxd[j] = result
                print("All futures completed")

                # Remove structures that failed the checks
                for j, atoms in reversed(list(enumerate(generated_structures))):
                    if atoms is None:
                        energy_unrlxd = np.delete(energy_unrlxd, j-num_structures_old)
                        energy_rlxd = np.delete(energy_rlxd, j-num_structures_old)
                        del generated_structures[j]
                        generator.remove_structure(j)
                num_structures_new = len(generated_structures) 

                # write the structures to files
                for i in range(num_structures_new - num_structures_old):
                    write(iterdir+f"POSCAR_{i}", generated_structures[num_structures_old + i])
                    print(f"Structure {i} energy per atom: {energy_rlxd[i]}")
                    # append energy per atom to the energies_unrlxd.txt file
                    with open("energies_unrlxd.txt", "a") as energy_file:
                        energy_file.write(f"{i+num_structures_old} {energy_unrlxd[i]}\n")
                    # append energy per atom to the energies_rlxd.txt file
                    with open("energies_rlxd.txt", "a") as energy_file:
                        energy_file.write(f"{i+num_structures_old} {energy_rlxd[i]}\n")

                # update the distribution functions
                print("Updating distributions")
                generator.distributions.update(generated_structures[num_structures_old:], from_host=False, deallocate_systems=False)

                # print the new distribution functions to a file
                print("Printing distributions")
                generator.distributions.write_dfs(iterdir+"distributions.txt")
                generator.distributions.write_2body(iterdir+"df2.txt")
                generator.distributions.write_3body(iterdir+"df3.txt")
                generator.distributions.write_4body(iterdir+"df4.txt")

                # update the number of structures generated
                num_structures_old = num_structures_new

    generator.distributions.write_gdfs("gdfs.txt")

    # Read energies from the file
    with open("energies_rlxd.txt", "r") as energy_file:
        energies = energy_file.readlines()

    # Parse and sort the energies
    energies = [line.strip().split() for line in energies]
    energies = sorted(energies, key=lambda x: float(x[1]))

    # Write the sorted energies back to the file
    with open("energies_ordered.txt", "w") as energy_file:
        for entry in energies:
            energy_file.write(f"{int(entry[0])} {float(entry[1])}\n")
