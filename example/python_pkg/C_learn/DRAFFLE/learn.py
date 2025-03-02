from chgnet.model import CHGNetCalculator
from ase.calculators.singlepoint import SinglePointCalculator
from raffle.generator import raffle_generator
from ase import build, Atoms
from ase.optimize import FIRE
from ase.io import write
import numpy as np
import os
from joblib import Parallel, delayed

import logging
logging.basicConfig(level=logging.DEBUG)

# Function to relax a structure
def process_structure(i, atoms, num_structures_old, calc_params, optimise_structure, iteration, calc):
    # Check if the structure has already been processed
    if i < num_structures_old:
        return None, None, None

    # calc = Vasp(**calc_params, label=f"struct{i}", directory=f"iteration{iteration}/struct{i}/", txt=f"stdout{i}.o")
    inew = i - num_structures_old
    atoms.calc = calc

    # Calculate and save the initial energy per atom
    energy_unrlxd = atoms.get_potential_energy() / len(atoms)
    print(f"Initial energy per atom: {energy_unrlxd}")

    # Optimise the structure if requested
    if optimise_structure:
        optimizer = FIRE(atoms, trajectory = f"traje{inew}.traj", logfile=f"optimisation{inew}.log")
        try:
            optimizer.run(fmax=0.05, steps=100)
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

    if abs(energy_rlxd - energy_unrlxd) > 10.0:
        print(f"Energy difference too large: {energy_rlxd} vs {energy_unrlxd}")
        return None, None, None

    return atoms, energy_unrlxd, energy_rlxd


if __name__ == "__main__":

    # set up the calculator
    calc_params = {}
    calc = CHGNetCalculator()

    # set up the hosts
    crystal_structures = [
        'orthorhombic', # 'hcp',
    ]

    hosts = []
    for crystal_structure in crystal_structures:
        print(f'Crystal structure: {crystal_structure}')
        for a in np.linspace(3.567, 3.567, num=1):
            b = a
            c = a
            atom = build.bulk(
                    name = 'C',
                    crystalstructure = crystal_structure,
                    a = a,
                    b = b,
                    c = c,
            )
            hosts.append(Atoms('C', positions=[(0, 0, 0)], cell=atom.get_cell(), pbc=True, calculator=calc))
            print(hosts[-1])

    # set the parameters for the generator
    optimise_structure = True
    mass = 12.011
    num_atoms = 7

    # Loop over random seeds
    for seed in range(20):
        print(f"Seed: {seed}")

        # set up file names
        energies_rlxd_filename = f"energies_rlxd_seed{seed}.txt"
        energies_unrlxd_filename = f"energies_unrlxd_seed{seed}.txt"
        generator = raffle_generator()
        generator.distributions.set_element_energies(
            {
                'C': 0.0
            }
        )
        # set energy scale
        generator.distributions.set_kBT(0.4)

        # set the distribution function widths (2-body, 3-body, 4-body)
        generator.distributions.set_width([0.025, np.pi/200.0, np.pi/200.0])

        # set the initial database
        initial_database = [Atoms('C', positions=[(0, 0, 0)], cell=[8, 8, 8], pbc=True)]
        initial_database[0].calc = calc
        generator.distributions.create(initial_database)

        # set the host
        generator.set_host(hosts[0])

        # check if the energies file exists, if not create it
        if os.path.exists(energies_rlxd_filename):
            with open(energies_rlxd_filename, "w") as energy_file:
                pass
        else:
            open(energies_rlxd_filename, "w").close()

        if os.path.exists(energies_unrlxd_filename):
            with open(energies_unrlxd_filename, "w") as energy_file:
                pass
        else:
            open(energies_unrlxd_filename, "w").close()

        # initialise the number of structures generated
        num_structures_old = 0
        unrlxd_structures = []
        rlxd_structures = []
        # start the iterations
        for iter in range(200):

            # generate the structures
            generator.generate(
                num_structures = 5,
                stoichiometry = { 'C': num_atoms },
                seed = seed*1000+iter,
                method_ratio = {"void": 0.5, "rand": 0.001, "walk": 0.5, "grow": 0.0, "min": 1.0},
                verbose = 0,
            )

            # print the number of structures generated
            print("Total number of structures generated: ", generator.num_structures)
            generated_structures = generator.get_structures(calc)
            num_structures_new = len(generated_structures)

            # check if directory iteration[iter] exists, if not create it
            iterdir = f"iteration{iter}/"
            if not os.path.exists(iterdir):
                os.makedirs(iterdir)
            generator.print_settings(iterdir+"generator_settings.txt")

            # set up list of energies
            energy_unrlxd = np.zeros(num_structures_new - num_structures_old)
            energy_rlxd = np.zeros(num_structures_new - num_structures_old)
            for i in range(num_structures_new - num_structures_old):
                write(iterdir+f"POSCAR_unrlxd_{i}", generated_structures[num_structures_old + i])
                print(f"Structure {i} energy per atom: {generated_structures[num_structures_old + i].get_potential_energy() / len(generated_structures[num_structures_old + i])}")
                unrlxd_structures.append(generated_structures[num_structures_old + i].copy())
                unrlxd_structures[-1].calc = SinglePointCalculator(
                    generated_structures[num_structures_old + i],
                    energy=generated_structures[num_structures_old + i].get_potential_energy(),
                    forces=generated_structures[num_structures_old + i].get_forces()
                )

            # Start parallel execution
            print("Starting parallel execution")
            results = Parallel(n_jobs=5)(
                delayed(process_structure)(i, generated_structures[i].copy(), num_structures_old, calc_params, optimise_structure, iteration=seed, calc=calc)
                for i in range(num_structures_old, num_structures_new)
            )

            # Wait for all futures to complete
            for j, result in enumerate(results):
                generated_structures[num_structures_old + j], energy_unrlxd[j], energy_rlxd[j] = result
                if generated_structures[num_structures_old + j] is None:
                    print("Structure failed the checks")
                    continue
                rlxd_structures.append(generated_structures[num_structures_old + j].copy())
                rlxd_structures[-1].calc = SinglePointCalculator(
                    generated_structures[j+num_structures_old],
                    energy=generated_structures[num_structures_old + j].get_potential_energy(),
                    forces=generated_structures[num_structures_old + j].get_forces()
                )
            print("All futures completed")

            # Remove structures that failed the checks
            for j, atoms in reversed(list(enumerate(generated_structures))):
                if atoms is None:
                    energy_unrlxd = np.delete(energy_unrlxd, j-num_structures_old)
                    energy_rlxd = np.delete(energy_rlxd, j-num_structures_old)
                    del generated_structures[j]
                    del unrlxd_structures[j]
                    del rlxd_structures[j]
                    generator.remove_structure(j)
            num_structures_new = len(generated_structures)

            # write the structures to files
            for i in range(num_structures_new - num_structures_old):
                write(iterdir+f"POSCAR_{i}", generated_structures[num_structures_old + i])
                print(f"Structure {i} energy per atom: {energy_rlxd[i]}")
                # append energy per atom to the 'energies_unrlxd_filename' file
                with open(energies_unrlxd_filename, "a") as energy_file:
                    energy_file.write(f"{i+num_structures_old} {energy_unrlxd[i]}\n")
                # append energy per atom to the 'energies_rlxd_filename' file
                with open(energies_rlxd_filename, "a") as energy_file:
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
            generator.distributions.deallocate_systems()

            # update the number of structures generated
            num_structures_old = num_structures_new

        # Write the final distribution functions to a file
        generator.distributions.write_gdfs(f"gdfs_seed{seed}.txt")

        # Read energies from the file
        with open(energies_rlxd_filename, "r") as energy_file:
            energies = energy_file.readlines()

        # Parse and sort the energies
        energies = [line.strip().split() for line in energies]
        energies = sorted(energies, key=lambda x: float(x[1]))

        # Write the sorted energies back to the file
        with open(f"sorted_{energies_rlxd_filename}", "w") as energy_file:
            for entry in energies:
                energy_file.write(f"{int(entry[0])} {float(entry[1])}\n")

        # Write the structures to files
        write(f"unrlxd_structures_seed{seed}.traj", unrlxd_structures)
        write(f"rlxd_structures_seed{seed}.traj", rlxd_structures)
        print("All generated and relaxed structures written")

    print("Learning complete")
