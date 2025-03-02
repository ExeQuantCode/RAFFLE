from mace.calculators import mace_mp
from ase.calculators.singlepoint import SinglePointCalculator
from raffle.generator import raffle_generator
from ase import build
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
    calc_params = { 'model': "../mace-mpa-0-medium.model" }
    calc = mace_mp(**calc_params)

    # set up the hosts
    hosts = []

    Si_bulk = build.bulk("Si", crystalstructure="diamond", a=5.43)
    Si_bulk.calc = calc
    Si_reference_energy = Si_bulk.get_potential_energy() / len(Si_bulk)
    Si_cubic = build.make_supercell(Si_bulk, [[-1, 1, 1], [1, -1, 1], [1, 1, -1]])
    Ge_bulk = build.bulk("Ge", crystalstructure="diamond", a=5.65)
    Ge_bulk.calc = calc
    Ge_cubic = build.make_supercell(Ge_bulk, [[-1, 1, 1], [1, -1, 1], [1, 1, -1]])
    Ge_reference_energy = Ge_bulk.get_potential_energy() / len(Ge_bulk)

    Si_supercell = build.make_supercell(Si_cubic, [[2, 0, 0], [0, 2, 0], [0, 0, 1]])
    Ge_supercell = build.make_supercell(Ge_cubic, [[2, 0, 0], [0, 2, 0], [0, 0, 1]])

    Si_surface = build.surface(Si_supercell, indices=(0, 0, 1), layers=2)
    Ge_surface = build.surface(Ge_supercell, indices=(0, 0, 1), layers=2)

    Si_slab = build.surface(Si_supercell, indices=(0, 0, 1), layers=2, vacuum=12, periodic=True)
    Si_slab.calc = calc
    Ge_slab = build.surface(Ge_supercell, indices=(0, 0, 1), layers=2, vacuum=12, periodic=True)
    Ge_slab.calc = calc

    host = build.stack(Si_surface, Ge_surface, axis=2, distance= 5.43/2 + 5.65/2)
    cell = host.get_cell()
    print("Host cell: ", host.get_cell())
    cell[2, 2] -= 3.8865 # (5.43 + 5.65) / 2 * 3/4
    host.set_cell(cell, scale_atoms=False)

    hosts.append(host)

    # set the parameters for the generator
    generator = raffle_generator()
    generator.distributions.set_element_energies(
        {
            'Si': Si_reference_energy,
            'Ge': Ge_reference_energy,
        }
    )

    # set energy scale
    generator.distributions.set_kBT(0.2)

    # set the distribution function widths (2-body, 3-body, 4-body)
    generator.distributions.set_width([0.04, np.pi/160.0, np.pi/160.0])

    # set the initial database
    initial_database = [Si_bulk, Ge_bulk]
    generator.distributions.create(initial_database, deallocate_systems=False)

    # set the parameters for the generator
    seed = 0

    # check if the energies file exists, if not create it
    energies_rlxd_filename = f"energies_rlxd_seed{seed}.txt"
    energies_unrlxd_filename = f"energies_unrlxd_seed{seed}.txt"
    
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
    iter = -1
    unrlxd_structures = []
    rlxd_structures = []
    num_structures_old = 0
    optimise_structure = True
    # start the iterations, loop over the hosts, then repeat the process X times
    for ival in range(40):
        for host in hosts:
            print("setting host")
            generator.set_host(host)
            generator.set_bounds([[0, 0, 0.34], [1, 1, 0.52]])

            iter += 1
            print(f"Iteration {iter}")

            # generate the structures
            generator.generate(
                num_structures = 5,
                stoichiometry = { 'Si': 16, 'Ge': 16 },
                seed = seed*1000+iter,
                method_ratio = {"void": 0.1, "rand": 0.01, "walk": 0.25, "grow": 0.25, "min": 1.0},
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
                # print(f"Structure {i} energy per unit area: {( generated_structures[num_structures_old + i].get_potential_energy() - Si_slab.get_potential_energy() - Ge_slab.get_potential_energy() ) / (2 * ( cell[0, 0] * cell[1, 1] ) )}")
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
                if j < num_structures_old:
                    continue
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