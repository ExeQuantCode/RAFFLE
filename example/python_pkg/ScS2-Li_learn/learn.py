from chgnet.model import CHGNetCalculator
# from mace.calculators import mace_mp
from raffle.generator import raffle_generator
from ase import build, Atoms
from ase.optimize import BFGS, FIRE
from ase.io import write, read
import numpy as np
import os
from concurrent.futures import ProcessPoolExecutor
from joblib import Parallel, delayed


def process_structure(i, atoms, num_structures_old, calc_params, optimise_structure, iteration, calc):
    if i < num_structures_old:
        return
    
    # calc = Vasp(**calc_params, label=f"struct{i}", directory=f"iteration{iteration}/struct{i}/", txt=f"stdout{i}.o")
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
            optimizer.run(fmax=0.02, steps=300)
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

# crystal_structures = [
#     'sc', 'fcc', 'bcc', 'hcp',
#     'diamond', 'zincblende', 'rocksalt', 'cesiumchloride',
#     'fluorite', 'wurtzite', 'tetragonal', 'orthorhombic',
#     'bct', 'rhombohedral', 'mcl'
# ]



if __name__ == "__main__":

    calc_params = {}
    calc = CHGNetCalculator()
    # calc = mace_mp(model="medium", dispersion=False, default_dtype="float32", device='cpu')
    
    host = read("../POSCAR_host")
    host.calc = calc
    host_reference_energy = host.get_potential_energy() / len(host)

    hosts = []
    for a in [0.8, 0.9, 1.0, 1.1, 1.2]:
        for c in [1.0, 1.1, 1.2, 1.3, 1.4]:
            atom = host.copy()
            atom.set_cell([a * host.cell[0], a * host.cell[1], c * host.cell[2]], scale_atoms=True)
            hosts.append(atom)
            hosts[-1].calc = calc
            print(hosts[-1])

    Li_reference = Atoms(
        "Li2",
        positions=[[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]],
        cell=[3.42, 3.42, 3.42],
        pbc=[True, True, True],
    )
    Li_reference.calc = calc
    Li_reference_energy = Li_reference.get_potential_energy() / len(Li_reference)

    generator = raffle_generator()
    generator.distributions.set_element_energies(
        {
            'Sc': 0.0,
            'S': 0.0,
            'Li': Li_reference_energy,
        }
    )
    # set energy scale
    generator.distributions.set_kBT(0.4)
    # set the distribution function widths (2-body, 3-body, 4-body)
    generator.distributions.set_width([0.025, np.pi/200.0, np.pi/200.0])


    initial_database = [host, Li_reference]

    generator.distributions.create(initial_database)

    seed = 0
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

    iter = -1
    unrlxd_structures = []
    rlxd_structures = []
    num_structures_old = 0
    optimise_structure = True
    for host in hosts:
        print("setting host")
        generator.set_host(host)
        volume = host.get_volume()
        for num_atoms in range(1, 3):

            iter += 1
            print(f"Seed: {seed}")
            generator.generate(
                num_structures = 5,
                stoichiometry = { 'Li': num_atoms },
                seed = seed*1000+iter,
                method_probab = {"void": 1.0, "rand": 0.5, "walk": 1.0, "grow": 0.0, "min": 1.0},
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
            
            print("Starting parallel execution")
            results = Parallel(n_jobs=5)(
                delayed(process_structure)(i, generated_structures[i], num_structures_old, calc_params, optimise_structure, iteration=seed, calc=calc)
                for i in range(num_structures_old, num_structures_new)
            )


            # Wait for all futures to complete
            for j, result in enumerate(results):
                generated_structures[j+num_structures_old], energy_unrlxd[j], energy_rlxd[j] = result
                if generated_structures[j+num_structures_old] is None:
                    print("Structure failed the checks")
                    continue
                rlxd_structures.append(generated_structures[j+num_structures_old].copy())
            print("All futures completed")


            # Wait for all futures to complete
            for j, result in enumerate(results):
                generated_structures[j+num_structures_old], energy_unrlxd[j], energy_rlxd[j] = result
                if generated_structures[j+num_structures_old] is None:
                    print("Structure failed the checks")
                    continue
                rlxd_structures.append(generated_structures[j+num_structures_old].copy())
            print("All futures completed")

            # Remove structures that failed the checks
            for j, atoms in reversed(list(enumerate(generated_structures))):
                if atoms is None:
                    energy_unrlxd = np.delete(energy_unrlxd, j-num_structures_old)
                    energy_rlxd = np.delete(energy_rlxd, j-num_structures_old)
                    del generated_structures[j]
                    # del unrlxd_structures[j]
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

    write(f"unrlxd_structures_seed{seed}.traj", unrlxd_structures)
    write(f"rlxd_structures_seed{seed}.traj", rlxd_structures)
    print("All generated and relaxed structures written")

    print("Learning complete")