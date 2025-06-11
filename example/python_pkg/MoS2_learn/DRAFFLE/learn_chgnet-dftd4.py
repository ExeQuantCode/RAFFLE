from dftd4.ase import DFTD4
from chgnet.model import CHGNetCalculator
from ase.calculators.mixing import SumCalculator
from ase.calculators.singlepoint import SinglePointCalculator
from raffle.generator import raffle_generator
from ase import build, Atoms
from ase.optimize import FIRE
from ase.io import write, read
import numpy as np
import os
from joblib import Parallel, delayed
from pathlib import Path

script_dir = Path(__file__).resolve().parent

import logging
logging.basicConfig(level=logging.DEBUG)

# Function to relax a structure
def process_structure(i, atoms, calc_params, optimise_structure, calc):
    # Check if the structure has already been processed
    
    # calc = Vasp(**calc_params, label=f"struct{i}", directory=f"iteration{iteration}/struct{i}/", txt=f"stdout{i}.o")
    atoms.calc = calc
    # positions_initial = atoms.get_positions()

    # Calculate and save the initial energy per atom
    energy_unrlxd = atoms.get_potential_energy() / len(atoms)
    print(f"Initial energy per atom: {energy_unrlxd}")

    # Optimise the structure if requested
    if optimise_structure:
        optimizer = FIRE(atoms, trajectory = f"traje{i}.traj", logfile=f"optimisation{i}.log")
        try:
            optimizer.run(fmax=0.02, steps=200)
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
    calc = SumCalculator( [ DFTD4(method='PBE'), CHGNetCalculator() ] )

    # set up the hosts
    hosts = []
    a = 3.19
    hosts.append(Atoms('Mo', positions=[(0, 0, 0)], cell=[
        [a * 0.5, a * -np.sqrt(3.0)/2.0, 0.0],
        [a * 0.5, a * np.sqrt(3.0)/2.0, 0.0],
        [0.0, 0.0, 13.1]
        ], pbc=True, calculator=calc))

    # set the parameters for the generator
    optimise_structure = True

    # loop over random seeds
    for seed in range(1):
        print(f"Seed: {seed}")
        energies_rlxd_filename = f"energies_rlxd_seed{seed}.txt"
        energies_unrlxd_filename = f"energies_unrlxd_seed{seed}.txt"
        generator = raffle_generator()
        generator.distributions.set_element_energies(
            {
                'Mo': -23.00869551/2,
                'S': -525.85160941/128
            }
        )
        # set energy scale
        generator.distributions.set_kBT(0.4)
        # set the distribution function widths (2-body, 3-body, 4-body)
        generator.distributions.set_width([0.025, np.pi/200.0, np.pi/200.0])

        # set the initial database
        Mo_bulk = read(script_dir/ ".." / "Mo.poscar")
        Mo_bulk.calc = calc
        S_bulk = read(script_dir/ ".." / "S.poscar")
        S_bulk.calc = calc
        initial_database = [Mo_bulk, S_bulk]
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
        generator.init_seed(put=seed)
        # start the iterations
        for iter in range(200):
            # generate the structures
            generated_structures = generator.generate(
                num_structures = 5,
                stoichiometry = { 'Mo': 1 , 'S': 4 },
                method_ratio = {"void": 0.1, "rand": 0.01, "walk": 0.25, "grow": 0.25, "min": 1.0},
                verbose = 0,
                calc = calc,
            )

            # print the number of structures generated
            print("Total number of structures generated: ", generator.num_structures)

            # check if directory iteration[iter] exists, if not create it
            iterdir = f"iteration{iter}/"
            if not os.path.exists(iterdir):
                os.makedirs(iterdir)
            generator.print_settings(iterdir+"generator_settings.txt")

            # set up list of energies
            for i in reversed(range(len(generated_structures))):
                if generated_structures[i] is None:
                    print(f"Structure {i} is None, skipping")
                    del generated_structures[i]
                    continue
                write(iterdir+f"POSCAR_unrlxd_{i}", generated_structures[i])
                print(f"Structure {i} energy per atom: {generated_structures[i].get_potential_energy() / len(generated_structures[i])}")

            
            # Wait for all futures to complete
            energies_unrlxd = []
            energies_rlxd = []
            for j in range(len(generated_structures)):
                rlxd_generated_structure, energy_unrlxd, energy_rlxd = process_structure(
                    j, generated_structures[j].copy(), calc_params, optimise_structure, calc=calc
                )
                if rlxd_generated_structure is None:
                    print("Structure failed the checks")
                    continue
                else:
                    unrlxd_structures.append(generated_structures[j].copy())
                    unrlxd_structures[-1].calc = SinglePointCalculator(
                        generated_structures[j],
                        energy=generated_structures[j].get_potential_energy(),
                        forces=generated_structures[j].get_forces()
                    )
                    rlxd_structures.append(rlxd_generated_structure.copy())
                    rlxd_structures[-1].calc = SinglePointCalculator(
                        rlxd_generated_structure,
                        energy=rlxd_generated_structure.get_potential_energy(),
                        forces=rlxd_generated_structure.get_forces()
                    )
                    energies_unrlxd.append(energy_unrlxd)
                    energies_rlxd.append(energy_rlxd)
            print("All futures completed")
            num_structures_new = len(energies_unrlxd)


            # write the structures to files
            for i in range(num_structures_new):
                write(iterdir+f"POSCAR_{i}", rlxd_structures[num_structures_old + i])
                print(f"Structure {i} energy per atom: {energies_rlxd[i]}")
                # append energy per atom to the 'energies_unrlxd_filename' file
                with open(energies_unrlxd_filename, "a") as energy_file:
                    energy_file.write(f"{i+num_structures_old} {energies_unrlxd[i]}\n")
                # append energy per atom to the 'energies_rlxd_filename' file
                with open(energies_rlxd_filename, "a") as energy_file:
                    energy_file.write(f"{i+num_structures_old} {energies_rlxd[i]}\n")

            # update the distribution functions
            print("Updating distributions")
            generator.distributions.update(generated_structures, from_host=False, deallocate_systems=False)

            # print the new distribution functions to a file
            print("Printing distributions")
            generator.distributions.write_dfs(iterdir+"distributions.txt")
            generator.distributions.write_2body(iterdir+"df2.txt")
            generator.distributions.write_3body(iterdir+"df3.txt")
            generator.distributions.write_4body(iterdir+"df4.txt")
            generator.distributions.deallocate_systems()

            # update the number of structures generated
            num_structures_old += num_structures_new

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