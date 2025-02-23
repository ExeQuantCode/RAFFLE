# from mace.calculators import mace_mp
# from ase.calculators.vasp import Vasp
from chgnet.model import CHGNetCalculator
from raffle.generator import raffle_generator
from ase import build, Atoms
from ase.optimize import BFGS, FIRE
from ase.io import write
from ase.visualize import view
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


def process_structure_with_queue(i, structure, num_old, calc_params, optimise_structure, iteration, queue):
    # Perform the computation
    result = process_structure(i, structure, num_old, calc_params, optimise_structure, iteration)
    queue.put(result)  # Report completion

def process_structure(i, atoms, num_structures_old, calc_params, optimise_structure, iteration):
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

# crystal_structures = [
#     'sc', 'fcc', 'bcc', 'hcp',
#     'diamond', 'zincblende', 'rocksalt', 'cesiumchloride',
#     'fluorite', 'wurtzite', 'tetragonal', 'orthorhombic',
#     'bct', 'rhombohedral', 'mcl'
# ]



if __name__ == "__main__":

    calc_params = {}
    calc = CHGNetCalculator()
    # calc_params = {
    #     "model": "large",
    #     "dispersion": False,
    #     "default_dtype": "float32",
    #     "device": 'cpu'
    # }
    # calc = mace_mp(**calc_params)
    # calc_params = {
    #     "command": "$HOME/DVASP/vasp.6.4.3/bin/vasp_std",
    #     # "label": "iteration",
    #     # "txt": "stdout.o",
    #     "xc": 'pbe',
    #     "setups": {'C': ''},
    #     "kpts": (3, 3, 3),
    #     "encut": 400,
    #     "istart": 0,
    #     "icharg": 0,
    # }
    # ## DON'T FORGET TO EXPORT THE VASP_PP_PATH
    # calc = Vasp(**calc_params, label="tmp", directory="tmp", txt="stdout.o")

    #crystal_structures = [
    #    'orthorhombic', 'diamond',
    #    'bct', 'sc',
    #    'fcc', 'bcc', 'hcp',
    #]

    """
    This needs to be changed to set the material cells that will be searched through.
    The crystal_structures list contains the crystal structures that will be used to generate the hosts.
    """
    crystal_structures = [
        'orthorhombic', 'hcp',
    ]

    """
    This needs to be changed to set the lattice constants that will be searched through.
    The lattice_constants list contains the lattice constants that will be used to generate the hosts.
    """
    lattice_constants = np.linspace(3.1, 5.4, num=6)

    """
    This needs to be changed to set the values of the method_ratio.
    """
    #void_val; rand_val; walk_val; grow_val; min_val
    ##############[  V,   R,   W,   G,   M]    
    method_val = [85.0, 15.0, 0.0, 0.0, 0.0] #change



    hosts = []
    for crystal_structure in crystal_structures:
        print(f'Crystal structure: {crystal_structure}')
        for a in lattice_constants:
            b = a
            c = a
            atom = build.bulk(
                    name = 'Al',
                    crystalstructure = crystal_structure,
                    a = a,
                    b = b,
                    c = c,
            )
            hosts.append(Atoms('Al', positions=[(0, 0, 0)], cell=atom.get_cell(), pbc=True, calculator=calc))
            # hosts[-1].set_pbc(True)
            # hosts[-1].calc = calc
            print(hosts[-1])

    print("number of hosts: ", len(hosts))

    optimise_structure = True
    mass = 26.9815385
    density = 1.61 # u/A^3
    # num_atoms = 7

    for seed in range(1):
        print(f"Seed: {seed}")
        energies_rlxd_filename = f"energies_rlxd_seed{seed}.txt"
        energies_unrlxd_filename = f"energies_unrlxd_seed{seed}.txt"
        generator = raffle_generator()
        generator.distributions.set_element_energies(
            {
                'Al': 0.0
            }
        )
        # set energy scale
        generator.distributions.set_kBT(0.4)
        # set the distribution function widths (2-body, 3-body, 4-body)
        generator.distributions.set_width([0.025, np.pi/200.0, np.pi/200.0])

        initial_database = [Atoms('Al', positions=[(0, 0, 0)], cell=[8, 8, 8], pbc=True)]
        initial_database[0].calc = calc

        generator.distributions.create(initial_database)

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


        num_structures_old = 0
        unrlxd_structures = []
        rlxd_structures = []
        iter2 = 0
        for iter in range(10):
            for host in hosts:
                generator.set_host(host)
                volume = host.get_volume()

                num_atoms = round(density * volume / mass) - 1
                if(num_atoms < 1):
                    continue
                iter2 += 1
                print(f"Volume: {volume}")
                print(f"Number of atoms: {num_atoms}")
    
                generator.generate(
                    num_structures = 5,
                    stoichiometry = { 'Al': num_atoms },
                    seed = seed*1000+iter,
                    method_ratio = {"void": method_val[0], "rand": method_val[1], "walk": method_val[2], "grow": method_val[3], "min": method_val[4]},
                    verbose = 0,
                )

                # print the number of structures generated
                print("Total number of structures generated: ", generator.num_structures)
                generated_structures = generator.get_structures(calc)
                num_structures_new = len(generated_structures)

                # check if directory iteration[iter] exists, if not create it
                iterdir = f"iteration{iter2}/"
                if not os.path.exists(iterdir):
                    os.makedirs(iterdir)
                generator.print_settings(iterdir+"generator_settings.txt")

                # set up list of energies
                energy_unrlxd = np.zeros(num_structures_new - num_structures_old)
                energy_rlxd = np.zeros(num_structures_new - num_structures_old)
                for i in range(num_structures_new - num_structures_old):
                    write(iterdir+f"POSCAR_unrlxd_{i}", generated_structures[num_structures_old + i])
                    print(f"Structure {i} energy per atom: {generated_structures[num_structures_old + i].get_potential_energy() / len(generated_structures[num_structures_old + i])}")
                    unrlxd_structures.append(deepcopy(generated_structures[num_structures_old + i]))
                
                # Start parallel execution
                print("Starting parallel execution")
                results = Parallel(n_jobs=5)(
                    delayed(process_structure)(i, deepcopy(generated_structures[i]), num_structures_old, calc_params, optimise_structure, iteration=seed)
                    for i in range(num_structures_old, num_structures_new)
                )

                # Wait for all futures to complete
                for j, result in enumerate(results):
                    generated_structures[j+num_structures_old], energy_unrlxd[j], energy_rlxd[j] = result
                    rlxd_structures.append(deepcopy(generated_structures[j+num_structures_old]))
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