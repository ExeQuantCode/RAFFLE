# from mace.calculators import mace_mp
# from ase.calculators.vasp import Vasp
from chgnet.model import CHGNetCalculator
from ase.calculators.singlepoint import SinglePointCalculator
from raffle.generator import raffle_generator
from ase import build, Atoms
from ase.optimize import FIRE
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
    method_dicts = [
        {"void":  1.0, "rand":  0.0, "walk":  0.0, "grow":  0.0, "min":  0.0},
        {"void":  0.0, "rand":  1.0, "walk":  0.0, "grow":  0.0, "min":  0.0},
        {"void":  0.0, "rand":  0.0, "walk":  1.0, "grow":  0.0, "min":  0.0},
        {"void":  0.0, "rand":  0.0, "walk":  0.0, "grow":  1.0, "min":  0.0},
        {"void":  0.0, "rand":  0.0, "walk":  0.0, "grow":  0.0, "min":  1.0},
        {"void": 85.0, "rand": 15.0, "walk":  0.0, "grow":  0.0, "min":  0.0},
        {"void":  0.0, "rand": 15.0, "walk": 85.0, "grow":  0.0, "min":  0.0},
        {"void":  0.0, "rand": 15.0, "walk":  0.0, "grow": 85.0, "min":  0.0},
        {"void":  0.0, "rand": 15.0, "walk":  0.0, "grow":  0.0, "min": 85.0},
    ]

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

    # get the current directory
    cwd = os.getcwd()
    print(f"Current working directory: {cwd}")

    seed = 0
    for methods in method_dicts:
        # move to the cwd
        os.chdir(cwd)
        print(f"Changed directory to {cwd}")

        # print the methods
        print(f"Method: {methods}")

        # make method directory
        method_dir = f"method_v{methods['void']}_r{methods['rand']}_w{methods['walk']}_g{methods['grow']}_m{methods['min']}/"
        if not os.path.exists(method_dir):
            os.makedirs(method_dir)
        os.chdir(method_dir)
        print(f"Changed directory to {method_dir}")
        
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

        # delete .traj files if they exist
        if os.path.exists(f"unrlxd_structures_seed{seed}.xyz"):
            os.remove(f"unrlxd_structures_seed{seed}.xyz")
        if os.path.exists(f"rlxd_structures_seed{seed}.xyz"):
            os.remove(f"rlxd_structures_seed{seed}.xyz")

        num_structures_tot = 0
        unrlxd_structures = []
        rlxd_structures = []
        iter2 = -1
        generator.init_seed(put=seed)
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
    
                generated_structures = generator.generate(
                    num_structures = 5,
                    stoichiometry = { 'Al': num_atoms },
                    method_ratio = methods,
                    verbose = 0,
                    calc = calc,
                )

                # print the number of structures generated
                print("Total number of structures generated: ", generator.num_structures)
                num_structures_new = len(generated_structures)

                # check if directory iteration[iter] exists, if not create it
                iterdir = f"iteration{iter2}/"
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

            
                # Start parallel execution
                print("Starting parallel execution")
                results = Parallel(n_jobs=5)(
                    delayed(process_structure)(i, generated_structures[i].copy(), calc_params, optimise_structure, calc=calc)
                    for i in range(len(generated_structures))
                )

                # Wait for all futures to complete
                energies_unrlxd = []
                energies_rlxd = []
                for j, result in enumerate(results):
                    rlxd_generated_structure, energy_unrlxd, energy_rlxd = result
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
                    write(iterdir+f"POSCAR_{i}", rlxd_structures[num_structures_tot + i])
                    print(f"Structure {i} energy per atom: {energies_rlxd[i]}")
                    # append energy per atom to the 'energies_unrlxd_filename' file
                    with open(energies_unrlxd_filename, "a") as energy_file:
                        energy_file.write(f"{i+num_structures_tot} {energies_unrlxd[i]}\n")
                    # append energy per atom to the 'energies_rlxd_filename' file
                    with open(energies_rlxd_filename, "a") as energy_file:
                        energy_file.write(f"{i+num_structures_tot} {energies_rlxd[i]}\n")

                # update the distribution functions
                print("Updating distributions")
                generator.distributions.update(rlxd_structures[num_structures_tot:], from_host=False, deallocate_systems=False)

                # print the new distribution functions to a file
                print("Printing distributions")
                generator.distributions.write_dfs(iterdir+"distributions.txt")
                generator.distributions.write_2body(iterdir+"df2.txt")
                generator.distributions.write_3body(iterdir+"df3.txt")
                generator.distributions.write_4body(iterdir+"df4.txt")
                generator.distributions.deallocate_systems()

                write(f"unrlxd_structures_seed{seed}.xyz", unrlxd_structures[num_structures_tot:], append=True)
                write(f"rlxd_structures_seed{seed}.xyz", rlxd_structures[num_structures_tot:], append=True)

                # update the number of structures generated
                num_structures_tot += num_structures_new


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