from chgnet.model import CHGNetCalculator
from raffle.generator import raffle_generator
from ase import build, Atoms
from ase.optimize import FIRE
from ase.io import write, read
import numpy as np
import os
from joblib import Parallel, delayed


# Function to relax a structure
def process_structure(i, atoms, num_structures_old, calc_params, optimise_structure, iteration):
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

    hosts = []
    # get all hosts in ../../data/C-MgO_hosts/ directory
    for filename in os.listdir("../../data/C-MgO_hosts/"):
        if filename.startswith("POSCAR"):
            atoms = read(f"../../data/C-MgO_hosts/{filename}")
            if not all(element == 'C' for element in set(atoms.get_chemical_symbols())):
                print(f"Skipping {filename} as it has non-carbon atoms")
                continue
            elif len(atoms) > 20:
                print(f"Skipping {filename} as it has more than 20 atoms (want to do smaller cells for now)")
                continue
            hosts.append(read(f"../../data/C-MgO_hosts/{filename}"))
            hosts[-1].calc = calc

    # set up element reference energies
    C_reference = Atoms('C4', positions=[
        [0.000000,    0.000000,    1.950768],
        [0.000000,    0.000000,    5.852305],
        [0.000000,    1.424491,    1.950768],
        [1.233646,    0.712246,    5.852305],
        ], cell=[[1.2336, -2.1367, 0.0], [1.2336, 2.1367, 0.0], [0.0, 0.0, 7.8031]], pbc=True)
    C_reference.calc = calc
    C_reference_energy = C_reference.get_potential_energy() / len(C_reference)
    Mg_reference = Atoms("Mg2", positions=[
            [0.000000,    1.831369,    1.285301],
            [1.586012,    0.915684,    3.855904],
        ], cell=[
            [1.58601,   -2.7471,    0.0000],
            [1.58601,    2.7471,    0.0000],
            [0.00000,    0.0001,    5.1412],
        ], pbc=True)
    Mg_reference.calc = calc
    Mg_reference_energy = Mg_reference.get_potential_energy() / len(Mg_reference)
    O_reference = Atoms("O2", positions=[
        [0.0, 0.0, 0.0],
        [1.23, 0.0, 0.0]], cell=[10, 10, 10], pbc=False)
    O_reference.calc = calc
    O_reference_energy = O_reference.get_potential_energy() / len(O_reference)

    generator = raffle_generator()
    generator.distributions.set_element_energies(
        {
            'C': C_reference_energy,
            'Mg': Mg_reference_energy,
            'O': O_reference_energy
        }
    )

    # set the covalent radii
    generator . distributions . set_bond_radii (
        {
            ('C', 'Mg'): 1.4,
            ('Mg', 'Mg'): 1.4,
            ('C', 'O'): 1.4,
            # ('C', 'Mg'): 1.15, # 2.3 = what Joe used https://github.com/ExeQuantCode/RAFFLE/blob/1cfdf2d220dce41f20ace4cce5a4135a37d1f1b1/chem.in
            # ('C', 'O'): 0.587,
            # ('Mg', 'Mg'): 1.592,
            # INITIALISE O2 BOND TO LARGER VALUE
        }
    )

    # set energy scale
    generator.distributions.set_kBT(0.4)

    # set the distribution function widths (2-body, 3-body, 4-body)
    generator.distributions.set_width([0.025, np.pi/200.0, np.pi/200.0])

    # generate MgO bulk
    MgO_atom = Atoms("Mg4O4", positions=[
            [0.0, 0.0, 0.0],
            [0.0, 2.097, 2.097],
            [2.097, 0.0, 2.097],
            [2.097, 2.097, 0.0],
            [0.0, 0.0, 2.097],
            [0.0, 2.097, 0.0],
            [2.097, 0.0, 0.0],
            [2.097, 2.097, 2.097],
        ], cell=[4.1940, 4.1940, 4.1940], pbc=True)
    MgO_atom.calc = calc

    # read C-Mg-O database, obtained from Materials Project
    print("Reading database")
    database = read("../database.xyz", index=":")

    # set up the initial database
    initial_database = []
    for i, atoms in enumerate(database):
        if atoms.calc is None:
            database.remove(atoms)
            continue
        # check if any other species in the system other than C, Mg, O
        if not all(element in ['C', 'Mg', 'O'] for element in set(atoms.get_chemical_symbols())):
            print(f"Skipping structure {i} as it has non-C-Mg-O atoms")
            continue
        atoms.calc = calc
        initial_database.append(atoms)
    initial_database.append(C_reference)
    initial_database.append(Mg_reference)
    initial_database.append(O_reference)
    initial_database.append(MgO_atom)
    generator.distributions.create(initial_database, deallocate_systems=False)

    # set parameters
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
    num_structures_old = 0
    unrlxd_structures = []
    rlxd_structures = []
    optimise_structure = True
    mass_C = 12.011
    mass_Mg = 24.305
    mass_O = 15.999
    # loop over all hosts
    for host in hosts:
        print("setting host")
        generator.set_host(host)
        volume = host.get_volume()
        num_atoms_C = len(host)
    
        host_positions = host.get_scaled_positions(wrap=True)
        host_height = np.linalg.norm(host.get_cell()[2, :])
        print(f"host_height: {host_height}")

        # get all unique z values
        z_values = np.unique(host_positions[:, 2])
        print(f"z_values: {z_values}")
        # return the two closes to 0.5 in the cell
        z_values = z_values[np.argsort(np.abs(z_values - 0.5))[:2]]
        print(f"z_values: {z_values}")
        # set the bounds for the generator
        z_min = min(z_values) + 1.0 / host_height
        z_max = max(z_values) - 1.0 / host_height
        generator.set_bounds([[0, 0, z_min], [1, 1, z_max]])

        # loop over range of MgO densities to intercalate
        for num_atoms in range(1, 50):

            num_atoms_Mg = num_atoms
            num_atoms_O = num_atoms
            density = ( mass_C * num_atoms_C + mass_Mg * num_atoms_Mg + mass_O * num_atoms_O ) / volume
            if density > 2.2:
                # print("Density too high:", density, "u/A^3")
                continue
            elif density < 1.0:
                # print("Density too low", density, "u/A^3")
                continue
            iter += 1

            # generate structures
            generator.generate(
                num_structures = 5,
                stoichiometry = { 'Mg': num_atoms_Mg, 'O': num_atoms_O },
                seed = seed*1000+iter,
                method_ratio = {"void": 1.0, "rand": 0.001, "walk": 1.0, "grow": 0.0, "min": 1.0},
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
            
            # Start parallel execution
            print("Starting parallel execution")
            results = Parallel(n_jobs=5)(
                delayed(process_structure)(i, generated_structures[i].copy(), num_structures_old, calc_params, optimise_structure, iteration=seed)
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