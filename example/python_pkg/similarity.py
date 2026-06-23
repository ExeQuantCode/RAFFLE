from ase.build import bulk, make_supercell
from ase.io import read
import numpy as np


from raffle import structure_similarity_rmsd


def main():

    # Example 1:
    # Primitive vs supercell comparison

    si_primitive = bulk("Si", cubic=True)

    P = np.diag([2, 2, 2])

    si_supercell = make_supercell(
        si_primitive,
        P,
    )

    rmsd, translation = structure_similarity_rmsd(
        si_primitive,
        si_supercell,
        rcut=6.0,
        return_translation=True,
    )

    print("Primitive vs supercell")
    print("----------------------")
    print(f"RMSD: {rmsd:.8f} Å")
    print(f"Translation: {translation}")

    # Example 2:
    # Compare structures loaded from files

    initial = read("initial.xyz")
    final = read("final.xyz")
    target = read("target.xyz")

    initial_rmsd, initial_vector = structure_similarity_rmsd(
        target,
        initial,
        rcut=6.0,
        return_translation=True,
    )

    final_rmsd, final_vector = structure_similarity_rmsd(
        target,
        final,
        rcut=6.0,
        return_translation=True,
    )



    print("\nFile-based comparison")
    print("---------------------")
    print(f"Initial RMSD: {initial_rmsd:.8f} Å")
    print(f"Final RMSD: {final_rmsd:.8f} Å")

    print(f"Initial translation: {initial_vector}")
    print(f"Final translation: {final_vector}")


if __name__ == "__main__":
    main()
