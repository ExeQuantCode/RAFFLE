import unittest
import os
from raffle.raffle import Geom_Rw, Raffle__Distribs_Container, Generator
from ase import Atoms

class TestGeomRw(unittest.TestCase):

    def test_species_type_initialization(self):
        species = Geom_Rw.species_type()
        self.assertIsNotNone(species._handle)

    def test_species_type_properties(self):
        species = Geom_Rw.species_type()
        species.mass = 12.0
        self.assertEqual(species.mass, 12.0)
        species.charge = 1.0
        self.assertEqual(species.charge, 1.0)
        species.radius = 0.5
        self.assertEqual(species.radius, 0.5)
        species.name = "C"
        self.assertEqual(species.name, b"C")
        species.num = 2
        self.assertEqual(species.num, 2)

    def test_basis_initialization(self):
        basis = Geom_Rw.basis()
        self.assertIsNotNone(basis._handle)

    def test_basis_properties(self):
        basis = Geom_Rw.basis()
        basis.nspec = 2
        self.assertEqual(basis.nspec, 2)
        basis.natom = 4
        self.assertEqual(basis.natom, 4)
        basis.energy = -10.0
        self.assertEqual(basis.energy, -10.0)
        basis.lat = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        self.assertTrue((basis.lat == [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]).all())
        basis.lcart = True
        self.assertTrue(basis.lcart)
        basis.pbc = [True, True, True]
        self.assertTrue((basis.pbc == [True, True, True]).all())
        basis.sysname = "Test"
        self.assertEqual(basis.sysname, b"Test")

    def test_basis_array_initialization(self):
        basis_array = Geom_Rw.basis_array()
        self.assertIsNotNone(basis_array._handle)

    def test_basis_array_allocate(self):
        basis_array = Geom_Rw.basis_array()
        basis_array.allocate(2)
        self.assertEqual(len(basis_array.items), 2)

class TestRaffleDistribsContainer(unittest.TestCase):

    def test_distribs_container_initialization(self):
        distribs_container = Raffle__Distribs_Container.distribs_container_type()
        self.assertIsNotNone(distribs_container._handle)

    def test_distribs_container_set_kBT(self):
        distribs_container = Raffle__Distribs_Container.distribs_container_type()
        distribs_container.set_kBT(0.025)
        self.assertTrue( abs(distribs_container.kBT - 0.025) < 1e-6 )

    def test_distribs_container_set_weight_method(self):
        distribs_container = Raffle__Distribs_Container.distribs_container_type()
        distribs_container.set_weight_method('formation_energy')
        self.assertFalse(distribs_container.weight_by_hull)
        distribs_container.set_weight_method('energy_above_hull')
        self.assertTrue(distribs_container.weight_by_hull)

    def test_distribs_container_set_width(self):
        distribs_container = Raffle__Distribs_Container.distribs_container_type()
        distribs_container.set_width([0.1, 0.2, 0.3])
        self.assertTrue( ( abs( distribs_container.width - [0.1, 0.2, 0.3] ) < 1e-6 ).all() )

    def test_distribs_container_set_sigma(self):
        distribs_container = Raffle__Distribs_Container.distribs_container_type()
        distribs_container.set_sigma([0.1, 0.2, 0.3])
        self.assertTrue( ( abs( distribs_container.sigma - [0.1, 0.2, 0.3] ) < 1e-6 ).all() )

    def test_distribs_container_set_cutoff_min(self):
        distribs_container = Raffle__Distribs_Container.distribs_container_type()
        distribs_container.set_cutoff_min([0.1, 0.2, 0.3])
        self.assertTrue( ( abs( distribs_container.cutoff_min - [0.1, 0.2, 0.3] ) < 1e-6 ).all() )

    def test_distribs_container_set_cutoff_max(self):
        distribs_container = Raffle__Distribs_Container.distribs_container_type()
        distribs_container.set_cutoff_max([0.1, 0.2, 0.3])
        self.assertTrue( ( abs( distribs_container.cutoff_max - [0.1, 0.2, 0.3] ) < 1e-6 ).all() )

    def test_distribs_container_set_radius_distance_tol(self):
        distribs_container = Raffle__Distribs_Container.distribs_container_type()
        distribs_container.set_radius_distance_tol([0.1, 0.2, 0.3, 0.4])
        self.assertTrue( ( abs(distribs_container.radius_distance_tol - [0.1, 0.2, 0.3, 0.4]) < 1e-6 ).all() )

class TestGenerator(unittest.TestCase):

    def test_generator_initialization(self):
        generator = Generator.raffle_generator()
        self.assertIsNotNone(generator._handle)

    def test_set_max_attempts(self):
        generator = Generator.raffle_generator()
        generator.set_max_attempts(10)
        self.assertEqual(generator.max_attempts, 10)

    def test_set_walk_step_size(self):
        generator = Generator.raffle_generator()
        generator.set_walk_step_size(coarse=0.5, fine=0.1)
        self.assertTrue( abs( generator.walk_step_size_coarse - 0.5 ) < 1e-6 )
        self.assertTrue( abs( generator.walk_step_size_fine - 0.1 ) < 1e-6 )

    def test_set_host(self):
        generator = Generator.raffle_generator()
        atoms = Atoms('H2', positions=[[0, 0, 0], [0, 0, 0.74]])
        generator.set_host(atoms)
        self.assertIsNotNone(generator.host)

    def test_set_grid(self):
        generator = Generator.raffle_generator()
        generator.set_grid(grid=[10, 10, 10], grid_offset=[0.1, 0.2, 0.3])
        self.assertTrue( ( generator.grid == [10, 10, 10] ).all() )
        self.assertTrue( ( abs( generator.grid_offset - [0.1, 0.2, 0.3] ) < 1e-6 ).all() )

        generator.set_grid(grid_spacing=0.1)
        self.assertTrue( abs( generator.grid_spacing - 0.1 ) < 1e-6 )

    def test_reset_grid(self):
        generator = Generator.raffle_generator()
        generator.set_grid(grid=[10, 10, 10], grid_offset=[0.0, 0.0, 0.0])
        generator.reset_grid()
        self.assertTrue((generator.grid == [0, 0, 0]).all())

    def test_set_bounds(self):
        generator = Generator.raffle_generator()
        bounds = [ [ 0.1, 0.2, 0.3 ], [0.4, 0.5, 0.6 ] ]
        generator.set_bounds(bounds)
        self.assertTrue( ( abs( generator.bounds - bounds ) < 1e-6 ).all() )

    def test_reset_bounds(self):
        generator = Generator.raffle_generator()
        bounds = [ [ 0.1, 0.2, 0.3 ], [0.4, 0.5, 0.6 ] ]
        generator.set_bounds(bounds)
        generator.reset_bounds()
        self.assertTrue( ( abs( generator.bounds - [ [ 0.0, 0.0, 0.0 ], [ 1.0, 1.0, 1.0 ] ] ) < 1e-6 ).all() )

    # def test_generate(self):
    #     generator = Generator.raffle_generator()
    #     atoms = Atoms('H2', positions=[[0, 0, 0], [0, 0, 0.74]], cell=[10, 10, 10])
    #     generator.set_host(atoms)
    #     generator.distributions.set_element_energies( {"H": 0.0} )
    #     stoichiometry = {"H": 2}
    #     generator.generate(num_structures=1, stoichiometry=stoichiometry)
    #     self.assertEqual(generator.num_structures, 1)

    def test_get_structures(self):
        generator = Generator.raffle_generator()
        atoms = Atoms('H2', positions=[[0, 0, 0], [0, 0, 0.74]], cell=[10, 10, 10])
        generator.set_structures([atoms])
        structures = generator.get_structures()
        self.assertEqual( len(structures), 1 )

    def test_set_structures(self):
        generator = Generator.raffle_generator()
        atoms = Atoms('H2', positions=[[0, 0, 0], [0, 0, 0.74]], cell=[10, 10, 10])
        generator.set_structures([atoms])
        self.assertEqual(generator.num_structures, 1)

    def test_remove_structure(self):
        generator = Generator.raffle_generator()
        atoms = []
        atoms.append(Atoms('H2', positions=[[0, 0, 0], [0, 0, 0.74]], cell=[10, 10, 10]))
        atoms.append(Atoms('H2', positions=[[0, 0, 0], [0, 0, 0.74]], cell=[10, 10, 10]))
        generator.set_structures(atoms)
        generator.remove_structure(0)
        self.assertEqual(generator.num_structures, 1)

    # def test_evaluate(self):
    #     generator = Generator.raffle_generator()
    #     atoms = Atoms('H2', positions=[[0, 0, 0], [0, 0, 0.74]], cell=[10, 10, 10])
    #     generator.set_host(atoms)
    #     generator.distributions.set_element_energies( {"H": 0.0} )
    #     generator.distributions.create([atoms])
    #     basis = Geom_Rw.basis(atoms)
    #     viability = generator.evaluate(basis)
    #     self.assertIsNotNone(viability)

    def test_print_settings(self):
        generator = Generator.raffle_generator()
        atoms = Atoms('H2', positions=[[0, 0, 0], [0, 0, 0.74]], cell=[10, 10, 10])
        generator.set_host(atoms)
        generator.distributions.set_element_energies( {"H": 0.0} )
        generator.distributions.create([atoms])

        # check if filename exists
        i = 0
        filename = ".raffle_generator_settings.txt"
        while os.path.isfile(filename):
            i += 1
            filename = ".raffle_generator_settings"+str(i)+".txt"

        generator.print_settings(filename)
        with open(filename, "r") as file:
            settings = file.read()
        self.assertTrue(len(settings) > 0)
        os.remove(filename)


if __name__ == '__main__':
    unittest.main()