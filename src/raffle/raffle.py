from __future__ import print_function, absolute_import, division
import raffle._raffle as _raffle
import f90wrap.runtime
import logging
import numpy
import warnings


class Geom_Rw(f90wrap.runtime.FortranModule):
    """
    Code for handling geometry read/write operations.

    This module provides the necessary functionality to read, write, and
    store atomic geometries.
    In this module, and all of the codebase, element and species are used
    interchangeably.

    Defined in ../src/lib/mod_geom_rw.f90

    .. note::
        It is recommended not to use this module directly, but to handle
        atom objects through the ASE interface.
        This is provided mostly for compatibility with the existing codebase
        and Fortran code.
    """
    @f90wrap.runtime.register_class("raffle.species_type")
    class species_type(f90wrap.runtime.FortranDerivedType):
        def __init__(self, handle=None):
            """
            Create a ``species_type`` object.

            Returns:
                species (species_type):
                    Object to be constructed
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _raffle.f90wrap_geom_rw__species_type_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result

        def __del__(self):
            """
            Destructor for class species_type


            Defined at ../src/lib/mod_geom_rw.f90 lines \
                26-32

            Parameters
            ----------
            this : species_type
            	Object to be destructed


            Automatically generated destructor for species_type
            """
            if self._alloc:
                _raffle.f90wrap_geom_rw__species_type_finalise(this=self._handle)

        @property
        def atom(self):
            """
            Derived type containing the atomic information of a crystal.
            """
            array_ndim, array_type, array_shape, array_handle = \
                _raffle.f90wrap_species_type__array__atom(self._handle)
            if array_handle in self._arrays:
                atom = self._arrays[array_handle]
            else:
                atom = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _raffle.f90wrap_species_type__array__atom)
                self._arrays[array_handle] = atom
            return atom

        @atom.setter
        def atom(self, atom):
            self.atom[...] = atom

        @property
        def mass(self):
            """
            The mass of the element.
            """
            return _raffle.f90wrap_species_type__get__mass(self._handle)

        @mass.setter
        def mass(self, mass):
            _raffle.f90wrap_species_type__set__mass(self._handle, mass)

        @property
        def charge(self):
            """
            The charge of the element.
            """
            return _raffle.f90wrap_species_type__get__charge(self._handle)

        @property
        def radius(self):
            """
            The radius of the element.
            """
            return _raffle.f90wrap_species_type__get__radius(self._handle)

        @radius.setter
        def radius(self, radius):
            _raffle.f90wrap_species_type__set__radius(self._handle, radius)

        @charge.setter
        def charge(self, charge):
            _raffle.f90wrap_species_type__set__charge(self._handle, charge)

        @property
        def name(self):
            """
            The symbol of the element.
            """
            return _raffle.f90wrap_species_type__get__name(self._handle)

        @name.setter
        def name(self, name):
            _raffle.f90wrap_species_type__set__name(self._handle, name)

        @property
        def num(self):
            """
            The number of atoms of this species/element.
            """
            return _raffle.f90wrap_species_type__get__num(self._handle)

        @num.setter
        def num(self, num):
            _raffle.f90wrap_species_type__set__num(self._handle, num)

        def __str__(self):
            ret = ['<species_type>{\n']
            ret.append('    atom : ')
            ret.append(repr(self.atom))
            ret.append(',\n    mass : ')
            ret.append(repr(self.mass))
            ret.append(',\n    charge : ')
            ret.append(repr(self.charge))
            ret.append(',\n    name : ')
            ret.append(repr(self.name))
            ret.append(',\n    num : ')
            ret.append(repr(self.num))
            ret.append('}')
            return ''.join(ret)

        _dt_array_initialisers = []


    @f90wrap.runtime.register_class("raffle.basis")
    class basis(f90wrap.runtime.FortranDerivedType):
        def __init__(self, atoms=None, handle=None):
            """
            Create a ``basis`` object.

            This object is used to store the atomic information of a crystal,
            including lattice and basis information.
            This is confusingly named as a crystal = lattice + basis.

            Returns:
                basis (basis):
                    Object to be constructed
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _raffle.f90wrap_geom_rw__basis_type_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result

            if atoms is not None:
                self.fromase(atoms)

        def __del__(self):
            """
            Destructor for class basis


            Defined at ../src/lib/mod_geom_rw.f90 lines \
                34-42

            Parameters
            ----------
            this : basis
            	Object to be destructed


            Automatically generated destructor for basis
            """
            if self._alloc:
                _raffle.f90wrap_geom_rw__basis_type_finalise(this=self._handle)

        def allocate_species(self, num_species=None, species_symbols=None, species_count=None, \
            positions=None):
            """
            Allocate memory for the species list.

            Parameters:
                num_species (int):
                    Number of species
                species_symbols (list of str):
                    List of species symbols
                species_count (list of int):
                    List of species counts
                atoms (list of float):
                    List of atomic positions
            """
            _raffle.f90wrap_geom_rw__allocate_species__binding__basis_type(this=self._handle, \
                num_species=num_species, species_symbols=species_symbols, species_count=species_count, \
                atoms=positions)

        def _init_array_spec(self):
            """
            Initialise the species array.
            """
            self.spec = f90wrap.runtime.FortranDerivedTypeArray(self,
                                            _raffle.f90wrap_basis_type__array_getitem__spec,
                                            _raffle.f90wrap_basis_type__array_setitem__spec,
                                            _raffle.f90wrap_basis_type__array_len__spec,
                                            """
            Element spec ftype=type(species_type) pytype=species_type


            Defined at ../src/lib/mod_geom_rw.f90 line 35

            """, Geom_Rw.species_type)
            return self.spec

        def toase(self, calculator=None):
            """
            Convert the basis object to an ASE Atoms object.

            Parameters:
                calculator (ASE Calculator):
                    ASE calculator object to be assigned to the Atoms object.
            """
            from ase import Atoms

            # Set the species list
            positions = []
            species_string = ""
            for i in range(self.nspec):
                for j in range(self.spec[i].num):
                    species_string += str(self.spec[i].name.decode()).strip()
                    positions.append(self.spec[i].atom[j])

            # Set the atoms
            if(self.lcart):
                atoms = Atoms(species_string, positions=positions, cell=self.lat, pbc=self.pbc)
            else:
                atoms = Atoms(species_string, scaled_positions=positions, cell=self.lat, pbc=self.pbc)

            if calculator is not None:
                atoms.calc = calculator
            return atoms

        def fromase(self, atoms, verbose=False):
            """
            Convert the ASE Atoms object to a basis object.

            Parameters:
                atoms (ASE Atoms):
                    ASE Atoms object to be converted.
                verbose (bool):
                    Boolean whether to print warnings.
            """
            from ase.calculators.singlepoint import SinglePointCalculator

            # Get the species symbols
            species_symbols = atoms.get_chemical_symbols()
            species_symbols_unique = sorted(set(species_symbols))

            # Set the number of species
            self.nspec = len(species_symbols_unique)

            # Set the number of atoms
            self.natom = len(atoms)

            # check if calculator is present
            if atoms.calc is None:
                if verbose:
                    print("WARNING: No calculator present, setting energy to 0.0")
                atoms.calc = SinglePointCalculator(atoms, energy=0.0)
            self.energy = atoms.get_potential_energy()

            # # Set the lattice vectors
            self.lat = numpy.reshape(atoms.get_cell().flatten(), [3,3], order='A')
            self.pbc = atoms.pbc

            # Set the system name
            self.sysname = atoms.get_chemical_formula()

            # Set the species list
            species_count = []
            atom_positions = []
            positions = atoms.get_scaled_positions()
            for species in species_symbols_unique:
                species_count.append(sum([1 for symbol in species_symbols if symbol == species]))
                for j, symbol in enumerate(species_symbols):
                    if symbol == species:
                        atom_positions.append(positions[j])

            # Allocate memory for the atom list
            self.lcart = False
            self.allocate_species(species_symbols=species_symbols_unique, species_count=species_count, positions=atom_positions)

        @property
        def nspec(self):
            """
            The number of species in the basis.
            """
            return _raffle.f90wrap_basis_type__get__nspec(self._handle)

        @nspec.setter
        def nspec(self, nspec):
            _raffle.f90wrap_basis_type__set__nspec(self._handle, nspec)

        @property
        def natom(self):
            """
            The number of atoms in the basis.
            """
            return _raffle.f90wrap_basis_type__get__natom(self._handle)

        @natom.setter
        def natom(self, natom):
            _raffle.f90wrap_basis_type__set__natom(self._handle, natom)

        @property
        def energy(self):
            """
            The energy associated with the basis (or crystal).
            """
            return _raffle.f90wrap_basis_type__get__energy(self._handle)

        @energy.setter
        def energy(self, energy):
            _raffle.f90wrap_basis_type__set__energy(self._handle, energy)

        @property
        def lat(self):
            """
            The lattice vectors of the basis.
            """
            array_ndim, array_type, array_shape, array_handle = \
                _raffle.f90wrap_basis_type__array__lat(self._handle)
            if array_handle in self._arrays:
                lat = self._arrays[array_handle]
            else:
                lat = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _raffle.f90wrap_basis_type__array__lat)
                self._arrays[array_handle] = lat
            return lat

        @lat.setter
        def lat(self, lat):
            self.lat[...] = lat

        @property
        def lcart(self):
            """
            Boolean whether the atomic positions are in cartesian coordinates.
            """
            return _raffle.f90wrap_basis_type__get__lcart(self._handle)

        @lcart.setter
        def lcart(self, lcart):
            _raffle.f90wrap_basis_type__set__lcart(self._handle, lcart)

        @property
        def pbc(self):
            """
            Boolean array indicating the periodic boundary conditions.
            """
            array_ndim, array_type, array_shape, array_handle = \
                _raffle.f90wrap_basis_type__array__pbc(self._handle)
            if array_handle in self._arrays:
                pbc = self._arrays[array_handle]
            else:
                pbc = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _raffle.f90wrap_basis_type__array__pbc)
                self._arrays[array_handle] = pbc
            return pbc

        @pbc.setter
        def pbc(self, pbc):
            self.pbc[...] = pbc

        @property
        def sysname(self):
            """
            The name of the system.
            """
            return _raffle.f90wrap_basis_type__get__sysname(self._handle)

        @sysname.setter
        def sysname(self, sysname):
            _raffle.f90wrap_basis_type__set__sysname(self._handle, sysname)

        def __str__(self):
            ret = ['<basis>{\n']
            ret.append('    nspec : ')
            ret.append(repr(self.nspec))
            ret.append(',\n    natom : ')
            ret.append(repr(self.natom))
            ret.append(',\n    energy : ')
            ret.append(repr(self.energy))
            ret.append(',\n    lat : ')
            ret.append(repr(self.lat))
            ret.append(',\n    lcart : ')
            ret.append(repr(self.lcart))
            ret.append(',\n    pbc : ')
            ret.append(repr(self.pbc))
            ret.append(',\n    sysname : ')
            ret.append(repr(self.sysname))
            ret.append('}')
            return ''.join(ret)

        _dt_array_initialisers = [_init_array_spec]



    @f90wrap.runtime.register_class("raffle.basis_array")
    class basis_array(f90wrap.runtime.FortranDerivedType):
        def __init__(self, atoms=None, handle=None):
            """
            Create a ``basis_array`` object.


            Returns:
                basis_array (basis_array):
                    Object to be constructed
            """

            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _raffle.f90wrap_geom_rw__basis_type_xnum_array_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result


            # check if atoms is an ASE Atoms object or a list of ASE Atoms objects
            if atoms:
                from ase import Atoms
                if isinstance(atoms, Atoms):
                    self.allocate(1)
                    self.items[0].fromase(atoms)
                elif isinstance(atoms, list):
                    self.allocate(len(atoms))
                    for i, atom in enumerate(atoms):
                        self.items[i].fromase(atom)

        def __del__(self):
            """
            Destructor for class basis_array


            Defined at ../src/lib/mod_generator.f90 lines \
                19-21

            Parameters
            ----------
            this : basis_array
            	Object to be destructed


            Automatically generated destructor for basis_array
            """
            if self._alloc:
                _raffle.f90wrap_geom_rw__basis_type_xnum_array_finalise(this=self._handle)

        def _init_array_items(self):
            """
            Initialise the items array.
            """
            self.items = f90wrap.runtime.FortranDerivedTypeArray(self,
                                            _raffle.f90wrap_basis_type_xnum_array__array_getitem__items,
                                            _raffle.f90wrap_basis_type_xnum_array__array_setitem__items,
                                            _raffle.f90wrap_basis_type_xnum_array__array_len__items,
                                            """
            Element items ftype=type(basis_type) pytype=basis


            Defined at  line 0

            """, Geom_Rw.basis)
            return self.items

        def toase(self):
            """
            Convert the basis_array object to a list of ASE Atoms objects.
            """

            # Set the species list
            atoms = []
            for i in range(len(self.items)):
                atoms.append(self.items[i].toase())
            return atoms

        def allocate(self, size):
            """
            Allocate the items array with the given size.

            Parameters:
                size (int):
                    Size of the items array
            """
            _raffle.f90wrap_basis_type_xnum_array__array_alloc__items(self._handle, num=size)

        def deallocate(self):
            """
            Deallocate the items array
            """
            _raffle.f90wrap_basis_type_xnum_array__array_dealloc__items(self._handle)

        _dt_array_initialisers = [_init_array_items]

    _dt_array_initialisers = []


geom_rw = Geom_Rw()

class Raffle__Distribs_Container(f90wrap.runtime.FortranModule):
    """
    Code for handling distribution functions.

    This module provides the necessary functionality to create, update, and
    store distribution functions.
    The distribution functions are used as descriptors for the atomic
    environments in a crystal.
    The generalised distribution function (GDF) is a generalised descriptor
    for the atomic configurations that each species can adopt.

    Defined in ../src/fortran/lib/mod_distribs_container.f90

    """
    @f90wrap.runtime.register_class("raffle.distribs_container_type")
    class distribs_container_type(f90wrap.runtime.FortranDerivedType):
        def __init__(self, handle=None):
            """
            Create a ``Distribs_Container_Type`` object.

            Returns:
                distribution_container (Distribs_Container_Type):
                    Object to be constructed

            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = \
                _raffle.f90wrap_raffle__dc__dc_type_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result

        def __del__(self):
            """
            Destructor for class Distribs_Container_Type


            Defined at ../fortran/lib/mod_distribs_container.f90 \
                lines 25-162

            Parameters
            ----------
            this : Distribs_Container_Type
            	Object to be destructed


            Automatically generated destructor for distribs_container_type
            """
            if self._alloc:
                _raffle.f90wrap_raffle__dc__dc_type_finalise(this=self._handle)

        def set_kBT(self, kBT):
            """
            Set the energy scale for the distribution functions.

            Parameters:
                kBT (float):
                    Energy scale for the distribution functions.

            """
            self.kBT = kBT

        def set_weight_method(self, method):
            """
            Set the weight method for combining the the distribution functions
            to form the generalised distribution function (GDF).

            Parameters:
                method (str):
                    Method to be used for weighting the distribution functions.
                    Allowed values are:
                    - 'formation_energy' or 'formation' or 'form' or 'e_form'
                    - 'energy_above_hull' or 'hull_distance' or 'hull' or 'distance' or 'convex_hull'

            """

            if method in ['empirical', 'formation_energy', 'formation', 'form', 'e_form']:
                self.weight_by_hull = False
            elif method in ['energy_above_hull', 'hull_distance', 'hull', 'distance', 'convex_hull']:
                self.weight_by_hull = True
            else:
                raise ValueError("Invalid weight method: {}".format(method))

        def set_width(self, width):
            """
            Set the distribution function widths.

            Parameters:
                width (list of float):
                    List of distribution function widths.
                    The first element is the 2-body distribution function width,
                    the second element is the 3-body distribution function width,
                    and the third element is the 4-body distribution function width.
            """
            _raffle.f90wrap_raffle__dc__set_width__binding__dc_type(this=self._handle, \
                width=width)

        def set_sigma(self, sigma):
            """
            Set the sigma values of the Gaussians used to
            build the distribution functions.

            Parameters:
                sigma (list of float):
                    List of sigma values.
                    The first element is the 2-body distribution function sigma,
                    the second element is the 3-body distribution function sigma,
                    and the third element is the 4-body distribution function sigma.
            """
            _raffle.f90wrap_raffle__dc__set_sigma__binding__dc_type(this=self._handle, \
                sigma=sigma)

        def set_cutoff_min(self, cutoff_min):
            """
            Set the minimum cutoff values for the distribution functions.

            Parameters:
                cutoff_min (list of float):
                    List of minimum cutoff values.
                    The first element is the 2-body distribution function minimum cutoff,
                    the second element is the 3-body distribution function minimum cutoff,
                    and the third element is the 4-body distribution function minimum cutoff.
            """
            _raffle.f90wrap_raffle__dc__set_cutoff_min__binding__dc_type(this=self._handle, \
                cutoff_min=cutoff_min)

        def set_cutoff_max(self, cutoff_max):
            """
            Set the maximum cutoff values for the distribution functions.

            Parameters:
                cutoff_min (list of float):
                    List of maximum cutoff values.
                    The first element is the 2-body distribution function maximum cutoff,
                    the second element is the 3-body distribution function maximum cutoff,
                    and the third element is the 4-body distribution function maximum cutoff.
            """
            _raffle.f90wrap_raffle__dc__set_cutoff_max__binding__dc_type(this=self._handle, \
                cutoff_max=cutoff_max)

        def set_radius_distance_tol(self, radius_distance_tol):
            """
            Set the radius distance tolerance for the distribution functions.

            The radius distance tolerance represents a multiplier to the bond radii to
            determine the cutoff distance for the distribution functions.

            Parameters:
                radius_distance_tol (list of float):
                    List of radius distance tolerance values.
                    The first two values are the lower and upper bounds for the
                    3-body distribution function radius distance tolerance.
                    The third and fourth values are the lower and upper bounds for the
                    4-body distribution function radius distance tolerance.
            """
            _raffle.f90wrap_raffle__dc__set_radius_distance_tol__binding__dc_type(this=self._handle, \
                radius_distance_tol=radius_distance_tol)

        def create(self, basis_list, energy_above_hull_list=None, deallocate_systems=True, verbose=None):
            """
            Create the distribution functions.

            Parameters:
                basis_list (basis_array or Atoms or list of Atoms):
                    List of atomic configurations to be used to create the distribution functions.
                energy_above_hull_list (list of float):
                    List of energy above hull values for the atomic configurations.
                deallocate_systems (bool):
                    Boolean whether to deallocate the atomic configurations after creating the distribution functions.
                verbose (int):
                    Verbosity level.
            """
            from ase import Atoms
            if isinstance(basis_list, Atoms):
                basis_list = geom_rw.basis_array(basis_list)
            elif isinstance(basis_list, list):
                if all([isinstance(basis, Atoms) for basis in basis_list]):
                    basis_list = geom_rw.basis_array(basis_list)

            _raffle.f90wrap_raffle__dc__create__binding__dc_type(this=self._handle, \
                basis_list=basis_list._handle, \
                energy_above_hull_list=energy_above_hull_list, \
                deallocate_systems=deallocate_systems, \
                verbose=verbose \
            )

        def update(self, basis_list, energy_above_hull_list=None, from_host=True, deallocate_systems=True, verbose=None):
            """
            Update the distribution functions.

            Parameters:
                basis_list (basis_array or Atoms or list of Atoms):
                    List of atomic configurations to be used to create the distribution functions.
                energy_above_hull_list (list of float):
                    List of energy above hull values for the atomic configurations.
                deallocate_systems (bool):
                    Boolean whether to deallocate the atomic configurations after creating the distribution functions.
                from_host (bool):
                    Boolean whether the provided basis_list is based on the host.
                verbose (int):
                    Verbosity level.
            """
            from ase import Atoms
            if isinstance(basis_list, Atoms):
                basis_list = geom_rw.basis_array(basis_list)
            elif isinstance(basis_list, list):
                if all([isinstance(basis, Atoms) for basis in basis_list]):
                    basis_list = geom_rw.basis_array(basis_list)


            _raffle.f90wrap_raffle__dc__update__binding__dc_type(this=self._handle, \
                basis_list=basis_list._handle, \
                energy_above_hull_list=energy_above_hull_list, \
                from_host=from_host, \
                deallocate_systems=deallocate_systems, \
                verbose=verbose \
            )

        def deallocate_systems(self):
            """
            Deallocate the atomic configurations.
            """
            _raffle.f90wrap_raffle__dc__deallocate_systems__binding__dc_type(this=self._handle)

        def add_basis(self, basis):
            """
            Add a basis to the distribution functions.

            It is not recommended to use this function directly, but to use the
            create or update functions instead.

            Parameters:
                basis (geom_rw.basis):
                    Basis object to be added to the distribution functions.

            """
            _raffle.f90wrap_raffle__dc__add_basis__binding__dc_type(this=self._handle, \
                basis=basis._handle)

        def set_element_energies(self, element_energies):
            """
            Set the element reference energies for the distribution functions.

            These energies are used to calculate the formation energies of the
            atomic configurations.

            Parameters:
                element_energies (dict):
                    Dictionary of element reference energies.
                    The keys are the element symbols and the values are the reference energies.
            """

            element_list = list(element_energies.keys())
            energies = [element_energies[element] for element in element_list]
            _raffle.f90wrap_raffle__dc__set_element_energies__binding__dc_type(this=self._handle, \
                elements=element_list, energies=energies)

        def get_element_energies(self):
            """
            Get the element reference energies for the distribution functions.

            Returns:
                element_energies (dict):
                    Dictionary of element reference energies.
                    The keys are the element symbols and the values are the reference energies
            """

            num_elements = _raffle.f90wrap_raffle__dc__get__num_elements(self._handle)
            elements = numpy.zeros((num_elements,), dtype='S3')
            energies = numpy.zeros((num_elements,), dtype=numpy.float32)

            _raffle.f90wrap_raffle__dc__get_element_energies_sm__binding__dc_type(this=self._handle, \
                elements=elements, energies=energies)

            # convert the fortran array to a python dictionary
            element_energies = {}
            for i, element in enumerate(elements):
                name = str(element.decode()).strip()
                element_energies[name] = energies[i]

            return element_energies

        def set_bond_info(self):
            """
            Allocate the bond information array.

            It is not recommended to use this function directly, but to use the
            set_bond_radius or set_bond_radii functions instead.
            """
            _raffle.f90wrap_raffle__dc__set_bond_info__binding__dc_type(this=self._handle)

        def set_bond_radius(self, radius_dict):
            """
            Set the bond radius for the distribution functions.

            Parameters:
                radius_dict (dict):
                    Dictionary of bond radii.
                    The keys are a tuple of the two element symbols and the values are the bond radii.
            """

            # convert radius_dict to elements and radius
            # radius_dict = {('C', 'C'): 1.5}
            elements = list(radius_dict.keys()[0])
            radius = radius_dict.values()[0]

            _raffle.f90wrap_raffle__dc__set_bond_radius__binding__dc_type(this=self._handle, \
                elements=elements, radius=radius)

        def set_bond_radii(self, radius_dict):
            """
            Set the bond radii for the distribution functions.

            Parameters:
                radius_dict (dict):
                    Dictionary of bond radii.
                    The keys are a tuple of the two element symbols and the values are the bond radii.
            """

            # convert radius_list to elements and radii
            # radius_list = {('C', 'C'): 1.5, ('C', 'H'): 1.1}
            elements = []
            radii = []
            for key, value in radius_dict.items():
                elements.append(list(key))
                radii.append(value)


            _raffle.f90wrap_raffle__dc__set_bond_radii__binding__dc_type(this=self._handle, \
                elements=elements, radii=radii)

        def get_bond_radii(self):
            """
            Get the bond radii for the distribution functions.

            Returns:
                bond_radii (dict):
                    Dictionary of bond radii.
                    The keys are a tuple of the two element symbols and the values are the bond radii.
            """

            num_elements = _raffle.f90wrap_raffle__dc__get__num_elements(self._handle)
            if num_elements == 0:
                return {}
            num_pairs = round(num_elements * ( num_elements + 1 ) / 2)
            elements = numpy.zeros((num_pairs,2,), dtype='S3', order='F')
            radii = numpy.zeros((num_pairs,), dtype=numpy.float32, order='F')

            _raffle.f90wrap_raffle__dc__get_bond_radii_staticmem__binding__dc_type(this=self._handle, \
                elements=elements, radii=radii)
            # _raffle.f90wrap_raffle__dc__get_bond_radii_staticmem__binding__dc_type(this=self._handle, \
            #     elements=elements, energies=energies)

            # convert the fortran array to a python dictionary
            bond_radii = {}
            for i, element in enumerate(elements):
                names = tuple([str(name.decode()).strip() for name in element])
                bond_radii[names] = radii[i]

            return bond_radii

        def initialise_gdfs(self):
            """
            Initialise the generalised distribution functions (GDFs).

            It is not recommended to use this function directly, but to use the
            create or update functions instead.
            """
            _raffle.f90wrap_raffle__dc__initialise_gdfs__binding__dc_type(this=self._handle)

        def evolve(self): #, system=None):
            """
            Evolve the distribution functions.

            It is not recommended to use this function directly, but to use the
            create or update functions instead.
            """
            _raffle.f90wrap_raffle__dc__evolve__binding__dc_type(this=self._handle)
            # _raffle.f90wrap_raffle__dc__evolve__binding__dc_type(this=self._handle, \
            #     system=None if system is None else system._handle)

        def write_gdfs(self, file):
            """
            Write the generalised distribution functions (GDFs) to a file.

            Parameters:
                file (str):
                    Name of file to write the GDFs to.
            """
            _raffle.f90wrap_raffle__dc__write_gdfs__binding__dc_type(this=self._handle, \
                file=file)

        def read_gdfs(self, file):
            """
            Read the generalised distribution functions (GDFs) from a file.

            Parameters:
                file (str):
                    Name of file to read the GDFs from.
            """
            _raffle.f90wrap_raffle__dc__read_gdfs__binding__dc_type(this=self._handle, \
                file=file)

        def write_dfs(self, file):
            """
            Write the distribution functions (DFs) associated with all
            allocated systems to a file.

            Parameters:
                file (str):
                    Name of file to write the DFs to.
            """
            _raffle.f90wrap_raffle__dc__write_dfs__binding__dc_type(this=self._handle, \
                file=file)

        def read_dfs(self, file):
            """
            Read the distribution functions (DFs) associated with a set of
            systems from a file.

            Parameters:
                file (str):
                    Name of file to read the DFs from.
            """
            _raffle.f90wrap_raffle__dc__read_dfs__binding__dc_type(this=self._handle, \
                file=file)

        def write_2body(self, file):
            """
            Write the 2-body generalised distribution functions (GDFs) to a file.

            Parameters:
                file (str):
                    Name of file to write the 2-body GDFs to.
            """
            _raffle.f90wrap_raffle__dc__write_2body__binding__dc_type(this=self._handle, \
                file=file)

        def write_3body(self, file):
            """
            Write the 3-body generalised distribution functions (GDFs) to a file.

            Parameters:
                file (str):
                    Name of file to write the 3-body GDFs to.
            """
            _raffle.f90wrap_raffle__dc__write_3body__binding__dc_type(this=self._handle, \
                file=file)

        def write_4body(self, file):
            """
            Write the 4-body generalised distribution functions (GDFs) to a file.

            Parameters:
                file (str):
                    Name of file to write the 4-body GDFs to.
            """
            _raffle.f90wrap_raffle__dc__write_4body__binding__dc_type(this=self._handle, \
                file=file)

        def _get_pair_index(self, species1, species2):
            """
            Get the index of the pair of species in the distribution functions.

            This is meant as an internal function and not likely to be used directly.

            Parameters:
                species1 (str):
                    Name of the first species
                species2 (str):
                    Name of the second species

            Returns:
                idx (int):
                    Index of the pair of species in the distribution functions.
            """
            idx = \
                _raffle.f90wrap_raffle__dc__get_pair_index__binding__dc_type(this=self._handle, \
                species1=species1, species2=species2)
            return idx

        def _get_bin(self, value, dim):
            """
            Get the bin index for a value in the distribution functions.

            This is meant as an internal function and not likely to be used directly.

            Parameters:
                value (float):
                    Value to get the bin index for.
                dim (int):
                    Dimension of the distribution function.
                    1 for 2-body, 2 for 3-body, and 3 for 4-body.

            Returns
            -------
                bin (int):
                    Bin index for the value in the distribution functions.

            """
            bin = \
                _raffle.f90wrap_raffle__dc__get_bin__binding__dc_type(this=self._handle, \
                value=value, dim=dim)
            return bin

        @property
        def num_evaluated(self):
            """
            Number of evaluated distribution functions.
            """
            return _raffle.f90wrap_distribs_container_type__get__num_evaluated(self._handle)

        @num_evaluated.setter
        def num_evaluated(self, num_evaluated):
            _raffle.f90wrap_distribs_container_type__set__num_evaluated(self._handle, \
                num_evaluated)

        @property
        def num_evaluated_allocated(self):
            """
            Number of evaluated distribution functions still allocated.
            """
            return \
                _raffle.f90wrap_distribs_container_type__get__num_evaluated_allocated(self._handle)

        @num_evaluated_allocated.setter
        def num_evaluated_allocated(self, num_evaluated_allocated):
            _raffle.f90wrap_distribs_container_type__set__num_evaluated_allocated(self._handle, \
                num_evaluated_allocated)

        @property
        def kBT(self):
            """
            Energy scale for the distribution functions.
            """
            return _raffle.f90wrap_distribs_container_type__get__kbt(self._handle)

        @kBT.setter
        def kBT(self, kBT):
            _raffle.f90wrap_distribs_container_type__set__kbt(self._handle, kBT)

        @property
        def weight_by_hull(self):
            """
            Boolean whether to weight the distribution functions by the energy above hull.
            """
            return \
                _raffle.f90wrap_distribs_container_type__get__weight_by_hull(self._handle)

        @weight_by_hull.setter
        def weight_by_hull(self, weight_by_hull):
            _raffle.f90wrap_distribs_container_type__set__weight_by_hull(self._handle, \
                weight_by_hull)

        @property
        def nbins(self):
            """
            Number of bins in the distribution functions.
            """
            array_ndim, array_type, array_shape, array_handle = \
                _raffle.f90wrap_distribs_container_type__array__nbins(self._handle)
            if array_handle in self._arrays:
                nbins = self._arrays[array_handle]
            else:
                nbins = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _raffle.f90wrap_distribs_container_type__array__nbins)
                self._arrays[array_handle] = nbins
            return nbins

        @nbins.setter
        def nbins(self, nbins):
            self.nbins[...] = nbins

        @property
        def sigma(self):
            """
            Sigma values for the Gaussians used by the
            2-, 3-, and 4-body distribution functions.
            """
            array_ndim, array_type, array_shape, array_handle = \
                _raffle.f90wrap_distribs_container_type__array__sigma(self._handle)
            if array_handle in self._arrays:
                sigma = self._arrays[array_handle]
            else:
                sigma = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _raffle.f90wrap_distribs_container_type__array__sigma)
                self._arrays[array_handle] = sigma
            return sigma

        @sigma.setter
        def sigma(self, sigma):
            self.sigma[...] = sigma

        @property
        def width(self):
            """
            Bin widths for the 2-, 3-, and 4-body distribution functions.
            """
            array_ndim, array_type, array_shape, array_handle = \
                _raffle.f90wrap_distribs_container_type__array__width(self._handle)
            if array_handle in self._arrays:
                width = self._arrays[array_handle]
            else:
                width = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _raffle.f90wrap_distribs_container_type__array__width)
                self._arrays[array_handle] = width
            return width

        @width.setter
        def width(self, width):
            self.width[...] = width

        @property
        def cutoff_min(self):
            """
            The lower cutoff values for the 2-, 3-, and 4-body distribution functions.
            """
            array_ndim, array_type, array_shape, array_handle = \
                _raffle.f90wrap_distribs_container_type__array__cutoff_min(self._handle)
            if array_handle in self._arrays:
                cutoff_min = self._arrays[array_handle]
            else:
                cutoff_min = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _raffle.f90wrap_distribs_container_type__array__cutoff_min)
                self._arrays[array_handle] = cutoff_min
            return cutoff_min

        @cutoff_min.setter
        def cutoff_min(self, cutoff_min):
            self.cutoff_min[...] = cutoff_min

        @property
        def cutoff_max(self):
            """
            The upper cutoff values for the 2-, 3-, and 4-body distribution functions.
            """
            array_ndim, array_type, array_shape, array_handle = \
                _raffle.f90wrap_distribs_container_type__array__cutoff_max(self._handle)
            if array_handle in self._arrays:
                cutoff_max = self._arrays[array_handle]
            else:
                cutoff_max = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _raffle.f90wrap_distribs_container_type__array__cutoff_max)
                self._arrays[array_handle] = cutoff_max
            return cutoff_max

        @cutoff_max.setter
        def cutoff_max(self, cutoff_max):
            self.cutoff_max[...] = cutoff_max

        @property
        def radius_distance_tol(self):
            """
            The radius distance tolerance for the distribution functions.

            The radius distance tolerance represents a multiplier to the bond radii to
            determine the cutoff distance for the distribution functions.

            The first two values are the lower and upper bounds for the
            3-body distribution function radius distance tolerance.
            The third and fourth values are the lower and upper bounds for the
            4-body distribution function radius distance tolerance.
            """
            array_ndim, array_type, array_shape, array_handle = \
                _raffle.f90wrap_distribs_container_type__array__radius_distance_tol(self._handle)
            if array_handle in self._arrays:
                radius_distance_tol = self._arrays[array_handle]
            else:
                radius_distance_tol = \
                    f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _raffle.f90wrap_distribs_container_type__array__radius_distance_tol)
                self._arrays[array_handle] = radius_distance_tol
            return radius_distance_tol

        @radius_distance_tol.setter
        def radius_distance_tol(self, radius_distance_tol):
            self.radius_distance_tol[...] = radius_distance_tol

        def __str__(self):
            ret = ['<distribs_container_type>{\n']
            ret.append('    num_evaluated : ')
            ret.append(repr(self.num_evaluated))
            ret.append(',\n    num_evaluated_allocated : ')
            ret.append(repr(self.num_evaluated_allocated))
            ret.append(',\n    kBT : ')
            ret.append(repr(self.kBT))
            ret.append(',\n    weight_by_hull : ')
            ret.append(repr(self.weight_by_hull))
            ret.append(',\n    nbins : ')
            ret.append(repr(self.nbins))
            ret.append(',\n    sigma : ')
            ret.append(repr(self.sigma))
            ret.append(',\n    width : ')
            ret.append(repr(self.width))
            ret.append(',\n    cutoff_min : ')
            ret.append(repr(self.cutoff_min))
            ret.append(',\n    cutoff_max : ')
            ret.append(repr(self.cutoff_max))
            ret.append(',\n    radius_distance_tol : ')
            ret.append(repr(self.radius_distance_tol))
            ret.append('}')
            return ''.join(ret)

        _dt_array_initialisers = []#_init_array_system]


    _dt_array_initialisers = []


raffle__distribs_container = Raffle__Distribs_Container()

class Generator(f90wrap.runtime.FortranModule):
    """
    Code for generating interface structures.

    The module handles converting Python objects to Fortran objects and vice versa.
    These include converting between dictionaries and Fortran derived types, and
    between ASE Atoms objects and Fortran derived types.

    Defined in ../src/lib/mod_generator.f90

    """
    @f90wrap.runtime.register_class("raffle.stoichiometry_type")
    class stoichiometry_type(f90wrap.runtime.FortranDerivedType):
        def __init__(self, dict=None, element=None, num=None, handle=None):
            """
            Object to store the stoichiometry of a structure.

            Returns:
                stoichiometry (stoichiometry_type):
                    Stoichiometry object
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _raffle.f90wrap_stoichiometry_type_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result

            if element:
                self.element = element
            if num:
                self.num = num


        def __del__(self):
            """
            Destructor for class Stoichiometry_Type


            Defined at ../src/lib/mod_generator.f90 lines \
                19-21

            Automatically generated destructor for stoichiometry_type
            """
            if self._alloc:
                _raffle.f90wrap_stoichiometry_type_finalise(this=self._handle)

        @property
        def element(self):
            """
            String representing an element symbol.
            """
            return _raffle.f90wrap_stoichiometry_type__get__element(self._handle)

        @element.setter
        def element(self, element):
            _raffle.f90wrap_stoichiometry_type__set__element(self._handle, element)

        @property
        def num(self):
            """
            Integer representing the number of atoms of the element.
            """
            return _raffle.f90wrap_stoichiometry_type__get__num(self._handle)

        @num.setter
        def num(self, num):
            _raffle.f90wrap_stoichiometry_type__set__num(self._handle, num)

        def __str__(self):
            ret = ['<stoichiometry_type>{\n']
            ret.append('    element : ')
            ret.append(repr(self.element))
            ret.append(',\n    num : ')
            ret.append(repr(self.num))
            ret.append('}')
            return ''.join(ret)

        _dt_array_initialisers = []


    @f90wrap.runtime.register_class("raffle.stoichiometry_array")
    class stoichiometry_array(f90wrap.runtime.FortranDerivedType):
        def __init__(self, dict=None, handle=None):
            """
            Array or list of stoichiometry objects.

            Returns:
                stoichiometry_array (stoichiometry_array):
                    Stoichiometry array object
            """

            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _raffle.f90wrap_generator__stoich_type_xnum_array_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
            if dict:
                num_elements = len(dict)
                elements = list(dict.keys())
                nums = list(dict.values())
                self.allocate(num_elements)
                for i in range(num_elements):
                    self.items[i].element = elements[i]
                    self.items[i].num = nums[i]

        def __del__(self):
            """
            Destructor for class Stoichiometry_Type


            Defined at ../src/lib/mod_generator.f90 lines \
                19-21

            Parameters
            ----------
            this : Stoichiometry_Type
            	Object to be destructed


            Automatically generated destructor for stoichiometry_type
            """
            if self._alloc:
                _raffle.f90wrap_generator__stoich_type_xnum_array_finalise(this=self._handle)

        def _init_array_items(self):
            self.items = f90wrap.runtime.FortranDerivedTypeArray(self,
                                            _raffle.f90wrap_stoich_type_xnum_array__array_getitem__items,
                                            _raffle.f90wrap_stoich_type_xnum_array__array_setitem__items,
                                            _raffle.f90wrap_stoich_type_xnum_array__array_len__items,
                                            """
            Element items ftype=type(stoichiometry_type) pytype=stoichiometry_type


            Defined at  line 0

            """, Generator.stoichiometry_type)
            return self.items

        def allocate(self, size):
            """
            Allocate the items array with the given size.

            Parameters:
                size (int):
                    Size of the items array
            """
            _raffle.f90wrap_stoich_type_xnum_array__array_alloc__items(self._handle, num=size)

        def deallocate(self):
            """
            Deallocate the items array.
            """
            _raffle.f90wrap_stoich_type_xnum_array__array_dealloc__items(self._handle)



        _dt_array_initialisers = [_init_array_items]


    @f90wrap.runtime.register_class("raffle.raffle_generator")
    class raffle_generator(f90wrap.runtime.FortranDerivedType):

        def __init__(self, handle=None):
            """
            Create a ``raffle_generator`` object.

            This object is used to generate structures using the RAFFLE method.
            The object has procedures to set the parameters for the generation
            and to generate the structures.
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _raffle.f90wrap_generator__raffle_generator_type_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result

        def __del__(self):
            """
            Destructor for class raffle_generator


            Defined at ../src/lib/mod_generator.f90 lines \
                23-34

            """
            if self._alloc:
                _raffle.f90wrap_generator__raffle_generator_type_finalise(this=self._handle)

        def set_max_attempts(self, max_attempts):
            """
            Set the walk-method maximum attempts parameter.
            This parameter determines the maximum number of attempts to generate a structure
            using the walk method before reverting to the void method.

            Parameters:
                max_attempts (int):
                    The maximum number of attempts to generate a structure using the walk method.
            """
            self.max_attempts = max_attempts

        def set_walk_step_size(self, coarse=None, fine=None):
            """
            Set the walk-method step size parameters.

            Parameters:
                coarse (float):
                    The coarse step size for the walk method.
                fine (float):
                    The fine step size for the walk method.
            """
            if coarse is not None:
                self.walk_step_size_coarse = coarse
            if fine is not None:
                self.walk_step_size_fine = fine

        def set_host(self, host):
            """
            Set the host structure for the generation.

            Parameters:
                host (ase.Atoms or geom_rw.basis):
                    The host structure for the generation.
            """
            from ase import Atoms
            # check if host is ase.Atoms object or a Fortran derived type basis_type
            if isinstance(host, Atoms):
                host = geom_rw.basis(atoms=host)

            _raffle.f90wrap_generator__set_host__binding__rgt(this=self._handle, \
                host=host._handle)

        def set_grid(self, grid=None, grid_spacing=None, grid_offset=None):
            """
            Set the grid parameters for the generation.

            Parameters:
                grid (list of int):
                    The number of grid points along each axis of the host.
                grid_spacing (float)
                    The spacing between grid points.
                grid_offset (list of float):
                    The offset of the grid from the origin.
            """
            _raffle.f90wrap_generator__set_grid__binding__raffle_generator_type(this=self._handle, \
                grid=grid, grid_spacing=grid_spacing, grid_offset=grid_offset)

        def reset_grid(self):
            """
            Reset the grid parameters to their default values.
            """
            _raffle.f90wrap_generator__reset_grid__binding__raffle_generator_type(this=self._handle)

        def set_bounds(self, bounds=None):
            """
            Set the bounding box for the generation.

            Parameters:
                bounds (list of list of float):
                    The bounding box within which to constrain placement of atoms.
                    In the form [[a_min, a_max], [b_min, b_max], [c_min, c_max]].
                    Values given in direct (crystal) coordinates, ranging from 0 to 1.
            """
            _raffle.f90wrap_generator__set_bounds__binding__rgt(this=self._handle, \
                bounds=bounds)

        def reset_bounds(self):
            """
            Reset the bounding box to the full host structure.
            """
            _raffle.f90wrap_generator__reset_bounds__binding__rgt(this=self._handle)

        def generate(
                self, num_structures, stoichiometry,
                method_ratio={"void": 0.0, "rand": 0.0, "walk": 0.0, "grow": 0.0, "min": 0.0},
                method_probab=None,
                seed=None, settings_out_file=None, verbose=0,
                calc=None
        ):
            """
            Generate structures using the RAFFLE method.

            Parameters:
                num_structures (int):
                    The number of structures to generate.
                stoichiometry (stoichiometry_array or dict):
                    The stoichiometry of the structures to generate.
                method_ratio (dict):
                    The ratio of using each method to generate a structure.
                method_probab (dict):
                    DEPRECATED - The ratio of using each method to generate a structure.
                seed (int):
                    The seed for the random number generator.
                settings_out_file (str):
                    The file to write the settings to.
                verbose (int):
                    The verbosity level for the generation.
                calc (ASE calculator):
                    The calculator to use for the generated structures.
                
            Returns:
                structures (geom_rw.basis_array):
                    The generated structures.
                exit_code (int):
                    The exit code from the generation.
            """

            exit_code = 0
            structures = None

            # check if method_ratio is provided, if so, use it only if method_ratio is not provided
            if method_probab is not None and method_ratio == {"void": 0.0, "rand": 0.0, "walk": 0.0, "grow": 0.0, "min": 0.0}:
                method_ratio = method_probab
                warnings.warn("method_probab is deprecated, use method_ratio instead", DeprecationWarning)
            elif method_probab is not None:
                warnings.warn("method_probab is deprecated, use method_ratio instead", DeprecationWarning)
                # break if both method_ratio and method_probab are provided
                raise ValueError("Both method_ratio and method_probab are provided, use only one (method_ratio)")
            method_ratio_list = []
            method_ratio_list.append(method_ratio.get("void", 0.0))
            method_ratio_list.append(method_ratio.get("rand", 0.0)) # or method_ratio.get("random", 0.0))
            method_ratio_list.append(method_ratio.get("walk", 0.0))
            method_ratio_list.append(method_ratio.get("grow", 0.0)) # or method_ratio.get("growth", 0.0))
            method_ratio_list.append(method_ratio.get("min", 0.0))  # or method_ratio.get("minimum", 0.0) or method_ratio.get("global", 0.0))

            # check if all values are 0.0, if so, set them to the default of all 1.0
            if all([val < 1E-6 for val in method_ratio_list]):
                method_ratio_list = [1.0, 0.1, 0.5, 0.5, 1.0]

            # if stoichiometry is a dictionary, convert it to a stoichiometry_array
            if isinstance(stoichiometry, dict):
                stoichiometry = Generator.stoichiometry_array(dict=stoichiometry)

            exit_code = _raffle.f90wrap_generator__generate__binding__rgt(
                this=self._handle,
                num_structures=num_structures,
                stoichiometry=stoichiometry._handle,
                method_ratio=method_ratio_list,
                settings_out_file=settings_out_file,
                seed=seed, verbose=verbose)
                
            structures = self.get_structures(calc)[-num_structures:]
            return structures, exit_code

        def get_structures(self, calculator=None):
            """
            Get the generated structures as a list of ASE Atoms objects.

            Parameters:
                calculator (ASE calculator):
                    The calculator to use for the generated structures.
            """
            atoms = []
            for structure in self.structures:
                atoms.append(structure.toase(calculator))
            return atoms

        def set_structures(self, structures):
            """
            Set the list of generated structures.

            Parameters:
                structures (list of geom_rw.basis or list of ase.Atoms):
                    The list of structures to set.
            """
            structures = geom_rw.basis_array(atoms=structures)
            _raffle.f90wrap_generator__set_structures__binding__rgt(this=self._handle, \
                structures=structures._handle)

        def remove_structure(self, index):
            """
            Remove the structure at the given indices from the generator.

            Parameters:
                index (int or list of int):
                    The indices of the structure to remove.
            """
            index_list = [index] if isinstance(index, int) else index
            index_list = [ i + 1 for i in index_list ]
            _raffle.f90wrap_generator__remove_structure__binding__rgt(this=self._handle, \
                index_bn=index_list)

        def evaluate(self, basis):
            """
            Evaluate the viability of the structures.

            WARNING: This function is not implemented yet.

            Parameters:
                basis (geom_rw.basis or Atoms):
                    The basis to use for the evaluation.

            Returns:
                viability (float):
                    The viability of the structures.
            """
            from ase import Atoms
            if isinstance(basis, Atoms):
                basis = geom_rw.basis(atoms=basis)

            viability = \
                _raffle.f90wrap_generator__evaluate__binding__rgt(this=self._handle, \
                basis=basis._handle)
            return viability

        def print_settings(self, file):
            """
            Print the settings for the generation to a file.

            Parameters:
                file (str):
                    Name of the file to write the settings to.
            """
            _raffle.f90wrap_generator__print_settings__binding__rgt(this=self._handle, \
                file=file)

        def read_settings(self, file):
            """
            Read the settings for the generation from a file.

            Parameters:
                file (str):
                    Name of the file to read the settings from.
            """
            _raffle.f90wrap_generator__read_settings__binding__rgt(this=self._handle, \
                file=file)

        @property
        def num_structures(self):
            """
            The number of generated structures currently stored in the generator.
            """
            return _raffle.f90wrap_raffle_generator_type__get__num_structures(self._handle)

        @num_structures.setter
        def num_structures(self, num_structures):
            _raffle.f90wrap_raffle_generator_type__set__num_structures(self._handle, \
                num_structures)

        @property
        def host(self):
            """
            The host structure for the generation.
            """
            host_handle = _raffle.f90wrap_raffle_generator_type__get__host(self._handle)
            if tuple(host_handle) in self._objs:
                host = self._objs[tuple(host_handle)]
            else:
                host = geom_rw.basis.from_handle(host_handle)
                self._objs[tuple(host_handle)] = host
            return host

        @host.setter
        def host(self, host):
            host = host._handle
            _raffle.f90wrap_raffle_generator_type__set__host(self._handle, host)

        @property
        def grid(self):
            """
            The grid parameters for the generation.
            """
            array_ndim, array_type, array_shape, array_handle = \
                _raffle.f90wrap_raffle_generator_type__array__grid(self._handle)
            if array_handle in self._arrays:
                grid = self._arrays[array_handle]
            else:
                grid = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _raffle.f90wrap_raffle_generator_type__array__grid)
                self._arrays[array_handle] = grid
            return grid

        @grid.setter
        def grid(self, grid):
            self.grid[...] = grid

        @property
        def grid_offset(self):
            """
            The offset of the grid from the origin.
            """
            array_ndim, array_type, array_shape, array_handle = \
                _raffle.f90wrap_raffle_generator_type__array__grid_offset(self._handle)
            if array_handle in self._arrays:
                grid_offset = self._arrays[array_handle]
            else:
                grid_offset = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _raffle.f90wrap_raffle_generator_type__array__grid_offset)
                self._arrays[array_handle] = grid_offset
            return grid_offset

        @grid_offset.setter
        def grid_offset(self, grid_offset):
            self.grid_offset[...] = grid_offset

        @property
        def grid_spacing(self):
            """
            The spacing between grid points.
            """
            return _raffle.f90wrap_raffle_generator_type__get__grid_spacing(self._handle)

        @grid_spacing.setter
        def grid_spacing(self, grid_spacing):
            _raffle.f90wrap_raffle_generator_type__set__grid_spacing(self._handle, \
                grid_spacing)

        @property
        def bounds(self):
            """
            The bounds in direct coordinates of the host cell for the generation.
            """
            array_ndim, array_type, array_shape, array_handle = \
                _raffle.f90wrap_raffle_generator_type__array__bounds(self._handle)
            if array_handle in self._arrays:
                bounds = self._arrays[array_handle]
            else:
                bounds = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _raffle.f90wrap_raffle_generator_type__array__bounds)
                self._arrays[array_handle] = bounds
            return bounds

        @bounds.setter
        def bounds(self, bounds):
            self.bounds[...] = bounds

        @property
        def distributions(self):
            """
            The container for the distribution functions used in the generation.
            """
            distributions_handle = \
                _raffle.f90wrap_raffle_generator_type__get__distributions(self._handle)
            if tuple(distributions_handle) in self._objs:
                distributions = self._objs[tuple(distributions_handle)]
            else:
                distributions = \
                    raffle__distribs_container.distribs_container_type.from_handle(distributions_handle)
                self._objs[tuple(distributions_handle)] = distributions
            return distributions

        @distributions.setter
        def distributions(self, distributions):
            distributions = distributions._handle
            _raffle.f90wrap_raffle_generator_type__set__distributions(self._handle, \
                distributions)

        @property
        def max_attempts(self):
            """
            The maximum number of attempts to generate a structure using the walk method.
            """
            return _raffle.f90wrap_raffle_generator_type__get__max_attempts(self._handle)

        @max_attempts.setter
        def max_attempts(self, max_attempts):
            _raffle.f90wrap_raffle_generator_type__set__max_attempts(self._handle, \
                max_attempts)

        @property
        def walk_step_size_coarse(self):
            """
            The coarse step size for the walk method.
            """
            return \
                _raffle.f90wrap_raffle_generator_type__get__walk_step_size_coarse(self._handle)

        @walk_step_size_coarse.setter
        def walk_step_size_coarse(self, walk_step_size_coarse):
            _raffle.f90wrap_raffle_generator_type__set__walk_step_size_coarse(self._handle, \
                walk_step_size_coarse)

        @property
        def walk_step_size_fine(self):
            """
            The fine step size for the walk method.
            """
            return \
                _raffle.f90wrap_raffle_generator_type__get__walk_step_size_fine(self._handle)

        @walk_step_size_fine.setter
        def walk_step_size_fine(self, walk_step_size_fine):
            _raffle.f90wrap_raffle_generator_type__set__walk_step_size_fine(self._handle, \
                walk_step_size_fine)

        @property
        def method_ratio(self):
            """
            The ratio of methods to employ to generate a structure.
            """
            array_ndim, array_type, array_shape, array_handle = \
                _raffle.f90wrap_raffle_generator_type__array__method_ratio(self._handle)
            if array_handle in self._arrays:
                method_ratio = self._arrays[array_handle]
            else:
                method_ratio = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _raffle.f90wrap_raffle_generator_type__array__method_ratio)
                self._arrays[array_handle] = method_ratio
            return method_ratio

        @method_ratio.setter
        def method_ratio(self, method_ratio):
            self.method_ratio[...] = method_ratio

        def _init_array_structures(self):
            """
            Initialise the structures array.

            It is not recommended to use this function directly. Use the `structures` property instead.
            """
            self.structures = f90wrap.runtime.FortranDerivedTypeArray(self,
                                            _raffle.f90wrap_raffle_generator_type__array_getitem__structures,
                                            _raffle.f90wrap_raffle_generator_type__array_setitem__structures,
                                            _raffle.f90wrap_raffle_generator_type__array_len__structures,
                                            """
            Element items ftype=type(basis_type) pytype=basis


            Defined at ../src/lib/mod_generator.f90 line \
                29

            """, Geom_Rw.basis)
            return self.structures

        def __str__(self):
            ret = ['<raffle_generator>{\n']
            ret.append('    num_structures : ')
            ret.append(repr(self.num_structures))
            ret.append(',\n    host : ')
            ret.append(repr(self.host))
            ret.append(',\n    grid : ')
            ret.append(repr(self.grid))
            ret.append(',\n    grid_offset : ')
            ret.append(repr(self.grid_offset))
            ret.append(',\n    grid_spacing : ')
            ret.append(repr(self.grid_spacing))
            ret.append(',\n    bounds : ')
            ret.append(repr(self.bounds))
            ret.append(',\n    distributions : ')
            ret.append(repr(self.distributions))
            ret.append(',\n    max_attempts : ')
            ret.append(repr(self.max_attempts))
            ret.append(',\n    method_ratio : ')
            ret.append(repr(self.method_ratio))
            ret.append(',\n    structures : ')
            ret.append(repr(self.structures))
            ret.append('}')
            return ''.join(ret)

        _dt_array_initialisers = [_init_array_structures]


    _dt_array_initialisers = []


generator = Generator()

class Raffle(f90wrap.runtime.FortranModule):
    """
    Module raffle


    Defined at ../src/raffle.f90 lines 1-4

    """
    pass
    _dt_array_initialisers = []


raffle = Raffle()
