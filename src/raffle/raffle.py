from __future__ import print_function, absolute_import, division
import raffle._raffle as _raffle
import f90wrap.runtime
import logging
import numpy

class Rw_Geom(f90wrap.runtime.FortranModule):
    """
    Module rw_geom
    
    
    Defined at ../src/lib/mod_rw_geom.f90 lines \
        13-968
    
    """
    @f90wrap.runtime.register_class("raffle.species_type")
    class species_type(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=species_type)
        
        
        Defined at ../src/lib/mod_rw_geom.f90 lines \
            26-32
        
        """
        def __init__(self, handle=None):
            """
            self = species_type()
            
            
            Defined at ../src/lib/mod_rw_geom.f90 lines \
                26-32
            
            
            Returns
            -------
            this : species_type
            	Object to be constructed
            
            
            Automatically generated constructor for species_type
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _raffle.f90wrap_rw_geom__species_type_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class species_type
            
            
            Defined at ../src/lib/mod_rw_geom.f90 lines \
                26-32
            
            Parameters
            ----------
            this : species_type
            	Object to be destructed
            
            
            Automatically generated destructor for species_type
            """
            if self._alloc:
                _raffle.f90wrap_rw_geom__species_type_finalise(this=self._handle)
        
        @property
        def atom(self):
            """
            Element atom ftype=real(real32) pytype=float
            
            
            Defined at ../src/lib/mod_rw_geom.f90 line 27
            
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
            Element mass ftype=real(real32) pytype=float
            
            
            Defined at ../src/lib/mod_rw_geom.f90 line 28
            
            """
            return _raffle.f90wrap_species_type__get__mass(self._handle)
        
        @mass.setter
        def mass(self, mass):
            _raffle.f90wrap_species_type__set__mass(self._handle, mass)
        
        @property
        def charge(self):
            """
            Element charge ftype=real(real32) pytype=float
            
            
            Defined at ../src/lib/mod_rw_geom.f90 line 29
            
            """
            return _raffle.f90wrap_species_type__get__charge(self._handle)
        
        @property
        def radius(self):
            """
            Element radius ftype=real(real32) pytype=float
            
            
            Defined at ../src/lib/mod_rw_geom.f90 line 29
            
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
            Element name ftype=character(len=3) pytype=str
            
            
            Defined at ../src/lib/mod_rw_geom.f90 line 30
            
            """
            return _raffle.f90wrap_species_type__get__name(self._handle)
        
        @name.setter
        def name(self, name):
            _raffle.f90wrap_species_type__set__name(self._handle, name)
        
        @property
        def num(self):
            """
            Element num ftype=integer  pytype=int
            
            
            Defined at ../src/lib/mod_rw_geom.f90 line 31
            
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
        """
        Type(name=basis)
        
        
        Defined at ../src/lib/mod_rw_geom.f90 lines \
            34-42
        
        """
        def __init__(self, atoms=None, handle=None):
            """
            self = basis()
            
            
            Defined at ../src/lib/mod_rw_geom.f90 lines \
                34-42
            
            
            Returns
            -------
            this : basis
            	Object to be constructed
            
            
            Automatically generated constructor for basis
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _raffle.f90wrap_rw_geom__basis_type_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result

            if atoms is not None:
                self.fromase(atoms)
        
        def __del__(self):
            """
            Destructor for class basis
            
            
            Defined at ../src/lib/mod_rw_geom.f90 lines \
                34-42
            
            Parameters
            ----------
            this : basis
            	Object to be destructed
            
            
            Automatically generated destructor for basis
            """
            if self._alloc:
                _raffle.f90wrap_rw_geom__basis_type_finalise(this=self._handle)
        
        def allocate_species(self, num_species=None, species_symbols=None, species_count=None, \
            positions=None):
            """
            allocate_species__binding__basis_type(self[, num_species, species_symbols, \
                species_count, atoms])
            
            
            Defined at ../src/lib/mod_rw_geom.f90 lines \
                47-74
            
            Parameters
            ----------
            this : unknown
            num_species : int
            species_symbols : str array
            species_count : int array
            atoms : float array
            
            """
            _raffle.f90wrap_rw_geom__allocate_species__binding__basis_type(this=self._handle, \
                num_species=num_species, species_symbols=species_symbols, species_count=species_count, \
                atoms=positions)
        
        def _init_array_spec(self):
            self.spec = f90wrap.runtime.FortranDerivedTypeArray(self,
                                            _raffle.f90wrap_basis_type__array_getitem__spec,
                                            _raffle.f90wrap_basis_type__array_setitem__spec,
                                            _raffle.f90wrap_basis_type__array_len__spec,
                                            """
            Element spec ftype=type(species_type) pytype=species_type
            
            
            Defined at ../src/lib/mod_rw_geom.f90 line 35
            
            """, Rw_Geom.species_type)
            return self.spec
        
        def toase(self, calculator=None):
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
        
        def fromase(self, atoms):
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
            Element nspec ftype=integer  pytype=int
            
            
            Defined at ../src/lib/mod_rw_geom.f90 line 36
            
            """
            return _raffle.f90wrap_basis_type__get__nspec(self._handle)
        
        @nspec.setter
        def nspec(self, nspec):
            _raffle.f90wrap_basis_type__set__nspec(self._handle, nspec)
        
        @property
        def natom(self):
            """
            Element natom ftype=integer  pytype=int
            
            
            Defined at ../src/lib/mod_rw_geom.f90 line 37
            
            """
            return _raffle.f90wrap_basis_type__get__natom(self._handle)
        
        @natom.setter
        def natom(self, natom):
            _raffle.f90wrap_basis_type__set__natom(self._handle, natom)
        
        @property
        def energy(self):
            """
            Element energy ftype=real(real32) pytype=float
            
            
            Defined at ../src/lib/mod_rw_geom.f90 line 38
            
            """
            return _raffle.f90wrap_basis_type__get__energy(self._handle)
        
        @energy.setter
        def energy(self, energy):
            _raffle.f90wrap_basis_type__set__energy(self._handle, energy)
        
        @property
        def lat(self):
            """
            Element lat ftype=real(real32) pytype=float
            
            
            Defined at ../src/lib/mod_rw_geom.f90 line 38
            
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
            Element lcart ftype=logical pytype=bool
            
            
            Defined at ../src/lib/mod_rw_geom.f90 line 39
            
            """
            return _raffle.f90wrap_basis_type__get__lcart(self._handle)
        
        @lcart.setter
        def lcart(self, lcart):
            _raffle.f90wrap_basis_type__set__lcart(self._handle, lcart)
        
        @property
        def pbc(self):
            """
            Element pbc ftype=logical pytype=bool
            
            
            Defined at ../src/lib/mod_rw_geom.f90 line 40
            
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
            Element sysname ftype=character(len=1024) pytype=str
            
            
            Defined at ../src/lib/mod_rw_geom.f90 line 41
            
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
        """
        Type(name=basis_array)
        
        
        Defined at ../src/lib/mod_generator.f90 lines \
            19-21
        
        """
        def __init__(self, atoms=None, handle=None):
            """
            self = basis_array()
            
            
            Defined at ../src/lib/mod_generator.f90 lines \
                19-21
            
            
            Returns
            -------
            this : basis_array
            	Object to be constructed
            
            
            Automatically generated constructor for basis_array
            """

            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _raffle.f90wrap_rw_geom__basis_type_xnum_array_initialise()
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
                _raffle.f90wrap_rw_geom__basis_type_xnum_array_finalise(this=self._handle)
        
        def _init_array_items(self):
            self.items = f90wrap.runtime.FortranDerivedTypeArray(self,
                                            _raffle.f90wrap_basis_type_xnum_array__array_getitem__items,
                                            _raffle.f90wrap_basis_type_xnum_array__array_setitem__items,
                                            _raffle.f90wrap_basis_type_xnum_array__array_len__items,
                                            """
            Element items ftype=type(basis_type) pytype=basis
            
            
            Defined at  line 0
            
            """, Rw_Geom.basis)
            return self.items
        
        def toase(self):

            # Set the species list
            atoms = []
            for i in range(len(self.items)):
                atoms.append(self.items[i].toase())
            return atoms

        def allocate(self, size):
            """
            Allocate the items array with the given size
            
            Parameters
            ----------
            self : basis_conatiner
            size : int
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
    

rw_geom = Rw_Geom()

class Raffle__Distribs_Container(f90wrap.runtime.FortranModule):
    """
    Module raffle__distribs_container
    
    
    Defined at ../fortran/lib/mod_distribs_container.f90 \
        lines 1-1839
    
    """
    @f90wrap.runtime.register_class("raffle.distribs_container_type")
    class distribs_container_type(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=distribs_container_type)
        
        
        Defined at \
            ../fortran/lib/mod_distribs_container.f90 \
            lines 25-162
        
        """
        def __init__(self, handle=None):
            """
            self = Distribs_Container_Type()
            
            
            Defined at ../fortran/lib/mod_distribs_container.f90 \
                lines 25-162
            
            
            Returns
            -------
            this : Distribs_Container_Type
            	Object to be constructed
            
            
            Automatically generated constructor for distribs_container_type
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
            Parameters
            ----------
            this : unknown
            kBT : float
            
            """
            self.kBT = kBT

        def set_weight_method(self, method):
            """
            Parameters
            ----------
            this : unknown
            method : str
            
            """
            # method can be 'formation_energy' or 'energy_above_hull'
            # allowed abbreviations for 'formation_energy':
            #  'empirical', 'formation', 'form', 'e_form'
            # allowed abbreviations for 'hull_distance':
            #  'hull_distance', 'hull', 'distance', 'convex_hull'

            if method in ['empirical', 'formation_energy', 'formation', 'form', 'e_form']:
                self.weight_by_hull = False
            elif method in ['energy_above_hull', 'hull_distance', 'hull', 'distance', 'convex_hull']:
                self.weight_by_hull = True
            else:
                raise ValueError("Invalid weight method: {}".format(method))

        def set_width(self, width):
            """
            set_width__binding__dc_type(self, width)
            
            
            Defined at ../fortran/lib/mod_distribs_container.f90 \
                lines 237-247
            
            Parameters
            ----------
            this : unknown
            width : float array
            
            """
            _raffle.f90wrap_raffle__dc__set_width__binding__dc_type(this=self._handle, \
                width=width)
        
        def set_sigma(self, sigma):
            """
            set_sigma__binding__dc_type(self, sigma)
            
            
            Defined at ../fortran/lib/mod_distribs_container.f90 \
                lines 251-261
            
            Parameters
            ----------
            this : unknown
            sigma : float array
            
            """
            _raffle.f90wrap_raffle__dc__set_sigma__binding__dc_type(this=self._handle, \
                sigma=sigma)
        
        def set_cutoff_min(self, cutoff_min):
            """
            set_cutoff_min__binding__dc_type(self, cutoff_min)
            
            
            Defined at ../fortran/lib/mod_distribs_container.f90 \
                lines 265-273
            
            Parameters
            ----------
            this : unknown
            cutoff_min : float array
            
            """
            _raffle.f90wrap_raffle__dc__set_cutoff_min__binding__dc_type(this=self._handle, \
                cutoff_min=cutoff_min)
        
        def set_cutoff_max(self, cutoff_max):
            """
            set_cutoff_max__binding__dc_type(self, cutoff_max)
            
            
            Defined at ../fortran/lib/mod_distribs_container.f90 \
                lines 277-285
            
            Parameters
            ----------
            this : unknown
            cutoff_max : float array
            
            """
            _raffle.f90wrap_raffle__dc__set_cutoff_max__binding__dc_type(this=self._handle, \
                cutoff_max=cutoff_max)
        
        def set_radius_distance_tol(self, radius_distance_tol):
            """
            set_radius_distance_tol__binding__dc_type(self, \
                radius_distance_tol)
            
            
            Defined at ../fortran/lib/mod_distribs_container.f90 \
                lines 289-297
            
            Parameters
            ----------
            this : unknown
            radius_distance_tol : float array
            
            """
            _raffle.f90wrap_raffle__dc__set_radius_distance_tol__binding__dc_type(this=self._handle, \
                radius_distance_tol=radius_distance_tol)
        
        def create(self, basis_list, energy_above_hull_list=None, deallocate_systems=True):
            """
            create__binding__dc_type(self, basis_list)

            Defined at ../fortran/lib/mod_distribs_container.f90 \
                lines 353-440

            Parameters
            ----------
            this : unknown
            basis_list : basis_array or Atoms or list of Atoms
            energy_above_hull_list : list of float
            deallocate_systems : bool

            """
            from ase import Atoms
            if isinstance(basis_list, Atoms):
                basis_list = rw_geom.basis_array(basis_list)
            elif isinstance(basis_list, list):
                if all([isinstance(basis, Atoms) for basis in basis_list]):
                    basis_list = rw_geom.basis_array(basis_list)

            _raffle.f90wrap_raffle__dc__create__binding__dc_type(this=self._handle, \
                basis_list=basis_list._handle, \
                energy_above_hull_list=energy_above_hull_list, \
                deallocate_systems=deallocate_systems \
            )
            
        def update(self, basis_list, energy_above_hull_list=None, from_host=True, deallocate_systems=True):
            """
            update__binding__dc_type(self, basis_list)

            Defined at ../fortran/lib/mod_distribs_container.f90 \
                445-503

            Parameters
            ----------
            this : unknown
            basis_list : basis_array or list of Atoms
            energy_above_hull_list : list of float
            from_host : bool
            deallocate_systems : bool

            """
            from ase import Atoms
            if isinstance(basis_list, Atoms):
                basis_list = rw_geom.basis_array(basis_list)
            elif isinstance(basis_list, list):
                if all([isinstance(basis, Atoms) for basis in basis_list]):
                    basis_list = rw_geom.basis_array(basis_list)
            

            _raffle.f90wrap_raffle__dc__update__binding__dc_type(this=self._handle, \
                basis_list=basis_list._handle, \
                energy_above_hull_list=energy_above_hull_list, \
                from_host=from_host, \
                deallocate_systems=deallocate_systems \
            )
            
        def deallocate_systems(self):
            """
            deallocate_systems__binding__dc_type(self)
            
            
            Defined at ../fortran/lib/mod_distribs_container.f90 \
                lines 497-506
            
            Parameters
            ----------
            this : unknown
            
            """
            _raffle.f90wrap_raffle__dc__deallocate_systems__binding__dc_type(this=self._handle)
        
        def add_basis(self, basis):
            """
            add_basis__binding__dc_type(self, basis)
            
            
            Defined at ../fortran/lib/mod_distribs_container.f90 \
                lines 776-797
            
            Parameters
            ----------
            this : unknown
            basis : Basis_Type
            
            """
            _raffle.f90wrap_raffle__dc__add_basis__binding__dc_type(this=self._handle, \
                basis=basis._handle)
        
        def set_element_energies(self, element_energies):
            """
            set_element_energies__binding__dc_type(self, element_energies)
            
            
            Defined at ../fortran/lib/mod_distribs_container.f90 \
                lines 944-958
            
            Parameters
            ----------
            this : unknown
            element_energies : dict
            """

            element_list = list(element_energies.keys())
            energies = [element_energies[element] for element in element_list]
            _raffle.f90wrap_raffle__dc__set_element_energies__binding__dc_type(this=self._handle, \
                elements=element_list, energies=energies)
        
        def get_element_energies(self):
            """
            get_element_energies_static__binding__dc_type(self, elements, \
                energies)
            
            
            Defined at ../fortran/lib/mod_distribs_container.f90 \
                lines 984-1004
            
            Parameters
            ----------
            this : unknown

            Returns
            -------
            element_energies : dict
            
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
            set_bond_info__binding__dc_type(self)
            
            
            Defined at ../fortran/lib/mod_distribs_container.f90 \
                lines 1008-1052
            
            Parameters
            ----------
            this : unknown
            
            ---------------------------------------------------------------------------
             allocate the bond information array
            ---------------------------------------------------------------------------
            """
            _raffle.f90wrap_raffle__dc__set_bond_info__binding__dc_type(this=self._handle)
        
        def set_bond_radius(self, radius_dict):
            """
            set_bond_radius__binding__dc_type(self, elements, radius)
            
            
            Defined at ../fortran/lib/mod_distribs_container.f90 \
                lines 1197-1247
            
            Parameters
            ----------
            this : unknown
            radius_dict : dict
            
            ---------------------------------------------------------------------------
             remove python formatting
            ---------------------------------------------------------------------------
            """
            
            # convert radius_dict to elements and radius
            # radius_dict = {('C', 'C'): 1.5}
            elements = list(radius_dict.keys()[0])
            radius = radius_dict.values()[0]

            _raffle.f90wrap_raffle__dc__set_bond_radius__binding__dc_type(this=self._handle, \
                elements=elements, radius=radius)
        
        def set_bond_radii(self, radius_dict):
            """
            set_bond_radii__binding__dc_type(self, elements, radii)
            
            
            Defined at ../fortran/lib/mod_distribs_container.f90 \
                lines 1251-1266
            
            Parameters
            ----------
            this : unknown
            radius_dict : dict
            
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
            get_bond_radii_staticmem__binding__dc_type(self, elements, \
                radii)
            
            
            Defined at ../fortran/lib/mod_distribs_container.f90 \
                lines 1292-1312
            
            Parameters
            ----------
            this : unknown
            
            Returns
            -------
            bond_radii : dict
            
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
        
        def initialise_distribs(self):
            """
            initialise_distribs__binding__dc_type(self)
            
            
            Defined at ../fortran/lib/mod_distribs_container.f90 \
                lines 1474-1493
            
            Parameters
            ----------
            this : unknown
            
            """
            _raffle.f90wrap_raffle__dc__initialise_distribs__binding__dc_type(this=self._handle)
        
        def evolve(self): #, system=None):
            """
            evolve__binding__dc_type(self)
            
            
            Defined at ../fortran/lib/mod_distribs_container.f90 \
                lines 1539-1838
            
            Parameters
            ----------
            this : unknown
            
            ---------------------------------------------------------------------------
             if present, add the system to the container
            ---------------------------------------------------------------------------
            """
            _raffle.f90wrap_raffle__dc__evolve__binding__dc_type(this=self._handle)
            # _raffle.f90wrap_raffle__dc__evolve__binding__dc_type(this=self._handle, \
            #     system=None if system is None else system._handle)
        
        def write(self, file):
            """
            write__binding__dc_type(self, file)
            
            
            Defined at ../fortran/lib/mod_distribs_container.f90 \
                lines 510-559
            
            Parameters
            ----------
            this : unknown
            file : str
            
            """
            _raffle.f90wrap_raffle__dc__write__binding__dc_type(this=self._handle, \
                file=file)
        
        def read(self, file):
            """
            read__binding__dc_type(self, file)
            
            
            Defined at ../fortran/lib/mod_distribs_container.f90 \
                lines 563-620
            
            Parameters
            ----------
            this : unknown
            file : str
            
            """
            _raffle.f90wrap_raffle__dc__read__binding__dc_type(this=self._handle, \
                file=file)
        
        def write_2body(self, file):
            """
            write_2body__binding__dc_type(self, file)
            
            
            Defined at ../fortran/lib/mod_distribs_container.f90 \
                lines 624-663
            
            Parameters
            ----------
            this : unknown
            file : str
            
            """
            _raffle.f90wrap_raffle__dc__write_2body__binding__dc_type(this=self._handle, \
                file=file)
        
        def write_3body(self, file):
            """
            write_3body__binding__dc_type(self, file)
            
            
            Defined at ../fortran/lib/mod_distribs_container.f90 \
                lines 667-689
            
            Parameters
            ----------
            this : unknown
            file : str
            
            """
            _raffle.f90wrap_raffle__dc__write_3body__binding__dc_type(this=self._handle, \
                file=file)
        
        def write_4body(self, file):
            """
            write_4body__binding__dc_type(self, file)
            
            
            Defined at ../fortran/lib/mod_distribs_container.f90 \
                lines 693-715
            
            Parameters
            ----------
            this : unknown
            file : str
            
            """
            _raffle.f90wrap_raffle__dc__write_4body__binding__dc_type(this=self._handle, \
                file=file)
        
        def get_pair_index(self, species1, species2):
            """
            idx = get_pair_index__binding__dc_type(self, species1, species2)
            
            
            Defined at ../fortran/lib/mod_distribs_container.f90 \
                lines 1408-1430
            
            Parameters
            ----------
            this : unknown
            species1 : str
            species2 : str
            
            Returns
            -------
            idx : int
            
            """
            idx = \
                _raffle.f90wrap_raffle__dc__get_pair_index__binding__dc_type(this=self._handle, \
                species1=species1, species2=species2)
            return idx
        
        def get_bin(self, value, dim):
            """
            bin = get_bin__binding__dc_type(self, value, dim)
            
            
            Defined at ../fortran/lib/mod_distribs_container.f90 \
                lines 1451-1470
            
            Parameters
            ----------
            this : unknown
            value : float
            dim : int
            
            Returns
            -------
            bin : int
            
            """
            bin = \
                _raffle.f90wrap_raffle__dc__get_bin__binding__dc_type(this=self._handle, \
                value=value, dim=dim)
            return bin
        
        @property
        def num_evaluated(self):
            """
            Element num_evaluated ftype=integer  pytype=int
            
            
            Defined at ../fortran/lib/mod_distribs_container.f90 \
                line 30
            
            """
            return _raffle.f90wrap_distribs_container_type__get__num_evaluated(self._handle)
        
        @num_evaluated.setter
        def num_evaluated(self, num_evaluated):
            _raffle.f90wrap_distribs_container_type__set__num_evaluated(self._handle, \
                num_evaluated)
        
        @property
        def num_evaluated_allocated(self):
            """
            Element num_evaluated_allocated ftype=integer  pytype=int
            
            
            Defined at ../fortran/lib/mod_distribs_container.f90 \
                line 32
            
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
            Element kBT ftype=real(real32) pytype=float
            
            
            Defined at ../fortran/lib/mod_distribs_container.f90 \
                line 34
            
            """
            return _raffle.f90wrap_distribs_container_type__get__kbt(self._handle)
        
        @kBT.setter
        def kBT(self, kBT):
            _raffle.f90wrap_distribs_container_type__set__kbt(self._handle, kBT)
        
        @property
        def weight_by_hull(self):
            """
            Element weight_by_hull ftype=logical pytype=bool
            
            
            Defined at ../fortran/lib/mod_distribs_container.f90 \
                line 36
            
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
            Element nbins ftype=integer pytype=int
            
            
            Defined at ../fortran/lib/mod_distribs_container.f90 \
                line 54
            
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
            Element sigma ftype=real(real32) pytype=float
            
            
            Defined at ../fortran/lib/mod_distribs_container.f90 \
                line 57
            
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
            Element width ftype=real(real32) pytype=float
            
            
            Defined at ../fortran/lib/mod_distribs_container.f90 \
                line 61
            
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
            Element cutoff_min ftype=real(real32) pytype=float
            
            
            Defined at ../fortran/lib/mod_distribs_container.f90 \
                line 64
            
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
            Element cutoff_max ftype=real(real32) pytype=float
            
            
            Defined at ../fortran/lib/mod_distribs_container.f90 \
                line 67
            
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
            Element radius_distance_tol ftype=real(real32) pytype=float
            
            
            Defined at ../fortran/lib/mod_distribs_container.f90 \
                line 70
            
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
        
        # @property
        # def total(self):
        #     """
        #     Element total ftype=type(distribs_base_type) pytype=Distribs_Base_Type
            
            
        #     Defined at ../src/lib/mod_distribs_container.f90 line 38
            
        #     """
        #     total_handle = _raffle.f90wrap_distribs_container_type__get__total(self._handle)
        #     if tuple(total_handle) in self._objs:
        #         total = self._objs[tuple(total_handle)]
        #     else:
        #         total = distribs.distribs_base_type.from_handle(total_handle)
        #         self._objs[tuple(total_handle)] = total
        #     return total
        
        # @total.setter
        # def total(self, total):
        #     total = total._handle
        #     _raffle.f90wrap_distribs_container_type__set__total(self._handle, total)
        
        # def _init_array_system(self):
        #     self.system = f90wrap.runtime.FortranDerivedTypeArray(self,
        #                                     _raffle.f90wrap_distribs_container_type__array_getitem__system,
        #                                     _raffle.f90wrap_distribs_container_type__array_setitem__system,
        #                                     _raffle.f90wrap_distribs_container_type__array_len__system,
        #                                     """
        #     Element system ftype=type(distribs_type) pytype=Distribs_Type
            
            
        #     Defined at ../src/lib/mod_distribs_container.f90 line 39
            
        #     """, Distribs.distribs_type)
        #     return self.system
        
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
            # ret.append(',\n    total : ')
            # ret.append(repr(self.total))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []#_init_array_system]
        
    
    _dt_array_initialisers = []
    

raffle__distribs_container = Raffle__Distribs_Container()

class Generator(f90wrap.runtime.FortranModule):
    """
    Module generator
    
    
    Defined at ../src/lib/mod_generator.f90 lines \
        1-286
    
    """
    @f90wrap.runtime.register_class("raffle.stoichiometry_type")
    class stoichiometry_type(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=stoichiometry_type)
        
        
        Defined at ../src/lib/mod_generator.f90 lines \
            19-21
        
        """
        def __init__(self, dict=None, element=None, num=None, handle=None):
            """
            self = Stoichiometry_Type()
            
            
            Defined at ../src/lib/mod_generator.f90 lines \
                19-21
            
            
            Returns
            -------
            this : Stoichiometry_Type
            	Object to be constructed
            
            
            Automatically generated constructor for stoichiometry_type
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
            
            Parameters
            ----------
            this : Stoichiometry_Type
            	Object to be destructed
            
            
            Automatically generated destructor for stoichiometry_type
            """
            if self._alloc:
                _raffle.f90wrap_stoichiometry_type_finalise(this=self._handle)
        
        @property
        def element(self):
            """
            Element element ftype=character(len=3) pytype=str
            
            
            Defined at ../src/lib/mod_generator.f90 line \
                20
            
            """
            return _raffle.f90wrap_stoichiometry_type__get__element(self._handle)
        
        @element.setter
        def element(self, element):
            _raffle.f90wrap_stoichiometry_type__set__element(self._handle, element)
        
        @property
        def num(self):
            """
            Element num ftype=integer  pytype=int
            
            
            Defined at ../src/lib/mod_generator.f90 line \
                21
            
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
        """
        Type(name=stoichiometry_array)
        
        
        Defined at ../src/lib/mod_generator.f90 lines \
            19-21
        
        """
        def __init__(self, dict=None, handle=None):
            """
            self = Stoichiometry_Type()
            
            
            Defined at ../src/lib/mod_generator.f90 lines \
                19-21
            
            
            Returns
            -------
            this : Stoichiometry_Type
            	Object to be constructed
            
            
            Automatically generated constructor for stoichiometry_type
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
            Allocate the items array with the given size
            
            Parameters
            ----------
            self : Stoichiometry_Type
            size : int
                Size of the items array
            """
            _raffle.f90wrap_stoich_type_xnum_array__array_alloc__items(self._handle, num=size)

        def deallocate(self):
            """
            Deallocate the items array
            """
            _raffle.f90wrap_stoich_type_xnum_array__array_dealloc__items(self._handle)

        

        _dt_array_initialisers = [_init_array_items]
        
    
    @f90wrap.runtime.register_class("raffle.raffle_generator")
    class raffle_generator(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=raffle_generator)
        
        
        Defined at ../src/lib/mod_generator.f90 lines \
            23-34
        
        """
        def __init__(self, handle=None):
            """
            self = raffle_generator()
            
            
            Defined at ../src/lib/mod_generator.f90 lines \
                23-34
            
            
            Returns
            -------
            this : raffle_generator
            	Object to be constructed
            
            
            Automatically generated constructor for raffle_generator
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _raffle.f90wrap_generator__raffle_generator_type_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class raffle_generator
            
            
            Defined at ../src/lib/mod_generator.f90 lines \
                23-34
            
            Parameters
            ----------
            this : raffle_generator
            	Object to be destructed
            
            
            Automatically generated destructor for raffle_generator
            """
            if self._alloc:
                _raffle.f90wrap_generator__raffle_generator_type_finalise(this=self._handle)

        def set_max_attempts(self, max_attempts):
            """
            Parameters
            ----------
            this : unknown
            max_attempts : integer
            
            """
            self.max_attempts = max_attempts

        def set_walk_step_size(self, coarse=None, fine=None):
            """
            Parameters
            ----------
            this : unknown
            coarse : float
            fine: float
             
            """
            if coarse is not None:
                self.walk_step_size_coarse = coarse
            if fine is not None:
                self.walk_step_size_fine = fine
        
        def set_host(self, host):
            """
            set_host__binding__raffle_generator(self, host)
            
            
            Defined at ../src/lib/mod_generator.f90 lines \
                99-108
            
            Parameters
            ----------
            this : unknown
            host : basis
            
            """
            from ase import Atoms
            # check if host is ase.Atoms object
            if isinstance(host, Atoms):
                host = rw_geom.basis(atoms=host)

            _raffle.f90wrap_generator__set_host__binding__rgt(this=self._handle, \
                host=host._handle)
        
        def set_grid(self, grid=None, grid_spacing=None, grid_offset=None):
            """
            set_grid__binding__raffle_generator(self[, grid, grid_spacing, \
                grid_offset])
            
            
            Defined at ../fortran/lib/mod_generator.f90 lines \
                142-173
            
            Parameters
            ----------
            this : unknown
            grid : int array
            grid_spacing : float
            grid_offset : float array
            
            """
            _raffle.f90wrap_generator__set_grid__binding__raffle_generator_type(this=self._handle, \
                grid=grid, grid_spacing=grid_spacing, grid_offset=grid_offset)
        
        def reset_grid(self):
            """
            reset_grid__binding__raffle_generator(self)
            
            
            Defined at ../fortran/lib/mod_generator.f90 lines \
                170-176
            
            Parameters
            ----------
            this : unknown
            
            """
            _raffle.f90wrap_generator__reset_grid__binding__raffle_generator_type(this=self._handle)
        
        def generate(self, num_structures, stoichiometry, method_probab={"void": 0.0, "rand": 0.0, "walk": 0.0, "grow": 0.0, "min": 0.0}, seed=None, verbose=0):
            """
            generate__binding__raffle_generator(self, num_structures, stoichiometry, method_probab)

            Defined at ../src/lib/mod_generator.f90 lines \
                76-84

            Parameters
            ----------
            this : unknown
            num_structures : int
            stoichiometry : stoichiometry_array
            method_probab : dict
            verbose : int

            """
            
            method_probab_list = []
            method_probab_list.append(method_probab.get("void", 0.0))
            method_probab_list.append(method_probab.get("rand", 0.0)) # or method_probab.get("random", 0.0))
            method_probab_list.append(method_probab.get("walk", 0.0))
            method_probab_list.append(method_probab.get("grow", 0.0)) # or method_probab.get("growth", 0.0))
            method_probab_list.append(method_probab.get("min", 0.0))  # or method_probab.get("minimum", 0.0) or method_probab.get("global", 0.0))
            
            # check if all values are 0.0, if so, set them to the default of all 1.0
            if all([probab < 1E-6 for probab in method_probab_list]):
                method_probab_list = [1.0, 0.1, 0.5, 0.5, 1.0]

            # if stoichiometry is a dictionary, convert it to a stoichiometry_array
            if isinstance(stoichiometry, dict):
                stoichiometry = Generator.stoichiometry_array(dict=stoichiometry)

            if seed is not None:
                _raffle.f90wrap_generator__generate__binding__rgt(
                    this=self._handle,
                    num_structures=num_structures,
                    stoichiometry=stoichiometry._handle,
                    method_probab=method_probab_list, seed=seed, verbose=verbose)
            else:
                _raffle.f90wrap_generator__generate__binding__rgt(
                    this=self._handle,
                    num_structures=num_structures,
                    stoichiometry=stoichiometry._handle,
                    method_probab=method_probab_list, verbose=verbose)
            
        def get_structures(self, calculator=None):
            atoms = []
            for structure in self.structures:
                atoms.append(structure.toase(calculator))
            return atoms


        def evaluate(self, basis):
            """
            viability = evaluate__binding__raffle_generator(self, basis)
            
            
            Defined at ../src/lib/mod_generator.f90 lines \
                311-322
            
            Parameters
            ----------
            this : unknown
            basis : basis
            
            Returns
            -------
            viability : float
            
            """
            viability = \
                _raffle.f90wrap_generator__evaluate__binding__rgt(this=self._handle, \
                basis=basis._handle)
            return viability
        
        @property
        def num_structures(self):
            """
            Element num_structures ftype=integer  pytype=int
            
            
            Defined at ../src/lib/mod_generator.f90 line \
                24
            
            """
            return _raffle.f90wrap_raffle_generator_type__get__num_structures(self._handle)
        
        @num_structures.setter
        def num_structures(self, num_structures):
            _raffle.f90wrap_raffle_generator_type__set__num_structures(self._handle, \
                num_structures)
        
        @property
        def host(self):
            """
            Element host ftype=type(basis_type) pytype=basis
            
            
            Defined at ../src/lib/mod_generator.f90 line \
                25
            
            """
            host_handle = _raffle.f90wrap_raffle_generator_type__get__host(self._handle)
            if tuple(host_handle) in self._objs:
                host = self._objs[tuple(host_handle)]
            else:
                host = rw_geom.basis.from_handle(host_handle)
                self._objs[tuple(host_handle)] = host
            return host
        
        @host.setter
        def host(self, host):
            host = host._handle
            _raffle.f90wrap_raffle_generator_type__set__host(self._handle, host)
        
        @property
        def grid(self):
            """
            Element grid ftype=integer pytype=int
            
            
            Defined at ../src/lib/mod_generator.f90 line \
                24
            
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
            Element grid_offset ftype=real pytype=float
            
            
            Defined at ../fortran/lib/mod_generator.f90 line \
                45
            
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
            Element grid_spacing ftype=real(real32) pytype=float
            
            
            Defined at ../fortran/lib/mod_generator.f90 line \
                47
            
            """
            return _raffle.f90wrap_raffle_generator_type__get__grid_spacing(self._handle)
        
        @grid_spacing.setter
        def grid_spacing(self, grid_spacing):
            _raffle.f90wrap_raffle_generator_type__set__grid_spacing(self._handle, \
                grid_spacing)
        
        @property
        def distributions(self):
            """
            Element distributions ftype=type(distribs_container_type) \
                pytype=Distribs_Container_Type
            
            
            Defined at ../fortran/lib/mod_generator.f90 line \
                54
            
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
            Element max_attempts ftype=integer  pytype=int
            
            
            Defined at ../src/lib/mod_generator.f90 line \
                24
            
            """
            return _raffle.f90wrap_raffle_generator_type__get__max_attempts(self._handle)
        
        @max_attempts.setter
        def max_attempts(self, max_attempts):
            _raffle.f90wrap_raffle_generator_type__set__max_attempts(self._handle, \
                max_attempts)

        @property
        def walk_step_size_coarse(self):
            """
            Element walk_step_size_coarse ftype=real(real12) pytype=float
            
            
            Defined at \
                /Users/nedtaylor/DCoding/DGit/raffle/src/fortran/lib/mod_generator.f90 line \
                60
            
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
            Element walk_step_size_fine ftype=real(real12) pytype=float
            
            
            Defined at \
                /Users/nedtaylor/DCoding/DGit/raffle/src/fortran/lib/mod_generator.f90 line \
                60
            
            """
            return \
                _raffle.f90wrap_raffle_generator_type__get__walk_step_size_fine(self._handle)
        
        @walk_step_size_fine.setter
        def walk_step_size_fine(self, walk_step_size_fine):
            _raffle.f90wrap_raffle_generator_type__set__walk_step_size_fine(self._handle, \
                walk_step_size_fine)
        
        @property
        def method_probab(self):
            """
            Element method_probab ftype=real(real32) pytype=float
            
            
            Defined at ../src/lib/mod_generator.f90 line \
                67
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _raffle.f90wrap_raffle_generator_type__array__method_probab(self._handle)
            if array_handle in self._arrays:
                method_probab = self._arrays[array_handle]
            else:
                method_probab = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _raffle.f90wrap_raffle_generator_type__array__method_probab)
                self._arrays[array_handle] = method_probab
            return method_probab
        
        @method_probab.setter
        def method_probab(self, method_probab):
            self.method_probab[...] = method_probab
        
        def _init_array_structures(self):
            self.structures = f90wrap.runtime.FortranDerivedTypeArray(self,
                                            _raffle.f90wrap_raffle_generator_type__array_getitem__structures,
                                            _raffle.f90wrap_raffle_generator_type__array_setitem__structures,
                                            _raffle.f90wrap_raffle_generator_type__array_len__structures,
                                            """
            Element items ftype=type(basis_type) pytype=basis
            
            
            Defined at ../src/lib/mod_generator.f90 line \
                29
            
            """, Rw_Geom.basis)
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
            ret.append(',\n    distributions : ')
            ret.append(repr(self.distributions))
            ret.append(',\n    max_attempts : ')
            ret.append(repr(self.max_attempts))
            ret.append(',\n    method_probab : ')
            ret.append(repr(self.method_probab))
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

