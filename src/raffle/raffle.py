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
            Element atom ftype=real(real12) pytype=float
            
            
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
            Element mass ftype=real(real12) pytype=float
            
            
            Defined at ../src/lib/mod_rw_geom.f90 line 28
            
            """
            return _raffle.f90wrap_species_type__get__mass(self._handle)
        
        @mass.setter
        def mass(self, mass):
            _raffle.f90wrap_species_type__set__mass(self._handle, mass)
        
        @property
        def charge(self):
            """
            Element charge ftype=real(real12) pytype=float
            
            
            Defined at ../src/lib/mod_rw_geom.f90 line 29
            
            """
            return _raffle.f90wrap_species_type__get__charge(self._handle)
        
        @property
        def radius(self):
            """
            Element radius ftype=real(real12) pytype=float
            
            
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
        
    
    @f90wrap.runtime.register_class("raffle.basis_type")
    class basis_type(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=basis_type)
        
        
        Defined at ../src/lib/mod_rw_geom.f90 lines \
            34-42
        
        """
        def __init__(self, atoms=None, handle=None):
            """
            self = basis_type()
            
            
            Defined at ../src/lib/mod_rw_geom.f90 lines \
                34-42
            
            
            Returns
            -------
            this : basis_type
            	Object to be constructed
            
            
            Automatically generated constructor for basis_type
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _raffle.f90wrap_rw_geom__basis_type_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result

            if atoms is not None:
                self.fromase(atoms)
        
        def __del__(self):
            """
            Destructor for class basis_type
            
            
            Defined at ../src/lib/mod_rw_geom.f90 lines \
                34-42
            
            Parameters
            ----------
            this : basis_type
            	Object to be destructed
            
            
            Automatically generated destructor for basis_type
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
        
        def init_array_spec(self):
            self.spec = f90wrap.runtime.FortranDerivedTypeArray(self,
                                            _raffle.f90wrap_basis_type__array_getitem__spec,
                                            _raffle.f90wrap_basis_type__array_setitem__spec,
                                            _raffle.f90wrap_basis_type__array_len__spec,
                                            """
            Element spec ftype=type(species_type) pytype=species_type
            
            
            Defined at ../src/lib/mod_rw_geom.f90 line 35
            
            """, Rw_Geom.species_type)
            return self.spec
        
        def toase(self):
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
            Element energy ftype=real(real12) pytype=float
            
            
            Defined at ../src/lib/mod_rw_geom.f90 line 38
            
            """
            return _raffle.f90wrap_basis_type__get__energy(self._handle)
        
        @energy.setter
        def energy(self, energy):
            _raffle.f90wrap_basis_type__set__energy(self._handle, energy)
        
        @property
        def lat(self):
            """
            Element lat ftype=real(real12) pytype=float
            
            
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
            ret = ['<basis_type>{\n']
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
        
        _dt_array_initialisers = [init_array_spec]
        


    @f90wrap.runtime.register_class("raffle.basis_type_xnum_array")
    class basis_type_xnum_array(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=basis_type_xnum_array)
        
        
        Defined at ../src/lib/mod_generator.f90 lines \
            19-21
        
        """
        def __init__(self, handle=None):
            """
            self = basis_Type()
            
            
            Defined at ../src/lib/mod_generator.f90 lines \
                19-21
            
            
            Returns
            -------
            this : basis_Type
            	Object to be constructed
            
            
            Automatically generated constructor for basis_type
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _raffle.f90wrap_rw_geom__basis_type_xnum_array_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class basis_Type
            
            
            Defined at ../src/lib/mod_generator.f90 lines \
                19-21
            
            Parameters
            ----------
            this : basis_type
            	Object to be destructed
            
            
            Automatically generated destructor for basis_type
            """
            if self._alloc:
                _raffle.f90wrap_rw_geom__basis_type_xnum_array_finalise(this=self._handle)
        
        def init_array_items(self):
            self.items = f90wrap.runtime.FortranDerivedTypeArray(self,
                                            _raffle.f90wrap_basis_type_xnum_array__array_getitem__items,
                                            _raffle.f90wrap_basis_type_xnum_array__array_setitem__items,
                                            _raffle.f90wrap_basis_type_xnum_array__array_len__items,
                                            """
            Element items ftype=type(basis_type) pytype=basis_type
            
            
            Defined at  line 0
            
            """, Rw_Geom.basis_type)
            return self.items
        
        def allocate(self, size):
            """
            Allocate the items array with the given size
            
            Parameters
            ----------
            self : basis_type
            size : int
                Size of the items array
            """
            _raffle.f90wrap_basis_type_xnum_array__array_alloc__items(self._handle, num=size)

        def deallocate(self):
            """
            Deallocate the items array
            """
            _raffle.f90wrap_basis_type_xnum_array__array_dealloc__items(self._handle)

        
        _dt_array_initialisers = [init_array_items]


    # @staticmethod
    # def geom_read(unit, lat, bas, length=None):
    #     """
    #     geom_read(unit, lat, bas[, length])
        
        
    #     Defined at ../src/lib/mod_rw_geom.f90 lines \
    #         79-111
        
    #     Parameters
    #     ----------
    #     unit : int
    #     lat : float array
    #     bas : basis_type
    #     length : int
        
    #     """
    #     _raffle.f90wrap_rw_geom__geom_read(unit=unit, lat=lat, bas=bas._handle, \
    #         length=length)
    
    # @staticmethod
    # def geom_write(unit, lat, bas):
    #     """
    #     geom_write(unit, lat, bas)
        
        
    #     Defined at ../src/lib/mod_rw_geom.f90 lines \
    #         117-139
        
    #     Parameters
    #     ----------
    #     unit : int
    #     lat : float array
    #     bas : basis_type
        
    #     """
    #     _raffle.f90wrap_rw_geom__geom_write(unit=unit, lat=lat, bas=bas._handle)
    
    # @staticmethod
    # def convert_bas(self, latconv):
    #     """
    #     outbas = convert_bas(self, latconv)
        
        
    #     Defined at ../src/lib/mod_rw_geom.f90 lines \
    #         821-840
        
    #     Parameters
    #     ----------
    #     inbas : basis_type
    #     latconv : float array
        
    #     Returns
    #     -------
    #     outbas : basis_type
        
    #     """
    #     outbas = _raffle.f90wrap_rw_geom__convert_bas(inbas=self._handle, \
    #         latconv=latconv)
    #     outbas = f90wrap.runtime.lookup_class("raffle.basis_type").from_handle(outbas, \
    #         alloc=True)
    #     return outbas
    
    # @staticmethod
    # def clone_bas(self, outbas, inlat=None, outlat=None, trans_dim=None):
    #     """
    #     clone_bas(self, outbas[, inlat, outlat, trans_dim])
        
        
    #     Defined at ../src/lib/mod_rw_geom.f90 lines \
    #         897-967
        
    #     Parameters
    #     ----------
    #     inbas : basis_type
    #     outbas : basis_type
    #     inlat : float array
    #     outlat : float array
    #     trans_dim : bool
        
    #     -----------------------------------------------------------------------------
    #      determines whether user wants output basis extra translational dimension
    #     -----------------------------------------------------------------------------
    #     """
    #     _raffle.f90wrap_rw_geom__clone_bas(inbas=self._handle, outbas=outbas._handle, \
    #         inlat=inlat, outlat=outlat, trans_dim=trans_dim)
    
    # @property
    # def igeom_input(self):
    #     """
    #     Element igeom_input ftype=integer  pytype=int
        
        
    #     Defined at ../src/lib/mod_rw_geom.f90 line 24
        
    #     """
    #     return _raffle.f90wrap_rw_geom__get__igeom_input()
    
    # @igeom_input.setter
    # def igeom_input(self, igeom_input):
    #     _raffle.f90wrap_rw_geom__set__igeom_input(igeom_input)
    
    # @property
    # def igeom_output(self):
    #     """
    #     Element igeom_output ftype=integer  pytype=int
        
        
    #     Defined at ../src/lib/mod_rw_geom.f90 line 24
        
    #     """
    #     return _raffle.f90wrap_rw_geom__get__igeom_output()
    
    # @igeom_output.setter
    # def igeom_output(self, igeom_output):
    #     _raffle.f90wrap_rw_geom__set__igeom_output(igeom_output)
    
    # def __str__(self):
    #     ret = ['<rw_geom>{\n']
    #     ret.append('    igeom_input : ')
    #     ret.append(repr(self.igeom_input))
    #     ret.append(',\n    igeom_output : ')
    #     ret.append(repr(self.igeom_output))
    #     ret.append('}')
    #     return ''.join(ret)
    
    _dt_array_initialisers = []
    

rw_geom = Rw_Geom()

class Evolver(f90wrap.runtime.FortranModule):
    """
    Module evolver
    
    
    Defined at ../src/lib/mod_evolver.f90 lines \
        1-1204
    
    """
    @f90wrap.runtime.register_class("raffle.gvector_base_type")
    class gvector_base_type(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=gvector_base_type)
        
        
        Defined at ../src/lib/mod_evolver.f90 lines \
            14-17
        
        """
        def __init__(self, handle=None):
            """
            self = Gvector_Base_Type()
            
            
            Defined at ../src/lib/mod_evolver.f90 lines \
                14-17
            
            
            Returns
            -------
            this : Gvector_Base_Type
            	Object to be constructed
            
            
            Automatically generated constructor for gvector_base_type
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _raffle.f90wrap_evolver__gvector_base_type_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Gvector_Base_Type
            
            
            Defined at ../src/lib/mod_evolver.f90 lines \
                14-17
            
            Parameters
            ----------
            this : Gvector_Base_Type
            	Object to be destructed
            
            
            Automatically generated destructor for gvector_base_type
            """
            if self._alloc:
                _raffle.f90wrap_evolver__gvector_base_type_finalise(this=self._handle)
        
        @property
        def df_2body(self):
            """
            Element df_2body ftype=real(real12) pytype=float
            
            
            Defined at ../src/lib/mod_evolver.f90 line 15
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _raffle.f90wrap_gvector_base_type__array__df_2body(self._handle)
            if array_handle in self._arrays:
                df_2body = self._arrays[array_handle]
            else:
                df_2body = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _raffle.f90wrap_gvector_base_type__array__df_2body)
                self._arrays[array_handle] = df_2body
            return df_2body
        
        @df_2body.setter
        def df_2body(self, df_2body):
            self.df_2body[...] = df_2body
        
        @property
        def df_3body(self):
            """
            Element df_3body ftype=real(real12) pytype=float
            
            
            Defined at ../src/lib/mod_evolver.f90 line 16
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _raffle.f90wrap_gvector_base_type__array__df_3body(self._handle)
            if array_handle in self._arrays:
                df_3body = self._arrays[array_handle]
            else:
                df_3body = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _raffle.f90wrap_gvector_base_type__array__df_3body)
                self._arrays[array_handle] = df_3body
            return df_3body
        
        @df_3body.setter
        def df_3body(self, df_3body):
            self.df_3body[...] = df_3body
        
        @property
        def df_4body(self):
            """
            Element df_4body ftype=real(real12) pytype=float
            
            
            Defined at ../src/lib/mod_evolver.f90 line 17
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _raffle.f90wrap_gvector_base_type__array__df_4body(self._handle)
            if array_handle in self._arrays:
                df_4body = self._arrays[array_handle]
            else:
                df_4body = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _raffle.f90wrap_gvector_base_type__array__df_4body)
                self._arrays[array_handle] = df_4body
            return df_4body
        
        @df_4body.setter
        def df_4body(self, df_4body):
            self.df_4body[...] = df_4body
        
        def __str__(self):
            ret = ['<gvector_base_type>{\n']
            ret.append('    df_2body : ')
            ret.append(repr(self.df_2body))
            ret.append(',\n    df_3body : ')
            ret.append(repr(self.df_3body))
            ret.append(',\n    df_4body : ')
            ret.append(repr(self.df_4body))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    @f90wrap.runtime.register_class("raffle.gvector_type")
    class gvector_type(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=gvector_type)
        
        
        Defined at ../src/lib/mod_evolver.f90 lines \
            19-25
        
        """
        def __init__(self, handle=None):
            """
            self = Gvector_Type()
            
            
            Defined at ../src/lib/mod_evolver.f90 lines \
                19-25
            
            
            Returns
            -------
            this : Gvector_Type
            	Object to be constructed
            
            
            Automatically generated constructor for gvector_type
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _raffle.f90wrap_evolver__gvector_type_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Gvector_Type
            
            
            Defined at ../src/lib/mod_evolver.f90 lines \
                19-25
            
            Parameters
            ----------
            this : Gvector_Type
            	Object to be destructed
            
            
            Automatically generated destructor for gvector_type
            """
            if self._alloc:
                _raffle.f90wrap_evolver__gvector_type_finalise(this=self._handle)
        
        def calculate(self, lattice, basis, nbins=None, width=None, sigma=None, \
            cutoff_min=None, cutoff_max=None, radius_distance_tol=None):
            """
            calculate__binding__gvector_type(self, lattice, basis[, nbins, width, sigma, \
                cutoff_min, cutoff_max])
            
            
            Defined at ../src/lib/mod_evolver.f90 lines \
                747-1135
            
            Parameters
            ----------
            this : unknown
            lattice : float array
            basis : Basis_Type
            nbins : int array
            width : float array
            sigma : float array
            cutoff_min : float array
            cutoff_max : float array
            radius_distance_tol : float array
            
            --------------------------------------------------------------------------
             initialise optional variables
            --------------------------------------------------------------------------
            """
            _raffle.f90wrap_evolver__calculate__binding__gvector_type(this=self._handle, \
                lattice=lattice, basis=basis._handle, nbins=nbins, width=width, sigma=sigma, \
                cutoff_min=cutoff_min, cutoff_max=cutoff_max, radius_distance_tol=radius_distance_tol)
        
        @property
        def num_atoms(self):
            """
            Element num_atoms ftype=integer  pytype=int
            
            
            Defined at ../src/lib/mod_evolver.f90 line 20
            
            """
            return _raffle.f90wrap_gvector_type__get__num_atoms(self._handle)
        
        @num_atoms.setter
        def num_atoms(self, num_atoms):
            _raffle.f90wrap_gvector_type__set__num_atoms(self._handle, num_atoms)
        
        @property
        def energy(self):
            """
            Element energy ftype=real(real12) pytype=float
            
            
            Defined at ../src/lib/mod_evolver.f90 line 21
            
            """
            return _raffle.f90wrap_gvector_type__get__energy(self._handle)
        
        @energy.setter
        def energy(self, energy):
            _raffle.f90wrap_gvector_type__set__energy(self._handle, energy)
        
        @property
        def stoichiometry(self):
            """
            Element stoichiometry ftype=integer pytype=int
            
            
            Defined at ../src/lib/mod_evolver.f90 line 22
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _raffle.f90wrap_gvector_type__array__stoichiometry(self._handle)
            if array_handle in self._arrays:
                stoichiometry = self._arrays[array_handle]
            else:
                stoichiometry = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _raffle.f90wrap_gvector_type__array__stoichiometry)
                self._arrays[array_handle] = stoichiometry
            return stoichiometry
        
        @stoichiometry.setter
        def stoichiometry(self, stoichiometry):
            self.stoichiometry[...] = stoichiometry
        
        @property
        def species(self):
            """
            Element species ftype=character(len=3) pytype=str
            
            
            Defined at ../src/lib/mod_evolver.f90 line 23
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _raffle.f90wrap_gvector_type__array__species(self._handle)
            if array_handle in self._arrays:
                species = self._arrays[array_handle]
            else:
                species = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _raffle.f90wrap_gvector_type__array__species)
                self._arrays[array_handle] = species
            return species
        
        @species.setter
        def species(self, species):
            self.species[...] = species
        
        def __str__(self):
            ret = ['<gvector_type>{\n']
            ret.append('    num_atoms : ')
            ret.append(repr(self.num_atoms))
            ret.append(',\n    energy : ')
            ret.append(repr(self.energy))
            ret.append(',\n    stoichiometry : ')
            ret.append(repr(self.stoichiometry))
            ret.append(',\n    species : ')
            ret.append(repr(self.species))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    @f90wrap.runtime.register_class("raffle.gvector_container_type")
    class gvector_container_type(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=gvector_container_type)
        
        
        Defined at ../src/lib/mod_evolver.f90 lines \
            30-62
        
        """
        def __init__(self, handle=None):
            """
            self = Gvector_Container_Type()
            
            
            Defined at ../src/lib/mod_evolver.f90 lines \
                30-62
            
            
            Returns
            -------
            this : Gvector_Container_Type
            	Object to be constructed
            
            
            Automatically generated constructor for gvector_container_type
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _raffle.f90wrap_evolver__gvector_container_type_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Gvector_Container_Type
            
            
            Defined at ../src/lib/mod_evolver.f90 lines \
                30-62
            
            Parameters
            ----------
            this : Gvector_Container_Type
            	Object to be destructed
            
            
            Automatically generated destructor for gvector_container_type
            """
            if self._alloc:
                _raffle.f90wrap_evolver__gvector_container_type_finalise(this=self._handle)
        
        def set_width(self, width):
            """
            set_width__binding__gvector_container_type(self, width)
            
            
            Defined at ../src/lib/mod_evolver.f90 lines \
                108-118
            
            Parameters
            ----------
            this : unknown
            width : float array
            
            """
            _raffle.f90wrap_evolver__set_width__binding__gvector_container_type(this=self._handle, \
                width=width)
        
        def set_sigma(self, sigma):
            """
            set_sigma__binding__gvector_container_type(self, sigma)
            
            
            Defined at ../src/lib/mod_evolver.f90 lines \
                120-130
            
            Parameters
            ----------
            this : unknown
            sigma : float array
            
            """
            _raffle.f90wrap_evolver__set_sigma__binding__gvector_container_type(this=self._handle, \
                sigma=sigma)
        
        def set_cutoff_min(self, cutoff_min):
            """
            set_cutoff_min__binding__gvector_container_type(self, cutoff_min)
            
            
            Defined at ../src/lib/mod_evolver.f90 lines \
                132-140
            
            Parameters
            ----------
            this : unknown
            cutoff_min : float array
            
            """
            _raffle.f90wrap_evolver__set_cutoff_min__binding__gvector_container7007(this=self._handle, \
                cutoff_min=cutoff_min)
        
        def set_cutoff_max(self, cutoff_max):
            """
            set_cutoff_max__binding__gvector_container_type(self, cutoff_max)
            
            
            Defined at ../src/lib/mod_evolver.f90 lines \
                142-150
            
            Parameters
            ----------
            this : unknown
            cutoff_max : float array
            
            """
            _raffle.f90wrap_evolver__set_cutoff_max__binding__gvector_container047c(this=self._handle, \
                cutoff_max=cutoff_max)
        
        def set_radius_distance_tol(self, radius_distance_tol):
            """
            set_radius_distance_tol__binding__gvector_container_type(self, radius_distance_tol)
            
            
            Defined at ../src/lib/mod_evolver.f90 lines \
                142-150
            
            Parameters
            ----------
            this : unknown
            cutoff_max : float array
            
            """
            _raffle.f90wrap_evolver__set_radius_distance_tol__binding__gvector_1dda(this=self._handle, \
                radius_distance_tol=radius_distance_tol)
        
        def create(self, basis_list, deallocate_systems=True):
            """
            create__binding__gvector_container_type(self, basis_list)

            Defined at ../src/lib/mod_evolver.f90 lines \
                152-162

            Parameters
            ----------
            this : unknown
            basis_list : Basis_Type array
            deallocate_systems : bool

            """
            _raffle.f90wrap_evolver__create__binding__gvector_container_type(this=self._handle, \
                basis_list=basis_list._handle, deallocate_systems=deallocate_systems)
            
        def update(self, basis_list, deallocate_systems=True):
            """
            update__binding__gvector_container_type(self, basis_list)

            Defined at ../src/lib/mod_evolver.f90 lines \
                152-162

            Parameters
            ----------
            this : unknown
            basis_list : Basis_Type array
            deallocate_systems : bool

            """
            _raffle.f90wrap_evolver__update__binding__gvector_container_type(this=self._handle, \
                basis_list=basis_list._handle, deallocate_systems=deallocate_systems)
            
        def deallocate_systems(self):
            """
            deallocate_systems__binding__gvector_container_type(self)
            
            
            Defined at /Users/nedtaylor/DCoding/DGit/raffle/src/fortran/lib/mod_evolver.f90 \
                lines 323-331
            
            Parameters
            ----------
            this : unknown
            
            """
            _raffle.f90wrap_evolver__deallocate_systems__binding__gvector_conta8f02(this=self._handle)
        
        def add_basis(self, lattice, basis):
            """
            add_basis__binding__gvector_container_type(self, lattice, basis)
            
            
            Defined at ../src/lib/mod_evolver.f90 lines \
                415-430
            
            Parameters
            ----------
            this : unknown
            lattice : float array
            basis : Basis_Type
            
            """
            _raffle.f90wrap_evolver__add_basis__binding__gvector_container_type(this=self._handle, \
                lattice=lattice, basis=basis._handle)
        
        # def set_element_info(self, element_file=None, element_list=None):
        #     """
        #     set_element_info__binding__gvector_container_type(self[, element_file, \
        #         element_list])
            
            
        #     Defined at ../src/lib/mod_evolver.f90 lines \
        #         436-466
            
        #     Parameters
        #     ----------
        #     this : unknown
        #     element_file : str
        #     element_list : str array
            
        #     --------------------------------------------------------------------------
        #      load the elements database
        #     --------------------------------------------------------------------------
        #     """
        #     _raffle.f90wrap_evolver__set_element_info__binding__gvector_containbcb0(this=self._handle, \
        #         element_file=element_file, element_list=element_list)

        def set_element_energies(self, element_energies):
            """
            set_element_energies__binding__gvector_container_type(self, element_energies)

            Defined at ../src/lib/mod_evolver.f90 lines \
                472-526
            
            Parameters
            ----------
            this : unknown
            element_energies : dict
            """

            element_list = list(element_energies.keys())
            energies = [element_energies[element] for element in element_list]
            _raffle.f90wrap_evolver__set_element_energies__binding__gvector_con0537(this=self._handle, \
                elements=element_list, energies=energies)

        def get_element_energies(self):
            """
            get_element_energies_static__binding__gvector_container_type(self, elements, \
                energies)
            
            
            Defined at /Users/nedtaylor/DCoding/DGit/raffle/src/fortran/lib/mod_evolver.f90 \
                lines 557-574
            
            Parameters
            ----------
            this : unknown

            Returns
            -------
            element_energies : dict
            
            """

            num_elements = _raffle.f90wrap_evolver__get__num_elements(self._handle)
            elements = numpy.zeros((num_elements,), dtype='S3')
            energies = numpy.zeros((num_elements,), dtype=numpy.float32)

            _raffle.f90wrap_evolver__get_element_energies_staticmem__binding__g4f53(this=self._handle, \
                elements=elements, energies=energies)
            
            # convert the fortran array to a python dictionary
            element_energies = {}
            for i, element in enumerate(elements):
                name = str(element.decode()).strip()
                element_energies[name] = energies[i]

            return element_energies

        def set_bond_info(self, bond_file=None):
            """
            set_bond_info__binding__gvector_container_type(self[, bond_file])
            
            
            Defined at ../src/lib/mod_evolver.f90 lines \
                472-526
            
            Parameters
            ----------
            this : unknown
            bond_file : str
            
            --------------------------------------------------------------------------
             load the element bonds database
            --------------------------------------------------------------------------
            """
            _raffle.f90wrap_evolver__set_bond_info__binding__gvector_container_type(this=self._handle, \
                bond_file=bond_file)
        
        def set_bond_radius(self, radius_dict):
            """
            set_bond_radius__binding__gvector_container_type(self, elements, radius)
            
            
            Defined at /Users/nedtaylor/DCoding/DGit/raffle/src/fortran/lib/mod_evolver.f90 \
                lines 711-757
            
            Parameters
            ----------
            this : unknown
            elements : str array
            radius : float
            
            ---------------------------------------------------------------------------
             remove python formatting
            ---------------------------------------------------------------------------
            """
            
            # convert radius_dict to elements and radius
            # radius_dict = {('C', 'C'): 1.5}
            elements = list(radius_dict.keys()[0])
            radius = radius_dict.values()[0]

            _raffle.f90wrap_evolver__set_bond_radius__binding__gvector_containe7df9(this=self._handle, \
                elements=elements, radius=radius)
        
        def set_bond_radii(self, radius_dict):
            """
            set_bond_radii__binding__gvector_container_type(self, elements, radii)
            
            
            Defined at /Users/nedtaylor/DCoding/DGit/raffle/src/fortran/lib/mod_evolver.f90 \
                lines 761-776
            
            Parameters
            ----------
            this : unknown
            elements : str array
            radii : float array
            
            """

            # convert radius_list to elements and radii
            # radius_list = {('C', 'C'): 1.5, ('C', 'H'): 1.1}
            elements = []
            radii = []
            for key, value in radius_dict.items():
                elements.append(list(key))
                radii.append(value)
               

            _raffle.f90wrap_evolver__set_bond_radii__binding__gvector_container83c5(this=self._handle, \
                elements=elements, radii=radii)
        
        def get_bond_radii(self):
            """
            get_bond_radii_staticmem__binding__gvector_container_type(self, elements, radii)
            
            
            Defined at /Users/nedtaylor/DCoding/DGit/raffle/src/fortran/lib/mod_evolver.f90 \
                lines 808-828
            
            Parameters
            ----------
            this : unknown
            elements : str array
            radii : float array
            
            Returns
            -------
            element_energies : dict
            
            """

            num_elements = _raffle.f90wrap_evolver__get__num_elements(self._handle)
            if num_elements == 0:
                return {}
            num_pairs = round(num_elements * ( num_elements + 1 ) / 2)
            elements = numpy.zeros((num_pairs,2,), dtype='S3', order='F')
            radii = numpy.zeros((num_pairs,), dtype=numpy.float32, order='F')

            _raffle.f90wrap_evolver__get_bond_radii_staticmem__binding__gvectord2e1(this=self._handle, \
                elements=elements, radii=radii)
            # _raffle.f90wrap_evolver__get_element_energies_staticmem__binding__g4f53(this=self._handle, \
            #     elements=elements, energies=energies)
            
            # convert the fortran array to a python dictionary
            bond_radii = {}
            for i, element in enumerate(elements):
                names = tuple([str(name.decode()).strip() for name in element])
                bond_radii[names] = radii[i]

            return bond_radii
        
        def set_best_energy(self):
            """
            set_best_energy__binding__gvector_container_type(self)
            
            
            Defined at ../src/lib/mod_evolver.f90 lines \
                532-554
            
            Parameters
            ----------
            this : unknown
            
            """
            _raffle.f90wrap_evolver__set_best_energy__binding__gvector_containe4680(this=self._handle)
        
        def initialise_gvectors(self):
            """
            initialise_gvectors__binding__gvector_container_type(self)
            
            
            Defined at ../src/lib/mod_evolver.f90 lines \
                600-630
            
            Parameters
            ----------
            this : unknown
            
            """
            _raffle.f90wrap_evolver__initialise_gvectors__binding__gvector_contc1f2(this=self._handle)
        
        def evolve(self, system=None):
            """
            evolve__binding__gvector_container_type(self[, system, \
                deallocate_systems_after_evolve])
            
            
            Defined at ../src/lib/mod_evolver.f90 lines \
                637-740
            
            Parameters
            ----------
            this : unknown
            system : Gvector_Type array
            deallocate_systems_after_evolve : bool
            
            --------------------------------------------------------------------------
             if present, set the deallocate flag
            --------------------------------------------------------------------------
            """
            _raffle.f90wrap_evolver__evolve__binding__gvector_container_type(this=self._handle, \
                system=None if system is None else system._handle)
        
        def write(self, file):
            """
            write__binding__gvector_container_type(self, file)
            
            
            Defined at ../src/lib/mod_evolver.f90 lines \
                182-210
            
            Parameters
            ----------
            this : unknown
            file : str
            
            """
            _raffle.f90wrap_evolver__write__binding__gvector_container_type(this=self._handle, \
                file=file)
        
        def read(self, file):
            """
            read__binding__gvector_container_type(self, file)
            
            
            Defined at ../src/lib/mod_evolver.f90 lines \
                216-260
            
            Parameters
            ----------
            this : unknown
            file : str
            
            """
            _raffle.f90wrap_evolver__read__binding__gvector_container_type(this=self._handle, \
                file=file)
        
        def write_2body(self, file):
            """
            write_2body__binding__gvector_container_type(self, file)
            
            
            Defined at ../src/lib/mod_evolver.f90 lines \
                266-295
            
            Parameters
            ----------
            this : unknown
            file : str
            
            """
            _raffle.f90wrap_evolver__write_2body__binding__gvector_container_type(this=self._handle, \
                file=file)
        
        def write_3body(self, file):
            """
            write_3body__binding__gvector_container_type(self, file)
            
            
            Defined at ../src/lib/mod_evolver.f90 lines \
                301-315
            
            Parameters
            ----------
            this : unknown
            file : str
            
            """
            _raffle.f90wrap_evolver__write_3body__binding__gvector_container_type(this=self._handle, \
                file=file)
        
        def write_4body(self, file):
            """
            write_4body__binding__gvector_container_type(self, file)
            
            
            Defined at ../src/lib/mod_evolver.f90 lines \
                321-335
            
            Parameters
            ----------
            this : unknown
            file : str
            
            """
            _raffle.f90wrap_evolver__write_4body__binding__gvector_container_type(this=self._handle, \
                file=file)
        
        def get_pair_index(self, species1, species2):
            """
            idx = get_pair_index__binding__gvector_container_type(self, species1, species2)
            
            
            Defined at ../src/lib/mod_evolver.f90 lines \
                560-575
            
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
                _raffle.f90wrap_evolver__get_pair_index__binding__gvector_container4618(this=self._handle, \
                species1=species1, species2=species2)
            return idx
        
        def get_bin(self, value, dim):
            """
            bin = get_bin__binding__gvector_container_type(self, value, dim)
            
            
            Defined at ../src/lib/mod_evolver.f90 lines \
                581-594
            
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
                _raffle.f90wrap_evolver__get_bin__binding__gvector_container_type(this=self._handle, \
                value=value, dim=dim)
            return bin
        
        @property
        def num_evaluated(self):
            """
            Element num_evaluated ftype=integer  pytype=int
            
            
            Defined at /Users/nedtaylor/DCoding/DGit/raffle/src/fortran/lib/mod_evolver.f90 \
                line 57
            
            """
            return _raffle.f90wrap_gvector_container_type__get__num_evaluated(self._handle)
        
        @num_evaluated.setter
        def num_evaluated(self, num_evaluated):
            _raffle.f90wrap_gvector_container_type__set__num_evaluated(self._handle, \
                num_evaluated)
        
        @property
        def num_evaluated_allocated(self):
            """
            Element num_evaluated_allocated ftype=integer  pytype=int
            
            
            Defined at /Users/nedtaylor/DCoding/DGit/raffle/src/fortran/lib/mod_evolver.f90 \
                line 59
            
            """
            return \
                _raffle.f90wrap_gvector_container_type__get__num_evaluated_allocated(self._handle)
        
        @num_evaluated_allocated.setter
        def num_evaluated_allocated(self, num_evaluated_allocated):
            _raffle.f90wrap_gvector_container_type__set__num_evaluated_allocated(self._handle, \
                num_evaluated_allocated)
        
        @property
        def best_system(self):
            """
            Element best_system ftype=integer  pytype=int
            
            
            Defined at ../src/lib/mod_evolver.f90 line 31
            
            """
            return _raffle.f90wrap_gvector_container_type__get__best_system(self._handle)
        
        @best_system.setter
        def best_system(self, best_system):
            _raffle.f90wrap_gvector_container_type__set__best_system(self._handle, \
                best_system)
        
        @property
        def best_energy(self):
            """
            Element best_energy ftype=real(real12) pytype=float
            
            
            Defined at ../src/lib/mod_evolver.f90 line 32
            
            """
            return _raffle.f90wrap_gvector_container_type__get__best_energy(self._handle)
        
        @best_energy.setter
        def best_energy(self, best_energy):
            _raffle.f90wrap_gvector_container_type__set__best_energy(self._handle, \
                best_energy)
        
        @property
        def nbins(self):
            """
            Element nbins ftype=integer pytype=int
            
            
            Defined at ../src/lib/mod_evolver.f90 line 33
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _raffle.f90wrap_gvector_container_type__array__nbins(self._handle)
            if array_handle in self._arrays:
                nbins = self._arrays[array_handle]
            else:
                nbins = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _raffle.f90wrap_gvector_container_type__array__nbins)
                self._arrays[array_handle] = nbins
            return nbins
        
        @nbins.setter
        def nbins(self, nbins):
            self.nbins[...] = nbins
        
        @property
        def sigma(self):
            """
            Element sigma ftype=real(real12) pytype=float
            
            
            Defined at ../src/lib/mod_evolver.f90 line 34
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _raffle.f90wrap_gvector_container_type__array__sigma(self._handle)
            if array_handle in self._arrays:
                sigma = self._arrays[array_handle]
            else:
                sigma = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _raffle.f90wrap_gvector_container_type__array__sigma)
                self._arrays[array_handle] = sigma
            return sigma
        
        @sigma.setter
        def sigma(self, sigma):
            self.sigma[...] = sigma
        
        @property
        def width(self):
            """
            Element width ftype=real(real12) pytype=float
            
            
            Defined at ../src/lib/mod_evolver.f90 line 35
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _raffle.f90wrap_gvector_container_type__array__width(self._handle)
            if array_handle in self._arrays:
                width = self._arrays[array_handle]
            else:
                width = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _raffle.f90wrap_gvector_container_type__array__width)
                self._arrays[array_handle] = width
            return width
        
        @width.setter
        def width(self, width):
            self.width[...] = width
        
        @property
        def cutoff_min(self):
            """
            Element cutoff_min ftype=real(real12) pytype=float
            
            
            Defined at ../src/lib/mod_evolver.f90 line 36
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _raffle.f90wrap_gvector_container_type__array__cutoff_min(self._handle)
            if array_handle in self._arrays:
                cutoff_min = self._arrays[array_handle]
            else:
                cutoff_min = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _raffle.f90wrap_gvector_container_type__array__cutoff_min)
                self._arrays[array_handle] = cutoff_min
            return cutoff_min
        
        @cutoff_min.setter
        def cutoff_min(self, cutoff_min):
            self.cutoff_min[...] = cutoff_min
        
        @property
        def cutoff_max(self):
            """
            Element cutoff_max ftype=real(real12) pytype=float
            
            
            Defined at ../src/lib/mod_evolver.f90 line 37
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _raffle.f90wrap_gvector_container_type__array__cutoff_max(self._handle)
            if array_handle in self._arrays:
                cutoff_max = self._arrays[array_handle]
            else:
                cutoff_max = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _raffle.f90wrap_gvector_container_type__array__cutoff_max)
                self._arrays[array_handle] = cutoff_max
            return cutoff_max
        
        @cutoff_max.setter
        def cutoff_max(self, cutoff_max):
            self.cutoff_max[...] = cutoff_max
        
        @property
        def radius_distance_tol(self):
            """
            Element radius_distance_tol ftype=real(real12) pytype=float
            
            
            Defined at /Users/nedtaylor/DCoding/DGit/raffle/src/fortran/lib/mod_evolver.f90 \
                line 81
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _raffle.f90wrap_gvector_container_type__array__radius_distance_tol(self._handle)
            if array_handle in self._arrays:
                radius_distance_tol = self._arrays[array_handle]
            else:
                radius_distance_tol = \
                    f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _raffle.f90wrap_gvector_container_type__array__radius_distance_tol)
                self._arrays[array_handle] = radius_distance_tol
            return radius_distance_tol
        
        @radius_distance_tol.setter
        def radius_distance_tol(self, radius_distance_tol):
            self.radius_distance_tol[...] = radius_distance_tol

        @property
        def total(self):
            """
            Element total ftype=type(gvector_base_type) pytype=Gvector_Base_Type
            
            
            Defined at ../src/lib/mod_evolver.f90 line 38
            
            """
            total_handle = _raffle.f90wrap_gvector_container_type__get__total(self._handle)
            if tuple(total_handle) in self._objs:
                total = self._objs[tuple(total_handle)]
            else:
                total = evolver.gvector_base_type.from_handle(total_handle)
                self._objs[tuple(total_handle)] = total
            return total
        
        @total.setter
        def total(self, total):
            total = total._handle
            _raffle.f90wrap_gvector_container_type__set__total(self._handle, total)
        
        def init_array_system(self):
            self.system = f90wrap.runtime.FortranDerivedTypeArray(self,
                                            _raffle.f90wrap_gvector_container_type__array_getitem__system,
                                            _raffle.f90wrap_gvector_container_type__array_setitem__system,
                                            _raffle.f90wrap_gvector_container_type__array_len__system,
                                            """
            Element system ftype=type(gvector_type) pytype=Gvector_Type
            
            
            Defined at ../src/lib/mod_evolver.f90 line 39
            
            """, Evolver.gvector_type)
            return self.system
        
        def __str__(self):
            ret = ['<gvector_container_type>{\n']
            ret.append('    best_system : ')
            ret.append(repr(self.best_system))
            ret.append(',\n    best_energy : ')
            ret.append(repr(self.best_energy))
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
            ret.append(',\n    total : ')
            ret.append(repr(self.total))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = [init_array_system]
        
    
    _dt_array_initialisers = []
    

evolver = Evolver()

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
        def __init__(self, element=None, num=None, handle=None):
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
        
    
    @f90wrap.runtime.register_class("raffle.stoichiometry_type_xnum_array")
    class stoichiometry_type_xnum_array(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=stoichiometry_type_xnum_array)
        
        
        Defined at ../src/lib/mod_generator.f90 lines \
            19-21
        
        """
        def __init__(self, handle=None):
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
        
        def init_array_items(self):
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

        

        _dt_array_initialisers = [init_array_items]
        
    
    @f90wrap.runtime.register_class("raffle.raffle_generator_type")
    class raffle_generator_type(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=raffle_generator_type)
        
        
        Defined at ../src/lib/mod_generator.f90 lines \
            23-34
        
        """
        def __init__(self, handle=None):
            """
            self = Raffle_Generator_Type()
            
            
            Defined at ../src/lib/mod_generator.f90 lines \
                23-34
            
            
            Returns
            -------
            this : Raffle_Generator_Type
            	Object to be constructed
            
            
            Automatically generated constructor for raffle_generator_type
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _raffle.f90wrap_generator__raffle_generator_type_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Raffle_Generator_Type
            
            
            Defined at ../src/lib/mod_generator.f90 lines \
                23-34
            
            Parameters
            ----------
            this : Raffle_Generator_Type
            	Object to be destructed
            
            
            Automatically generated destructor for raffle_generator_type
            """
            if self._alloc:
                _raffle.f90wrap_generator__raffle_generator_type_finalise(this=self._handle)

        def set_host(self, host):
            """
            set_host__binding__raffle_generator_type(self, host)
            
            
            Defined at ../src/lib/mod_generator.f90 lines \
                99-108
            
            Parameters
            ----------
            this : unknown
            host : basis_type
            
            """
            _raffle.f90wrap_generator__set_host__binding__rgt(this=self._handle, \
                host=host._handle)
        
        def set_grid(self, grid=None, grid_spacing=None):
            """
            set_grid__binding__raffle_generator_type(self[, grid, grid_spacing])
            
            
            Defined at \
                /Users/nedtaylor/DCoding/DGit/raffle/src/fortran/lib/mod_generator.f90 lines \
                139-166
            
            Parameters
            ----------
            this : unknown
            grid : int array
            grid_spacing : float
            
            """
            _raffle.f90wrap_generator__set_grid__binding__raffle_generator_type(this=self._handle, \
                grid=grid, grid_spacing=grid_spacing)
        
        def reset_grid(self):
            """
            reset_grid__binding__raffle_generator_type(self)
            
            
            Defined at \
                /Users/nedtaylor/DCoding/DGit/raffle/src/fortran/lib/mod_generator.f90 lines \
                170-176
            
            Parameters
            ----------
            this : unknown
            
            """
            _raffle.f90wrap_generator__reset_grid__binding__raffle_generator_type(this=self._handle)
        
        def generate(self, num_structures, stoichiometry, method_probab={"void": 1.0, "walk": 1.0, "min": 1.0}, seed=None, verbose=0):
            """
            generate__binding__raffle_generator_type(self, num_structures, stoichiometry, method_probab)

            Defined at ../src/lib/mod_generator.f90 lines \
                76-84

            Parameters
            ----------
            this : unknown
            num_structures : int
            stoichiometry : stoichiometry_type_xnum_array
            method_probab : list of float
            verbose : int

            """
            
            method_probab_list = []
            method_probab_list.append(method_probab.get("void", 1.0))
            method_probab_list.append(method_probab.get("walk", 1.0))
            method_probab_list.append(method_probab.get("min", 1.0))

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
            
        def get_structures(self):
            """
            structures = get_structures__binding__raffle_generator_type(self)
            
            
            Defined at ../src/lib/mod_generator.f90 lines \
                86-97
            
            Parameters
            ----------
            this : unknown
            
            """
            structures = _raffle.f90wrap_generator__get_structures__binding__rgt(this=self._handle)
            return structures
        
        def evaluate(self, basis):
            """
            viability = evaluate__binding__raffle_generator_type(self, basis)
            
            
            Defined at ../src/lib/mod_generator.f90 lines \
                311-322
            
            Parameters
            ----------
            this : unknown
            basis : basis_type
            
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
            Element host ftype=type(basis_type) pytype=basis_type
            
            
            Defined at ../src/lib/mod_generator.f90 line \
                25
            
            """
            host_handle = _raffle.f90wrap_raffle_generator_type__get__host(self._handle)
            if tuple(host_handle) in self._objs:
                host = self._objs[tuple(host_handle)]
            else:
                host = rw_geom.basis_type.from_handle(host_handle)
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
        def distributions(self):
            """
            Element distributions ftype=type(gvector_container_type) \
                pytype=Gvector_Container_Type
            
            
            Defined at ../src/lib/mod_generator.f90 line \
                27
            
            """
            distributions_handle = \
                _raffle.f90wrap_raffle_generator_type__get__distributions(self._handle)
            if tuple(distributions_handle) in self._objs:
                distributions = self._objs[tuple(distributions_handle)]
            else:
                distributions = evolver.gvector_container_type.from_handle(distributions_handle)
                self._objs[tuple(distributions_handle)] = distributions
            return distributions
        
        @distributions.setter
        def distributions(self, distributions):
            distributions = distributions._handle
            _raffle.f90wrap_raffle_generator_type__set__distributions(self._handle, \
                distributions)
        
        @property
        def method_probab(self):
            """
            Element method_probab ftype=real(real12) pytype=float
            
            
            Defined at ../src/lib/mod_generator.f90 line \
                28
            
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
        
        def init_array_structures(self):
            self.structures = f90wrap.runtime.FortranDerivedTypeArray(self,
                                            _raffle.f90wrap_raffle_generator_type__array_getitem__structures,
                                            _raffle.f90wrap_raffle_generator_type__array_setitem__structures,
                                            _raffle.f90wrap_raffle_generator_type__array_len__structures,
                                            """
            Element items ftype=type(basis_type) pytype=basis_type
            
            
            Defined at ../src/lib/mod_generator.f90 line \
                29
            
            """, Rw_Geom.basis_type)
            return self.structures

        def __str__(self):
            ret = ['<raffle_generator_type>{\n']
            ret.append('    num_structures : ')
            ret.append(repr(self.num_structures))
            ret.append(',\n    host : ')
            ret.append(repr(self.host))
            ret.append(',\n    grid : ')
            ret.append(repr(self.grid))
            ret.append(',\n    grid_spacing : ')
            ret.append(repr(self.grid_spacing))
            ret.append(',\n    distributions : ')
            ret.append(repr(self.distributions))
            ret.append(',\n    method_probab : ')
            ret.append(repr(self.method_probab))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = [init_array_structures]
        
    
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

