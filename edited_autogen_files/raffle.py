from __future__ import print_function, absolute_import, division
import _raffle
import f90wrap.runtime
import logging
import numpy
from ase import Atoms

class Rw_Geom(f90wrap.runtime.FortranModule):
    """
    Module rw_geom
    
    
    Defined at ../src/lib/mod_rw_geom.f90 lines \
        13-968
    
    """
    @f90wrap.runtime.register_class("raffle.spec_type")
    class spec_type(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=spec_type)
        
        
        Defined at ../src/lib/mod_rw_geom.f90 lines \
            26-32
        
        """
        def __init__(self, handle=None):
            """
            self = Spec_Type()
            
            
            Defined at ../src/lib/mod_rw_geom.f90 lines \
                26-32
            
            
            Returns
            -------
            this : Spec_Type
            	Object to be constructed
            
            
            Automatically generated constructor for spec_type
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _raffle.f90wrap_rw_geom__spec_type_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class Spec_Type
            
            
            Defined at ../src/lib/mod_rw_geom.f90 lines \
                26-32
            
            Parameters
            ----------
            this : Spec_Type
            	Object to be destructed
            
            
            Automatically generated destructor for spec_type
            """
            if self._alloc:
                _raffle.f90wrap_rw_geom__spec_type_finalise(this=self._handle)
        
        @property
        def atom(self):
            """
            Element atom ftype=real(real12) pytype=float
            
            
            Defined at ../src/lib/mod_rw_geom.f90 line 27
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _raffle.f90wrap_spec_type__array__atom(self._handle)
            if array_handle in self._arrays:
                atom = self._arrays[array_handle]
            else:
                atom = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _raffle.f90wrap_spec_type__array__atom)
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
            return _raffle.f90wrap_spec_type__get__mass(self._handle)
        
        @mass.setter
        def mass(self, mass):
            _raffle.f90wrap_spec_type__set__mass(self._handle, mass)
        
        @property
        def charge(self):
            """
            Element charge ftype=real(real12) pytype=float
            
            
            Defined at ../src/lib/mod_rw_geom.f90 line 29
            
            """
            return _raffle.f90wrap_spec_type__get__charge(self._handle)
        
        @charge.setter
        def charge(self, charge):
            _raffle.f90wrap_spec_type__set__charge(self._handle, charge)
        
        @property
        def name(self):
            """
            Element name ftype=character(len=3) pytype=str
            
            
            Defined at ../src/lib/mod_rw_geom.f90 line 30
            
            """
            return _raffle.f90wrap_spec_type__get__name(self._handle)
        
        @name.setter
        def name(self, name):
            _raffle.f90wrap_spec_type__set__name(self._handle, name)
        
        @property
        def num(self):
            """
            Element num ftype=integer  pytype=int
            
            
            Defined at ../src/lib/mod_rw_geom.f90 line 31
            
            """
            return _raffle.f90wrap_spec_type__get__num(self._handle)
        
        @num.setter
        def num(self, num):
            _raffle.f90wrap_spec_type__set__num(self._handle, num)
        
        def __str__(self):
            ret = ['<spec_type>{\n']
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
        
    
    @f90wrap.runtime.register_class("raffle.bas_type")
    class bas_type(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=bas_type)
        
        
        Defined at ../src/lib/mod_rw_geom.f90 lines \
            34-42
        
        """
        def __init__(self, atoms=None, handle=None):
            """
            self = bas_type()
            
            
            Defined at ../src/lib/mod_rw_geom.f90 lines \
                34-42
            
            
            Returns
            -------
            this : bas_type
            	Object to be constructed
            
            
            Automatically generated constructor for bas_type
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _raffle.f90wrap_rw_geom__bas_type_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result

            if atoms is not None:
                self.fromase(atoms)
        
        def __del__(self):
            """
            Destructor for class bas_type
            
            
            Defined at ../src/lib/mod_rw_geom.f90 lines \
                34-42
            
            Parameters
            ----------
            this : bas_type
            	Object to be destructed
            
            
            Automatically generated destructor for bas_type
            """
            if self._alloc:
                _raffle.f90wrap_rw_geom__bas_type_finalise(this=self._handle)
        
        def allocate_species(self, num_species=None, species_symbols=None, species_count=None, \
            atoms=None):
            """
            allocate_species__binding__bas_type(self[, num_species, species_symbols, \
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
            _raffle.f90wrap_rw_geom__allocate_species__binding__bas_type(this=self._handle, \
                num_species=num_species, species_symbols=species_symbols, species_count=species_count, \
                atoms=atoms)
        
        def init_array_spec(self):
            self.spec = f90wrap.runtime.FortranDerivedTypeArray(self,
                                            _raffle.f90wrap_bas_type__array_getitem__spec,
                                            _raffle.f90wrap_bas_type__array_setitem__spec,
                                            _raffle.f90wrap_bas_type__array_len__spec,
                                            """
            Element spec ftype=type(spec_type) pytype=Spec_Type
            
            
            Defined at ../src/lib/mod_rw_geom.f90 line 35
            
            """, Rw_Geom.spec_type)
            return self.spec
        
        def toase(self):

            # Set the species list
            positions = []
            species_string = ""
            for i in range(self.nspec):
                for j in range(self.spec[i].num):
                    species_string += str(self.spec[i].name.decode()).strip()
                    positions.append(self.spec[i].atom[j])
            
            # Set the atoms
            atoms = Atoms(species_string, positions)
            atoms.set_pbc(self.pbc)

            # Set the lattice vectors
            atoms.set_cell(bas.lat)

            return atoms
        
        def fromase(self, atoms):
            
            # Get the species symbols
            species_symbols = atoms.get_chemical_symbols()
            species_symbols_unique = sorted(set(species_symbols))

            # Set the number of species
            self.nspec = len(species_symbols_unique)
            
            # Set the number of atoms
            self.natom = len(atoms)
            
            # Set the energy
            # self.energy = atoms.get_total_energy()
            
            # # Set the lattice vectors
            self.lat = numpy.reshape(atoms.get_cell().flatten(), [3,3], order='F')
            self.pbc = atoms.pbc
            
            # Set the system name
            self.sysname = atoms.get_chemical_formula()
            
            # Set the species list
            species_count = []
            atom_positions = []
            positions = atoms.get_positions()
            for species in species_symbols_unique:
                species_count.append(sum([1 for symbol in species_symbols if symbol == species]))
                for j, symbol in enumerate(species_symbols):
                    if symbol == species:
                        atom_positions.append(positions[j])
            
            # Allocate memory for the atom list
            self.allocate_species(species_symbols=species_symbols_unique, species_count=species_count, atoms=atom_positions)

        @property
        def nspec(self):
            """
            Element nspec ftype=integer  pytype=int
            
            
            Defined at ../src/lib/mod_rw_geom.f90 line 36
            
            """
            return _raffle.f90wrap_bas_type__get__nspec(self._handle)
        
        @nspec.setter
        def nspec(self, nspec):
            _raffle.f90wrap_bas_type__set__nspec(self._handle, nspec)
        
        @property
        def natom(self):
            """
            Element natom ftype=integer  pytype=int
            
            
            Defined at ../src/lib/mod_rw_geom.f90 line 37
            
            """
            return _raffle.f90wrap_bas_type__get__natom(self._handle)
        
        @natom.setter
        def natom(self, natom):
            _raffle.f90wrap_bas_type__set__natom(self._handle, natom)
        
        @property
        def energy(self):
            """
            Element energy ftype=real(real12) pytype=float
            
            
            Defined at ../src/lib/mod_rw_geom.f90 line 38
            
            """
            return _raffle.f90wrap_bas_type__get__energy(self._handle)
        
        @energy.setter
        def energy(self, energy):
            _raffle.f90wrap_bas_type__set__energy(self._handle, energy)
        
        @property
        def lat(self):
            """
            Element lat ftype=real(real12) pytype=float
            
            
            Defined at ../src/lib/mod_rw_geom.f90 line 38
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _raffle.f90wrap_bas_type__array__lat(self._handle)
            if array_handle in self._arrays:
                lat = self._arrays[array_handle]
            else:
                lat = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _raffle.f90wrap_bas_type__array__lat)
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
            return _raffle.f90wrap_bas_type__get__lcart(self._handle)
        
        @lcart.setter
        def lcart(self, lcart):
            _raffle.f90wrap_bas_type__set__lcart(self._handle, lcart)
        
        @property
        def pbc(self):
            """
            Element pbc ftype=logical pytype=bool
            
            
            Defined at ../src/lib/mod_rw_geom.f90 line 40
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _raffle.f90wrap_bas_type__array__pbc(self._handle)
            if array_handle in self._arrays:
                pbc = self._arrays[array_handle]
            else:
                pbc = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _raffle.f90wrap_bas_type__array__pbc)
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
            return _raffle.f90wrap_bas_type__get__sysname(self._handle)

        @sysname.setter
        def sysname(self, sysname):
            _raffle.f90wrap_bas_type__set__sysname(self._handle, sysname)
                
        def __str__(self):
            ret = ['<bas_type>{\n']
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
        


    @f90wrap.runtime.register_class("raffle.bas_type_xnum_array")
    class bas_type_xnum_array(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=bas_type_xnum_array)
        
        
        Defined at ../src/lib/mod_generator.f90 lines \
            19-21
        
        """
        def __init__(self, handle=None):
            """
            self = bas_Type()
            
            
            Defined at ../src/lib/mod_generator.f90 lines \
                19-21
            
            
            Returns
            -------
            this : bas_Type
            	Object to be constructed
            
            
            Automatically generated constructor for bas_type
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = _raffle.f90wrap_rw_geom__bas_type_xnum_array_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result
        
        def __del__(self):
            """
            Destructor for class bas_Type
            
            
            Defined at ../src/lib/mod_generator.f90 lines \
                19-21
            
            Parameters
            ----------
            this : bas_type
            	Object to be destructed
            
            
            Automatically generated destructor for bas_type
            """
            if self._alloc:
                _raffle.f90wrap_rw_geom__bas_type_xnum_array_finalise(this=self._handle)
        
        def init_array_items(self):
            self.items = f90wrap.runtime.FortranDerivedTypeArray(self,
                                            _raffle.f90wrap_bas_type_xnum_array__array_getitem__items,
                                            _raffle.f90wrap_bas_type_xnum_array__array_setitem__items,
                                            _raffle.f90wrap_bas_type_xnum_array__array_len__items,
                                            """
            Element items ftype=type(bas_type) pytype=bas_type
            
            
            Defined at  line 0
            
            """, Rw_Geom.bas_type)
            return self.items
        
        def allocate(self, size):
            """
            Allocate the items array with the given size
            
            Parameters
            ----------
            self : bas_type
            size : int
                Size of the items array
            """
            _raffle.f90wrap_bas_type_xnum_array__array_alloc__items(self._handle, num=size)

        def deallocate(self):
            """
            Deallocate the items array
            """
            _raffle.f90wrap_bas_type_xnum_array__array_dealloc__items(self._handle)

        
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
    #     bas : bas_type
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
    #     bas : bas_type
        
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
    #     inbas : bas_type
    #     latconv : float array
        
    #     Returns
    #     -------
    #     outbas : bas_type
        
    #     """
    #     outbas = _raffle.f90wrap_rw_geom__convert_bas(inbas=self._handle, \
    #         latconv=latconv)
    #     outbas = f90wrap.runtime.lookup_class("raffle.bas_type").from_handle(outbas, \
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
    #     inbas : bas_type
    #     outbas : bas_type
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
            result = _raffle.f90wrap_stoichiometry_type_initialise()
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
            host : bas_type
            
            """
            _raffle.f90wrap_generator__set_host__binding__rgt(this=self._handle, \
                host=host._handle)
        
        def generate(self, num_structures, stoichiometry, method_probab=[1.0, 1.0, 1.0]):
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

            """
            
            _raffle.f90wrap_generator__generate__binding__rgt(
                this=self._handle,
                num_structures=num_structures,
                stoichiometry=stoichiometry._handle,
                method_probab=method_probab)
        
        def evaluate(self, basis):
            """
            viability = evaluate__binding__raffle_generator_type(self, basis)
            
            
            Defined at ../src/lib/mod_generator.f90 lines \
                311-322
            
            Parameters
            ----------
            this : unknown
            basis : bas_type
            
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
            Element host ftype=type(bas_type) pytype=bas_type
            
            
            Defined at ../src/lib/mod_generator.f90 line \
                25
            
            """
            host_handle = _raffle.f90wrap_raffle_generator_type__get__host(self._handle)
            if tuple(host_handle) in self._objs:
                host = self._objs[tuple(host_handle)]
            else:
                host = rw_geom.bas_type.from_handle(host_handle)
                self._objs[tuple(host_handle)] = host
            return host
        
        @host.setter
        def host(self, host):
            host = host._handle
            _raffle.f90wrap_raffle_generator_type__set__host(self._handle, host)
        
        @property
        def bins(self):
            """
            Element bins ftype=integer pytype=int
            
            
            Defined at ../src/lib/mod_generator.f90 line \
                24
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _raffle.f90wrap_raffle_generator_type__array__bins(self._handle)
            if array_handle in self._arrays:
                bins = self._arrays[array_handle]
            else:
                bins = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _raffle.f90wrap_raffle_generator_type__array__bins)
                self._arrays[array_handle] = bins
            return bins
        
        @bins.setter
        def bins(self, bins):
            self.bins[...] = bins
        
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
            Element items ftype=type(bas_type) pytype=bas_type
            
            
            Defined at ../src/lib/mod_generator.f90 line \
                29
            
            """, Rw_Geom.bas_type)
            return self.structures

        def __str__(self):
            ret = ['<raffle_generator_type>{\n']
            ret.append('    num_structures : ')
            ret.append(repr(self.num_structures))
            ret.append(',\n    host : ')
            ret.append(repr(self.host))
            ret.append(',\n    bins : ')
            ret.append(repr(self.bins))
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

