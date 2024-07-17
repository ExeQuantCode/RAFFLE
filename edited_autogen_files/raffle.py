from __future__ import print_function, absolute_import, division
import _raffle
import f90wrap.runtime
import logging
import numpy

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
            result = _raffle.f90wrap_generator__stoichiometry_type_initialise()
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
                _raffle.f90wrap_generator__stoichiometry_type_finalise(this=self._handle)
        
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
            Element items ftype=type(test_type) pytype=Test_Type
            
            
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
        
        def print_hello(self):
            """
            print_hello__binding__raffle_generator_type(self)
            
            
            Defined at ../src/lib/mod_generator.f90 lines \
                69-74
            
            Parameters
            ----------
            this : unknown
            
            """
            _raffle.f90wrap_generator__print_hello__binding__raffle_generator_type(this=self._handle)

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
                method_probab=method_probab)#, n0=len(method_probab))
        
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
        def lattice_host(self):
            """
            Element lattice_host ftype=real(real12) pytype=float
            
            
            Defined at ../src/lib/mod_generator.f90 line \
                25
            
            """
            array_ndim, array_type, array_shape, array_handle = \
                _raffle.f90wrap_raffle_generator_type__array__lattice_host(self._handle)
            if array_handle in self._arrays:
                lattice_host = self._arrays[array_handle]
            else:
                lattice_host = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _raffle.f90wrap_raffle_generator_type__array__lattice_host)
                self._arrays[array_handle] = lattice_host
            return lattice_host
        
        @lattice_host.setter
        def lattice_host(self, lattice_host):
            self.lattice_host[...] = lattice_host
        
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
        
        def __str__(self):
            ret = ['<raffle_generator_type>{\n']
            ret.append('    bins : ')
            ret.append(repr(self.bins))
            ret.append(',\n    lattice_host : ')
            ret.append(repr(self.lattice_host))
            ret.append(',\n    method_probab : ')
            ret.append(repr(self.method_probab))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
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

