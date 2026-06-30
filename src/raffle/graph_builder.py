from __future__ import print_function, absolute_import, division
import raffle._raffle as _raffle
import f90wrap.runtime
import logging
import numpy

class Graph_Builder(f90wrap.runtime.FortranModule):
    """
    Module raffle__graph_builder
    """
    @f90wrap.runtime.register_class("raffle.topology")
    class topology(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=topology)
        """
        def __init__(self, handle=None):
            """
            self = Topology_Type()

            Returns
            -------
            this : Topology_Type
            	Object to be constructed


            Automatically generated constructor for topology
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = \
                _raffle.f90wrap_raffle__graph_builder__topology_type_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result

        def __del__(self):
            """
            Destructor for class Topology_Type

            Parameters
            ----------
            this : Topology_Type
            	Object to be destructed


            Automatically generated destructor for topology
            """
            if self._alloc:
                _raffle.f90wrap_raffle__graph_builder__topology_type_finalise(this=self._handle)

        def allocate_arrays(self, num_atoms, num_pairs, num_triplets, \
            num_quadruplets):
            """
            allocate_arrays__binding__topology_type(self, num_atoms, num_pairs, \
                num_triplets, num_quadruplets)

            Parameters
            ----------
            this : Topology_Type
            num_atoms : int
            num_pairs : int
            num_triplets : int
            num_quadruplets : int

            """
            _raffle.f90wrap_raffle__graph_builder__allocate_arrays__binding__toda88(this=self._handle, \
                num_atoms=num_atoms, num_pairs=num_pairs, \
                num_triplets=num_triplets, num_quadruplets=num_quadruplets)

        def finalize(self):
            """
            finalize__binding__topology_type(self)

            Parameters
            ----------
            this : Topology_Type

            """
            _raffle.f90wrap_raffle__graph_builder__finalize__binding__topology_type(this=self._handle)

        @property
        def symbols(self):
            """
            Element symbols ftype=character(len=3) pytype=str
            """
            # get the number of species
            num_species = _raffle.f90wrap_topology_type__get__num_species(self._handle)
            symbols = numpy.zeros((num_species,), dtype='S3')
            _raffle.f90wrap_topology_type__get_symbols(self._handle, symbols=symbols)

            for i in range(num_species):
                symbols[i] = symbols[i].decode('utf-8').strip()
            return symbols

        @symbols.setter
        def symbols(self, symbols):
            _raffle.f90wrap_topology_type__set_symbols(self._handle, symbols=symbols)

        @property
        def species_index(self):
            """
            Element species_index ftype=integer pytype=int
            """
            array_ndim, array_type, array_shape, array_handle = \
                _raffle.f90wrap_topology_type__array__species_index(self._handle)
            if array_handle in self._arrays:
                species_index = self._arrays[array_handle]
            else:
                species_index = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _raffle.f90wrap_topology_type__array__species_index)
                self._arrays[array_handle] = species_index
            return species_index

        @species_index.setter
        def species_index(self, species_index):
            self.species_index[...] = species_index

        @property
        def atomic_numbers(self):
            """
            Element atomic_numbers ftype=real(real32) pytype=float
            """
            array_ndim, array_type, array_shape, array_handle = \
                _raffle.f90wrap_topology_type__array__atomic_numbers(self._handle)
            if array_handle in self._arrays:
                atomic_numbers = self._arrays[array_handle]
            else:
                atomic_numbers = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _raffle.f90wrap_topology_type__array__atomic_numbers)
                self._arrays[array_handle] = atomic_numbers
            return atomic_numbers

        @atomic_numbers.setter
        def atomic_numbers(self, atomic_numbers):
            self.atomic_numbers[...] = atomic_numbers

        @property
        def covalent_radii(self):
            """
            Element covalent_radii ftype=real(real32) pytype=float
            """
            array_ndim, array_type, array_shape, array_handle = \
                _raffle.f90wrap_topology_type__array__covalent_radii(self._handle)
            if array_handle in self._arrays:
                covalent_radii = self._arrays[array_handle]
            else:
                covalent_radii = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _raffle.f90wrap_topology_type__array__covalent_radii)
                self._arrays[array_handle] = covalent_radii
            return covalent_radii

        @covalent_radii.setter
        def covalent_radii(self, covalent_radii):
            self.covalent_radii[...] = covalent_radii

        @property
        def pair_image_shift(self):
            """
            Element pair_image_shift ftype=integer pytype=int
            """
            array_ndim, array_type, array_shape, array_handle = \
                _raffle.f90wrap_topology_type__array__pair_image_shift(self._handle)
            if array_handle in self._arrays:
                pair_image_shift = self._arrays[array_handle]
            else:
                pair_image_shift = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _raffle.f90wrap_topology_type__array__pair_image_shift)
                self._arrays[array_handle] = pair_image_shift
            return pair_image_shift

        @pair_image_shift.setter
        def pair_image_shift(self, pair_image_shift):
            self.pair_image_shift[...] = pair_image_shift

        @property
        def pair_target_species_index(self):
            """
            Element pair_target_species_index ftype=integer pytype=int
            """
            array_ndim, array_type, array_shape, array_handle = \
                _raffle.f90wrap_topology_type__array__pair_target_species_index(self._handle)
            if array_handle in self._arrays:
                pair_target_species_index = self._arrays[array_handle]
            else:
                pair_target_species_index = \
                    f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _raffle.f90wrap_topology_type__array__pair_target_species_index)
                self._arrays[array_handle] = pair_target_species_index
            return pair_target_species_index

        @pair_target_species_index.setter
        def pair_target_species_index(self, pair_target_species_index):
            self.pair_target_species_index[...] = pair_target_species_index

        @property
        def pair_index(self):
            """
            Element pair_index ftype=integer pytype=int
            """
            array_ndim, array_type, array_shape, array_handle = \
                _raffle.f90wrap_topology_type__array__pair_index(self._handle)
            if array_handle in self._arrays:
                pair_index = self._arrays[array_handle]
            else:
                pair_index = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _raffle.f90wrap_topology_type__array__pair_index)
                self._arrays[array_handle] = pair_index
            return pair_index

        @pair_index.setter
        def pair_index(self, pair_index):
            self.pair_index[...] = pair_index

        @property
        def pair_type_index(self):
            """
            Element pair_type_index ftype=integer pytype=int
            """
            array_ndim, array_type, array_shape, array_handle = \
                _raffle.f90wrap_topology_type__array__pair_type_index(self._handle)
            if array_handle in self._arrays:
                pair_type_index = self._arrays[array_handle]
            else:
                pair_type_index = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _raffle.f90wrap_topology_type__array__pair_type_index)
                self._arrays[array_handle] = pair_type_index
            return pair_type_index

        @pair_type_index.setter
        def pair_type_index(self, pair_type_index):
            self.pair_type_index[...] = pair_type_index

        @property
        def pair_cutoff_weight_3body(self):
            """
            Element pair_cutoff_weight_3body ftype=real(real32) pytype=float
            """
            array_ndim, array_type, array_shape, array_handle = \
                _raffle.f90wrap_topology_type__array__pair_cutoff_weight_3body(self._handle)
            if array_handle in self._arrays:
                pair_cutoff_weight_3body = self._arrays[array_handle]
            else:
                pair_cutoff_weight_3body = \
                    f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _raffle.f90wrap_topology_type__array__pair_cutoff_weight_3body)
                self._arrays[array_handle] = pair_cutoff_weight_3body
            return pair_cutoff_weight_3body

        @pair_cutoff_weight_3body.setter
        def pair_cutoff_weight_3body(self, pair_cutoff_weight_3body):
            self.pair_cutoff_weight_3body[...] = pair_cutoff_weight_3body

        @property
        def pair_cutoff_weight_4body(self):
            """
            Element pair_cutoff_weight_4body ftype=real(real32) pytype=float
            """
            array_ndim, array_type, array_shape, array_handle = \
                _raffle.f90wrap_topology_type__array__pair_cutoff_weight_4body(self._handle)
            if array_handle in self._arrays:
                pair_cutoff_weight_4body = self._arrays[array_handle]
            else:
                pair_cutoff_weight_4body = \
                    f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _raffle.f90wrap_topology_type__array__pair_cutoff_weight_4body)
                self._arrays[array_handle] = pair_cutoff_weight_4body
            return pair_cutoff_weight_4body

        @pair_cutoff_weight_4body.setter
        def pair_cutoff_weight_4body(self, pair_cutoff_weight_4body):
            self.pair_cutoff_weight_4body[...] = pair_cutoff_weight_4body

        @property
        def triplet_species_index(self):
            """
            Element triplet_species_index ftype=integer pytype=int
            """
            array_ndim, array_type, array_shape, array_handle = \
                _raffle.f90wrap_topology_type__array__triplet_species_index(self._handle)
            if array_handle in self._arrays:
                triplet_species_index = self._arrays[array_handle]
            else:
                triplet_species_index = \
                    f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _raffle.f90wrap_topology_type__array__triplet_species_index)
                self._arrays[array_handle] = triplet_species_index
            return triplet_species_index

        @triplet_species_index.setter
        def triplet_species_index(self, triplet_species_index):
            self.triplet_species_index[...] = triplet_species_index

        @property
        def triplet_index(self):
            """
            Element triplet_index ftype=integer pytype=int
            """
            array_ndim, array_type, array_shape, array_handle = \
                _raffle.f90wrap_topology_type__array__triplet_index(self._handle)
            if array_handle in self._arrays:
                triplet_index = self._arrays[array_handle]
            else:
                triplet_index = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _raffle.f90wrap_topology_type__array__triplet_index)
                self._arrays[array_handle] = triplet_index
            return triplet_index

        @triplet_index.setter
        def triplet_index(self, triplet_index):
            self.triplet_index[...] = triplet_index

        @property
        def triplet_pair_ids(self):
            """
            Element triplet_pair_ids ftype=integer pytype=int
            """
            array_ndim, array_type, array_shape, array_handle = \
                _raffle.f90wrap_topology_type__array__triplet_pair_ids(self._handle)
            if array_handle in self._arrays:
                triplet_pair_ids = self._arrays[array_handle]
            else:
                triplet_pair_ids = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _raffle.f90wrap_topology_type__array__triplet_pair_ids)
                self._arrays[array_handle] = triplet_pair_ids
            return triplet_pair_ids

        @triplet_pair_ids.setter
        def triplet_pair_ids(self, triplet_pair_ids):
            self.triplet_pair_ids[...] = triplet_pair_ids

        @property
        def triplet_center_index(self):
            """
            Element triplet_center_index ftype=integer pytype=int
            """
            array_ndim, array_type, array_shape, array_handle = \
                _raffle.f90wrap_topology_type__array__triplet_center_index(self._handle)
            if array_handle in self._arrays:
                triplet_center_index = self._arrays[array_handle]
            else:
                triplet_center_index = \
                    f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _raffle.f90wrap_topology_type__array__triplet_center_index)
                self._arrays[array_handle] = triplet_center_index
            return triplet_center_index

        @triplet_center_index.setter
        def triplet_center_index(self, triplet_center_index):
            self.triplet_center_index[...] = triplet_center_index

        @property
        def quadruplet_pair_ids(self):
            """
            Element quadruplet_pair_ids ftype=integer pytype=int
            """
            array_ndim, array_type, array_shape, array_handle = \
                _raffle.f90wrap_topology_type__array__quadruplet_pair_ids(self._handle)
            if array_handle in self._arrays:
                quadruplet_pair_ids = self._arrays[array_handle]
            else:
                quadruplet_pair_ids = \
                    f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _raffle.f90wrap_topology_type__array__quadruplet_pair_ids)
                self._arrays[array_handle] = quadruplet_pair_ids
            return quadruplet_pair_ids

        @quadruplet_pair_ids.setter
        def quadruplet_pair_ids(self, quadruplet_pair_ids):
            self.quadruplet_pair_ids[...] = quadruplet_pair_ids

        @property
        def quadruplet_species_index(self):
            """
            Element quadruplet_species_index ftype=integer pytype=int
            """
            array_ndim, array_type, array_shape, array_handle = \
                _raffle.f90wrap_topology_type__array__quadruplet_species_index(self._handle)
            if array_handle in self._arrays:
                quadruplet_species_index = self._arrays[array_handle]
            else:
                quadruplet_species_index = \
                    f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _raffle.f90wrap_topology_type__array__quadruplet_species_index)
                self._arrays[array_handle] = quadruplet_species_index
            return quadruplet_species_index

        @quadruplet_species_index.setter
        def quadruplet_species_index(self, quadruplet_species_index):
            self.quadruplet_species_index[...] = quadruplet_species_index

        @property
        def num_atoms(self):
            """
            Element num_atoms ftype=integer  pytype=int
            """
            return \
                _raffle.f90wrap_topology_type__get__num_atoms(self._handle)

        @num_atoms.setter
        def num_atoms(self, num_atoms):
            _raffle.f90wrap_topology_type__set__num_atoms(self._handle, \
                num_atoms)

        @property
        def num_pairs(self):
            """
            Element num_pairs ftype=integer  pytype=int
            """
            return \
                _raffle.f90wrap_topology_type__get__num_pairs(self._handle)

        @num_pairs.setter
        def num_pairs(self, num_pairs):
            _raffle.f90wrap_topology_type__set__num_pairs(self._handle, \
                num_pairs)

        @property
        def num_triplets(self):
            """
            Element num_triplets ftype=integer  pytype=int
            """
            return \
                _raffle.f90wrap_topology_type__get__num_triplets(self._handle)

        @num_triplets.setter
        def num_triplets(self, num_triplets):
            _raffle.f90wrap_topology_type__set__num_triplets(self._handle, \
                num_triplets)

        @property
        def num_quadruplets(self):
            """
            Element num_quadruplets ftype=integer  pytype=int
            """
            return \
                _raffle.f90wrap_topology_type__get__num_quadruplets(self._handle)

        @num_quadruplets.setter
        def num_quadruplets(self, num_quadruplets):
            _raffle.f90wrap_topology_type__set__num_quadruplets(self._handle, \
                num_quadruplets)

        def __str__(self):
            ret = ['<topology>{\n']
            ret.append('    symbols : ')
            ret.append(repr(self.symbols))
            ret.append(',\n    species_index : ')
            ret.append(repr(self.species_index))
            ret.append(',\n    atomic_numbers : ')
            ret.append(repr(self.atomic_numbers))
            ret.append(',\n    covalent_radii : ')
            ret.append(repr(self.covalent_radii))
            ret.append(',\n    pair_image_shift : ')
            ret.append(repr(self.pair_image_shift))
            ret.append(',\n    pair_target_species_index : ')
            ret.append(repr(self.pair_target_species_index))
            ret.append(',\n    pair_index : ')
            ret.append(repr(self.pair_index))
            ret.append(',\n    pair_type_index : ')
            ret.append(repr(self.pair_type_index))
            ret.append(',\n    pair_cutoff_weight_3body : ')
            ret.append(repr(self.pair_cutoff_weight_3body))
            ret.append(',\n    pair_cutoff_weight_4body : ')
            ret.append(repr(self.pair_cutoff_weight_4body))
            ret.append(',\n    triplet_index : ')
            ret.append(repr(self.triplet_index))
            ret.append(',\n    triplet_pair_ids : ')
            ret.append(repr(self.triplet_pair_ids))
            ret.append(',\n    triplet_center_index : ')
            ret.append(repr(self.triplet_center_index))
            ret.append(',\n    quadruplet_pair_ids : ')
            ret.append(repr(self.quadruplet_pair_ids))
            ret.append(',\n    quadruplet_species_index : ')
            ret.append(repr(self.quadruplet_species_index))
            ret.append(',\n    num_atoms : ')
            ret.append(repr(self.num_atoms))
            ret.append(',\n    num_pairs : ')
            ret.append(repr(self.num_pairs))
            ret.append(',\n    num_triplets : ')
            ret.append(repr(self.num_triplets))
            ret.append(',\n    num_quadruplets : ')
            ret.append(repr(self.num_quadruplets))
            ret.append('}')
            return ''.join(ret)

        _dt_array_initialisers = []


    @f90wrap.runtime.register_class("raffle.graph_tensors")
    class graph_tensors(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=graph_tensors)
        """
        def __init__(self, handle=None):
            """
            self = Graph_Tensors_Type()

            Returns
            -------
            this : Graph_Tensors_Type
            	Object to be constructed


            Automatically generated constructor for graph_tensors
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            result = \
                _raffle.f90wrap_raffle__graph_builder__graph_tensors_type_initialise()
            self._handle = result[0] if isinstance(result, tuple) else result

        def __del__(self):
            """
            Destructor for class Graph_Tensors_Type

            Parameters
            ----------
            this : Graph_Tensors_Type
            	Object to be destructed


            Automatically generated destructor for graph_tensors
            """
            if self._alloc:
                _raffle.f90wrap_raffle__graph_builder__graph_tensors_type_finalise(this=self._handle)

        @property
        def global_features(self):
            """
            Element global_features ftype=real(real32) pytype=float
            """
            array_ndim, array_type, array_shape, array_handle = \
                _raffle.f90wrap_graph_tensors_type__array__global_features(self._handle)
            if array_handle in self._arrays:
                global_features = self._arrays[array_handle]
            else:
                global_features = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _raffle.f90wrap_graph_tensors_type__array__global_features)
                self._arrays[array_handle] = global_features
            return global_features

        @global_features.setter
        def global_features(self, global_features):
            self.global_features[...] = global_features

        @property
        def atom_node_features(self):
            """
            Element atom_node_features ftype=real(real32) pytype=float
            """
            array_ndim, array_type, array_shape, array_handle = \
                _raffle.f90wrap_graph_tensors_type__array__atom_node_features(self._handle)
            if array_handle in self._arrays:
                atom_node_features = self._arrays[array_handle]
            else:
                atom_node_features = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _raffle.f90wrap_graph_tensors_type__array__atom_node_features)
                self._arrays[array_handle] = atom_node_features
            return atom_node_features

        @atom_node_features.setter
        def atom_node_features(self, atom_node_features):
            self.atom_node_features[...] = atom_node_features

        @property
        def pair_node_features(self):
            """
            Element pair_node_features ftype=real(real32) pytype=float
            """
            array_ndim, array_type, array_shape, array_handle = \
                _raffle.f90wrap_graph_tensors_type__array__pair_node_features(self._handle)
            if array_handle in self._arrays:
                pair_node_features = self._arrays[array_handle]
            else:
                pair_node_features = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _raffle.f90wrap_graph_tensors_type__array__pair_node_features)
                self._arrays[array_handle] = pair_node_features
            return pair_node_features

        @pair_node_features.setter
        def pair_node_features(self, pair_node_features):
            self.pair_node_features[...] = pair_node_features

        @property
        def atom_edge_index(self):
            """
            Element atom_edge_index ftype=integer pytype=int
            """
            array_ndim, array_type, array_shape, array_handle = \
                _raffle.f90wrap_graph_tensors_type__array__atom_edge_index(self._handle)
            if array_handle in self._arrays:
                atom_edge_index = self._arrays[array_handle]
            else:
                atom_edge_index = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _raffle.f90wrap_graph_tensors_type__array__atom_edge_index)
                self._arrays[array_handle] = atom_edge_index
            return atom_edge_index

        @atom_edge_index.setter
        def atom_edge_index(self, atom_edge_index):
            self.atom_edge_index[...] = atom_edge_index

        @property
        def pair_edge_index(self):
            """
            Element pair_edge_index ftype=integer pytype=int
            """
            array_ndim, array_type, array_shape, array_handle = \
                _raffle.f90wrap_graph_tensors_type__array__pair_edge_index(self._handle)
            if array_handle in self._arrays:
                pair_edge_index = self._arrays[array_handle]
            else:
                pair_edge_index = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _raffle.f90wrap_graph_tensors_type__array__pair_edge_index)
                self._arrays[array_handle] = pair_edge_index
            return pair_edge_index

        @pair_edge_index.setter
        def pair_edge_index(self, pair_edge_index):
            self.pair_edge_index[...] = pair_edge_index

        @property
        def pair_index(self):
            """
            Element pair_index ftype=integer pytype=int
            """
            array_ndim, array_type, array_shape, array_handle = \
                _raffle.f90wrap_graph_tensors_type__array__pair_index(self._handle)
            if array_handle in self._arrays:
                pair_index = self._arrays[array_handle]
            else:
                pair_index = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _raffle.f90wrap_graph_tensors_type__array__pair_index)
                self._arrays[array_handle] = pair_index
            return pair_index

        @pair_index.setter
        def pair_index(self, pair_index):
            self.pair_index[...] = pair_index

        @property
        def atom_edge_attr(self):
            """
            Element atom_edge_attr ftype=real(real32) pytype=float
            """
            array_ndim, array_type, array_shape, array_handle = \
                _raffle.f90wrap_graph_tensors_type__array__atom_edge_attr(self._handle)
            if array_handle in self._arrays:
                atom_edge_attr = self._arrays[array_handle]
            else:
                atom_edge_attr = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _raffle.f90wrap_graph_tensors_type__array__atom_edge_attr)
                self._arrays[array_handle] = atom_edge_attr
            return atom_edge_attr

        @atom_edge_attr.setter
        def atom_edge_attr(self, atom_edge_attr):
            self.atom_edge_attr[...] = atom_edge_attr

        @property
        def pair_edge_attr(self):
            """
            Element pair_edge_attr ftype=real(real32) pytype=float
            """
            array_ndim, array_type, array_shape, array_handle = \
                _raffle.f90wrap_graph_tensors_type__array__pair_edge_attr(self._handle)
            if array_handle in self._arrays:
                pair_edge_attr = self._arrays[array_handle]
            else:
                pair_edge_attr = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _raffle.f90wrap_graph_tensors_type__array__pair_edge_attr)
                self._arrays[array_handle] = pair_edge_attr
            return pair_edge_attr

        @pair_edge_attr.setter
        def pair_edge_attr(self, pair_edge_attr):
            self.pair_edge_attr[...] = pair_edge_attr

        @property
        def atom_edge_weight(self):
            """
            Element atom_edge_weight ftype=real(real32) pytype=float
            """
            array_ndim, array_type, array_shape, array_handle = \
                _raffle.f90wrap_graph_tensors_type__array__atom_edge_weight(self._handle)
            if array_handle in self._arrays:
                atom_edge_weight = self._arrays[array_handle]
            else:
                atom_edge_weight = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _raffle.f90wrap_graph_tensors_type__array__atom_edge_weight)
                self._arrays[array_handle] = atom_edge_weight
            return atom_edge_weight

        @atom_edge_weight.setter
        def atom_edge_weight(self, atom_edge_weight):
            self.atom_edge_weight[...] = atom_edge_weight

        @property
        def pair_edge_weight(self):
            """
            Element pair_edge_weight ftype=real(real32) pytype=float
            """
            array_ndim, array_type, array_shape, array_handle = \
                _raffle.f90wrap_graph_tensors_type__array__pair_edge_weight(self._handle)
            if array_handle in self._arrays:
                pair_edge_weight = self._arrays[array_handle]
            else:
                pair_edge_weight = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _raffle.f90wrap_graph_tensors_type__array__pair_edge_weight)
                self._arrays[array_handle] = pair_edge_weight
            return pair_edge_weight

        @pair_edge_weight.setter
        def pair_edge_weight(self, pair_edge_weight):
            self.pair_edge_weight[...] = pair_edge_weight

        @property
        def pair_hyperedge_index(self):
            """
            Element pair_hyperedge_index ftype=integer pytype=int
            """
            array_ndim, array_type, array_shape, array_handle = \
                _raffle.f90wrap_graph_tensors_type__array__pair_hyperedge_index(self._handle)
            if array_handle in self._arrays:
                pair_hyperedge_index = self._arrays[array_handle]
            else:
                pair_hyperedge_index = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _raffle.f90wrap_graph_tensors_type__array__pair_hyperedge_index)
                self._arrays[array_handle] = pair_hyperedge_index
            return pair_hyperedge_index

        @pair_hyperedge_index.setter
        def pair_hyperedge_index(self, pair_hyperedge_index):
            self.pair_hyperedge_index[...] = pair_hyperedge_index

        @property
        def pair_hyperedge_weight(self):
            """
            Element pair_hyperedge_weight ftype=real(real32) pytype=float
            """
            array_ndim, array_type, array_shape, array_handle = \
                _raffle.f90wrap_graph_tensors_type__array__pair_hyperedge_weight(self._handle)
            if array_handle in self._arrays:
                pair_hyperedge_weight = self._arrays[array_handle]
            else:
                pair_hyperedge_weight = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _raffle.f90wrap_graph_tensors_type__array__pair_hyperedge_weight)
                self._arrays[array_handle] = pair_hyperedge_weight
            return pair_hyperedge_weight

        @pair_hyperedge_weight.setter
        def pair_hyperedge_weight(self, pair_hyperedge_weight):
            self.pair_hyperedge_weight[...] = pair_hyperedge_weight

        @property
        def pair_hyperedge_attr(self):
            """
            Element pair_hyperedge_attr ftype=real(real32) pytype=float
            """
            array_ndim, array_type, array_shape, array_handle = \
                _raffle.f90wrap_graph_tensors_type__array__pair_hyperedge_attr(self._handle)
            if array_handle in self._arrays:
                pair_hyperedge_attr = self._arrays[array_handle]
            else:
                pair_hyperedge_attr = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                        self._handle,
                                        _raffle.f90wrap_graph_tensors_type__array__pair_hyperedge_attr)
                self._arrays[array_handle] = pair_hyperedge_attr
            return pair_hyperedge_attr

        @pair_hyperedge_attr.setter
        def pair_hyperedge_attr(self, pair_hyperedge_attr):
            self.pair_hyperedge_attr[...] = pair_hyperedge_attr

        def __str__(self):
            ret = ['<graph_tensors>{\n']
            ret.append('    global_features : ')
            ret.append(repr(self.global_features))
            ret.append(',\n    atom_node_features : ')
            ret.append(repr(self.atom_node_features))
            ret.append(',\n    pair_node_features : ')
            ret.append(repr(self.pair_node_features))
            ret.append(',\n    atom_edge_index : ')
            ret.append(repr(self.atom_edge_index))
            ret.append(',\n    pair_edge_index : ')
            ret.append(repr(self.pair_edge_index))
            ret.append(',\n    pair_index : ')
            ret.append(repr(self.pair_index))
            ret.append(',\n    atom_edge_attr : ')
            ret.append(repr(self.atom_edge_attr))
            ret.append(',\n    pair_edge_attr : ')
            ret.append(repr(self.pair_edge_attr))
            ret.append(',\n    atom_edge_weight : ')
            ret.append(repr(self.atom_edge_weight))
            ret.append(',\n    pair_edge_weight : ')
            ret.append(repr(self.pair_edge_weight))
            ret.append(',\n    pair_hyperedge_index : ')
            ret.append(repr(self.pair_hyperedge_index))
            ret.append(',\n    pair_hyperedge_weight : ')
            ret.append(repr(self.pair_hyperedge_weight))
            ret.append(',\n    pair_hyperedge_attr : ')
            ret.append(repr(self.pair_hyperedge_attr))
            ret.append('}')
            return ''.join(ret)

        _dt_array_initialisers = []

graph_builder = Graph_Builder()
