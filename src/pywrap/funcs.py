from ase import Atoms
import numpy as np
from raffle import raffle 
from ase.io import read

def get_stoichiometry(stoichiometry_dict):
    symbols=[key for key in stoichiometry_dict]
    configs = []
    atoms=[]
    for i in range(len(stoichiometry_dict[symbols[0]])):
        symbol_ase=""
        for symbol in symbols:
            symbol_ase+=(f'{symbol}{stoichiometry_dict[symbol][i]}')
        atoms.append(Atoms(symbol_ase))
    return(atoms)

def check_bondlengths(trial, atoms): 
    for atom in atoms: 
        if get_distance(atom, trial) < 0.1: 
            return False 
    return True 


def get_distance(a,b):
    return np.linalg.norm(a.position-b.position)

def _initialise_generator():
    generator = raffle.generator.raffle_generator_type()
    return generator

def _populate_database(**kwargs):
    database_path=kwargs.get("database")
    if database_path is not None:
        database = read(kwargs.get("database"),index=":")
    else:
        #database = read("../../examples/example_files/database/database.xyz", index=":")
        from ase.build import bulk
        from ase.build import fcc111
        database=[bulk("Au"),fcc111("Au",size=(2,2,1))]

    database_basis = raffle.rw_geom.bas_type_xnum_array()
    database_basis.allocate(len(database))
    for i, atoms in enumerate(database):
        database_basis.items[i].fromase(atoms)

    return database_basis




def RAFFLE(host_structures,atoms,energies_dict,iterations,**kwargs): 
    

    database_basis=_populate_database()

    for host in host_structures:
        generator=_initialise_generator()
        generator.set_host(raffle.rw_geom.bas_type(host.structure))
        generator.distributions.set_element_energies(energies_dict)
        generator.distributions.create(database_basis)
        generator.distributions.get_element_energies()

        generator.bins=[50,50,50]
        stoich_list = raffle.generator.stoichiometry_type_xnum_array()
        list_stoich=[]
        print(atoms)
        k=0
        for i, structure in enumerate(atoms):
            for atom in structure: 
                list_stoich.append(atom.symbol)
                
        stoich_list.allocate(len(list_stoich))
        for i, item in enumerate(list_stoich):
            stoich_list.items[i].element=atom.symbol
            stoich_list.items[i].num=1
        
        generator.generate(num_structures=iterations,stoichiometry=stoich_list,method_probab=[0.0,0.0,1.0])
        
        for structure in generator.structures:
            
            ase_structure=structure.toase()
            ase_structure.calc=host.structure.calc
            
            host.add_child(ase_structure)
        



def rss(hosts,atoms,iterations):
    
    """
    Execute RAFFLE call
    """

    _RSS(hosts,atoms,iterations)
    
def _RSS(hosts,atoms,iterations): 
    for parent_id, host in enumerate(hosts): 
        for i in range(iterations):
            for master_group in atoms:
                atom_group=master_group
                for atom in atom_group:
                    for tries in range(1000):
                        test_atom = atom
                        test_atom.position = np.matmul(np.random.rand(1,3),host.structure.cell)
                        if check_bondlengths(test_atom,host.structure): 
                            atom = test_atom 
                            break
            host.add_child(host.structure+atom_group)            
          
