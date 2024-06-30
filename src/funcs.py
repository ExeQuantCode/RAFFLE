from ase import Atoms
import numpy as np


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


def raffle(hosts,atoms,iterations):
    
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
          
