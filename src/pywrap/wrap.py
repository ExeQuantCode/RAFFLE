import matplotlib.pyplot as plt

from ase import Atoms
from ase.io import read, write
from ase.build import surface
from host import Host
from ase.calculators.emt import EMT
from argparse import ArgumentParser
import numpy as np
import funcs 
from raffle import raffle

parser = ArgumentParser()
parser.add_argument('-i', '--process_id', type=int, default=0)
parser.add_argument('-t', '--total_submissions', type=int, default=1)
args=parser.parse_args()

np.random.seed(0)

#############################################################################
#############################################################################
"""
Here, we initialise the training database, if one is provided/desired.
Not implemented yet, not sure if we want to keep this in the fortran code
"""

database=[]
#############################################################################
#############################################################################
"""
kwargs defined for the calculation. 
"""

kwargs={"volume_method" : "manual", "manual_value" : "10"}

    
#############################################################################
#############################################################################
"""
calculator can be any ASE calculator, including VASP. EMT for example 
"""
calculator=EMT()

#############################################################################
#############################################################################

"""
Here, we can either read in or call ARTEMIS directly. I create a simple surface as an example
"""
created_by_artemis=[]
host_structures=[]
for i in range(1):
    structure=surface('Au', (2,2,1), 5)
    structure.center(vacuum = 10 ,axis = 2)
    created_by_artemis.append(structure)

for i,structures in enumerate(created_by_artemis): 
    host=Host(structures,struct_id=i,**kwargs)
    host.add_calculator(calculator)
    host_structures.append(host)

#############################################################################
#############################################################################
"""
Stoichoimetry dict can be populated by a function, but here it is initialised simply
"""
stoichiometry_dict={"Au":[15]}
iterations=1

#############################################################################
#############################################################################

"""
Create Distributions
"""

try: 
    read("example_distribution.traj",index=":")

except: 
    from ase.build import bulk 
    from ase.build import fcc111


    database=[bulk("Au"),bulk("Au")]
    for structure in database:

        structure.calc = EMT()
        structure.get_potential_energy()
    write("example_distribution.traj",database)

"""
Calculate parallelisation over array job
"""




atoms=funcs.get_stoichiometry(stoichiometry_dict)


energies_dict={'Au': -5.0}
#raffle.Generator.generate(host_structures,atoms,iterations)
#funcs.rss(host_structures,atoms,iterations)

funcs.RAFFLE(host_structures,atoms,energies_dict,iterations,**kwargs)

for host in host_structures:
    host.calculate_children(calculator)
    
    for child in host.children:
        for atom in child.structure: 
            print(atom.position)
    host.dump_children()






x=[]
y=[]
for i in range(iterations): 
    traj=read(f'../tests/{i}.traj')
    x.append(i)
    y.append(traj.get_potential_energy())

plt.plot(x,y)
plt.show()










