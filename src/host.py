from structure_class import Structure
from child import Child
from ase.io import write, Trajectory

class Host(Structure): 

    """
    Base class for the host object 
    
    This object contains any functions which one might want to execute on the host object i.e get strain. 
    These functions could contain calls to the f90 implementation 
    """
    def __init__(self, atoms, struct_id, **kwargs):
        super(Host, self).__init__(atoms)
        self.struct_id = struct_id
        self.get_search_volume(**kwargs)
        self.children = []
        

    def get_id(self):
        return self.struct_id

    def get_host_energy(self):
        if self.structure.calc is not None: 
            return structure.get_potential_energy()
        else: 
            print("No Calculator for Host structure")
            exit()

    def get_search_volume(self, **kwargs):
        if "volume_method" in kwargs:
            method=_determine_volume_method(**kwargs)
        else:
            method="manual"

        exec(f'self._{method}_calc(**kwargs)')

    def _lom_calc(self,**kwargs):
        self.PNI_error()
        
    def _manual_calc(self,**kwargs):
        volume=kwargs.get("manual_value")
        if volume is None: 
            self.dict_error("manual_value")
        self.volume = volume

    def _cell_calc(self,**kwargs):
        self.PNI_error()
        


    def add_child(self,atoms):
        child=Child(atoms)
        self.children.append(child)


    def calculate_children(self,calc,path="tests/"): 
        calculated_children=[]
        energies = []
        for i,child in enumerate(self.children): 
            if child.structure.calc is not calc:
                child.add_calculator(self.structure.calc)
            child.get_energy()
            energies.append(child.get_energy)
            traj=Trajectory(f'{path}{i}.traj',mode="w")
            traj.write(child.structure)
        return energies

    def dump_children(self,path="tests/"):
        self.calculate_children(path)
        self.children=[]
        
    
def _determine_volume_method(**kwargs):
    for key, value in kwargs.items(): 
        if key == "volume_method": 
            for val in ["manual","cell_based","lom"]: 
                if value == val: 
                    return value
    




