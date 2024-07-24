from class_method import MethodClass

class Structure(MethodClass):
    
    def __init__(self,structure):
        super().__init__()
        self.structure = structure

    def add_calculator(self,calc): 
        self.structure.calc=calc 
        

    def get_energy(self):
        if self.structure.calc is not None:
            self.structure.get_potential_energy()
        else:
            print("No Calculator for Host structure")
            exit()
