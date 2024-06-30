class MethodClass():
    def __init__(self):
        pass
        
    def PNI_error(self): 
        print("Property not implemented")
        print("TERMINATING!")
        exit() 
        
    def dict_error(self, key):
        print(f'dict error, key {key} unspecified') 
        print("TERMINATING!")
        exit() 
       
