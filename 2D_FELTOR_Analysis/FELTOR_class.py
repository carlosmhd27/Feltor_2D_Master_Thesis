import json
from subprocess import call
from numpy      import copy

from IPython.display import clear_output

data_example = {
                "model" : "HW",
                "modified": 1,
                "n" :  3,
                "Nx" : 120,
                "Ny" : 120,
                "dt" : 0.01,
                "n_out"  : 3,
                "Nx_out" : 60,
                "Ny_out" : 60,
                "itstp"  : 10,
                "maxout" : 80,
                "eps_pol"   : 1e-6,
                "eps_time"  : 1e-10,
                "stages" : 3,
                "curvature"  : 0.015,
                "dens_prof"  : -0.015,
                "adiabatic"  : 0.0001,
                "nu_perp"  : 5e-3,
                "amplitude"  : 0.5,
                "sigma"  : 10,
                "posX"  : 0.5,
                "posY"  : 0.5,
                "lx"  : 200,
                "ly"  : 200,
                "bc_x"  : "DIR",
                "bc_y"  : "PER",
                "rows": 1,
                "cols": 2,
                "width": 1000,
                "height": 500}

class FELTOR_Cpp ():

    def __init__(self, PATH_FELTOR = '', FELTOR_file = 'convection_hpc', 
                 input_file = 'input.json', output_file = 'output.nc', data = None):
        
        self.PATH_FELTOR = PATH_FELTOR
        self.FELTOR_file = FELTOR_file
        self.input_file  = input_file
        self.output_file = output_file
        
        if data != None:
            self.data = data
        else:
            self.data = data_example.copy()
            
    def save_input (self):
        with open(self.input_file, 'w') as outfile:
            json.dump(self.data, outfile, indent = 1)
    
    def run(self, save = True):
        if save:
            self.save_input()
        program = self.PATH_FELTOR + self.FELTOR_file
        inp     = self.input_file
        out     = self.output_file
#         call([program, inp, out])
        return program + ' ' + inp + ' ' + out
    
    
