from netCDF4           import Dataset
from scipy.io          import netcdf
from scipy.integrate   import simps
from boutdata          import collect
from numpy             import tile, array, copy, ndim, shape
from matplotlib.pyplot import imshow, show, plot, colorbar, title, savefig

import numpy as np
import matplotlib.pyplot as plt
import json


class Analyse ():
    
    
    def __init__ (self, File_name):
        self.Data  = Dataset(File_name, 'r', format="NETCDF4")
        self.input = json.loads(self.Data.inputfile)
        
        self.x    = array(self.Data['x'][:])
        self.Nx   = self.input['Nx'] + self.input['Nx_out'] ## Nx_out are the points out of the grid

        self.y    = array(self.Data['y'][:])
        self.Ny   = self.input['Ny'] + self.input['Ny_out']
        
        self.ions = array(self.Data['ions'][:])
        self.CM()

    def CM(self):

        x_vec, y_vec = copy(self.x), copy(self.y)
        rho          = copy(self.ions)
        
        X_mat = tile(x_vec.reshape(self.Nx, 1), self.Ny)
        Y_mat = tile(y_vec.reshape(self.Ny, 1), self.Nx).transpose()

        n_total = np.sum(rho, axis = (1,2))
        
        self.X_CM = np.sum(X_mat * rho, axis = (1,2)) / n_total 
        self.Y_CM = np.sum(Y_mat * rho, axis = (1,2)) / n_total
     
    def imshow_variable(self, variable_name, time_step, set_color_bar = False, save = None, **kwargs):
        
        parameter = self.Data[variable_name][time_step]
        imshow(parameter, **kwargs)
        if set_color_bar:
            colorbar()
        if save != None and type(save) == str:
            savefig(save)
        show();
    
    def integrate(self, variable, dim_integral = 2, axis = 1):
        
        if type(variable) == str:
            var = copy(self.Data[variable][:])
        else:
            var = copy(var)
            
        dim_var = ndim(var)
        if dim_var == 3 and dim_integral == 2:
            Integral = simps(var, self.x, axis = axis)
            Integral = simps(Integral, self.y, axis = axis)
            
        else:
            if dim_var == 1:
                axis = 0
            Integral = simps(var, self.x, axis = axis)
        
        return Integral
