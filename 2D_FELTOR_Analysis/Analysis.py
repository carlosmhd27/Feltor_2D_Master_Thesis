from netCDF4           import Dataset
from json              import loads
from scipy.integrate   import simps
from numpy             import tile, array, copy, ndim, shape, arange, roll, sqrt, conjugate
from numpy             import correlate, average, empty, array_equal, squeeze, transpose
from matplotlib.pyplot import pcolormesh, show, plot, colorbar, title, savefig

import numpy as np


class Analyse ():
    ''' 
    Class to open and analyse the output data from the 2D FELTOR simulations
    to run it one must introduce the name of the output File
    '''
    
    def __init__ (self, File_name, Access_Mode = 'r', parallel = False):
        '''
        Open the data and extract some important parameters
        Also calculate the Center of Mass
        '''
        self.Data  = Dataset(File_name, Access_Mode, format="NETCDF4", parallel = parallel)
        self.input = loads(self.Data.inputfile)
        
        self.x  = array(self.Data['x'][:])
        self.y  = array(self.Data['y'][:])
        self.lx, self.ly  = self.input['lx'], self.input['ly']
        
        self.time    = array(self.Data['time'][:])
        self.dt      = self.time[1] - self.time[0]
        
        self.Nx = self.input['n_out'] * self.input['Nx_out'] 
        self.Ny = self.input['n_out'] * self.input['Ny_out']
        self.nt = len(self.time)
        
        self.ions      = array(self.Data['ions'][:])
        self.potential = array(self.Data['potential'][:])
        self.vorticity = array(self.Data['vorticity'][:])
        
        self.Mass      = self.integrate('ions') / (self.lx * self.ly)
        self.Potential = self.integrate('potential') / (self.lx * self.ly)
        self.Mass_err  = (self.Mass[0] - self.Mass) / self.Mass[0]
#         self.CM(); self.V_CM()

    def CM(self):
        '''
        Function to calculate the Center of Mass of the system for all the time steps
        The Center of Mass is a 2D vector, for that reason we obtain 2 coordinate, X and Y,
        separately
        '''
        
        X_mat = tile(self.x.reshape(self.Nx, 1), self.Ny)
        Y_mat = tile(self.y.reshape(self.Ny, 1), self.Nx).transpose()
        
        self.X_CM = self.integrate(X_mat * self.ions) / self.Mass 
        self.Y_CM = self.integrate(Y_mat * self.ions) / self.Mass
        
    def V_CM(self):
        '''
        Function to calculate the velocity of the Center of Mass of the system 
        for all the time steps
        The velocidty of the Center of Mass is a 2D vector, for that reason we
        obtain 2 coordinate, X and Y, separately.
        In this case it is calculated as V = (x_i - x_i-1) / (t_i - t_i-1)
        For that reason V(0) = (0, 0), this prevent as from having an arbitrary
        speed for the las position.
        '''
        
        
        D_x, D_y = self.X_CM - roll(self.X_CM, 1), self.Y_CM - roll(self.Y_CM, 1)
        D_t      = self.time - roll(self.time, 1)
        
        self.V_CM_x,     self.V_CM_y    = D_x / D_t, D_y / D_t
        self.V_CM_x[0] = self.V_CM_y[0] = 0    
    
    def integrate(self, variable, dim_integral = 2, axis = 1, typ = 't'):
        '''
        Function to integrate a variable over the grid, this variable could be
        a vector or a matrix, or even been integrated over only on dimension or two
        '''
        
        if type(variable) == str:
            Integral = copy(self.Data[variable][:])
        else:
            Integral = copy(variable)
            
        dim_var = ndim(Integral)
        
        if dim_var == 3:
            indep_vars = [self.time, self.x, self.y]  
        elif dim_var == 2:
            indep_vars = [self.x, self.y] 
        else:
            indep_vars = [self.time] if typ == 't' else [self.y] if typ == 'y' else [self.x]
        
        if dim_var == 3 and dim_integral == 2:
            Integral = simps(Integral, indep_vars[axis], axis = axis)
            if axis == 2:
                axis = 1
            Integral = simps(Integral, indep_vars[axis + 1], axis = axis)
            
        else:
            if dim_var == 1:
                axis  = 0 
            Integral = simps(Integral, indep_vars[axis], axis = axis)
        
        return Integral
    
    def pcolormesh_variable(self, variable_name, time_step, set_color_bar = False, save = None, **kwargs):
        '''
        Easy function to plot 2D variables such ions or vorticity
        '''
        
        if type(variable) == str:
            parameter = copy(self.Data[variable][time_step])
        else:
            parameter = copy(variable)
            
        pcolormesh(self.x, self.y, parameter, **kwargs)
        if set_color_bar:
            colorbar()
        if save != None and type(save) == str:
            savefig(save)
        show();
    
    def plot_parameter(self, variable, x = None, save = None, **kwargs):
        '''
        Easy function to plot time dependance scalar variables such CM or the total Mass
        '''
        if type(variable) == str:
            parameter = copy(self.Data[variable][time_step])
        else:
            parameter = copy(variable)
        
        if ndim(parameter) != 1:
            raise ValueError('The dimension of the variable are different than 1, try to use imshow.')
        
        if type(x) == type(None):
            if len(parameter) == self.tm_stps:
                x = self.time
            else:
                x = arange(len(parameter))
                
        plot(x, parameter, **kwargs)
        
        if save != None and type(save) == str:
            savefig(save)
        show();
   
    
    def c_corr_dt (self, f, g, time_units = 1):
        '''
        cross-correlation time-delay
        '''

        assert array_equal(shape(f), shape(g)), 'The dimensions of f and g should be equal'
        
        ## + 1 because we wanna do [-time_units, time_units]
        steps  = int(time_units / self.dt) + 1
        corr   = empty((2 * steps + 1, self.Nx, self.Ny))
        nt     = len(f)
        g_conj = conjugate(g)

        corr  = empty((2 * steps + 1, self.Nx, self.Ny))
        for i in range(-steps, steps + 1):
            corr[steps + i, :, :] = np.sum((f * roll(g_conj, i, axis=0))[steps : -steps], axis=0) / (nt - 2 * steps)

        corr = self.integrate(corr, 1, axis = 2) / self.ly

        return corr

    def c_corr_sp (self, f, g, x0=None, y0=None, integrate = False):
        '''
        Spatial cross-correlation, x and y should not be integrated
        '''
        
        assert array_equal(shape(f), shape(g)), 'The dimensions of f and g should be equal'


        if type(x0) != type(None) and type(y0) == type(None):
            f_cp, g_cp     = copy(f), conjugate(g)
            z, axis        = copy(x0), 0
            transpose_back = False

        elif type(x0) == type(None) and type(y0) != type(None):
            f_cp, g_cp     = transpose(f, (0, 2, 1)), conjugate(transpose(g, (0, 2, 1)))
            z, axis        = copy(y0), 1
            transpose_back = True

        else:
            raise TypeError('One of the variables should be None, the other and int or an array of ints')

        nt, nx, ny = shape(f_cp)
        mz         = shape(z)[0]
        
        if mz == 0:
            z = [z]
            
        f_cp, g_cp = f_cp[:, z, :], g_cp[:, z, :]
                
        steps = int(ny / 2)
        fcg = empty((nt, mz, 2 * steps + 1))
        for i in range(-steps, steps + 1):
            fcg[:,:, steps + i] = self.integrate((f_cp * roll(g_cp,i,axis=2)),1,axis=2) 
            fcg[:,:, steps + i] /= self.ly

        fcg = average(fcg, axis = 0)

        if transpose_back:
            fcg = transpose(fcg) 

        if integrate:
            l   = [self.lx, self.ly][axis]
            fcg = self.integrate(fcg, axis = axis) / l

        return squeeze(fcg)
    
    
class Essential():
    """ Empty class, used to create ad hoc objects """
    def __init__ (self, File_name):
        self.Data  = Dataset(File_name, format="NETCDF4")
        
        self.x  = copy(array(self.Data['x'][:]))
        self.y  = copy(array(self.Data['y'][:]))
        
        self.time    = copy(array(self.Data['time'][:]))
        
        self.ions      = copy(array(self.Data['ions'][:]))
        self.potential = copy(array(self.Data['potential'][:]))
        self.vorticity = copy(array(self.Data['vorticity'][:]))
        
        self.Data.close()
        
        def c_corr_dt (self, f, g, time_units = 1):
        '''
        cross-correlation time-delay
        '''

        assert array_equal(shape(f), shape(g)), 'The dimensions of f and g should be equal'
        
        ## + 1 because we wanna do [-time_units, time_units]
        steps  = int(time_units / self.dt) + 1
        corr   = empty((2 * steps + 1, self.Nx, self.Ny))
        nt     = len(f)
        g_conj = conjugate(g)

        corr  = empty((2 * steps + 1, self.Nx, self.Ny))
        for i in range(-steps, steps + 1):
            corr[steps + i, :, :] = np.sum((f * roll(g_conj, i, axis=0))[steps : -steps], axis=0) / (nt - 2 * steps)

        corr = self.integrate(corr, 1, axis = 2) / self.ly

        return corr

    def c_corr_sp (self, f, g, x0=None, y0=None, integrate = False):
        '''
        Spatial cross-correlation, x and y should not be integrated
        '''
        
        assert array_equal(shape(f), shape(g)), 'The dimensions of f and g should be equal'


        if type(x0) != type(None) and type(y0) == type(None):
            f_cp, g_cp     = copy(f), conjugate(g)
            z, axis        = copy(x0), 0
            transpose_back = False

        elif type(x0) == type(None) and type(y0) != type(None):
            f_cp, g_cp     = transpose(f, (0, 2, 1)), conjugate(transpose(g, (0, 2, 1)))
            z, axis        = copy(y0), 1
            transpose_back = True

        else:
            raise TypeError('One of the variables should be None, the other and int or an array of ints')

        nt, nx, ny = shape(f_cp)
        mz         = shape(z)[0]
        
        if mz == 0:
            z = [z]
            
        f_cp, g_cp = f_cp[:, z, :], g_cp[:, z, :]
                
        steps = int(ny / 2)
        fcg = empty((nt, mz, 2 * steps + 1))
        for i in range(-steps, steps + 1):
            fcg[:,:, steps + i] = self.integrate((f_cp * roll(g_cp,i,axis=2)),1,axis=2) 
            fcg[:,:, steps + i] /= self.ly

        fcg = average(fcg, axis = 0)

        if transpose_back:
            fcg = transpose(fcg) 

        if integrate:
            l   = [self.lx, self.ly][axis]
            fcg = self.integrate(fcg, axis = axis) / l

        return squeeze(fcg)

    def integrate(self, variable, dim_integral = 2, axis = 1):
        '''
        Function to integrate a variable over the grid, this variable could be
        a vector or a matrix, or even been integrated over only on dimension or two
        '''
        
        if type(variable) == str:
            var = copy(self.Data[variable][:])
        else:
            var = copy(variable)
            
        dim_var = ndim(var)
        if dim_var == 3 and dim_integral == 2:
            Integral = simps(var, self.x, axis = axis)
            Integral = simps(Integral, self.y, axis = axis)
            
        else:
            if dim_var == 1:
                axis  = 0
            Integral = simps(var, [self.x, self.y][axis], axis = axis)
        
        return Integral
