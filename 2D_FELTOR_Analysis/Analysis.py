from netCDF4           import Dataset
from json              import loads
from scipy.integrate   import simps
from numpy             import tile, array, copy, ndim, shape, arange, roll, sqrt 
from numpy             import correlate, average, empty, array_equal, squeeze, transpose
from matplotlib.pyplot import pcolormesh, show, plot, colorbar, title, savefig

import numpy as np

class Analyse ():
    ''' 
    Class to ope and analyse the output data from the 2D FELTOR simulations
    to run it one must introduce the name of the output File
    '''
    
    def __init__ (self, File_name):
        '''
        Open the data and extract some important parameters
        Also calculate the Center of Mass
        '''
        self.Data  = Dataset(File_name, 'r', format="NETCDF4")
        self.input = loads(self.Data.inputfile)
        
        self.x  = array(self.Data['x'][:])
        self.y  = array(self.Data['y'][:])
        ## Nx_out are the points out of the grid
        self.Nx = self.input['n_out'] * self.input['Nx_out'] 
        self.Ny = self.input['n_out'] * self.input['Ny_out']
        
        self.time    = array(self.Data['time'][:])
        self.tm_stps = len(self.time)
        
        self.ions      = array(self.Data['ions'][:])
        self.potential = array(self.Data['potential'][:])
        self.vorticity = array(self.Data['vorticity'][:])
        
        self.Mass      = self.integrate('ions')
        self.Potential = self.integrate('potential')
        self.Mass_err  = (self.Mass[0] - self.Mass) / self.Mass[0]
        self.CM(); self.V_CM()

    def CM(self):
        '''
        Function to calculate the Center of Mass of the system for all the time steps
        The Center of Mass is a 2D vector, for that reason we obtain 2 coordinate, X and Y,
        separately
        '''
        
        x_vec, y_vec = copy(self.x), copy(self.y)
        rho          = copy(self.ions)
        
        X_mat = tile(x_vec.reshape(self.Nx, 1), self.Ny)
        Y_mat = tile(y_vec.reshape(self.Ny, 1), self.Nx).transpose()

        n_total = np.sum(rho, axis = (1,2))
        
        self.X_CM = np.sum(X_mat * rho, axis = (1,2)) / n_total 
        self.Y_CM = np.sum(Y_mat * rho, axis = (1,2)) / n_total
        
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
        
        X, Y, t = copy(self.X_CM), copy(self.Y_CM), copy(self.time)
        D_x, D_y, D_t = X - roll(X, 1), Y - roll(Y, 1), t - roll(t, 1)
        self.V_CM_x, self.V_CM_y      = D_x / D_t, D_y / D_t
        self.V_CM_x[0] = self.V_CM_y[0] = 0
        
    def c_corr_dt (self, f, g):
        '''
        cross-correlation time-delay
        '''
        if not array_equal(shape(f), shape(g)):
            raise Exception('The dimensions of f and g should be equal')
        
        nt, nx, ny = shape(f)
        corr = empty((nt, nx, ny))
    
        for ix in range(nx):
            for iy in range(ny):
                corr [:, ix, iy] = correlate(f[:, ix, iy], g[:, ix, iy], 'same')
            
        corr = self.integrate(corr) / ((self.x[-1] - self.x[0]) * (self.y[-1] - self.y[0]))

        return corr
    
    def c_corr_sp (self, f, g, x0=None, y0=None, integrate = False, axis = 1):
        '''
        Spatial cross-correlation, x and y should not be integrated
        '''

        if not array_equal(shape(f), shape(g)):
            raise Exception('The dimensions of f and g should be equal')

        if type(x0) != type(None) and type(y0) == type(None):
            f_cp, g_cp     = copy(f), copy(g)
            z, axis        = copy(x0), 0
            transpose_back = False
            
        elif type(x0) == None and type(y0) != None:
            f_cp, g_cp     = copy(transpose(f, (0, 2, 1))), copy(transpose(g, (0, 2, 1)))
            z, axis        = copy(y0), 1
            transpose_back = True

        else:
            raise TypeError('One of the variables should be None, the other and int or an array of ints')

        nt, nx, ny = shape(f_cp)
        mz         = shape(z)

        if len(mz) == 0:
            mz, z = 1, [z]
        else:
            mz = mz[0]

        fcg = empty((nt, mz, ny))
        
        for i in range(nt):
            for j in range(mz):
                fcg[i, j, :] = correlate(f_cp[i, z[j], :], g_cp[i, z[j], :], 'same')

        fcg = average(fcg, axis = 0)

        if transpose_back:
            fcg = transpose(fcg) 
        
        if integrate:
            b   = [self.x, self.y][axis][-1]
            a   = [self.x, self.y][axis][0]
            fcg = self.integrate(fcg, axis = axis) / (b - a)

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
    
        
#     def imshow_variable(self, variable_name, time_step, set_color_bar = False, save = None, **kwargs):
#         '''
#         Easy function to plot 2D variables such ions or vorticity
#         '''
        
#         if type(variable) == str:
#             parameter = copy(self.Data[variable][time_step])
#         else:
#             parameter = copy(variable)
            
#         imshow(parameter, **kwargs)
#         if set_color_bar:
#             colorbar()
#         if save != None and type(save) == str:
#             savefig(save)
#         show();
    
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
