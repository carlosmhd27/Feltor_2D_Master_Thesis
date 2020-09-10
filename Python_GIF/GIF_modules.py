import os
import matplotlib
matplotlib.use("Agg")
import sys

sys.path.insert(1, '../2D_FELTOR_Analysis/')
from Analysis import Analyse
from matplotlib.pyplot    import axis, savefig, pcolormesh, imshow, show, figure
from matplotlib.pyplot    import subplot, title, plot, colorbar

import numpy                as np
import matplotlib.animation as animation
import matplotlib.pyplot    as plt

ftsz_title = 15
ftsz_label = 12


def Analyzed (File_name):
    '''
    A function to get the parameters we want to plot, CM, error in Mass, Velocity of the
    CM, Potential, density and vorticity-
    '''
    
    Analitics = Analyse(File_name)
    Analitics.Mass      = Analitics.integrate('ions')
    Analitics.Potential = Analitics.integrate('potential')
    Analitics.Mass_err  = (Analitics.Mass[0] - Analitics.Mass) / Analitics.Mass[0]
    
    return Analitics

def animate (position, Analytics, ax):
    '''
    A function to update each frame
    '''
        
    ax[0, 0].plot(Analytics.X_CM[:position+1], Analytics.Y_CM[:position+1], 'o', color='tab:blue');
    ax[0, 0].plot(Analytics.X_CM[position], Analytics.Y_CM[position], 'o', color = 'red')

    ax[0, 1].plot(Analytics.time[:position+1], Analytics.Mass_err[:position+1], color='tab:blue'); 
#     ax[0, 1].plot(Analytics.time[point], Analytics.Mass_err[point], 'o', color = 'red')
    
    ax[0, 2].plot(Analytics.V_CM_x[:position+1], Analytics.V_CM_y[:position+1], 'o', color='tab:blue');
    ax[0, 2].plot(Analytics.V_CM_x[position], Analytics.V_CM_y[position], 'o', color = 'red')

    ax[1, 0].pcolormesh(Analytics.x, Analytics.y, Analytics.potential[position], cmap = 'hot', shading = 'gouraud');
    
    ax[1, 1].pcolormesh(Analytics.x, Analytics.y, Analytics.ions[position], cmap = 'hot', shading = 'gouraud');

    ax[1, 2].pcolormesh(Analytics.x, Analytics.y, Analytics.Data['vorticity'][position], cmap = 'hot', shading = 'gouraud');
    return ax


def init(Analytics, ax, fig, model):
    '''
    Initialization of the plots, set the limits in the graphics, labels, titles and scales
    '''
    
    fig.suptitle(model, fontsize=20)
    
    ax[0, 0].set_xlabel('CM X axis', fontsize = ftsz_label)
    ax[0, 0].set_ylabel('CM Y axis', fontsize = ftsz_label)
    ax[0, 0].set_xlim(np.amin(Analytics.X_CM) - 0.5, np.amax(Analytics.X_CM) + 0.5)
    ax[0, 0].set_ylim(np.amin(Analytics.Y_CM) - 0.5, np.amax(Analytics.Y_CM) + 0.5)
    ax[0, 0].set_title('Center of mass', fontsize = ftsz_title)
    
    ax[0, 1].set_xlabel('Time', fontsize = ftsz_label)
    ax[0, 1].set_ylabel('Error in total Mass', fontsize = ftsz_label)
    ax[0, 1].set_xlim(np.amin(Analytics.time) - 0.5, np.amax(Analytics.time) + 0.5)
    ax[0, 1].set_ylim(np.amin(Analytics.Mass_err), np.amax(Analytics.Mass_err) + 0.01)
    ax[0, 1].set_xscale('symlog')
    ax[0, 1].set_title('Error in total Mass', fontsize = ftsz_title)
    
    ax[0, 2].set_xlabel('V_CM X axis', fontsize = ftsz_label)
    ax[0, 2].set_ylabel('V_CM Y axis', fontsize = ftsz_label)
    ax[0, 2].set_xlim(np.amin(Analytics.V_CM_x) - 0.5, np.amax(Analytics.V_CM_x) + 0.5)
    ax[0, 2].set_ylim(np.amin(Analytics.V_CM_y) - 0.5, np.amax(Analytics.V_CM_y) + 0.5)
    ax[0, 2].set_title('Velocity of the CM', fontsize = ftsz_title)
    
    ax[1, 0].set_title('Potential', fontsize = ftsz_title)
    ax[1, 1].set_title('Ions density', fontsize = ftsz_title)
    ax[1, 2].set_title('Vorticity', fontsize = ftsz_title)
    
#     im1 = ax[1, 0].pcolormesh(Analytics.x, Analytics.y, Analytics.potential[0], cmap = 'plasma', shading = 'gouraud');
#     fig.colorbar(im1,ax=ax[1, 0])
    
#     im2 = ax[1, 1].pcolormesh(Analytics.x, Analytics.y, Analytics.ions[0], cmap = 'plasma', shading = 'gouraud');
#     fig.colorbar(im2,ax=ax[1, 1])
    
#     im3 = ax[1, 2].pcolormesh(Analytics.x, Analytics.y, Analytics.Data['vorticity'][0], cmap = 'plasma', shading = 'gouraud');
#     fig.colorbar(im3,ax=ax[1, 2])

    return ax