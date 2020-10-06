import matplotlib
matplotlib.use("Agg")

import sys
sys.path.i1nsert(1, '../2D_FELTOR_Analysis/')

from Analysis import Analyse
from matplotlib.pyplot    import axis, savefig, pcolormesh, imshow, show, figure
from matplotlib.pyplot    import subplot, title, plot, colorbar
from numpy    import amax, amin

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
    
    Analitics.int_vort     = Analitics.integrate('vorticity')
    Analitics.int_vort_sqr = Analitics.integrate(Analitics.vorticity ** 2)
    
    return Analitics

def animate (position, Analytics, ax):
    '''
    A function to update each frame
    '''
    
    ### Plots
    ## Integrated Vorticity
    ax[0, 0].plot(Analytics.time[:position+1], Analytics.int_vort[:position+1], color='tab:blue'); 
    
    ## Mass
    ax[0, 1].plot(Analytics.time[:position+1], Analytics.Mass[:position+1], color='tab:blue'); 
    
    ## Integrated Vorticity square ** 2
    ax[0, 2].plot(Analytics.time[:position+1], Analytics.int_vort_sqr[:position+1], color='tab:blue'); 

    ## Pcolormesh
    ## Potential
    ax[1, 0].pcolormesh(Analytics.x, Analytics.y, Analytics.potential[position], cmap = 'hot', shading = 'gouraud');
    
    ## Density
    ax[1, 1].pcolormesh(Analytics.x, Analytics.y, Analytics.ions[position], cmap = 'hot', shading = 'gouraud');

    ## Vorticity
    ax[1, 2].pcolormesh(Analytics.x, Analytics.y, Analytics.vorticity[position], cmap = 'hot', shading = 'gouraud');
    return ax


def init(Analytics, ax, fig, model):
    '''
    Initialization of the plots, set the limits in the graphics, labels, titles and scales
    '''
    
    fig.suptitle(model, fontsize=20)
    
    ## Integrated Vorticity
    ax[0, 0].set_xlabel('Time', fontsize = ftsz_label)
    ax[0, 0].set_ylabel(r'Integrated $\omega$', fontsize = ftsz_label)
    ax[0, 0].set_xlim(amin(Analytics.time) - 0.5, amax(Analytics.time) + 0.5)
    ax[0, 0].set_ylim(amin(Analytics.int_vort) - 0.5, amax(Analytics.int_vort) + 0.5)
    ax[0, 0].set_title('Integrated Vorticity', fontsize = ftsz_title)
    
    ## Total Mass
    ax[0, 1].set_xlabel('Time', fontsize = ftsz_label)
    ax[0, 1].set_ylabel('Error in total Mass', fontsize = ftsz_label)
    ax[0, 1].set_xlim(amin(Analytics.time) - 0.5, amax(Analytics.time) + 0.5)
    ax[0, 1].set_ylim(amin(Analytics.Mass), amax(Analytics.Mass) + 0.01)
    ax[0, 1].set_xscale('symlog')
    ax[0, 1].set_title('Mass', fontsize = ftsz_title)
    
    ## Integrated Vorticity square ** 2
    ax[0, 2].set_xlabel('Time', fontsize = ftsz_label)
    ax[0, 2].set_ylabel(r'Integrated $\omega^2$', fontsize = ftsz_label)
    ax[0, 2].set_xlim(amin(Analytics.time) - 0.5, amax(Analytics.time) + 0.5)
    ax[0, 2].set_ylim(amin(Analytics.int_vort_sqr) - 0.5, amax(Analytics.int_vort_sqr) + 0.5)
    ax[0, 2].set_title('Integrated Vorticity Square', fontsize = ftsz_title)
    
    ## Pcolormesh
    ## Potential, density, Vorticity
    ax[1, 0].set_title('Potential', fontsize = ftsz_title)
    ax[1, 1].set_title('Ions density', fontsize = ftsz_title)
    ax[1, 2].set_title('Vorticity', fontsize = ftsz_title)

    return ax
