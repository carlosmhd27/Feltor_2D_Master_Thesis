import matplotlib
matplotlib.use("Agg")

import sys
sys.path.i1nsert(1, '../2D_FELTOR_Analysis/')

from Analysis          import Analyse
from matplotlib.pyplot import axis, savefig, pcolormesh, imshow, show, figure
from matplotlib.pyplot import subplot, title, plot, colorbar
from numpy             import amax, amin

import matplotlib.animation as animation
import matplotlib.pyplot    as plt

ftsz_title = 15
ftsz_label = 12
cbar       = []

min_time,         max_time         = 0, 0
min_int_vort,     max_int_vort     = 0, 0
min_Mass,         max_Mass         = 0, 0
min_int_vort_sqr, max_int_vort_sqr = 0, 0

def Analyzed (File_name):
    '''
    A function to get the parameters we want to plot, CM, error in Mass, Velocity of the
    CM, Potential, density and vorticity-
    '''
    
    global min_time, min_int_vort, min_Mass, min_int_vort_sqr
    global max_time, max_int_vort, max_Mass, max_int_vort_sqr
    
    Analytics = Analyse(File_name)
    
    Analytics.int_vort     = Analytics.integrate('vorticity')
    Analytics.int_vort_sqr = Analytics.integrate(Analytics.vorticity ** 2)
    
    min_time         = amin(Analytics.time);         max_time         = amax(Analytics.time)
    min_int_vort     = amin(Analytics.int_vort);     max_int_vort     = amax(Analytics.int_vort)
    min_Mass         = amin(Analytics.Mass);         max_Mass         = amax(Analytics.Mass)
    min_int_vort_sqr = amin(Analytics.int_vort_sqr); max_int_vort_sqr = amax(Analytics.int_vort_sqr)
    
    return Analitics

def animate (position, Analytics, ax, fig):
    '''
    A function to update each frame
    '''
    
    for i in range(3):
        for j in range(2):
            ax[j, i].clear()
    ax = init(Analytics, ax, fig, None, suptitle = False)

    ### Plots
    ## Center of mass    
    ## Integrated Vorticity
    ax[0, 0].plot(Analytics.time[:position+1], Analytics.int_vort[:position+1], color='tab:blue'); 
    
    ## Mass
    ax[0, 1].plot(Analytics.time[:position+1], Analytics.Mass[:position+1], color='tab:blue'); 

    
    ## Integrated Vorticity square ** 2
    ax[0, 2].plot(Analytics.time[:position+1], Analytics.int_vort_sqr[:position+1], color='tab:blue'); 

    ## The colormesh get updated when colling init
    ## This happens because we wanna put the colorbar from the biegining
 
    return ax


def init(Analytics, ax, fig, model, pst = 0, suptitle = True):
    '''
    Initialization of the plots, set the limits in the graphics, labels, titles and scales
    '''
    
    global cbar
    
    if suptitle:
        fig.suptitle(model, fontsize=20)
    
    ## Integrated Vorticity
    ax[0, 0].set_xlabel('Time', fontsize = ftsz_label)
    ax[0, 0].set_ylabel(r'Integrated $\omega$', fontsize = ftsz_label)
    ax[0, 0].set_xlim(min_time - 0.5, max_time + 0.5)
    ax[0, 0].set_ylim(min_int_vort - 0.5, max_int_vort + 0.5)
    ax[0, 0].set_title('Integrated Vorticity', fontsize = ftsz_title)
    
    ## Total Mass
    ax[0, 1].set_xlabel('Time', fontsize = ftsz_label)
    ax[0, 1].set_ylabel('Error in total Mass', fontsize = ftsz_label)
    ax[0, 1].set_xlim(min_time - 0.5, max_time + 0.5)
    ax[0, 1].set_ylim(min_Mass, max_Mass + 0.01)
#     ax[0, 1].set_xscale('symlog')
    ax[0, 1].set_title('Mass', fontsize = ftsz_title)
    
    ## Integrated Vorticity square ** 2
    ax[0, 2].set_xlabel('Time', fontsize = ftsz_label)
    ax[0, 2].set_ylabel(r'Integrated $\omega^2$', fontsize = ftsz_label)
    ax[0, 1].set_xlim(min_time - 0.5, max_time + 0.5)
    ax[0, 0].set_ylim(min_int_vort_sqr - 0.5, max_int_vort_sqr + 0.5)
    ax[0, 2].set_title('Integrated Vorticity Square', fontsize = ftsz_title)
    
    ## Pcolormesh
    ## Potential, density, Vorticity
    ax[1, 0].set_title('Potential', fontsize = ftsz_title)
    ax[1, 1].set_title('Ions density', fontsize = ftsz_title)
    ax[1, 2].set_title('Vorticity', fontsize = ftsz_title)
    
    im1 = ax[1, 0].pcolormesh(Analytics.x, Analytics.y, Analytics.potential[pst], cmap = 'hot', shading = 'gouraud');

    
    im2 = ax[1, 1].pcolormesh(Analytics.x, Analytics.y, Analytics.ions[pst], cmap = 'hot', shading = 'gouraud');
    
    im3 = ax[1, 2].pcolormesh(Analytics.x, Analytics.y, Analytics.vorticity[pst], cmap = 'hot', shading = 'gouraud');
    
    if len(cbar) < 3:
        cbar.append(fig.colorbar(im1,ax=ax[1, 0]))
        cbar.append(fig.colorbar(im2,ax=ax[1, 1]))
        cbar.append(fig.colorbar(im3,ax=ax[1, 2]))
    else:
        cbar[0].mappable.set_clim(vmin = amin(Analytics.potential[pst]), 
                                  vmax = amax(Analytics.potential[pst]))
        cbar[1].mappable.set_clim(vmin = amin(Analytics.ions[pst]),
                                  vmax = amax(Analytics.ions[pst]))
        cbar[2].mappable.set_clim(vmin = amin(Analytics.vorticity[pst]),
                                  vmax = amax(Analytics.vorticity[pst])

    return ax


#=====================================================================================
# For contour with color mesh
#=====================================================================================

def init_mesh_contour(Analytics, ax, fig, model):
    '''
    Initialization of the plots, set the limits in the graphics, labels, titles and scales
    '''
    
    fig.suptitle(model, fontsize=20)
    
    ## Pcolormesh
    ## Potential, density, Vorticity
    ax[0].set_title('Potential over Ions density', fontsize = ftsz_title)
    ax[1].set_title('Vorticity', fontsize = ftsz_title)
    
    return ax

                                  def animate_mesh_contour (position, Analytics, ax):
    '''
    A function to update each frame
    '''
    
    global cbar
    
    ## Pcolormesh
    ## Potentia
    
    ax[0].clear()
    ax[1].clear()
    ax = init_mesh_contour(Analytics, ax, None, None, suptitle = False)

    ax[0].contour(Analytics.x, Analytics.y, Analytics.potential[position], cmap = 'hot')
    
    ## Density
    im1 = ax[0].pcolormesh(Analytics.x, Analytics.y, Analytics.ions[position], cmap = 'hot', shading = 'gouraud')

    ## Vorticity
    im2 = ax[1].pcolormesh(Analytics.x, Analytics.y, Analytics.vorticity[position], cmap = 'hot', shading = 'gouraud')
    if len(cbar) < 2:
        cbar.append(fig.colorbar(im1,ax=ax[0]))
        cbar.append(fig.colorbar(im2,ax=ax[1]))
    else:
        cbar[0].mappable.set_clim(vmin = amin(Analytics.ions[pst]),
                                  vmax = amax(Analytics.ions[pst]))
                                  
        cbar[1].mappable.set_clim(vmin = amin(Analytics.vorticity[pst]),
                                  vmax = amax(Analytics.vorticity[pst]))
                                  
    return ax
