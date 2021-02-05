from sys import path
path.insert(1, '../2D_FELTOR_Analysis/')

from Analysis          import Analyse
from numpy             import amax, amin
from matplotlib        import use
use("Agg")


ftsz_title_big = 27
ftsz_title     = 23
ftsz_label     = 20
models_name    = {'IC': 'Interchange', 'HW_mod': 'Modified Hasegawa Wakatani', 'HW_ord': 'Ordinary Hasegawa Wakatani', "IC_HW_mod": 'Interchange + Modified HW', 
                  "IC_HW_ord": 'Interchange + Ordinary HW', 'HW_mod_IC': 'Modified HW + Interchange', 'HW_ord_IC':'Ordinary HW + Interchange'}

cbar       = []

def Analyzed (File_name):
    '''
    A function to get the parameters we want to plot, CM, error in Mass, Velocity of the
    CM, Potential, density and vorticity-
    '''
    
    Analytics = Analyse(File_name)
    
    Analytics.int_vort     = Analytics.integrate('vorticity')
    Analytics.int_vort_sqr = Analytics.integrate(Analytics.vorticity ** 2)

    return Analytics

def animate (position, Analytics, ax, fig):
    '''
    A function to update each frame
    '''
    
    for i in range(3):
        for j in range(2):
            ax[j, i].clear()
    ax = init(Analytics, ax, fig, None, pst = position, suptitle = False)
    time_pos = Analytics.time[position]

    # Plots
    ## Radial Velocity
    ax[0, 0].plot(time_pos, Analytics.V_r[position], 'o', color='tab:red')
    ax[0, 0].annotate(round(time_pos, 1), (time_pos, Analytics.V_r[position]), fontsize = ftsz_label)

    ## Mass
    ax[0, 1].plot(time_pos, Analytics.Mass[position], 'o',  color='tab:red')
    ax[0, 1].annotate(round(time_pos, 1), (time_pos, Analytics.Mass[position]), fontsize = ftsz_label)

    ## Integrated Vorticity square ** 2
    ax[0, 2].plot(time_pos, Analytics.int_vort_sqr[position], 'o', color='tab:red')
    ax[0, 2].annotate(round(time_pos, 1), (time_pos, Analytics.int_vort_sqr[position]), fontsize = ftsz_label)


    ## The colormesh get updated when calling init
    ## This happens because we wanna put the colorbar from the begining
 
    return ax


def init(Analytics, ax, fig, model, extra = '', pst = 0, suptitle = True):
    '''
    Initialization of the plots, set the limits in the graphics, labels, titles and scales
    '''
    
    global cbar
    
    if suptitle:
        fig.suptitle(models_name[model] + '\n' + r'lx = {} ly = {} $\nu$ = {} dt = {:1.4f}'.
                     format(Analytics.lx, Analytics.ly,
                     Analytics.input['nu_perp'], Analytics.dt),
                     fontsize = ftsz_title_big)
                     
    ## Average Radial Velocity
    ax[0, 0].set_xlabel(r't [$\omega_{ci}^{-1}$]', fontsize = ftsz_label)
    ax[0, 0].set_ylabel(r'Average $V_r$ [$\omega_{ci}\rho_s$]', fontsize = ftsz_label) 
    ax[0, 0].plot(Analytics.time, Analytics.V_r, color='tab:blue')
    ax[0, 0].set_title('Average Radial transport', fontsize = ftsz_title)
    
    ## Total Mass
    ax[0, 1].set_xlabel(r't [$\omega_{ci}^{-1}]$', fontsize = ftsz_label)
    ax[0, 1].set_ylabel(r'n [$N$]', fontsize = ftsz_label)
    ax[0, 1].plot(Analytics.time, Analytics.Mass, color='tab:blue')
    ax[0, 1].set_title('Average density', fontsize = ftsz_title)
    
    ## Integrated Vorticity square ** 2
    ax[0, 2].set_xlabel(r't [$\omega_{ci}^{-1}$]', fontsize = ftsz_label)
    ax[0, 2].set_ylabel(r'$\omega^2$ [$\omega_{ci}^{2}$]', fontsize = ftsz_label)
    ax[0, 2].plot(Analytics.time, Analytics.int_vort_sqr, color='tab:blue')
    ax[0, 2].set_title('Average Vorticity Square', fontsize = ftsz_title)
    
    ## Pcolormesh
    ## Potential, density, Vorticity
    ax[1, 0].set_title(r'Potential [$e / T_e$]', fontsize = ftsz_title)
    ax[1, 1].set_title(r'Ions density [$N$]', fontsize = ftsz_title)
    ax[1, 2].set_title(r'Vorticity [$\omega_{ci}$]', fontsize = ftsz_title)
    
    im1 = ax[1, 0].pcolormesh(Analytics.x, Analytics.y, Analytics.potential[pst], cmap = 'hot', shading = 'gouraud')
    im2 = ax[1, 1].pcolormesh(Analytics.x, Analytics.y, Analytics.ions[pst], cmap = 'hot', shading = 'gouraud')
    im3 = ax[1, 2].pcolormesh(Analytics.x, Analytics.y, Analytics.vorticity[pst], cmap = 'hot', shading = 'gouraud')
    
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
                                  vmax = amax(Analytics.vorticity[pst]))

    return ax


#=====================================================================================
# For contour with color mesh
#=====================================================================================

def init_mesh_contour(Analytics, ax, fig, model = None, extra = '',  suptitle = True):
    '''
    Initialization of the plots, set the limits in the graphics, labels, titles and scales
    '''
    
    if suptitle:
        fig.suptitle(models_name[model] + '\n' + r'lx = {} ly = {} $\nu$ = {} dt = {:1.4f}'.
                     format(Analytics.lx, Analytics.ly,
                     Analytics.input['nu_perp'], Analytics.dt),
                     fontsize = ftsz_title_big)
    
    ## Pcolormesh
    ## Radial Velocity, Potential, density, Vorticity
    ax[0].set_title('Radial Velocity', fontsize = ftsz_title)
    ax[1].set_title('Potential over Ions density', fontsize = ftsz_title)
    ax[2].set_title('Vorticity', fontsize = ftsz_title)
    
    return ax

def animate_mesh_contour (position, Analytics, ax, fig):
    '''
    A function to update each frame
    '''
    
    global cbar
        
    ax[0].clear()
    ax[1].clear()
    ax[2].clear()

    ## Contour
    # Potential
    ax = init_mesh_contour(Analytics, ax, fig, suptitle = False)
    ax[1].contour(Analytics.x, Analytics.y, Analytics.potential[position], cmap = 'hot')
    
    ## Pcolormesh
    # Radial Velocity
    im0 = ax[0].pcolormesh(Analytics.x, Analytics.y, Analytics.v_r[position], cmap = 'hot', shading = 'gouraud')
    
    ## Density
    im1 = ax[1].pcolormesh(Analytics.x, Analytics.y, Analytics.ions[position], cmap = 'hot', shading = 'gouraud')

    ## Vorticity
    im2 = ax[2].pcolormesh(Analytics.x, Analytics.y, Analytics.vorticity[position], cmap = 'hot', shading = 'gouraud')
    if len(cbar) < 2:
        cbar.append(fig.colorbar(im0, ax = ax[0]))
        cbar.append(fig.colorbar(im1, ax = ax[1]))
        cbar.append(fig.colorbar(im2, ax = ax[2]))
    else:
        cbar[0].mappable.set_clim(vmin = amin(Analytics.v_r[position]),
                                  vmax = amax(Analytics.v_r[position]))
        
        cbar[1].mappable.set_clim(vmin = amin(Analytics.ions[position]),
                                  vmax = amax(Analytics.ions[position]))
                                  
        cbar[2].mappable.set_clim(vmin = amin(Analytics.vorticity[position]),
                                  vmax = amax(Analytics.vorticity[position]))
                                  
    return ax
