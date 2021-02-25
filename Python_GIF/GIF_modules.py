from sys import path
path.insert(1, '/marconi/home/userexternal/crodrigu/Plasma/prod/2D_FELTOR_Analysis/')

from Analysis          import Analyse
from numpy             import amax, amin
from matplotlib        import use
use("Agg")

ftsz_title_big = 27
ftsz_title     = 23
ftsz_label     = 20
models_name    = {'IC': 'Interchange', 'HW_mod': 'Modified Hasegawa Wakatani', 'HW_ord': 'Ordinary Hasegawa Wakatani', "IC_HW_mod": 'Interchange + Modified HW', "IC_HW_ord": 'Interchange + Ordinary HW', 'HW_mod_IC': 'Modified HW + Interchange', 'HW_ord_IC':'Ordinary HW + Interchange'}
cbar           = []

min_time,         max_time         = 0, 0
min_int_vort,     max_int_vort     = 0, 0
min_int_vrad,     max_int_vrad     = 0, 0
min_Mass,         max_Mass         = 0, 0
min_int_vort_sqr, max_int_vort_sqr = 0, 0
min_ions,         max_ions         = 0, 0
min_potential,    max_potential    = 0, 0
min_vorticity,    max_vorticity    = 0, 0
min_v_r,          max_v_r          = 0, 0
def Analyzed (File_name):
    '''
    A function to get the parameters we want to plot, CM, error in Mass, Velocity of the
    CM, Potential, density and vorticity-
    '''
    
#    global min_time, min_int_vort, min_Mass, min_int_vort_sqr
#    global max_time, max_int_vort, max_Mass, max_int_vort_sqr
    global min_ions, min_potential, min_vorticity, min_v_r
    global max_ions, max_potential, max_vorticity, max_v_r
    
    Analytics = Analyse(File_name)
    
    Analytics.int_vort     = Analytics.integrate('vorticity') / (Analytics.lx * Analytics.ly)
    Analytics.int_vort_sqr = Analytics.integrate(Analytics.vorticity ** 2) / (Analytics.lx * Analytics.ly)
    Analytics.V_r          = Analytics.integrate(Analytics.v_r) / (Analytics.lx * Analytics.ly)

#    min_time         = amin(Analytics.time);         max_time         = amax(Analytics.time)
#    min_int_vrad     = amin(Analytics.V_r);          max_int_vrad     = amax(Analytics.V_r)     
#    min_Mass         = amin(Analytics.Mass);         max_Mass         = amax(Analytics.Mass)
#    min_int_vort_sqr = amin(Analytics.int_vort_sqr); max_int_vort_sqr = amax(Analytics.int_vort_sqr)

    min_ions      = amin(Analytics.ions);      max_ions      = amax(Analytics.ions)
    min_potential = amin(Analytics.potential); max_potential = amax(Analytics.potential)
    min_vorticity = amin(Analytics.vorticity); max_vorticity = amax(Analytics.vorticity)
    min_v_r       = amin(Analytics.v_r);       max_v_r       = amax(Analytics.v_r)
    
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
#    ax[0, 0].set_ylabel(r'$\omega$ [$\omega_{ci}$]', fontsize = ftsz_label)
    ax[0, 0].set_ylabel(r'Average $V_r$ [$\omega_{ci}\rho_s$]', fontsize = ftsz_label) #  $\omega$ [$\omega_{ci}$]
    ax[0, 0].plot(Analytics.time, Analytics.V_r, color='tab:blue');
#    ax[0, 0].set_xlim(min_time - 0.5, max_time + 0.5)
#    ax[0, 0].set_ylim(min_int_vort - 0.5, max_iint_vort + 0.5)
#    ax[0, 0].set_title('Average Vorticity', fontsize = ftsz_title)
    ax[0, 0].set_title('Average Radial transport', fontsize = ftsz_title) # Vorticity
    
    ## Total Mass
    ax[0, 1].set_xlabel(r't [$\omega_{ci}^{-1}]$', fontsize = ftsz_label)
    ax[0, 1].set_ylabel(r'n [$N$]', fontsize = ftsz_label)
    ax[0, 1].plot(Analytics.time, Analytics.Mass, color='tab:blue')
#    ax[0, 1].set_xlim(min_time - 0.5, max_time + 0.5)
#    ax[0, 1].set_ylim(min_Mass, max_Mass + 0.01)
#     ax[0, 1].set_xscale('symlog')
    ax[0, 1].set_title('Average density', fontsize = ftsz_title)
    
    ## Integrated Vorticity square ** 2
    ax[0, 2].set_xlabel(r't [$\omega_{ci}^{-1}$]', fontsize = ftsz_label)
    ax[0, 2].set_ylabel(r'$\omega^2$ [$\omega_{ci}^{2}$]', fontsize = ftsz_label)
    ax[0, 2].plot(Analytics.time, Analytics.int_vort_sqr, color='tab:blue')
#    ax[0, 2].set_xlim(min_time - 0.5, max_time + 0.5)
#    ax[0, 2].set_ylim(min_int_vort_sqr - 0.5, max_int_vort_sqr + 0.5)
    ax[0, 2].set_title('Average Vorticity Square', fontsize = ftsz_title)
    
    ## Pcolormesh
    ## Potential, density, Vorticity
    ax[1, 0].set_title(r'Potential [$e / T_e$]', fontsize = ftsz_title)
    ax[1, 1].set_title(r'Ions density [$N$]', fontsize = ftsz_title)
    ax[1, 2].set_title(r'Vorticity [$\omega_{ci}$]', fontsize = ftsz_title)
    
    im1 = ax[1, 0].pcolormesh(Analytics.x, Analytics.y, Analytics.potential[pst], cmap = 'hot', shading = 'gouraud');
    
    im2 = ax[1, 1].pcolormesh(Analytics.x, Analytics.y, Analytics.ions[pst], cmap = 'hot', shading = 'gouraud');
    
    im3 = ax[1, 2].pcolormesh(Analytics.x, Analytics.y, Analytics.vorticity[pst], cmap = 'hot', shading = 'gouraud');
    
    if len(cbar) < 3:
        cbar.append(fig.colorbar(im1,ax=ax[1, 0]))
        cbar.append(fig.colorbar(im2,ax=ax[1, 1]))
        cbar.append(fig.colorbar(im3,ax=ax[1, 2]))
    else:
        cbar[0].mappable.set_clim(vmin = min_potential, 
                                  vmax = max_potential)
        cbar[1].mappable.set_clim(vmin = min_ions,
                                  vmax = max_ions)
        cbar[2].mappable.set_clim(vmin = min_vorticity,
                                  vmax = max_vorticity)

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
    ## Potential, density, Vorticity
    ax[0].set_title('Radial Velocity', fontsize = ftsz_title)
    ax[1].set_title('Potential over Ions density', fontsize = ftsz_title)
    ax[2].set_title('Vorticity', fontsize = ftsz_title)
    
    return ax

def animate_mesh_contour (position, Analytics, ax, fig):
    '''
    A function to update each frame
    '''
    
    global cbar
    
    ## Pcolormesh
    ## Potentia
    
    # Contour
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
        cbar[0].mappable.set_clim(vmin = min_v_r,
                                  vmax = max_v_r)
        cbar[1].mappable.set_clim(vmin = min_ions,
                                  vmax = max_ions)
        cbar[2].mappable.set_clim(vmin = min_vorticity,
                                  vmax = max_vorticity)

    return ax
