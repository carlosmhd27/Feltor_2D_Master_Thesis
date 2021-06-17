from sys import path
path.insert(1, '/m100/home/userexternal/crodrigu/Plasma/Feltor_2D_Master_Thesis/2D_FELTOR_Analysis/')

from Analysis          import Analyse
from numpy             import amax, amin, absolute, log, gradient
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
amp_int_vrad                       = 0
min_Mass,         max_Mass         = 0, 0
min_int_vort_sqr, max_int_vort_sqr = 0, 0
min_ions,         max_ions         = 0, 0
min_potential,    max_potential    = 0, 0
min_vorticity,    max_vorticity    = 0, 0
min_v_r,          max_v_r          = 0, 0

amp_ions, amp_potential, amp_vorticity = 0, 0, 0

def Analyzed (File_name):
    '''
    A function to get the parameters we want to plot, CM, error in Mass, Velocity of the
    CM, Potential, density and vorticity-
    '''

#    global min_time, min_int_vort, min_Mass, min_int_vort_sqr
#    global max_time, max_int_vort, max_Mass, max_int_vort_sqr
    global min_ions, min_potential, min_vorticity, min_v_r
    global max_ions, max_potential, max_vorticity, max_v_r
    global amp_ions, max_potential, amp_vorticity, amp_int_vrad

    Analytics = Analyse(File_name, input_model = True, dimensions = True, fields = False, get_everything = False, integrate_fields = False)

    Analytics.v_r          = -gradient(Analytics.Data['potential'][:], Analytics.y, axis = 1).transpose()
    Analytics.V_r          = Analytics.integrate(Analytics.v_r) / (Analytics.lx * Analytics.ly)
    Analytics.Mass         = Analitics.integrate(Analytics.Data['ions'][:]) / (Analytics.lx * Analytics.ly)
    Analytics.int_vort_sqr = Analytics.integrate(Analytics.Data['vorticity'][:] ** 2) / (Analytics.lx * Analytics.ly)

    min_ions      = amin(Analytics.Data['ions'][:]);      max_ions      = amax(Analytics.Data['ions'][:])
    min_potential = amin(Analytics.Data['potential'][:]); max_potential = amax(Analytics.Data['potential'][:])
    min_vorticity = amin(Analytics.Data['vorticity'][:]); max_vorticity = amax(Analytics.Data['vorticity'][:])
    min_v_r       = amin(Analytics.v_r);           max_v_r       = amax(Analytics.v_r)

    amp_Mass      = max(absolute(min_ions), absolute(min_ions))
    amp_potential = max(absolute(min_potential), absolute(max_potential))
    amp_vorticity = max(absolute(min_vorticity), absolute(max_vorticity))
    amp_int_vrad  = max(absolute(amax(Analytics.V_r)), absolute(amin(Analytics.V_r)))


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
    ## The colormesh get updated when calling init, but if there is some Probes,
    ## or boundaries, we will make them visible
    ## This happens because we wanna put the colorbar from the begining

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




    return ax


def init(Analytics, ax, fig, model, extra = '', pst = 0, suptitle = True, colormap = 'seismic', colormap_n = 'hot', log_n = True):
                                                                                    #'hot'
    '''
    Initialization of the plots, set the limits in the graphics, labels, titles and scales
    '''

    global cbar
    if min_ions < 0:
        log_n = False

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
    ax[0, 0].set_ylim(- amp_int_vrad - 0.1 * amp_int_vrad, amp_int_vrad + 0.1 * amp_int_vrad)
#    ax[0, 0].set_title('Average Vorticity', fontsize = ftsz_title)
    ax[0, 0].set_title('Average Radial transport', fontsize = ftsz_title) # Vorticity

    ## Total Mass
    ax[0, 1].set_xlabel(r't [$\omega_{ci}^{-1}]$', fontsize = ftsz_label)
    ax[0, 1].set_ylabel(r'n [$N$]', fontsize = ftsz_label)
    ax[0, 1].plot(Analytics.time, Analytics.Mass, color='tab:blue')
#    ax[0, 1].set_xlim(min_time - 0.5, max_time + 0.5)
    # ax[0, 1].set_ylim(- amp_Mass - 0.1 * amp_Mass, amp_Mass + 0.01)
#     ax[0, 1].set_xscale('symlog')
    ax[0, 1].set_title('Average density', fontsize = ftsz_title)

    ## Integrated Vorticity square ** 2
    ax[0, 2].set_xlabel(r't [$\omega_{ci}^{-1}$]', fontsize = ftsz_label)
    ax[0, 2].set_ylabel(r'$\omega^2$ [$\omega_{ci}^{2}$]', fontsize = ftsz_label)
    ax[0, 2].plot(Analytics.time, Analytics.int_vort_sqr, color='tab:blue')
#    ax[0, 2].set_xlim(min_time - 0.5, max_time + 0.5)
    # ax[0, 2].set_ylim(- amp_int_vort_sqr - 0.5, amp_int_vort_sqr + 0.5)
    ax[0, 2].set_title('Average Vorticity Square', fontsize = ftsz_title)

    ## Pcolormesh
    ## Potential, density, Vorticity
    n_title = r'Logarithm of the Density [$\log{N}$]' if log_n else r'Density [$N$]'

    ax[1, 0].set_title(r'Potential [$e / T_e$]', fontsize = ftsz_title)
    ax[1, 1].set_title(n_title, fontsize = ftsz_title)
    ax[1, 2].set_title(r'Vorticity [$\omega_{ci}$]', fontsize = ftsz_title)

    for i in range(len(ax[1])):
        ax[1, i].vlines(Analytics.input['x_a'] * Analytics.lx, Analytics.y[0], Analytics.y[-1], color = 'black')
        ax[1, i].vlines(Analytics.input['x_b'] * Analytics.lx, Analytics.y[0], Analytics.y[-1], color = 'black')

        if Analytics.input['x_c'] < 1:
            ax[1, i].vlines(Analytics.input['x_c'] * Analytics.lx, Analytics.y[0], Analytics.y[-1], color = 'black')

    n_ions = log(Analytics.Data['ions'][pst]) if log_n else Analytics.Data['ions'][pst]

    im1 = ax[1, 0].pcolormesh(Analytics.x, Analytics.y, Analytics.Data['potential'][pst], cmap = colormap, shading = 'gouraud');

    im2 = ax[1, 1].pcolormesh(Analytics.x, Analytics.y, n_ions, cmap = colormap_n, shading = 'gouraud');

    im3 = ax[1, 2].pcolormesh(Analytics.x, Analytics.y, Analytics.Data['vorticity'][pst], cmap = colormap, shading = 'gouraud');

    if len(cbar) < 3:
        cbar.append(fig.colorbar(im1,ax=ax[1, 0]))
        cbar.append(fig.colorbar(im2,ax=ax[1, 1]))
        cbar.append(fig.colorbar(im3,ax=ax[1, 2]))
    else:

        min_n = log(min_ions) if log_n else min_ions
        max_n = log(max_ions) if log_n else max_ions

        cbar[0].mappable.set_clim(vmin = -amp_potential,
                                  vmax =  amp_potential)
        cbar[1].mappable.set_clim(vmin =  min_n,
                                  vmax =  max_n)
        cbar[2].mappable.set_clim(vmin = -amp_vorticity,
                                  vmax =  amp_vorticity)

    return ax

#=====================================================================================
# Plotting Fluxes
#=====================================================================================

def Flux_plot(Analitics, ax, fig, model = None, extra = '',  suptitle = True):
    '''
    Plot the Average Flux, it should be kind of constant in the inner region between
    the source and the LCFS
    '''
    hlf   = Analitics.nt // 2
    Gamma = Analitics.integrate(Analitics.ions[hlf:] * Analitics.v_r[hlf:], dim_integral = 1, axis = -1, typ = 'y') / Analitics.ly
    n_yt  = Analitics.integrate(Analitics.ions[hlf:],                       dim_integral = 1, axis = -1, typ = 'y') / Analitics.ly
    Gamma = Analitics.integrate(Gamma, typ = 't', axis = 0, dim_integral = 1, indep_vars = [Analitics.time[hlf:]]) / (Analitics.time[-1] - Analitics.time[hlf])
    n_yt  = Analitics.integrate(n_yt,  typ = 't', axis = 0, dim_integral = 1, indep_vars = [Analitics.time[hlf:]]) / (Analitics.time[-1] - Analitics.time[hlf])

    if suptitle:
        fig.suptitle(models_name[model] + '\n' + r'lx = {} ly = {} $\nu$ = {} dt = {:1.3f}'.
                     format(Analitics.lx, Analitics.ly,
                     Analitics.input['nu_perp'], Analitics.dt),
                     fontsize = ftsz_title_big)

    ax[0].plot(Analitics.x, Gamma)
    ax[1].plot(Analitics.x, n_yt )

    ax[0].set_xlabel(r'y [$\rho_{s}$]', fontsize = ftsz_label)
    ax[0].set_ylabel(r'$\langle n v_r\rangle_{yt}$ [$n\omega_{ci}\rho_{s}$]', fontsize = ftsz_label) #  $\omega$ [$\omega_{ci}$]
    ax[0].set_title('Average Flux', fontsize = ftsz_title)

    ax[1].set_xlabel(r'y [$\rho_{s}$]', fontsize = ftsz_label)
    ax[1].set_ylabel(r'$\langle n \rangle_{yt}$ [$n$]', fontsize = ftsz_label) #  $\omega$ [$\omega_{ci}$]
    ax[1].set_title('Average Density', fontsize = ftsz_title)
    return fig, ax

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
    ax[1].contour(Analytics.x, Analytics.y, Analytics.potential[position].transpose((0, 2, 1)), cmap = 'hot')

    ## Pcolormesh
    # Radial Velocity
    im0 = ax[0].pcolormesh(Analytics.x, Analytics.y, Analytics.v_r[position].transpose((0, 2, 1)), cmap = 'hot', shading = 'gouraud')

    ## Density
    im1 = ax[1].pcolormesh(Analytics.x, Analytics.y, Analytics.ions[position].transpose((0, 2, 1)), cmap = 'hot', shading = 'gouraud')

    ## Vorticity
    im2 = ax[2].pcolormesh(Analytics.x, Analytics.y, Analytics.vorticity[position].transpose((0, 2, 1)), cmap = 'hot', shading = 'gouraud')
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
