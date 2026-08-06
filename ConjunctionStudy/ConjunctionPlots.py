#%%

import pyspedas
import numpy as np
import matplotlib.pyplot as plt
import pytplot
import matplotlib.cm as cm
from pytplot import get_data, store_data,timespan
from pyspedas import tplot
from pyspedas import tinterpol
from scipy.stats import binned_statistic_2d
from matplotlib.colors import ListedColormap
import cdflib
import matplotlib.colors as colors

me = 9.1094e-31 #kg
mi = 1837*me
mu0 = 1.2566370e-06  #;m kg / C^2
mu_0 = mu0
eps0 = 8.85e-12   # C^2/Nm^2
e = 1.602e-19 #C
Z = 1 # 1 for H+, 2 for He2+
gamma = 5/3
kb = 1.380649e-23

def stability_condition(betapar, a, b, beta0):
    denom = (betapar - beta0)**b
    Tperppar = 1 + a/denom
    return Tperppar

def brazilplot(z, x_label='', y_label='', z_label='', title='',mask = allbeta_par>0):    
    mirror_params = [0.77, 0.76, -0.016]  # [a, b, beta0]
    firehose_params = [-1.4, 1.0, -0.11]  # [a, b, beta0]
    cyclotron_params = [0.45, 0.4, -0.0004]  # [a, b, beta0]
    parfirehose_params = [-0.47, 0.53, 0.59]  # [a, b, beta0]

    beta_par_range = np.logspace(-3, 3, 100)  # Adjust range and resolution as needed
    mirror_threshold = stability_condition(beta_par_range, *mirror_params)
    firehose_threshold = stability_condition(beta_par_range, *firehose_params)
    cyclotron_threshold = stability_condition(beta_par_range, *cyclotron_params)
    parfirehose_threshold = stability_condition(beta_par_range, *parfirehose_params)

    plt.figure(figsize=(10, 8))
    sc = plt.scatter(allbeta_par[mask], 1/allTparperp[mask], c='k', cmap='plasma', s=0.01)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(title)
    plt.xlim(0.001, 1000)  # Adjust as needed
    plt.ylim(0.1, 10)  # Adjust as needed
    plt.axhline(y=1, color='k', linestyle='-')
    plt.axvline(x=1, color='k', linestyle='-')
    plt.xscale('log')
    plt.yscale('log')
    plt.plot(beta_par_range, mirror_threshold, label='Mirror', color='k',linestyle='dotted', linewidth=2)
    plt.plot(beta_par_range, firehose_threshold, label='Oblique FH', color='k', linestyle='--', linewidth=2)
    plt.plot(beta_par_range, cyclotron_threshold, label='IC', color='k',linestyle='dotted', linewidth=1)
    plt.plot(beta_par_range, parfirehose_threshold, label='Parallel FH', color='k', linestyle='--', linewidth=1)
    plt.legend()
    plt.show()

def brazilhist(z, x_label=r'$\beta_\parallel$', y_label=r'$T_{\perp}/T_{\parallel}$',
     z_label='', title='',nbins=50,mincount=10,vmin=1,vmax=10000, mask = allTparperp>0,
     count=False,scale = 'linear',fitline=False,xrange = [-3,3],yrange = [-1,1],cmap='jet'):    
    mirror_params = [0.77, 0.76, -0.016]  # [a, b, beta0]
    firehose_params = [-1.4, 1.0, -0.11]  # [a, b, beta0]
    cyclotron_params = [0.45, 0.4, -0.0004]  # [a, b, beta0]
    parfirehose_params = [-0.47, 0.53, 0.59]  # [a, b, beta0]

    beta_par_range = np.logspace(xrange[0], xrange[1], 100)  # Adjust range and resolution as needed
    mirror_threshold = stability_condition(beta_par_range, *mirror_params)
    firehose_threshold = stability_condition(beta_par_range, *firehose_params)
    cyclotron_threshold = stability_condition(beta_par_range, *cyclotron_params)
    parfirehose_threshold = stability_condition(beta_par_range, *parfirehose_params)

    # plt.figure(figsize=(10, 8))
    # sc = plt.scatter(allbeta_par[mask], 1/allTparperp[mask], c=z[mask], cmap='plasma', s=0.01)
    x = allbeta_par[mask]
    y = 1/allTparperp[mask]
    z = z[mask]
    valmask = np.isfinite(x) & np.isfinite(y) & np.isfinite(z) 
    x_valid = x[valmask]
    y_valid = y[valmask]
    z_valid = z[valmask]

    x_fit = np.log10(x_valid)
    y_fit = np.log10(y_valid)
    m,c = np.polyfit(x_fit, y_fit, 1)


    nbins = nbins
    mincount=mincount
    # xbins = np.linspace(0.001, 100, nbins)
    # ybins = np.linspace(0.1, 10, nbins)
    xbins = np.logspace(xrange[0], xrange[1], nbins)
    ybins = np.logspace(yrange[0], yrange[1], nbins)

    xline = np.logspace(-3, 3, 200)
    yline = 10**(m*np.log10(xline) + c)

    if count==True:
        hist, xedges, yedges, _= binned_statistic_2d(x_valid, y_valid, z_valid, statistic='count', bins=[xbins, ybins])
    else:
        hist, xedges, yedges, _= binned_statistic_2d(x_valid, y_valid, z_valid, statistic='mean', bins=[xbins, ybins])

    counts, _, _, _= binned_statistic_2d(x_valid, y_valid, z_valid, statistic='count', bins=[xbins, ybins])
    zero_mask = counts <= mincount
    hist_masked = np.ma.masked_where(counts <= mincount, hist)

    fig, ax = plt.subplots(figsize=(10, 8))

    # if np.any(zero_mask):
        # grey_data = np.where(zero_mask, 1.0, np.nan)
        # ax.pcolormesh(xedges, yedges, grey_data.T, cmap='gray', vmin=0, vmax=1)
    if scale == 'log':
        norm = colors.LogNorm(vmin=vmin, vmax=vmax)
    else:
        norm = colors.Normalize(vmin=vmin, vmax=vmax)

    im = ax.pcolormesh(xedges, yedges, hist_masked.T, cmap=cmap, norm=norm)

    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label(z_label,fontdict={'size': 15})
    if count==True:
        cbar.set_label('Count',fontdict={'size': 15})



    if fitline==True:
        ax.plot(xline, yline, 'r-', linewidth=3,label=rf'Best fit: $\beta_\parallel^{{{m:.2f}}}$',color='b',linestyle='-')
    else:
        pass
  

    # plt.colorbar(im, label=z_label)
    ax.set_facecolor('lightgray')
    ax.set_xlabel(x_label, fontdict={'size': 15})
    ax.set_ylabel(y_label, fontdict={'size': 15})
    ax.set_title(title, fontdict={'size': 15})
    ax.set_xlim(10**(xrange[0]), 10**(xrange[1]))  # Adjust as needed
    ax.set_ylim(10**(-1), 10**(1))  # Adjust as needed
    ax.axhline(y=1, color='k', linestyle='-')
    ax.axvline(x=1, color='k', linestyle='-')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.plot(beta_par_range, mirror_threshold, label='Mirror', color='k',linestyle='dotted', linewidth=2)
    ax.plot(beta_par_range, firehose_threshold, label='Oblique FH', color='k', linestyle='--', linewidth=2)
    ax.plot(beta_par_range, cyclotron_threshold, label='IC', color='k',linestyle='dotted', linewidth=1)
    ax.plot(beta_par_range, parfirehose_threshold, label='Parallel FH', color='k', linestyle='--', linewidth=1)
    ax.legend()
    plt.show()
    return fig, ax


# def brazil_multi(x_label=r'$\beta_\parallel$',count = False, y_label=r'$T_{\perp}/T_{\parallel}$',
#      title='',nbins=50,mincount=10,vmin=1,vmax=10000, mask = allTparperp>0,
#      scale = 'linear',fitline=False,xrange = [-3,3],yrange = [-1,1],cmap='jet'):    
#     mirror_params = [0.77, 0.76, -0.016]  # [a, b, beta0]
#     firehose_params = [-1.4, 1.0, -0.11]  # [a, b, beta0]
#     cyclotron_params = [0.45, 0.4, -0.0004]  # [a, b, beta0]
#     parfirehose_params = [-0.47, 0.53, 0.59]  # [a, b, beta0]

#     beta_par_range = np.logspace(xrange[0], xrange[1], 100)
#     mirror_threshold = stability_condition(beta_par_range, *mirror_params)
#     firehose_threshold = stability_condition(beta_par_range, *firehose_params)
#     cyclotron_threshold = stability_condition(beta_par_range, *cyclotron_params)
#     parfirehose_threshold = stability_condition(beta_par_range, *parfirehose_params)

#     x = allbeta_par[mask]
#     y = 1/allTparperp[mask]
#     valmask = np.isfinite(x) & np.isfinite(y)
#     x_valid = x[valmask]
#     y_valid = y[valmask]

#     x_fit = np.log10(x_valid)
#     y_fit = np.log10(y_valid)
#     m,c = np.polyfit(x_fit, y_fit, 1)

#     xbins = np.logspace(xrange[0], xrange[1], nbins)
#     ybins = np.logspace(yrange[0], yrange[1], nbins)

#     xline = np.logspace(-3, 3, 200)
#     yline = 10**(m*np.log10(xline) + c)

#     counts, xedges, yedges, _ = binned_statistic_2d(x_valid, y_valid, None, statistic='count', bins=[xbins, ybins])
#     hist_masked = np.ma.masked_where(counts <= mincount, counts)

#     fig, ax = plt.subplots(figsize=(10, 8))

#     if scale == 'log':
#         norm = colors.LogNorm(vmin=vmin, vmax=vmax)
#     else:
#         norm = colors.Normalize(vmin=vmin, vmax=vmax)



#     if count==True:
#         im = ax.pcolormesh(xedges, yedges, hist_masked.T, cmap=cmap, norm=norm)
#         cbar = plt.colorbar(im, ax=ax)
#         cbar.set_label('Count', fontdict={'size': 15})
#     else:
#         presence = np.ma.masked_where(hist_masked.mask, np.ones_like(hist_masked))
#         # im = ax.pcolormesh(xedges, yedges, presence.T, cmap=colors.ListedColormap(['g']))

#         xcenters = 0.5 * (xedges[:-1] + xedges[1:])
#         ycenters = 0.5 * (yedges[:-1] + yedges[1:])
#         contour_levels = np.logspace(np.log10(mincount), np.log10(counts.max()), 5)
#         cs = ax.contour(xcenters, ycenters, counts.T, levels=contour_levels,
#                         colors='r', linewidths=1)
#         ax.clabel(cs, inline=True, fontsize=8)
        
    

#     ax.set_facecolor('lightgray')
#     ax.set_xlabel(x_label, fontdict={'size': 15})
#     ax.set_ylabel(y_label, fontdict={'size': 15})
#     ax.set_title(title, fontdict={'size': 15})
#     ax.set_xlim(10**(xrange[0]), 10**(xrange[1]))
#     ax.set_ylim(10**(-1), 10**(1))
#     ax.axhline(y=1, color='k', linestyle='-')
#     ax.axvline(x=1, color='k', linestyle='-')
#     ax.set_xscale('log')
#     ax.set_yscale('log')
#     ax.plot(beta_par_range, mirror_threshold, label='Mirror', color='k',linestyle='dotted', linewidth=2)
#     ax.plot(beta_par_range, firehose_threshold, label='Oblique FH', color='k', linestyle='--', linewidth=2)
#     ax.plot(beta_par_range, cyclotron_threshold, label='IC', color='k',linestyle='dotted', linewidth=1)
#     ax.plot(beta_par_range, parfirehose_threshold, label='Parallel FH', color='k', linestyle='--', linewidth=1)
#     ax.legend()
#     plt.show()
#     return fig, ax





def brazil_multi(betapar_psp=None, Tparperp_psp=None,
                  betapar_solo=None, Tparperp_solo=None,
                  betapar_wind=None, Tparperp_wind=None,
                  mask_psp=None, mask_solo=None, mask_wind=None,
                  colors_dict={'psp': 'red', 'solo': 'green', 'wind': 'blueviolet'},
                  x_label=r'$\beta_\parallel$', y_label=r'$T_{\perp}/T_{\parallel}$',
                  title='', nbins=50, mincount=10,
                  fitline=False, xrange=[-3, 3], yrange=[-1, 1],
                  contour_nlevels=5, norm=False,
                  plot_type='contour', scatter_size=1, scatter_alpha=0.3,
                  scatter_max_points=None, mean_marker_size=120):

    mirror_params = [0.77, 0.76, -0.016]  # [a, b, beta0]
    firehose_params = [-1.4, 1.0, -0.11]  # [a, b, beta0]
    cyclotron_params = [0.45, 0.4, -0.0004]  # [a, b, beta0]
    parfirehose_params = [-0.47, 0.53, 0.59]  # [a, b, beta0]

    beta_par_range = np.logspace(xrange[0], xrange[1], 100)
    mirror_threshold = stability_condition(beta_par_range, *mirror_params)
    firehose_threshold = stability_condition(beta_par_range, *firehose_params)
    cyclotron_threshold = stability_condition(beta_par_range, *cyclotron_params)
    parfirehose_threshold = stability_condition(beta_par_range, *parfirehose_params)

    xbins = np.logspace(xrange[0], xrange[1], nbins)
    ybins = np.logspace(yrange[0], yrange[1], nbins)
    xcenters = 0.5 * (xbins[:-1] + xbins[1:])
    ycenters = 0.5 * (ybins[:-1] + ybins[1:])

    spacecraft_data = {
        'psp': (betapar_psp, Tparperp_psp, mask_psp),
        'solo': (betapar_solo, Tparperp_solo, mask_solo),
        'wind': (betapar_wind, Tparperp_wind, mask_wind),
    }

    fig, ax = plt.subplots(figsize=(10, 8))

    for sc, (betapar, Tparperp, mask) in spacecraft_data.items():
        if betapar is None or Tparperp is None:
            continue  # skip spacecraft the user didn't provide

        if mask is None:
            mask = Tparperp > 0

        x = betapar[mask]
        y = 1 / Tparperp[mask]
        valmask = np.isfinite(x) & np.isfinite(y)
        x_valid = x[valmask]
        y_valid = y[valmask]

        color = colors_dict.get(sc, 'k')

        if plot_type == 'scatter':
            x_plot, y_plot = x_valid, y_valid
            if scatter_max_points is not None and x_plot.size > scatter_max_points:
                idx = np.random.choice(x_plot.size, scatter_max_points, replace=False)
                x_plot, y_plot = x_plot[idx], y_plot[idx]

            ax.scatter(x_plot, y_plot, s=scatter_size, alpha=scatter_alpha,
                       color=color, edgecolors='none')

            # mark the mean position of this spacecraft's points
            x_mean = np.mean(x_valid)
            y_mean = np.mean(y_valid)
            ax.scatter(x_mean, y_mean, marker='^', color=color, s=2*mean_marker_size,
                       edgecolors='k', linewidths=1.5, zorder=10, label=sc.upper())

        elif plot_type == 'means_only':
            # just the mean marker, no contour or scatter cloud
            x_mean = np.mean(x_valid)
            y_mean = np.mean(y_valid)
            ax.scatter(x_mean, y_mean, marker='^', color=color, s=4*mean_marker_size,
                       edgecolors='k', linewidths=1.5, zorder=10, label=sc.upper())

        else:  # plot_type == 'contour'
            counts, xedges, yedges, _ = binned_statistic_2d(
                x_valid, y_valid, None, statistic='count', bins=[xbins, ybins])

            if counts.max() <= mincount:
                continue  # nothing to contour for this spacecraft

            if norm:
                total = x_valid.size  # total valid points for this spacecraft
                plot_counts = counts / total
                plot_mincount = mincount / total
                plot_max = counts.max() / total
            else:
                plot_counts = counts
                plot_mincount = mincount
                plot_max = counts.max()

            contour_levels = np.logspace(np.log10(plot_mincount), np.log10(plot_max), contour_nlevels)
            cs = ax.contour(xcenters, ycenters, plot_counts.T, levels=contour_levels,
                             colors=color, linewidths=0.6, alpha=0.4)

            # mark the counts-weighted mean bin location (centroid of the binned distribution)
            xgrid, ygrid = np.meshgrid(xcenters, ycenters, indexing='ij')
            total_counts = counts.sum()
            if total_counts > 0:
                x_mean = np.sum(xgrid * counts) / total_counts
                y_mean = np.sum(ygrid * counts) / total_counts
                ax.scatter(x_mean, y_mean, marker='^', color=color, s=2*mean_marker_size,
                           edgecolors=color, linewidths=1.5, zorder=10, label=sc.upper())

    ax.set_facecolor('lightgray')
    ax.set_xlabel(x_label, fontdict={'size': 15})
    ax.set_ylabel(y_label, fontdict={'size': 15})
    ax.set_title(title, fontdict={'size': 15})
    ax.set_xlim(10**(xrange[0]), 10**(xrange[1]))
    ax.set_ylim(10**(-1), 10**(1))
    ax.axhline(y=1, color='k', linestyle='-')
    ax.axvline(x=1, color='k', linestyle='-')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.plot(beta_par_range, mirror_threshold, label='Mirror', color='k', linestyle='dotted', linewidth=2)
    ax.plot(beta_par_range, firehose_threshold, label='Oblique FH', color='k', linestyle='--', linewidth=2)
    ax.plot(beta_par_range, cyclotron_threshold, label='IC', color='k', linestyle='dotted', linewidth=1)
    ax.plot(beta_par_range, parfirehose_threshold, label='Parallel FH', color='k', linestyle='--', linewidth=1)
    ax.legend(fontsize=15)
    plt.show()
    return fig, ax
# brazilplot(1e3*allbeta_par, r'$\beta_\parallel$', r'$T_{\perp}/T_{\parallel}$', r'$v_x$'

# nratio = abs(alln)/np.nanmean(abs(alln))
# Tratio = abs(allT)/np.nanmean(abs(allT))
# vratio = abs(allvmags)/np.nanmean(abs(allvmags))
# nTratio = nratio*Tratio

# bulkratio = abs(alln)*allvmags**2/np.nanmean(abs(alln)*allvmags**2)

# brazilhist(1e-3*abs(allvmags), z_label = r'$v \ (km/s)$',fitline=False,nbins=80,vmin=200,vmax=300,count=False,scale='linear')

# brazilhist(abs(allT)/np.nanmean(abs(allT)), z_label = r'$T/\langle T \rangle$',fitline=False,nbins=80,vmin=0,vmax=200,count=True,scale='linear')
# brazilhist(abs(allvmags)/np.nanmean(abs(allvmags)), z_label = r'$v/\langle v \rangle$',fitline=False,nbins=80,vmin=0,vmax=200,count=False,scale='linear')
# brazilhist(nratio*Tratio, z_label = r'$nT/\langle nT \rangle$',fitline=False,nbins=80,vmin=0,vmax=200,count=True,scale='linear')
# brazilhist(bulkratio, z_label = r'$nv^2/\langle nv^2 \rangle$',fitline=False,nbins=80,vmin=0,vmax=200,count=True,scale='linear')



# brazilhist(1.602e19*abs(allT), z_label = r'$T \ (eV)$',fitline=False,nbins=60,vmin=00,vmax=200,count=False,scale='linear')
# brazilhist(1.602e19*abs(allT), z_label = r'$T \ (eV)$',fitline=False,nbins=100,vmin=0,vmax=250,count=False,scale='linear',cmap='plasma')
# brazilhist(1e-6*abs(alln), z_label = r'$n \ (cm^{-3})$',fitline=False,nbins=100,vmin=1e3,vmax=4e3,count=False,scale='linear',cmap='plasma')
# brazilhist(np.log10(allmachs), z_label = r'$Log(Ma) \ (cm^{-3})$',fitline=False,nbins=100,vmin=-0.2,vmax=0.2,count=False,scale='linear',cmap='seismic')

# brazilhist(allT, z_label = r'$T \ (eV)$',fitline=False,nbins=100,vmin=0,vmax=1000,count=True,scale='linear',cmap='plasma')
# brazil_multi(fitline=False,nbins=100,vmin=0,vmax=1000,count=False,scale='linear',cmap='plasma')
# brazilhist(1e-6*abs(alln), z_label = r'$n \ (cm^{-3})$',fitline=False,nbins=100,vmin=1e3,vmax=4e3,count=False,scale='linear',cmap='plasma')
# brazilhist(1e-6*abs(alln), z_label = r'$n \ (cm^{-3})$',fitline=False,nbins=60,vmin=0,vmax=10,count=False,scale='linear')


#Here's synthetic test data for all three spacecraft, positioned in different regions of beta_par space with different characteristic anisotropies (mimicking the real physics — PSP near the Sun sees lower beta and larger anisotropies, WIND at 1 AU sees higher beta and more isotropic plasma):

# python
np.random.seed(42)

def make_test_data(n, beta_center_log10, beta_spread, aniso_center_log10, aniso_spread):
    """Generate synthetic betapar and Tparperp arrays clustered around a
    given region of (beta_par, Tperp/Tpar) space."""
    log_beta = np.random.normal(beta_center_log10, beta_spread, n)
    betapar = 10**log_beta

    # y-axis in the plot is Tperp/Tpar = 1/Tparperp, so build in that space
    log_aniso = np.random.normal(aniso_center_log10, aniso_spread, n)
    Tperp_over_Tpar = 10**log_aniso
    Tparperp = 1 / Tperp_over_Tpar  # since brazil_multi computes y = 1/Tparperp

    return betapar, Tparperp

# PSP: near the Sun -> lower beta_par, larger anisotropies (both mirror- and
# firehose-leaning), clustered around beta_par ~ 0.3
betapar_psp, Tparperp_psp = make_test_data(
    n=20000, beta_center_log10=-0.5, beta_spread=0.25,
    aniso_center_log10=0.25, aniso_spread=0.15)

# Solo: intermediate distance -> beta_par ~ 1, moderate anisotropy
betapar_solo, Tparperp_solo = make_test_data(
    n=20000, beta_center_log10=0.0, beta_spread=0.25,
    aniso_center_log10=0.0, aniso_spread=0.12)

# WIND: 1 AU -> higher beta_par, closer to isotropic (Tperp/Tpar near 1)
betapar_wind, Tparperp_wind = make_test_data(
    n=20000, beta_center_log10=0.6, beta_spread=0.3,
    aniso_center_log10=-0.15, aniso_spread=0.15)

# masks (all positive Tparperp already, but matching the function's default pattern)
mask_psp = Tparperp_psp > 0
mask_solo = Tparperp_solo > 0
mask_wind = Tparperp_wind > 0







#%%

fig, ax = brazil_multi(betapar_psp=betaparspsp_preceding, Tparperp_psp=Tparperppsp_preceding,
                        betapar_solo=betaparssolo_preceding, Tparperp_solo=Tparperpsolo_preceding,
                        betapar_wind=betaparswind_preceding, Tparperp_wind=Tparperpwind_preceding,
                        contour_nlevels=10,title = "Preceding Wind",plot_type='contour',
                        mean_marker_size=35)

#%%
fig, ax = brazil_multi(betapar_psp=betaparspsp_trailing, Tparperp_psp=Tparperppsp_trailing,
                        betapar_solo=betaparssolo_trailing, Tparperp_solo=Tparperpsolo_trailing,
                        betapar_wind=betaparswind_trailing, Tparperp_wind=Tparperpwind_trailing,
                        contour_nlevels=10,title = "Trailing Wind",plot_type='contour',
                        mean_marker_size=35)
#%%
fig, ax = brazil_multi(betapar_psp=betaparspsp_shoulder, Tparperp_psp=Tparperppsp_shoulder,
                        betapar_solo=betaparssolo_shoulder, Tparperp_solo=Tparperpsolo_shoulder,
                        betapar_wind=betaparswind_shoulder, Tparperp_wind=Tparperpwind_shoulder,
                        contour_nlevels=10,title = "Shoulder Wind",plot_type='contour',
                        mean_marker_size=35)

#%%
fig, ax = brazil_multi(betapar_psp=betaparspsp_fast, Tparperp_psp=Tparperppsp_fast,
                        betapar_solo=betaparssolo_fast, Tparperp_solo=Tparperpsolo_fast,
                        betapar_wind=betaparswind_fast, Tparperp_wind=Tparperpwind_fast,
                        contour_nlevels=10,title = "Fast Wind",plot_type='contour',
                        mean_marker_size=35)






# #%%
# fig, ax = brazil_multi(betapar_psp=betaparspsp_slow, Tparperp_psp=Tparperppsp_slow,
#                         betapar_solo=betaparssolo_slow, Tparperp_solo=Tparperpsolo_slow,
#                         betapar_wind=betaparswind_slow, Tparperp_wind=Tparperpwind_slow,
#                         contour_nlevels=10,title = "Slow Wind",plot_type='contour',
#                         mean_marker_size=35)

# #%%

# fig, ax = brazil_multi(betapar_psp=betaparspsp, Tparperp_psp=Tparperppsp,
#                         betapar_solo=betaparssolo, Tparperp_solo=Tparperpsolo,
#                         betapar_wind=betaparswind, Tparperp_wind=Tparperpwind,
#                         contour_nlevels=10,title = "Full Intervals",plot_type='contour',
#                         mean_marker_size=35)
#%%

# brazilhist(abs(), z_label = r'$n_iT_i \ (J/cm^3)$',fitline=False,nbins=80,vmin=0,vmax=1e-15,count=False,scale='linear')
# brazilhist(1e-6*0.5*abs(alln)*mi*allvmags**2, z_label = r'$\frac{1}{2}nm_i |v_i|^2 \ (J/cm^3)$',fitline=False,nbins=80,vmin=0,vmax=5e-14,count=False,scale='linear')
# brazilhist((1e-6*0.5*allBmags**2)/mu0, z_label = r'$\frac{|B|^2}{2\mu_0} \ (J/cm^3)$',fitline=False,nbins=80,vmin=0,vmax=1e-15,count=False,scale='linear')
# Ratio = (1e-6*abs(allT)*abs(alln)) / (1e-6*0.5*abs(alln)*mi*allvmags**2)
# brazilhist(Ratio, z_label = r'$\frac{n_iT_i}{\frac{1}{2}nm_i |v_i|^2}$',fitline=False,nbins=80,vmin=0,vmax=0.03,count=False,scale='linear')
# # brazilhist((1e-6*0.5*allBmags**2)/mu0, z_label = r'$\frac{1}{2}n|B|^2 \ (J/cm^3)$',fitline=False,nbins=80,vmin=0,vmax=1e-15,count=False,scale='linear')



x = np.array([1,2,3])

def get_vars(type = 'fast'):
    rarray = np.array([0.05/0.45,0.3/0.45,1/0.45])

    if type == 'fast':
        narray = 1e6*np.array([npsp_fast,nsolo_fast,nwind_fast])
        n_stdarray = 1e6*np.array([nstdpsp_fast,nstdsolo_fast,nstdwind_fast])
        Tarray = 1.602e-19*np.array([Tpsp_fast,Tsolo_fast,Twind_fast])
        T_stdarray = 1.602e-19*np.array([Tstdpsp_fast,Tstdsolo_fast,Tstdwind_fast])
        varray = 1e3*np.array([vpsp_fast,vsolo_fast,vwind_fast])
        v_stdarray = 1e3*np.array([vstdpsp_fast,vstdsolo_fast,vstdwind_fast])
        Barray = 1e-9*np.array([Bpsp_fast,Bsolo_fast,Bwind_fast])
        B_stdarray = 1e-9*np.array([Bstdpsp_fast,Bstdsolo_fast,Bstdwind_fast])
        dB_array = 1e-9*np.array([dBpsp_fast,dBsolo_fast,dBwind_fast])
        dB_stdarray = 1e-9*np.array([dBstdpsp_fast,dBstdsolo_fast,dBstdwind_fast])
    elif type == 'slow':
        narray = 1e6*np.array([npsp_slow,nsolo_slow,nwind_slow])
        n_stdarray = 1e6*np.array([nstdpsp_slow,nstdsolo_slow,nstdwind_slow])
        Tarray = 1.602e-19*np.array([Tpsp_slow,Tsolo_slow,Twind_slow])
        T_stdarray = 1.602e-19*np.array([Tstdpsp_slow,Tstdsolo_slow,Tstdwind_slow])
        varray = 1e3*np.array([vpsp_slow,vsolo_slow,vwind_slow])
        v_stdarray = 1e3*np.array([vstdpsp_slow,vstdsolo_slow,vstdwind_slow])
        Barray = 1e-9*np.array([Bpsp_slow,Bsolo_slow,Bwind_slow])
        B_stdarray = 1e-9*np.array([Bstdpsp_slow,Bstdsolo_slow,Bstdwind_slow])
        dB_array = 1e-9*np.array([dBpsp_slow,dBsolo_slow,dBwind_slow])
        dB_stdarray = 1e-9*np.array([dBstdpsp_slow,dBstdsolo_slow,dBstdwind_slow])
    elif type == 'ir':
        narray = 1e6*np.array([npsp,nsolo,nwind])
        n_stdarray = 1e6*np.array([nstdpsp,nstdsolo,nstdwind])
        Tarray = 1.602e-19*np.array([Tpsp,Tsolo,Twind])
        T_stdarray = 1.602e-19*np.array([Tstdpsp,Tstdsolo,Tstdwind])
        varray = 1e3*np.array([vpsp,vsolo,vwind])
        v_stdarray = 1e3*np.array([vstdpsp,vstdsolo,vstdwind])
        Barray = 1e-9*np.array([Bpsp,Bsolo,Bwind])
        B_stdarray = 1e-9*np.array([Bstdpsp,Bstdsolo,Bstdwind])
        dB_array = 1e-9*np.array([dBpsp,dBsolo,dBwind])
        dB_stdarray = 1e-9*np.array([dBstdpsp,dBstdsolo,dBstdwind])
    else:
        print('No var type selected')
    
  
    narray_norm = (rarray**2)*narray
    n_stdarray_norm = (rarray**2)*n_stdarray
    Tarray_norm = (rarray**(4/3))*Tarray
    T_stdarray_norm = (rarray**(4/3))*T_stdarray
    varray_norm = varray
    v_stdarray_norm = v_stdarray
    Barray_norm = (rarray**2)*Barray
    B_stdarray_norm = (rarray**2)*B_stdarray
    dB_array_norm = (rarray**2)*dB_array
    dB_stdarray_norm = (rarray**2)*dB_stdarray

    return narray_norm,Tarray_norm,varray_norm,Barray_norm,dB_array_norm, n_stdarray_norm, T_stdarray_norm, v_stdarray_norm, B_stdarray_norm, dB_stdarray_norm



nfast,Tfast,vfast,Bfast,dBfast,nfast_std,Tfast_std,vfast_std,Bfast_std,dBfast_std = get_vars('fast')
nslow,Tslow,vslow,Bslow,dBslow,nslow_std,Tslow_std,vslow_std,Bslow_std,dBslow_std = get_vars('slow')
nir,Tir,vir,Bir,dBir,nir_std,Tir_std,vir_std,Bir_std,dBir_std = get_vars('ir')

n = nfast + nslow
T = (Tfast + Tslow)/2
bulkfast = 0.5*nfast*mi*vfast**2
bulkslow = 0.5*nslow*mi*vslow**2
bulkir = 0.5*nir*mi*vir**2
thermfast = (3/2)*kb*nfast*Tfast
thermslow = (3/2)*kb*nslow*Tslow
thermir = (3/2)*kb*nir*Tir
magfast = (0.5*Bfast**2)/(mu0)
magslow = (0.5*Bslow**2)/(mu0)
magir = (0.5*Bir**2)/(mu0)
magmean = (magfast + magslow)/2
deltamagfast = (0.5*dBfast**2)/(mu0)
deltamagslow = (0.5*dBslow**2)/(mu0)
deltamagir = (0.5*dBir**2)/(mu0)
plasmafast = thermfast + bulkfast
plasmaslow = thermslow + bulkslow
plasmair = thermir + bulkir

plasma = plasmafast + plasmaslow
fieldfast = magfast + deltamagfast
fieldslow = magslow + deltamagslow
fieldir = magir + deltamagir
field = fieldfast + fieldslow
# %%

# Number Density
plt.errorbar(x,nfast,yerr=nfast_std,fmt='o-',color = 'r',label = r'$n^{\prime}_{fast}$',capsize=5)
plt.errorbar(x,nslow,yerr=nslow_std,fmt='o-',color='b',label = r'$n^{\prime}_{slow}$',capsize=5)
plt.errorbar(x,nir,yerr=nir_std,fmt='o-',color='purple',label = r'$n^{\prime}_{ir}$',capsize=5)
plt.plot(x,(nslow+nfast)/2,color='g',marker='s',linestyle='dashed',label = r'$n^{\prime}_{avg}$')
plt.title('Adjusted Number Density',fontsize=15)
plt.axvline(x=1,linestyle='dashed',color='k',linewidth=0.2)
plt.axvline(x=2,linestyle='dashed',color='k',linewidth=0.2)
plt.axvline(x=3,linestyle='dashed',color='k',linewidth=0.2)
plt.xticks(x,['PSP \n $\sim0.05 \  AU$','SOLO \n $\sim0.3 \  AU$','WIND \n $\sim1 \ AU$	'],fontsize=12)
plt.legend(fontsize=15)

#%%
# Temperature
plt.errorbar(x,Tfast,yerr=Tfast_std,fmt='o-',color = 'r',label = r'$T^{\prime}_{fast}$',capsize=5)
plt.errorbar(x,Tslow,yerr=Tslow_std,fmt='o-',color='b',label = r'$T^{\prime}_{slow}$',capsize=5)
plt.errorbar(x,Tir,yerr=Tir_std,fmt='o-',color='purple',label = r'$T^{\prime}_{ir}$',capsize=5)
plt.plot(x,(Tslow+Tfast)/2,color='g',marker='s',linestyle='dashed',label = r'$T^{\prime}_{avg}$')
plt.title('Adjusted Temperature',fontsize=15)
plt.axvline(x=1,linestyle='dashed',color='k',linewidth=0.2)
plt.axvline(x=2,linestyle='dashed',color='k',linewidth=0.2)
plt.axvline(x=3,linestyle='dashed',color='k',linewidth=0.2)
plt.xticks(x,['PSP \n $\sim0.05 \  AU$','SOLO \n $\sim0.3 \  AU$','WIND \n $\sim1 \ AU$	'],fontsize=12)
plt.legend(fontsize=15)



#%%
# Velocity
plt.errorbar(x,1e-3*vfast,yerr=1e-3*vfast_std,fmt='o-',color = 'r',label = r'$v_{fast}$',capsize=5)
plt.errorbar(x,1e-3*vslow,yerr=1e-3*vslow_std,fmt='o-',color='b',label = r'$v_{slow}$',capsize=5)
plt.errorbar(x,1e-3*vir,yerr=1e-3*vir_std,fmt='o-',color='purple',label = r'$v_{ir}$',capsize=5)
plt.errorbar(x,1e-3*(vslow+vfast)/2,color = 'g',marker='s',linestyle='dashed',label = r'$v_{avg}$')

# plt.ylim(200,800)
plt.title('Velocity [km/s]',fontsize=15)
plt.axvline(x=1,linestyle='dashed',color='k',linewidth=0.2)
plt.axvline(x=2,linestyle='dashed',color='k',linewidth=0.2)
plt.axvline(x=3,linestyle='dashed',color='k',linewidth=0.2)
plt.ylim(200,800)
plt.yticks(fontsize=12)
plt.xticks(x,['PSP \n $\sim0.05 \  AU$','SOLO \n $\sim0.3 \  AU$','WIND \n $\sim1 \ AU$	'],fontsize=12)
plt.legend(fontsize=15)

#%%
# Magnetic Field Strength
plt.errorbar(x,Bfast,yerr=Bfast_std,fmt='o-',color = 'r',label = r'$|B|_{fast}$',capsize=5)
plt.errorbar(x,Bslow,yerr=Bslow_std,fmt='o-',color='b',label = r'$|B|_{slow}$',capsize=5)
plt.errorbar(x,Bir,yerr=Bir_std,fmt='o-',color='purple',label = r'$|B|_{ir}$',capsize=5)
plt.errorbar(x,(Bslow+Bfast)/2,color = 'g',marker='s',linestyle='dashed',label = r'$|B|_{avg}$')
plt.title('Adjusted Magnetic Field Strength [T]')
plt.ylim(0,2.5e-9)
plt.axvline(x=1,linestyle='dashed',color='k',linewidth=0.2)
plt.axvline(x=2,linestyle='dashed',color='k',linewidth=0.2)
plt.axvline(x=3,linestyle='dashed',color='k',linewidth=0.2)
plt.yticks(fontsize=12)
plt.xticks(x,['PSP \n $\sim0.05 \  AU$','SOLO \n $\sim0.3 \  AU$','WIND \n $\sim1 \ AU$	'],fontsize=12)
plt.legend(fontsize=15)

#%%
# Magnetic Fluctuation Strength
plt.errorbar(x,dBfast,yerr=dBfast_std,fmt='o-',color = 'r',label = r'$\delta B_{fast}$',capsize=5)
plt.errorbar(x,dBslow,yerr=dBslow_std,fmt='o-',color='b',label = r'$\delta B_{slow}$',capsize=5)
plt.errorbar(x,dBir,yerr=dBir_std,fmt='o-',color='purple',label = r'$\delta B_{ir}$',capsize=5)
plt.errorbar(x,(dBslow+dBfast)/2,color = 'g',marker='s',linestyle='dashed',label = r'$\delta B_{avg}$')
plt.ylim(0,2.5e-9)
plt.title('Adjusted Magnetic Fluctuation Magnitude [T]')
plt.axvline(x=1,linestyle='dashed',color='k',linewidth=0.2)
plt.axvline(x=2,linestyle='dashed',color='k',linewidth=0.2)
plt.axvline(x=3,linestyle='dashed',color='k',linewidth=0.2)
plt.yticks(fontsize=12)
plt.xticks(x,['PSP \n $\sim0.05 \  AU$','SOLO \n $\sim0.3 \  AU$','WIND \n $\sim1 \ AU$	'],fontsize=12)
plt.legend(fontsize=15)



#%%
# Magnetic Fluctuation Ratio
fastratio_err = (dBfast/Bfast)*np.sqrt((dBfast_std/dBfast)**2 + (Bfast_std/Bfast)**2)
slowsratio_err = (dBslow/Bslow)*np.sqrt((dBslow_std/dBslow)**2 + (Bslow_std/Bslow)**2)
irratio_err = (dBir/Bir)*np.sqrt((dBir_std/dBir)**2 + (Bir_std/Bir)**2)

plt.errorbar(x,dBfast/Bfast,yerr=fastratio_err,fmt='o-',color = 'r',label = r'$|\delta B/B|_{fast}$',capsize=5)
plt.errorbar(x,dBslow/Bslow,yerr=slowsratio_err,fmt='o-',color='b',label = r'$|\delta B/B|_{slow}$',capsize=5)
plt.errorbar(x,dBir/Bir,yerr=irratio_err,fmt='o-',color='purple',label = r'$|\delta B/B|_{ir}$',capsize=5)
plt.plot(x,(dBslow/Bslow+dBfast/Bfast)/2,color='g',marker='s',linestyle='dashed',label = r'$|\delta B/B|_{avg}$')
plt.title('$\delta B/B$',fontsize=15)
plt.axvline(x=1,linestyle='dashed',color='k',linewidth=0.2)
plt.axvline(x=2,linestyle='dashed',color='k',linewidth=0.2)
plt.axvline(x=3,linestyle='dashed',color='k',linewidth=0.2)
plt.yticks(fontsize=12)
plt.xticks(x,['PSP \n $\sim0.05 \  AU$','SOLO \n $\sim0.3 \  AU$','WIND \n $\sim1 \ AU$	'],fontsize=12)
plt.legend(fontsize=15)



#%%
# Bulk Energy Density
bulkfast_err = bulkfast*np.sqrt((nfast_std/nfast)**2 + 2*(vfast_std/vfast)**2)
bulkslow_err = bulkslow*np.sqrt((nslow_std/nslow)**2 + 2*(vslow_std/vslow)**2)
bulkir_err = bulkir*np.sqrt((nir_std/nir)**2 + 2*(vir_std/vir)**2)

plt.errorbar(x,bulkfast,yerr=bulkfast_err,fmt='o-',color = 'r',label = r'Fast Wind',capsize=5)
plt.errorbar(x,bulkslow,yerr=bulkslow_err,fmt='o-',color='b',label = r'Slow Wind',capsize=5)
plt.errorbar(x,bulkir,yerr=bulkir_err,fmt='o-',color='purple',label = r'IR',capsize=5)
plt.plot(x,(bulkfast+bulkslow)/2,color='g',marker='s',linestyle='dashed',label = r'Mean')
plt.title('Adjusted Bulk Energy Density [J/m^3]',fontsize=15)
plt.axvline(x=1,linestyle='dashed',color='k',linewidth=0.2)
plt.axvline(x=2,linestyle='dashed',color='k',linewidth=0.2)
plt.axvline(x=3,linestyle='dashed',color='k',linewidth=0.2)
plt.xticks(x,['PSP','SOLO','WIND'])
plt.legend(fontsize=15)

#%%
# Thermal Energy Density
thermfast_err = thermfast*np.sqrt((nfast_std/nfast)**2 + (Tfast_std/Tfast)**2)
thermslow_err = thermslow*np.sqrt((nslow_std/nslow)**2 + (Tslow_std/Tslow)**2)
thermir_err = thermir*np.sqrt((nir_std/nir)**2 + (Tir_std/Tir)**2)

plt.errorbar(x,thermfast,yerr=thermfast_err,fmt='o-',color = 'r',label = r'Fast Wind',capsize=5)
plt.errorbar(x,thermslow,yerr=thermslow_err,fmt='o-',color='b',label = r'Slow Wind',capsize=5)
plt.errorbar(x,thermir,yerr=thermir_err,fmt='o-',color='purple',label = r'IR',capsize=5)
plt.plot(x,(thermfast+thermslow)/2,color='g',marker='s',linestyle='dashed',label = r'Mean')
plt.title('Adjusted Thermal Energy Density',fontsize=15)
plt.axvline(x=1,linestyle='dashed',color='k',linewidth=0.2)
plt.axvline(x=2,linestyle='dashed',color='k',linewidth=0.2)
plt.axvline(x=3,linestyle='dashed',color='k',linewidth=0.2)
plt.xticks(x,['PSP','SOLO','WIND'])
plt.legend(fontsize=15)

# %%
# Magnetic Energy Density
magfast_err = magfast*np.sqrt(2*(Bfast_std/Bfast)**2)
magslow_err = magslow*np.sqrt(2*(Bslow_std/Bslow)**2)
magir_err = magir*np.sqrt(2*(Bir_std/Bir)**2)
plt.errorbar(x,magfast,yerr=magfast_err,fmt='o-',color = 'r',label = r'Fast Wind',capsize=5)
plt.errorbar(x,magslow,yerr=magslow_err,fmt='o-',color='b',label = r'Slow Wind',capsize=5)
plt.errorbar(x,magir,yerr=magir_err,fmt='o-',color='purple',label = r'IR',capsize=5)
plt.plot(x,(magfast+magslow)/2,color='g',marker='s',linestyle='dashed',label = r'Mean')
plt.title('Adjusted Magnetic Energy Density')
plt.ylim(0,2e-12)
plt.axvline(x=1,linestyle='dashed',color='k',linewidth=0.2)
plt.axvline(x=2,linestyle='dashed',color='k',linewidth=0.2)
plt.axvline(x=3,linestyle='dashed',color='k',linewidth=0.2)
plt.xticks(x,['PSP','SOLO','WIND'])
plt.xlabel('Spacecraft')
plt.legend(fontsize=15)

#%%
#Magnetic Fluctuation Energy Density
deltamagfast_err = deltamagfast*np.sqrt(2*(dBfast_std/dBfast)**2)
deltamagslow_err = deltamagslow*np.sqrt(2*(dBslow_std/dBslow)**2)
deltamagir_err = deltamagir*np.sqrt(2*(dBir_std/dBir)**2)
plt.errorbar(x,deltamagfast,yerr=deltamagfast_err,fmt='o-',color = 'r',label = r'Fast Wind',capsize=5)
plt.errorbar(x,deltamagslow,yerr=deltamagslow_err,fmt='o-',color='b',label = r'Slow Wind',capsize=5)
plt.errorbar(x,deltamagir,yerr=deltamagir_err,fmt='o-',color='purple',label = r'IR',capsize=5)
plt.plot(x,(deltamagfast+deltamagslow)/2,color='g',marker='s',linestyle='dashed',label = r'Mean')
plt.title('Magnetic Fluctuation Energy Density [J/m^3]')
plt.axvline(x=1,linestyle='dashed',color='k',linewidth=0.2)
plt.axvline(x=2,linestyle='dashed',color='k',linewidth=0.2)
plt.axvline(x=3,linestyle='dashed',color='k',linewidth=0.2)
plt.xticks(x,['PSP','SOLO','WIND'])
plt.xlabel('Spacecraft')
plt.legend(fontsize=15)
# %%
