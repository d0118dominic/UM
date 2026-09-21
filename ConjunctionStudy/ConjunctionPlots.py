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

# def brazilhist(z, x_label=r'$\beta_\parallel$', y_label=r'$T_{\perp}/T_{\parallel}$',
#      z_label='', title='',nbins=50,mincount=10,vmin=1,vmax=10000, mask = allTparperp>0,
#      count=False,scale = 'linear',fitline=False,xrange = [-3,3],yrange = [-1,1],cmap='jet'):    
#     mirror_params = [0.77, 0.76, -0.016]  # [a, b, beta0]
#     firehose_params = [-1.4, 1.0, -0.11]  # [a, b, beta0]
#     cyclotron_params = [0.45, 0.4, -0.0004]  # [a, b, beta0]
#     parfirehose_params = [-0.47, 0.53, 0.59]  # [a, b, beta0]

#     beta_par_range = np.logspace(xrange[0], xrange[1], 100)  # Adjust range and resolution as needed
#     mirror_threshold = stability_condition(beta_par_range, *mirror_params)
#     firehose_threshold = stability_condition(beta_par_range, *firehose_params)
#     cyclotron_threshold = stability_condition(beta_par_range, *cyclotron_params)
#     parfirehose_threshold = stability_condition(beta_par_range, *parfirehose_params)

#     # plt.figure(figsize=(10, 8))
#     # sc = plt.scatter(allbeta_par[mask], 1/allTparperp[mask], c=z[mask], cmap='plasma', s=0.01)
#     x = allbeta_par[mask]
#     y = 1/allTparperp[mask]
#     z = z[mask]
#     valmask = np.isfinite(x) & np.isfinite(y) & np.isfinite(z) 
#     x_valid = x[valmask]
#     y_valid = y[valmask]
#     z_valid = z[valmask]

#     x_fit = np.log10(x_valid)
#     y_fit = np.log10(y_valid)
#     m,c = np.polyfit(x_fit, y_fit, 1)


#     nbins = nbins
#     mincount=mincount
#     # xbins = np.linspace(0.001, 100, nbins)
#     # ybins = np.linspace(0.1, 10, nbins)
#     xbins = np.logspace(xrange[0], xrange[1], nbins)
#     ybins = np.logspace(yrange[0], yrange[1], nbins)

#     xline = np.logspace(-3, 3, 200)
#     yline = 10**(m*np.log10(xline) + c)

#     if count==True:
#         hist, xedges, yedges, _= binned_statistic_2d(x_valid, y_valid, z_valid, statistic='count', bins=[xbins, ybins])
#     else:
#         hist, xedges, yedges, _= binned_statistic_2d(x_valid, y_valid, z_valid, statistic='mean', bins=[xbins, ybins])

#     counts, _, _, _= binned_statistic_2d(x_valid, y_valid, z_valid, statistic='count', bins=[xbins, ybins])
#     zero_mask = counts <= mincount
#     hist_masked = np.ma.masked_where(counts <= mincount, hist)

#     fig, ax = plt.subplots(figsize=(10, 8))

#     # if np.any(zero_mask):
#         # grey_data = np.where(zero_mask, 1.0, np.nan)
#         # ax.pcolormesh(xedges, yedges, grey_data.T, cmap='gray', vmin=0, vmax=1)
#     if scale == 'log':
#         norm = colors.LogNorm(vmin=vmin, vmax=vmax)
#     else:
#         norm = colors.Normalize(vmin=vmin, vmax=vmax)

#     im = ax.pcolormesh(xedges, yedges, hist_masked.T, cmap=cmap, norm=norm)

#     cbar = plt.colorbar(im, ax=ax)
#     cbar.set_label(z_label,fontdict={'size': 15})
#     if count==True:
#         cbar.set_label('Count',fontdict={'size': 15})



#     if fitline==True:
#         ax.plot(xline, yline, 'r-', linewidth=3,label=rf'Best fit: $\beta_\parallel^{{{m:.2f}}}$',color='b',linestyle='-')
#     else:
#         pass
  

#     # plt.colorbar(im, label=z_label)
#     ax.set_facecolor('lightgray')
#     ax.set_xlabel(x_label, fontdict={'size': 15})
#     ax.set_ylabel(y_label, fontdict={'size': 15})
#     ax.set_title(title, fontdict={'size': 15})
#     ax.set_xlim(10**(xrange[0]), 10**(xrange[1]))  # Adjust as needed
#     ax.set_ylim(10**(-1), 10**(1))  # Adjust as needed
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

def brazilhist(z, x_label=r'$\beta_\parallel$', y_label=r'$T_{\perp}/T_{\parallel}$',
     z_label='', title='',nbins=50,mincount=10,vmin=1,vmax=10000, mask=None,
     count=False,scale = 'linear',fitline=False,xrange = [-3,3],yrange = [-1,1],cmap='jet'):    
    mirror_params = [0.77, 0.76, -0.016]  # [a, b, beta0]
    firehose_params = [-1.4, 1.0, -0.11]  # [a, b, beta0]
    cyclotron_params = [0.45, 0.4, -0.0004]  # [a, b, beta0]
    parfirehose_params = [-0.47, 0.53, 0.59]  # [a, b, beta0]

    if mask is None:
        mask = np.ones(allTparperp.shape, dtype=bool)

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
                  scatter_max_points=None, mean_marker_size=120,
                  cgl_line=True, cgl_exponent=-1):

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

    xline_cgl = np.logspace(xrange[0], xrange[1], 100)

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

            if cgl_line:
                yline_cgl = y_mean * (xline_cgl / x_mean) ** cgl_exponent
                ax.plot(xline_cgl, yline_cgl, color=color, linestyle=':', linewidth=1.5, alpha=0.8, zorder=9)

        elif plot_type == 'means_only':
            # just the mean marker, no contour or scatter cloud
            x_mean = np.mean(x_valid)
            y_mean = np.mean(y_valid)
            ax.scatter(x_mean, y_mean, marker='^', color=color, s=4*mean_marker_size,
                       edgecolors='k', linewidths=1.5, zorder=10, label=sc.upper())

            if cgl_line:
                yline_cgl = y_mean * (xline_cgl / x_mean) ** cgl_exponent
                ax.plot(xline_cgl, yline_cgl, color=color, linestyle=':', linewidth=1.5, alpha=0.8, zorder=9)

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

                if cgl_line:
                    yline_cgl = y_mean * (xline_cgl / x_mean) ** cgl_exponent
                    ax.plot(xline_cgl, yline_cgl, color=color, linestyle='-', linewidth=1, alpha=0.8, zorder=9)

    ax.set_facecolor('lightgray')
    ax.set_xlabel(x_label, fontdict={'size': 20})
    ax.set_ylabel(y_label, fontdict={'size': 20})
    ax.set_title(title, fontdict={'size': 20})
    ax.set_xlim(10**(xrange[0]), 10**(xrange[1]))
    ax.set_ylim(10**(-1), 10**(1))
    ax.axhline(y=1, color='k', linestyle='-')
    ax.axvline(x=1, color='k', linestyle='-')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.tick_params(labelsize=15)
    ax.plot(beta_par_range, mirror_threshold, label='Mirror', color='k', linestyle='dotted', linewidth=2)
    ax.plot(beta_par_range, firehose_threshold, label='Oblique FH', color='k', linestyle='--', linewidth=2)
    ax.plot(beta_par_range, cyclotron_threshold, label='IC', color='k', linestyle='dotted', linewidth=1)
    ax.plot(beta_par_range, parfirehose_threshold, label='Parallel FH', color='k', linestyle='--', linewidth=1)
    ax.legend(fontsize=17,framealpha=1)
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



brazilhist(1.602e19*abs(allT), z_label = r'$T \ (eV)$',fitline=False,nbins=200,vmin=00,vmax=200,count=False,scale='linear')
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
                        contour_nlevels=10,title = "Leading Slow Wind",plot_type='contour',
                        mean_marker_size=35,cgl_exponent=-0.55)

#%%
fig, ax = brazil_multi(betapar_psp=betaparspsp_trailing, Tparperp_psp=Tparperppsp_trailing,
                        betapar_solo=betaparssolo_trailing, Tparperp_solo=Tparperpsolo_trailing,
                        betapar_wind=betaparswind_trailing, Tparperp_wind=Tparperpwind_trailing,
                        contour_nlevels=10,title = "Trailing Slow Wind",plot_type='contour',
                        mean_marker_size=35)
#%%
fig, ax = brazil_multi(betapar_psp=betaparspsp_shoulder, Tparperp_psp=Tparperppsp_shoulder,
                        betapar_solo=betaparssolo_shoulder, Tparperp_solo=Tparperpsolo_shoulder,
                        betapar_wind=betaparswind_shoulder, Tparperp_wind=Tparperpwind_shoulder,
                        contour_nlevels=10,title = "Compression Region",plot_type='contour',
                        mean_marker_size=35)

#%%
fig, ax = brazil_multi(betapar_psp=betaparspsp_fast, Tparperp_psp=Tparperppsp_fast,
                        betapar_solo=betaparssolo_fast, Tparperp_solo=Tparperpsolo_fast,
                        betapar_wind=betaparswind_fast, Tparperp_wind=Tparperpwind_fast,
                        contour_nlevels=10,title = "Fast CH Wind",plot_type='contour',
                        mean_marker_size=35,fitline=True,cgl_exponent=-1)

#%%
fig, ax = brazil_multi(betapar_psp=betaparspsp_rarefaction, Tparperp_psp=Tparperppsp_rarefaction,
                        betapar_solo=betaparssolo_rarefaction, Tparperp_solo=Tparperpsolo_rarefaction,
                        betapar_wind=betaparswind_rarefaction, Tparperp_wind=Tparperpwind_rarefaction,
                        contour_nlevels=10,title = "Rarefaction Region",plot_type='contour',
                        mean_marker_size=35)

#%%

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
        charray = np.array([chpsp_fast,chsolo_fast,chwind_fast])
        ch_stdarray = np.array([chstdpsp_fast,chstdsolo_fast,chstdwind_fast])
        rearray = np.array([repsp_fast,resolo_fast,rewind_fast])
        re_stdarray = np.array([restdpsp_fast,restdsolo_fast,restdwind_fast])
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
        charray = np.array([chpsp_slow,chsolo_slow,chwind_slow])
        ch_stdarray = np.array([chstdpsp_slow,chstdsolo_slow,chstdwind_slow])
        rearray = np.array([repsp_slow,resolo_slow,rewind_slow])
        re_stdarray = np.array([restdpsp_slow,restdsolo_slow,restdwind_slow])
    elif type == 'preceding':
        narray = 1e6*np.array([npsp_preceding,nsolo_preceding,nwind_preceding])
        n_stdarray = 1e6*np.array([nstdpsp_preceding,nstdsolo_preceding,nstdwind_preceding])
        Tarray = 1.602e-19*np.array([Tpsp_preceding,Tsolo_preceding,Twind_preceding])
        T_stdarray = 1.602e-19*np.array([Tstdpsp_preceding,Tstdsolo_preceding,Tstdwind_preceding])
        varray = 1e3*np.array([vpsp_preceding,vsolo_preceding,vwind_preceding])
        v_stdarray = 1e3*np.array([vstdpsp_preceding,vstdsolo_preceding,vstdwind_preceding])
        Barray = 1e-9*np.array([Bpsp_preceding,Bsolo_preceding,Bwind_preceding])
        B_stdarray = 1e-9*np.array([Bstdpsp_preceding,Bstdsolo_preceding,Bstdwind_preceding])
        dB_array = 1e-9*np.array([dBpsp_preceding,dBsolo_preceding,dBwind_preceding])
        dB_stdarray = 1e-9*np.array([dBstdpsp_preceding,dBstdsolo_preceding,dBstdwind_preceding])
        charray = np.array([chpsp_preceding,chsolo_preceding,chwind_preceding])
        ch_stdarray = np.array([chstdpsp_preceding,chstdsolo_preceding,chstdwind_preceding])
        rearray = np.array([repsp_preceding,resolo_preceding,rewind_preceding])
        re_stdarray = np.array([restdpsp_preceding,restdsolo_preceding,restdwind_preceding])
    elif type == 'trailing':
        narray = 1e6*np.array([npsp_trailing,nsolo_trailing,nwind_trailing])
        n_stdarray = 1e6*np.array([nstdpsp_trailing,nstdsolo_trailing,nstdwind_trailing])
        Tarray = 1.602e-19*np.array([Tpsp_trailing,Tsolo_trailing,Twind_trailing])
        T_stdarray = 1.602e-19*np.array([Tstdpsp_trailing,Tstdsolo_trailing,Tstdwind_trailing])
        varray = 1e3*np.array([vpsp_trailing,vsolo_trailing,vwind_trailing])
        v_stdarray = 1e3*np.array([vstdpsp_trailing,vstdsolo_trailing,vstdwind_trailing])
        Barray = 1e-9*np.array([Bpsp_trailing,Bsolo_trailing,Bwind_trailing])
        B_stdarray = 1e-9*np.array([Bstdpsp_trailing,Bstdsolo_trailing,Bstdwind_trailing])
        dB_array = 1e-9*np.array([dBpsp_trailing,dBsolo_trailing,dBwind_trailing])
        dB_stdarray = 1e-9*np.array([dBstdpsp_trailing,dBstdsolo_trailing,dBstdwind_trailing])
        charray = np.array([chpsp_trailing,chsolo_trailing,chwind_trailing])
        ch_stdarray = np.array([chstdpsp_trailing,chstdsolo_trailing,chstdwind_trailing])
        rearray = np.array([repsp_trailing,resolo_trailing,rewind_trailing])
        re_stdarray = np.array([restdpsp_trailing,restdsolo_trailing,restdwind_trailing])
    elif type == 'shoulder':
        narray = 1e6*np.array([npsp_shoulder,nsolo_shoulder,nwind_shoulder])
        n_stdarray = 1e6*np.array([nstdpsp_shoulder,nstdsolo_shoulder,nstdwind_shoulder])
        Tarray = 1.602e-19*np.array([Tpsp_shoulder,Tsolo_shoulder,Twind_shoulder])
        T_stdarray = 1.602e-19*np.array([Tstdpsp_shoulder,Tstdsolo_shoulder,Tstdwind_shoulder])
        varray = 1e3*np.array([vpsp_shoulder,vsolo_shoulder,vwind_shoulder])
        v_stdarray = 1e3*np.array([vstdpsp_shoulder,vstdsolo_shoulder,vstdwind_shoulder])
        Barray = 1e-9*np.array([Bpsp_shoulder,Bsolo_shoulder,Bwind_shoulder])
        B_stdarray = 1e-9*np.array([Bstdpsp_shoulder,Bstdsolo_shoulder,Bstdwind_shoulder])
        dB_array = 1e-9*np.array([dBpsp_shoulder,dBsolo_shoulder,dBwind_shoulder])
        dB_stdarray = 1e-9*np.array([dBstdpsp_shoulder,dBstdsolo_shoulder,dBstdwind_shoulder])
        charray = np.array([chpsp_shoulder,chsolo_shoulder,chwind_shoulder])
        ch_stdarray = np.array([chstdpsp_shoulder,chstdsolo_shoulder,chstdwind_shoulder])
        rearray = np.array([repsp_shoulder,resolo_shoulder,rewind_shoulder])
        re_stdarray = np.array([restdpsp_shoulder,restdsolo_shoulder,restdwind_shoulder])
    elif type == 'rarefaction':
        narray = 1e6*np.array([npsp_rarefaction,nsolo_rarefaction,nwind_rarefaction])
        n_stdarray = 1e6*np.array([nstdpsp_rarefaction,nstdsolo_rarefaction,nstdwind_rarefaction])
        Tarray = 1.602e-19*np.array([Tpsp_rarefaction,Tsolo_rarefaction,Twind_rarefaction])
        T_stdarray = 1.602e-19*np.array([Tstdpsp_rarefaction,Tstdsolo_rarefaction,Tstdwind_rarefaction])
        varray = 1e3*np.array([vpsp_rarefaction,vsolo_rarefaction,vwind_rarefaction])
        v_stdarray = 1e3*np.array([vstdpsp_rarefaction,vstdsolo_rarefaction,vstdwind_rarefaction])
        Barray = 1e-9*np.array([Bpsp_rarefaction,Bsolo_rarefaction,Bwind_rarefaction])
        B_stdarray = 1e-9*np.array([Bstdpsp_rarefaction,Bstdsolo_rarefaction,Bstdwind_rarefaction])
        dB_array = 1e-9*np.array([dBpsp_rarefaction,dBsolo_rarefaction,dBwind_rarefaction])
        dB_stdarray = 1e-9*np.array([dBstdpsp_rarefaction,dBstdsolo_rarefaction,dBstdwind_rarefaction])
        charray = np.array([chpsp_rarefaction,chsolo_rarefaction,chwind_rarefaction])
        ch_stdarray = np.array([chstdpsp_rarefaction,chstdsolo_rarefaction,chstdwind_rarefaction])
        rearray = np.array([repsp_rarefaction,resolo_rarefaction,rewind_rarefaction])
        re_stdarray = np.array([restdpsp_rarefaction,restdsolo_rarefaction,restdwind_rarefaction])
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
        charray = np.array([chpsp,chsolo,chwind])
        ch_stdarray = np.array([chstdpsp,chstdsolo,chstdwind])
        rearray = np.array([repsp,resolo,rewind])
        re_stdarray = np.array([restdpsp,restdsolo,restdwind])
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

    return narray_norm,Tarray_norm,varray_norm,Barray_norm,dB_array_norm,charray,rearray, n_stdarray_norm, T_stdarray_norm, v_stdarray_norm, B_stdarray_norm, dB_stdarray_norm, ch_stdarray, re_stdarray



nfast,Tfast,vfast,Bfast,dBfast,chfast,refast,nfast_std,Tfast_std,vfast_std,Bfast_std,dBfast_std,chfast_std,refast_std = get_vars('fast')
npreceding,Tpreceding,vpreceding,Bpreceding,dBpreceding,chpreceding,repreceding,npreceding_std,Tpreceding_std,vpreceding_std,Bpreceding_std,dBpreceding_std,chpreceding_std,repreceding_std = get_vars('preceding')
ntrailing,Ttrailing,vtrailing,Btrailing,dBtrailing,chtrailing,retrailing,ntrailing_std,Ttrailing_std,vtrailing_std,Btrailing_std,dBtrailing_std,chtrailing_std,retrailing_std = get_vars('trailing')
nshoulder,Tshoulder,vshoulder,Bshoulder,dBshoulder,chshoulder,reshoulder,nshoulder_std,Tshoulder_std,vshoulder_std,Bshoulder_std,dBshoulder_std,chshoulder_std,reshoulder_std = get_vars('shoulder')
nrarefaction,Trarefaction,vrarefaction,Brarefaction,dBrarefaction,chrarefaction,rerarefaction,nrarefaction_std,Trarefaction_std,vrarefaction_std,Brarefaction_std,dBrarefaction_std,chrarefaction_std,rerarefaction_std = get_vars('rarefaction')
# nslow,Tslow,vslow,Bslow,dBslow,nslow_std,Tslow_std,vslow_std,Bslow_std,dBslow_std = get_vars('slow')
# nir,Tir,vir,Bir,dBir,nir_std,Tir_std,vir_std,Bir_std,dBir_std = get_vars('ir')
# n = nfast + nslow
# T = (Tfast + Tslow)/2

bulkfast = 0.5*nfast*mi*vfast**2
bulkpreceding = 0.5*npreceding*mi*vpreceding**2
bulktrailing = 0.5*ntrailing*mi*vtrailing**2
bulkshoulder = 0.5*nshoulder*mi*vshoulder**2
bulkrarefaction = 0.5*nrarefaction*mi*vrarefaction**2

thermfast = (3/2)*kb*nfast*Tfast
thermpreceding = (3/2)*kb*npreceding*Tpreceding
thermtrailing = (3/2)*kb*ntrailing*Ttrailing
thermshoulder = (3/2)*kb*nshoulder*Tshoulder
thermrarefaction = (3/2)*kb*nrarefaction*Trarefaction

magfast = (0.5*Bfast**2)/(mu0)
magpreceding = (0.5*Bpreceding**2)/(mu0)
magtrailing = (0.5*Btrailing**2)/(mu0)
magshoulder = (0.5*Bshoulder**2)/(mu0)
magrarefaction = (0.5*Brarefaction**2)/(mu0)
# magmean = (magfast + magpreceding + magtrailing + magshoulder + magrarefaction)/5

deltamagfast = (0.5*dBfast**2)/(mu0)
deltamagpreceding = (0.5*dBpreceding**2)/(mu0)
deltamagtrailing = (0.5*dBtrailing**2)/(mu0)
deltamagshoulder = (0.5*dBshoulder**2)/(mu0)
deltamagrarefaction = (0.5*dBrarefaction**2)/(mu0)

plasmafast = thermfast + bulkfast
plasmapreceding = thermpreceding + bulkpreceding
plasmatrailing = thermtrailing + bulktrailing
plasmashoulder = thermshoulder + bulkshoulder
plasmararefaction = thermrarefaction + bulkrarefaction

plasma = plasmafast + plasmapreceding + plasmatrailing + plasmashoulder + plasmararefaction

fieldfast = magfast + deltamagfast
fieldpreceding = magpreceding + deltamagpreceding
fieldtrailing = magtrailing + deltamagtrailing
fieldshoulder = magshoulder + deltamagshoulder
fieldrarefaction = magrarefaction + deltamagrarefaction

field = fieldfast + fieldpreceding + fieldtrailing + fieldshoulder + fieldrarefaction
# %%


# Number Density
plt.errorbar(x,npreceding,yerr=npreceding_std,fmt='o-',color='darkblue',label = r'$n^{\prime}_{leading}$',capsize=5)
plt.errorbar(x,nshoulder,yerr=nshoulder_std,fmt='o-',color='purple',label = r'$n^{\prime}_{compression}$',capsize=5)
plt.errorbar(x,nfast,yerr=nfast_std,fmt='o-',color = 'r',label = r'$n^{\prime}_{fast}$',capsize=5)
plt.errorbar(x,nrarefaction,yerr=nrarefaction_std,fmt='o-',color='green',label = r'$n^{\prime}_{rarefaction}$',capsize=5)
plt.errorbar(x,ntrailing,yerr=ntrailing_std,fmt='o-',color='skyblue',label = r'$n^{\prime}_{trailing}$',capsize=5)
# plt.plot(x,(nfast+npreceding+ntrailing+nshoulder+nrarefaction)/5,color='g',marker='s',linestyle='dashed',label = r'$n^{\prime}_{avg}$')
plt.title('Adjusted Number Density',fontsize=15)
plt.axvline(x=1,linestyle='dashed',color='k',linewidth=0.2)
plt.axvline(x=2,linestyle='dashed',color='k',linewidth=0.2)
plt.axvline(x=3,linestyle='dashed',color='k',linewidth=0.2)
plt.xticks(x,['PSP \n $\sim0.05 \  AU$','SOLO \n $\sim0.3 \  AU$','WIND \n $\sim1 \ AU$	'],fontsize=12)
# plt.legend(fontsize=12)

#%%
# Temperature
plt.errorbar(x,Tpreceding,yerr=Tpreceding_std,fmt='o-',color='darkblue',label = r'$T^{\prime}_{leading}$',capsize=5)
plt.errorbar(x,Tshoulder,yerr=Tshoulder_std,fmt='o-',color='purple',label = r'$T^{\prime}_{compression}$',capsize=5)
plt.errorbar(x,Tfast,yerr=Tfast_std,fmt='o-',color = 'r',label = r'$T^{\prime}_{fast}$',capsize=5)
plt.errorbar(x,Trarefaction,yerr=Trarefaction_std,fmt='o-',color='green',label = r'$T^{\prime}_{rarefaction}$',capsize=5)
plt.errorbar(x,Ttrailing,yerr=Ttrailing_std,fmt='o-',color='skyblue',label = r'$T^{\prime}_{trailing}$',capsize=5)
# plt.plot(x,(Tfast+Tpreceding+Ttrailing+Tshoulder+Trarefaction)/5,color='g',marker='s',linestyle='dashed',label = r'$T^{\prime}_{avg}$')
plt.title('Adjusted Temperature',fontsize=15)
plt.yscale('log')
plt.axvline(x=1,linestyle='dashed',color='k',linewidth=0.2)
plt.axvline(x=2,linestyle='dashed',color='k',linewidth=0.2)
plt.axvline(x=3,linestyle='dashed',color='k',linewidth=0.2)
plt.xticks(x,['PSP \n $\sim0.05 \  AU$','SOLO \n $\sim0.3 \  AU$','WIND \n $\sim1 \ AU$	'],fontsize=12)
# plt.legend(fontsize=12)



#%%
# Velocity
plt.errorbar(x,1e-3*vpreceding,yerr=1e-3*vpreceding_std,fmt='o-',color='darkblue',label = r'$v_{leading}$',capsize=5)
plt.errorbar(x,1e-3*vshoulder,yerr=1e-3*vshoulder_std,fmt='o-',color='purple',label = r'$v_{compression}$',capsize=5)
plt.errorbar(x,1e-3*vfast,yerr=1e-3*vfast_std,fmt='o-',color = 'r',label = r'$v_{fast}$',capsize=5)
plt.errorbar(x,1e-3*vrarefaction,yerr=1e-3*vrarefaction_std,fmt='o-',color='green',label = r'$v_{rarefaction}$',capsize=5)
plt.errorbar(x,1e-3*vtrailing,yerr=1e-3*vtrailing_std,fmt='o-',color='skyblue',label = r'$v_{trailing}$',capsize=5)
# plt.errorbar(x,1e-3*(vfast+vpreceding+vtrailing+vshoulder+vrarefaction)/5,color = 'g',marker='s',linestyle='dashed',label = r'$v_{avg}$')

# plt.ylim(200,800)
plt.title('Velocity [km/s]',fontsize=15)
plt.axvline(x=1,linestyle='dashed',color='k',linewidth=0.2)
plt.axvline(x=2,linestyle='dashed',color='k',linewidth=0.2)
plt.axvline(x=3,linestyle='dashed',color='k',linewidth=0.2)
plt.ylim(150,900)
plt.yticks(fontsize=12)
plt.xticks(x,['PSP \n $\sim0.05 \  AU$','SOLO \n $\sim0.3 \  AU$','WIND \n $\sim1 \ AU$	'],fontsize=12)
# plt.legend(fontsize=12)

#%%
# Magnetic Field Strength
plt.errorbar(x,Bpreceding,yerr=Bpreceding_std,fmt='o-',color='darkblue',label = r'$|B|_{leading}$',capsize=5)
plt.errorbar(x,Bshoulder,yerr=Bshoulder_std,fmt='o-',color='purple',label = r'$|B|_{compression}$',capsize=5)
plt.errorbar(x,Bfast,yerr=Bfast_std,fmt='o-',color = 'r',label = r'$|B|_{fast}$',capsize=5)
plt.errorbar(x,Brarefaction,yerr=Brarefaction_std,fmt='o-',color='green',label = r'$|B|_{rarefaction}$',capsize=5)
plt.errorbar(x,Btrailing,yerr=Btrailing_std,fmt='o-',color='skyblue',label = r'$|B|_{trailing}$',capsize=5)
# plt.errorbar(x,(Bfast+Bpreceding+Btrailing+Bshoulder+Brarefaction)/5,color = 'g',marker='s',linestyle='dashed',label = r'$|B|_{avg}$')
plt.title('Adjusted Magnetic Field Strength',fontsize=15)
# plt.ylim(0,2.5e-9)
plt.axvline(x=1,linestyle='dashed',color='k',linewidth=0.2)
plt.axvline(x=2,linestyle='dashed',color='k',linewidth=0.2)
plt.axvline(x=3,linestyle='dashed',color='k',linewidth=0.2)
plt.yticks(fontsize=12)
plt.yscale('log')
plt.xticks(x,['PSP \n $\sim0.05 \  AU$','SOLO \n $\sim0.3 \  AU$','WIND \n $\sim1 \ AU$	'],fontsize=12)
# plt.legend(fontsize=12)

#%%
# Magnetic Fluctuation Strength
# plt.errorbar(x,dBfast,yerr=dBfast_std,fmt='o-',color = 'r',label = r'$\delta B_{fast}$',capsize=5)
# plt.errorbar(x,dBpreceding,yerr=dBpreceding_std,fmt='o-',color='darkblue',label = r'$\delta B_{leading}$',capsize=5)
# plt.errorbar(x,dBtrailing,yerr=dBtrailing_std,fmt='o-',color='skyblue',label = r'$\delta B_{trailing}$',capsize=5)
# plt.errorbar(x,dBshoulder,yerr=dBshoulder_std,fmt='o-',color='purple',label = r'$\delta B_{compression}$',capsize=5)
# plt.errorbar(x,dBrarefaction,yerr=dBrarefaction_std,fmt='o-',color='green',label = r'$\delta B_{rarefaction}$',capsize=5)
# # plt.errorbar(x,(dBfast+dBpreceding+dBtrailing+dBshoulder+dBrarefaction)/5,color = 'g',marker='s',linestyle='dashed',label = r'$\delta B_{avg}$')
# # plt.ylim(0,2.5e-9)
# plt.title('Adjusted Magnetic Fluctuation Magnitude [T]')
# plt.axvline(x=1,linestyle='dashed',color='k',linewidth=0.2)
# plt.axvline(x=2,linestyle='dashed',color='k',linewidth=0.2)
# plt.axvline(x=3,linestyle='dashed',color='k',linewidth=0.2)
# plt.yticks(fontsize=12)
# plt.xticks(x,['PSP \n $\sim0.05 \  AU$','SOLO \n $\sim0.3 \  AU$','WIND \n $\sim1 \ AU$	'],fontsize=12)
# plt.legend(fontsize=12)

#%%
# Magnetic Fluctuation Ratio
# fastratio_err = (dBfast/Bfast)*np.sqrt((dBfast_std/dBfast)**2 + (Bfast_std/Bfast)**2)
# slowsratio_err = (dBslow/Bslow)*np.sqrt((dBslow_std/dBslow)**2 + (Bslow_std/Bslow)**2)
# irratio_err = (dBir/Bir)*np.sqrt((dBir_std/dBir)**2 + (Bir_std/Bir)**2)

# plt.errorbar(x,dBfast/Bfast,yerr=fastratio_err,fmt='o-',color = 'r',label = r'$|\delta B/B|_{fast}$',capsize=5)
# plt.errorbar(x,dBslow/Bslow,yerr=slowsratio_err,fmt='o-',color='b',label = r'$|\delta B/B|_{slow}$',capsize=5)
# plt.errorbar(x,dBir/Bir,yerr=irratio_err,fmt='o-',color='purple',label = r'$|\delta B/B|_{ir}$',capsize=5)
# plt.plot(x,(dBslow/Bslow+dBfast/Bfast)/2,color='g',marker='s',linestyle='dashed',label = r'$|\delta B/B|_{avg}$')
# plt.title('$\delta B/B$',fontsize=15)
# plt.axvline(x=1,linestyle='dashed',color='k',linewidth=0.2)
# plt.axvline(x=2,linestyle='dashed',color='k',linewidth=0.2)
# plt.axvline(x=3,linestyle='dashed',color='k',linewidth=0.2)
# plt.yticks(fontsize=12)
# plt.xticks(x,['PSP \n $\sim0.05 \  AU$','SOLO \n $\sim0.3 \  AU$','WIND \n $\sim1 \ AU$	'],fontsize=12)
# plt.legend(fontsize=15)



#%%
# Bulk Energy Density
bulkfast_err = bulkfast*np.sqrt((nfast_std/nfast)**2 + 2*(vfast_std/vfast)**2)
bulkpreceding_err = bulkpreceding*np.sqrt((npreceding_std/npreceding)**2 + 2*(vpreceding_std/vpreceding)**2)
bulktrailing_err = bulktrailing*np.sqrt((ntrailing_std/ntrailing)**2 + 2*(vtrailing_std/vtrailing)**2)
bulkshoulder_err = bulkshoulder*np.sqrt((nshoulder_std/nshoulder)**2 + 2*(vshoulder_std/vshoulder)**2)
bulkrarefaction_err = bulkrarefaction*np.sqrt((nrarefaction_std/nrarefaction)**2 + 2*(vrarefaction_std/vrarefaction)**2)

plt.errorbar(x,bulkpreceding,yerr=bulkpreceding_err,fmt='o-',color='darkblue',label = r'Leading',capsize=5)
plt.errorbar(x,bulkshoulder,yerr=bulkshoulder_err,fmt='o-',color='purple',label = r'Compression',capsize=5)
plt.errorbar(x,bulkfast,yerr=bulkfast_err,fmt='o-',color = 'r',label = r'Fast Wind',capsize=5)
plt.errorbar(x,bulkrarefaction,yerr=bulkrarefaction_err,fmt='o-',color='green',label = r'Rarefaction',capsize=5)
plt.errorbar(x,bulktrailing,yerr=bulktrailing_err,fmt='o-',color='skyblue',label = r'Trailing',capsize=5)
# plt.plot(x,(bulkfast+bulkpreceding+bulktrailing+bulkshoulder+bulkrarefaction)/5,color='g',marker='s',linestyle='dashed',label = r'Mean')
plt.title('Adjusted Bulk Energy Density',fontsize=15)
plt.axvline(x=1,linestyle='dashed',color='k',linewidth=0.2)
plt.axvline(x=2,linestyle='dashed',color='k',linewidth=0.2)
plt.axvline(x=3,linestyle='dashed',color='k',linewidth=0.2)
plt.yscale('log')
plt.xticks(x,['PSP \n $\sim0.05 \  AU$','SOLO \n $\sim0.3 \  AU$','WIND \n $\sim1 \ AU$	'],fontsize=12)
# plt.legend(fontsize=12)

#%%
# Thermal Energy Density
thermfast_err = thermfast*np.sqrt((nfast_std/nfast)**2 + (Tfast_std/Tfast)**2)
thermpreceding_err = thermpreceding*np.sqrt((npreceding_std/npreceding)**2 + (Tpreceding_std/Tpreceding)**2)
thermtrailing_err = thermtrailing*np.sqrt((ntrailing_std/ntrailing)**2 + (Ttrailing_std/Ttrailing)**2)
thermshoulder_err = thermshoulder*np.sqrt((nshoulder_std/nshoulder)**2 + (Tshoulder_std/Tshoulder)**2)
thermrarefaction_err = thermrarefaction*np.sqrt((nrarefaction_std/nrarefaction)**2 + (Trarefaction_std/Trarefaction)**2)

plt.errorbar(x,thermpreceding,yerr=thermpreceding_err,fmt='o-',color='darkblue',label = r'Leading Slow Wind',capsize=5)
plt.errorbar(x,thermshoulder,yerr=thermshoulder_err,fmt='o-',color='purple',label = r'Compression Region',capsize=5)
plt.errorbar(x,thermfast,yerr=thermfast_err,fmt='o-',color = 'r',label = r'Fast CH Wind',capsize=5)
plt.errorbar(x,thermrarefaction,yerr=thermrarefaction_err,fmt='o-',color='green',label = r'Rarefaction Region',capsize=5)
plt.errorbar(x,thermtrailing,yerr=thermtrailing_err,fmt='o-',color='skyblue',label = r'Trailing Slow Wind',capsize=5)
# plt.plot(x,(thermfast+thermpreceding+thermtrailing+thermshoulder+thermrarefaction)/5,color='g',marker='s',linestyle='dashed',label = r'Mean')
plt.title('Adjusted Thermal Energy Density',fontsize=15)
plt.axvline(x=1,linestyle='dashed',color='k',linewidth=0.2)
plt.axvline(x=2,linestyle='dashed',color='k',linewidth=0.2)
plt.axvline(x=3,linestyle='dashed',color='k',linewidth=0.2)
plt.yscale('linear')
plt.xticks(x,['PSP \n $\sim0.05 \  AU$','SOLO \n $\sim0.3 \  AU$','WIND \n $\sim1 \ AU$	'],fontsize=12)
# plt.legend(fontsize=12)

# %%
# Magnetic Energy Density
magfast_err = magfast*np.sqrt(2*(Bfast_std/Bfast)**2)
magpreceding_err = magpreceding*np.sqrt(2*(Bpreceding_std/Bpreceding)**2)
magtrailing_err = magtrailing*np.sqrt(2*(Btrailing_std/Btrailing)**2)
magshoulder_err = magshoulder*np.sqrt(2*(Bshoulder_std/Bshoulder)**2)
magrarefaction_err = magrarefaction*np.sqrt(2*(Brarefaction_std/Brarefaction)**2)
plt.errorbar(x,magpreceding,yerr=magpreceding_err,fmt='o-',color='darkblue',label = r'Leading',capsize=5)
plt.errorbar(x,magshoulder,yerr=magshoulder_err,fmt='o-',color='purple',label = r'Compression',capsize=5)
plt.errorbar(x,magfast,yerr=magfast_err,fmt='o-',color = 'r',label = r'Fast',capsize=5)
plt.errorbar(x,magrarefaction,yerr=magrarefaction_err,fmt='o-',color='green',label = r'Rarefaction',capsize=5)
plt.errorbar(x,magtrailing,yerr=magtrailing_err,fmt='o-',color='skyblue',label = r'Trailing',capsize=5)
# plt.plot(x,(magfast+magpreceding+magtrailing+magshoulder+magrarefaction)/5,color='g',marker='s',linestyle='dashed',label = r'Mean')
plt.title('Adjusted Magnetic Energy Density',fontsize=15)
# plt.ylim(0,2e-12)
plt.axvline(x=1,linestyle='dashed',color='k',linewidth=0.2)
plt.axvline(x=2,linestyle='dashed',color='k',linewidth=0.2)
plt.axvline(x=3,linestyle='dashed',color='k',linewidth=0.2)
plt.yscale('linear')
plt.xticks(x,['PSP \n $\sim0.05 \  AU$','SOLO \n $\sim0.3 \  AU$','WIND \n $\sim1 \ AU$	'],fontsize=12)
# plt.legend(fontsize=20)

#%%
#Magnetic Fluctuation Energy Density
# deltamagfast_err = deltamagfast*np.sqrt(2*(dBfast_std/dBfast)**2)
# deltamagpreceding_err = deltamagpreceding*np.sqrt(2*(dBpreceding_std/dBpreceding)**2)
# deltamagtrailing_err = deltamagtrailing*np.sqrt(2*(dBtrailing_std/dBtrailing)**2)
# deltamagshoulder_err = deltamagshoulder*np.sqrt(2*(dBshoulder_std/dBshoulder)**2)
# deltamagrarefaction_err = deltamagrarefaction*np.sqrt(2*(dBrarefaction_std/dBrarefaction)**2)
# plt.errorbar(x,deltamagfast,yerr=deltamagfast_err,fmt='o-',color = 'r',label = r'Fast Wind',capsize=5)
# plt.errorbar(x,deltamagpreceding,yerr=deltamagpreceding_err,fmt='o-',color='skyblue',label = r'Preceding',capsize=5)
# plt.errorbar(x,deltamagtrailing,yerr=deltamagtrailing_err,fmt='o-',color='darkblue',label = r'Trailing',capsize=5)
# plt.errorbar(x,deltamagshoulder,yerr=deltamagshoulder_err,fmt='o-',color='purple',label = r'Shoulder',capsize=5)
# plt.errorbar(x,deltamagrarefaction,yerr=deltamagrarefaction_err,fmt='o-',color='green',label = r'Rarefaction',capsize=5)
# # plt.plot(x,(deltamagfast+deltamagpreceding+deltamagtrailing+deltamagshoulder+deltamagrarefaction)/5,color='g',marker='s',linestyle='dashed',label = r'Mean')
# plt.title('Magnetic Fluctuation Energy Density [J/m^3]')
# plt.axvline(x=1,linestyle='dashed',color='k',linewidth=0.2)
# plt.axvline(x=2,linestyle='dashed',color='k',linewidth=0.2)
# plt.axvline(x=3,linestyle='dashed',color='k',linewidth=0.2)
# plt.xticks(x,['PSP','SOLO','WIND'])
# plt.xlabel('Spacecraft')
# plt.legend(fontsize=15)
# %%
