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
     z_label='', title='',nbins=50,mincount=10,vmin=None,vmax=None, mask=None,
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

    # auto-scale color range from the data itself if not explicitly given
    if vmin is None or vmax is None:
        hist_for_range = hist_masked.compressed()  # ignores masked (low-count) bins
        if hist_for_range.size == 0:
            hist_for_range = z_valid  # fallback if everything got masked out
        if vmin is None:
            vmin = np.nanpercentile(hist_for_range, 2)
        if vmax is None:
            vmax = np.nanpercentile(hist_for_range, 98)

    fig, ax = plt.subplots(figsize=(10, 8))

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
    ax.plot(beta_par_range, mirror_threshold, label='Mirror', color='k',linestyle='dotted', linewidth=2)
    ax.plot(beta_par_range, firehose_threshold, label='Oblique FH', color='k', linestyle='--', linewidth=2)
    ax.plot(beta_par_range, cyclotron_threshold, label='IC', color='k',linestyle='dotted', linewidth=1)
    ax.plot(beta_par_range, parfirehose_threshold, label='Parallel FH', color='k', linestyle='--', linewidth=1)
    ax.legend()
    plt.show()
    return fig, ax







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
                  cgl_line=True, cgl_exponent=-1,
                  meanfit=False):

    # From Hellinger (2011)
    # mirror_params = [0.77, 0.76, -0.016]  # [a, b, beta0]
    # firehose_params = [-1.4, 1.0, -0.11]  # [a, b, beta0]
    # cyclotron_params = [0.43, 0.42, -0.0004]  # [a, b, beta0]
    # parfirehose_params = [-0.47, 0.53, 0.59]  # [a, b, beta0]

    mirror_params = [0.77, 0.76, -0.016]  # [a, b, beta0]
    firehose_params = [-1.4, 1.0, -0.11]  # [a, b, beta0]
    cyclotron_params = [0.43, 0.42, -0.0004]  # [a, b, beta0]
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

    mean_points = []  # collects (x_mean, y_mean) for each spacecraft actually plotted

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
            mean_points.append((x_mean, y_mean))

            if cgl_line:
                yline_cgl = y_mean * (xline_cgl / x_mean) ** cgl_exponent
                ax.plot(xline_cgl, yline_cgl, color=color, linestyle=':', linewidth=1.5, alpha=0.8, zorder=9)

        elif plot_type == 'means_only':
            # just the mean marker, no contour or scatter cloud
            x_mean = np.mean(x_valid)
            y_mean = np.mean(y_valid)
            ax.scatter(x_mean, y_mean, marker='^', color=color, s=4*mean_marker_size,
                       edgecolors='k', linewidths=1.5, zorder=10, label=sc.upper())
            mean_points.append((x_mean, y_mean))

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
                mean_points.append((x_mean, y_mean))

                if cgl_line:
                    yline_cgl = y_mean * (xline_cgl / x_mean) ** cgl_exponent
                    ax.plot(xline_cgl, yline_cgl, color=color, linestyle='-', linewidth=1, alpha=0.8, zorder=9)

    if meanfit and len(mean_points) >= 2:
        mean_points_arr = np.array(mean_points)
        mx_fit = np.log10(mean_points_arr[:, 0])
        my_fit = np.log10(mean_points_arr[:, 1])
        m_mean, c_mean = np.polyfit(mx_fit, my_fit, 1)  # equally weighted, one point per spacecraft

        xline_mean = np.logspace(xrange[0], xrange[1], 100)
        yline_mean = 10**(m_mean * np.log10(xline_mean) + c_mean)
        ax.plot(xline_mean, yline_mean, color='blue', linestyle='-', linewidth=1.5, zorder=11,
                label=rf'Means fit: $\beta_\parallel^{{{m_mean:.2f}}}$')

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
    ax.legend(fontsize=15,framealpha=1,loc='upper right')
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



# brazilhist(1.602e19*abs(allT), z_label = r'$T \ (eV)$',fitline=False,nbins=200,vmin=00,vmax=200,count=False,scale='linear')
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
                        mean_marker_size=35,cgl_exponent=-0.55,cgl_line=True,meanfit=False)

#%%
fig, ax = brazil_multi(betapar_psp=betaparspsp_trailing, Tparperp_psp=Tparperppsp_trailing,
                        betapar_solo=betaparssolo_trailing, Tparperp_solo=Tparperpsolo_trailing,
                        betapar_wind=betaparswind_trailing, Tparperp_wind=Tparperpwind_trailing,
                        contour_nlevels=10,title = "Trailing Slow Wind",plot_type='contour',
                        mean_marker_size=35,cgl_line=False,meanfit=True)
#%%
fig, ax = brazil_multi(betapar_psp=betaparspsp_shoulder, Tparperp_psp=Tparperppsp_shoulder,
                        betapar_solo=betaparssolo_shoulder, Tparperp_solo=Tparperpsolo_shoulder,
                        betapar_wind=betaparswind_shoulder, Tparperp_wind=Tparperpwind_shoulder,
                        contour_nlevels=10,title = "Compression Region",plot_type='contour',
                        mean_marker_size=35,cgl_line=False)

#%%
fig, ax = brazil_multi(betapar_psp=betaparspsp_fast, Tparperp_psp=Tparperppsp_fast,
                        betapar_solo=betaparssolo_fast, Tparperp_solo=Tparperpsolo_fast,
                        betapar_wind=betaparswind_fast, Tparperp_wind=Tparperpwind_fast,
                        contour_nlevels=10,title = "Fast CH Wind",plot_type='contour',
                        mean_marker_size=35,fitline=False,cgl_exponent=-1,meanfit=True)

#%%
fig, ax = brazil_multi(betapar_psp=betaparspsp_rarefaction, Tparperp_psp=Tparperppsp_rarefaction,
                        betapar_solo=betaparssolo_rarefaction, Tparperp_solo=Tparperpsolo_rarefaction,
                        betapar_wind=betaparswind_rarefaction, Tparperp_wind=Tparperpwind_rarefaction,
                        contour_nlevels=10,title = "Rarefaction Region",plot_type='contour',
                        mean_marker_size=35,cgl_line=False,cgl_exponent=-1,meanfit=False)

#%%

x = np.array([1,2,3])
def get_vars(type = 'fast'):
    # rarray = np.array([0.05/0.45,0.3/0.45,1/0.45])
    rarray = np.array([0.05/0.05,0.3/0.05,1/0.05])
    # rarray = np.ones(3)
    # rarray = np.array([1,6.3,21.2])

    if type == 'fast':
        narray = 1e6*np.array([npsp_fast,nsolo_fast,nwind_fast])
        n_stdarray = 1e6*np.array([nstdpsp_fast,nstdsolo_fast,nstdwind_fast])
        n_pctarray = 1e6*np.array([npctpsp_fast,npctsolo_fast,npctwind_fast])
        Tarray = 1.602e-19*np.array([Tpsp_fast,Tsolo_fast,Twind_fast])
        T_stdarray = 1.602e-19*np.array([Tstdpsp_fast,Tstdsolo_fast,Tstdwind_fast])
        T_pctarray = 1.602e-19*np.array([Tpctpsp_fast,Tpctsolo_fast,Tpctwind_fast])
        varray = 1e3*np.array([vpsp_fast,vsolo_fast,vwind_fast])
        v_stdarray = 1e3*np.array([vstdpsp_fast,vstdsolo_fast,vstdwind_fast])
        v_pctarray = 1e3*np.array([vpctpsp_fast,vpctsolo_fast,vpctwind_fast])
        Barray = 1e-9*np.array([Bpsp_fast,Bsolo_fast,Bwind_fast])
        B_stdarray = 1e-9*np.array([Bstdpsp_fast,Bstdsolo_fast,Bstdwind_fast])
        B_pctarray = 1e-9*np.array([Bpctpsp_fast,Bpctsolo_fast,Bpctwind_fast])
        dB_array = 1e-9*np.array([dBpsp_fast,dBsolo_fast,dBwind_fast])
        dB_stdarray = 1e-9*np.array([dBstdpsp_fast,dBstdsolo_fast,dBstdwind_fast])
        dB_pctarray = 1e-9*np.array([dBpctpsp_fast,dBpctsolo_fast,dBpctwind_fast])
        charray = np.array([chpsp_fast,chsolo_fast,chwind_fast])
        ch_stdarray = np.array([chstdpsp_fast,chstdsolo_fast,chstdwind_fast])
        ch_pctarray = np.array([chpctpsp_fast,chpctsolo_fast,chpctwind_fast])
        rearray = np.array([repsp_fast,resolo_fast,rewind_fast])
        re_stdarray = np.array([restdpsp_fast,restdsolo_fast,restdwind_fast])
        re_pctarray = np.array([repctpsp_fast,repctsolo_fast,repctwind_fast])
        Cpararray = np.array([Cparpsp_fast,Cparsolo_fast,Cparwind_fast])
        Cpar_stdarray = np.array([Cparstdpsp_fast,Cparstdsolo_fast,Cparstdwind_fast])
        Cpar_pctarray = np.array([Cparpctpsp_fast,Cparpctsolo_fast,Cparpctwind_fast])
        Cperparray = np.array([Cperppsp_fast,Cperpsolo_fast,Cperpwind_fast])
        Cperp_stdarray = np.array([Cperpstdpsp_fast,Cperpstdsolo_fast,Cperpstdwind_fast])
        Cperp_pctarray = np.array([Cperppctpsp_fast,Cperppctsolo_fast,Cperppctwind_fast])
    elif type == 'slow':
        narray = 1e6*np.array([npsp_slow,nsolo_slow,nwind_slow])
        n_stdarray = 1e6*np.array([nstdpsp_slow,nstdsolo_slow,nstdwind_slow])
        n_pctarray = 1e6*np.array([npctpsp_slow,npctsolo_slow,npctwind_slow])
        Tarray = 1.602e-19*np.array([Tpsp_slow,Tsolo_slow,Twind_slow])
        T_stdarray = 1.602e-19*np.array([Tstdpsp_slow,Tstdsolo_slow,Tstdwind_slow])
        T_pctarray = 1.602e-19*np.array([Tpctpsp_slow,Tpctsolo_slow,Tpctwind_slow])
        varray = 1e3*np.array([vpsp_slow,vsolo_slow,vwind_slow])
        v_stdarray = 1e3*np.array([vstdpsp_slow,vstdsolo_slow,vstdwind_slow])
        v_pctarray = 1e3*np.array([vpctpsp_slow,vpctsolo_slow,vpctwind_slow])
        Barray = 1e-9*np.array([Bpsp_slow,Bsolo_slow,Bwind_slow])
        B_stdarray = 1e-9*np.array([Bstdpsp_slow,Bstdsolo_slow,Bstdwind_slow])
        B_pctarray = 1e-9*np.array([Bpctpsp_slow,Bpctsolo_slow,Bpctwind_slow])
        dB_array = 1e-9*np.array([dBpsp_slow,dBsolo_slow,dBwind_slow])
        dB_stdarray = 1e-9*np.array([dBstdpsp_slow,dBstdsolo_slow,dBstdwind_slow])
        dB_pctarray = 1e-9*np.array([dBpctpsp_slow,dBpctsolo_slow,dBpctwind_slow])
        charray = np.array([chpsp_slow,chsolo_slow,chwind_slow])
        ch_stdarray = np.array([chstdpsp_slow,chstdsolo_slow,chstdwind_slow])
        ch_pctarray = np.array([chpctpsp_slow,chpctsolo_slow,chpctwind_slow])
        rearray = np.array([repsp_slow,resolo_slow,rewind_slow])
        re_stdarray = np.array([restdpsp_slow,restdsolo_slow,restdwind_slow])
        re_pctarray = np.array([repctpsp_slow,repctsolo_slow,repctwind_slow])
        Cpararray = np.array([Cparpsp_slow,Cparsolo_slow,Cparwind_slow])
        Cpar_stdarray = np.array([Cparstdpsp_slow,Cparstdsolo_slow,Cparstdwind_slow])
        Cpar_pctarray = np.array([Cparpctpsp_slow,Cparpctsolo_slow,Cparpctwind_slow])
        Cperparray = np.array([Cperppsp_slow,Cperpsolo_slow,Cperpwind_slow])
        Cperp_stdarray = np.array([Cperpstdpsp_slow,Cperpstdsolo_slow,Cperpstdwind_slow])
        Cperp_pctarray = np.array([Cperppctpsp_slow,Cperppctsolo_slow,Cperppctwind_slow])
    elif type == 'preceding':
        narray = 1e6*np.array([npsp_preceding,nsolo_preceding,nwind_preceding])
        n_stdarray = 1e6*np.array([nstdpsp_preceding,nstdsolo_preceding,nstdwind_preceding])
        n_pctarray = 1e6*np.array([npctpsp_preceding,npctsolo_preceding,npctwind_preceding])
        Tarray = 1.602e-19*np.array([Tpsp_preceding,Tsolo_preceding,Twind_preceding])
        T_stdarray = 1.602e-19*np.array([Tstdpsp_preceding,Tstdsolo_preceding,Tstdwind_preceding])
        T_pctarray = 1.602e-19*np.array([Tpctpsp_preceding,Tpctsolo_preceding,Tpctwind_preceding])
        varray = 1e3*np.array([vpsp_preceding,vsolo_preceding,vwind_preceding])
        v_stdarray = 1e3*np.array([vstdpsp_preceding,vstdsolo_preceding,vstdwind_preceding])
        v_pctarray = 1e3*np.array([vpctpsp_preceding,vpctsolo_preceding,vpctwind_preceding])
        Barray = 1e-9*np.array([Bpsp_preceding,Bsolo_preceding,Bwind_preceding])
        B_stdarray = 1e-9*np.array([Bstdpsp_preceding,Bstdsolo_preceding,Bstdwind_preceding])
        B_pctarray = 1e-9*np.array([Bpctpsp_preceding,Bpctsolo_preceding,Bpctwind_preceding])
        dB_array = 1e-9*np.array([dBpsp_preceding,dBsolo_preceding,dBwind_preceding])
        dB_stdarray = 1e-9*np.array([dBstdpsp_preceding,dBstdsolo_preceding,dBstdwind_preceding])
        dB_pctarray = 1e-9*np.array([dBpctpsp_preceding,dBpctsolo_preceding,dBpctwind_preceding])
        charray = np.array([chpsp_preceding,chsolo_preceding,chwind_preceding])
        ch_stdarray = np.array([chstdpsp_preceding,chstdsolo_preceding,chstdwind_preceding])
        ch_pctarray = np.array([chpctpsp_preceding,chpctsolo_preceding,chpctwind_preceding])
        rearray = np.array([repsp_preceding,resolo_preceding,rewind_preceding])
        re_stdarray = np.array([restdpsp_preceding,restdsolo_preceding,restdwind_preceding])
        re_pctarray = np.array([repctpsp_preceding,repctsolo_preceding,repctwind_preceding])
        Cpararray = np.array([Cparpsp_preceding,Cparsolo_preceding,Cparwind_preceding])
        Cpar_stdarray = np.array([Cparstdpsp_preceding,Cparstdsolo_preceding,Cparstdwind_preceding])
        Cpar_pctarray = np.array([Cparpctpsp_preceding,Cparpctsolo_preceding,Cparpctwind_preceding])
        Cperparray = np.array([Cperppsp_preceding,Cperpsolo_preceding,Cperpwind_preceding])
        Cperp_stdarray = np.array([Cperpstdpsp_preceding,Cperpstdsolo_preceding,Cperpstdwind_preceding])
        Cperp_pctarray = np.array([Cperppctpsp_preceding,Cperppctsolo_preceding,Cperppctwind_preceding])
    elif type == 'trailing':
        narray = 1e6*np.array([npsp_trailing,nsolo_trailing,nwind_trailing])
        n_stdarray = 1e6*np.array([nstdpsp_trailing,nstdsolo_trailing,nstdwind_trailing])
        n_pctarray = 1e6*np.array([npctpsp_trailing,npctsolo_trailing,npctwind_trailing])
        Tarray = 1.602e-19*np.array([Tpsp_trailing,Tsolo_trailing,Twind_trailing])
        T_stdarray = 1.602e-19*np.array([Tstdpsp_trailing,Tstdsolo_trailing,Tstdwind_trailing])
        T_pctarray = 1.602e-19*np.array([Tpctpsp_trailing,Tpctsolo_trailing,Tpctwind_trailing])
        varray = 1e3*np.array([vpsp_trailing,vsolo_trailing,vwind_trailing])
        v_stdarray = 1e3*np.array([vstdpsp_trailing,vstdsolo_trailing,vstdwind_trailing])
        v_pctarray = 1e3*np.array([vpctpsp_trailing,vpctsolo_trailing,vpctwind_trailing])
        Barray = 1e-9*np.array([Bpsp_trailing,Bsolo_trailing,Bwind_trailing])
        B_stdarray = 1e-9*np.array([Bstdpsp_trailing,Bstdsolo_trailing,Bstdwind_trailing])
        B_pctarray = 1e-9*np.array([Bpctpsp_trailing,Bpctsolo_trailing,Bpctwind_trailing])
        dB_array = 1e-9*np.array([dBpsp_trailing,dBsolo_trailing,dBwind_trailing])
        dB_stdarray = 1e-9*np.array([dBstdpsp_trailing,dBstdsolo_trailing,dBstdwind_trailing])
        dB_pctarray = 1e-9*np.array([dBpctpsp_trailing,dBpctsolo_trailing,dBpctwind_trailing])
        charray = np.array([chpsp_trailing,chsolo_trailing,chwind_trailing])
        ch_stdarray = np.array([chstdpsp_trailing,chstdsolo_trailing,chstdwind_trailing])
        ch_pctarray = np.array([chpctpsp_trailing,chpctsolo_trailing,chpctwind_trailing])
        rearray = np.array([repsp_trailing,resolo_trailing,rewind_trailing])
        re_stdarray = np.array([restdpsp_trailing,restdsolo_trailing,restdwind_trailing])
        re_pctarray = np.array([repctpsp_trailing,repctsolo_trailing,repctwind_trailing])
        Cpararray = np.array([Cparpsp_trailing,Cparsolo_trailing,Cparwind_trailing])
        Cpar_stdarray = np.array([Cparstdpsp_trailing,Cparstdsolo_trailing,Cparstdwind_trailing])
        Cpar_pctarray = np.array([Cparpctpsp_trailing,Cparpctsolo_trailing,Cparpctwind_trailing])
        Cperparray = np.array([Cperppsp_trailing,Cperpsolo_trailing,Cperpwind_trailing])
        Cperp_stdarray = np.array([Cperpstdpsp_trailing,Cperpstdsolo_trailing,Cperpstdwind_trailing])
        Cperp_pctarray = np.array([Cperppctpsp_trailing,Cperppctsolo_trailing,Cperppctwind_trailing])
    elif type == 'shoulder':
        narray = 1e6*np.array([npsp_shoulder,nsolo_shoulder,nwind_shoulder])
        n_stdarray = 1e6*np.array([nstdpsp_shoulder,nstdsolo_shoulder,nstdwind_shoulder])
        n_pctarray = 1e6*np.array([npctpsp_shoulder,npctsolo_shoulder,npctwind_shoulder])
        Tarray = 1.602e-19*np.array([Tpsp_shoulder,Tsolo_shoulder,Twind_shoulder])
        T_stdarray = 1.602e-19*np.array([Tstdpsp_shoulder,Tstdsolo_shoulder,Tstdwind_shoulder])
        T_pctarray = 1.602e-19*np.array([Tpctpsp_shoulder,Tpctsolo_shoulder,Tpctwind_shoulder])
        varray = 1e3*np.array([vpsp_shoulder,vsolo_shoulder,vwind_shoulder])
        v_stdarray = 1e3*np.array([vstdpsp_shoulder,vstdsolo_shoulder,vstdwind_shoulder])
        v_pctarray = 1e3*np.array([vpctpsp_shoulder,vpctsolo_shoulder,vpctwind_shoulder])
        Barray = 1e-9*np.array([Bpsp_shoulder,Bsolo_shoulder,Bwind_shoulder])
        B_stdarray = 1e-9*np.array([Bstdpsp_shoulder,Bstdsolo_shoulder,Bstdwind_shoulder])
        B_pctarray = 1e-9*np.array([Bpctpsp_shoulder,Bpctsolo_shoulder,Bpctwind_shoulder])
        dB_array = 1e-9*np.array([dBpsp_shoulder,dBsolo_shoulder,dBwind_shoulder])
        dB_stdarray = 1e-9*np.array([dBstdpsp_shoulder,dBstdsolo_shoulder,dBstdwind_shoulder])
        dB_pctarray = 1e-9*np.array([dBpctpsp_shoulder,dBpctsolo_shoulder,dBpctwind_shoulder])
        charray = np.array([chpsp_shoulder,chsolo_shoulder,chwind_shoulder])
        ch_stdarray = np.array([chstdpsp_shoulder,chstdsolo_shoulder,chstdwind_shoulder])
        ch_pctarray = np.array([chpctpsp_shoulder,chpctsolo_shoulder,chpctwind_shoulder])
        rearray = np.array([repsp_shoulder,resolo_shoulder,rewind_shoulder])
        re_stdarray = np.array([restdpsp_shoulder,restdsolo_shoulder,restdwind_shoulder])
        re_pctarray = np.array([repctpsp_shoulder,repctsolo_shoulder,repctwind_shoulder])
        Cpararray = np.array([Cparpsp_shoulder,Cparsolo_shoulder,Cparwind_shoulder])
        Cpar_stdarray = np.array([Cparstdpsp_shoulder,Cparstdsolo_shoulder,Cparstdwind_shoulder])
        Cpar_pctarray = np.array([Cparpctpsp_shoulder,Cparpctsolo_shoulder,Cparpctwind_shoulder])
        Cperparray = np.array([Cperppsp_shoulder,Cperpsolo_shoulder,Cperpwind_shoulder])
        Cperp_stdarray = np.array([Cperpstdpsp_shoulder,Cperpstdsolo_shoulder,Cperpstdwind_shoulder])
        Cperp_pctarray = np.array([Cperppctpsp_shoulder,Cperppctsolo_shoulder,Cperppctwind_shoulder])
    elif type == 'rarefaction':
        narray = 1e6*np.array([npsp_rarefaction,nsolo_rarefaction,nwind_rarefaction])
        n_stdarray = 1e6*np.array([nstdpsp_rarefaction,nstdsolo_rarefaction,nstdwind_rarefaction])
        n_pctarray = 1e6*np.array([npctpsp_rarefaction,npctsolo_rarefaction,npctwind_rarefaction])
        Tarray = 1.602e-19*np.array([Tpsp_rarefaction,Tsolo_rarefaction,Twind_rarefaction])
        T_stdarray = 1.602e-19*np.array([Tstdpsp_rarefaction,Tstdsolo_rarefaction,Tstdwind_rarefaction])
        T_pctarray = 1.602e-19*np.array([Tpctpsp_rarefaction,Tpctsolo_rarefaction,Tpctwind_rarefaction])
        varray = 1e3*np.array([vpsp_rarefaction,vsolo_rarefaction,vwind_rarefaction])
        v_stdarray = 1e3*np.array([vstdpsp_rarefaction,vstdsolo_rarefaction,vstdwind_rarefaction])
        v_pctarray = 1e3*np.array([vpctpsp_rarefaction,vpctsolo_rarefaction,vpctwind_rarefaction])
        Barray = 1e-9*np.array([Bpsp_rarefaction,Bsolo_rarefaction,Bwind_rarefaction])
        B_stdarray = 1e-9*np.array([Bstdpsp_rarefaction,Bstdsolo_rarefaction,Bstdwind_rarefaction])
        B_pctarray = 1e-9*np.array([Bpctpsp_rarefaction,Bpctsolo_rarefaction,Bpctwind_rarefaction])
        dB_array = 1e-9*np.array([dBpsp_rarefaction,dBsolo_rarefaction,dBwind_rarefaction])
        dB_stdarray = 1e-9*np.array([dBstdpsp_rarefaction,dBstdsolo_rarefaction,dBstdwind_rarefaction])
        dB_pctarray = 1e-9*np.array([dBpctpsp_rarefaction,dBpctsolo_rarefaction,dBpctwind_rarefaction])
        charray = np.array([chpsp_rarefaction,chsolo_rarefaction,chwind_rarefaction])
        ch_stdarray = np.array([chstdpsp_rarefaction,chstdsolo_rarefaction,chstdwind_rarefaction])
        ch_pctarray = np.array([chpctpsp_rarefaction,chpctsolo_rarefaction,chpctwind_rarefaction])
        rearray = np.array([repsp_rarefaction,resolo_rarefaction,rewind_rarefaction])
        re_stdarray = np.array([restdpsp_rarefaction,restdsolo_rarefaction,restdwind_rarefaction])
        re_pctarray = np.array([repctpsp_rarefaction,repctsolo_rarefaction,repctwind_rarefaction])
        Cpararray = np.array([Cparpsp_rarefaction,Cparsolo_rarefaction,Cparwind_rarefaction])
        Cpar_stdarray = np.array([Cparstdpsp_rarefaction,Cparstdsolo_rarefaction,Cparstdwind_rarefaction])
        Cpar_pctarray = np.array([Cparpctpsp_rarefaction,Cparpctsolo_rarefaction,Cparpctwind_rarefaction])
        Cperparray = np.array([Cperppsp_rarefaction,Cperpsolo_rarefaction,Cperpwind_rarefaction])
        Cperp_stdarray = np.array([Cperpstdpsp_rarefaction,Cperpstdsolo_rarefaction,Cperpstdwind_rarefaction])
        Cperp_pctarray = np.array([Cperppctpsp_rarefaction,Cperppctsolo_rarefaction,Cperppctwind_rarefaction])
    elif type == 'ir':
        narray = 1e6*np.array([npsp,nsolo,nwind])
        n_stdarray = 1e6*np.array([nstdpsp,nstdsolo,nstdwind])
        n_pctarray = 1e6*np.array([npctpsp,npctsolo,npctwind])
        Tarray = 1.602e-19*np.array([Tpsp,Tsolo,Twind])
        T_stdarray = 1.602e-19*np.array([Tstdpsp,Tstdsolo,Tstdwind])
        T_pctarray = 1.602e-19*np.array([Tpctpsp,Tpctsolo,Tpctwind])
        varray = 1e3*np.array([vpsp,vsolo,vwind])
        v_stdarray = 1e3*np.array([vstdpsp,vstdsolo,vstdwind])
        v_pctarray = 1e3*np.array([vpctpsp,vpctsolo,vpctwind])
        Barray = 1e-9*np.array([Bpsp,Bsolo,Bwind])
        B_stdarray = 1e-9*np.array([Bstdpsp,Bstdsolo,Bstdwind])
        B_pctarray = 1e-9*np.array([Bpctpsp,Bpctsolo,Bpctwind])
        dB_array = 1e-9*np.array([dBpsp,dBsolo,dBwind])
        dB_stdarray = 1e-9*np.array([dBstdpsp,dBstdsolo,dBstdwind])
        dB_pctarray = 1e-9*np.array([dBpctpsp,dBpctsolo,dBpctwind])
        charray = np.array([chpsp,chsolo,chwind])
        ch_stdarray = np.array([chstdpsp,chstdsolo,chstdwind])
        ch_pctarray = np.array([chpctpsp,chpctsolo,chpctwind])
        rearray = np.array([repsp,resolo,rewind])
        re_stdarray = np.array([restdpsp,restdsolo,restdwind])
        re_pctarray = np.array([repctpsp,repctsolo,repctwind])
        Cpararray = np.array([Cparpsp,Cparsolo,Cparwind])
        Cpar_stdarray = np.array([Cparstdpsp,Cparstdsolo,Cparstdwind])
        Cpar_pctarray = np.array([Cparpctpsp,Cparpctsolo,Cparpctwind])
        Cperparray = np.array([Cperppsp,Cperpsolo,Cperpwind])
        Cperp_stdarray = np.array([Cperpstdpsp,Cperpstdsolo,Cperpstdwind])
        Cperp_pctarray = np.array([Cperppctpsp,Cperppctsolo,Cperppctwind])
    else:
        print('No var type selected')
    
  
    narray_norm = (rarray**2)*narray
    n_stdarray_norm = (rarray**2)*n_stdarray
    n_pctarray_norm = (rarray[:,None]**2)*n_pctarray
    Tarray_norm = (rarray**(4/3))*Tarray
    T_stdarray_norm = (rarray**(4/3))*T_stdarray
    T_pctarray_norm = (rarray[:,None]**(4/3))*T_pctarray
    varray_norm = varray
    v_stdarray_norm = v_stdarray
    v_pctarray_norm = v_pctarray
    # Barray_norm = (rarray**2)*Barray
    Barray_norm = (rarray**1)*Barray
    # B_stdarray_norm = (rarray**2)*B_stdarray
    B_stdarray_norm = (rarray**1)*B_stdarray
    # B_pctarray_norm = (rarray[:,None]**2)*B_pctarray
    B_pctarray_norm = (rarray[:,None]**1)*B_pctarray
    dB_array_norm = (rarray**2)*dB_array
    dB_stdarray_norm = (rarray**2)*dB_stdarray
    dB_pctarray_norm = (rarray[:,None]**2)*dB_pctarray
    # ch, re, Cpar, Cperp (and their std/pct) are left unscaled, matching their mean arrays
    
    # Temporary
    # Barray_norm = (rarray**1)*Barray
    # B_stdarray_norm = (rarray**1)*B_stdarray
    # B_pctarray_norm = (rarray[:,None]**1)*B_pctarray
    # dB_array_norm = (rarray**1)*dB_array
    # dB_stdarray_norm = (rarray**1)*dB_stdarray
    # dB_pctarray_norm = (rarray[:,None]**1)*dB_pctarray
    
    return (narray_norm,Tarray_norm,varray_norm,Barray_norm,dB_array_norm,charray,rearray,Cpararray,Cperparray,
            n_stdarray_norm, T_stdarray_norm, v_stdarray_norm, B_stdarray_norm, dB_stdarray_norm, ch_stdarray, re_stdarray, Cpar_stdarray, Cperp_stdarray,
            n_pctarray_norm, T_pctarray_norm, v_pctarray_norm, B_pctarray_norm, dB_pctarray_norm, ch_pctarray, re_pctarray, Cpar_pctarray, Cperp_pctarray)



(nfast,Tfast,vfast,Bfast,dBfast,chfast,refast,Cparfast,Cperpfast,
 nfast_std,Tfast_std,vfast_std,Bfast_std,dBfast_std,chfast_std,refast_std,Cparfast_std,Cperpfast_std,
 nfast_pct,Tfast_pct,vfast_pct,Bfast_pct,dBfast_pct,chfast_pct,refast_pct,Cparfast_pct,Cperpfast_pct) = get_vars('fast')

(npreceding,Tpreceding,vpreceding,Bpreceding,dBpreceding,chpreceding,repreceding,Cparpreceding,Cperppreceding,
 npreceding_std,Tpreceding_std,vpreceding_std,Bpreceding_std,dBpreceding_std,chpreceding_std,repreceding_std,Cparpreceding_std,Cperppreceding_std,
 npreceding_pct,Tpreceding_pct,vpreceding_pct,Bpreceding_pct,dBpreceding_pct,chpreceding_pct,repreceding_pct,Cparpreceding_pct,Cperppreceding_pct) = get_vars('preceding')

(ntrailing,Ttrailing,vtrailing,Btrailing,dBtrailing,chtrailing,retrailing,Cpartrailing,Cperptrailing,
 ntrailing_std,Ttrailing_std,vtrailing_std,Btrailing_std,dBtrailing_std,chtrailing_std,retrailing_std,Cpartrailing_std,Cperptrailing_std,
 ntrailing_pct,Ttrailing_pct,vtrailing_pct,Btrailing_pct,dBtrailing_pct,chtrailing_pct,retrailing_pct,Cpartrailing_pct,Cperptrailing_pct) = get_vars('trailing')

(nshoulder,Tshoulder,vshoulder,Bshoulder,dBshoulder,chshoulder,reshoulder,Cparshoulder,Cperpshoulder,
 nshoulder_std,Tshoulder_std,vshoulder_std,Bshoulder_std,dBshoulder_std,chshoulder_std,reshoulder_std,Cparshoulder_std,Cperpshoulder_std,
 nshoulder_pct,Tshoulder_pct,vshoulder_pct,Bshoulder_pct,dBshoulder_pct,chshoulder_pct,reshoulder_pct,Cparshoulder_pct,Cperpshoulder_pct) = get_vars('shoulder')

(nrarefaction,Trarefaction,vrarefaction,Brarefaction,dBrarefaction,chrarefaction,rerarefaction,Cparrarefaction,Cperprarefaction,
 nrarefaction_std,Trarefaction_std,vrarefaction_std,Brarefaction_std,dBrarefaction_std,chrarefaction_std,rerarefaction_std,Cparrarefaction_std,Cperprarefaction_std,
 nrarefaction_pct,Trarefaction_pct,vrarefaction_pct,Brarefaction_pct,dBrarefaction_pct,chrarefaction_pct,rerarefaction_pct,Cparrarefaction_pct,Cperprarefaction_pct) = get_vars('rarefaction')




#%%
## ALL THE NEWER VERSIONS USING MEDIANS AND PCTS

def pct_to_yerr(mean_array, pct_array):
    """
    mean_array: shape (N,) central values
    pct_array:  shape (N,2) [lower,upper] bounds
    returns:    shape (2,N) [lower_err, upper_err] for yerr=
    """
    lower_err = mean_array - pct_array[:,0]
    upper_err = pct_array[:,1] - mean_array
    return np.array([lower_err, upper_err])


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
#%%

import numpy as np
import matplotlib.pyplot as plt


xr = np.array([0.05, 0.3, 1.0])
S = 1e-15  # one scale factor for both points and errorbars

def flux_err(n, n_pct, v, v_pct):
    n, v = np.asarray(n, float), np.asarray(v, float)
    en = np.asarray(pct_to_yerr(n, n_pct))   # (2, N)
    ev = np.asarray(pct_to_yerr(v, v_pct))   # (2, N)
    f = n * v
    return f * S, f * np.sqrt((en / n)**2 + (ev / v)**2) * S

rois = [('Leading',     npreceding,   npreceding_pct,   vpreceding,   vpreceding_pct,   'darkblue'),
        ('Compression', nshoulder,    nshoulder_pct,    vshoulder,    vshoulder_pct,    'purple'),
        ('Fast',        nfast,        nfast_pct,        vfast,        vfast_pct,        'r'),
        ('Rarefaction', nrarefaction, nrarefaction_pct, vrarefaction, vrarefaction_pct, 'green'),
        ('Trailing',    ntrailing,    ntrailing_pct,    vtrailing,    vtrailing_pct,    'skyblue')]

fs, errs = [], []
for name, n, npct, v, vpct, c in rois:
    f, yerr = flux_err(n, npct, v, vpct)
    fs.append(f)
    errs.append(yerr)
    # faded individual ROIs
    plt.errorbar(xr, f, yerr=yerr, fmt='o-', color=c, capsize=5,
                 alpha=0.3, label=rf"{name}")

fs   = np.array(fs)       # (5, 3)
errs = np.array(errs)     # (5, 2, 3)

# average of the ROI medians, with propagated (asymmetric) error
fmean    = fs.mean(axis=0)                                  # (3,)
mean_err = np.sqrt((errs**2).sum(axis=0)) / len(rois)       # (2, 3)

plt.errorbar(xr, fmean, yerr=mean_err, fmt='o-', color='k', capsize=6,
             lw=2.5, ms=5, elinewidth=2, zorder=10,
             label='avg of ROI medians')

for r in xr: plt.axvline(x=r, linestyle='dashed', color='k', linewidth=0.2)
plt.ylabel(r"$n'v_R \ (\times 10^{15} m^{-2} s^{-1})$", fontsize=15)
plt.xlim(0, 1.1)
plt.yticks(fontsize=12)
plt.xticks(xr, ['PSP \n $\\sim0.05$ AU', 'SOLO \n $\\sim0.3$ AU', 'WIND \n $\\sim1$ AU'], fontsize=12)
plt.legend()

#%%
import numpy as np
import matplotlib.pyplot as plt

xr = np.array([0.05, 0.3, 1.0])
f = 1e9  # T -> nT

rois = [('leading',     Bpreceding,   Bpreceding_pct,   'darkblue'),
        ('compression', Bshoulder,    Bshoulder_pct,    'purple'),
        ('fast',        Bfast,        Bfast_pct,        'r'),
        ('rarefaction', Brarefaction, Brarefaction_pct, 'green'),
        ('trailing',    Btrailing,    Btrailing_pct,    'skyblue')]

Bs, errs = [], []
for name, B, Bpct, c in rois:
    B = np.asarray(B, float)
    yerr = f * np.asarray(pct_to_yerr(B, Bpct))     # (2, 3)
    Bs.append(f * B)
    errs.append(yerr)
    # faded individual ROIs
    plt.errorbar(xr, f * B, yerr=yerr, fmt='o-', color=c, capsize=5,
                 alpha=0.3, label=rf'${{name}}$')

Bs   = np.array(Bs)       # (5, 3)
errs = np.array(errs)     # (5, 2, 3)

# average of the ROI medians, with propagated (asymmetric) error
Bavg     = Bs.mean(axis=0)
Bavg_err = np.sqrt((errs**2).sum(axis=0)) / len(rois)   # (2, 3)

plt.errorbar(xr, Bavg, yerr=Bavg_err, fmt='o', linestyle='-', color='k',
             capsize=6, lw=2.5, ms=7, elinewidth=2, zorder=10,
             label='avg of ROI medians')

plt.ylabel(r"$B'$ [nT]", fontsize=15)
for r in xr: plt.axvline(x=r, linestyle='dashed', color='k', linewidth=0.2)
plt.xlim(0, 1.1)
plt.yscale('linear')
plt.yticks(fontsize=12)
plt.xticks(xr, ['PSP \n $\\sim0.05$ AU', 'SOLO \n $\\sim0.3$ AU', 'WIND \n $\\sim1$ AU'], fontsize=12)
# plt.legend(fontsize=12)
plt.yscale('linear')



#%%
#%%
# NEW UPDATED:  PROPER SCALING OF X-AXIS
spreceding = Tpreceding/(npreceding**(2/3))
sshoulder = Tshoulder/(nshoulder**(2/3))
sfast = Tfast/(nfast**(2/3))
srarefaction = Trarefaction/(nrarefaction**(2/3))
strailing = Ttrailing/(ntrailing**(2/3))

def s_with_err(T, T_pct, n, n_pct, scale=1.0):
    """T, n in SI (K or eV, m^-3). Returns scaled s = T/n^(2/3) and asymmetric yerr (2,N)."""
    N = len(T)
    eT = np.broadcast_to(np.asarray(pct_to_yerr(T, T_pct)), (2, N))  # [minus, plus]
    en = np.broadcast_to(np.asarray(pct_to_yerr(n, n_pct)), (2, N))
    Tlo, Thi = T - eT[0], T + eT[1]
    nlo, nhi = n - en[0], n + en[1]

    s = T / n**(2/3)
    up   = np.hypot(np.log(Thi/T), (2/3)*np.log(n/nlo))   # high T pairs with low n
    down = np.hypot(np.log(T/Tlo), (2/3)*np.log(nhi/n))   # low T pairs with high n
    yerr = np.array([s*(1 - np.exp(-down)), s*(np.exp(up) - 1)])
    return scale*s, scale*yerr

pops = [
    (Tpreceding,   Tpreceding_pct,   npreceding,   npreceding_pct,   'darkblue', 'leading'),
    (Tshoulder,    Tshoulder_pct,    nshoulder,    nshoulder_pct,    'purple',   'compression'),
    (Tfast,        Tfast_pct,        nfast,        nfast_pct,        'r',        'fast'),
    (Trarefaction, Trarefaction_pct, nrarefaction, nrarefaction_pct, 'green',    'rarefaction'),
    (Ttrailing,    Ttrailing_pct,    ntrailing,    ntrailing_pct,    'skyblue',  'trailing'),
]

for T_, Tp, n_, npct, col, lab in pops:
    y, yerr = s_with_err(T_, Tp, n_, npct, scale=f)
    plt.errorbar(xr, y, yerr=yerr, fmt='o-', color=col,
                 label=rf'${{}}_{{{lab}}}$', capsize=5)

plt.ylabel(r'$T/n^{2/3} \ [eV \ cm^2]$', fontsize=15)
for r in xr: plt.axvline(x=r, linestyle='dashed', color='k', linewidth=0.2)
plt.xlim(0,1.1)
plt.yscale('linear')
plt.yticks(fontsize=12)
plt.xticks(xr, ['PSP \n $\sim0.05 \  AU$','SOLO \n $\sim0.3 \  AU$','WIND \n $\sim1 \ AU$'], fontsize=12)
# plt.legend()


#%%
# Number Density
navg =(npreceding+nshoulder+nfast+ntrailing+nrarefaction)/5 

plt.errorbar(xr,1e-6*npreceding,yerr=1e-6*pct_to_yerr(npreceding,npreceding_pct),fmt='o-',color='darkblue',label = r'$n^{\prime}_{leading}$',capsize=5)
plt.errorbar(xr,1e-6*nshoulder,yerr=1e-6*pct_to_yerr(nshoulder,nshoulder_pct),fmt='o-',color='purple',label = r'$n^{\prime}_{compression}$',capsize=5)
plt.errorbar(xr,1e-6*nfast,yerr=1e-6*pct_to_yerr(nfast,nfast_pct),fmt='o-',color = 'r',label = r'$n^{\prime}_{fast}$',capsize=5)
plt.errorbar(xr,1e-6*nrarefaction,yerr=1e-6*pct_to_yerr(nrarefaction,nrarefaction_pct),fmt='o-',color='green',label = r'$n^{\prime}_{rarefaction}$',capsize=5)
plt.errorbar(xr,1e-6*ntrailing,yerr=1e-6*pct_to_yerr(ntrailing,ntrailing_pct),fmt='o-',color='skyblue',label = r'$n^{\prime}_{trailing}$',capsize=5)
# plt.plot(xr,1e-6*navg,color='k',linestyle='dashed')
plt.ylabel('$n\'$ [$cm^{-3}$]',fontsize=15)
for r in xr: plt.axvline(x=r,linestyle='dashed',color='k',linewidth=0.2)
plt.xlim(0,1.1)
plt.yticks(fontsize=12)
# plt.yscale('log')
plt.xticks(xr,['PSP \n $\sim0.05 \  AU$','SOLO \n $\sim0.3 \  AU$','WIND \n $\sim1 \ AU$	'],fontsize=12)
# plt.legend(fontsize=12)

#%%
# Temperature
f = (1.602e-19)**-1
plt.errorbar(xr,f*Tpreceding,yerr=f*pct_to_yerr(Tpreceding,Tpreceding_pct),fmt='o-',color='darkblue',label = r'$T^{\prime}_{leading}$',capsize=5)
plt.errorbar(xr,f*Tshoulder,yerr=f*pct_to_yerr(Tshoulder,Tshoulder_pct),fmt='o-',color='purple',label = r'$T^{\prime}_{compression}$',capsize=5)
plt.errorbar(xr,f*Tfast,yerr=f*pct_to_yerr(Tfast,Tfast_pct),fmt='o-',color = 'r',label = r'$T^{\prime}_{fast}$',capsize=5)
plt.errorbar(xr,f*Trarefaction,yerr=f*pct_to_yerr(Trarefaction,Trarefaction_pct),fmt='o-',color='green',label = r'$T^{\prime}_{rarefaction}$',capsize=5)
plt.errorbar(xr,f*Ttrailing,yerr=f*pct_to_yerr(Ttrailing,Ttrailing_pct),fmt='o-',color='skyblue',label = r'$T^{\prime}_{trailing}$',capsize=5)
plt.ylabel('$T\'$ [$eV$]',fontsize=15)
plt.yscale('linear')
for r in xr: plt.axvline(x=r,linestyle='dashed',color='k',linewidth=0.2)
plt.xlim(0,1.1)
plt.yticks(fontsize=12)
plt.xticks(xr,['PSP \n $\sim0.05 \  AU$','SOLO \n $\sim0.3 \  AU$','WIND \n $\sim1 \ AU$	'],fontsize=12)
# plt.legend(fontsize=12)




#%%
# Velocity
vavg = (vpreceding+vshoulder+vfast+vrarefaction+vtrailing)/5
plt.errorbar(xr,1e-3*vpreceding,yerr=1e-3*pct_to_yerr(vpreceding,vpreceding_pct),fmt='o-',color='darkblue',label = r'$v_{leading}$',capsize=5)
plt.errorbar(xr,1e-3*vshoulder,yerr=1e-3*pct_to_yerr(vshoulder,vshoulder_pct),fmt='o-',color='purple',label = r'$v_{compression}$',capsize=5)
plt.errorbar(xr,1e-3*vfast,yerr=1e-3*pct_to_yerr(vfast,vfast_pct),fmt='o-',color = 'r',label = r'$v_{fast}$',capsize=5)
plt.errorbar(xr,1e-3*vrarefaction,yerr=1e-3*pct_to_yerr(vrarefaction,vrarefaction_pct),fmt='o-',color='green',label = r'$v_{rarefaction}$',capsize=5)
plt.errorbar(xr,1e-3*vtrailing,yerr=1e-3*pct_to_yerr(vtrailing,vtrailing_pct),fmt='o-',color='skyblue',label = r'$v_{trailing}$',capsize=5)
plt.ylabel('$v_R$ [$km/s$]',fontsize=15)
for r in xr: plt.axvline(x=r,linestyle='dashed',color='k',linewidth=0.2)
plt.xlim(0,1.1)
plt.ylim(150,900)
plt.yticks(fontsize=12)
plt.xticks(xr,['PSP \n $\sim0.05 \  AU$','SOLO \n $\sim0.3 \  AU$','WIND \n $\sim1 \ AU$	'],fontsize=12)
# plt.legend(fontsize=12)


#%%
# Magnetic Field Strength
f = 1e9
Bavg = (Bpreceding+Bshoulder+Bfast+Brarefaction+Btrailing)/5
plt.errorbar(xr,f*Bpreceding,yerr=f*pct_to_yerr(Bpreceding,Bpreceding_pct),fmt='o-',color='darkblue',label = r'$|B|_{leading}$',capsize=5)
plt.errorbar(xr,f*Bshoulder,yerr=f*pct_to_yerr(Bshoulder,Bshoulder_pct),fmt='o-',color='purple',label = r'$|B|_{compression}$',capsize=5)
plt.errorbar(xr,f*Bfast,yerr=f*pct_to_yerr(Bfast,Bfast_pct),fmt='o-',color = 'r',label = r'$|B|_{fast}$',capsize=5)
plt.errorbar(xr,f*Brarefaction,yerr=f*pct_to_yerr(Brarefaction,Brarefaction_pct),fmt='o-',color='green',label = r'$|B|_{rarefaction}$',capsize=5)
plt.errorbar(xr,f*Btrailing,yerr=f*pct_to_yerr(Btrailing,Btrailing_pct),fmt='o-',color='skyblue',label = r'$|B|_{trailing}$',capsize=5)
# plt.plot(xr,f*Bavg, linestyle='dotted',color='k',linewidth=3)
# plt.ylabel('$B\'$ [$nT$]',fontsize=15)
plt.ylabel('$B_R\'$ [$nT$]',fontsize=15)
for r in xr: plt.axvline(x=r,linestyle='dashed',color='k',linewidth=0.2)
plt.xlim(0,1.1)
plt.yticks(fontsize=12)
plt.yscale('linear')
plt.yticks(fontsize=12)
plt.xticks(xr,['PSP \n $\sim0.05 \  AU$','SOLO \n $\sim0.3 \  AU$','WIND \n $\sim1 \ AU$	'],fontsize=12)
# plt.legend(fontsize=12)

#%%

#%%
# Bulk Energy Density (interval propagation: bulk = 0.5*n*mi*v^2, monotonic increasing in n and v)

f=(6.242e18)*(1e-6)*(1e-6) #MeV/cm^3
bulkfast_pct = np.stack([0.5*nfast_pct[:,0]*mi*vfast_pct[:,0]**2, 0.5*nfast_pct[:,1]*mi*vfast_pct[:,1]**2],axis=1)
bulkpreceding_pct = np.stack([0.5*npreceding_pct[:,0]*mi*vpreceding_pct[:,0]**2, 0.5*npreceding_pct[:,1]*mi*vpreceding_pct[:,1]**2],axis=1)
bulktrailing_pct = np.stack([0.5*ntrailing_pct[:,0]*mi*vtrailing_pct[:,0]**2, 0.5*ntrailing_pct[:,1]*mi*vtrailing_pct[:,1]**2],axis=1)
bulkshoulder_pct = np.stack([0.5*nshoulder_pct[:,0]*mi*vshoulder_pct[:,0]**2, 0.5*nshoulder_pct[:,1]*mi*vshoulder_pct[:,1]**2],axis=1)
bulkrarefaction_pct = np.stack([0.5*nrarefaction_pct[:,0]*mi*vrarefaction_pct[:,0]**2, 0.5*nrarefaction_pct[:,1]*mi*vrarefaction_pct[:,1]**2],axis=1)

plt.errorbar(xr,f*bulkpreceding,yerr=f*pct_to_yerr(bulkpreceding,bulkpreceding_pct),fmt='o-',color='darkblue',label = r'Leading',capsize=5)
plt.errorbar(xr,f*bulkshoulder,yerr=f*pct_to_yerr(bulkshoulder,bulkshoulder_pct),fmt='o-',color='purple',label = r'Compression',capsize=5)
plt.errorbar(xr,f*bulkfast,yerr=f*pct_to_yerr(bulkfast,bulkfast_pct),fmt='o-',color = 'r',label = r'Fast Wind',capsize=5)
plt.errorbar(xr,f*bulkrarefaction,yerr=f*pct_to_yerr(bulkrarefaction,bulkrarefaction_pct),fmt='o-',color='green',label = r'Rarefaction',capsize=5)
plt.errorbar(xr,f*bulktrailing,yerr=f*pct_to_yerr(bulktrailing,bulktrailing_pct),fmt='o-',color='skyblue',label = r'Trailing',capsize=5)
plt.ylabel('$u_{k}\'$ [$MeV \ cm^{-3}$]',fontsize=15)
for r in xr: plt.axvline(x=r,linestyle='dashed',color='k',linewidth=0.2)
plt.xlim(0,1.1)
plt.yscale('linear')
plt.yticks(fontsize=12)
plt.xticks(xr,['PSP \n $\sim0.05 \  AU$','SOLO \n $\sim0.3 \  AU$','WIND \n $\sim1 \ AU$	'],fontsize=12)
# plt.legend(fontsize=12)

#%%
# Thermal Energy Density (interval propagation: therm = 1.5*kb*n*T, monotonic increasing in n and T)
thermfast_pct = np.stack([1.5*kb*nfast_pct[:,0]*Tfast_pct[:,0], 1.5*kb*nfast_pct[:,1]*Tfast_pct[:,1]],axis=1)
thermpreceding_pct = np.stack([1.5*kb*npreceding_pct[:,0]*Tpreceding_pct[:,0], 1.5*kb*npreceding_pct[:,1]*Tpreceding_pct[:,1]],axis=1)
thermtrailing_pct = np.stack([1.5*kb*ntrailing_pct[:,0]*Ttrailing_pct[:,0], 1.5*kb*ntrailing_pct[:,1]*Ttrailing_pct[:,1]],axis=1)
thermshoulder_pct = np.stack([1.5*kb*nshoulder_pct[:,0]*Tshoulder_pct[:,0], 1.5*kb*nshoulder_pct[:,1]*Tshoulder_pct[:,1]],axis=1)
thermrarefaction_pct = np.stack([1.5*kb*nrarefaction_pct[:,0]*Trarefaction_pct[:,0], 1.5*kb*nrarefaction_pct[:,1]*Trarefaction_pct[:,1]],axis=1)

plt.errorbar(xr,f*thermpreceding,yerr=f*pct_to_yerr(thermpreceding,thermpreceding_pct),fmt='o-',color='darkblue',label = r'Leading Slow Wind',capsize=5)
plt.errorbar(xr,f*thermshoulder,yerr=f*pct_to_yerr(thermshoulder,thermshoulder_pct),fmt='o-',color='purple',label = r'Compression Region',capsize=5)
plt.errorbar(xr,f*thermfast,yerr=f*pct_to_yerr(thermfast,thermfast_pct),fmt='o-',color = 'r',label = r'Fast CH Wind',capsize=5)
plt.errorbar(xr,f*thermrarefaction,yerr=f*pct_to_yerr(thermrarefaction,thermrarefaction_pct),fmt='o-',color='green',label = r'Rarefaction Region',capsize=5)
plt.errorbar(xr,f*thermtrailing,yerr=f*pct_to_yerr(thermtrailing,thermtrailing_pct),fmt='o-',color='skyblue',label = r'Trailing Slow Wind',capsize=5)
plt.ylabel('$u_{th}\'$ [$MeV \ cm^{-3}$]',fontsize=15)
for r in xr: plt.axvline(x=r,linestyle='dashed',color='k',linewidth=0.2)
plt.xlim(0,1.1)
plt.yscale('linear')
plt.yticks(fontsize=12)
plt.xticks(xr,['PSP \n $\sim0.05 \  AU$','SOLO \n $\sim0.3 \  AU$','WIND \n $\sim1 \ AU$	'],fontsize=12)
# plt.legend(fontsize=12)

# %%
# Magnetic Energy Density (interval propagation: mag = 0.5*B^2/mu0, monotonic increasing in B for B>0)
magfast_pct = np.stack([0.5*Bfast_pct[:,0]**2/mu0, 0.5*Bfast_pct[:,1]**2/mu0],axis=1)
magpreceding_pct = np.stack([0.5*Bpreceding_pct[:,0]**2/mu0, 0.5*Bpreceding_pct[:,1]**2/mu0],axis=1)
magtrailing_pct = np.stack([0.5*Btrailing_pct[:,0]**2/mu0, 0.5*Btrailing_pct[:,1]**2/mu0],axis=1)
magshoulder_pct = np.stack([0.5*Bshoulder_pct[:,0]**2/mu0, 0.5*Bshoulder_pct[:,1]**2/mu0],axis=1)
magrarefaction_pct = np.stack([0.5*Brarefaction_pct[:,0]**2/mu0, 0.5*Brarefaction_pct[:,1]**2/mu0],axis=1)

plt.errorbar(xr,f*magpreceding,yerr=f*pct_to_yerr(magpreceding,magpreceding_pct),fmt='o-',color='darkblue',label = r'Leading',capsize=5)
plt.errorbar(xr,f*magshoulder,yerr=f*pct_to_yerr(magshoulder,magshoulder_pct),fmt='o-',color='purple',label = r'Compression',capsize=5)
plt.errorbar(xr,f*magfast,yerr=f*pct_to_yerr(magfast,magfast_pct),fmt='o-',color = 'r',label = r'Fast',capsize=5)
plt.errorbar(xr,f*magrarefaction,yerr=f*pct_to_yerr(magrarefaction,magrarefaction_pct),fmt='o-',color='green',label = r'Rarefaction',capsize=5)
plt.errorbar(xr,f*magtrailing,yerr=f*pct_to_yerr(magtrailing,magtrailing_pct),fmt='o-',color='skyblue',label = r'Trailing',capsize=5)
plt.ylabel('$u_{m}\'$ [$MeV \ cm^{-3}$]',fontsize=15)
for r in xr: plt.axvline(x=r,linestyle='dashed',color='k',linewidth=0.2)
plt.xlim(0,1.1)
plt.yscale('linear')
plt.yticks(fontsize=12)
plt.xticks(xr,['PSP \n $\sim0.05 \  AU$','SOLO \n $\sim0.3 \  AU$','WIND \n $\sim1 \ AU$	'],fontsize=12)
# plt.legend(fontsize=20)

#%%

# %%
#%%

# Cperp
f = 1e-9 / kb
plt.errorbar(xr,f*Cperppreceding,yerr=f*pct_to_yerr(Cperppreceding,Cperppreceding_pct),fmt='o-',color='darkblue',label = r'$Leading$',capsize=5)
plt.errorbar(xr,f*Cperpshoulder,yerr=f*pct_to_yerr(Cperpshoulder,Cperpshoulder_pct),fmt='o-',color='purple',label = r'$Compression$',capsize=5)
plt.errorbar(xr,f*Cperpfast,yerr=f*pct_to_yerr(Cperpfast,Cperpfast_pct),fmt='o-',color = 'r',label = r'$Fast$',capsize=5)
plt.errorbar(xr,f*Cperprarefaction,yerr=f*pct_to_yerr(Cperprarefaction,Cperprarefaction_pct),fmt='o-',color='green',label = r'$Rarefaction$',capsize=5)
plt.errorbar(xr,f*Cperptrailing,yerr=f*pct_to_yerr(Cperptrailing,Cperptrailing_pct),fmt='o-',color='skyblue',label = r'$Trailing$',capsize=5)
plt.ylabel(r'$C_\perp \ [K/nT]$',fontsize=15)
for r in xr: plt.axvline(x=r,linestyle='dashed',color='k',linewidth=0.2)
plt.xlim(0,1.1)
plt.yscale('log')
plt.yticks(fontsize=12)
plt.xticks(xr,['PSP \n $\sim0.05 \  AU$','SOLO \n $\sim0.3 \  AU$','WIND \n $\sim1 \ AU$	'],fontsize=12)
plt.legend(fontsize=12)
#%%
# Cpar
f = 1e30 / kb
plt.errorbar(xr,f*Cparpreceding,yerr=f*pct_to_yerr(Cparpreceding,Cparpreceding_pct),fmt='o-',color='darkblue',label = r'$Leading$',capsize=5)
plt.errorbar(xr,f*Cparshoulder,yerr=f*pct_to_yerr(Cparshoulder,Cparshoulder_pct),fmt='o-',color='purple',label = r'$Compression$',capsize=5)
plt.errorbar(xr,f*Cparfast,yerr=f*pct_to_yerr(Cparfast,Cparfast_pct),fmt='o-',color = 'r',label = r'$Fast$',capsize=5)
plt.errorbar(xr,f*Cparrarefaction,yerr=f*pct_to_yerr(Cparrarefaction,Cparrarefaction_pct),fmt='o-',color='green',label = r'$Rarefaction$',capsize=5)
plt.errorbar(xr,f*Cpartrailing,yerr=f*pct_to_yerr(Cpartrailing,Cpartrailing_pct),fmt='o-',color='skyblue',label = r'$Trailing$',capsize=5)
plt.ylabel(r'$C_\parallel \ [K \ nT^2 \  cm^6]$',fontsize=15)
for r in xr: plt.axvline(x=r,linestyle='dashed',color='k',linewidth=0.2)
plt.xlim(0,1.1)
plt.yscale('log')
plt.yticks(fontsize=12)
plt.xticks(xr,['PSP \n $\sim0.05 \  AU$','SOLO \n $\sim0.3 \  AU$','WIND \n $\sim1 \ AU$	'],fontsize=12)
# plt.legend(fontsize=12,framealpha=1)











#%%
# Cperp normalized to PSP (dividing by a positive constant preserves monotonicity, so bounds divide too)
plt.errorbar(xr,Cperppreceding/Cperppreceding[0],yerr=pct_to_yerr(Cperppreceding/Cperppreceding[0],Cperppreceding_pct/Cperppreceding[0]),fmt='o-',color='darkblue',label = r'$Leading$',capsize=5)
plt.errorbar(xr,Cperpshoulder/Cperpshoulder[0],yerr=pct_to_yerr(Cperpshoulder/Cperpshoulder[0],Cperpshoulder_pct/Cperpshoulder[0]),fmt='o-',color='purple',label = r'$Compression$',capsize=5)
plt.errorbar(xr,Cperpfast/Cperpfast[0],yerr=pct_to_yerr(Cperpfast/Cperpfast[0],Cperpfast_pct/Cperpfast[0]),fmt='o-',color = 'r',label = r'$Fast$',capsize=5)
plt.errorbar(xr,Cperprarefaction/Cperprarefaction[0],yerr=pct_to_yerr(Cperprarefaction/Cperprarefaction[0],Cperprarefaction_pct/Cperprarefaction[0]),fmt='o-',color='green',label = r'$Rarefaction$',capsize=5)
plt.errorbar(xr,Cperptrailing/Cperptrailing[0],yerr=pct_to_yerr(Cperptrailing/Cperptrailing[0],Cperptrailing_pct/Cperptrailing[0]),fmt='o-',color='skyblue',label = r'$Trailing$',capsize=5)
plt.ylabel('$C_\perp/C_{\perp psp}$',fontsize=15)
plt.axhline(y=1,color='k',linestyle='dashed')
for r in xr: plt.axvline(x=r,linestyle='dashed',color='k',linewidth=0.2)
plt.xlim(0,1.1)
plt.yscale('log')
plt.ylim(1e-1,10**2.2)
plt.yticks(fontsize=12)
plt.xticks(xr,['PSP \n $\sim0.05 \  AU$','SOLO \n $\sim0.3 \  AU$','WIND \n $\sim1 \ AU$	'],fontsize=12)
# plt.legend(fontsize=12)
#%%
# Cpar normalized to PSP
plt.errorbar(xr,Cparpreceding/Cparpreceding[0],yerr=pct_to_yerr(Cparpreceding/Cparpreceding[0],Cparpreceding_pct/Cparpreceding[0]),fmt='o-',color='darkblue',label = r'$Leading$',capsize=5)
plt.errorbar(xr,Cparshoulder/Cparshoulder[0],yerr=pct_to_yerr(Cparshoulder/Cparshoulder[0],Cparshoulder_pct/Cparshoulder[0]),fmt='o-',color='purple',label = r'$Compression$',capsize=5)
plt.errorbar(xr,Cparfast/Cparfast[0],yerr=pct_to_yerr(Cparfast/Cparfast[0],Cparfast_pct/Cparfast[0]),fmt='o-',color = 'r',label = r'$Fast$',capsize=5)
plt.errorbar(xr,Cparrarefaction/Cparrarefaction[0],yerr=pct_to_yerr(Cparrarefaction/Cparrarefaction[0],Cparrarefaction_pct/Cparrarefaction[0]),fmt='o-',color='green',label = r'$Rarefaction$',capsize=5)
plt.errorbar(xr,Cpartrailing/Cpartrailing[0],yerr=pct_to_yerr(Cpartrailing/Cpartrailing[0],Cpartrailing_pct/Cpartrailing[0]),fmt='o-',color='skyblue',label = r'$Trailing$',capsize=5)
plt.ylabel('$C_\parallel/C_{\parallel psp}$',fontsize=15)
plt.axhline(y=1,color='k',linestyle='dashed')
for r in xr: plt.axvline(x=r,linestyle='dashed',color='k',linewidth=0.2)
plt.xlim(0,1.1)
plt.yscale('log')
plt.ylim(1e-1,10**2.2)
plt.yticks(fontsize=12)
plt.xticks(xr,['PSP \n $\sim0.05 \  AU$','SOLO \n $\sim0.3 \  AU$','WIND \n $\sim1 \ AU$	'],fontsize=12)
# plt.legend(fontsize=12,framealpha=1)


# %%
