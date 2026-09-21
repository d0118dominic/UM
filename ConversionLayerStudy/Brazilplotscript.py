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


def brazil_multi(betapar0=None, Tparperp0=None,
                  betapar1=None, Tparperp1=None,
                  betapar2=None, Tparperp2=None,
                  betapar3=None, Tparperp3=None,
                  betapar4=None, Tparperp4=None,
                  mask=None,colors_dict={'0': 'grey', '1': 'darkblue', '2': 'limegreen', '3': 'orange', '4': 'red'},
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

    data = {
        '0': (betapar0, Tparperp0, mask),
        '1': (betapar1, Tparperp1, mask),
        '2': (betapar2, Tparperp2, mask),
        '3': (betapar3, Tparperp3, mask),
        '4': (betapar4, Tparperp4, mask),
    }

    fig, ax = plt.subplots(figsize=(10, 8))

    for sc, (betapar, Tparperp, mask) in data.items():
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

dva_norm0 = abs(dva[0]/vameans[0])
dva_norm1 = abs(dva[1]/vameans[1])
dva_norm2 = abs(dva[2]/vameans[2])
dva_norm3 = abs(dva[3]/vameans[3])



lperp_over_lpar0 = lperp_over_lpar[0]
lperp_over_lpar1 = lperp_over_lpar[1]
lperp_over_lpar2 = lperp_over_lpar[2]
lperp_over_lpar3 = lperp_over_lpar[3]
lperp_over_lpar4 = lperp_over_lpar[4]

threshold=0.05


def multimask(var,threshold,lcond = None):
    mask = []
    for i in range(len(var)):
        mask.append((var[i]>=threshold))
    return mask


def multimask(var,threshold,lcond = None):
    mask = []
    if lcond == 'perp':
        for i in range(len(var)):
            mask.append((var[i]>=threshold) & (lperp_over_lpar[i]>=1))
    elif lcond == 'par':
        for i in range(len(var)):
            mask.append((var[i]>=threshold) & (lperp_over_lpar[i]<=1))
    else:
        for i in range(len(var)):
            mask.append((var[i]>=threshold))
    return mask


dva_norm = np.array(dva)/np.array(vameans)
dvth_norm = np.array(dvth)/np.array(vthmeans)
mask = multimask(np.abs(dva_norm),0.0,lcond=None)
mask = multimask(np.abs(dvth_norm),0.05,lcond='par')


fig, ax = brazil_multi(betapar0=betapar[0][mask[0]], Tparperp0=Tparperp[0][mask[0]],
                        betapar1=betapar[1][mask[1]], Tparperp1=Tparperp[1][mask[1]],
                        betapar2=betapar[2][mask[2]], Tparperp2=Tparperp[2][mask[2]],
                        betapar3=betapar[3][mask[3]], Tparperp3=Tparperp[3][mask[3]],
                        # betapar4=betapar[4], Tparperp4=Tparperp[4],
                        contour_nlevels=10,title = "dva gradients > $5\%$ ",plot_type='contour',
                        mean_marker_size=100)

#%%
fig, ax = brazil_multi(betapar0=betapar[0][dva_norm0>=threshold], Tparperp0=Tparperp[0][dva_norm0>=threshold],
                        betapar1=betapar[1][dva_norm1>=threshold], Tparperp1=Tparperp[1][dva_norm1>=threshold],
                        betapar2=betapar[2][dva_norm2>=threshold], Tparperp2=Tparperp[2][dva_norm2>=threshold],
                        betapar3=betapar[3][dva_norm3>=threshold], Tparperp3=Tparperp[3][dva_norm3>=threshold],
                        # betapar4=betapar[4], Tparperp4=Tparperp[4],
                        contour_nlevels=10,title = "dva gradients > $5\%$ ",plot_type='contour',
                        mean_marker_size=100)


#%%
fig, ax = brazil_multi(betapar0=betapar[0][lperp_over_lpar0<=1], Tparperp0=Tparperp[0][lperp_over_lpar0<=1],
                        betapar1=betapar[1][lperp_over_lpar1<=1], Tparperp1=Tparperp[1][lperp_over_lpar1<=1],
                        betapar2=betapar[2][lperp_over_lpar2<=1], Tparperp2=Tparperp[2][lperp_over_lpar2<=1],
                        betapar3=betapar[3][lperp_over_lpar3<=1], Tparperp3=Tparperp[3][lperp_over_lpar3<=1],
                        # betapar4=betapar[4], Tparperp4=Tparperp[4],
                        contour_nlevels=10,title = "$\parallel$ gradients",plot_type='contour',
                        mean_marker_size=100)


fig, ax = brazil_multi(betapar0=betapar[0][lperp_over_lpar0>=1], Tparperp0=Tparperp[0][lperp_over_lpar0>=1],
                        betapar1=betapar[1][lperp_over_lpar1>=1], Tparperp1=Tparperp[1][lperp_over_lpar1>=1],
                        betapar2=betapar[2][lperp_over_lpar2>=1], Tparperp2=Tparperp[2][lperp_over_lpar2>=1],
                        betapar3=betapar[3][lperp_over_lpar3>=1], Tparperp3=Tparperp[3][lperp_over_lpar3>=1],
                        # betapar4=betapar[4], Tparperp4=Tparperp[4],
                        contour_nlevels=10,title = " $\perp$ gradients ",plot_type='contour',
                        mean_marker_size=100)