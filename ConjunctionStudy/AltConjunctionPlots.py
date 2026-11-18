
#%%
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import binned_statistic_2d
from IPython.display import Image
from matplotlib.animation import FuncAnimation, PillowWriter

# Same colors as the v_{region} line plots
REGION_COLORS = {
    'leading':     'darkblue',
    'compression': 'purple',
    'fast':        'r',
    'rarefaction': 'green',
    'trailing':    'black',
}

def brazil_regions(region_data, sc_name='',
                   colors_dict=REGION_COLORS,
                   x_label=r'$\beta_\parallel$', y_label=r'$T_{\perp}/T_{\parallel}$',
                   title=None, nbins=50, mincount=10,
                   xrange=[-3, 3], yrange=[-1, 1],
                   contour_nlevels=5, norm=False,
                   mean_marker_size=120,
                   cgl_line=True, cgl_exponent=-1,
                   ax=None):
    """
    Brazil plot for ONE spacecraft with several solar wind regions overplotted.

    region_data : dict {region_name: (betapar, Tparperp, mask)}
                  mask may be None (defaults to Tparperp > 0).
                  Tparperp is inverted internally (y = 1/Tparperp), as in brazil_multi.
    """
    # From Hellinger (2011)
    mirror_params      = [0.77, 0.76, -0.016]
    firehose_params    = [-1.4, 1.0, -0.11]
    cyclotron_params   = [0.43, 0.42, -0.0004]
    parfirehose_params = [-0.47, 0.53, 0.59]

    beta_par_range = np.logspace(xrange[0], xrange[1], 100)
    mirror_threshold      = stability_condition(beta_par_range, *mirror_params)
    firehose_threshold    = stability_condition(beta_par_range, *firehose_params)
    cyclotron_threshold   = stability_condition(beta_par_range, *cyclotron_params)
    parfirehose_threshold = stability_condition(beta_par_range, *parfirehose_params)

    xbins = np.logspace(xrange[0], xrange[1], nbins)
    ybins = np.logspace(yrange[0], yrange[1], nbins)
    xcenters = 0.5 * (xbins[:-1] + xbins[1:])
    ycenters = 0.5 * (ybins[:-1] + ybins[1:])
    xline_cgl = np.logspace(xrange[0], xrange[1], 100)
    xgrid, ygrid = np.meshgrid(xcenters, ycenters, indexing='ij')

    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 8))
    else:
        fig = ax.figure

    for region, (betapar, Tparperp, mask) in region_data.items():
        if betapar is None or Tparperp is None:
            continue
        if mask is None:
            mask = Tparperp > 0

        x = betapar[mask]
        y = 1 / Tparperp[mask]
        valmask = np.isfinite(x) & np.isfinite(y)
        x_valid, y_valid = x[valmask], y[valmask]

        color = colors_dict.get(region, 'k')

        counts, _, _, _ = binned_statistic_2d(
            x_valid, y_valid, None, statistic='count', bins=[xbins, ybins])

        if counts.max() <= mincount:
            continue

        if norm:
            total = x_valid.size
            plot_counts  = counts / total
            plot_mincount = mincount / total
            plot_max     = counts.max() / total
        else:
            plot_counts  = counts
            plot_mincount = mincount
            plot_max     = counts.max()

        levels = np.logspace(np.log10(plot_mincount), np.log10(plot_max), contour_nlevels)
        ax.contour(xcenters, ycenters, plot_counts.T, levels=levels,
                   colors=color, linewidths=0.8, alpha=0.4)

        # counts-weighted centroid of the binned distribution
        total_counts = counts.sum()
        if total_counts > 0:
            x_mean = np.sum(xgrid * counts) / total_counts
            y_mean = np.sum(ygrid * counts) / total_counts
            ax.scatter(x_mean, y_mean, marker='^', color=color, s=2*mean_marker_size,
                       edgecolors=color, linewidths=1.5, zorder=10, label=region.capitalize())

            if cgl_line:
                yline_cgl = y_mean * (xline_cgl / x_mean) ** cgl_exponent
                ax.plot(xline_cgl, yline_cgl, color=color, linestyle='-',
                        linewidth=1, alpha=0.8, zorder=9)

    ax.set_facecolor('lightgray')
    ax.set_xlabel(x_label, fontdict={'size': 20})
    ax.set_ylabel(y_label, fontdict={'size': 20})
    ax.set_title(sc_name if title is None else title, fontdict={'size': 20})
    ax.set_xlim(10**(xrange[0]), 10**(xrange[1]))
    ax.set_ylim(10**(-1), 10**(1))
    ax.axhline(y=1, color='k', linestyle='-')
    ax.axvline(x=1, color='k', linestyle='-')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.tick_params(labelsize=15)
    ax.plot(beta_par_range, mirror_threshold,      label='Mirror',       color='k', linestyle='dotted', linewidth=2)
    ax.plot(beta_par_range, firehose_threshold,    label='Oblique FH',   color='k', linestyle='--',     linewidth=2)
    ax.plot(beta_par_range, cyclotron_threshold,   label='IC',           color='k', linestyle='dotted', linewidth=1)
    ax.plot(beta_par_range, parfirehose_threshold, label='Parallel FH',  color='k', linestyle='--',     linewidth=1)
    ax.legend(fontsize=13, framealpha=1, loc='upper right')
    return fig, ax


def brazil_by_spacecraft(sc_region_data, one_figure=False, **kwargs):
    """
    sc_region_data : dict {'PSP': {region: (beta, Tparperp, mask), ...},
                           'SOLO': {...}, 'WIND': {...}}
    one_figure=True puts the spacecraft side by side in a single 1xN figure;
    otherwise each spacecraft gets its own figure.
    Returns a dict {spacecraft: (fig, ax)}.
    """
    out = {}
    if one_figure:
        n = len(sc_region_data)
        fig, axes = plt.subplots(1, n, figsize=(9*n, 8), sharey=True)
        axes = np.atleast_1d(axes)
        for ax, (sc, rd) in zip(axes, sc_region_data.items()):
            brazil_regions(rd, sc_name=sc, ax=ax, **kwargs)
            out[sc] = (fig, ax)
        for ax in axes[1:]:
            ax.set_ylabel('')
        fig.tight_layout()
    else:
        for sc, rd in sc_region_data.items():
            out[sc] = brazil_regions(rd, sc_name=sc, **kwargs)
    plt.show()
    return out
# %%
sc_region_data = {
    'PSP': {
        'leading':     (betaparspsp_preceding,   Tparperppsp_preceding,   None),
        'compression': (betaparspsp_shoulder,    Tparperppsp_shoulder,    None),
        'fast':        (betaparspsp_fast,        Tparperppsp_fast,        None),
        'rarefaction': (betaparspsp_rarefaction, Tparperppsp_rarefaction, None),
        'trailing':    (betaparspsp_trailing,    Tparperppsp_trailing,    None),
    },
    'SOLO': {
        'leading':     (betaparssolo_preceding,   Tparperpsolo_preceding,   None),
        'compression': (betaparssolo_shoulder,    Tparperpsolo_shoulder,    None),
        'fast':        (betaparssolo_fast,        Tparperpsolo_fast,        None),
        'rarefaction': (betaparssolo_rarefaction, Tparperpsolo_rarefaction, None),
        'trailing':    (betaparssolo_trailing,    Tparperpsolo_trailing,    None),
    },
    'WIND': {
        'leading':     (betaparswind_preceding,   Tparperpwind_preceding,   None),
        'compression': (betaparswind_shoulder,    Tparperpwind_shoulder,    None),
        'fast':        (betaparswind_fast,        Tparperpwind_fast,        None),
        'rarefaction': (betaparswind_rarefaction, Tparperpwind_rarefaction, None),
        'trailing':    (betaparswind_trailing,    Tparperpwind_trailing,    None),
    },
}

figs = brazil_by_spacecraft(sc_region_data,
                            one_figure=False,      # True for a 1x3 side-by-side figure
                            contour_nlevels=10,
                            mean_marker_size=35,
                            cgl_line=False,
                            cgl_exponent=-1,norm=False)
# %%
keep = ['fast', 'compression', 'rarefaction']   # fast CH wind + the two slow wind bins

sc_region_data_sub = {
    sc: {r: v for r, v in rd.items() if r in keep}
    for sc, rd in sc_region_data.items()
}

figs = brazil_by_spacecraft(sc_region_data_sub,
                            one_figure=False,
                            contour_nlevels=10,
                            mean_marker_size=35,
                            cgl_line=False,
                            cgl_exponent=-1,norm=False)
# %%

def brazil_movie(sc_region_data, order=('PSP', 'SOLO', 'WIND'),
                 labels=None, filename='brazil_evolution.gif',
                 fps=1, dpi=150, **kwargs):
    """
    sc_region_data : dict {'PSP': {region: (beta, Tparperp, mask)}, ...}
    order          : frame order (one frame per spacecraft)
    labels         : optional dict of frame titles, e.g. {'PSP': 'PSP (~0.05 AU)'}
    kwargs         : passed through to brazil_regions
    """
    if labels is None:
        labels = {'PSP': 'PSP (~0.05 AU)',
                  'SOLO': 'SOLO (~0.3 AU)',
                  'WIND': 'WIND (~1 AU)'}
    
    # if labels is None:
    #     labels = {'PSP': 'Below $R_A$',
    #               'SOLO': 'Above $R_A$',
    #               'WIND': 'WIND (~1 AU)'}

    fig, ax = plt.subplots(figsize=(10, 8))

    def update(i):
        sc = order[i]
        ax.clear()                      # wipe the previous frame
        brazil_regions(sc_region_data[sc], sc_name=labels.get(sc, sc),
                       ax=ax, **kwargs)
        return ax,

    anim = FuncAnimation(fig, update, frames=len(order), blit=False)
    anim.save(filename, writer=PillowWriter(fps=fps), dpi=dpi)
    plt.close(fig)
    return anim



#%%
anim = brazil_movie(sc_region_data_sub,
                    filename='brazil_fast_slow.gif',
                    fps=2,                 # 1 s per frame; try 0.7 to 1.5
                    contour_nlevels=10,
                    mean_marker_size=35,
                    cgl_line=False,
                    cgl_exponent=-1)

# Image(filename='brazil_slow.gif')
Image(filename='brazil_fast.gif')
# %%
# SIMPLE ANIMATION EXAMPLE
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

# 1. Make one figure and one axes. Every frame is drawn on these.
fig, ax = plt.subplots()

# 2. Your data: one entry per frame.
x = np.linspace(0, 10, 200)
curves = [np.sin(x), np.sin(2*x), np.sin(3*x)]
titles = ['Frame 1: sin(x)', 'Frame 2: sin(2x)', 'Frame 3: sin(3x)']

# 3. A function that draws frame i.
def draw_frame(i):
    ax.clear()                      # erase the previous frame
    ax.plot(x, curves[i])           # draw this frame's data
    ax.set_title(titles[i])
    ax.set_ylim(-1.2, 1.2)          # re-set limits, since clear() resets them

# 4. Tell matplotlib to call draw_frame for i = 0, 1, 2.
anim = FuncAnimation(fig, draw_frame, frames=3)

# 5. Save it.
anim.save('basic_movie.gif', writer=PillowWriter(fps=1))
Image(filename='basic_movie.gif')
# %%


# ALTERNATE WITH TRAILS







import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import binned_statistic_2d
from matplotlib.animation import FuncAnimation, PillowWriter

# ---------------------------------------------------------------
# Helper (skipped automatically if you already defined your own)
# ---------------------------------------------------------------
if 'stability_condition' not in globals():
    def stability_condition(beta, a, b, beta0):
        return 1.0 + a / (beta - beta0) ** b

# Same colors as the v_{region} line plots
REGION_COLORS = {
    'leading':     'darkblue',
    'compression': 'purple',
    'fast':        'r',
    'rarefaction': 'green',
    'trailing':    'dodgerblue',
}

# ---------------------------------------------------------------
# One spacecraft, several regions
# ---------------------------------------------------------------
def brazil_regions(region_data, sc_name='',
                   colors_dict=REGION_COLORS,
                   x_label=r'$\beta_\parallel$', y_label=r'$T_{\perp}/T_{\parallel}$',
                   title=None, nbins=50, mincount=10,
                   xrange=[-3, 2], yrange=[-1, 1],
                   contour_nlevels=5, norm=False,
                   mean_marker_size=120,
                   cgl_line=True, cgl_exponent=-1,
                   ax=None, return_centroids=False):
    """
    region_data : {region_name: (betapar, Tparperp, mask)}; mask may be None.
                  Tparperp is inverted internally (y = 1/Tparperp), as in brazil_multi.
    """
    # From Hellinger (2011)
    mirror_params      = [0.77, 0.76, -0.016]
    firehose_params    = [-1.4, 1.0, -0.11]
    cyclotron_params   = [0.43, 0.42, -0.0004]
    parfirehose_params = [-0.47, 0.53, 0.59]

    beta_par_range = np.logspace(xrange[0], xrange[1], 100)
    mirror_threshold      = stability_condition(beta_par_range, *mirror_params)
    firehose_threshold    = stability_condition(beta_par_range, *firehose_params)
    cyclotron_threshold   = stability_condition(beta_par_range, *cyclotron_params)
    parfirehose_threshold = stability_condition(beta_par_range, *parfirehose_params)

    xbins = np.logspace(xrange[0], xrange[1], nbins)
    ybins = np.logspace(yrange[0], yrange[1], nbins)
    xcenters = 0.5 * (xbins[:-1] + xbins[1:])
    ycenters = 0.5 * (ybins[:-1] + ybins[1:])
    xline_cgl = np.logspace(xrange[0], xrange[1], 100)
    xgrid, ygrid = np.meshgrid(xcenters, ycenters, indexing='ij')

    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 8))
    else:
        fig = ax.figure

    centroids = {}

    for region, (betapar, Tparperp, mask) in region_data.items():
        if betapar is None or Tparperp is None:
            continue
        if mask is None:
            mask = Tparperp > 0

        x = betapar[mask]
        y = 1 / Tparperp[mask]
        valmask = np.isfinite(x) & np.isfinite(y)
        x_valid, y_valid = x[valmask], y[valmask]

        color = colors_dict.get(region, 'k')

        counts, _, _, _ = binned_statistic_2d(
            x_valid, y_valid, None, statistic='count', bins=[xbins, ybins])

        if counts.max() <= mincount:
            continue

        if norm:
            total = x_valid.size
            plot_counts   = counts / total
            plot_mincount = mincount / total
            plot_max      = counts.max() / total
        else:
            plot_counts   = counts
            plot_mincount = mincount
            plot_max      = counts.max()

        levels = np.logspace(np.log10(plot_mincount), np.log10(plot_max), contour_nlevels)
        ax.contour(xcenters, ycenters, plot_counts.T, levels=levels,
                   colors=color, linewidths=0.7, alpha=1)

        # counts-weighted centroid of the binned distribution
        total_counts = counts.sum()
        if total_counts > 0:
            x_mean = np.sum(xgrid * counts) / total_counts
            y_mean = np.sum(ygrid * counts) / total_counts
            centroids[region] = (x_mean, y_mean)
            ax.scatter(x_mean, y_mean, marker='^', color=color, s=2*mean_marker_size,
                       edgecolors=color, linewidths=1.5, zorder=10,
                       label=region.capitalize())

            if cgl_line:
                yline_cgl = y_mean * (xline_cgl / x_mean) ** cgl_exponent
                ax.plot(xline_cgl, yline_cgl, color=color, linestyle='-',
                        linewidth=1, alpha=0.8, zorder=9)

    ax.set_facecolor('lightgray')
    ax.set_xlabel(x_label, fontdict={'size': 20})
    ax.set_ylabel(y_label, fontdict={'size': 20})
    ax.set_title(sc_name if title is None else title, fontdict={'size': 25})
    ax.set_xlim(10**(xrange[0]), 10**(xrange[1]))
    ax.set_ylim(10**(-1), 10**(1))
    ax.axhline(y=1, color='k', linestyle='-')
    ax.axvline(x=1, color='k', linestyle='-')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.tick_params(labelsize=15)
    ax.plot(beta_par_range, mirror_threshold,      label='Mirror',      color='k', linestyle='dotted', linewidth=2)
    ax.plot(beta_par_range, firehose_threshold,    label='Oblique FH',  color='k', linestyle='--',     linewidth=2)
    ax.plot(beta_par_range, cyclotron_threshold,   label='IC',          color='k', linestyle='dotted', linewidth=1)
    ax.plot(beta_par_range, parfirehose_threshold, label='Parallel FH', color='k', linestyle='--',     linewidth=1)
    ax.legend(fontsize=13, framealpha=1, loc='upper right')

    if return_centroids:
        return fig, ax, centroids
    return fig, ax


# ---------------------------------------------------------------
# Static: one plot per spacecraft (or 1xN panel)
# ---------------------------------------------------------------
def brazil_by_spacecraft(sc_region_data, one_figure=False, **kwargs):
    out = {}
    if one_figure:
        n = len(sc_region_data)
        fig, axes = plt.subplots(1, n, figsize=(9*n, 8), sharey=True)
        axes = np.atleast_1d(axes)
        for ax, (sc, rd) in zip(axes, sc_region_data.items()):
            brazil_regions(rd, sc_name=sc, ax=ax, **kwargs)
            out[sc] = (fig, ax)
        for ax in axes[1:]:
            ax.set_ylabel('')
        fig.tight_layout()
    else:
        for sc, rd in sc_region_data.items():
            out[sc] = brazil_regions(rd, sc_name=sc, **kwargs)
    plt.show()
    return out


# ---------------------------------------------------------------
# Movie: PSP -> SOLO -> WIND with growing centroid trails
# ---------------------------------------------------------------
def brazil_movie_trails(sc_region_data, order=('PSP', 'SOLO', 'WIND'),
                        labels=None, filename='brazil_trails.gif',
                        fps=1, dpi=150, colors_dict=REGION_COLORS,
                        trail_lw=1.0, trail_alpha=0.9, **kwargs):
    if labels is None:
        labels = {'PSP': 'PSP (~0.05 AU)',
                  'SOLO': 'SOLO (~0.3 AU)',
                  'WIND': 'WIND (~1 AU)'}

    # Pass 1: compute centroids for every spacecraft on a scratch figure
    cents = {}
    for sc in dict.fromkeys(order):                 # unique, order-preserving
        f, _, c = brazil_regions(sc_region_data[sc], colors_dict=colors_dict,
                                 return_centroids=True, **kwargs)
        plt.close(f)
        cents[sc] = c

    # Pass 2: animate
    fig, ax = plt.subplots(figsize=(10, 8))

    def update(i):
        ax.clear()
        sc = order[i]
        brazil_regions(sc_region_data[sc], sc_name=labels.get(sc, sc),
                       colors_dict=colors_dict, ax=ax, **kwargs)

        visited = list(dict.fromkeys(order[:i + 1]))   # spacecraft seen so far
        regions = {r for s in visited for r in cents[s]}
        for r in regions:
            pts = [cents[s][r] for s in visited if r in cents[s]]
            xs, ys = zip(*pts)
            color = colors_dict.get(r, 'k')
            if len(pts) > 1:
                ax.plot(xs, ys, '-', color=color, lw=trail_lw,
                        alpha=trail_alpha, zorder=8)
                # hollow triangles mark earlier positions
                ax.scatter(xs[:-1], ys[:-1], marker='^', facecolors='none',
                           edgecolors=color, s=70, linewidths=1.5, zorder=9)

    anim = FuncAnimation(fig, update, frames=len(order), blit=False)
    anim.save(filename, writer=PillowWriter(fps=fps), dpi=dpi)
    plt.close(fig)
    return anim


# ---------------------------------------------------------------
# Your data, organized by spacecraft -> region
# ---------------------------------------------------------------
sc_region_data = {
    'PSP': {
        'leading':     (betaparspsp_preceding,   Tparperppsp_preceding,   None),
        'compression': (betaparspsp_shoulder,    Tparperppsp_shoulder,    None),
        'fast':        (betaparspsp_fast,        Tparperppsp_fast,        None),
        'rarefaction': (betaparspsp_rarefaction, Tparperppsp_rarefaction, None),
        'trailing':    (betaparspsp_trailing,    Tparperppsp_trailing,    None),
    },
    'SOLO': {
        'leading':     (betaparssolo_preceding,   Tparperpsolo_preceding,   None),
        'compression': (betaparssolo_shoulder,    Tparperpsolo_shoulder,    None),
        'fast':        (betaparssolo_fast,        Tparperpsolo_fast,        None),
        'rarefaction': (betaparssolo_rarefaction, Tparperpsolo_rarefaction, None),
        'trailing':    (betaparssolo_trailing,    Tparperpsolo_trailing,    None),
    },
    'WIND': {
        'leading':     (betaparswind_preceding,   Tparperpwind_preceding,   None),
        'compression': (betaparswind_shoulder,    Tparperpwind_shoulder,    None),
        'fast':        (betaparswind_fast,        Tparperpwind_fast,        None),
        'rarefaction': (betaparswind_rarefaction, Tparperpwind_rarefaction, None),
        'trailing':    (betaparswind_trailing,    Tparperpwind_trailing,    None),
    },
}

# Fast CH wind + the two slow wind bins only
keep = ['fast']
sc_region_data_sub = {
    sc: {r: v for r, v in rd.items() if r in keep}
    for sc, rd in sc_region_data.items()
}

# ---------------------------------------------------------------
# Make the movie (last frame held an extra beat)
# ---------------------------------------------------------------
anim = brazil_movie_trails(sc_region_data_sub,
                           order=('PSP', 'SOLO', 'WIND', 'WIND'),
                           filename='brazil_trails.gif',
                           fps=0.6,
                           contour_nlevels=10,
                           mean_marker_size=35,
                           cgl_line=False,labels = {'PSP': 'Below $CL \ (M_A<1)$','SOLO': 'Above $CL$ (R_A>1)','WIND': 'L1'})

# Optional: static panels instead of / in addition to the movie
# figs = brazil_by_spacecraft(sc_region_data_sub, one_figure=True,
#                             contour_nlevels=10, mean_marker_size=35, cgl_line=False)

# Optional: view inline in a notebook
# from IPython.display import Image
# Image(filename='brazil_fast_slow_trails.gif')

anim.save('brazil_fastonly.gif', writer=PillowWriter(fps=0.6))
Image(filename='brazil_fastonly.gif')

# %%
