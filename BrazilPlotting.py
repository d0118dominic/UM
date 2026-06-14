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

def brazilhist(z, x_label=r'$\beta_\parallel$', y_label=r'$T_{\parallel}/T_{\perp}$',
     z_label='', title='',nbins=50,mincount=10,vmin=1,vmax=10000, mask = allTparperp>0,count=False,scale = 'linear',fitline=False):    
    mirror_params = [0.77, 0.76, -0.016]  # [a, b, beta0]
    firehose_params = [-1.4, 1.0, -0.11]  # [a, b, beta0]
    cyclotron_params = [0.45, 0.4, -0.0004]  # [a, b, beta0]
    parfirehose_params = [-0.47, 0.53, 0.59]  # [a, b, beta0]

    beta_par_range = np.logspace(-3, 3, 100)  # Adjust range and resolution as needed
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

    # x_fit = np.log10(x_valid)
    # y_fit = np.log10(y_valid)
    # m,c = np.polyfit(x_fit, y_fit, 1)


    nbins = nbins
    mincount=mincount
    # xbins = np.linspace(0.001, 100, nbins)
    # ybins = np.linspace(0.1, 10, nbins)
    xbins = np.logspace(-3, 3, nbins)
    ybins = np.logspace(-1, 1, nbins)

    # xline = np.logspace(-3, 3, 200)
    # yline = 10**(m*np.log10(xline) + c)

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

    im = ax.pcolormesh(xedges, yedges, hist_masked.T, cmap='jet', norm=norm)

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
    ax.set_xlim(1e-3, 10**(3))  # Adjust as needed
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





# brazilplot(1e3*allbeta_par, r'$\beta_\parallel$', r'$T_{\perp}/T_{\parallel}$', r'$v_x$'

nratio = abs(alln)/np.nanmean(abs(alln))
Tratio = abs(allT)/np.nanmean(abs(allT))
vratio = abs(allvmags)/np.nanmean(abs(allvmags))
nTratio = nratio*Tratio

bulkratio = abs(alln)*allvmags**2/np.nanmean(abs(alln)*allvmags**2)

brazilhist(abs(alln)/np.nanmean(abs(alln)), z_label = r'$n/\langle n \rangle$',fitline=False,nbins=80,vmin=0,vmax=2,count=False,scale='linear')
brazilhist(abs(allT)/np.nanmean(abs(allT)), z_label = r'$T/\langle T \rangle$',fitline=False,nbins=80,vmin=0,vmax=2,count=False,scale='linear')
brazilhist(abs(allvmags)/np.nanmean(abs(allvmags)), z_label = r'$v/\langle v \rangle$',fitline=False,nbins=80,vmin=0,vmax=2,count=False,scale='linear')
brazilhist(nratio*Tratio, z_label = r'$nT/\langle nT \rangle$',fitline=False,nbins=80,vmin=0,vmax=2,count=False,scale='linear')
brazilhist(bulkratio, z_label = r'$nv^2/\langle nv^2 \rangle$',fitline=False,nbins=80,vmin=0,vmax=2,count=False,scale='linear')


#%%
# brazilhist(1.602e19*abs(allT), z_label = r'$T \ (eV)$',fitline=False,nbins=60,vmin=00,vmax=200,count=False,scale='linear')
brazilhist(1e-3*abs(allvmags), z_label = r'$v \ (km/s)$',fitline=False,nbins=60,vmin=400,vmax=800,count=False,scale='linear')
# brazilhist(1e-6*abs(alln), z_label = r'$n \ (cm^{-3})$',fitline=False,nbins=60,vmin=0,vmax=10,count=False,scale='linear')

#%%

brazilhist(1e-6*abs(allT)*abs(alln), z_label = r'$n_iT_i \ (J/cm^3)$',fitline=False,nbins=80,vmin=0,vmax=1e-15,count=False,scale='linear')
brazilhist(1e-6*0.5*abs(alln)*mi*allvmags**2, z_label = r'$\frac{1}{2}nm_i |v_i|^2 \ (J/cm^3)$',fitline=False,nbins=80,vmin=0,vmax=5e-14,count=False,scale='linear')
brazilhist((1e-6*0.5*allBmags**2)/mu0, z_label = r'$\frac{|B|^2}{2\mu_0} \ (J/cm^3)$',fitline=False,nbins=80,vmin=0,vmax=1e-15,count=False,scale='linear')
Ratio = (1e-6*abs(allT)*abs(alln)) / (1e-6*0.5*abs(alln)*mi*allvmags**2)
brazilhist(Ratio, z_label = r'$\frac{n_iT_i}{\frac{1}{2}nm_i |v_i|^2}$',fitline=False,nbins=80,vmin=0,vmax=0.03,count=False,scale='linear')
# brazilhist((1e-6*0.5*allBmags**2)/mu0, z_label = r'$\frac{1}{2}n|B|^2 \ (J/cm^3)$',fitline=False,nbins=80,vmin=0,vmax=1e-15,count=False,scale='linear')
# %%
brazilhist(allpositions, z_label = r'$R/R_s$',nbins=80,vmin=0,vmax=40,count=False)
# %%
brazilhist(abs(allcrosshelicity), z_label = r'$n$',nbins=200,vmin = 0,vmax=1,count=False)
# %%
