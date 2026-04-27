#%%

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import binned_statistic_2d
from matplotlib.colors import ListedColormap
import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.colors as mcolors


me = 9.1094e-31 #kg
mi = 1837*me
mu0 = 1.2566370e-06  #;m kg / C^2
eps0 = 8.85e-12   # C^2/Nm^2

e = 1.602e-19 #C
Z = 1 # 1 for H+, 2 for He2+
gamma = 5/3
kb = 1.380649e-23

def stability_condition(betapar, a, b, beta0):
    denom = (betapar - beta0)**b
    Tperppar = 1 + a/denom
    return Tperppar

# Define parameters for different instabilities
mirror_params = [0.77, 0.76, -0.016]  # [a, b, beta0]
firehose_params = [-1.4, 1.0, -0.11]  # [a, b, beta0]
cyclotron_params = [0.45, 0.4, -0.0004]  # [a, b, beta0]
parfirehose_params = [-0.47, 0.53, 0.59]  # [a, b, beta0]

# Create beta_parallel array (logarithmically spaced for smoother curves)
betapar = np.logspace(-3, 2, 10000)  # From 0.001 to 1000

# Calculate T_perp/T_par for both instabilities
Tperppar_m = stability_condition(betapar, mirror_params[0], mirror_params[1], mirror_params[2])
Tperppar_f = stability_condition(betapar, firehose_params[0], firehose_params[1], firehose_params[2])
Tperppar_pf = stability_condition(betapar, parfirehose_params[0], parfirehose_params[1], parfirehose_params[2])
Tperppar_c = stability_condition(betapar, cyclotron_params[0], cyclotron_params[1], cyclotron_params[2])

log_betapar = np.log10(betapar)
log_Tperppar_m = np.log10(Tperppar_m)
log_Tperppar_f = np.log10(Tperppar_f)
log_Tperppar_pf = np.log10(Tperppar_pf)
log_Tperppar_c = np.log10(Tperppar_c)



# Create the plot
fig, ax = plt.subplots(figsize=(8, 5))

# Plot the stability conditions
# Plot the stability condition curves
ax.plot(log_betapar, log_Tperppar_m, 'w-', linewidth=3, 
            label='Mirror Instability', linestyle='dotted',color='k')
ax.plot(log_betapar, log_Tperppar_c, 'w-', linewidth=1.5, 
        label='IC Instability', linestyle='dotted',color='k')
ax.plot(log_betapar, log_Tperppar_f, 'w-', linewidth=3, 
        label='Oblique Firehose Instability', linestyle='--',color='k')
ax.plot(log_betapar, log_Tperppar_pf, 'w-', linewidth=1.5, 
        label='Parallel Firehose Instability', linestyle='--',color='k')

# Add reference lines at T_perp/T_par = 1 and beta_par = 1
ax.axhline(y=0, color='k', linestyle='--', linewidth=1, alpha=0.5)
ax.axvline(x=0, color='k', linestyle='--', linewidth=1, alpha=0.5)

# Set logarithmic scales
# ax.set_xscale('log')
# ax.set_yscale('log')

# Set axis limits
ax.set_xlim(-3, 2)
ax.set_ylim(-1, 1)  

ax.set_xlabel(r'$Log_{10}(\beta_{\parallel})$', fontsize=14)
ax.set_ylabel(r'$Log_{10}(T_{\perp}/T_{\parallel})$', fontsize=14)
ax.set_title('Plasma Instability Conditions', fontsize=16)

# Add legend
ax.legend(loc='best', fontsize=10)

# Improve layout
plt.tight_layout()
plt.show()
#%%




#%%
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from scipy.stats import binned_statistic_2d


def stability_condition(betapar, a, b, beta0):
    """
    Calculate T_perp/T_par stability boundary based on beta_parallel and instability parameters.
    From Bale et al. 2009: T_perp/T_par = 1 + a/(beta_par - beta0)^b
    """
    denom = (betapar - beta0)**b
    Tperppar = 1 + a/denom
    return Tperppar

def cgl_condition(betapar, betapar0,Tperppar0):
    # Condition 1: Tperp/B = const
    # Condition 2: TparB^2/n = const
    return Tperppar



def Instability_plot(allbeta_par, allTparperp, alldB_norm_mag):
    """
    Plot 2D histogram of plasma data with stability condition boundaries overlaid.
    """
    # Prepare data for histogram

    innerlim = 0
    outerlim = 60
    mask = ((allpositions <= outerlim) & (allpositions>=innerlim))
#     mask = ((allvmags >= 200e3) & (allvmags<=300e3))
    #mask = ((allvmags >= 100e3) & (allvmags<=200e3) & 
        #(allpositions <= outerlim) & (allpositions>=innerlim))
    mask = ((allpositions>=0)&(allpositions<=70)&(allvmags>=000e3)&(allvmags<600e3) & (allmachs>=0))

    x = np.log10(allbeta_par[mask])
    y = np.log10(1/allTparperp[mask])
    z = abs(alldv_norm_mag[mask] * allvmags[mask] / allva[mask])
    z = np.log10(allmachs[mask])
#     z = r_norm[mask]**(2)*allTperp[mask]/1.602e-19 # variable for colorbar
#     z = r_norm[mask]**(2)*1/allTparperp[mask] # variable for colorbar
    #z = 1000*allentropy[mask]/(1.602e-19) #puts it back in ev/m^2 # variable for colorbar
    #z = 1000*alldvpar_norm[mask] #puts it back in ev/m^2 # variable for colorbar
#     z = abs(allcrosshelicity[mask]) # variable for colorbar
    #z = allresidenergy[mask] # variable for colorbar


    # z = abs(alldv_norm_mag)
    # z = (r_norm**(2))*alln/1e6 # variable for colorbar
    # z = (r_norm**10/3)*(alln/1e6)*(allT/1.602e-19) # variable for colorbar
    #z = (1e9*allBmags**2)*(allTpar/1.602e-19)/((alln/1e6)**2) # variable for colorbar
    #z = (r_norm**(4/3))*(allTperp/1.602e-19) # variable for colorbar
    # z = 1000*allentropy/(1.602e-19)
    # z = (r_norm**2)*allBmags**2
    #z = (r_norm**(-4/3))*1000*allentropy/(1.602e-19) #puts it back in ev/m^2 # variable for colorbar
    #z = abs(allTperp)/abs(allT)
    #z = ((1e9*alldBmag)**2)/(alln/1e6)
    #z = allKr/allSr
    #z = ((1e9*allBmags)**2)/(alln/1e6) #Mag energy per part
    #z = allangles[mask]

    # Filter out invalid values (inf, -inf, nan)
    valid_mask = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    x_valid = x[valid_mask]
    y_valid = y[valid_mask]
    z_valid = z[valid_mask]
    
    print(f"Filtered out {np.sum(~valid_mask)} invalid points out of {len(x)} total")
    
    # Define bins
    bins = 90
    mincount = 100
    x_bins = np.linspace(-4, 4, bins)
    y_bins = np.linspace(-1, 1, bins)
    
    # Calculate the mean/max of 'z' for each bin
    mean_z, x_edge, y_edge, bin_number = binned_statistic_2d(
        x_valid, y_valid, z_valid, statistic='mean', bins=[x_bins, y_bins]  # Now use 'mean'
    )
    
    # mean_z, x_edge, y_edge, bin_number = binned_statistic_2d(
    #     x_valid, y_valid, z_valid, statistic=lambda v:np.nanstd(v)/np.sqrt(np.sum(~np.isnan(v))), bins=[x_bins, y_bins]  # Now use 'mean'
    # )
    # Calculate counts to identify zero bins
    count_z, _, _, _ = binned_statistic_2d(
        x_valid, y_valid, z_valid, statistic='count', bins=[x_bins, y_bins]
    )
    
    # Create a masked array where zero-count bins are masked
    mean_z_masked = np.ma.masked_where(count_z <=mincount, mean_z)

    # Also mask any remaining invalid values
    mean_z_masked = np.ma.masked_invalid(mean_z_masked)
    
    # Create figure
    fig, ax = plt.subplots(figsize=(10, 7))
    
    # First, plot the zero-count bins in light grey
    zero_mask = count_z <= mincount
    if np.any(zero_mask):
        zero_array = np.full_like(mean_z, 1.0)
        zero_array[~zero_mask] = np.nan
        
        ax.imshow(zero_array.T, origin='lower', 
                extent=[x_edge[0], x_edge[-1], y_edge[0], y_edge[-1]], 
                cmap='gray', aspect='auto', vmin=0.8, vmax=1.2, alpha=0.4)
    


    norm = colors.LogNorm(vmin=1, vmax=1e-7)
    norm = colors.LogNorm(vmin=1e-1, vmax=1e1)
    disccolors = ['b','w','r']
    bounds = [-0.5,-0.2,0.2,0.5]
    cmap = mcolors.ListedColormap(disccolors)
    discnorm = mcolors.BoundaryNorm(bounds,cmap.N)
    # Plot the 2D histogram
    im = ax.imshow(mean_z_masked.T, origin='lower', 
                extent=[x_edge[0], x_edge[-1], y_edge[0], y_edge[-1]], 
                cmap=cmap, aspect='auto',norm = discnorm)
    # Add colorbar
    cbar = plt.colorbar(im, ax=ax)
    cbar.ax.tick_params(labelsize=15)
    cbar.set_label(r'$\delta v / v_a$', fontsize=15)
    # cbar.set_label(r'$\delta B$', fontsize=12)
    cbar.set_label(r'$\theta_{def}$', fontsize=15)
    # cbar.set_label(r'$\theta_{def}$', fontsize=12)
    #cbar.set_label(r'$R^{10/3} \  nT \  (ev \ cm^{-1})$', fontsize=12)
    cbar.set_label(r'$R^{4/3} \ T_{} \ (ev)$', fontsize=15)
    #cbar.set_label(r'$V_{p} \ (km/s)$', fontsize=12)
    #cbar.set_label(r'$R^2 nT $', fontsize=12)
    #cbar.set_label(r'$R^{-4/3} \ S_p  \ (keV/m^2)$', fontsize=12)
    #cbar.set_label(r'$Log_{10}(M_a)$', fontsize=12)
#     cbar.set_label(r'$R/R_{sun}$', fontsize=12)
    # cbar.set_label(r'$Log_{10}(M_a)$', fontsize=15)
    cbar.set_label(r'$|\sigma_c|$', fontsize=17)
    
    # Now overlay the stability condition curves
    # Define parameters for different instabilities
    # These all come from Hellinger 06, and are used in Bale 09
    mirror_params = [0.77, 0.76, -0.016]  # [a, b, beta0]
    firehose_params = [-1.4, 1.0, -0.11]  # [a, b, beta0]
    cyclotron_params = [0.45, 0.4, -0.0004]  # [a, b, beta0]
    parfirehose_params = [-0.47, 0.53, 0.59]  # [a, b, beta0]
    # cyclotron_params = [0.62, 0.41, -0.0002]  # [a, b, beta0] #from coello-guzman manuscript
    # parfirehose_params = [-0.74, 0.37, 0.91]  # [a, b, beta0] #from corllo-guzman manuscript
    
    
    
    # Create beta_parallel array (logarithmically spaced)
    betapar = np.logspace(-3, 4, 10000)  # From 0.001 to 100
    
    # Calculate T_perp/T_par for both instabilities
    Tperppar_m = stability_condition(betapar, mirror_params[0], mirror_params[1], mirror_params[2])
    Tperppar_f = stability_condition(betapar, firehose_params[0], firehose_params[1], firehose_params[2])
    Tperppar_pf = stability_condition(betapar, parfirehose_params[0], parfirehose_params[1], parfirehose_params[2])
    Tperppar_c = stability_condition(betapar, cyclotron_params[0], cyclotron_params[1], cyclotron_params[2])
    
    # Convert to log10 space to match the histogram axes
    log_betapar = np.log10(betapar)
    log_Tperppar_m = np.log10(Tperppar_m)
    log_Tperppar_f = np.log10(Tperppar_f)
    log_Tperppar_pf = np.log10(Tperppar_pf)
    log_Tperppar_c = np.log10(Tperppar_c)
    
    # Plot the stability condition curves
    ax.plot(log_betapar, log_Tperppar_m, 'w-', linewidth=3, 
            label='Mirror Instability', linestyle='dotted',color='k')
    ax.plot(log_betapar, log_Tperppar_c, 'w-', linewidth=1.5, 
             label='IC Instability', linestyle='dotted',color='k')
    ax.plot(log_betapar, log_Tperppar_f, 'w-', linewidth=3, 
            label='Oblique Firehose Instability', linestyle='--',color='k')
    ax.plot(log_betapar, log_Tperppar_pf, 'w-', linewidth=1.5, 
            label='Parallel Firehose Instability', linestyle='--',color='k')
    # # Lines for T_par,perp = T/2 conditions
    # ax.axhline(y=0.4,color='w',linewidth=5,linestyle='-')
    # ax.axhline(y=-0.6,color='w',linewidth=5,linestyle='-')
    
    # Add reference lines
    ax.axhline(y=0, color='k', linewidth=1.5)
    ax.axvline(x=0, color='k', linewidth=1.5)
    
    # Set limits and labels
    ax.set_xlim(-3, 2)
    ax.set_ylim(-1, 1)
    ax.tick_params(labelsize=15)
    ax.set_xlabel(r'$Log_{10}(\beta_\parallel)$', fontsize=15)
    ax.set_ylabel(r'$Log_{10}(T_\perp/T_\parallel)$', fontsize=15)
    # ax.set_xlabel(r'$Log_{10}(\ \frac{r^2}{\langle r \rangle^2} \ \beta_\parallel \ )$', fontsize=12)
    # ax.set_ylabel(r'$Log_{10}(\ \frac{\langle r \rangle^2}{r^2} \ (T_\perp/T_\parallel) \ \ )$', fontsize=12)
    
    # Add legend with white background for visibility
    ax.legend(loc='best', fontsize=11, facecolor='white', framealpha=0.8)
    
    plt.tight_layout()
    plt.show()
    
    return fig, ax

# Example usage (you would call this with your actual data):
fig, ax = Instability_plot(allbeta_par, allTparperp, alldB_norm_mag)

# %%

def stability_condition(betapar, a, b, beta0):
    """
    Calculate T_perp/T_par stability boundary based on beta_parallel and instability parameters.
    From Bale et al. 2009: T_perp/T_par = 1 + a/(beta_par - beta0)^b
    """
    denom = (betapar - beta0)**b
    Tperppar = 1 + a/denom
    return Tperppar


def cgl_curve(beta_par, beta_par_0=1.0, Tperppar_0=1.0,exponent=1.0):
   
    ratio = (beta_par_0 / beta_par)**exponent
    Tperppar = Tperppar_0 * ratio
    return Tperppar


def Instability_plot(allbeta_par, allTparperp, alldB_norm_mag):
    """
    Plot 2D histogram of plasma data with stability condition boundaries overlaid.
    """
    # Prepare data for histogram

    innerlim = 0
    outerlim = 50
    mask = ((allpositions <= outerlim) & (allpositions>=innerlim))
#     mask = ((allvmags >= 200e3) & (allvmags<=300e3))
    #mask = ((allvmags >= 100e3) & (allvmags<=200e3) & 
        #(allpositions <= outerlim) & (allpositions>=innerlim))
#     mask = ((allvmags>=400e3) & (allpositions>=55) & (allpositions<=700))

    mask = ((allpositions>=25)&(allpositions<=40)&(allvmags>=400e3)&(allvmags<=600e3))
    mask = ((np.log10(allmachs)>=-1) & (np.log10(allmachs)<=-0.2) & (allpositions<=40) &(allvmags>=300e3) & (allvmags<=600e3))
    mask1 = ((allangles>=90) & (np.log10(allmachs)>=-0.2) & (np.log10(allmachs)<=0.2) & (allpositions<=40) &(allvmags>=400e3) & (allvmags<=600e3))
#     mask1 = ((allangles>=90) & (allpositions>=25) & (allpositions<=40) &(allvmags>=00e3) & (allvmags<=600e3))
    mask2 = ((allangles>=90) & (allpositions>=30) & (allpositions<=70) &(allvmags>=00e3) & (allvmags<=600e3))

    x = np.log10(allbeta_par[mask])
    y = np.log10(1/allTparperp[mask])
    # z = abs(alldv_norm_mag[mask] * allvmags[mask] / allva[mask])
    z = allpositions[mask]
    # z = r_norm[mask]**(0)*allT[mask]/1.602e-19 # variable for colorbar
    # z = r_norm[mask]**(2)*allTperp[mask]/allTpar[mask]# variable for colorbar
#     z = np.log10(allmachs[mask])
#     z = abs(allcrosshelicity[mask])
    # z = (r_norm[mask]**2)*allSmags[mask]
    # z = (r_norm[mask]**2)*allKr[mask]
    # z = allmachs[mask]

    # Filter out invalid values (inf, -inf, nan)
    valid_mask = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    x_valid = x[valid_mask]
    y_valid = y[valid_mask]
    z_valid = z[valid_mask]
    
    print(f"Filtered out {np.sum(~valid_mask)} invalid points out of {len(x)} total")
    
    # Define bins
    bins = 90
    mincount=0
    x_bins = np.linspace(-3, 1, bins)
    y_bins = np.linspace(-1, 1, bins)
    
    # Calculate the mean/max of 'z' for each bin
    mean_z, x_edge, y_edge, bin_number = binned_statistic_2d(
        x_valid, y_valid, z_valid, statistic='count', bins=[x_bins, y_bins]  # Now use 'mean'


    )
    


    # Calculate counts to identify zero bins
    count_z, _, _, _ = binned_statistic_2d(
        x_valid, y_valid, z_valid, statistic='count', bins=[x_bins, y_bins]
    )
    
    # Create a masked array where zero-count bins are masked
    mean_z_masked = np.ma.masked_where(count_z <=mincount, mean_z)

    # Also mask any remaining invalid values
    mean_z_masked = np.ma.masked_invalid(mean_z_masked)
    
    # Create figure
    fig, ax = plt.subplots(figsize=(10, 7.5))
    
    # First, plot the zero-count bins in light grey
    zero_mask = count_z <= mincount
    if np.any(zero_mask):
        zero_array = np.full_like(mean_z, 1.0)
        zero_array[~zero_mask] = np.nan
        
        ax.imshow(zero_array.T, origin='lower', 
                extent=[x_edge[0], x_edge[-1], y_edge[0], y_edge[-1]], 
                cmap='gray', aspect='auto', vmin=0.8, vmax=1.2, alpha=0.4)
    


#     norm = colors.LogNorm(vmin=1e-8, vmax=1e-6)
# #     norm = colors.LogNorm(vmin=1e7, vmax=1e10)
#     norm = colors.LogNorm(vmin=1e2, vmax=4e3)

    disccolors = ['b','w','r']
    bounds = [-0.5,-0.2,0.2,0.5]
    cmap = mcolors.ListedColormap(disccolors)
    discnorm = mcolors.BoundaryNorm(bounds,cmap.N)
    norm = colors.LogNorm(vmin=1e-1, vmax=1e4)
    # norm = colors.LogNorm(vmin=1e-2, vmax=1e-1)
    # norm = colors.LogNorm(vmin=1e2 , vmax=1e3)
    # norm = colors.Normalize(vmin=-1, vmax=1)
    im = ax.imshow(mean_z_masked.T, origin='lower', 
                extent=[x_edge[0], x_edge[-1], y_edge[0], y_edge[-1]], 
                # cmap=cmap, aspect='auto',norm=discnorm)
                cmap='jet', aspect='auto',norm=norm)
    # Add colorbar
    cbar = plt.colorbar(im, ax=ax)
    cbar.ax.tick_params(labelsize=17)
    cbar.set_label(r'$\delta v / v_a$', fontsize=15)
    # cbar.set_label(r'$\delta B$', fontsize=12)
    cbar.set_label(r'$\theta_{def}$', fontsize=15)
    cbar.set_label(r'$\langle\theta \rangle$', fontsize=15)
    #cbar.set_label(r'$R^{10/3} \  nT \  (ev \ cm^{-1})$', fontsize=12)
    # cbar.set_label(r'$  \mathscr{R}^{4/3} T \ (eV)$', fontsize=17)
    cbar.set_label(r'$T_{\alpha} \ (eV)$', fontsize=17)
    cbar.set_label(r'$V_{\alpha} \ (km/s)$', fontsize=17)
    # cbar.set_label(r'$dB \  Energy \ per \ particle \ (eV) $', fontsize=12)
    # cbar.set_label(r'$R^{4/3} \ S_p  \ (keV/m^2)$', fontsize=12)
    cbar.set_label(r'$Log_{10}(M_a)$', fontsize=17)
    # cbar.set_label(r'$\langle\theta\rangle$', fontsize=18)
    # cbar.set_label(r'$ R/R_s$', fontsize=15)
    # cbar.set_label(r'$R/R_s$', fontsize=17)
    # cbar.set_label(r'$\mathscr{R}^2 \ |S|$', fontsize=17)
    cbar.set_label(r'Counts', fontsize=17)
#     cbar.set_label(r'|$\sigma_c$|', fontsize=17)

    # Now overlay the stability condition curves
    # Define parameters for different instabilities
    # These all come from Hellinger 06, and are used in Bale 09
    mirror_params = [0.77, 0.76, -0.016]  # [a, b, beta0]
    firehose_params = [-1.4, 1.0, -0.11]  # [a, b, beta0]
    cyclotron_params = [0.45, 0.4, -0.0004]  # [a, b, beta0]
    parfirehose_params = [-0.47, 0.53, 0.59]  # [a, b, beta0]
    # cyclotron_params = [0.62, 0.41, -0.0002]  # [a, b, beta0] #from coello-guzman manuscript
    # parfirehose_params = [-0.74, 0.37, 0.91]  # [a, b, beta0] #from corllo-guzman manuscript
    
    
    
    # Create beta_parallel array (logarithmically spaced)
    betapar = np.logspace(-3, 1, 10000)  # From 0.001 to 100
    
    # Calculate T_perp/T_par for both instabilities
    Tperppar_m = stability_condition(betapar, mirror_params[0], mirror_params[1], mirror_params[2])
    Tperppar_f = stability_condition(betapar, firehose_params[0], firehose_params[1], firehose_params[2])
    Tperppar_pf = stability_condition(betapar, parfirehose_params[0], parfirehose_params[1], parfirehose_params[2])
    Tperppar_c = stability_condition(betapar, cyclotron_params[0], cyclotron_params[1], cyclotron_params[2])
    
    # Calculate CGL double adiabatic curves
    beta_par_0 = 0.02  # Reference beta_parallel
    Tperppar_0 = 1  # Reference T_perp/T_par (isotropic starting point)
    exp1 = 1
    exp2 = 0.5
    exp3 = 0.6
    exp4 = 1

    beta_par_0 = 10**-0.7  # Reference beta_parallel
    Tperppar_0 = 10**0.45  # Reference T_perp/T_par (isotropic starting point)
    exp1 = 1
#     exp1 = 0.4
#     exp2 = 0.5
#     exp3 = 0.6
#     exp4 = 1
   
    Tperppar_cgl1 = cgl_curve(betapar, beta_par_0, Tperppar_0,exp1)
    # Tperppar_cgl2 = cgl_curve(betapar, 10**-0., 1,1.3)
    Tperppar_cgl3 = cgl_curve(beta_par, 1, Tperppar_0,1)
    
    Tperppar_cgl2 = cgl_curve(betapar, 10**-1, Tperppar_0,1)
    Tperppar_cgl3 = cgl_curve(betapar,10**-0.5, Tperppar_0,1)
    Tperppar_cgl4 = cgl_curve(betapar,10**0, Tperppar_0,1)
    
#     Tperppar_cgl1 = cgl_curve(betapar, 10**-1.6,Tperppar_0,0.55)
    Tperppar_cgl2 = cgl_curve(betapar, 10**-1.6,Tperppar_0,0.55)
    # Tperppar_cgl1 = cgl_curve(betapar, 10**-0.65,Tperppar_0,0.55)
#     Tperppar_cgl2 = cgl_curve(betapar, beta_par_0, 10**0.1,0.2)
    
#     Tperppar_cgl4 = cgl_curve(betapar,beta_par_0, 10**0.1,0.3) #non-CGL
    
    # Convert to log10 space to match the histogram axes
    log_betapar = np.log10(betapar)
    log_Tperppar_m = np.log10(Tperppar_m)
    log_Tperppar_f = np.log10(Tperppar_f)
    log_Tperppar_pf = np.log10(Tperppar_pf)
    log_Tperppar_c = np.log10(Tperppar_c)
    log_Tperppar_cgl1 = np.log10(Tperppar_cgl1)
    log_Tperppar_cgl2 = np.log10(Tperppar_cgl2)
    log_Tperppar_cgl3 = np.log10(Tperppar_cgl3)
    log_Tperppar_cgl4 = np.log10(Tperppar_cgl4)
    
    # Plot the stability condition curves
    ax.plot(log_betapar, log_Tperppar_m, 'w-', linewidth=3, 
            label='Mirror', linestyle='dotted',color='k')
    ax.plot(log_betapar, log_Tperppar_c, 'w-', linewidth=1.5, 
             label='IC', linestyle='dotted',color='k')
    ax.plot(log_betapar, log_Tperppar_f, 'w-', linewidth=3, 
            label='Oblique FH', linestyle='--',color='k')
    ax.plot(log_betapar, log_Tperppar_pf, 'w-', linewidth=1.5, 
            label='Parallel FH', linestyle='--',color='k')
    
#     Plot CGL double adiabatic curves
    ax.plot(log_betapar, log_Tperppar_cgl1, linewidth=1, linestyle='-', color='blue',label=r'$\beta^{-1}_\parallel$')
 
#     ax.plot(log_betapar, log_Tperppar_cgl2, linewidth=1, 
#             label=r'$\beta_{\parallel}^{-0.5}$', linestyle='-', color='red')
    
    ax.plot(log_betapar, log_Tperppar_cgl2, linewidth=1, linestyle='-', color='r',label=r'$\beta^{-0.55}_\parallel$')
    # ax.plot(log_betapar, log_Tperppar_cgl3, linewidth=1, linestyle='-', color='b')
    # ax.plot(log_betapar, log_Tperppar_cgl4, linewidth=1, linestyle='-', color='b', label=r'$\beta^{-1}_\parallel$')

    # Add reference lines
    ax.axhline(y=0, color='k', linewidth=1.5)
    ax.axvline(x=0, color='k', linewidth=1.5)
    
    # Set limits and labels
    ax.set_xlim(-3, 1)
    ax.set_ylim(-1, 1)
    ax.tick_params(labelsize=17)
    ax.set_xlabel(r'$Log_{10}(\beta_\parallel)$', fontsize=17)
    ax.set_ylabel(r'$Log_{10}(T_\perp/T_\parallel)$', fontsize=17)
    # ax.set_xlabel(r'$Log_{10}(\ \frac{r^2}{\langle r \rangle^2} \ \beta_\parallel \ )$', fontsize=12)
    # ax.set_ylabel(r'$Log_{10}(\ \frac{\langle r \rangle^2}{r^2} \ (T_\perp/T_\parallel) \ \ )$', fontsize=12)
    
    # Add legend with white background for visibility



    # Scatterplot of small subset
    a = np.log10(allbeta_par)
    b = np.log10(1/allTparperp)
    sc = ax.scatter(a[mask1], (b[mask1]), marker='o', c='k', s=15,label=r'$\theta \geq 90^o$')
    x_fit = a[mask1]
    y_fit = b[mask1]
    # Clean data
    valid = np.isfinite(x_fit) & np.isfinite(y_fit)
    x_fit = x_fit[valid]
    y_fit = y_fit[valid]
    m, c = np.polyfit(x_fit, y_fit, 1)
    # Create fit line
    x_line = np.linspace(-3, 1, 200)
    y_line = m * x_line + c
    # Plot fit
    # ax.plot(x_line, y_line, 'r-', linewidth=3,
    #         label=rf'Best fit: $\beta^{{{m:.2f}}}$',color='b',linestyle='-')

    # ax.plot(log_betapar, log_Tperppar_cgl1, linewidth=1, linestyle='-', color='k',label=r'$\beta^{-1}_\parallel$')
    # ax.plot(log_betapar, log_Tperppar_cgl2, linewidth=1, linestyle='-', color='k',label=r'$\beta^{-0.2}_\parallel$')

    
    ax.legend(loc='upper right', fontsize=15, facecolor='white', framealpha=0.8)
    plt.tight_layout()
    plt.show()
    
    return fig, ax

fig, ax = Instability_plot(allTparperp, allTparperp, allTparperp)
# sc = plt.scatter(x[mask1],(y[mask1]),marker='*',c='r',s=3)



# %%

# %%

# %%

distlist = np.array([10,15,20,25,30,35,40,45,50,55,60])

wind1 = 00e3
wind2 = 600e3
msize=10

mask1 = ((allvmags>=wind1)& (allvmags<wind2)& (allpositions>=10) & (allpositions<=30))
mask2 = ((allvmags>=wind1)&(allvmags<=wind2)&(allpositions>=30) & (allpositions<=50))
mask3 = ((allvmags>=wind1)&(allvmags<=wind2)&(allpositions>=50) & (allpositions<=70))
# mask4 = ((allvmags>=wind1)&(allvmags<=wind2)&(allpositions>=50) & (allpositions<=70))
# mask5 = ((allpositions>=50) & (allpositions<=60))
masks = [mask1,mask2,mask3]

 # variable for colorbar
ch = []
re = []
ntpar = []
ntperp = []
nt = []
exp = 4/3
for i in range(len(masks)):
        ch.append(np.nanmean(abs(allcrosshelicity[masks[i]])))
        re.append(np.nanmean(allresidenergy[masks[i]]))
        ntpar.append(np.nanmean(r_norm[masks[i]]**(exp)*allTpar[masks[i]]/1.602e-19))
        ntperp.append(np.nanmean(r_norm[masks[i]]**(exp)*allTperp[masks[i]]/1.602e-19))
        nt.append(np.nanmean(r_norm[masks[i]]**(exp)*allT[masks[i]]/1.602e-19))

ch = np.array(ch)
re = np.array(re)
ntpar = np.array(ntpar)
ntperp = np.array(ntperp)
nt = np.array(nt)

#normalized
plt.plot(nt,marker='o',linestyle=' ',label = r'$\langle R^{4/3} T\rangle$',ms=msize)
plt.plot(ntpar,marker='o',linestyle=' ',label = r'$\langle R^{4/3} T_\parallel\rangle$',ms=msize)
plt.plot(ntperp,marker='o',linestyle=' ',label = r'$\langle R^{4/3} T_\perp\rangle$',ms=msize)
plt.legend(fontsize=15)
plt.ylim(0,210)
plt.xticks([0, 1, 2], ['$10-30 \ R_S$', '$30-50 \ R_S$', '$50-70 \ R_S$'])
# plt.xticks([0, 1, 2, 3], ['$10-15 \ R_S$', '$15-20 \ R_S$', '$20-25 \ R_S$', '$25-30 \ R_S$'])



fig, ax1 = plt.subplots()

ax1.plot(ch, marker='o',linestyle='', label=r'$\langle|\sigma_c|\rangle$',color='r',markersize=msize)
ax1.set_ylabel(r'$\langle|\sigma_c|\rangle$', color='r')
ax1.set_ylim(0.4,1)
ax1.tick_params(axis='y', labelcolor='r')

ax2 = ax1.twinx()
ax2.plot(re, marker='o',linestyle=' ', label=r'$\langle\sigma_r \rangle$',color='b',markersize=msize)
ax2.set_ylabel(r'$\langle\sigma_r \rangle$', color='b')
ax2.set_ylim(-0.3,0)
ax2.tick_params(axis='y', labelcolor='b')
plt.xticks([0, 1, 2], ['$10-30 \ R_S$', '$30-50 \ R_S$', '$50-70 \ R_S$'])

fig.tight_layout()
plt.show()

######



mask1 = ((allvmags>=wind1)& (allvmags<wind2)& (allpositions>=10) & (allpositions<=15))
mask2 = ((allvmags>=wind1)&(allvmags<=wind2)&(allpositions>=15) & (allpositions<=20))
mask3 = ((allvmags>=wind1)&(allvmags<=wind2)&(allpositions>=20) & (allpositions<=25))
mask4 = ((allvmags>=wind1)&(allvmags<=wind2)&(allpositions>=25) & (allpositions<=30))
# mask5 = ((allpositions>=50) & (allpositions<=60))
masks = [mask1,mask2,mask3,mask4]

 # variable for colorbar
ch = []
re = []
ntpar = []
ntperp = []
nt = []

for i in range(len(masks)):
        ch.append(np.nanmean(abs(allcrosshelicity[masks[i]])))
        re.append(np.nanmean(allresidenergy[masks[i]]))
        ntpar.append(np.nanmean(r_norm[masks[i]]**(exp)*allTpar[masks[i]]/1.602e-19))
        ntperp.append(np.nanmean(r_norm[masks[i]]**(exp)*allTperp[masks[i]]/1.602e-19))
        nt.append(np.nanmean(r_norm[masks[i]]**(exp)*allT[masks[i]]/1.602e-19))

ch = np.array(ch)
re = np.array(re)
ntpar = np.array(ntpar)
ntperp = np.array(ntperp)
nt = np.array(nt)
plt.plot(nt,marker='o',linestyle=' ',label = r'$\langle R^{4/3} T\rangle$',ms=msize)
plt.plot(ntpar,marker='o',linestyle=' ',label = r'$\langle R^{4/3} T_\parallel\rangle$',ms=msize)
plt.plot(ntperp,marker='o',linestyle=' ',label = r'$\langle R^{4/3} T_\perp\rangle$',ms=msize)
plt.legend(fontsize=15)
plt.ylim(0,100)
plt.xticks([0, 1, 2,3], ['$10-15 \ R_S$', '$15-20 \ R_S$', '$20-25 \ R_S$','$25-30 \ R_S$'])
# plt.xticks([0, 1, 2, 3], ['$10-15 \ R_S$', '$15-20 \ R_S$', '$20-25 \ R_S$', '$25-30 \ R_S$'])



fig, ax1 = plt.subplots()

ax1.plot(ch, marker='o',linestyle=' ', label=r'$\langle|\sigma_c|\rangle$',color='r',markersize=10)
ax1.set_ylabel(r'$\langle|\sigma_c|\rangle$', color='r')
ax1.set_ylim(0.4,1)
ax1.tick_params(axis='y', labelcolor='r')

ax2 = ax1.twinx()
ax2.plot(re, marker='o',linestyle=' ', label=r'$\langle\sigma_r \rangle$',color='b',markersize=msize)
ax2.set_ylabel(r'$\langle\sigma_r \rangle$', color='b')
ax2.set_ylim(-0.3,0)
ax2.tick_params(axis='y', labelcolor='b')
plt.xticks([0, 1, 2,3], ['$10-15 \ R_S$', '$15-20 \ R_S$', '$20-25 \ R_S$','$25-30 \ R_S$'])

fig.tight_layout()
plt.show()





# %%


#VERSION WITH ERROR BARS (std)
# distlist = np.array([10,15,20,25,30,35,40,45,50,55,60])

# wind1 = 400e3
# wind2 = 600e3
# msize=5

# mask1 = ((allvmags>=wind1)& (allvmags<wind2)& (allpositions>=10) & (allpositions<=30))
# mask2 = ((allvmags>=wind1)&(allvmags<=wind2)&(allpositions>=30) & (allpositions<=50))
# mask3 = ((allvmags>=wind1)&(allvmags<=wind2)&(allpositions>=50) & (allpositions<=70))
# masks = [mask1,mask2,mask3]

# ch, ch_err = [], []
# re, re_err = [], []
# ntpar, ntpar_err = [], []
# ntperp, ntperp_err = [], []
# nt, nt_err = [], []
# exp = 4/3

# for i in range(len(masks)):
#     ch.append(np.nanmean(abs(allcrosshelicity[masks[i]])))
#     ch_err.append(np.nanstd(abs(allcrosshelicity[masks[i]])))

#     re.append(np.nanmean(allresidenergy[masks[i]]))
#     re_err.append(np.nanstd(allresidenergy[masks[i]]))

#     vals_par = r_norm[masks[i]]**(exp)*allTpar[masks[i]]/1.602e-19
#     ntpar.append(np.nanmean(vals_par));   ntpar_err.append(np.nanstd(vals_par))

#     vals_perp = r_norm[masks[i]]**(exp)*allTperp[masks[i]]/1.602e-19
#     ntperp.append(np.nanmean(vals_perp)); ntperp_err.append(np.nanstd(vals_perp))

#     vals = r_norm[masks[i]]**(exp)*allT[masks[i]]/1.602e-19
#     nt.append(np.nanmean(vals));          nt_err.append(np.nanstd(vals))

# ch, ch_err = np.array(ch), np.array(ch_err)
# re, re_err = np.array(re), np.array(re_err)
# ntpar, ntpar_err = np.array(ntpar), np.array(ntpar_err)
# ntperp, ntperp_err = np.array(ntperp), np.array(ntperp_err)
# nt, nt_err = np.array(nt), np.array(nt_err)

# # --- Normalized temperature plot (wide bins) ---
# # plt.errorbar(range(len(nt)),    nt,    yerr=nt_err,    marker='o', linestyle=' ', label=r'$\langle R^{4/3} T\rangle$',          ms=msize, capsize=4)
# plt.errorbar(range(len(ntpar)), ntpar, yerr=ntpar_err, marker='o', linestyle=' ', label=r'$\langle R^{4/3} T_\parallel\rangle$', ms=msize, capsize=4)
# plt.errorbar(range(len(ntperp)),ntperp,yerr=ntperp_err,marker='o', linestyle=' ', label=r'$\langle R^{4/3} T_\perp\rangle$',    ms=msize, capsize=4)
# plt.legend(fontsize=15)
# plt.ylim(0,250)
# plt.xticks([0, 1, 2], ['$10-30 \ R_S$', '$30-50 \ R_S$', '$50-70 \ R_S$'])

# # --- Cross-helicity / residual energy plot (wide bins) ---
# fig, ax1 = plt.subplots()
# ax1.errorbar(range(len(ch)), ch, yerr=ch_err, marker='o', linestyle='', color='r',
#              markersize=msize, capsize=4, label=r'$\langle|\sigma_c|\rangle$')
# ax1.set_ylabel(r'$\langle|\sigma_c|\rangle$', color='r')
# ax1.set_ylim(0,2)
# ax1.tick_params(axis='y', labelcolor='r')

# ax2 = ax1.twinx()
# ax2.errorbar(range(len(re)), re, yerr=re_err, marker='o', linestyle=' ', color='b',
#              markersize=msize, capsize=4, label=r'$\langle\sigma_r \rangle$')
# ax2.set_ylabel(r'$\langle\sigma_r \rangle$', color='b')
# ax2.set_ylim(-1,1)
# ax2.tick_params(axis='y', labelcolor='b')
# plt.xticks([0, 1, 2], ['$10-30 \ R_S$', '$30-50 \ R_S$', '$50-70 \ R_S$'])
# fig.tight_layout()
# plt.show()

# ######

# mask1 = ((allvmags>=wind1)& (allvmags<wind2)& (allpositions>=10) & (allpositions<=15))
# mask2 = ((allvmags>=wind1)&(allvmags<=wind2)&(allpositions>=15) & (allpositions<=20))
# mask3 = ((allvmags>=wind1)&(allvmags<=wind2)&(allpositions>=20) & (allpositions<=25))
# mask4 = ((allvmags>=wind1)&(allvmags<=wind2)&(allpositions>=25) & (allpositions<=30))
# masks = [mask1,mask2,mask3,mask4]

# ch, ch_err = [], []
# re, re_err = [], []
# ntpar, ntpar_err = [], []
# ntperp, ntperp_err = [], []
# nt, nt_err = [], []

# for i in range(len(masks)):
#     ch.append(np.nanmean(abs(allcrosshelicity[masks[i]])))
#     ch_err.append(np.nanstd(abs(allcrosshelicity[masks[i]])))

#     re.append(np.nanmean(allresidenergy[masks[i]]))
#     re_err.append(np.nanstd(allresidenergy[masks[i]]))

#     vals_par = r_norm[masks[i]]**(exp)*allTpar[masks[i]]/1.602e-19
#     ntpar.append(np.nanmean(vals_par));   ntpar_err.append(np.nanstd(vals_par))

#     vals_perp = r_norm[masks[i]]**(exp)*allTperp[masks[i]]/1.602e-19
#     ntperp.append(np.nanmean(vals_perp)); ntperp_err.append(np.nanstd(vals_perp))

#     vals = r_norm[masks[i]]**(exp)*allT[masks[i]]/1.602e-19
#     nt.append(np.nanmean(vals));          nt_err.append(np.nanstd(vals))

# ch, ch_err = np.array(ch), np.array(ch_err)
# re, re_err = np.array(re), np.array(re_err)
# ntpar, ntpar_err = np.array(ntpar), np.array(ntpar_err)
# ntperp, ntperp_err = np.array(ntperp), np.array(ntperp_err)
# nt, nt_err = np.array(nt), np.array(nt_err)

# # --- Normalized temperature plot (narrow bins) ---
# plt.errorbar(range(len(nt)),    nt,    yerr=nt_err,    marker='o', linestyle=' ', label=r'$\langle R^{4/3} T\rangle$',          ms=msize, capsize=4)
# plt.errorbar(range(len(ntpar)), ntpar, yerr=ntpar_err, marker='o', linestyle=' ', label=r'$\langle R^{4/3} T_\parallel\rangle$', ms=msize, capsize=4)
# plt.errorbar(range(len(ntperp)),ntperp,yerr=ntperp_err,marker='o', linestyle=' ', label=r'$\langle R^{4/3} T_\perp\rangle$',    ms=msize, capsize=4)
# plt.legend(fontsize=15)
# plt.ylim(0,200)
# plt.xticks([0, 1, 2, 3], ['$10-15 \ R_S$', '$15-20 \ R_S$', '$20-25 \ R_S$', '$25-30 \ R_S$'])

# # --- Cross-helicity / residual energy plot (narrow bins) ---
# fig, ax1 = plt.subplots()
# ax1.errorbar(range(len(ch)), ch, yerr=ch_err, marker='o', linestyle=' ', color='r',
#              markersize=msize, capsize=4, label=r'$\langle|\sigma_c|\rangle$')
# ax1.set_ylabel(r'$\langle|\sigma_c|\rangle$', color='r')
# ax1.set_ylim(0.1,2)
# ax1.tick_params(axis='y', labelcolor='r')

# ax2 = ax1.twinx()
# ax2.errorbar(range(len(re)), re, yerr=re_err, marker='o', linestyle=' ', color='b',
#              markersize=msize, capsize=4, label=r'$\langle\sigma_r \rangle$')
# ax2.set_ylabel(r'$\langle\sigma_r \rangle$', color='b')
# ax2.set_ylim(-1,1)
# ax2.tick_params(axis='y', labelcolor='b')
# plt.xticks([0, 1, 2, 3], ['$10-15 \ R_S$', '$15-20 \ R_S$', '$20-25 \ R_S$', '$25-30 \ R_S$'])
# fig.tight_layout()
# plt.show()








# %%
# %%

# %%

#Scatterplotting

mirror_params = [0.77, 0.76, -0.016]  # [a, b, beta0]
firehose_params = [-1.4, 1.0, -0.11]  # [a, b, beta0]
cyclotron_params = [0.45, 0.4, -0.0004]  # [a, b, beta0]
parfirehose_params = [-0.47, 0.53, 0.59]  # [a, b, beta0]
# Create beta_parallel array (logarithmically spaced)
betapar = np.logspace(-1.3, 1, 10000)  # From 0.001 to 100
Tperppar_m = stability_condition(betapar, mirror_params[0], mirror_params[1], mirror_params[2])
Tperppar_f = stability_condition(betapar, firehose_params[0], firehose_params[1], firehose_params[2])
Tperppar_pf = stability_condition(betapar, parfirehose_params[0], parfirehose_params[1], parfirehose_params[2])
Tperppar_c = stability_condition(betapar, cyclotron_params[0], cyclotron_params[1], cyclotron_params[2])
# Calculate T_perp/T_par for both instabilities
log_betapar = np.log10(betapar)
log_Tperppar_m = np.log10(Tperppar_m)
log_Tperppar_f = np.log10(Tperppar_f)
log_Tperppar_pf = np.log10(Tperppar_pf)
log_Tperppar_c = np.log10(Tperppar_c)

# Calculate CGL double adiabatic curves
beta_par_0 = 10**-1.1  # Reference beta_parallel
Tperppar_0 = 10**0.45  # Reference T_perp/T_par (isotropic starting point)
exp1 = 1
exp2 = 0.3
exp3 = 0.6

Tperppar_cgl1 = cgl_curve(betapar, 10**-0.7, Tperppar_0,exp1)
Tperppar_cgl2 = cgl_curve(betapar, 10**0.45, 10**0.15,0.2)
Tperppar_cgl3 = cgl_curve(betapar, beta_par_0, 10**0.15,exp3)

# Convert to log10 space to match the histogram axes
log_betapar = np.log10(betapar)
log_Tperppar_m = np.log10(Tperppar_m)
log_Tperppar_f = np.log10(Tperppar_f)
log_Tperppar_pf = np.log10(Tperppar_pf)
log_Tperppar_c = np.log10(Tperppar_c)
log_Tperppar_cgl1 = np.log10(Tperppar_cgl1)
log_Tperppar_cgl2 = np.log10(Tperppar_cgl2)
log_Tperppar_cgl3 = np.log10(Tperppar_cgl3)
    
# Plot the stability condition curves
plt.plot(log_betapar, log_Tperppar_m, 'w-', linewidth=3, 
        label='Mirror Instability', linestyle='dotted',color='k')
plt.plot(log_betapar, log_Tperppar_c, 'w-', linewidth=1.5, 
        label='IC Instability', linestyle='dotted',color='k')
plt.plot(log_betapar, log_Tperppar_f, 'w-', linewidth=3, 
        label='Oblique FH', linestyle='--',color='k')
plt.plot(log_betapar, log_Tperppar_pf, 'w-', linewidth=1.5, 
        label='Parallel FH', linestyle='--',color='k')


# Plot CGL double adiabatic curves


# plt.plot(log_betapar, log_Tperppar_cgl3, linewidth=1, 
#         label=r'$\beta^{-0.6}$', linestyle='-', color='green')

mask = ((allangles>=0) & (allpositions>=10) & (allpositions<=70) &(allvmags>=00e3) & (allvmags<=600e3))
mask1 = ((allangles>=90) & (allpositions>=10) & (allpositions<=30) &(allvmags>=00e3) & (allvmags<=600e3))
mask2 = ((allangles>=90) & (allpositions>30) & (allpositions<=50) &(allvmags>=00e3) & (allvmags<=600e3))
#mask2 = ((allvmags>=350e3))
mask3 = ((allangles>=90) & (allmachs>=0.6) & (allmachs<=1.6))
mask4 = ((allangles>=90) & (allmachs>=1.6))
x = np.log10(allbeta_par)
y = np.log10(1/allTparperp)
z = alldv_par_norm
#z = 1e-3*allvmags
z = np.log10(allmachs)
# z = allmachs
# z = allangles


# disccolors = ['k','b', 'g', 'y', 'orange','red']
# bounds= [10, 25, 40, 55, 60,70]
disccolors = ['b','r']
bounds= [10, 30, 70]

cmap = mcolors.ListedColormap(disccolors)
discnorm = mcolors.BoundaryNorm(bounds, cmap.N)

norm = colors.LogNorm(vmin=0.5e-1, vmax=0.5e1)
normlin = colors.Normalize(vmin=90,vmax=130)


# z = r_norm**(0)*allT/1.602e-19 # variable for colorbar





# Only 10-30 Rs
# sc = plt.scatter(x[mask],(y[mask]),marker='*',c='g',s=3)
# plt.title(r'Distribution of Magnetic Deflections with $\theta \geq 90^o $',fontsize=15)

# Only 10-30 Rs
# sc = plt.scatter(x[mask1],(y[mask1]),marker='*',c='r',s=3)
# plt.plot(log_betapar, log_Tperppar_cgl2, linewidth=1, 
#         label=r'$\beta^{-0.2}$', linestyle='-', color='red')
# plt.title(r'Distribution of Magnetic Deflections with $\theta \geq 90^o $' + '\n ($R = 10-30 \ R_s$)',fontsize=15)

# Only 30-70 Rs
# sc = plt.scatter(x[mask2],(y[mask2]),marker='*',c='b',s=1)
# plt.plot(log_betapar, log_Tperppar_cgl1, linewidth=1, 
#          label=r'$\beta^{-1}$ (CGL)', linestyle='-', color='blue')
# plt.title(r'Distribution of Magnetic Deflections with $\theta \geq 90^o $' + '\n ($R = 30-50 \ R_s$)',fontsize=15)




# z = allmachs
# sc = plt.scatter(x[mask1],(y[mask1]),marker='*',c=z[mask1],cmap='jet',vmin=0.6,vmax=3,s=10)
# sc = plt.scatter(x[mask2],(y[mask2]),marker='*',c=z[mask2],cmap='jet',vmin=1,vmax=5,s=10)
# cbar = plt.colorbar(sc)
# cbar.set_label(r'$M_a$', fontsize=15)


# z = allpositions
# # sc = plt.scatter(x[mask1],(y[mask1]),marker='*',c=z[mask1],cmap='jet',vmin=10,vmax=30,s=10)
# sc = plt.scatter(x[mask2],(y[mask2]),marker='*',c=z[mask2],cmap='jet',vmin=30,vmax=70,s=10)
# plt.xticks(fontsize=15)
# plt.yticks(fontsize=15)
# cbar = plt.colorbar(sc)
# cbar.set_label(r'$R/R_s$', fontsize=15)
disccolors = ['r','cyan','c','darkcyan']
bounds = [-0.2,0.2,0.3,0.4,0.5]
bounds = [-0.2,0.2,0.4,0.6,0.8,1]
bounds = [-0.2,0.2,0.3,0.4,0.5]
cmap = mcolors.ListedColormap(disccolors)
discnorm = mcolors.BoundaryNorm(bounds,cmap.N)
# Plot the 2D histogram

# sc = plt.scatter(x[mask1],(y[mask1]),marker='*',c='b',vmin=0,vmax=1,s=10)
# sc = plt.scatter(x[mask1],(y[mask1]),marker='*',c=z[mask1],cmap='jet',vmin=0,vmax=1,s=10)
sc = plt.scatter(x[mask1],(y[mask1]),marker='*',c=z[mask1],cmap=cmap,norm=discnorm,s=15)
# sc = plt.scatter(x[mask2],(y[mask2]),marker='*',c=z[mask2],cmap='jet',vmin=0.5,vmax=1,s=10)
# plt.xticks(fontsize=15)
# plt.yticks(fontsize=15)
cbar = plt.colorbar(sc)
cbar.set_label(r'$\log_{10}(M_a)$', fontsize=15)


# sc = plt.scatter(x[mask],(y[mask]),marker='*',c=z[mask],cmap='jet',vmin=10,vmax=70,s=4)
# sc = plt.scatter(x[mask2],(y[mask2]),marker='*',c=z[mask2],cmap='plasma',vmin=30,vmax=70,s=4)
# sc = plt.scatter(x[mask1],(y[mask1]),marker='*',c=z[mask1],cmap=cmap,norm=colors,s=2)
# sc = plt.scatter(x[mask1],(y[mask1]),marker='*',c=allangles[mask1],cmap=cmap,norm=discnorm)

# Add colorbar



# cbar.ax.tick_params(labelsize=15)
plt.ylim(-1,1)
plt.xlim(-1.3,1)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.xlabel(r'$Log_{10}(\beta_\parallel)$', fontsize=15)
plt.ylabel(r'$Log_{10}(T_\perp/T_\parallel)$', fontsize=15)
plt.axhline(y=0,color='k',alpha=0.5)
plt.axvline(x=0,color='k',alpha=0.5)
plt.legend(fontsize=11)
# %





#%%


#Scatterplotting

mirror_params = [0.77, 0.76, -0.016]  # [a, b, beta0]
firehose_params = [-1.4, 1.0, -0.11]  # [a, b, beta0]
cyclotron_params = [0.45, 0.4, -0.0004]  # [a, b, beta0]
parfirehose_params = [-0.47, 0.53, 0.59]  # [a, b, beta0]
# Create beta_parallel array (logarithmically spaced)
betapar = np.logspace(-2, 1, 10000)  # From 0.001 to 100
Tperppar_m = stability_condition(betapar, mirror_params[0], mirror_params[1], mirror_params[2])
Tperppar_f = stability_condition(betapar, firehose_params[0], firehose_params[1], firehose_params[2])
Tperppar_pf = stability_condition(betapar, parfirehose_params[0], parfirehose_params[1], parfirehose_params[2])
Tperppar_c = stability_condition(betapar, cyclotron_params[0], cyclotron_params[1], cyclotron_params[2])
# Calculate T_perp/T_par for both instabilities
log_betapar = np.log10(betapar)
log_Tperppar_m = np.log10(Tperppar_m)
log_Tperppar_f = np.log10(Tperppar_f)
log_Tperppar_pf = np.log10(Tperppar_pf)
log_Tperppar_c = np.log10(Tperppar_c)

# Calculate CGL double adiabatic curves
beta_par_0 = 10**-1.1  # Reference beta_parallel
Tperppar_0 = 10**0.45  # Reference T_perp/T_par (isotropic starting point)
exp1 = 1
exp2 = 0.3
exp3 = 0.6

Tperppar_cgl1 = cgl_curve(betapar, 10**-0.7, Tperppar_0,exp1)
Tperppar_cgl2 = cgl_curve(betapar, 10**0.45, 10**-0.2,0.3)
Tperppar_cgl3 = cgl_curve(betapar, beta_par_0, 10**0.15,exp3)

# Convert to log10 space to match the histogram axes
log_betapar = np.log10(betapar)
log_Tperppar_m = np.log10(Tperppar_m)
log_Tperppar_f = np.log10(Tperppar_f)
log_Tperppar_pf = np.log10(Tperppar_pf)
log_Tperppar_c = np.log10(Tperppar_c)
log_Tperppar_cgl1 = np.log10(Tperppar_cgl1)
log_Tperppar_cgl2 = np.log10(Tperppar_cgl2)
log_Tperppar_cgl3 = np.log10(Tperppar_cgl3)
    
# Plot the stability condition curves
plt.plot(log_betapar, log_Tperppar_m, 'w-', linewidth=3, 
        label='Mirror Instability', linestyle='dotted',color='k')
plt.plot(log_betapar, log_Tperppar_c, 'w-', linewidth=1.5, 
        label='IC Instability', linestyle='dotted',color='k')
plt.plot(log_betapar, log_Tperppar_f, 'w-', linewidth=3, 
        label='Oblique FH', linestyle='--',color='k')
plt.plot(log_betapar, log_Tperppar_pf, 'w-', linewidth=1.5, 
        label='Parallel FH', linestyle='--',color='k')


# Plot CGL double adiabatic curves


# plt.plot(log_betapar, log_Tperppar_cgl3, linewidth=1, 
#         label=r'$\beta^{-0.6}$', linestyle='-', color='green')

mask = ((allangles>=0) & (allpositions>=10) & (allpositions<=70) &(allvmags>=00e3) & (allvmags<=600e3))
mask1 = ((allangles>=90) & (allpositions>=10) & (allpositions<=30) &(allvmags>=00e3) & (allvmags<=600e3))
mask2 = ((allangles>=90) & (allpositions>=30) & (allpositions<=70) &(allvmags>=00e3) & (allvmags<=600e3))
# mask2 = ((allangles>=90) & (allpositions>=50) & (allpositions<=70) &(allvmags>=00e3) & (allvmags<=600e3))
#mask2 = ((allvmags>=350e3))
#mask3 = ((allangles>=70) & (allangles<80))
x = np.log10(allbeta_par)
y = np.log10(1/allTparperp)
z = alldv_par_norm
#z = 1e-3*allvmags
z = abs(allcrosshelicity)
# z = allmachs
# z = allangles


# disccolors = ['k','b', 'g', 'y', 'orange','red']
# bounds= [10, 25, 40, 55, 60,70]
disccolors = ['b','r']
bounds= [10, 30, 70]

cmap = mcolors.ListedColormap(disccolors)
discnorm = mcolors.BoundaryNorm(bounds, cmap.N)

norm = colors.LogNorm(vmin=0.5e-1, vmax=0.5e1)
normlin = colors.Normalize(vmin=90,vmax=130)


# z = r_norm**(0)*allT/1.602e-19 # variable for colorbar





# Only 10-30 Rs
# sc = plt.scatter(x[mask],(y[mask]),marker='*',c='g',s=3)
# plt.title(r'Distribution of Magnetic Deflections with $\theta \geq 90^o $',fontsize=15)

# Only 10-30 Rs
# sc = plt.scatter(x[mask1],(y[mask1]),marker='*',c='r',s=3)
# plt.plot(log_betapar, log_Tperppar_cgl2, linewidth=1, 
#         label=r'$\beta^{-0.3}$', linestyle='-', color='red')
# plt.title(r'Distribution of Magnetic Deflections with $\theta \geq 90^o $' + '\n ($R = 10-30 \ R_s$)',fontsize=15)

# Only 30-70 Rs
sc = plt.scatter(x[mask2],(y[mask2]),marker='*',c='b',s=1)

x_fit = x[mask2]
y_fit = y[mask2]
x_fit = x_fit[np.isfinite(x_fit) & np.isfinite(y_fit)]
y_fit = y_fit[np.isfinite(x_fit) & np.isfinite(y_fit)]
m,c = np.polyfit(x_fit,y_fit,1)
xline = np.linspace(-2,1,10000)
yline = m*xline + c

plt.plot(xline, yline, 'r--', linewidth=1)

plt.plot(log_betapar, log_Tperppar_cgl1, linewidth=1, 
         label=r'$\beta^{-1}$ (CGL)', linestyle='-', color='blue')
plt.title(r'Distribution of Magnetic Deflections with $\theta \geq 90^o $' + '\n ($R = 30-70 \ R_s$)',fontsize=15)




# plt.plot(log_betapar, log_Tperppar_cgl1, linewidth=1, 
#          label=r'$\beta^{-1}$ (CGL)', linestyle='-', color='blue')
# plt.title(r'Distribution of Magnetic Deflections with $\theta \geq 90^o $' + '\n ($R = 30-70 \ R_s$)',fontsize=15)




# z = allmachs
# sc = plt.scatter(x[mask1],(y[mask1]),marker='*',c=z[mask1],cmap='jet',vmin=0.6,vmax=3,s=10)
# sc = plt.scatter(x[mask2],(y[mask2]),marker='*',c=z[mask2],cmap='jet',vmin=1,vmax=5,s=10)
# cbar = plt.colorbar(sc)
# cbar.set_label(r'$M_a$', fontsize=15)


# z = allpositions
# # sc = plt.scatter(x[mask1],(y[mask1]),marker='*',c=z[mask1],cmap='jet',vmin=10,vmax=30,s=10)
# sc = plt.scatter(x[mask2],(y[mask2]),marker='*',c=z[mask2],cmap='jet',vmin=30,vmax=70,s=10)
# plt.xticks(fontsize=15)
# plt.yticks(fontsize=15)
# cbar = plt.colorbar(sc)
# cbar.set_label(r'$R/R_s$', fontsize=15)


# sc = plt.scatter(x[mask1],(y[mask1]),marker='*',c=z[mask1],cmap='jet',vmin=0,vmax=1,s=10)
# sc = plt.scatter(x[mask2],(y[mask2]),marker='*',c=z[mask2],cmap='jet',vmin=0.5,vmax=1,s=10)
# plt.xticks(fontsize=15)
# plt.yticks(fontsize=15)
# cbar = plt.colorbar(sc)
# cbar.set_label(r'$|\sigma_c|$', fontsize=15)


# sc = plt.scatter(x[mask],(y[mask]),marker='*',c=z[mask],cmap='jet',vmin=10,vmax=70,s=4)
# sc = plt.scatter(x[mask2],(y[mask2]),marker='*',c=z[mask2],cmap='plasma',vmin=30,vmax=70,s=4)
# sc = plt.scatter(x[mask1],(y[mask1]),marker='*',c=z[mask1],cmap=cmap,norm=colors,s=2)
# sc = plt.scatter(x[mask1],(y[mask1]),marker='*',c=allangles[mask1],cmap=cmap,norm=discnorm)

# Add colorbar



# cbar.ax.tick_params(labelsize=15)
plt.ylim(-1,1)
plt.xlim(-2,1)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.xlabel(r'$Log_{10}(\beta_\parallel)$', fontsize=15)
plt.ylabel(r'$Log_{10}(T_\perp/T_\parallel)$', fontsize=15)
plt.axhline(y=0,color='k',alpha=0.5)
plt.axvline(x=0,color='k',alpha=0.5)
plt.legend(fontsize=11)

# %%



#Scatterplotting (trying to fit stuff)

mirror_params = [0.77, 0.76, -0.016]  # [a, b, beta0]
firehose_params = [-1.4, 1.0, -0.11]  # [a, b, beta0]
cyclotron_params = [0.45, 0.4, -0.0004]  # [a, b, beta0]
parfirehose_params = [-0.47, 0.53, 0.59]  # [a, b, beta0]
# Create beta_parallel array (logarithmically spaced)
betapar = np.logspace(-2, 1, 10000)  # From 0.001 to 100
Tperppar_m = stability_condition(betapar, mirror_params[0], mirror_params[1], mirror_params[2])
Tperppar_f = stability_condition(betapar, firehose_params[0], firehose_params[1], firehose_params[2])
Tperppar_pf = stability_condition(betapar, parfirehose_params[0], parfirehose_params[1], parfirehose_params[2])
Tperppar_c = stability_condition(betapar, cyclotron_params[0], cyclotron_params[1], cyclotron_params[2])
# Calculate T_perp/T_par for both instabilities
log_betapar = np.log10(betapar)
log_Tperppar_m = np.log10(Tperppar_m)
log_Tperppar_f = np.log10(Tperppar_f)
log_Tperppar_pf = np.log10(Tperppar_pf)
log_Tperppar_c = np.log10(Tperppar_c)

# Calculate CGL double adiabatic curves
beta_par_0 = 10**-1.1  # Reference beta_parallel
Tperppar_0 = 10**0.45  # Reference T_perp/T_par (isotropic starting point)
exp1 = 1
exp2 = 0.3
exp3 = 0.6

Tperppar_cgl1 = cgl_curve(betapar, 10**-0.7, Tperppar_0,exp1)
Tperppar_cgl2 = cgl_curve(betapar, 10**0.45, 10**-0.2,0.3)
Tperppar_cgl3 = cgl_curve(betapar, beta_par_0, 10**0.15,exp3)

# Convert to log10 space to match the histogram axes
log_betapar = np.log10(betapar)
log_Tperppar_m = np.log10(Tperppar_m)
log_Tperppar_f = np.log10(Tperppar_f)
log_Tperppar_pf = np.log10(Tperppar_pf)
log_Tperppar_c = np.log10(Tperppar_c)
log_Tperppar_cgl1 = np.log10(Tperppar_cgl1)
log_Tperppar_cgl2 = np.log10(Tperppar_cgl2)
log_Tperppar_cgl3 = np.log10(Tperppar_cgl3)
    
# Plot the stability condition curves
plt.plot(log_betapar, log_Tperppar_m, 'w-', linewidth=3, 
        label='Mirror Instability', linestyle='dotted',color='k')
plt.plot(log_betapar, log_Tperppar_c, 'w-', linewidth=1.5, 
        label='IC Instability', linestyle='dotted',color='k')
plt.plot(log_betapar, log_Tperppar_f, 'w-', linewidth=3, 
        label='Oblique FH', linestyle='--',color='k')
plt.plot(log_betapar, log_Tperppar_pf, 'w-', linewidth=1.5, 
        label='Parallel FH', linestyle='--',color='k')


# Plot CGL double adiabatic curves


# plt.plot(log_betapar, log_Tperppar_cgl3, linewidth=1, 
#         label=r'$\beta^{-0.6}$', linestyle='-', color='green')

mask = ((allangles>=90) & (allpositions>=10) & (allpositions<=70) &(allvmags>=00e3) & (allvmags<=600e3))
mask1 = ((allangles>=90) & (allpositions>=10) & (allpositions<=30) &(allvmags>=00e3) & (allvmags<=600e3))
mask2 = ((allangles>=90) & (allpositions>=30) & (allpositions<=60) &(allvmags>=00e3) & (allvmags<=600e3))
# mask2 = ((allangles>=90) & (allpositions>=50) & (allpositions<=70) &(allvmags>=00e3) & (allvmags<=600e3))
#mask2 = ((allvmags>=350e3))
#mask3 = ((allangles>=70) & (allangles<80))
x = np.log10(allbeta_par)
y = np.log10(1/allTparperp)
z = alldv_par_norm
#z = 1e-3*allvmags
z = abs(allcrosshelicity)
# z = allmachs
# z = allangles


# disccolors = ['k','b', 'g', 'y', 'orange','red']
# bounds= [10, 25, 40, 55, 60,70]
disccolors = ['b','r']
bounds= [10, 30, 70]

cmap = mcolors.ListedColormap(disccolors)
discnorm = mcolors.BoundaryNorm(bounds, cmap.N)

norm = colors.LogNorm(vmin=0.5e-1, vmax=0.5e1)
normlin = colors.Normalize(vmin=90,vmax=130)


# z = r_norm**(0)*allT/1.602e-19 # variable for colorbar





# Only 10-30 Rs
# sc = plt.scatter(x[mask],(y[mask]),marker='*',c='g',s=3)
# plt.title(r'Distribution of Magnetic Deflections with $\theta \geq 90^o $',fontsize=15)

# Only 10-30 Rs
# sc = plt.scatter(x[mask1],(y[mask1]),marker='*',c='r',s=3)
# plt.plot(log_betapar, log_Tperppar_cgl2, linewidth=1, 
#         label=r'$\beta^{-0.3}$', linestyle='-', color='red')
# plt.title(r'Distribution of Magnetic Deflections with $\theta \geq 90^o $' + '\n ($R = 10-30 \ R_s$)',fontsize=15)

sc = plt.scatter(x[mask1], (y[mask1]), marker='*', c='b', s=1)
x_fit = x[mask1]
y_fit = y[mask1]
# Clean data
valid = np.isfinite(x_fit) & np.isfinite(y_fit)
x_fit = x_fit[valid]
y_fit = y_fit[valid]
m, c = np.polyfit(x_fit, y_fit, 1)
# Create fit line
x_line = np.linspace(min(x_fit), max(x_fit), 200)
y_line = m * x_line + c
# Plot fit
plt.plot(x_line, y_line, 'r-', linewidth=2,
         label=f'Best fit (slope = {m:.2f})')
# --- END BLOCK ---
# Title (unchanged)
plt.title(r'Distribution of Magnetic Deflections with $\theta \geq 90^o $' 
          + '\n ($R = 10-30 \ R_s$)', fontsize=15)
plt.legend()





# sc = plt.scatter(x[mask], (y[mask]), marker='*', c='b', s=1)

# # --- ADD THIS BLOCK (fit only mask2 data) ---
# x_fit = x[mask]
# y_fit = y[mask]

# # Clean data
# valid = np.isfinite(x_fit) & np.isfinite(y_fit)
# x_fit = x_fit[valid]
# y_fit = y_fit[valid]

# # Linear fit in log-log space
# m, c = np.polyfit(x_fit, y_fit, 1)

# # Create fit line
# x_line = np.linspace(min(x_fit), max(x_fit), 200)
# y_line = m * x_line + c

# # Plot fit
# plt.plot(x_line, y_line, 'r-', linewidth=2,
#          label=f'Best fit (slope = {m:.2f})')
# # --- END BLOCK ---
# # Title (unchanged)
# plt.title(r'Distribution of Magnetic Deflections with $\theta \geq 90^o $' 
#           + '\n ($R = 30-70 \ R_s$)', fontsize=15)

# plt.legend()






# plt.plot(log_betapar, log_Tperppar_cgl1, linewidth=1, 
#          label=r'$\beta^{-1}$ (CGL)', linestyle='-', color='blue')
# plt.title(r'Distribution of Magnetic Deflections with $\theta \geq 90^o $' + '\n ($R = 30-70 \ R_s$)',fontsize=15)




# z = allmachs
# sc = plt.scatter(x[mask1],(y[mask1]),marker='*',c=z[mask1],cmap='jet',vmin=0.6,vmax=3,s=10)
# sc = plt.scatter(x[mask2],(y[mask2]),marker='*',c=z[mask2],cmap='jet',vmin=1,vmax=5,s=10)
# cbar = plt.colorbar(sc)
# cbar.set_label(r'$M_a$', fontsize=15)


# z = allpositions
# # sc = plt.scatter(x[mask1],(y[mask1]),marker='*',c=z[mask1],cmap='jet',vmin=10,vmax=30,s=10)
# sc = plt.scatter(x[mask2],(y[mask2]),marker='*',c=z[mask2],cmap='jet',vmin=30,vmax=70,s=10)
# plt.xticks(fontsize=15)
# plt.yticks(fontsize=15)
# cbar = plt.colorbar(sc)
# cbar.set_label(r'$R/R_s$', fontsize=15)


# sc = plt.scatter(x[mask1],(y[mask1]),marker='*',c=z[mask1],cmap='jet',vmin=0,vmax=1,s=10)
# sc = plt.scatter(x[mask2],(y[mask2]),marker='*',c=z[mask2],cmap='jet',vmin=0.5,vmax=1,s=10)
# plt.xticks(fontsize=15)
# plt.yticks(fontsize=15)
# cbar = plt.colorbar(sc)
# cbar.set_label(r'$|\sigma_c|$', fontsize=15)


# sc = plt.scatter(x[mask],(y[mask]),marker='*',c=z[mask],cmap='jet',vmin=10,vmax=70,s=4)
# sc = plt.scatter(x[mask2],(y[mask2]),marker='*',c=z[mask2],cmap='plasma',vmin=30,vmax=70,s=4)
# sc = plt.scatter(x[mask1],(y[mask1]),marker='*',c=z[mask1],cmap=cmap,norm=colors,s=2)
# sc = plt.scatter(x[mask1],(y[mask1]),marker='*',c=allangles[mask1],cmap=cmap,norm=discnorm)

# Add colorbar



# cbar.ax.tick_params(labelsize=15)
plt.ylim(-1,1)
plt.xlim(-2,1)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.xlabel(r'$Log_{10}(\beta_\parallel)$', fontsize=15)
plt.ylabel(r'$Log_{10}(T_\perp/T_\parallel)$', fontsize=15)
plt.axhline(y=0,color='k',alpha=0.5)
plt.axvline(x=0,color='k',alpha=0.5)
plt.legend(fontsize=11)