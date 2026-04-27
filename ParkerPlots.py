#%%

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import binned_statistic_2d
from matplotlib.colors import ListedColormap
import matplotlib.cm as cm
import matplotlib.colors as colors

def overview_plot():
    # Deflection Angle vs mach number and position
    bins = 30
    fig,ax = plt.subplots(4,1,figsize=(8,16))

    count=1
    # allmachs = filtered_mach

    ax[0].scatter(np.log10(allmachs),allangles,s=0.1)
    ax[0].set_xlim(-0.8,0.8)
    ax[0].set_ylim(0,180)
    ax[0].axhline(y=90,color='k')
    ax[0].axvline(x=0,color='k')
    ax[0].set_ylabel("Deflection Angle (degrees)")
    ax[0].set_xlabel("Log10(Ma)")

    ax[1].hist2d(np.log10(allmachs),allangles,range=[[-1,1.0],[0,180]],bins=bins,cmin=count,cmap='viridis')
    ax[1].set_xlabel('Log10(Ma)')
    ax[1].set_ylabel('Deflection Angle (degrees)')
    ax[1].set_facecolor('k')
    ax[1].axhline(y=90,color='w')
    ax[1].axvline(x=0,color='w')
    ax[1].set_xlim(-0.8,0.8)
    ax[1].set_ylim(0,180)

    ax[2].scatter(allpositions,allangles,s=0.1)
    ax[2].set_xlim(0,45)
    ax[2].set_ylim(0,180)
    ax[2].axhline(y=90,color='k')
    ax[2].axvline(x=0,color='k')
    ax[2].set_ylabel("Deflection Angle (degrees)")
    ax[2].set_xlabel("R (Soalr Radii)")

    ax[3].hist2d(allpositions,allangles,range=[[10,45],[0,180]],bins=bins,cmin=count,cmap='viridis')
    ax[3].set_xlabel('R (Solar Radii)')
    ax[3].set_ylabel('Deflection Angle (degrees)')
    ax[3].set_facecolor('k')
    ax[3].axhline(y=90,color='w')
    ax[3].axvline(x=0,color='w')
    ax[3].set_xlim(0,45)
    ax[3].set_ylim(0,180)
    return
overview_plot()

#%%
## Mach number vs Deflection Angle Distribution
def angledist_plot():
    x = np.log10(filtered_mach) # x axis
    y = allangles # y axis
    z = alldv_perp_norm  # variable for colorbar

    # Define bins
    bins=30
    x_bins = np.linspace(-0.5, 0.5, bins)
    y_bins = np.linspace(0, 180, bins)

    # Calculate the mean of 'z' for each bin
    mean_z, x_edge, y_edge, bin_number = binned_statistic_2d(x, y, z, statistic='count', bins=[x_bins, y_bins])

    # Calculate counts to identify zero bins
    count_z, _, _, _ = binned_statistic_2d(x, y, z, statistic='count', bins=[x_bins, y_bins])
    # Create a masked array where zero-count bins are masked
    mean_z_masked = np.ma.masked_where(count_z == 0, mean_z)

    # Plot the result
    fig, ax = plt.subplots(figsize=(10,7))

    # First, plot the zero-count bins in light grey
    zero_mask = count_z == 0
    if np.any(zero_mask):
        # Create an array filled with a dummy value for zero bins
        zero_array = np.full_like(mean_z, 1.0)  # Use middle value of your range
        zero_array[~zero_mask] = np.nan  # Set non-zero bins to NaN so they're transparent
        
        ax.imshow(zero_array.T, origin='lower', 
                extent=[x_edge[0], x_edge[-1], y_edge[0], y_edge[-1]], 
                cmap='gray', aspect='auto', vmin=0.8, vmax=1.2, alpha=0.4)

    # Then plot the actual data on top
    norm = colors.LogNorm(vmin=1,vmax=1e4)
    im = ax.imshow(mean_z_masked.T, origin='lower', 
                extent=[x_edge[0], x_edge[-1], y_edge[0], y_edge[-1]], 
                cmap='jet', aspect='auto', norm=norm)


    cbar = plt.colorbar(im, label='counts')
    cbar.ax.tick_params(labelsize=15)
    cbar.set_label('counts',fontsize=15)
    plt.xlim(-0.4, 0.5)
    plt.ylim(10, 180)
    plt.xlabel(r'$\log_{10}(M_a)$',fontsize=15)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.ylabel(r'$\theta$ [degrees]',fontsize=15)
    plt.axhline(y=90, color='k')
    plt.axvline(x=0, color='k')
    plt.show()
    return
angledist_plot()


#%%
#Changing up axis

# def angledist_plot():
#     x = filtered_mach # x axis
#     y = allangles # y axis
#     z = alldv_perp_norm  # variable for colorbar

#     # Define bins
#     bins=40
#     x_bins = np.logspace(-0.4,0.5, bins)
#     y_bins = np.linspace(0, 180, bins)

#     # Calculate the mean of 'z' for each bin
#     mean_z, x_edge, y_edge, bin_number = binned_statistic_2d(x, y, z, statistic='count', bins=[x_bins, y_bins])
#     # mean_z, x_edge, y_edge, bin_number = binned_statistic_2d(x, y, z, statistic=lambda v: np.nanmean(v), bins=[x_bins, y_bins])
    
    
#     # Calculate counts to identify zero bins
#     count_z, _, _, _ = binned_statistic_2d(x, y, z, statistic='count', bins=[x_bins, y_bins])
#     # Create a masked array where zero-count bins are masked
#     mean_z_masked = np.ma.masked_where(count_z == 0, mean_z)

#     # Plot the result
#     fig, ax = plt.subplots(figsize=(10,7))


#     norm = colors.LogNorm(vmin=1,vmax=1e4)

#     X, Y = np.meshgrid(x_edge, y_edge)
#     cmap = plt.cm.plasma.copy()
#     cmap.set_bad(color='lightgrey')  # masked cells → grey instead of white
#     im = ax.pcolormesh(X, Y, mean_z_masked.T, cmap=cmap,norm=norm)  # switch to pcolormesh
#     ax.set_xscale('log')  # now log scale works


#     cbar = plt.colorbar(im, label='counts')
#     cbar.ax.tick_params(labelsize=12)
#     cbar.set_label('counts',fontsize=12)
#     plt.xlim(10**-0.4, 10**0.5)
#     plt.ylim(10, 180)
#     plt.xlabel(r'$M_a$',fontsize=12)
#     plt.xticks(fontsize=12)
#     plt.yticks(fontsize=12)
#     plt.ylabel(r'$\theta_{def}$ (degrees)',fontsize=12)
#     plt.axhline(y=90, color='k')
#     plt.axvline(x=1, color='k')
#     plt.show()
#     return
# angledist_plot()











#%%
## Position vs Deflection Angle Distribution
def angledistpos_plot():
    x = allpositions # x axis
    y = allangles # y axis
    z = alldv_norm_mag  # variable for colorbar

    # Define bins

    bins=30
    x_bins = np.linspace(9, 40, bins)
    y_bins = np.linspace(0, 180, bins)

    # Calculate the mean of 'z' for each bin
    mean_z, x_edge, y_edge, bin_number = binned_statistic_2d(x, y, z, statistic='count', bins=[x_bins, y_bins])

    # Calculate counts to identify zero bins
    count_z, _, _, _ = binned_statistic_2d(x, y, z, statistic='count', bins=[x_bins, y_bins])
    # Create a masked array where zero-count bins are masked
    mean_z_masked = np.ma.masked_where(count_z == 0, mean_z)

    # Plot the result
    fig, ax = plt.subplots(figsize=(10,7))

    # First, plot the zero-count bins in light grey
    zero_mask = count_z == 0
    if np.any(zero_mask):
        # Create an array filled with a dummy value for zero bins
        zero_array = np.full_like(mean_z, 1.0)  # Use middle value of your range
        zero_array[~zero_mask] = np.nan  # Set non-zero bins to NaN so they're transparent
        
        ax.imshow(zero_array.T, origin='lower', 
                extent=[x_edge[0], x_edge[-1], y_edge[0], y_edge[-1]], 
                cmap='gray', aspect='auto', vmin=0.8, vmax=1.2, alpha=0.4)

    # Then plot the actual data on top
    norm = colors.LogNorm(vmin=1e0,vmax=1e4)
    im = ax.imshow(mean_z_masked.T, origin='lower', 
                extent=[x_edge[0], x_edge[-1], y_edge[0], y_edge[-1]], 
                cmap='jet', aspect='auto', norm=norm)

    cbar = plt.colorbar(im, label='counts')
    cbar.ax.tick_params(labelsize=15)
    cbar.set_label('counts',fontsize=15)
    plt.xlim(9, 40)
    plt.ylim(10, 180)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.xlabel(r'$R/R_{sun}$',fontsize=15)
    plt.ylabel(r'$\theta$ [degrees]',fontsize=15)
    plt.axhline(y=90, color='k')
    plt.axvline(x=0, color='k')
    plt.show()
    return

angledistpos_plot()
#%%
## Mach number vs Deflection Angle: dv/v Values 
def param_plot(z,label,limits = [0,2],cmap='seismic'):
    x = np.log10(filtered_mach) # x axis
    y = allangles # y axis
    #z = alldv_norm_mag # variable for colorbar
    #z = alldv_perp_norm # variable for colorbar
    #z = alldv_norm_mag*allvmags/allva # variable for colorbar
    z = 1000*allentropy/(1.602e-19) #puts it back in ev/m^2 # variable for colorbar

    # Define bins
    bins=30
    x_bins = np.linspace(-0.5, 0.5, bins)
    y_bins = np.linspace(0, 180, bins)

    # Calculate the mean of 'z' for each bin
    mean_z, x_edge, y_edge, bin_number = binned_statistic_2d(x, y, z, statistic='median', bins=[x_bins, y_bins])
    # mean_z, x_edge, y_edge, bin_number = binned_statistic_2d(x, y, z, statistic='std', bins=[x_bins, y_bins])
    mean_z, x_edge, y_edge, bin_number = binned_statistic_2d(x, y, z, statistic=lambda v: np.nanmedian(v), bins=[x_bins, y_bins])
    # mean_z, x_edge, y_edge, bin_number = binned_statistic_2d(x, y, z, statistic=lambda v: np.nanstd(v)/np.sqrt(np.sum(~np.isnan(v))), bins=[x_bins, y_bins])

    # Calculate counts to identify zero bins
    count_z, _, _, _ = binned_statistic_2d(x, y, z, statistic='count', bins=[x_bins, y_bins])
    # Create a masked array where zero-count bins are masked
    mean_z_masked = np.ma.masked_where(count_z == 0, mean_z)

    # Plot the result
    fig, ax = plt.subplots()

    # First, plot the zero-count bins in light grey
    zero_mask = count_z == 0
    if np.any(zero_mask):
        # Create an array filled with a dummy value for zero bins
        zero_array = np.full_like(mean_z, 1.0)  # Use middle value of your range
        zero_array[~zero_mask] = np.nan  # Set non-zero bins to NaN so they're transparent
        
        ax.imshow(zero_array.T, origin='lower', 
                extent=[x_edge[0], x_edge[-1], y_edge[0], y_edge[-1]], 
                cmap='gray', aspect='auto', vmin=0.8, vmax=1.2, alpha=0.4)

    # Then plot the actual data on top

    norm = colors.LogNorm(vmin=1,vmax=1e4)

    im = ax.imshow(mean_z_masked.T, origin='lower', 
                extent=[x_edge[0], x_edge[-1], y_edge[0], y_edge[-1]], 
                #cmap='nipy_spectral', aspect='auto', vmin=0., vmax=2)
                cmap=cmap, aspect='auto', vmin=limits[0],vmax=limits[1])
    
   
    cbar = plt.colorbar(im, label='|dv/va|')
    cbar.ax.tick_params(labelsize=15)
    # cbar.set_label(r'$(\delta v)_\parallel/ \ v_a \ Bin \ Counts$',fontsize=12)
    cbar.set_label(label,fontsize=15)
    # cbar.set_label(r'$\delta v/v$',fontsize=12)
    plt.xlim(-0.5, 0.5)
    plt.ylim(10, 180)

    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.xlabel(r'$\log_{10}(M_a)$',fontsize=15)
    plt.ylabel(r'$\theta $  [degrees]',fontsize=15)
    plt.axhline(y=90, color='k')
    plt.axvline(x=0, color='k')
    plt.show()
    return

# param_plot(abs(allpositions),label = r'$R/R_s$',limits=[10,40],cmap='jet')






param_plot(allcrosshelicity,label = r'$V_p \  (km/s)$',limits=[0,0.1],cmap='jet')

#%%
# Median Values
param_plot(abs(alldv_norm_mag),label = r'$(\delta v)/ \ v $')
param_plot(abs(alldv_norm_mag*allvmags/allva),label = r'$(\delta v)/ \ v_a$')
param_plot(abs(alldv_par_norm),label = r'$(\delta v)_\parallel/ \ v$',limits=[0,1],cmap='jet')
param_plot((r_norm**(0))*abs(allT)/1.602e-19,label = r'$ T \ (eV)$',limits=[0,200],cmap='jet')
#%%
# Standard Deviations
param_plot(abs(alldv_norm_mag),label = r'SEM of $(\delta v)/ \ v $',cmap = 'jet',limits = [0,0.3])
param_plot(abs(alldv_norm_mag*allvmags/allva),label = r'SEM of $(\delta v)/ \ v_a$',cmap = 'jet',limits=[0,0.3])
param_plot(abs(alldv_par_norm),label = r'SEM of $(\delta v)_\parallel/ \ v$',cmap='jet',limits=[0,0.3])
param_plot(abs(alldv_perp_norm),label = r'SEM of $(\delta v)_\perp/ \ v $',cmap='jet',limits=[0,0.3])


# %%
#Energy flix ratio
def param_logplot(z,label,limits=[0,2],cmap='seismic'):
    x = np.log10(filtered_mach) # x axis
    y = allangles # y axis
    #z = alldv_norm_mag # variable for colorbar
    #z = alldv_perp_norm # variable for colorbar
    #z = alldv_norm_mag*allvmags/allva # variable for colorbar

    # Define bins
    bins=30
    x_bins = np.linspace(-0.5, 0.5, bins)
    y_bins = np.linspace(0, 180, bins)

    # Calculate the mean of 'z' for each bin
    # mean_z, x_edge, y_edge, bin_number = binned_statistic_2d(x, y, z, statistic='med', bins=[x_bins, y_bins])
    # mean_z, x_edge, y_edge, bin_number = binned_statistic_2d(x, y, z, statistic='std', bins=[x_bins, y_bins])
    mean_z, x_edge, y_edge, bin_number = binned_statistic_2d(x, y, z, statistic=lambda v: np.nanmedian(v), bins=[x_bins, y_bins])
    # mean_z, x_edge, y_edge, bin_number = binned_statistic_2d(x, y, z, statistic=lambda v: np.nanstd(v)/np.sqrt(np.sum(~np.isnan(v))), bins=[x_bins, y_bins])

    # Calculate counts to identify zero bins
    count_z, _, _, _ = binned_statistic_2d(x, y, z, statistic='count', bins=[x_bins, y_bins])
    # Create a masked array where zero-count bins are masked
    mean_z_masked = np.ma.masked_where(count_z == 0, mean_z)

    # Plot the result
    fig, ax = plt.subplots()

    # First, plot the zero-count bins in light grey
    zero_mask = count_z == 0
    
    if np.any(zero_mask):
        # Create an array filled with a dummy value for zero bins
        zero_array = np.full_like(mean_z, 1.0)  # Use middle value of your range
        zero_array[~zero_mask] = np.nan  # Set non-zero bins to NaN so they're transparent
        
        ax.imshow(zero_array.T, origin='lower', 
                extent=[x_edge[0], x_edge[-1], y_edge[0], y_edge[-1]], 
                cmap='gray', aspect='auto', vmin=0.8, vmax=1.2, alpha=0.4)

    # Then plot the actual data on top

    norm = colors.LogNorm(vmin=limits[0],vmax=limits[1])
    cmap = plt.cm.seismic.copy()
    cmap.set_bad(color='lightgrey')  # masked cells → grey instead of white
    im = ax.imshow(mean_z_masked.T, origin='lower', 
                extent=[x_edge[0], x_edge[-1], y_edge[0], y_edge[-1]], 
                #cmap='nipy_spectral', aspect='auto', vmin=0., vmax=2)
                cmap=cmap, aspect='auto', norm=norm)
    

    cbar = plt.colorbar(im, label='|dv/va|')
    cbar.ax.tick_params(labelsize=15)
    # cbar.set_label(r'$(\delta v)_\parallel/ \ v_a \ Bin \ Counts$',fontsize=12)
    cbar.set_label(label,fontsize=15)
    # cbar.set_label(r'$\delta v/v$',fontsize=12)
    plt.xlim(-0.5, 0.5)
    plt.ylim(10, 180)
    plt.axvline(np.log10(np.sqrt(2)), color='k', linestyle='-.',linewidth=3)

    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.xlabel(r'$\log_{10}(M_a)$',fontsize=15)
    plt.ylabel(r'$\theta$ (degrees)',fontsize=15)
    plt.axhline(y=90, color='k')
    plt.axvline(x=0, color='k')
    plt.show()
    return

param_logplot(abs(allKr/allSr),label = r'$|K_R/S_R|$',limits = [0.1,10])
# param_logplot(abs(allKr/allSr),label = r'SEM of $|K_R/S_R|$',limits = [0.01,100])










#%%

def fluxratio_plot():
    x = np.log10(filtered_mach) # x axis
    y = allangles # y axis
    z = allKr/allSr # variable for colorbar
    #z = alldKmag/alldSmag # variable for colorbar
    #z1 = (allvmags*allBmags**2)/mu0 # variable for colorbar

    # z2 = 0.5*alln*mi*allvmags**3 # variable for colorbar
    #z = z
    # Define bins
    bins=40
    x_bins = np.linspace(-0.6, 0.6, bins)
    y_bins = np.linspace(0, 180, bins)

    # Calculate the mean of 'z' for each bin
    mean_z, x_edge, y_edge, bin_number = binned_statistic_2d(x, y, z, statistic='median', bins=[x_bins, y_bins])
    mean_z, x_edge, y_edge, bin_number = binned_statistic_2d(x, y, z, statistic=lambda v: np.nanmedian(v), bins=[x_bins, y_bins])

    # Calculate counts to identify zero bins
    count_z, _, _, _ = binned_statistic_2d(x, y, z, statistic='count', bins=[x_bins, y_bins])
    # Create a masked array where zero-count bins are masked
    mean_z_masked = np.ma.masked_where(count_z == 0, mean_z)

    # Plot the result
    fig, ax = plt.subplots()

    # First, plot the zero-count bins in light grey
    zero_mask = count_z == 0
    if np.any(zero_mask):
        # Create an array filled with a dummy value for zero bins
        zero_array = np.full_like(mean_z, 1.0)  # Use middle value of your range
        zero_array[~zero_mask] = np.nan  # Set non-zero bins to NaN so they're transparent
        
        ax.imshow(zero_array.T, origin='lower', 
                extent=[x_edge[0], x_edge[-1], y_edge[0], y_edge[-1]], 
                cmap='gray', aspect='auto', vmin=0.8, vmax=1.2, alpha=0.4)


    # Create logarithmic normalization with range 0.1 to 10, centered at 1
    norm = colors.LogNorm(vmin=0.1, vmax=10)
    im = ax.imshow(mean_z_masked.T, origin='lower', 
                extent=[x_edge[0], x_edge[-1], y_edge[0], y_edge[-1]], 
                cmap='seismic', aspect='auto', norm=norm)




    cbar = plt.colorbar(im, label='|dv/v|')
    cbar.ax.tick_params(labelsize=12)
    cbar.set_label(r'$|K_R/S_R|$',fontsize=12)
    plt.xlim(-0.6, 0.5)
    plt.ylim(10, 180)

    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.xlabel('$\log_{10}(M_a)$',fontsize=15)
    plt.ylabel(r'$\theta$ (degrees)',fontsize=15)
    plt.axhline(y=90, color='k')
    plt.axvline(x=0, color='k')
    plt.show()
    return
fluxratio_plot()
#%%

#x = np.log10(filtered_mach) # x axis
#y = allangles # y axis
# z = (allKr/allSr)#*(1/norm_factor) # variable for colorbar


# # Define bins
# bins=40
# x_bins = np.linspace(-0.6, 0.4, bins)
# y_bins = np.linspace(0, 180, bins)

# # Calculate the mean of 'z' for each bin
# mean_z, x_edge, y_edge, bin_number = binned_statistic_2d(x, y, z, statistic='mean', bins=[x_bins, y_bins])

# # Calculate counts to identify zero bins
# count_z, _, _, _ = binned_statistic_2d(x, y, z, statistic='count', bins=[x_bins, y_bins])
# # Create a masked array where zero-count bins are masked
# mean_z_masked = np.ma.masked_where(count_z == 0, mean_z)

# # Plot the result
# fig, ax = plt.subplots()

# # First, plot the zero-count bins in light grey
# zero_mask = count_z == 0
# if np.any(zero_mask):
#     # Create an array filled with a dummy value for zero bins
#     zero_array = np.full_like(mean_z, 1.0)  # Use middle value of your range
#     zero_array[~zero_mask] = np.nan  # Set non-zero bins to NaN so they're transparent
    
#     ax.imshow(zero_array.T, origin='lower', 
#               extent=[x_edge[0], x_edge[-1], y_edge[0], y_edge[-1]], 
#               cmap='gray', aspect='auto', vmin=0.8, vmax=1.2, alpha=0.4)

# # Create logarithmic normalization with range 0.1 to 10, centered at 1
# norm = colors.LogNorm(vmin=0.1, vmax=10)

# im = ax.imshow(mean_z_masked.T, origin='lower', 
#                extent=[x_edge[0], x_edge[-1], y_edge[0], y_edge[-1]], 
#                cmap='seismic', aspect='auto', norm=norm)

# cbar = plt.colorbar(im, label='|Kr/Sr|')
# cbar.ax.tick_params(labelsize=12)
# cbar.set_label('|Kr/Sr| * (r/<r>)',fontsize=12)
# #cbar.set_label('|Kr/Sr| ',fontsize=12)

# # Set custom colorbar ticks to highlight the center at 1
# cbar.set_ticks([0.1, 1, 10])
# cbar.set_ticklabels(['0.1',  '1',  '10'])

# plt.xlim(-0.6, 0.4)
# plt.ylim(10, 180)

# plt.xticks(fontsize=12)
# plt.yticks(fontsize=12)
# plt.xlabel('Log10(Ma)',fontsize=12)
# plt.ylabel('Deflection Angle (degrees)',fontsize=12)
# plt.axhline(y=90, color='k')
# plt.axvline(x=0, color='k')
# plt.show()





#%%

#version with logscale colorbar

# x = np.log10(allmachs) # x axis
# y = allangles # y axis
# z = (allKr/allSr)#*(norm_factor) # variable for colorbar
# # z = (allKr)#*(norm_factor) # variable for colorbar
# #z = alldSmag

# # Define bins
# bins=40
# x_bins = np.linspace(-0.5, 0.5, bins)
# y_bins = np.linspace(0, 180, bins)

# # Calculate the mean of 'z' for each bin
# mean_z, x_edge, y_edge, bin_number = binned_statistic_2d(x, y, z, statistic='mean', bins=[x_bins, y_bins])

# # Calculate counts to identify zero bins
# count_z, _, _, _ = binned_statistic_2d(x, y, z, statistic='count', bins=[x_bins, y_bins])
# # Create a masked array where zero-count bins are masked
# mean_z_masked = np.ma.masked_where(count_z == 0, mean_z)

# # Plot the result
# fig, ax = plt.subplots()

# # First, plot the zero-count bins in light grey
# zero_mask = count_z == 0
# if np.any(zero_mask):
#     # Create an array filled with a dummy value for zero bins
#     zero_array = np.full_like(mean_z, 1.0)  # Use middle value of your range
#     zero_array[~zero_mask] = np.nan  # Set non-zero bins to NaN so they're transparent
    
#     ax.imshow(zero_array.T, origin='lower', 
#               extent=[x_edge[0], x_edge[-1], y_edge[0], y_edge[-1]], 
#               cmap='gray', aspect='auto', vmin=0.8, vmax=1.2, alpha=0.4)

# # Create logarithmic normalization with range 0.1 to 10, centered at 1
# norm = colors.LogNorm(vmin=0.1, vmax=10)

# im = ax.imshow(mean_z_masked.T, origin='lower', 
#                extent=[x_edge[0], x_edge[-1], y_edge[0], y_edge[-1]], 
#                cmap='seismic', aspect='auto', norm=norm)

# cbar = plt.colorbar(im, label='|Kr/Sr|')
# cbar.ax.tick_params(labelsize=12)
# cbar.set_label('|Kr/Sr|',fontsize=12)

# # Set custom colorbar ticks to highlight the center at 1
# cbar.set_ticks([0.1, 1, 10])
# cbar.set_ticklabels(['0.1',  '1',  '10'])

# plt.xlim(-0.5, 0.5)
# plt.ylim(10, 180)

# plt.xticks(fontsize=12)
# plt.yticks(fontsize=12)
# plt.xlabel('Log10(Ma)',fontsize=12)
# plt.ylabel('Deflection Angle (degrees)',fontsize=12)
# plt.axhline(y=90, color='k')
# plt.axvline(x=0, color='k')
# plt.show()

#%%

def angledist_plot():
    x = np.log10(allbeta_par) # x axis
    y = allangles # y axis
    z = alldv_perp_norm  # variable for colorbar

    # Define bins
    bins=40
    x_bins = np.linspace(-0.6, 0.5, bins)
    y_bins = np.linspace(0, 180, bins)

    # Calculate the mean of 'z' for each bin
    mean_z, x_edge, y_edge, bin_number = binned_statistic_2d(x, y, z, statistic='count', bins=[x_bins, y_bins])

    # Calculate counts to identify zero bins
    count_z, _, _, _ = binned_statistic_2d(x, y, z, statistic='count', bins=[x_bins, y_bins])
    # Create a masked array where zero-count bins are masked
    mean_z_masked = np.ma.masked_where(count_z == 0, mean_z)

    # Plot the result
    fig, ax = plt.subplots(figsize=(10,7))

    # First, plot the zero-count bins in light grey
    zero_mask = count_z == 0
    if np.any(zero_mask):
        # Create an array filled with a dummy value for zero bins
        zero_array = np.full_like(mean_z, 1.0)  # Use middle value of your range
        zero_array[~zero_mask] = np.nan  # Set non-zero bins to NaN so they're transparent
        
        ax.imshow(zero_array.T, origin='lower', 
                extent=[x_edge[0], x_edge[-1], y_edge[0], y_edge[-1]], 
                cmap='gray', aspect='auto', vmin=0.8, vmax=1.2, alpha=0.4)

    # Then plot the actual data on top
    im = ax.imshow(mean_z_masked.T, origin='lower', 
                extent=[x_edge[0], x_edge[-1], y_edge[0], y_edge[-1]], 
                cmap='plasma', aspect='auto', vmin=0., vmax=1e3)


    cbar = plt.colorbar(im, label='counts')
    cbar.ax.tick_params(labelsize=12)
    cbar.set_label('counts',fontsize=12)
    plt.xlim(-0.6, 0.5)
    plt.ylim(10, 180)
    plt.xlabel('Log10(Ma)',fontsize=12)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.ylabel('Deflection Angle (degrees)',fontsize=12)
    plt.axhline(y=90, color='k')
    plt.axvline(x=0, color='k')
    plt.show()
    return
angledist_plot()


#%%

def Tparperp_plot():
    import matplotlib.colors as colors
    x = np.log10(filtered_mach) # x axis
    y = allangles # y axis
    z = allTparperp # variable for colorbar

    # Define bins
    bins=40
    x_bins = np.linspace(-0.6, 0.6, bins)
    y_bins = np.linspace(0, 180, bins)

    # Calculate the mean of 'z' for each bin
    mean_z, x_edge, y_edge, bin_number = binned_statistic_2d(x, y, z, statistic='mean', bins=[x_bins, y_bins])

    # Calculate counts to identify zero bins
    count_z, _, _, _ = binned_statistic_2d(x, y, z, statistic='count', bins=[x_bins, y_bins])
    # Create a masked array where zero-count bins are masked
    mean_z_masked = np.ma.masked_where(count_z == 0, mean_z)

    # Plot the result
    fig, ax = plt.subplots()

    # First, plot the zero-count bins in light grey
    zero_mask = count_z == 0
    if np.any(zero_mask):
        # Create an array filled with a dummy value for zero bins
        zero_array = np.full_like(mean_z, 1.0)  # Use middle value of your range
        zero_array[~zero_mask] = np.nan  # Set non-zero bins to NaN so they're transparent
        
        ax.imshow(zero_array.T, origin='lower', 
                extent=[x_edge[0], x_edge[-1], y_edge[0], y_edge[-1]], 
                cmap='gray', aspect='auto', vmin=0.8, vmax=1.2, alpha=0.4)

    # Create logarithmic normalization with range 0.1 to 10, centered at 1
    norm = colors.LogNorm(vmin=0.1, vmax=10)

    im = ax.imshow(mean_z_masked.T, origin='lower', 
                extent=[x_edge[0], x_edge[-1], y_edge[0], y_edge[-1]], 
                cmap='seismic', aspect='auto', norm=norm)

    cbar = plt.colorbar(im, label='|Kr/Sr|')
    cbar.ax.tick_params(labelsize=12)
    cbar.set_label('|T_par/T_perp|',fontsize=12)

    # Set custom colorbar ticks to highlight the center at 1
    cbar.set_ticks([0.1, 1, 10])
    cbar.set_ticklabels(['0.1',  '1',  '10'])

    plt.xlim(-0.6, 0.6)
    plt.ylim(10, 180)

    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.xlabel('Log10(Ma)',fontsize=12)
    plt.ylabel('Deflection Angle (degrees)',fontsize=12)
    plt.axhline(y=90, color='k')
    plt.axvline(x=0, color='k')
    plt.show()
    return

Tparperp_plot()



#%%

## Mach number vs Deflection Angle: dv/v Values 
x = np.log10(filtered_mach) # x axis
y = allangles # y axis
#z = (alln/1e6)*norm_factor # variable for colorbar
#z = abs(allSr-allSmags) # variable for colorbar
#z = (1/mu0)*allvmags*allBmags**2
#z=allvr*(allBt**2 + allBn**2)/mu0
z = alldv_norm_mag # variable for colorbar
z = allentropy

# Define bins
bins=40
x_bins = np.linspace(-0.6, 0.6, bins)
y_bins = np.linspace(0, 180, bins)

# Calculate the mean of 'z' for each bin
mean_z, x_edge, y_edge, bin_number = binned_statistic_2d(x, y, z, statistic='mean', bins=[x_bins, y_bins])

# Calculate counts to identify zero bins
count_z, _, _, _ = binned_statistic_2d(x, y, z, statistic='count', bins=[x_bins, y_bins])
# Create a masked array where zero-count bins are masked
mean_z_masked = np.ma.masked_where(count_z ==0, mean_z)
mean_z_masked = np.ma.masked_where(count_z <=0, mean_z)

# Plot the result
fig, ax = plt.subplots()

# First, plot the zero-count bins in light grey
zero_mask = count_z == 0
if np.any(zero_mask):
    # Create an array filled with a dummy value for zero bins
    zero_array = np.full_like(mean_z, 1.0)  # Use middle value of your range
    zero_array[~zero_mask] = np.nan  # Set non-zero bins to NaN so they're transparent
    
    ax.imshow(zero_array.T, origin='lower', 
              extent=[x_edge[0], x_edge[-1], y_edge[0], y_edge[-1]], 
              cmap='gray', aspect='auto', vmin=0.8, vmax=1.2, alpha=0.4)



im = ax.imshow(mean_z_masked.T, origin='lower', 
               extent=[x_edge[0], x_edge[-1], y_edge[0], y_edge[-1]], 
               cmap='seismic', aspect='auto', norm=norm)



# Then plot the actual data on top
im = ax.imshow(mean_z_masked.T, origin='lower', 
               extent=[x_edge[0], x_edge[-1], y_edge[0], y_edge[-1]], 
               cmap='seismic', aspect='auto', vmin=0, vmax= 2)

cbar = plt.colorbar(im, label='Temp')
cbar.ax.tick_params(labelsize=12)
cbar.set_label('dv/v',fontsize=12)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.xlim(-0.6, 0.6)
plt.ylim(10, 180)
plt.xlabel('Log10(Ma)',fontsize=12)
plt.ylabel('Deflection Angle (degrees)',fontsize=12)
plt.axhline(y=90, color='k')
plt.axvline(x=0, color='k')
plt.show()

#%%

## Mach number vs Deflection Angle: Specific Entropy Values

def sp_plot():
    x = np.log10(filtered_mach) # x axis
    y = allangles # y axis
    # z = alldv_perp_norm/alldv_norm_mag # variable for colorbar
    z = 1000*allentropy/(1.602e-19) #puts it back in ev/m^2 # variable for colorbar
    #z = alldv_perp_norm # variable for colorbar

    # Define bins
    bins=40
    x_bins = np.linspace(-0.6, 0.6, bins)
    y_bins = np.linspace(0, 180, bins)

    # Calculate the mean of 'z' for each bin
    mean_z, x_edge, y_edge, bin_number = binned_statistic_2d(x, y, z, statistic='mean', bins=[x_bins, y_bins])

    # Calculate counts to identify zero bins
    count_z, _, _, _ = binned_statistic_2d(x, y, z, statistic='count', bins=[x_bins, y_bins])
    # Create a masked array where zero-count bins are masked
    mean_z_masked = np.ma.masked_where(count_z == 0, mean_z)

    # Plot the result
    fig, ax = plt.subplots()

    # First, plot the zero-count bins in light grey
    zero_mask = count_z == 0
    if np.any(zero_mask):
        # Create an array filled with a dummy value for zero bins
        zero_array = np.full_like(mean_z, 1.0)  # Use middle value of your range
        zero_array[~zero_mask] = np.nan  # Set non-zero bins to NaN so they're transparent
        
        ax.imshow(zero_array.T, origin='lower', 
                extent=[x_edge[0], x_edge[-1], y_edge[0], y_edge[-1]], 
                cmap='gray', aspect='auto', vmin=0.8, vmax=1.2, alpha=0.4)

    # Then plot the actual data on top
    im = ax.imshow(mean_z_masked.T, origin='lower', 
                extent=[x_edge[0], x_edge[-1], y_edge[0], y_edge[-1]], 
                cmap='plasma', aspect='auto', vmin=0., vmax=0.3)

    plt.colorbar(im, label='Tn^(-2/3)')
    plt.xlim(-0.6, 0.6)
    plt.ylim(15, 160)
    plt.xlabel('Log10(Ma)')
    plt.ylabel('Deflection Angle (degrees)')
    plt.axhline(y=90, color='k')
    plt.axvline(x=0, color='k')
    plt.show()
    return
sp_plot()
#%%
print(np.nanstd(abs(allcrosshelicity[above30 & SBs])))
print(np.nanstd(abs(allcrosshelicity[above30 & wind])))
print(np.nanstd(abs(allcrosshelicity[below30 & SBs])))
print(np.nanstd(abs(allcrosshelicity[below30 & wind])))