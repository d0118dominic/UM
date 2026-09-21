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

me = 9.1094e-31 #kg
mi = 1837*me
mu0 = 1.2566370e-06  #;m kg / C^s
eps0 = 8.85e-12   # C^2/Nm^2

e = 1.602e-19 #C
Z = 1 # 1 for H+, 2 for He2+
gamma = 5/3
kb = 1.380649e-23



def reform(var):
	if not isinstance(var[1][0],np.ndarray):
		newvar = np.zeros(len(var[0]))
	elif isinstance(var[1][0][0],np.ndarray):
		newvar = np.zeros([len(var[0]),len(var[1][0]),len(var[1][0][0])])
	elif isinstance(var[1][0],np.ndarray):		
		newvar = np.zeros([len(var[0]),len(var[1][0])])
	for i in range(len(var[0])-1):
		newvar[i] = var[1][i]
	return newvar

# takes the timerange (trange) in [start,stop] format, and returns the number of minutes
def duration(trange):
    from datetime import datetime as dt
    start = dt.strptime(trange[0], '%Y-%m-%d/%H:%M')
    stop = dt.strptime(trange[1], '%Y-%m-%d/%H:%M')
    duration = stop-start
    duration_s = duration.total_seconds()
    duration_m = duration_s/60 #minutes
    duration_h = duration_m/60 #hours
    duration_d = duration_h/24 #days
    return duration_s, duration_m, duration_h,duration_d


# takes an array of some length and a timerange in [start,stop] format and determines how many elements of the array corresponds to 1 minute 
# (so an array of length 120 representing an hour of data would return a value of 2, because 2 data points of that array represent a minute)
def interval(timeax,trange):
	steps = len(timeax)
	seconds,minutes,hours,days = duration(trange)
	# second, minute, hour, day = np.ceil(steps/seconds),np.ceil(steps/minutes),np.ceil(steps/hours),np.ceil(steps/days)
	second, minute, hour, day = steps/seconds,steps/minutes,steps/hours,steps/days
	return second,minute,hour,day

# takes a timeseries and compute the mean (basically smooth out small fluctuations over some interval)
def get_mean(var, size):
    var = np.asarray(var, dtype=float)
    n = len(var)

    if size < 1:
        raise ValueError(f"Window size must be >= 1, got {size}")

    valid = ~np.isnan(var)
    var_filled = np.where(valid, var, 0.0)
    box = np.ones(size)

    # Use 'full' convolution and manually slice a centered length-n window,
    # so the output length is ALWAYS n, even if size > n.
    full_sum = np.convolve(var_filled, box, mode='full')
    full_count = np.convolve(valid.astype(float), box, mode='full')

    start = (size - 1) // 2
    sum_smoothed = full_sum[start:start + n]
    count_smoothed = full_count[start:start + n]

    with np.errstate(invalid='ignore', divide='ignore'):
        smoothed_var = sum_smoothed / count_smoothed
    smoothed_var[count_smoothed == 0] = np.nan

    return smoothed_var

# computes mean of each vector component (the n_minutes & minuteinterval stuff will be explained later)
def get_vecmean(vec,interval):   # Vector mean
	vec1_mean = get_mean(vec[:,0],interval)
	vec2_mean = get_mean(vec[:,1],interval)
	vec3_mean = get_mean(vec[:,2],interval)
	vec_mean = np.zeros_like(vec)
	for i in range(len(vec_mean)): vec_mean[i] = np.array([vec1_mean[i],vec2_mean[i],vec3_mean[i]])
	return vec_mean


def get_vecs (t_windows,vec):
    #vx = V[:, 0]
    vec_array = []
    # protect against going out of bounds from the array 
    for i in range(len(t_windows)):
        vec_small = get_vecmean(vec, int(np.round(t_windows[i]*secondint)))
        vec_array.append(vec_small)
    return np.array(vec_array)

def get_pm(B):
	pm = 0.5*(mu0**-1)*np.linalg.norm(B)**2
	return pm
def get_pth(n,T):
	pth = n*T
	return pth

def get_vals (t_windows,var):
    #vx = V[:, 0]
    var_array = []
    # protect against going out of bounds from the array 
    for i in range(len(t_windows)):
        var_small = get_mean(var, int(np.round(t_windows[i]*secondint)))
        var_array.append(var_small)
    return var_array

def get_scaldeltas (t_windows,scalvar):
    #vx = V[:, 0]
    delta_var_array,mean_var_array = [],[]
    # protect against going out of bounds from the array 
    for i in range(len(t_windows) - 1):

        # These are the time windows <B>_ts - <B>_tL / <B>_tL 
        T_small = t_windows[i]
        T_large = t_windows[i+1]
        # Compute moving-average magnetic field values using
        # the small and large time windows
        var_small = get_mean(scalvar, int(np.round(T_small*secondint)))
        var_large = get_mean(scalvar, int(np.round(T_large*secondint)))

        delta_var = (var_small - var_large)
        delta_var_array.append(delta_var)
        mean_var_array.append(var_large)

    return delta_var_array, mean_var_array

def get_vecdeltas (t_windows, vecvar):
    #vx = V[:, 0]
    delta_var_array,mean_var_array = [],[]
    # protect against going out of bounds from the array
	 
    for i in range(len(t_windows) - 1):
        # These are the time windows <B>_ts - <B>_tL / <B>_tL 
        T_small = t_windows[i]
        T_large = t_windows[i+1]
        # Compute moving-average magnetic field values using
        # the small and large time windows
        var_small = get_vecmean(vecvar, int(np.round(T_small*secondint)))
        var_large = get_vecmean(vecvar, int(np.round(T_large*secondint)))

        delta_var = (var_small - var_large)
        delta_var_array.append(delta_var)
        mean_var_array.append(var_large)

    return np.array(delta_var_array), np.array(mean_var_array)	
# given a magnetic field, density, and mass (usually use mi), this calculates the Alfven velocity
def get_va(B,n,m):
	va = np.linalg.norm(B)/np.sqrt(mu0*n*m)
	return va

def get_vth(T,m):
	vth = np.sqrt(kb*T/m)
	return vth
def filter_deflections(var,threshold):
	# newvar = np.array([])
	newvar = []
	for i in range(len(var)):
		if (abs(var[i])>=threshold):
			newvar.append(var[i]) 
		else: 
			pass
	return np.array(newvar)

def sw_km_to_Re(Vsw):
	#time 
	earth_rad_km = 6371
	Re_min = (Vsw * 60)/ earth_rad_km
	#distances 
	return Re_min

def get_parperps(n,T,B):  #B ant T tensor coord systems need to match for this
	trace = T[0] + T[1] + T[2]
	term1 = (T[0]*B[0]**2 + T[1]*B[1]**2 + T[2]*B[2]**2)/(np.linalg.norm(B)**2)
	term2 = 2*(T[3]*B[0]*B[1] + T[4]*B[0]*B[2] + T[5]*B[1]*B[2])/(np.linalg.norm(B)**2)
	Tpar = term1+term2
	Tperp=0.5*(trace-Tpar)
	Ppar = n*Tpar
	Pperp = n*Tperp
	return Tpar,Tperp,Ppar,Pperp

# Gets the angle between two vectors (radian version)
def get_angle(vec,meanvec):
	term1 = np.dot(vec,meanvec)
	term2 = np.dot(np.linalg.norm(vec),np.linalg.norm(meanvec))
	term3 = term1/term2
	term4 = np.arccos(term3) 
	return term4

# This gets all the angles across t_windows
def get_angles_t(t_windows, vecvar):
    angle_array = []

    for i in range(len(t_windows)-1):

        angle = np.zeros(len(vecvar))

        T_small = t_windows[i]
        T_large = t_windows[i+1]

        vec_small = get_vecmean(vecvar, int(T_small * minuteint))
        vec_large = get_vecmean(vecvar, int(T_large * minuteint))

        for j in range(len(vec_small)):
            angle[j] = get_angle(vec_small[j], vec_large[j])

        angle_array.append(angle)

    return angle_array


def get_counts(data_arr,threshold):
	counts_list = []
	#threshold = 0.05
	for i in range(len(d_windows) - 1):
		mask = data_arr[i]>=threshold
		number1 = len(data_arr[i][mask])
		counts_list.append(number1)
	return counts_list


# to use this cal the midpoint of the arr get the req counts within the data and choose which 
# edges you want to use 
def aggregate_counts_by_bounds(source_centers, source_counts, target_edges):
	aggregated = []
	for i in range(len(target_edges)-1):
		left = target_edges[i]
		right = target_edges[i+1]

		# Find source-bin centers inside the current interval

		mask = (source_centers >= left) & (source_centers < right)

		# Add the counts from those selected bins
		aggregated.append(int(np.sum(np.array(source_counts)[mask])))
	return np.array(aggregated)

def get_angles_t(t_windows, vecvar):
    angle_array = []
    for i in range(len(t_windows) - 1):
        angle = np.zeros(len(timeax))
        T_small = t_windows[i]
        T_large = t_windows[i+1]
        vec_small = get_vecmean(vecvar, int(T_small * minuteint))
        vec_large = get_vecmean(vecvar, int(T_large * minuteint))
        for j in range(len(vec_small)):
            angle[j] = get_angle(vec_small[j], vec_large[j])
        angle_array.append(angle)
    return angle_array

def y_error_bars(var,counts,p):
	y_error = []
	for var in counts:
		result = var * p
		error = y_error.append(result)
	return error 

recent_perihelia = [['2024-09-29/00:00', '2024-10-01/00:00'], #E21
					['2024-06-29/00:00', '2024-07-01/12:00'], #E20
					['2024-03-29/00:00','2024-03-31/00:00'], #E19
					['2023-12-24/00:00','2024-01-03/00:00'], #E18
					['2023-09-27/00:00','2023-09-29/00:00'], #E17
					['2023-06-19/00:00','2023-06-25/00:00']] #E16


# Specific Gradients for individual look?
event18 = ['2023-12-28/14:00','2023-12-28/16:00'], #E18
event20 = ['2024-06-29/00:00','2024-06-29/14:00'], #E20



event = ['2022-03-21/00:00','2022-12-14/00:00'], #E20
event = [recent_perihelia[3]]




# Alfven crossings (<2 hr)
# short_crossings = [['2021-11-21/21:00','2021-11-21/22:00'],
# 				   ['2021-08-10/00:15', '2021-08-10/00:45'], # Encounter 9 (some sub-Alfvenic)
#                    ['2023-12-29/01:30','2023-12-29/03:00'], # E18
#                    ['2024-03-29/22:00','2024-03-29/23:30'],  # E19 
#                    ['2024-06-29/11:00', '2024-06-29/13:00']] # E20


# eventlist = ['2023-12-28/00:00','2023-12-30/00:00']
# eventlist = [['202-21/00:00', '2024-06-24/00:00']]
eventlist = [['2025-03-21/00:00','2025-03-25/00:00']]
['2023-12-24/00:00','2024-01-03/00:00'], #E18
eventlist = [recent_perihelia[3]]

# Choose number of seconds for base, factors of 10 up until on order of a day or so
base_sec = 12
# t_windows = [base_sec,10*base_sec,100*base_sec, 1e3*base_sec]
t_windows = [base_sec,10*base_sec,100*base_sec, 1e3*base_sec, 1e4*base_sec]#,1e4*base_sec]
# 12 sec, 2 min, 20 min, 200 min (3.33 hrs), 33.3 hrs  








#%%
scalar_names = ['alln', 'allT', 'allTpar','allTperp',
                'allBmags','allvamags','allvmags','allrpsp']

vector_names = ['allv', 'allB', 'allTpar','allva','allvelocity_psp','allposition_psp']

for name in scalar_names:
    globals()[name] = np.empty(0)

for name in vector_names:
    globals()[name] = np.empty((0,3))




for i in range(len(eventlist)):
    trange = eventlist[i]
    # trange = event

    Bfld_vars = pyspedas.projects.psp.fields(trange=trange, level='l2', time_clip=True,no_update=True)
    swp_vars = pyspedas.projects.psp.spi(trange=trange,level='l3',get_support_data=True,time_clip=True,no_update=True)
    qtn_vars = pyspedas.projects.psp.fields(trange=trange,level='l3',datatype='sqtn_rfs_V1V2',time_clip=True,no_update=True)
    # alph_vars = pyspedas.projects.psp.spi(trange=trange,level='l3',datatype='sf0a_l3_mom',time_clip=True)
    # voltages_vars = pyspedas.projects.psp.fields(trange=trange, datatype='dfb_wf_dvdc', level='l2',time_clip=True)
    #On DC datatype: 'sqn_rfs_V1V2 has some kind of electron density & core temp, but looks weird
    # spc_vars = pyspedas.projects.psp.spc(trange=trange, datatype='l2', level='l2')

    # Reform all data to simple arrays and convert to SI units 
    B_name = 'psp_fld_l2_mag_RTN'
    vi_name = 'psp_spi_VEL_RTN_SUN'
    Bxyz_name = 'psp_spi_MAGF_INST'
    vxyz_name = 'psp_spi_VEL_INST'
    TiTensor_name = 'psp_spi_T_TENSOR_INST'
    Ti_name = 'psp_spi_TEMP'
    # ni_name = 'psp_spi_DENS'
    ni_name = 'electron_density'
    phivals_name = 'psp_spi_PHI_VALS'
    ephi_name = 'psp_spi_EFLUX_VS_PHI'
    # voltages_name = 'psp_fld_l2_dfb_wf_dVdc_sc'

    # Spacecraft Variables
    position_name = 'psp_spi_SUN_DIST'
    velocity_name = 'psp_spi_SC_VEL_RTN_SUN'

    interpvar_name = vi_name
    timeax = pytplot.get_data(interpvar_name).times
    secondint,minuteint,hourint,dayint = interval(timeax,trange)

    # Handling troublesome qtn indices
    for name in [B_name, Bxyz_name, vxyz_name, vi_name, Ti_name, ni_name, TiTensor_name, position_name]:
        data = pytplot.get_data(name)
        if data is None:
            print(f"WARNING: no data found for {name}")
            continue
        times = data.times
        _, unique_idx = np.unique(times, return_index=True)
        if len(unique_idx) < len(times):
            print(f"Removing {len(times)-len(unique_idx)} duplicate timestamps from {name}")
            pytplot.store_data(name, data={'x': times[unique_idx], 'y': data.y[unique_idx]})

    ##%%
    tinterpol(B_name,interpvar_name,newname='B')
    tinterpol(Bxyz_name,interpvar_name,newname='Bxyz')
    tinterpol(vxyz_name,interpvar_name,newname='vxyz')
    tinterpol(vi_name,interpvar_name,newname='vi')
    tinterpol(Ti_name,interpvar_name,newname='Ti')
    tinterpol(ni_name,interpvar_name,newname='ni')
    tinterpol(TiTensor_name,interpvar_name,newname='TiTensor')
    tinterpol(phivals_name,interpvar_name,newname='phivals')
    tinterpol(ephi_name,interpvar_name,newname='ephi')
    # tinterpol(voltages_name,interpvar_name,newname='voltages')
    tinterpol(position_name,interpvar_name,newname='position')
    tinterpol(velocity_name,interpvar_name,newname='velocity')

    Bvecs = 1e-9*reform(pytplot.get_data('B'))
    Bxyz = 1e-9*reform(pytplot.get_data('Bxyz'))
    vxyz = 1e3*reform(pytplot.get_data('vxyz'))
    vivecs = 1e3*reform(pytplot.get_data('vi'))
    ni = 1e6*reform(pytplot.get_data('ni'))
    Ti = 1.602e-19*reform(pytplot.get_data('Ti'))
    TiTensor = 1.602e-19*reform(pytplot.get_data('TiTensor')) #Comes in xyz
    phis = reform(pytplot.get_data('phivals'))
    ephi = reform(pytplot.get_data('ephi')).T
    # PiTensor = ni*TiTensor
    # voltages = reform(get_data('voltages'))
    position = reform(get_data('position'))/695700 # Solar radii
    velocity = reform(get_data('velocity')) # km/s
    vpsp,rpsp = velocity,position
    
    vshift = vpsp - vivecs # Difference between parker and bulk proton speed
    
    # va = np.zeros_like(ni)
    # for j in range(len(timeax)):
    #     va[j] = get_va(Bvecs[j],ni[j], mi)
   # Define a bunch of variables here which have no averaging or means. 
   # va, theta_PTB, energy fluxes, etc. 
   # After defining, they can be fed in to the delta functions 
    
    ma,va,vth,vsw,theta_ptb = np.zeros_like(ni),np.zeros_like(ni),np.zeros_like(ni),np.zeros_like(ni),np.zeros_like(ni)
    Tpar,Tperp,Ppar,Pperp,beta = np.zeros_like(ni),np.zeros_like(ni),np.zeros_like(ni),np.zeros_like(ni),np.zeros_like(ni)
    Pmag,Pth,beta,betapar,Tparperp = np.zeros_like(ni), np.zeros_like(ni),np.zeros_like(ni),np.zeros_like(ni),np.zeros_like(ni)
    
    for j in range(len(ni)):
        va[j] = get_va(Bvecs[j],ni[j],mi)
        vth[j] = get_vth(Ti[j],mi) 
        theta_ptb[j] = get_angle(vshift[j],Bvecs[j])
        ma[j] = np.linalg.norm(vivecs[j])/va[j]
        Pmag[j] = get_pm(Bvecs[j])
        Pth[j] = get_pth(ni[j],Ti[j])
        beta[j] = Pth[j]/Pmag[j]
        Tpar[j],Tperp[j],Ppar[j],Pperp[j] = get_parperps(ni[j],TiTensor[j],Bxyz[j])
        Tparperp[j] = Tpar[j]/Tperp[j]
        betapar[j] = Ppar[j]/Pmag[j]
        # Temp components
    vsw = vivecs[:,0]

    # Alginment parameter of shifted psp velocity.  
    # Ranges from 0 (orth to B) to 1 (aligned with B) 
    lpar_over_l = [abs(np.cos(theta)) for theta in theta_ptb]
    lperp_over_l = [abs(np.sin(theta)) for theta in theta_ptb]
    lperp_over_lpar = [abs(np.tan(theta)) for theta in theta_ptb]    
    
    # theta = np.arccos(brnorm)*180/np.pi
    ### OR ....

    # use the vector means to derive mor complex quantities after (maybe more efficient?)
    
    # Get Deltas and Means (vectors)
    # dBvecs,Bmeans = get_vecdeltas(t_windows,Bvecs)
    # dvvecs,vmeans = get_vecdeltas(t_windows,vivecs)
    # dvpsp,vpspmeans = get_vecdeltas(t_windows,vpsp)
    # dvshift,vshiftmeans = get_vecdeltas(t_windows,vshift)
    # temp anisotropies
    # Enthalpy Flux density vectors
    # Ion Agyrotropy??

   # Get Deltas and Means (Scalars)
    dn,nmeans = get_scaldeltas(t_windows,ni)
    dT,Tmeans = get_scaldeltas(t_windows,Ti)
    dva,vameans = get_scaldeltas(t_windows,va)
    dma,mameans = get_scaldeltas(t_windows,ma)
    dvth,vthmeans = get_scaldeltas(t_windows,vth)
    # # dvsw,vswmeans = get_scaldeltas(t_windows,vsw)

    lperp_over_l = get_vals(t_windows,lperp_over_l)
    lpar_over_l = get_vals(t_windows,lpar_over_l)
    lperp_over_lpar = get_vals(t_windows,lperp_over_lpar)

    v = get_vecs(t_windows,vshift) 
    Tpar = get_vals(t_windows,Tpar)
    Tperp = get_vals(t_windows,Tperp)
    Tparperp = get_vals(t_windows,Tparperp)
    betapar = get_vals(t_windows,betapar)
    # Something with Alphas? (idk their resolution tho)


    # Here, l correctly corresponds to small avg
    # Since we're using the t_windows starting from zero
    # vmeans already represents the large averages.
    l = np.zeros([len(t_windows[:-1]),len(ni)])
    for t in range(len(t_windows[:-1])):
        for j in range(len(ni)):
            l[t][j] = np.linalg.norm(v[t][j])*t_windows[t]



# What next??  What am I looking for?



    # l = [*t for t in t_windows[:-1]]
    # Means
    # Bvecs_mean = get_vecmean(Bvecs,int(np.round(n_sec*secondint)))
     
        # ma_mean = get_mean(v_mag/va,minutes*meaninterval)
# %
# 
#%%
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

def plot_vs_time(timeax, var, ylabel='', title='', ax=None):
    """Plot var against timeax (unix epoch seconds), with a proper datetime x-axis."""
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 4))

    time_dt = pd.to_datetime(timeax, unit='s')
    ax.plot(time_dt, var)

    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
    ax.xaxis.set_major_locator(mdates.AutoDateLocator())
    ax.set_xlabel('Time (UTC)')
    ax.set_ylabel(ylabel)
    ax.set_title(title)

    return ax

# Call it:
ax = plot_vs_time(timeax, Tpar[0], ylabel='n$_i$ [m$^{-3}$]', title='Ion Density')
plt.gcf().autofmt_xdate()
plt.show()


#%%

def plot_vs_time(timeax, var, ylabel='', title='', ax=None):
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 4))

    time_dt = pd.to_datetime(timeax, unit='s')
    ax.plot(time_dt, var)

    ax.xaxis.set_major_locator(mdates.AutoDateLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
    ax.set_xlabel('Time (UTC)')
    ax.set_ylabel(ylabel)
    ax.set_title(title)

    return ax

# Call it:
ax = plot_vs_time(timeax, Tpar[0], ylabel='n$_i$ [m$^{-3}$]', title='Ion Density')
plt.gcf().autofmt_xdate()
plt.show()
#%%

def plot_vs_time(timeax, vars_dict, ylabel='', title='', ax=None,
                  colors=('darkblue', 'limegreen', 'orange', 'red')):
    """Plot one or more variables against timeax (unix epoch seconds) on the same axes.
    
    vars_dict: dict mapping {label: array}, e.g. {'n_i': ni, 'n_e': ne}
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 7))

    time_dt = pd.to_datetime(timeax, unit='s')

    for i, (label, var) in enumerate(vars_dict.items()):
        ax.plot(time_dt, var, label=label, color=colors[i % len(colors)])

    ax.xaxis.set_major_locator(mdates.AutoDateLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
    ax.set_xlabel('Time (UTC)')
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend()

    return ax

# Call it:
ax = plot_vs_time(timeax, {'Small Scales (~34 $Mm$)': mameans[0], 'Medium Scales (~0.5 $R_s$)': mameans[1], 'Large Scales (~5 $R_s$)': mameans[2], 'Very Large Scales (~50 $R_s$)': mameans[3] },
                   ylabel='Alfvénic Mach Number ($M_A$)', title='Encounter 18')
ax.axhline(y=1, color='k')
ax.axhline(y=0.6, color='k',linestyle='dashed')
ax.axhline(y=1.6, color='k',linestyle='dashed')
ax.set_ylim([0.2, 10])
ax.set_xlabel('Date')
ax.set_yscale('log')
plt.gcf().autofmt_xdate()
plt.show()

#%%


#%%
ax = plot_vs_time(timeax, {'Small Scales (~34 $Mm$)': lperp_over_lpar[1], 'Medium Scales (~0.5 $R_s$)': lperp_over_lpar[2], 'Large Scales (~5 $R_s$)': lperp_over_lpar[3], 'Very Large Scales (~50 $R_s$)': lperp_over_lpar[4] },
                   ylabel='$\ell_\perp/\ell_\parallel$', title='Encounter 18')
ax.axhline(y=1, color='k')
# ax.axhline(y=0.6, color='k',linestyle='dashed')
# ax.axhline(y=1.6, color='k',linestyle='dashed')
ax.set_ylim([0.01, 100])
ax.set_xlabel('Date')
ax.set_yscale('log')
plt.gcf().autofmt_xdate()
plt.show()
#%%
ax = plot_vs_time(timeax, {'Small Scales (~34 $Mm$)': vthmeans[0], 'Medium Scales (~0.5 $R_s$)': vthmeans[1], 'Large Scales (~5 $R_s$)': vthmeans[2], 'Very Large Scales (~50 $R_s$)': vthmeans[3] },
                   ylabel='Thermal Velocity ($v_{th}$)', title='Encounter 18')
# ax.axhline(y=1, color='k')
# ax.axhline(y=0.6, color='k',linestyle='dashed')
# ax.axhline(y=1.6, color='k',linestyle='dashed')
# ax.set_ylim([0.2, 10])
ax.set_xlabel('Date')
ax.set_yscale('log')
plt.gcf().autofmt_xdate()
plt.show()
#%%
ax = plot_vs_time(timeax, {'Small Scales (~34 $Mm$)': 1/Tparperp[1], 'Medium Scales (~0.5 $R_s$)': 1/Tparperp[2], 'Large Scales (~5 $R_s$)': 1/Tparperp[3], 'Very Large Scales (~50 $R_s$)': 1/Tparperp[4] },
                   ylabel='$T_\perp/T_\parallel$', title='Encounter 18')
ax.axhline(y=1, color='k')
# ax.axhline(y=0.6, color='k',linestyle='dashed')
# ax.axhline(y=1.6, color='k',linestyle='dashed')
# ax.set_ylim([0.2, 10])
ax.set_xlabel('Date')
ax.set_yscale('log')
plt.ylim(0.1,10)
plt.gcf().autofmt_xdate()
plt.show()

#%%
ax = plot_vs_time(timeax, {'Small Scales (~34 $Mm$)': Tmeans[0], 'Medium Scales (~0.5 $R_s$)': Tmeans[1], 'Large Scales (~5 $R_s$)': Tmeans[2], 'Very Large Scales (~50 $R_s$)': Tmeans[3] },
                   ylabel='Temperature ($T$)', title='Encounter 18')
# ax.axhline(y=1, color='k')
# ax.axhline(y=0.6, color='k',linestyle='dashed')
# ax.axhline(y=1.6, color='k',linestyle='dashed')
# ax.set_ylim([0.2, 10])
ax.set_xlabel('Date')
ax.set_yscale('log')
plt.gcf().autofmt_xdate()
plt.show()


#%%

ax = plot_vs_time(timeax, {'Very Small Scales (~3.5 $Mm$)': abs(dva[0]/vameans[0]), 'Small Scales (~34 $Mm$)': abs(dva[1]/vameans[1]), 'Medium Scales (~0.5 $R_s$)': abs(dva[2]/vameans[2]), 'Large Scales (~5 $R_s$)': abs(dva[3]/vameans[3]) },
                   ylabel='$\delta v_A/v_A$', title='Encounter 18',colors=('grey','darkblue', 'limegreen', 'orange', 'red'))
ax.set_xlabel('Date')
ax.set_yscale('linear')
ax.set_ylim([0, 0.6])
plt.gcf().autofmt_xdate()
plt.show()

#%%

ax = plot_vs_time(timeax, {'Very Small Scales (~3.5 $Mm$)': abs(dvth[0]/vthmeans[0]), 'Small Scales (~34 $Mm$)': abs(dvth[1]/vthmeans[1]), 'Medium Scales (~0.5 $R_s$)': abs(dvth[2]/vthmeans[2]), 'Large Scales (~5 $R_s$)': abs(dvth[3]/vthmeans[3]) },
                   ylabel='$\delta v_{th}/v_{th}$', title='Encounter 18',colors=('grey','darkblue', 'limegreen', 'orange', 'red'))
ax.set_xlabel('Date')
ax.set_yscale('linear')
ax.set_ylim([0, 0.6])
plt.gcf().autofmt_xdate()
plt.show()
#%%

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import binned_statistic_2d
from matplotlib.colors import ListedColormap
import matplotlib.cm as cm
import matplotlib.colors as colors

# plt.hist2d(Tparperp[0][dvva[0]>=0.05],histtype='step',bins=100)
plt.hist2d(mameans[3][dvva[0]>=0.05],dvva[0][dvva[0]>0.05],bins=100)
# plt.hist(Tparperp[1][dvva[1]>=0.05],histtype='step',bins=100)
# plt.hist(Tparperp[2][dvva[2]>=0.05],histtype='step',bins=100)
# plt.hist(Tparperp[3][dvva[3]>=0.05],histtype='step',bins=100)
plt.yscale('log')
plt.xlim(0.1,1.6)













