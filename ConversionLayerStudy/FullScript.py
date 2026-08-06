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

# Gets the angle between two vectors
def get_angle(vec,meanvec):
	term1 = np.dot(vec,meanvec)
	term2 = np.dot(np.linalg.norm(vec),np.linalg.norm(meanvec))
	term3 = term1/term2
	term4 = np.arccos(term3)*180/np.pi
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
#%%
recent_perihelia = [['2024-09-29/00:00', '2024-10-01/00:00'], #E21
					['2024-06-29/00:00', '2024-07-01/12:00'], #E20
					['2024-03-29/00:00','2024-03-31/00:00'], #E19
					['2023-12-28/00:00','2023-12-30/00:00'], #E18
					['2023-09-27/00:00','2023-09-29/00:00'], #E17
					['2023-06-20/00:00','2023-06-24/06:00']] #E16


# Specific Gradients for individual look
event18 = ['2023-12-28/14:00','2023-12-28/16:00'], #E18
event20 = ['2024-06-29/09:00','2024-06-29/14:00'], #E20


# Alfven crossings (<2 hr)
# short_crossings = [['2021-11-21/21:00','2021-11-21/22:00'],
# 				   ['2021-08-10/00:15', '2021-08-10/00:45'], # Encounter 9 (some sub-Alfvenic)
#                    ['2023-12-29/01:30','2023-12-29/03:00'], # E18
#                    ['2024-03-29/22:00','2024-03-29/23:30'],  # E19 
#                    ['2024-06-29/11:00', '2024-06-29/13:00']] # E20


eventlist = [recent_perihelia[0]]
# eventlist = [['202-21/00:00', '2024-06-24/00:00']]

# Choose number of seconds for base, factors of 10 up until on order of a day or so
base_sec = 12
# t_windows = [base_sec,10*base_sec,100*base_sec, 1e3*base_sec]
t_windows = [10*base_sec,100*base_sec, 1e3*base_sec,1e4*base_sec]
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
    


    # va = np.zeros_like(ni)
    # for j in range(len(timeax)):
    #     va[j] = get_va(Bvecs[j],ni[j], mi)
   # Define a bunch of variables here which have no averaging or means. 
   # va, theta_PTB, energy fluxes, etc. 
   # After defining, they can be fed in to the delta functions 
    
    va,vth,vsw = np.zeros_like(ni),np.zeros_like(va),np.zeros_like(va)
    for j in range(len(ni)):
        va[j] = get_va(Bvecs[j],ni[j],mi)
        vth[j] = get_vth(Ti[j],mi)
    vsw = vivecs[:,0]
    ### OR ....

    # use the vector means to derive mor complex quantities after (maybe more efficient?)
    
    # Get Deltas and Means (vectors)
    # dBvecs,Bmeans = get_vecdeltas(t_windows,Bvecs)
    # dvvecs,vmeans = get_vecdeltas(t_windows,vivecs)
    # dvpsp,vpspmeans = get_vecdeltas(t_windows,vpsp)

    # Get Deltas and Means (Scalars)
    #dn,nmeans = get_scaldeltas(t_windows,ni)
    #dT,Tmeans = get_scaldeltas(t_windows,Ti)
    #dva,vameans = get_scaldeltas(t_windows,va)
    #dvth,vthmeans = get_scaldeltas(t_windows,vth)
    #dvsw,vswmeans = get_scaldeltas(t_windows,vsw)
 

    # Means
    # Bvecs_mean = get_vecmean(Bvecs,int(np.round(n_sec*secondint)))
     
        # ma_mean = get_mean(v_mag/va,minutes*meaninterval)
# %%