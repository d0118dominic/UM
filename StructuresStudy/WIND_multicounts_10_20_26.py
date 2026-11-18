
#%%
import pyspedas
import numpy as np
import matplotlib.pyplot as plt
import sys
import math 
from pyspedas.projects.kyoto import dst
from scipy.spatial import cKDTree

from pyspedas import get_data
from pyspedas import store_data
from pyspedas import timespan
from pyspedas import tplot_options

from pyspedas import tplot
from pyspedas import tinterpol

# List of useful physical constants
me = 9.1094e-31 #kg                                
mi = 1837*me
mu0 = 1.2566370e-06  #;m kg / C^2
eps0 = 8.85e-12   # C^2/Nm^2
e = 1.602e-19 #C
Z = 1 # 1 for H+, 2 for He2+
gamma = 5/3
kb = 1.380649e-23
Vsw = 400
#electron mass (kg)
#mi = 1837 * me       # ion mass ≈ proton mass (kg)
#mu0 = 1.2566370e-06  # permeability of free space / magnetic constant (H/m)
#eps0 = 8.85e-12      # permittivity of free space / electric constant (F/m)
#e = 1.602e-19        # elementary charge (C)
#Z = 1                # ion charge number (1 for H+, 2 for He2+, etc.)
#gamma = 5/3          # adiabatic index (monatomic ideal gas)
#kb = 1.380649e-23    # Boltzmann constant (J/K)



# List of useful physical constants
me = 9.1094e-31 #kg                                
mi = 1837*me
mu0 = 1.2566370e-06  #;m kg / C^2
eps0 = 8.85e-12   # C^2/Nm^2
e = 1.602e-19 #C
Z = 1 # 1 for H+, 2 for He2+
gamma = 5/3
kb = 1.380649e-23
Vsw = 400
#electron mass (kg)
#mi = 1837 * me       # ion mass ≈ proton mass (kg)
#mu0 = 1.2566370e-06  # permeability of free space / magnetic constant (H/m)
#eps0 = 8.85e-12      # permittivity of free space / electric constant (F/m)
#e = 1.602e-19        # elementary charge (C)
#Z = 1                # ion charge number (1 for H+, 2 for He2+, etc.)
#gamma = 5/3          # adiabatic index (monatomic ideal gas)
#kb = 1.380649e-23    # Boltzmann constant (J/K)

#%%


# takes a tplot variable and turns it into a simpler numpy array
def reform(var):
	if not isinstance(var[1][0],np.ndarray):
		newvar = np.zeros(len(var[0]))
	elif isinstance(var[1][0][0],np.ndarray):
		newvar = np.zeros([len(var[0]),len(var[1][0]),len(var[1][0][0])])
	elif isinstance(var[1][0],np.ndarray):		
		newvar = np.zeros([len(var[0]),len(var[1][0])])
	for i in range(len(var[0])-1):
		newvar[i] = var[1][i]
		# np.asarray(var[1], dtype=float)
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
	return duration_m


# takes an array of some length and a timerange in [start,stop] format and determines how many elements of the array corresponds to 1 minute 
# (so an array of length 120 representing an hour of data would return a value of 2, because 2 data points of that array represent a minute)
def minute_int(timeax,trange):
	steps = len(timeax)
	minutes = duration(trange)
	minute = np.ceil(steps/minutes)
	return int(minute)

# takes a timeseries and compute the mean (basically smooth out small fluctuations over some interval)
def get_mean(var,int):   
	box = np.ones(int)/int
	smoothed_var = np.convolve(var,box,mode='same')
	return smoothed_var

def get_mean(var, int):
    var = np.asarray(var, dtype=float)
    n = len(var)
    size = int  # (shadows the built-in int() inside this function — see note below)

    valid = ~np.isnan(var)
    var_filled = np.where(valid, var, 0.0)
    box = np.ones(size)

    # Use 'full' convolution + centered slice so output length always equals len(var),
    # even if the window is larger than the array itself.
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

# calculates the difference between 2 vectors, returns the difference vector and the normalized difference vector 
# (typically I was calculating the difference between a variable and its mean, thus the name)
def get_delta(vec,vec_mean):
	dvec = vec - vec_mean
	dvec_norm = dvec/np.linalg.norm(vec_mean)
	return dvec,dvec_norm

# calculates the difference between 2 scalars, returns the difference and normalized difference
def get_deltascalar(var,var_mean):
	dvar = var - var_mean
	dvar_norm = dvar/var_mean
	return dvar,dvar_norm

# given a magnetic field, density, and mass (usually use mi), this calculates the Alfven velocity
def get_va(B,n,m):
	va = np.linalg.norm(B)/np.sqrt(mu0*n*m)
	return va


def filter_deflections(var,threshold):
	# newvar = np.array([])
	newvar = []
	for i in range(len(var)):
		if (abs(var[i])>=threshold):
			newvar.append(var[i]) 
		else: 
			pass
	return np.array(newvar)

def get_bin_edge_labels(bins, bin_width):
	labels = []
	for left_edge in bins:
		right_edge = left_edge + bin_width
		labels.append(f'{left_edge}-{right_edge}')
	return labels

def get_nideltas (t_windows,ni):
	#vx = V[:, 0]
	delta_nivar_array = []
	# protect against going out of bounds from the array 
	for i in range(len(t_windows) - 1):

		# These are the time windows <B>_ts - <B>_tL / <B>_tL 
		T_small = t_windows[i]
		T_large = t_windows[i+1]
		# Compute moving-average magnetic field values using
		# the small and large time windows
		var_small = get_mean(ni, T_small)
		var_large = get_mean(ni, T_large)

		delta_var = (var_small - var_large) / var_large
		delta_nivar_array.append(delta_var)

	return np.abs(np.array(delta_nivar_array))


# helper fun for creating a array of ['0-200' .................] ~ building strings for the 
# x axis 
def get_span_labels(edges):
	labels = []
	for i in range(len(edges)-1):
		labels.append(f'{int(edges[i])}-{int(edges[i+1])}')
	return labels
#def Convert_sw_RE():
# Combine counts from smaller source bins into larger target intervals.
# For each large interval, find which source-bin centers fall inside it,
# then add those counts together and store the total.
#source_centers ~ [50,150,250,350,450,550,650,750,850,950]
#counts ~ of sw structures within the bins 
#target ~ [0, 200, 400, 600, 800, 1000]
def sw_km_to_Re(Vsw):
	#time 
	earth_rad_km = 6371
	Re_min = (Vsw * 60)/ earth_rad_km
	#distances 
	return Re_min
def Vsw_Tw (t_windows):
	#vx = V[:, 0]
	delta_B_array = []
	# protect against going out of bounds from the array 
	for i in range(len(t_windows) - 1):

		B_small = t_windows[i]
		B_large = t_windows[i+1]

		delta_B = (B_small - B_large) / B_large
		delta_B_array.append((B_small, B_large))

	return delta_B_array
#%%
def sw_km_to_Re(Vsw):	#time 
	earth_rad_km = 6371
	Re_min = (Vsw * 60)/ earth_rad_km
	#distances 

	return Re_min

def get_scaldeltas (t_windows,scalvar):
	#vx = V[:, 0]
	delta_var_array = []
	# protect against going out of bounds from the array 
	for i in range(len(t_windows) - 1):

		# These are the time windows <B>_ts - <B>_tL / <B>_tL 
		T_small = t_windows[i]
		T_large = t_windows[i+1]
		# Compute moving-average magnetic field values using
		# the small and large time windows
		var_small = get_mean(scalvar, T_small*minuteint)
		var_large = get_mean(scalvar, T_large*minuteint)

		delta_var = (var_small - var_large) / var_large
		delta_var_array.append(delta_var)

	return np.abs(np.array(delta_var_array))


def get_vecdeltas (t_windows, vecvar):
	#vx = V[:, 0]
	delta_var_array = []
	# protect against going out of bounds from the array
	 
	for i in range(len(t_windows) - 1):
		# These are the time windows <B>_ts - <B>_tL / <B>_tL 
		T_small = t_windows[i]
		T_large = t_windows[i+1]
		# Compute moving-average magnetic field values using
		# the small and large time windows
		var_small = get_vecmean(vecvar, T_small*minuteint)
		var_large = get_vecmean(vecvar, T_large*minuteint)

		delta_var = (var_small - var_large) / var_large
		delta_var_array.append(delta_var)

	return np.abs(delta_var_array)	
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

def get_counts_mask(data_arr,data_arr2,threshold):
	counts_list = []
	
	for i in range(len(t_windows) - 1):
		mask = (data_arr[i]>=threshold) & (data_arr2[i]>=threshold)
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


# New Functions

# Simple smoothing function to get list of averaged variables
def smooth_var (t_windows,scalvar):
	#vx = V[:, 0]
	var_array = []
	# protect against going out of bounds from the array 
	for i in range(len(t_windows)):
		T = t_windows[i]
		var = get_mean(scalvar, T*minuteint)
		var_array.append(var)
	return var_array

#Time derivative function.  Not used atm but may be later.  
def ddt_var (var):
	deltavar,deltavar_norm = np.zeros_like(var),np.zeros_like(var)
	for i in range(len(var)-1):
		deltavar[i] = np.abs(np.gradient(var[i]))
		deltavar_norm[i] = deltavar[i]/((t_windows[i+1]-t_windows[i])*var[i])
	return deltavar_norm

# The same basic delta function as before
# Computes the normalized difference between small and large mean
def delta_var (var):
	deltavar,deltavar_norm = np.zeros_like(var),np.zeros_like(var)
	for i in range(len(var)-1):
		deltavar[i] = np.abs(var[i]-var[i+1])
		deltavar_norm[i] = deltavar[i]/var[i+1]
	return deltavar_norm

# Takes a variable and a threshold
# If var>=threshold, sets value = 1, sets to zero otherwise
def make_binary (var,th):
	newvar = np.zeros_like(var)
	for i in range(len(var)):
		if var[i]>=th:
			newvar[i] = 1
		else:
			newvar[i] = 0
	return newvar


#(meant to handle the binary array)
# Determines if theres a change from 1 to 0 or 0 to 1
# This is designed to only identify the edges of structures
# the 1/2 factor is because abs(gradient of [0,1]) becomes [0.5,0.5]
# So each 'edge' has 2 nonzero elements
def get_struc_edges(var,label = '??'):
	newvar,counts = np.zeros_like(var),np.zeros(len(var))
	# print('\n Number of structure edges for '+label+':')
	# print('\n')
	for i in range(len(var)):
		newvar[i] = abs(np.gradient(var[i]))
		# newvar[i] = abs(np.diff(var[i]))
		counts[i] = len(newvar[i][newvar[i]!=0])
		# print(len(newvar[i][newvar[i]>0]))
	return newvar,counts/2

#%%

trange=['2019-01-01/00:00', '2019-01-30/00:00']
# trange=['2019-01-20/00:00', '2019-01-30/00:00']
swe_vars = pyspedas.projects.wind.swe(trange=trange,no_update=True,addmaster=False)
# denstiy 
mfi_vars = pyspedas.projects.wind.mfi(trange=trange,no_update=True,addmaster=False)
dst_vars = pyspedas.projects.kyoto.dst(trange=trange)
# tplot('kyoto_dst')
# 2. Load SYM-H (High-resolution Dst index from Kyoto)
# This loads 1-minute resolution sym_h data into pytplot
# pyspedas.kyoto.dst(trange=trange, datatype='symh')
# dst(trange=trange, datatype='symh')
# Define variable names for convenience
B_name='BGSE'
ni_name='N_elec'
vi_name = 'U_eGSE'
pos_name = 'DIST'
dst_name = 'kyoto_dst'

# Extract SYM-H values
# symh_data = p.get_data('sym_h')
# symh_time = symh_data.times
# symh_values = symh_data.y


# Interopolate variables to same cadence
interpvar_name = vi_name # Choose which variable to interpolate everything to 
interpdata = get_data(interpvar_name) # Get the interp var data
timeax = interpdata.times # Get the time axis of the interp var
minuteint = minute_int(timeax,trange) # Get the number of data points equal to one minute
tinterpol(B_name,interpvar_name,newname='B')
tinterpol(ni_name,interpvar_name,newname='ni') 
tinterpol(vi_name,interpvar_name,newname='vi') 
tinterpol(dst_name,interpvar_name,newname='dst') 

#%%
Bvecs = 1e-9*reform(get_data('B')) # Converting from nanotesla to Tesla
vivecs = 1e3*reform(get_data('vi')) #converting km/s to m/s
ni = 1e6*reform(get_data('ni')) # Converting from cm^-3 to m^-3
dst = 1e-9*reform(get_data('dst'))
#%%
# Get B and v magnitudes and Alfven speed
Bmag = np.zeros(len(timeax)) # this just creates the array of the correct length
vmag = np.zeros(len(timeax)) 
va = np.zeros(len(timeax))
# dstmag = 


for i in range(len(timeax)):
	Bmag[i] = np.linalg.norm(Bvecs[i]) 
	vmag[i] = np.linalg.norm(vivecs[i]) 
	va[i] = get_va(Bvecs[i],ni[i],mi)
	  

# Store as tplot variables
store_data('Bvecs', data = {'x':timeax,'y':1e9*Bvecs}) # with a factor to convert Tesla to nanotesla 
store_data('Bmag', data = {'x':timeax,'y':1e9*Bmag}) # with a factor to convert Tesla to nanotesla 
store_data('vmag', data = {'x':timeax,'y':1e-3*vmag}) 
store_data('va', data = {'x':timeax,'y':1e-3*va}) # with a factor to convert m/s to km/s


# %%
#Windspeeds in km/s
windspeed = 1e-3*np.nanmean(vmag)
windspeed_High = 1e-3*np.percentile(vmag,75) # Changed high speed to 75th percentile (instead of 90th)
windspeed_Low = 1e-3*np.percentile(vmag,25) # Changed low speed to 25th percentile (instead of 10th)

#Windspeeds in Re/min
sw_high = sw_km_to_Re(windspeed_High)#RE per min
sw_low = sw_km_to_Re(windspeed_Low)
sw_mean = sw_km_to_Re(windspeed)

print(r"Average Windspeed: ", windspeed, r'km/s')
print(r"Low: ", windspeed_Low, r'km/s')
print(r"High: ", windspeed_High, r'km/s')
#%%

# Set v to either sw_mean, sw_low, or sw_high

v_low = sw_low
v_high = sw_high
v_mean = sw_mean
# Define d_windows in Re/min
d_windows = np.array([4,10,100,1000,10000])

x = np.arange(len(d_windows)-1)


def make_d_windows():
	n = 101
	a = np.zeros(n)
	for i in range(len(a)):
		
		a[i] = 100*i
	a[0] = 4
	return a

def make_labels(d_windows):
    return [f"{d_windows[i]}-{d_windows[i+1]} $R_E$" for i in range(len(d_windows) - 1)]

# d_windows = make_d_windows()

# Get t_windows (list of integers)
t_windows_low,t_windows_high,t_windows_mean = [],[],[]
for i in d_windows: 
	t_windows_low.append(int(np.round(i/v_low)))
	t_windows_high.append(int(np.round(i/v_high)))
	t_windows_mean.append(int(np.round(i/v_mean)))

print("High Speed Wind: ",round(windspeed_High), 'km/s')
print("Low Speed Wind: ", round(windspeed_Low),'km/s \n')
print("Mean Wind Speed: ", round(windspeed), 'km/s')

t_windows=t_windows_mean

#%%

# 
def getvars(threshold,windtype='mean'):
	threshold=threshold
	# Just gets the smoothed variables on the different timescales
	if windtype=='mean':
		b = smooth_var(t_windows,Bmag)
		v = smooth_var(t_windows,vmag)
		n = smooth_var(t_windows,ni)
		# dst = smooth_var(t_windows,dst)
	if windtype=='fast':
		b = smooth_var(t_windows_high,Bmag)
		v = smooth_var(t_windows_high,vmag)
		n = smooth_var(t_windows_high,ni)
		# dst = smooth_var(t_windows_high,dst)
	if windtype=='slow':
		b = smooth_var(t_windows_low,Bmag)
		v = smooth_var(t_windows_low,vmag)
		n = smooth_var(t_windows_low,ni)
		# dst = smooth_var(t_windows_low,dst)
	# Returns the normalized delta (there's also a ddt function which is a time derivative) 
	# These arrays 
	db = delta_var(b)[:-1]
	dv = delta_var(v)[:-1]
	dn = delta_var(n)[:-1]
	# ddst = delta_var(dst)[:-1]

	# Sets all values above a threshold = 1, zero everywhere else
	# Just flags what is and isn't a structure'
	dbflags,dvflags,dnflags = np.zeros_like(db),np.zeros_like(dv),np.zeros_like(dn)
	# ddstflags = np.zeros_like(ddst)
	# print('For a threshold of '+str(int(threshold*100))+'%')
	for i in range(len(db)):
		dbflags[i] = make_binary(db[i],threshold)
		dvflags[i] = make_binary(dv[i],threshold)
		dnflags[i] = make_binary(dn[i],threshold)
		# ddstflags[i] = make_binary(ddst[i],threshold)

	# Make an array that looks for where the flagged structures begin or end
	# If the gradient>0, then it is counted
	# This should only catch the edges of structures 
	db_edges,bcounts = get_struc_edges(dbflags,label='db')
	# print('B Structures:')
	# print(bcounts/2)
	dv_edges,vcounts = get_struc_edges(dvflags,label='dv')
	# print('\nV Structures:')
	# print(vcounts/2)
	dn_edges,ncounts = get_struc_edges(dnflags,label ='dn')
	# print('\nN Structures:')
	# print(ncounts/2)
	# ddst_edges,ddstcounts = get_struc_edges(ddstflags,label = 'ddst')
	return threshold,b,v,n,db,dv,dn,dbflags,dvflags,dnflags,bcounts/2,vcounts/2,ncounts/2,2*db_edges,2*dv_edges,2*dn_edges


threshold=0.1
threshold_1 = 0.05
threshold,b_fast,v_fast,n_fast,db_fast,dv_fast,dn_fast,dbflags_fast,dvflags_fast,dnflags_fast,bcounts_fast,vcounts_fast,ncounts_fast,dbedges_fast,dvedges_fast,dnedges_fast = getvars(threshold=threshold,windtype='fast')
threshold,b_slow,v_slow,n_slow,db_slow,dv_slow,dn_slow,dbflags_slow,dvflags_slow,dnflags_slow,bcounts_slow,vcounts_slow,ncounts_slow,dbedges_slow,dvedges_slow,dnedges_slow = getvars(threshold=threshold,windtype='slow')
threshold,b,v,n,db,dv,dn,dbflags,dvflags,dnflags,bcounts,vcounts,ncounts,dbedges,dvedges,dnedges= getvars(threshold=threshold,windtype='mean')

threshold_1,b_1,v_1,n_1,db_1,dv_1,dn_1,dbflags_1,dvflags_1,dnflags_1,bcounts_1,vcounts_1,ncounts_1,dbedges_1,dvedges_1,dnedges_1= getvars(threshold=threshold_1,windtype='mean')
#%%

# def get_dst(windtype='mean'):
# if windtype=='mean':
ddst = smooth_var(t_windows,dst)
# ddst_fast = smooth_var(t_windows_high,dst)
# ddst_slow = smooth_var(t_windows_low,dst)

ddst = delta_var(ddst)[:-1]
# ddst_fast = delta_var(ddst_fast)[:-1]
# ddst_slow = delta_var(ddst_slow)[:-1]

def get_ddstflags(threshold=threshold):
	ddst = smooth_var(t_windows,dst)
	ddst = delta_var(ddst)[:-1]
	ddstflags = np.zeros_like(ddst)
	# print('For a threshold of '+str(int(threshold*100))+'%')
	for i in range(len(db)):
		ddstflags[i] = make_binary(ddst[i],threshold)
	return ddstflags

ddstflags = get_ddstflags()	

#%%



# Get a list of all 'events', including the start,stop times of each one
def get_startstops(flagvar):
	flagvar[0],flagvar[-1]=0,0 # Set first & last element to zero to handle diff issue
	diff = np.diff(flagvar) # sets 1 at first element of flag, -1 at last element

	# Get an array of events [start1,stop1], [start2,stop2], etc.
	start, stop = np.where(diff==1)[0]+1, np.where(diff==-1)[0]
	eventarr = np.array([start,stop]).T
	return eventarr

def get_ddstevents():
	ddst0events = get_startstops(ddstflags[0])
	ddst1events = get_startstops(ddstflags[1])
	ddst2events = get_startstops(ddstflags[2])
	ddst3events = get_startstops(ddstflags[3])
	ddstevents = [ddst0events,ddst1events,ddst2events,ddst3events]
	return ddstevents

# Make a function that can be called to return event lists 
def get_events(windtype='mean'):
	if windtype=='mean':
		db0events = get_startstops(dbflags[0])
		db1events = get_startstops(dbflags[1])
		db2events = get_startstops(dbflags[2])
		db3events = get_startstops(dbflags[3])
		dbevents = [db0events,db1events,db2events,db3events]

		dv0events = get_startstops(dvflags[0])
		dv1events = get_startstops(dvflags[1])
		dv2events = get_startstops(dvflags[2])
		dv3events = get_startstops(dvflags[3])
		dvevents = [dv0events,dv1events,dv2events,dv3events]

		dn0events = get_startstops(dnflags[0])
		dn1events = get_startstops(dnflags[1])
		dn2events = get_startstops(dnflags[2])
		dn3events = get_startstops(dnflags[3])
		dnevents = [dn0events,dn1events,dn2events,dn3events]
		return dbevents,dvevents,dnevents
	elif windtype=='fast':
		db0events_fast = get_startstops(dbflags_fast[0])
		db1events_fast = get_startstops(dbflags_fast[1])
		db2events_fast = get_startstops(dbflags_fast[2])
		db3events_fast = get_startstops(dbflags_fast[3])
		dbevents_fast = [db0events_fast,db1events_fast,db2events_fast,db3events_fast]
		
		dv0events_fast = get_startstops(dvflags_fast[0])
		dv1events_fast = get_startstops(dvflags_fast[1])
		dv2events_fast = get_startstops(dvflags_fast[2])
		dv3events_fast = get_startstops(dvflags_fast[3])
		dvevents_fast = [dv0events_fast,dv1events_fast,dv2events_fast,dv3events_fast]
		
		dn0events_fast = get_startstops(dnflags_fast[0])
		dn1events_fast = get_startstops(dnflags_fast[1])
		dn2events_fast = get_startstops(dnflags_fast[2])
		dn3events_fast = get_startstops(dnflags_fast[3])
		dnevents_fast = [dn0events_fast,dn1events_fast,dn2events_fast,dn3events_fast]
		return dbevents_fast,dvevents_fast,dnevents_fast

	elif windtype =='slow':
			db0events_slow = get_startstops(dbflags_slow[0])
			db1events_slow = get_startstops(dbflags_slow[1])
			db2events_slow = get_startstops(dbflags_slow[2])
			db3events_slow = get_startstops(dbflags_slow[3])
			dbevents_slow = [db0events_slow,db1events_slow,db2events_slow,db3events_slow]
			
			dv0events_slow = get_startstops(dvflags_slow[0])
			dv1events_slow = get_startstops(dvflags_slow[1])
			dv2events_slow = get_startstops(dvflags_slow[2])
			dv3events_slow = get_startstops(dvflags_slow[3])
			dvevents_slow = [dv0events_slow,dv1events_slow,dv2events_slow,dv3events_slow]
			
			dn0events_slow = get_startstops(dnflags_slow[0])
			dn1events_slow = get_startstops(dnflags_slow[1])
			dn2events_slow = get_startstops(dnflags_slow[2])
			dn3events_slow= get_startstops(dnflags_slow[3])
			dnevents_slow = [dn0events_slow,dn1events_slow,dn2events_slow,dn3events_slow]



			return dbevents_slow,dvevents_slow,dnevents_slow


	
dbevents,dvevents,dnevents = get_events(windtype='mean')
dbevents_fast,dvevents_fast,dnevents_fast = get_events(windtype='fast')
dbevents_slow,dvevents_slow,dnevents_slow = get_events(windtype='slow')

ddstevents = get_ddstevents()


#%%

ddst0events = get_startstops(ddstflags[0])
ddst1events = get_startstops(ddstflags[1])
ddst2events = get_startstops(ddstflags[2])
ddst3events = get_startstops(ddstflags[3])
ddstevents = [ddst0events,ddst1events,ddst2events,ddst3events]
ddstcounts = [len(i) for i in ddstevents]


db0events = get_startstops(dbflags[0])
db1events = get_startstops(dbflags[1])
db2events = get_startstops(dbflags[2])
db3events = get_startstops(dbflags[3])
dbevents = [db0events,db1events,db2events,db3events]
bcounts = [len(i) for i in dbevents]

dv0events = get_startstops(dvflags[0])
dv1events = get_startstops(dvflags[1])
dv2events = get_startstops(dvflags[2])
dv3events = get_startstops(dvflags[3])
dvevents = [dv0events,dv1events,dv2events,dv3events]
vcounts = [len(i) for i in dvevents]

dn0events = get_startstops(dnflags[0])
dn1events = get_startstops(dnflags[1])
dn2events = get_startstops(dnflags[2])
dn3events = get_startstops(dnflags[3])
dnevents = [dn0events,dn1events,dn2events,dn3events]
ncounts = [len(i) for i in dnevents]
## %5
db0events_1 = get_startstops(dbflags_1[0])
db1events_1 = get_startstops(dbflags_1[1])
db2events_1 = get_startstops(dbflags_1[2])
db3events_1 = get_startstops(dbflags_1[3])
dbevents_1 = [db0events_1,db1events_1,db2events_1,db3events_1]
bcounts_1 = [len(i) for i in dbevents_1]

dv0events_1 = get_startstops(dvflags_1[0])
dv1events_1 = get_startstops(dvflags_1[1])
dv2events_1 = get_startstops(dvflags_1[2])
dv3events_1 = get_startstops(dvflags_1[3])
dvevents_1 = [dv0events_1,dv1events_1,dv2events_1,dv3events_1]
vcounts_1 = [len(i) for i in dvevents_1]

dn0events_1 = get_startstops(dnflags_1[0])
dn1events_1 = get_startstops(dnflags_1[1])
dn2events_1 = get_startstops(dnflags_1[2])
dn3events_1 = get_startstops(dnflags_1[3])
dnevents_1 = [dn0events_1,dn1events_1,dn2events_1,dn3events_1]
ncounts_1 = [len(i) for i in dnevents_1]


###############
db0events_fast = get_startstops(dbflags_fast[0])
db1events_fast = get_startstops(dbflags_fast[1])
db2events_fast = get_startstops(dbflags_fast[2])
db3events_fast = get_startstops(dbflags_fast[3])
dbevents_fast = [db0events_fast,db1events_fast,db2events_fast,db3events_fast]
bcounts_fast = [len(i) for i in dbevents_fast]

dv0events_fast = get_startstops(dvflags_fast[0])
dv1events_fast = get_startstops(dvflags_fast[1])
dv2events_fast = get_startstops(dvflags_fast[2])
dv3events_fast = get_startstops(dvflags_fast[3])
dvevents_fast = [dv0events_fast,dv1events_fast,dv2events_fast,dv3events_fast]
vcounts_fast = [len(i) for i in dvevents_fast]

dn0events_fast = get_startstops(dnflags_fast[0])
dn1events_fast = get_startstops(dnflags_fast[1])
dn2events_fast= get_startstops(dnflags_fast[2])
dn3events_fast = get_startstops(dnflags_fast[3])
dnevents_fast = [dn0events_fast,dn1events_fast,dn2events_fast,dn3events_fast]
ncounts_fast = [len(i) for i in dnevents_fast]


db0events_slow = get_startstops(dbflags_slow[0])
db1events_slow = get_startstops(dbflags_slow[1])
db2events_slow = get_startstops(dbflags_slow[2])
db3events_slow = get_startstops(dbflags_slow[3])
dbevents_slow = [db0events_slow,db1events_slow,db2events_slow,db3events_slow]
bcounts_slow = [len(i) for i in dbevents_slow]

dv0events_slow = get_startstops(dvflags_slow[0])
dv1events_slow = get_startstops(dvflags_slow[1])
dv2events_slow = get_startstops(dvflags_slow[2])
dv3events_slow = get_startstops(dvflags_slow[3])
dvevents_slow = [dv0events_slow,dv1events_slow,dv2events_slow,dv3events_slow]
vcounts_slow = [len(i) for i in dvevents_slow]

dn0events_slow = get_startstops(dnflags_slow[0])
dn1events_slow = get_startstops(dnflags_slow[1])
dn2events_slow = get_startstops(dnflags_slow[2])
dn3events_slow = get_startstops(dnflags_slow[3])
dnevents_slow = [dn0events_slow,dn1events_slow,dn2events_slow,dv3events_slow]
ncounts_slow = [len(i) for i in dnevents_slow]


#%%

Transpose_arr1= np.array([bcounts_fast,bcounts_slow]).T


maxcounts = [int(np.max(Transpose_arr1[i])) for i in range (len(bcounts))]
mincounts = [int(np.min(Transpose_arr1[i])) for i in range (len(bcounts))]

#nanmax skips the NaN in the Trasnpose arr  finds the max
#nanmin skips the NaN in the arr but finds the min 

span = np.array(maxcounts) - np.array(mincounts) 
err = span/2
center = maxcounts-err

n_Transpose_arr1= np.array([ncounts_fast,ncounts_slow]).T

n_maxcounts = [int(np.max(n_Transpose_arr1[i])) for i in range (len(bcounts))]
n_mincounts = [int(np.min(n_Transpose_arr1[i])) for i in range (len(bcounts))]

#nanmax skips the NaN in the Trasnpose arr  finds the max
#nanmin skips the NaN in the arr but finds the min 

n_span = np.array(n_maxcounts) - np.array(n_mincounts) 
n_err = n_span/2
n_center = n_maxcounts-n_err



v_Transpose_arr1= np.array([vcounts_fast,vcounts_slow]).T

v_maxcounts = [int(np.max(v_Transpose_arr1[i])) for i in range (len(bcounts))]
v_mincounts = [int(np.min(v_Transpose_arr1[i])) for i in range (len(bcounts))]

#nanmax skips the NaN in the Trasnpose arr  finds the max
#nanmin skips the NaN in the arr but finds the min 

v_span = np.array(v_maxcounts) - np.array(v_mincounts) 
v_err = v_span/2
v_center = v_maxcounts-v_err



# Plot of single-variable strucures
print('B Structures: \n' + str(bcounts) + '\n')
print('V Structures: \n' + str(vcounts)+'\n')
print('N Structures: \n' + str(ncounts)+'\n')
# x = np.arange(len(bcounts))
labels = np.array(make_labels(d_windows)) if 'd_windows' in locals() else [f"W{i+1}" for i in range(len(d_windows))]
width = 0.25

plt.bar(x,ncounts_1, label=fr'N B V Structures at {threshold_1 *100:.0f}%', width=width, edgecolor='k', linestyle='--', fill=False)
plt.bar(x-width, bcounts_1, width=width, edgecolor='k', linestyle='--', fill=False)
plt.bar(x+width, vcounts_1, width=width, edgecolor='k', linestyle='--', fill=False)

plt.bar(x-width,bcounts,label = 'B Structures',width=width,color='b')
plt.bar(x,ncounts,label = 'N Structures',width=width,color='g')
plt.bar(x+width,vcounts,label = 'V Structures',width=width,color='red')
plt.errorbar(x-width,center,yerr=err,fmt='none',ecolor='k',elinewidth=2, capsize=4)
plt.errorbar(x,n_center,yerr= n_err,fmt='none',ecolor='k',elinewidth=2, capsize=4)
plt.errorbar(x+width,v_center,yerr=v_err,fmt='none',ecolor='k',elinewidth=2, capsize=4)
plt.xticks(x,labels)
plt.ylabel('Counts')
plt.ylim(1e0,1e4)
plt.title('Single-Variable Structures at '+str(int(threshold*100))+'% threshold')
plt.legend()
plt.yscale('log')



# %%
from scipy.spatial import cKDTree


def tolerance(arr1, arr2, tol=0):
    tree = cKDTree(arr1)
    # nearest neighbor in arr1 for every row of arr2, Chebyshev distance (p=inf)
    # means BOTH components must be within tol
    dist, _ = tree.query(arr2, k=1, p=np.inf,
                         distance_upper_bound=np.nextafter(tol, np.inf))
    return arr2[dist <= tol]


# Tolerance for each is based on the length of avg window
# Since t_windows is in minutes, need to multiply by minuteint
t = [minuteint*i for i in t_windows]

# COMMENT:  I think we need a different t for the tolerances on the fast and slow lists
# Luckily, those t_windows arrays are already defined (as low and high for slow and fast)
t_slow = [minuteint*i for i in t_windows_low]
t_fast = [minuteint*i for i in t_windows_high]


bn0events = tolerance(db0events,dn0events,tol=int(t[0]/2))
bn1events = tolerance(db1events,dn1events,tol=int(t[1]/2))
bn2events = tolerance(db2events,dn2events,tol=int(t[2]/2))
bn3events = tolerance(db3events,dn3events,tol=int(t[3]/2))

bn0events_fast = tolerance(db0events_fast,dn0events_fast,tol=int(t_fast[0]/2))
bn1events_fast = tolerance(db1events_fast,dn1events_fast,tol=int(t_fast[1]/2))
bn2events_fast = tolerance(db2events_fast,dn2events_fast,tol=int(t_fast[2]/2))
bn3events_fast = tolerance(db3events_fast,dn3events_fast,tol=int(t_fast[3]/2))

bn0events_slow = tolerance(db0events_slow ,dn0events_slow ,tol=int(t_slow[0]/2))
bn1events_slow = tolerance(db1events_slow ,dn1events_slow ,tol=int(t_slow[1]/2))
bn2events_slow  = tolerance(db2events_slow ,dn2events_slow ,tol=int(t_slow[2]/2))
bn3events_slow  = tolerance(db3events_slow ,dn3events_slow ,tol=int(t_slow[3]/2))

bv0events = tolerance(db0events,dv0events,tol=int(t[0]/2))
bv1events = tolerance(db1events,dv1events,tol=int(t[1]/2))
bv2events = tolerance(db2events,dv2events,tol=int(t[2]/2))
bv3events = tolerance(db3events,dv3events,tol=int(t[3]/2))

bv0events_fast = tolerance(db0events_fast,dv0events_fast,tol=int(t_fast[0]/2))
bv1events_fast = tolerance(db1events_fast,dv1events_fast,tol=int(t_fast[1]/2))
bv2events_fast = tolerance(db2events_fast,dv2events_fast,tol=int(t_fast[2]/2))
bv3events_fast = tolerance(db3events_fast,dv3events_fast,tol=int(t_fast[3]/2))

bv0events_slow = tolerance(db0events_slow ,dv0events_slow ,tol=int(t_slow[0]/2))
bv1events_slow = tolerance(db1events_slow ,dv1events_slow ,tol=int(t_slow[1]/2))
bv2events_slow  = tolerance(db2events_slow ,dv2events_slow ,tol=int(t_slow[2]/2))
bv3events_slow  = tolerance(db3events_slow ,dv3events_slow ,tol=int(t_slow[3]/2))

# VN events (NEW)
vn0events = tolerance(dv0events,dn0events,tol=int(t[0]/2))
vn1events = tolerance(dv1events,dn1events,tol=int(t[1]/2))
vn2events = tolerance(dv2events,dn2events,tol=int(t[2]/2))
vn3events = tolerance(dv3events,dn3events,tol=int(t[3]/2))

vn0events_fast = tolerance(dv0events_fast,dn0events_fast,tol=int(t_fast[0]/2))
vn1events_fast = tolerance(dv1events_fast,dn1events_fast,tol=int(t_fast[1]/2))
vn2events_fast = tolerance(dv2events_fast,dn2events_fast,tol=int(t_fast[2]/2))
vn3events_fast = tolerance(dv3events_fast,dn3events_fast,tol=int(t_fast[3]/2))

vn0events_slow = tolerance(dv0events_slow ,dn0events_slow ,tol=int(t_slow[0]/2))
vn1events_slow = tolerance(dv1events_slow ,dn1events_slow ,tol=int(t_slow[1]/2))
vn2events_slow  = tolerance(dv2events_slow ,dn2events_slow ,tol=int(t_slow[2]/2))
vn3events_slow  = tolerance(dv3events_slow ,dn3events_slow ,tol=int(t_slow[3]/2))

# AND events 
bvn0events = tolerance(dv0events,bn0events,tol=int(t[0]/2))
bvn1events = tolerance(dv1events,bn1events,tol=int(t[1]/2))
bvn2events = tolerance(dv2events,bn2events,tol=int(t[2]/2))
bvn3events = tolerance(dv3events,bn3events,tol=int(t[3]/2))

bvn0events_fast = tolerance(dv0events_fast,bn0events_fast,tol=int(t_fast[0]/2))
bvn1events_fast = tolerance(dv1events_fast,bn1events_fast,tol=int(t_fast[1]/2))
bvn2events_fast = tolerance(dv2events_fast,bn2events_fast,tol=int(t_fast[2]/2))
bvn3events_fast = tolerance(dv3events_fast,bn3events_fast,tol=int(t_fast[3]/2))

bvn0events_slow = tolerance(dv0events_slow,bn0events_slow,tol=int(t_slow[0]/2))
bvn1events_slow = tolerance(dv1events_slow,bn1events_slow,tol=int(t_slow[1]/2))
bvn2events_slow = tolerance(dv2events_slow,bn2events_slow,tol=int(t_slow[2]/2))
bvn3events_slow = tolerance(dv3events_slow,bn3events_slow,tol=int(t_slow[3]/2))


#5%
bn0events_1 = tolerance(db0events_1,dn0events_1,tol=int(t[0]/2))
bn1events_1 = tolerance(db1events_1,dn1events_1,tol=int(t[1]/2))
bn2events_1 = tolerance(db2events_1,dn2events_1,tol=int(t[2]/2))
bn3events_1 = tolerance(db3events_1,dn3events_1,tol=int(t[3]/2))

bv0events_1 = tolerance(db0events_1,dv0events_1,tol=int(t[0]/2))
bv1events_1 = tolerance(db1events_1,dv1events_1,tol=int(t[1]/2))
bv2events_1 = tolerance(db2events_1,dv2events_1,tol=int(t[2]/2))
bv3events_1 = tolerance(db3events_1,dv3events_1,tol=int(t[3]/2))

# VN events, 5% (NEW)
vn0events_1 = tolerance(dv0events_1,dn0events_1,tol=int(t[0]/2))
vn1events_1 = tolerance(dv1events_1,dn1events_1,tol=int(t[1]/2))
vn2events_1 = tolerance(dv2events_1,dn2events_1,tol=int(t[2]/2))
vn3events_1 = tolerance(dv3events_1,dn3events_1,tol=int(t[3]/2))

bvn0events_1 = tolerance(dv0events_1,bn0events_1,tol=int(t[0]/2))
bvn1events_1 = tolerance(dv1events_1,bn1events_1,tol=int(t[1]/2))
bvn2events_1 = tolerance(dv2events_1,bn2events_1,tol=int(t[2]/2))
bvn3events_1 = tolerance(dv3events_1,bn3events_1,tol=int(t[3]/2))




Transpose_arr1= np.array([bcounts_fast,bcounts_slow]).T

maxcounts = [int(np.max(Transpose_arr1[i])) for i in range (len(bcounts))]
mincounts = [int(np.min(Transpose_arr1[i])) for i in range (len(bcounts))]

#nanmax skips the NaN in the Trasnpose arr  finds the max
#nanmin skips the NaN in the arr but finds the min 

span = np.array(maxcounts) - np.array(mincounts) 
err = span/2
center = maxcounts-err

#5%
bnevents_tol_1 = [len(bn0events_1),len(bn1events_1),len(bn2events_1),len(bn3events_1)]
bvevents_tol_1 = [len(bv0events_1),len(bv1events_1),len(bv2events_1),len(bv3events_1)]
vnevents_tol_1 = [len(vn0events_1),len(vn1events_1),len(vn2events_1),len(vn3events_1)]   # NEW
bvnevents_tol_1 = [len(bvn0events_1),len(bvn1events_1),len(bvn2events_1),len(bvn3events_1)]

bvnevents_tol = [len(bvn0events),len(bvn1events),len(bvn2events),len(bvn3events)]
bvnevents_fast_tol = [len(bvn0events_fast),len(bvn1events_fast),len(bvn2events_fast),len(bvn3events_fast)]
bvnevents_slow_tol = [len(bvn0events_slow),len(bvn1events_slow),len(bvn2events_slow),len(bvn3events_slow)]

bnevents_tol = [len(bn0events),len(bn1events),len(bn2events),len(bn3events)]
bvevents_tol = [len(bv0events),len(bv1events),len(bv2events),len(bv3events)]
vnevents_tol = [len(vn0events),len(vn1events),len(vn2events),len(vn3events)]             # NEW
bvnevents_tol = [len(bvn0events),len(bvn1events),len(bvn2events),len(bvn3events)]
# notolerance = [len(notol0),len(notol1),len(notol2),len(notol3)]

bnevents_fast_tol = [len(bn0events_fast),len(bn1events_fast),len(bn2events_fast),len(bn3events_fast)]
bvevents_fast_tol = [len(bv0events_fast),len(bv1events_fast),len(bv2events_fast),len(bv3events_fast)]
vnevents_fast_tol = [len(vn0events_fast),len(vn1events_fast),len(vn2events_fast),len(vn3events_fast)]   # NEW
bvnevents_fast_tol = [len(bvn0events_fast),len(bvn1events_fast),len(bvn2events_fast),len(bvn3events_fast)]
# notolerance_fast = [len(notol0_fast),len(notol1_fast),len(notol2_fast),len(notol3_fast)]

bnevents_slow_tol = [len(bn0events_slow),len(bn1events_slow),len(bn2events_slow),len(bn3events_slow)]
bvevents_slow_tol = [len(bv0events_slow),len(bv1events_slow),len(bv2events_slow),len(bv3events_slow)]
vnevents_slow_tol = [len(vn0events_slow),len(vn1events_slow),len(vn2events_slow),len(vn3events_slow)]   # NEW
bvnevents_slow_tol = [len(bvn0events_slow),len(bvn1events_slow),len(bvn2events_slow),len(bvn3events_slow)]
# notolerance_slow = [len(notol0_slow),len(notol1_slow),len(notol2_slow),len(notol3_slow)]

# print('Event counts (exact match): ' + str(notolerance))
print('\nEvent counts (scale-dependent tolerance): ' + str(bnevents_tol))


Transpose_arr1_tol= np.array([bnevents_fast_tol,bnevents_slow_tol]).T

maxcounts_tol = [int(np.max(Transpose_arr1_tol[i])) for i in range (len(bnevents_tol))]
mincounts_tol = [int(np.min(Transpose_arr1_tol[i])) for i in range (len(bnevents_tol))]

#nanmax skips the NaN in the Trasnpose arr  finds the max
#nanmin skips the NaN in the arr but finds the min 

span_tol = np.array(maxcounts_tol) - np.array(mincounts_tol) 
err_tol = span_tol/2
center_tol = maxcounts_tol-err_tol

Transpose_bv_tol= np.array([bvevents_fast_tol,bvevents_slow_tol]).T

maxcountsbv_tol = [int(np.max(Transpose_bv_tol[i])) for i in range (len(bnevents_tol))]
mincountsbv_tol = [int(np.min(Transpose_bv_tol[i])) for i in range (len(bnevents_tol))]

#nanmax skips the NaN in the Trasnpose arr  finds the max
#nanmin skips the NaN in the arr but finds the min 

span_tolbv = np.array(maxcountsbv_tol) - np.array(mincountsbv_tol) 
err_tolbv = span_tolbv/2
center_tolbv = maxcountsbv_tol-err_tolbv

# NEW: VN fast/slow spread for error bars
Transpose_vn_tol= np.array([vnevents_fast_tol,vnevents_slow_tol]).T

maxcountsvn_tol = [int(np.max(Transpose_vn_tol[i])) for i in range (len(vnevents_tol))]
mincountsvn_tol = [int(np.min(Transpose_vn_tol[i])) for i in range (len(vnevents_tol))]

span_tolvn = np.array(maxcountsvn_tol) - np.array(mincountsvn_tol) 
err_tolvn = span_tolvn/2
center_tolvn = maxcountsvn_tol-err_tolvn

Transpose_bvn_tol= np.array([bvnevents_fast_tol,bvnevents_slow_tol]).T

maxcountsbvn_tol = [int(np.max(Transpose_bvn_tol[i])) for i in range (len(bvnevents_tol))]
mincountsbvn_tol = [int(np.min(Transpose_bvn_tol[i])) for i in range (len(bvnevents_tol))]

span_tolbvn = np.array(maxcountsbvn_tol) - np.array(mincountsbvn_tol)
err_tolbvn = span_tolbvn/2
center_tolbvn = maxcountsbvn_tol-err_tolbvn



width = 0.25   # narrowed from 0.25 to fit four bars per group
print(timeax[bvn0events[:,0]])  # start times
print(timeax[bvn0events[:,1]])  # stop times
print(sum(bvnevents_tol))
print(sum(bvnevents_fast_tol))
print(sum(bvnevents_slow_tol))
print(sum(vnevents_tol))        # NEW
print(sum(vnevents_fast_tol))   # NEW
print(sum(vnevents_slow_tol))   # NEW
# width = 0.2

# plt.bar(x,bvevents_tol_1, label=fr'N B V Structures at {threshold_1 *100:.0f}%', width=width, edgecolor='k', linestyle='--', fill=False)
# plt.bar(x-width, bcounts_1, width=width, edgecolor='k', linestyle='--', fill=False)
# plt.bar(x+width, vcounts_1, width=width, edgecolor='k', linestyle='--', fill=False)



plt.bar(x-width,bvevents_tol_1, label=fr'N B V Structures at {threshold_1 *100:.0f}%', width=width, edgecolor='k', linestyle='--', fill=False)
plt.bar(x, bvnevents_tol_1, width=width, edgecolor='k', linestyle='--', fill=False)
plt.bar(x+width, bnevents_tol_1, width=width, edgecolor='k', linestyle='--', fill=False)


plt.bar(x-width,bvevents_tol,label = fr'BV ',width=width,color='blue')
plt.bar(x,bvnevents_tol,label = 'BVN',width=width,color='k')
plt.bar(x+width,bnevents_tol,label = 'BN',width=width,color='r')
plt.ylabel('Counts')
plt.errorbar(x+width,center_tol,yerr=err_tol,fmt='none',ecolor='k',elinewidth=2, capsize=4)
plt.errorbar(x-width,center_tolbv,yerr=err_tolbv,fmt='none',ecolor='k',elinewidth=2, capsize=4)
plt.errorbar(x,center_tolbvn,yerr=err_tolbvn,fmt='none',ecolor='k',elinewidth=2, capsize=4)
plt.ylim(1e-1,1e4)
plt.yscale('log')
plt.xticks(x,labels)
plt.title('Multi-Variable Structure Analysis '+str(int(threshold*100))+'% threshold ')
plt.legend()
plt.yscale('log')
#%%

# Plots with vn as well
width=0.2
plt.bar(x-1.5*width,bvevents_tol_1, label=fr'N B V Structures at {threshold_1 *100:.0f}%', width=width, edgecolor='k', linestyle='--', fill=False)
plt.bar(x-0.5*width, bvnevents_tol_1, width=width, edgecolor='k', linestyle='--', fill=False)
plt.bar(x+0.5*width, bnevents_tol_1, width=width, edgecolor='k', linestyle='--', fill=False)
plt.bar(x+1.5*width, vnevents_tol_1, width=width, edgecolor='k', linestyle='--', fill=False)   # NEW


plt.bar(x-1.5*width,bvevents_tol,label = fr'BV ',width=width,color='blue')
plt.bar(x-0.5*width,bvnevents_tol,label = 'BVN',width=width,color='k')
plt.bar(x+0.5*width,bnevents_tol,label = 'BN',width=width,color='r')
plt.bar(x+1.5*width,vnevents_tol,label = 'VN',width=width,color='g')                           # NEW
plt.ylabel('Counts')
plt.errorbar(x+0.5*width,center_tol,yerr=err_tol,fmt='none',ecolor='k',elinewidth=2, capsize=4)
plt.errorbar(x-1.5*width,center_tolbv,yerr=err_tolbv,fmt='none',ecolor='k',elinewidth=2, capsize=4)
plt.errorbar(x-0.5*width,center_tolbvn,yerr=err_tolbvn,fmt='none',ecolor='k',elinewidth=2, capsize=4)
plt.errorbar(x+1.5*width,center_tolvn,yerr=err_tolvn,fmt='none',ecolor='k',elinewidth=2, capsize=4)   # NEW
plt.ylim(1e-1,1e4)
plt.xticks(x,labels)
plt.title('Multi-Variable Structure Analysis '+str(int(threshold*100))+'% threshold ')
plt.legend()
plt.yscale('log')
# %%
# Getting the full theta variables (defined at every element of the timerange)
theta_b = get_angles_t(t_windows,Bvecs)
thetab0,thetab1,thetab2,thetab3 = theta_b[0],theta_b[1],theta_b[2],theta_b[3]
theta_v = get_angles_t(t_windows,vivecs)
thetav0,thetav1,thetav2,thetav3 = theta_v[0],theta_v[1],theta_v[2],theta_v[3]

# Function to return a list of the mean angle for each event in a given event list
def event_angles(theta,eventlist):
    angles = []  # define empty angle list  
    # For each event, find mean angle & add to angle list
    for event in eventlist:     
        angles.append(np.nanmean(theta[event[0]:event[1]]))
    return angles


# Magnetic angles (bangles) for all db events
db0_bangles = event_angles(thetab0,db0events)
db1_bangles = event_angles(thetab1,db1events)
db2_bangles = event_angles(thetab2,db2events)
db3_bangles = event_angles(thetab3,db3events)

# Velocity angles (vangles) for all db events
db0_vangles = event_angles(thetav0,db0events)
db1_vangles = event_angles(thetav1,db1events)
db2_vangles = event_angles(thetav2,db2events)
db3_vangles = event_angles(thetav3,db3events)

# Magnetic angles for all dv events
dv0_bangles = event_angles(thetab0, dv0events)
dv1_bangles = event_angles(thetab1, dv1events)
dv2_bangles = event_angles(thetab2, dv2events)
dv3_bangles = event_angles(thetab3, dv3events)

# Velocity angles for all dv events
dv0_vangles = event_angles(thetav0, dv0events)
dv1_vangles = event_angles(thetav1, dv1events)
dv2_vangles = event_angles(thetav2, dv2events)
dv3_vangles = event_angles(thetav3, dv3events)
# Plotting
scale_labels = ['4-10', '10-100', '100-1k', '1k-10k']
nbins = 10
#%%
nbins = 30
##θ_B during δV events 
plt.hist(dv0_bangles,bins=nbins,histtype='step',label=rf'$ \theta_B \ (4-10 \ R_E)$')
plt.hist(dv1_bangles,bins=nbins,histtype='step',label = rf'$ \theta_B \ (10-100 \ R_E)$')
plt.hist(dv2_bangles,bins=nbins,histtype='step',label=rf'$ \theta_B \ (100-1k \ R_E)$')
plt.hist(dv3_bangles,bins=nbins,histtype='step',label=rf'$ \theta_B \ (1k-10k \ R_E)$')
plt.title(rf'$\theta_B$ for all events with $\delta V$ > {100*threshold:.0f}%')
plt.xlabel(r'$\theta_B \ (deg)$'); plt.ylabel('Counts')
plt.yscale('log')
plt.legend()
plt.figure()
#%%
##θ_V during (averaged over) δB events
plt.hist(db0_vangles,bins=nbins,histtype='step',label=rf'$ \theta_B \ (4-10 \ R_E)$')
plt.hist(db1_vangles,bins=nbins,histtype='step',label = rf'$ \theta_B \ (10-100 \ R_E)$')
plt.hist(db2_vangles,bins=nbins,histtype='step',label=rf'$ \theta_B \ (100-1k \ R_E)$')
plt.hist(db3_vangles,bins=nbins,histtype='step',label=rf'$ \theta_B \ (1k-10k \ R_E)$')
plt.title(rf'$\theta_V$ for all events with $\delta B$ >{100*threshold:.0f}%')
plt.xlabel(r'$\theta_B \ (deg)$'); plt.ylabel('Counts')
plt.yscale('log')
plt.legend()

#%%
#θ_B during δB events
nbins=20 
plt.hist(db0_bangles,bins=nbins,histtype='step',label=rf'$ \theta_B \ (4-10 \ R_E)$')
plt.hist(db1_bangles,bins=nbins,histtype='step',label = rf'$ \theta_B \ (10-100 \ R_E)$')
plt.hist(db2_bangles,bins=nbins,histtype='step',label=rf'$ \theta_B \ (100-1k \ R_E)$')
plt.hist(db3_bangles,bins=nbins,histtype='step',label=rf'$ \theta_B \ (1k-10k \ R_E)$')
plt.title(rf'$\theta_B$ for all events with $\delta B$ > {100*threshold:.0f}%')
plt.legend()
plt.xlabel(rf'$\theta_V \ (deg)$')
plt.ylabel(rf'Counts')
plt.yscale('log')
# %%
#θ_V during δV events
plt.hist(dv0_vangles,bins=nbins,histtype='step',label=rf'$ \theta_v \ (4-10 \ R_E)$')
plt.hist(dv1_vangles,bins=nbins,histtype='step',label = rf'$ \theta_v \ (10-100 \ R_E)$')
plt.hist(dv2_vangles,bins=nbins,histtype='step',label=rf'$ \theta_v \ (100-1k \ R_E)$')
plt.hist(dv3_vangles,bins=nbins,histtype='step',label=rf'$ \theta_v \ (1k-10k \ R_E)$')
plt.title(rf'$\theta_V$ for all events with $\delta V$ > {100*threshold:.0f}%')
plt.legend()
plt.xlabel(rf'$\theta_V \ (deg)$')
plt.ylabel(rf'Counts')
plt.yscale('log')

# %%
