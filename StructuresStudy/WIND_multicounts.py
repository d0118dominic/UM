
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
swe_vars = pyspedas.projects.wind.swe(trange=trange,no_update=True)
# denstiy 
mfi_vars = pyspedas.projects.wind.mfi(trange=trange,no_update=True)

# 2. Load SYM-H (High-resolution Dst index from Kyoto)
# This loads 1-minute resolution sym_h data into pytplot
# pyspedas.kyoto.dst(trange=trange, datatype='symh')
# dst(trange=trange, datatype='symh')
# Define variable names for convenience
B_name='BGSE'
ni_name='N_elec'
vi_name = 'U_eGSE'
pos_name = 'DIST'

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


Bvecs = 1e-9*reform(get_data('B')) # Converting from nanotesla to Tesla
vivecs = 1e3*reform(get_data('vi')) #converting km/s to m/s
ni = 1e6*reform(get_data('ni')) # Converting from cm^-3 to m^-3

#%%
# Get B and v magnitudes and Alfven speed
Bmag = np.zeros(len(timeax)) # this just creates the array of the correct length
vmag = np.zeros(len(timeax)) 
va = np.zeros(len(timeax))

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
windspeed_High = 1e-3*np.percentile(vmag,90)
windspeed_Low = 1e-3*np.percentile(vmag,10)

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
	if windtype=='fast':
		b = smooth_var(t_windows_high,Bmag)
		v = smooth_var(t_windows_high,vmag)
		n = smooth_var(t_windows_high,ni)
	if windtype=='slow':
		b = smooth_var(t_windows_low,Bmag)
		v = smooth_var(t_windows_low,vmag)
		n = smooth_var(t_windows_low,ni)
	# Returns the normalized delta (there's also a ddt function which is a time derivative) 
	# These arrays 
	db = delta_var(b)[:-1]
	dv = delta_var(v)[:-1]
	dn = delta_var(n)[:-1]

	# Sets all values above a threshold = 1, zero everywhere else
	# Just flags what is and isn't a structure'
	dbflags,dvflags,dnflags = np.zeros_like(db),np.zeros_like(dv),np.zeros_like(dn)
	print('For a threshold of '+str(int(threshold*100))+'%')
	for i in range(len(db)):
		dbflags[i] = make_binary(db[i],threshold)
		dvflags[i] = make_binary(dv[i],threshold)
		dnflags[i] = make_binary(dn[i],threshold)

	# Make an array that looks for where the flagged structures begin or end
	# If the gradient>0, then it is counted
	# This should only catch the edges of structures 
	db_edges,bcounts = get_struc_edges(dbflags,label='db')
	print('B Structures:')
	print(bcounts/2)
	dv_edges,vcounts = get_struc_edges(dvflags,label='dv')
	print('\nV Structures:')
	print(vcounts/2)
	dn_edges,ncounts = get_struc_edges(dnflags,label ='dn')
	print('\nN Structures:')
	print(ncounts/2)
	return threshold,b,v,n,db,dv,dn,dbflags,dvflags,dnflags,bcounts/2,vcounts/2,ncounts/2,2*db_edges,2*dv_edges,2*dn_edges


threshold=0.1
threshold,b_fast,v_fast,n_fast,db_fast,dv_fast,dn_fast,dbflags_fast,dvflags_fast,dnflags_fast,bcounts_fast,vcounts_fast,ncounts_fast,dbedges_fast,dvedges_fast,dnedges_fast = getvars(threshold=threshold,windtype='fast')
threshold,b_slow,v_slow,n_slow,db_slow,dv_slow,dn_slow,dbflags_slow,dvflags_slow,dnflags_slow,bcounts_slow,vcounts_slow,ncounts_slow,dbedges_slow,dvedges_slow,dnedges_slow = getvars(threshold=threshold,windtype='slow')
threshold,b,v,n,db,dv,dn,dbflags,dvflags,dnflags,bcounts,vcounts,ncounts,dbedges,dvedges,dnedges= getvars(threshold=threshold,windtype='mean')
#%%
# Get a list of all 'events', including the start,stop times of each one
def get_startstops(flagvar):
	flagvar[0],flagvar[-1]=0,0 # Set first & last element to zero to handle diff issue
	diff = np.diff(flagvar) # sets 1 at first element of flag, -1 at last element

	# Get an array of events [start1,stop1], [start2,stop2], etc.
	start, stop = np.where(diff==1)[0]+1, np.where(diff==-1)[0]
	eventarr = np.array([start,stop]).T
	return eventarr

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
	
	# Include other conditions if windtype is fast or slow

	return dbevents,dvevents,dnevents

dbevents,dvevents,dnevents = get_events(windtype='mean')
# dbevents_fast,dvevents_fast,dnevents_fast = get_events(windtype='fast')
# dbevents_slow,dvevents_slow,dnevents_slow = get_events(windtype='slow')


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

# Go through all events, check if it satisfies a multi-condition
# If it does, add it to appropriate list, return lengths of lists
def multistructs(scale=0):
	bvstructs = []
	bnstructs = []
	vnstructs=[]
	bvnstructs = []
	events = dbevents[scale]
	for i in range(len(events)):
		condition_bn = np.all(dnflags[scale][events[i,0]:events[i,1]+1]==1)
		condition_bv = np.all(dvflags[scale][events[i,0]:events[i,1]+1]==1)
		if  condition_bv == True:
			bvstructs.append(events[i])
			print('BV: ' + str(condition_bv))
		if condition_bn == True:
			bnstructs.append(events[i])
			print('BN: ' + str(condition_bn))
		if ((condition_bv==True)&(condition_bn==True)):
			bvnstructs.append(events[i])
			print('Both!!')

	events = dvevents[scale]
	for i in range(len(events)):
		condition_vn = np.all(dnflags[scale][events[i,0]:events[i,1]+1]==1)
		condition_vb = np.all(dbflags[scale][events[i,0]:events[i,1]+1]==1)
		if  condition_vn == True:
			vnstructs.append(events[i])
			print('VN: ' + str(condition_vn))
		if condition_vb == True:
			bvstructs.append(events[i])
			print('BV: ' + str(condition_vb))
		if ((condition_vn==True)&(condition_vb==True)):
			bvnstructs.append(events[i])
			print('Both!!')

	events = dnevents[scale]
	for i in range(len(events)):
		condition_nv = np.all(dvflags[scale][events[i,0]:events[i,1]+1]==1)
		condition_nb = np.all(dbflags[scale][events[i,0]:events[i,1]+1]==1)
		if  condition_nv == True:
			vnstructs.append(events[i])
			print('VN: ' + str(condition_nv))
		if condition_nb == True:
			bnstructs.append(events[i])
			print('BN: ' + str(condition_nb))
		if ((condition_nv==True)&(condition_nb==True)):
			bvnstructs.append(events[i])
			print('Both!!')
	return len(bvstructs),len(bnstructs),len(vnstructs),len(bvnstructs)
# 	# return bvstructs,bnstructs,vnstructs,bvnstructs

# Get the multi-condition structure counts for each scale
bv0counts,bn0counts,vn0counts,bvn0counts = multistructs(0)
bv1counts,bn1counts,vn1counts,bvn1counts = multistructs(1)
bv2counts,bn2counts,vn2counts,bvn2counts = multistructs(2)
bv3counts,bn3counts,vn3counts,bvn3counts = multistructs(3)

# Collect into convenient count lists for plotting 
# bvcounts = [bv0counts,bv1counts,bv2counts,bv3counts]
# bncounts = [bn0counts,bn1counts,bn2counts,bn3counts]
# vncounts = [vn0counts,vn1counts,vn2counts,vn3counts]
# bvncounts = [bvn0counts,bvn1counts,bvn2counts,bvn3counts]
#%%

# Plot of single-variable strucures
print('B Structures: \n' + str(bcounts) + '\n')
print('V Structures: \n' + str(vcounts)+'\n')
print('N Structures: \n' + str(ncounts)+'\n')
# x = np.arange(len(bcounts))
labels = np.array(make_labels(d_windows)) if 'd_windows' in locals() else [f"W{i+1}" for i in range(len(d_windows))]
width = 0.25
plt.bar(x-width,bcounts,label = 'B Structures',width=width,color='b')
plt.bar(x,ncounts,label = 'N Structures',width=width,color='g')
plt.bar(x+width,vcounts,label = 'V Structures',width=width,color='red')
plt.xticks(x,labels)
plt.ylabel('Counts')
plt.ylim(1e0,1e4)
plt.title('Single-Variable Structures at '+str(int(threshold*100))+'% threshold')
plt.legend()
plt.yscale('log')


#%%

# Plot of multi-variable structures
print('BV Structures: \n' + str(bvcounts) + '\n')
print('BN Structures: \n' + str(bncounts)+'\n')
print('VN Structures: \n' + str(vncounts)+'\n')
print('BVN Structures: \n' + str(bvncounts)+'\n')


width = 0.25
plt.bar(x-width,bvncounts,label = 'BVN Structures',width=width,color='k')
plt.bar(x-width/2,bvcounts,label = 'BV Structures',width=width,color='purple')
plt.bar(x,bncounts,label = 'BN Structures',width=width,color = 'teal')
plt.bar(x+width/2,vncounts,label = 'VN Structures',width=width,color='orange')
plt.ylabel('Counts')
plt.ylim(1e0,1e4)
plt.yscale('log')
plt.xticks(x,labels)
plt.title('Multi-Variable Structures at '+str(int(threshold*100))+'% threshold')
plt.legend()
plt.yscale('log')
# %%

# %%
import numpy as np
from scipy.spatial import cKDTree

# Test Arrays
# a = np.array([[1,2],[8,9],[10,12],[15,18],[19,26],[80,90]])
# b = np.array([[0,3],[10,12],[15,18],[18,27]])

# Honestly, I don't understand all the details of how it works
# But it does and its quick
# It has a built-in tolerance  
# Worth giving the function a different name if we decide this is the method to go with 

def func(arr1, arr2, tol=0):
    tree = cKDTree(arr1)
    # nearest neighbor in arr1 for every row of arr2, Chebyshev distance (p=inf)
    # means BOTH components must be within tol
    dist, _ = tree.query(arr2, k=1, p=np.inf,
                         distance_upper_bound=np.nextafter(tol, np.inf))
    return arr2[dist <= tol]


# matches_notolerance = func(a,b,tol=0)
# matches_smalltolerance = func(a,b,tol=1)
# print('Exact Matches:\n' + str(matches_notolerance))
# print('\n')
# print('Close Matches:\n' + str(matches_smalltolerance))

#%%

# Tolerance for each is based on the length of avg window
# Since t_windows is in minutes, need to multiply by minuteint
t = [minuteint*i for i in t_windows]
test0 = func(db0events,dn0events,tol=int(t[0]/2))
test1 = func(db1events,dn1events,tol=int(t[1]/2))
test2 = func(db2events,dn2events,tol=int(t[2]/2))
test3 = func(db3events,dn3events,tol=int(t[3]/2))



withtolerance = [len(test0),len(test1),len(test2),len(test3)]

print('Event counts (exact match): ' + str(notolerance))
print('\nEvent counts (scale-dependent tolerance): ' + str(withtolerance))


width = 0.25
plt.bar(x,notolerance,label = 'Exact Match',width=width,color='k')
plt.bar(x+width,withtolerance,label = 'Scale-dependent Tolerance',width=width,color='r')
# plt.bar(x+2*width,bncounts,label = 'Previous Method',width=width,color='b')
plt.ylabel('Counts')
plt.ylim(1e-1,1e4)
plt.yscale('log')
plt.xticks(x,labels)
plt.title('BN Structures at '+str(int(threshold*100))+'% threshold (event matching method)')
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


# Plotting
nbins=20 
plt.hist(db0_bangles,bins=nbins,histtype='step',label=rf'$ \theta_B \ (4-10 \ R_E)$')
plt.hist(db1_bangles,bins=nbins,histtype='step',label = rf'$ \theta_B \ (10-100 \ R_E)$')
plt.hist(db2_bangles,bins=nbins,histtype='step',label=rf'$ \theta_B \ (100-1k \ R_E)$')
plt.hist(db3_bangles,bins=nbins,histtype='step',label=rf'$ \theta_B \ (1k-10k \ R_E)$')
plt.title(rf'$\theta_B$ for all events with $\delta B$ > {100*threshold:.0f}%')
plt.legend()
plt.xlabel(rf'$\theta_B \ (deg)$')
plt.ylabel(rf'Counts')
plt.yscale('log')
# %%
