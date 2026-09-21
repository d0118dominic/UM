
#%%
import pyspedas
import numpy as np
import matplotlib.pyplot as plt
import sys
import math 
from pyspedas.projects.kyoto import dst

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
def get_struc_edges(var,label = '??'):
	newvar,counts = np.zeros_like(var),np.zeros(len(var))
	# print('\n Number of structure edges for '+label+':')
	# print('\n')
	for i in range(len(var)):
		newvar[i] = abs(np.gradient(var[i]))
		counts[i] = len(newvar[i][newvar[i]!=0])
		# print(len(newvar[i][newvar[i]>0]))
	return newvar,counts

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
# d_windows = np.array([4,100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000])
# d_windows = np.array([4,100,200,300,400,500,600,700,800,900,1000])#,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000])
# d_windows = np.array([1000,2000,3000,4000,5000,6000,7000,8000,9000,10000])#,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000])
# d_windows = np.array([10000,20000,30000,40000,50000,60000,70000,80000,90000,100000])#,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000])


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

def multicounts(type='bvn'):
	flags = np.zeros_like(dbflags)
	for j in range(len(flags)):
		if type=='bvn':
			for i in range(len(flags[0])):
				if ((dbflags[j][i]>0)&(dvflags[j][i]>0)&(dnflags[j][i]>0)):
					flags[j][i]=1
				else:
					flags[j][i]=0
			# edges,counts = get_struc_edges(flags)
		if type=='bv':
			for i in range(len(flags[0])):
				if ((dbflags[j][i]>0)&(dvflags[j][i]>0)):
					flags[j][i]=1
				else:
					flags[j][i]=0
			# edges,counts = get_struc_edges(flags)
		if type=='vn':
			for i in range(len(flags[0])):
				if ((dvflags[j][i]>0)&(dnflags[j][i]>0)):
					flags[j][i]=1
				else:
					flags[j][i]=0
			# edges,counts = get_struc_edges(flags)
		if type=='bn':
			for i in range(len(flags[0])):
				if ((dbflags[j][i]>0)&(dnflags[j][i]>0)):
					flags[j][i]=1
				else:
					flags[j][i]=0
	edges,counts = get_struc_edges(flags)
	return counts/2

#%%
threshold=0.1
# threshold,b_fast,v_fast,n_fast,db_fast,dv_fast,dn_fast,dbflags_fast,dvflags_fast,dnflags_fast,bcounts_fast,vcounts_fast,ncounts_fast,dbedges_fast,dvedges_fast,dnedges_fast = getvars(threshold=threshold,windtype='fast')
# threshold,b_slow,v_slow,n_slow,db_slow,dv_slow,dn_slow,dbflags_slow,dvflags_slow,dnflags_slow,bcounts_slow,vcounts_slow,ncounts_slow,dbedges_slow,dvedges_slow,dnedges_slow = getvars(threshold=threshold,windtype='slow')
threshold,b,v,n,db,dv,dn,dbflags,dvflags,dnflags,bcounts,vcounts,ncounts,dbedges,dvedges,dnedges= getvars(threshold=threshold,windtype='mean')
#%%
bvncounts = multicounts(type='bvn')
bvcounts = multicounts(type='bv')
vncounts = multicounts(type='vn')
bncounts = multicounts(type='bn')
#%%
print('BVN Structures')
print(bvncounts)
print('\n BV Structures')
print(bvcounts)
print('\n VN Structures')
print(vncounts)
print('\n BN Structures')
print(bncounts)
#%%
x = np.arange(len(bcounts))
labels = np.array(make_labels(d_windows)) if 'd_windows' in locals() else [f"W{i+1}" for i in range(len(d_windows))]
width = 0.25
plt.bar(x-width,bcounts,label = 'B Structures',width=width,color='b')
plt.bar(x,vcounts,label = 'V Structures',width=width,color='red')
plt.bar(x+width,ncounts,label = 'N Structures',width=width,color='g')
# plt.bar(x+1.5*width,bncounts,label = type.upper()+' Structures',width=width,color='k')
plt.xticks(x,labels)
plt.ylabel('Counts')
plt.ylim(1e0,1e5)
plt.title('Number of structure edges at '+str(int(threshold*100))+'% threshold')
plt.legend()
plt.yscale('log')


#%%
# width = 0.25
plt.bar(x-width,bvncounts,label = 'BVN Structures',width=width,color='k')
plt.bar(x-width/2,bvcounts,label = 'BV Structures',width=width,color='purple')
plt.bar(x,bncounts,label = 'BN Structures',width=width,color = 'teal')
plt.bar(x+width/2,vncounts,label = 'VN Structures',width=width,color='orange')
plt.ylabel('Counts')
plt.ylim(1e-1,1e5)
plt.yscale('log')
# plt.xticks(x,labels)
plt.title('Number of structure edges at '+str(int(threshold*100))+'% threshold')
plt.legend()
plt.yscale('log')
#%%





def hilight_structures(var,dvar,flags,threshold=threshold,scale=0,ylabel='Var'):
	fig,ax = plt.subplots(2,1,figsize=(12,6))
	ax[0].plot(var[scale],label='')
	ax[1].plot(dvar[scale],label='')
	ax[0].fill_between(range(len(timeax)), max(var[scale]), where=(flags[scale]>0),lw=0.2, alpha=0.3,color='r')
	ax[1].fill_between(range(len(timeax)), 1, where=(flags[scale]>0),lw=0.2, alpha=0.3,color='r')
	ax[1].axhline(y=threshold, color='k', linestyle='--', label='Threshold')
	ax[0].set_ylim(0,max(var[scale]))
	ax[1].set_ylim(0,1)
	ax[0].set_ylabel(ylabel)
	ax[1].set_ylabel('$\delta$'+ylabel+' / '+ylabel)
	return

hilight_structures(b,db,dbflags,threshold=threshold,scale=0,ylabel='B')
hilight_structures(v,dv,dvflags,threshold=threshold,scale=0,ylabel='V')
hilight_structures(n,dn,dnflags,threshold=threshold,scale=0,ylabel='N')




#%%










#%%



















# def make_multibinary (th,type = 'bvn'):
# 	var = db[0]
# 	newvar = np.zeros_like(db)
# 	for j in range(len(db)):
# 		if type=='bvn':
# 			for i in range(len(var)):
# 				if ((db[j][i]>=th)&(dv[j][i]>=th)&(dn[j][i]>=th)):
# 					newvar[j][i] = 1
# 				else:
# 					newvar[j][i] = 0
# 		if type=='bv':
# 			for i in range(len(var)):
# 				if ((db[j][i]>=th)&(dv[j][i]>=th)):
# 					newvar[j][i] = 1
# 				else:
# 					newvar[j][i] = 0
# 		if type=='vn':
# 			for i in range(len(var)):
# 				if ((dv[j][i]>=th)&(dn[j][i]>=th)):
# 					newvar[j][i] = 1
# 				else:
# 					newvar[j][i] = 0
# 		if type=='bn':
# 			for i in range(len(var)):
# 				if ((db[j][i]>=th)&(dn[j][i]>=th)):
# 					newvar[j][i] = 1
# 				else:
# 					newvar[j][i] = 0
# 	return newvar








#%%


# # Get the Delta arrays with the t_windows vals
dBmag_mean = get_scaldeltas(t_windows_mean,Bmag)
dni_mean = get_scaldeltas(t_windows_mean, ni)
dvmag_mean = get_scaldeltas(t_windows_mean , vmag)
dBvecs_mean = get_vecdeltas(t_windows_mean, Bvecs)
dvvecs_mean = get_vecdeltas(t_windows_mean, vivecs)
dB_angles_mean = get_angles_t(t_windows_mean, Bvecs)
dv_angles_mean = get_angles_t(t_windows_mean, vivecs) 


dBmag_high = get_scaldeltas(t_windows_high,Bmag)
dni_high = get_scaldeltas(t_windows_high, ni)
dvmag_high = get_scaldeltas(t_windows_high , vmag)
dBvecs_high = get_vecdeltas(t_windows_high, Bvecs)
dvvecs_high = get_vecdeltas(t_windows_high, vivecs)
dB_angles_high = get_angles_t(t_windows_high, Bvecs)
dv_angles_high = get_angles_t(t_windows_high, vivecs) 


dBmag_low = get_scaldeltas(t_windows_low,Bmag)
dni_low = get_scaldeltas(t_windows_low, ni)
dvmag_low = get_scaldeltas(t_windows_low , vmag)
dBvecs_low = get_vecdeltas(t_windows_low, Bvecs)
dvvecs_low = get_vecdeltas(t_windows_low, vivecs)
dB_angles_low = get_angles_t(t_windows_low, Bvecs)
dv_angles_low = get_angles_t(t_windows_low, vivecs) 



#%%



# Function and examples to get count arrays for multiple/correlated conditions
# As written, it only works with the delta_mean arrays.
# Function could easily be altered to filter based on delta_high or delta_low arrays too
# For that, I suggest making alternate functions to do it



# Set a percent threshold 
threshold = 0.5

# New function to filter a variable based on one to three thresholds on (B,V and/or N)
# If not specified, defaults to no threshold (technically threshold=0)
def get_counts_multi(data_arr,Bthreshold=0,Vthreshold=0,Nthreshold=0):
	counts_list = []
	#threshold = 0.05
	for i in range(len(d_windows) - 1):
		# Set masks for B,V,and N variables
		mask_B = (dBmag_mean[i]>=Bthreshold)
		mask_V = (dvmag_mean[i]>=Vthreshold)
		mask_N = (dni_mean[i]>=Nthreshold)
		count = len(data_arr[i][mask_B & mask_V & mask_N])
		counts_list.append(count)
	return counts_list

# This is what the highcounts version would look like 
def get_highcounts_multi(data_arr,Bthreshold=0,Vthreshold=0,Nthreshold=0):
	counts_list = []
	#threshold = 0.05
	for i in range(len(t_windows) - 1):
		# Set masks for B,V,and N variables
		mask_B = (dBmag_high[i]>=Bthreshold)
		mask_V = (dvmag_high[i]>=Vthreshold)
		mask_N = (dni_high[i]>=Nthreshold)
		count = len(data_arr[i][mask_B & mask_V & mask_N])
		counts_list.append(count)
	return counts_list

def get_lowcounts_multi(data_arr,Bthreshold=0,Vthreshold=0,Nthreshold=0):
	counts_list = []
	#threshold = 0.05
	for i in range(len(t_windows) - 1):
		# Set masks for B,V,and N variables
		mask_B = (dBmag_low[i]>=Bthreshold)
		mask_V = (dvmag_low[i]>=Vthreshold)
		mask_N = (dni_low[i]>=Nthreshold)
		count = len(data_arr[i][mask_B & mask_V & mask_N])
		counts_list.append(count)
	return counts_list
#### Commented out examples ######
# Example: Get list of counts in each scale range where dB/B AND dV/v >= threshold 
# It actually doesn't matter which data_arr_mean you use here
# The multi-condition filter would get the same results, you can check yourself if you want
# counts_db = get_counts_multi(dBmag_mean,Bthreshold=threshold,Vthreshold=threshold)
# counts_dv = get_counts_multi(dvmag_mean,Bthreshold=threshold,Vthreshold=threshold)

# Example: Only filter for dB/B>= threshold
# this would be equivalent to get_counts(dBmag_mean,threshold)
# counts_db = get_counts_multi(dBmag_mean,Bthreshold=threshold)

# Similar Example: Only Filter for dV/V >= threshold
# remember, it doesn't matter which mean variable is in the first argument here
# what matters is the filter condition(s) onlyx
# counts_dv = get_counts_multi(dBmag_mean,Vthreshold=threshold)


# Example: Customize multiple different conditions
BTH,VTH,NTH = 0.05,0.05,0.05
counts = get_counts_multi(dBmag_mean,Bthreshold=BTH,Vthreshold = VTH,Nthreshold=NTH)
counts_high = get_highcounts_multi(dBmag_high,Bthreshold=BTH,Vthreshold = VTH,Nthreshold=NTH)
counts_low = get_lowcounts_multi(dBmag_low,Bthreshold=BTH,Vthreshold = VTH,Nthreshold=NTH)

labels = np.array(make_labels(d_windows)) if 'd_windows' in locals() else [f"W{i+1}" for i in range(len(d_windows))]

plt.bar(x,counts)
plt.xticks(x, labels)
plt.xlabel(fr"Scale Size ($R_E$)", fontsize=15)
plt.yscale('log')
plt.title(fr'$\delta B \geq$ {BTH*100}% | $\delta V \geq$ {VTH*100}% | $\delta N \geq$ {NTH*100}%',fontsize=15)
plt.ylabel('Structure Count',fontsize=15)

#%%


#%%

BTH,VTH,NTH = 0.05,0.00,0.00

nullmeancounts = np.array(get_counts_multi(dBmag_mean,0,0,0))
nullfastcounts = np.array(get_counts_multi(dBmag_high,0,0,0))
nullslowcounts = np.array(get_counts_multi(dBmag_low,0,0,0))

meancounts = np.array(get_counts_multi(dBmag_mean,BTH, VTH, NTH))
slowcounts = np.array(get_lowcounts_multi(dBmag_low,BTH, VTH, NTH))
fastcounts = np.array(get_highcounts_multi(dBmag_high,BTH, VTH, NTH))

Transpose_arr = np.array([fastcounts,slowcounts]).T
maxcounts = np.array([np.nanmax(Transpose_arr[i]) for i in range(len(meancounts))])
mincounts = np.array([np.nanmin(Transpose_arr[i]) for i in range(len(meancounts))])
span = maxcounts-mincounts
err = span/2
center = maxcounts-err

print('Low Speed Counts:',counts_low)
print('High Speed Counts:',counts_high)
print('Counts:',counts)

# error_label = f"Error Bars Mean: (Low: -{avg_low_err:.2f}, Up: +{avg_up_err:.2f})"


labels = np.array(make_labels(d_windows)) if 'd_windows' in locals() else [f"W{i+1}" for i in range(len(d_windows))]
width = 0.35
plt.figure(figsize=(14,10))
# plt.bar(x,maxcounts , label="Max Counts")
# plt.bar(x,100*(meancounts/nullmeancounts) , label="Percentage of Structures")
plt.bar(x,meancounts)
plt.errorbar(x ,center, yerr=err, fmt='none', ecolor='navy',
              elinewidth=2, capsize=4)

# plt.bar(x,mincounts , label="Min Counts")
plt.xticks(x, labels)
plt.xlabel(fr"Scale Size ($R_E$)", fontsize=15)
plt.yscale('log')
plt.legend()
plt.title(fr'$\delta B \geq$ {BTH*100}% | $\delta V \geq$ {VTH*100}% | $\delta N \geq$ {NTH*100}%',fontsize=15)
plt.ylabel('Structure Count',fontsize=15)
# plt.ylim((1000,30000))















#%%
TH=0.05
#LOW
counts_db_low      = get_counts(dBmag_low, threshold)
counts_dv_low      = get_counts(dvmag_low, threshold)
counts_dni_low    = get_counts(dni_low, threshold)
counts_dBvecs_low  = get_counts(dBvecs_low, threshold)
counts_dvvecs_low  = get_counts(dvvecs_low, threshold)
print("counts_db_low: " , counts_db_low)
print("counts_dv_low: " , counts_dv_low)
print("counts_dni_low: " , counts_dni_low)
print("counts_dbvec_low: " , counts_dBvecs_low)
print("counts_dvvec_low: " , counts_dvvecs_low)

#HIGH

counts_db_high      = get_counts(dBmag_high, threshold)
counts_dv_high      = get_counts(dvmag_high, threshold)
counts_dni_high     = get_counts(dni_high, threshold)
counts_dBvecs_high  = get_counts(dBvecs_high, threshold)
counts_dvvecs_high  = get_counts(dvvecs_high, threshold)

print("counts_db_high: " , counts_db_high,'\n','\n','\n')
print("counts_dv_high: " , counts_dv_high)
print("counts_dni_high: " , counts_dni_high)
print("counts_dbvec_high:" , counts_dBvecs_high)
print("counts_dvvec_high: " , counts_dvvecs_high,'\n','\n','\n')
#MEAN

counts_db      = get_counts(dBmag_mean, threshold)
counts_dv      = get_counts(dvmag_mean, threshold)
counts_dni     = get_counts(dni_mean, threshold)
counts_dBvecs  = get_counts(dBvecs_mean, threshold)
counts_dvvecs  = get_counts(dvvecs_mean, threshold)
print("counts_db_mean: " , counts_db)
print("counts_dv_mean: " , counts_dv)
print("counts_dni_mean :" , counts_dni)
print("counts_dbvec_mean: " , counts_dBvecs)
print("counts_dvvec_mean :" , counts_dvvecs,'\n','\n','\n')


target_edges = d_windows 
x_midpoints = (d_windows[:-1] + d_windows[1:]) / 2

labels = labels = make_labels(d_windows) if 'd_windows' in locals() else [f"W{i+1}" for i in range(10)]

# x = np.arange(10)
#vmag 


high_vmag =  np.array(counts_dv_high)
low_vmag =   np.array(counts_dv_low)
mean_vmag =   np.array(counts_dv)

lower_error_vmag = np.abs(mean_vmag - low_vmag)
upper_error_vmag = np.abs(high_vmag - mean_vmag)
asy_error_vmag =  [lower_error_vmag,upper_error_vmag]



#dni
high_dni =  np.array(counts_dni_high)
low_dni =  np.array(counts_dni_low)
mean_dni =    np.array(counts_dni)

lower_error_dni = np.abs(mean_dni - low_dni)
upper_error_dni =  np.abs(high_dni - mean_dni)
asy_error_dni =  [lower_error_dni,upper_error_dni]


#just for db 
high_db =  np.array(counts_db_high)
low_db =  np.array(counts_db_low)
mean_db =  np.array(counts_db)


lower_error = np.abs(mean_db - low_db)
upper_error =  np.abs(high_db - mean_db)
asy_error_db = [lower_error,upper_error]

width = 0.30
plt.figure(figsize=(14, 10))

# DB 
plt.bar(x - width, mean_db, width=width, color="royalblue", label="sw: mean_db")
plt.errorbar(x - width, mean_db, yerr=asy_error_db, fmt='none', ecolor='navy',
			  elinewidth=2, capsize=4,label="Error Bars")
# DNI 
plt.bar(x, mean_dni, width=width, color="crimson", label="sw: mean_dni")
plt.errorbar(x, mean_dni, yerr=asy_error_dni, fmt='none', 
			 ecolor='darkred', elinewidth=2, capsize=4)
# VMAG 
plt.bar(x + width, mean_vmag, width=width, color="pink", label="sw: mean_vmag")
plt.errorbar(x + width, mean_vmag, yerr=asy_error_vmag, fmt='none',
			  ecolor='hotpink', elinewidth=2, capsize=4)
plt.xticks(x, labels,fontsize=13)
plt.xlabel(fr"Scale Size ($R_E$) Threshold: {TH * 100}", fontsize=25)
plt.ylabel("Counts",fontsize=25)
plt.legend(fontsize=25)
plt.tight_layout()
plt.yscale('log')
plt.yticks(fontsize=20)
plt.xticks(fontsize=14)
plt.ylim(1,1e5)
plt.show()
#%%


angle_threshold = 5 # degrees
mask0 = (dvmag[0] >= threshold) & (dBmag[0]>=threshold) & (dB_angles[0]>=angle_threshold)

# Plot Fluctuations
plt.plot(dBmag[0], label=r'$\delta B/B$ (Magnetic fluctuation)')
plt.plot(dvmag[0], label=r'$\delta v/v$ (Velocity fluctuation)')
plt.xlabel('Time Index')
plt.title(fr'Magnetic and Velocity Fluctuations ($\theta > {angle_threshold}^\circ$)')
plt.axhline(y=threshold, color='r', linestyle='-',label = r'Fluctuation threshold')
plt.legend()
plt.ylim(0,1)

print("Number of points above threshold: ", len(dBmag[0][mask0]))


#%%
d_windows = make_d_windows()

# d_windows = np.array([4,10,100,150,200,250,500,600,700,750,1000])
labels = make_labels(d_windows) 
labels = np.array(labels)
x = np.arange(len(labels))


threshold = 0.05
threshold_correlation = threshold
#get_counts_mask()
Bscal = get_scaldeltas(t_windows_test, Bmag)
dni = get_nideltas(t_windows_test, ni)
Velocity_mag = get_scaldeltas(t_windows_test, vmag)




#%%


Muiti_Threshold = 0.05

print(len(dBmag_mean[0][(dBmag_mean[0]>=0.05)&()]))

#%%


counts_db_mask = get_counts(dBmag_mean, threshold)
counts_dv_mask  = get_counts(dvmag_mean, threshold)
counts_dni_mask   = get_counts(dni_mean, threshold)
counts_dBvecs_mask = get_counts(dBvecs_mean, threshold)
counts_dvvecs_mask  = get_counts(dvvecs_mean, threshold)

db_counts   = np.array(counts_db_mask)
dv_counts   = np.array(counts_dv_mask)
dni_counts  = np.array(counts_dni_mask)


bins = np.logspace(np.log10(Muiti_Threshold), np.log10(1.0), 15)


def sum(arr):
	total = 0
	for i in range (len(arr)):
		total = total + arr[i]

	return total


mask_Bmag_Vmag = []
mask_Bmag_dni = []
mask_vmag_Dni = []

for i in range(len(dBmag_mean)):
	mask_b_v = (dBmag_mean[i] >= 
	Muiti_Threshold) & (dvmag_mean[i] >= Muiti_Threshold)


	mask_b_d = (dBmag_mean[i] >= 
	Muiti_Threshold) & (dni_mean[i] >= Muiti_Threshold)

	mask_v_d = (dvmag_mean[i] >= 
	Muiti_Threshold) & (dni_mean[i] >= Muiti_Threshold)


	
	mask_Bmag_Vmag.append(mask_b_v)
	mask_Bmag_dni.append(mask_b_d)
	mask_vmag_Dni.append(mask_v_d)

#%%
#test cell
print(len(mask_Bmag_Vmag))
print(len(mask_vmag_Dni)) 
print(len(mask_Bmag_dni))

j = sum(mask_Bmag_Vmag)
e = sum(mask_vmag_Dni)
r = sum(mask_Bmag_dni)

print(len(j))
print(len(e))
print(len(r))
#%%

# if(mask_Bmag_Vmag.sum() == mask_Bmag_dni.sum()):


# 	print(fr"ERROR : vals {mask_Bmag_Vmag.sum()} AND bdni {mask_Bmag_dni.sum()}")
#%%
# for i in range(len(dvmag_mean)):


# 	mask_b_v_d = (dBmag_mean[i] >= 
# 		Muiti_Threshold) & (dni_mean[i] >= Muiti_Threshold) & (dvmag_mean[i] >= Muiti_Threshold)
# 	mask_Bmag_Vmag_dni.append(mask_Bmag_Vmag_dni)



# n = mask_Bmag_Vmag_dni 


# plt.hist(dBmag_mean[0][n[0]], histtype='step', color='navy', bins=bins, label=fr'4 - 10 $R_E$ (counts {(n[0].sum())})')
# plt.hist(dBmag_mean[1][n[1]], histtype='step', color='red', bins=bins, label=fr'10 - 100 $R_E$ (counts {(n[1].sum())})')
# plt.hist(dBmag_mean[2][n[2]], histtype='step'  ,bins=bins,label=fr'100 - 1000 $R_E$ (counts {(n[2].sum())})')
# plt.hist(dBmag_mean[3][n[3]], histtype='step'  ,bins=bins,label=fr'1k - 10k $R_E$ (counts {(n[3].sum())})')
# plt.ylabel("Test", fontsize='16')
# #\geq$ Latex >=
# plt.xlabel(fr"$\Delta$ B/B AND Desnity $\geq$ {Muiti_Threshold * 100}")
# plt.ylim(1,1e5)
# plt.yscale('log')
# plt.legend()
# plt.show()



#%%
m = mask_Bmag_Vmag
k = mask_Bmag_dni
n = mask_vmag_Dni

percent_bins = 10
bins = percent_bins 


plt.hist(dBmag_mean[0][m[0]], histtype='step', color='navy', bins=bins, label=fr'4 - 10 $R_E$ (counts {(m[0].sum())})')
plt.hist(dBmag_mean[1][m[1]], histtype='step', color='blue', bins=bins, label=fr'10 - 100 $R_E$ (counts {(m[1].sum())})')
plt.hist(dBmag_mean[2][m[2]], histtype='step' , color = 'brown' ,bins=bins,label=fr'100 - 1000 $R_E$ (counts {(m[2].sum())})')
plt.hist(dBmag_mean[3][m[3]], histtype='step' , color = 'green' ,bins=bins,label=fr'1000 - 10000 $R_E$ (counts {(m[3].sum())})')
plt.ylabel("Counts above Threshold ", fontsize='14')
#\geq$ Latex >=
plt.xlabel(fr"$\Delta$ B/B AND Vmag $\geq$ {Muiti_Threshold * 100}")
plt.ylim(1,1e5)
plt.xlim(0.05,1)
plt.yscale('log')
plt.legend()
plt.show()


plt.hist(dvmag_mean[0][n[0]], histtype='step', color='navy', bins=bins, label=fr'4 - 10 $R_E$ (counts {(n[0].sum())})')
plt.hist(dvmag_mean[1][n[1]], histtype='step', color='blue', bins=bins, label=fr'10 - 100 $R_E$ (counts {(n[1].sum())})')
plt.hist(dvmag_mean[2][n[2]], histtype='step' , color = 'brown' ,bins=bins,label=fr'100 - 1000 $R_E$ (counts {(n[2].sum())})')
plt.hist(dvmag_mean[3][n[3]], histtype='step' , color = 'green' ,bins=bins,label=fr'1000 - 10000 $R_E$ (counts {(n[3].sum())})')
plt.ylabel("Counts above Threshold", fontsize='14')
#\geq$ Latex >=
plt.xlabel(fr"$\Delta$ V/V AND DNI $\geq$ {Muiti_Threshold * 100}")
plt.ylim(1,1e5)
plt.xlim(0.05,1)
plt.yscale('log')
plt.legend()
plt.show()



plt.hist(dBmag_mean[0][k[0]], histtype='step', color='navy', bins=bins, label=fr'4 - 10 $R_E$ (counts {(k[0].sum())})')
plt.hist(dBmag_mean[1][k[1]], histtype='step', color='red', bins=bins, label=fr'10 - 100 $R_E$ (counts {(k[1].sum())})')
plt.hist(dBmag_mean[2][k[2]], histtype='step'  ,bins=bins,label=fr'100 - 1000 $R_E$ (counts {(k[2].sum())})')
plt.hist(dBmag_mean[3][k[3]], histtype='step'  ,bins=bins,label=fr'1k - 10k $R_E$ (counts {(k[3].sum())})')
plt.hist(dBmag_mean[4][k[4]], histtype='step'  ,bins=bins,label=fr'1k - 10k $R_E$ (counts {(k[3].sum())})')
plt.ylabel("Counts above Threshold", fontsize='14')
#\geq$ Latex >=
plt.xlabel(fr"$\Delta$ B/B AND Desnity $\geq$ {Muiti_Threshold * 100}")
plt.ylim(1,1e5)
plt.xlim(0.05,1)
plt.yscale('log')
plt.legend()
plt.show()






#%%%

#%%

# ERROR_TEST_B_V = m.sum()
# ERROR_TEST_B_D= k.sum()

# if( ERROR_TEST_B_V == ERROR_TEST_B_D):
# 	print(fr"ERROR Values : {(ERROR_TEST_B_D)} AND {ERROR_TEST_B_V}")

# else:
# 	print(fr"they are not the same , Values : {ERROR_TEST_B_D} AND {ERROR_TEST_B_V}")
	
#%%
#Multiple histograms for multiple scale ranges
# This method should probably altered if there are a lot of scale bins
# Example here uses 4 

# Adjust Thresholds on dB/B and dv/v.  
# If you don't want any threshold on one, set it to zero
Bthreshold = 0.05
vthreshold = 0.05

# Set number of bins (50 seems fine)
bins = 50

# Set your mask to filter the data  
# Right now its set up to exclude all data that doesn't meet both thresholds
mask = [(dBmag_mean[i] >= Bthreshold) & (dvmag_mean[i] >= vthreshold) for i in range(len(dBmag))]

mask1 = [(dBmag_mean[i] >= Bthreshold) & (dni_mean[i] >= vthreshold) for i in range(len(dBmag))]

mask2 = [(dvmag_mean[i] >= Bthreshold) & (dni_mean[i] >= vthreshold) for i in range(len(dBmag))]

# Histograms, not filled in so we can see all of them 
# But change anything around as is convenient
plt.hist(dB_angles[0][mask[0]],bins=bins,histtype='step',color='r', label = '4-10 $R_E$')
plt.hist(dB_angles[1][mask[1]],bins=bins,histtype='step',color = 'g', label = '10-100 $R_E$')
plt.hist(dB_angles[2][mask[2]],bins=bins,histtype='step',color='b', label = '100-1k $R_E$')
plt.hist(dB_angles[3][mask[3]],bins=bins,histtype='step',color='violet', label = '1k-10k $R_E$')
plt.xlabel(r'$\theta_B$ (degrees)' + f' \n Threshold: $\delta B/B \geq$ {Bthreshold*100}% & $\delta v/v \geq$ {vthreshold*100}%')
plt.ylabel('Counts')
plt.yscale('log')
plt.ylim(1,1e4)
plt.legend(loc='upper right')
#%%

plt.hist(dB_angles[0][mask1[0]],bins=bins,histtype='step',color='r', label = '4-10 $R_E$')
plt.hist(dB_angles[1][mask1[1]],bins=bins,histtype='step',color = 'g', label = '10-100 $R_E$')
plt.hist(dB_angles[2][mask1[2]],bins=bins,histtype='step',color='b', label = '100-1k $R_E$')
plt.hist(dB_angles[3][mask1[3]],bins=bins,histtype='step',color='violet', label = '1k-10k $R_E$')
plt.xlabel(r'$\theta_B$ (degrees)' + f' \n Threshold: $\delta B/B \geq$ {Bthreshold*100}% & $\delta dni \geq$ {vthreshold*100}%')
plt.ylabel('Counts')
plt.yscale('log')
plt.ylim(1,1e4)
plt.legend(loc='upper right')

#%%

plt.hist(dv_angles[0][mask2[0]],bins=bins,histtype='step',color='r', label = '4-10 $R_E$')
plt.hist(dv_angles[1][mask2[1]],bins=bins,histtype='step',color = 'g', label = '10-100 $R_E$')
plt.hist(dv_angles[2][mask2[2]],bins=bins,histtype='step',color='b', label = '100-1k $R_E$')
plt.hist(dv_angles[3][mask2[3]],bins=bins,histtype='step',color='violet', label = '1k-10k $R_E$')
plt.xlabel(r'$\theta_B$ (degrees)' + f' \n Threshold: $\delta v/v \geq$ {Bthreshold*100}% & $\delta dni \geq$ {vthreshold*100}%')
plt.ylabel('Counts')
plt.yscale('log')
plt.ylim(1,1e4)
plt.xlim(0,25)
plt.legend(loc='upper right')

#%%


#%%
Bscal_0 = Bscal[0]

dni_0 = dni[0]

vmag_0 = Velocity_mag[0] 

mask_bmag_dni = (Bscal_0 >= threshold) & (dni_0 >= threshold)
mask_vmag_dni = (vmag_0 >= threshold) & (dni_0 >= threshold) 
mask_vmag_bmag = (vmag_0 >= threshold) & (Bscal_0 >= threshold)

#Keep values where BOTH Bmag AND dni meet threshold

clean_Bscal_A = Bscal_0[mask_bmag_dni] 
clean_dni_A   = dni_0[mask_bmag_dni]

print(clean_Bscal_A)
print(clean_dni_A)

arranged_clean_Bmag_A = np.arange(len(clean_Bscal_A))
arranged_clean_dni_A = np.arange(len(clean_dni_A))

print("Bmag  Length : ", len(arranged_clean_Bmag_A))
#Keep values where BOTH vmag AND dni meet threshold
clean_vmag_B  = vmag_0[mask_vmag_dni]
clean_dni_B   = dni_0[mask_vmag_dni]

arranged_clean_vmag_B = np.arange(len(clean_vmag_B))
arranged_clean_dni_B = np.arange(len(clean_dni_B))

# Keep values where BOTH vmag AND Bmag meet threshold
clean_vmag_C  = vmag_0[mask_vmag_dni]
clean_Bscal_C = Bscal_0[mask_vmag_bmag]

arranged_clean_vmag_C = np.arange(len(clean_vmag_C))
arranged_clean_bscal_C = np.arange(len(clean_Bscal_C ))


x_filtered = np.arange(len(clean_Bscal_A))	#254

#%%
width1 = 4
Move_Bar_Left = 0.05
Move_Bar_Right = 0.05
plt.figure(figsize=(12, 5))
plt.bar(arranged_clean_Bmag_A - Move_Bar_Left , clean_Bscal_A,
width=width1, alpha=0.7, color="blue", label=r"$\Delta B_{scal}$")
# Map your labels evenly across the 254 
tick_positions = np.linspace(0, len(clean_Bscal_A) - 1, len(labels))
plt.xticks(tick_positions, labels, rotation=45, ha='right')
plt.xlabel(fr"Scale Size ($R_E$) ")
plt.ylabel("Counts", fontsize=12)
plt.title(fr" Mag Field >= {threshold_correlation * 100} AND  dni >= {threshold_correlation * 100} ", fontsize=12)
plt.yscale('log')
plt.xlim(-1)
plt.grid(True)
plt.legend(fontsize=10)
# plt.tight_layout()
plt.show()
#%%

plt.figure(figsize=(12, 5))
# catgory , value
plt.bar(arranged_clean_vmag_B - Move_Bar_Left, clean_vmag_B,
width=width1, alpha=0.7, color="blue", label=r"$\Delta V_{mag}$",edgecolor='black')
plt.bar(arranged_clean_dni_B + Move_Bar_Right, clean_dni_B ,width=width1,alpha=0.9, color = "pink", label=r"$\Delta DNI$",edgecolor='black')
# Map your labels evenly across the 254 
tick_positions = np.linspace(0, len(clean_vmag_B) - 1, len(labels))
plt.xticks(tick_positions, labels, rotation=45, ha='right')
plt.xlabel(fr"Scale Size ($R_E$)", fontsize=12)
plt.ylabel("Counts", fontsize=12)
plt.yscale('log')
plt.title(fr"vmag >= {threshold_correlation * 100} AND  dni >= {threshold_correlation * 100} ", fontsize=12)
plt.xlim(-1)
plt.grid(True , alpha = 0.3)
plt.legend(fontsize=10)
plt.tight_layout()
plt.show()

#%%

plt.figure(figsize=(12, 5))

# catgory , value
plt.bar(arranged_clean_vmag_C + Move_Bar_Right, clean_vmag_C,
width=width1, alpha=0.7, color="blue", label=r"$\Delta V_{mag}$",edgecolor='black')
plt.bar(arranged_clean_bscal_C - Move_Bar_Left , clean_Bscal_C , width=width1,alpha=0.7, color = "pink", label=r"$\Delta Magnetic     Field$",edgecolor='black')
# # Map your labels evenly across the 254 
tick_positions = np.linspace(0, len(clean_vmag_B) - 1, len(labels))
plt.xticks(tick_positions, labels, rotation=45, ha='right')
plt.xlabel(fr"Scale Size ($R_E$)", fontsize=12)
plt.ylabel("Counts", fontsize=12)
plt.title(fr" vmag >= {threshold_correlation * 100} AND  dni >= {threshold_correlation * 100} ", fontsize=12)
plt.yscale('log')
plt.xlim(-1)
plt.grid(True , alpha = 0.3)
plt.legend(fontsize=10)
plt.tight_layout()
plt.show()


#%%
# plt.figure(figsize=(14, 10))

# # DB 
# plt.bar(x - width, mean_db, width=width, color="royalblue", label="sw: mean_db")
# plt.errorbar(x - width, mean_db, yerr=asy_error_db, fmt='none', ecolor='navy',
# 			  elinewidth=2, capsize=4,label="Error Bars")

# # DNI 
# plt.bar(x, mean_dni, width=width, color="crimson", label="sw: mean_dni")
# plt.errorbar(x, mean_dni, yerr=asy_error_dni, fmt='none', 
# 			 ecolor='darkred', elinewidth=2, capsize=4)

# # VMAG 
# plt.bar(x + width, mean_vmag, width=width, color="pink", label="sw: mean_vmag")
# plt.errorbar(x + width, mean_vmag, yerr=asy_error_vmag, fmt='none',
# 			  ecolor='hotpink', elinewidth=2, capsize=4)


# plt.xticks(x, labels,fontsize=13)
# plt.xlabel(fr"Scale Size ($R_E$) Threshold: {TH * 100}", fontsize=25)
# plt.ylabel("Counts",fontsize=25)
# plt.legend(fontsize=25)
# plt.tight_layout()
# plt.yscale('log')
# plt.yscale('log')
# plt.xlim()
# plt.ylim(1,100000)
# plt.show()

#%%
width = 0.2 
plt.figure()
plt.bar(clean_vmag_B, clean_dni_B, width=width, alpha=0.5,color="black")
plt.xlabel(fr"    Velocity Mag x-axis : vmag >= {threshold_correlation} AND Density >= {threshold_correlation}")
plt.ylabel('DNI')
plt.title('Vmag vs DNI (Filtered)')
plt.grid(True)
plt.show()
#%%
plt.figure()
plt.scatter(clean_vmag_C, clean_Bscal_C, alpha=0.5,color='black')
plt.xlabel(fr"Velocity Mag x-axis : vmag >= {threshold_correlation * 100} AND Magnetic Field >= {threshold_correlation * 100}")
plt.ylabel('Bscal')
plt.legend()
plt.grid(True)
plt.title('Vmag vs Bscal (Filtered)')
plt.show()
#%%
# d_windows = make_d_windows()
# # d_windows = np.array([4,10,100,150,200,250,500,600,700,750,1000])
# labels = make_labels(d_windows) 
# labels = np.array(labels)
# x = np.arange(len(labels))

# plt.scatter(new_Flitered_vector1, new_Flitered_dni1_v ,s=10,color='blue',
# 			label=fr"Scale 1 N={N0 + N3}")
# plt.scatter(new_Flitered_vector0, new_Flitered_dni0_v ,s=10,color='black',alpha= 0.7,
# 			label=fr"Scale 2 N={N1 + N4}")
# # plt.xticks(x, labels)
# plt.xlabel(fr'$\delta v/v$ Threshold : {threshold * 100}%')
# plt.ylabel(r'$\delta n/n$')
# plt.plot(Bscal_0,Bscal_0,color='r', label= 'Ref line')
# plt.legend() 
# plt.grid(True)
# plt.title("Correlated plot of V/V ~ ni/ni")
# plt.ylim(0,1)
# plt.xlim(0,1)

    #%%
# width = 0.2
# plt.figure(figsize=(9, 5))
# # x - 0.2 shifting bar to the left
# plt.bar(x - width ,counts_db_binned, width=0.3, 
# color="skyblue", edgecolor="black", label="Sw |B| Structures")
# plt.bar(x , counts_dv_binned, width=0.3,color="navy",
# 		edgecolor="black",label="Sw |V| Structures")
# plt.bar(x+width, counts_dni_binned, width=0.3,color="pink",
# 		edgecolor="black",label="Sw ni Structures")
# plt.xticks(x, labels)
# plt.title("Solar Wind Structures Above A Set Threshold")
# plt.xlabel(fr"Scale Size ($R_E$) Choosen Fliter: {threshold * 100}" "%" "\n" 
# 		   "Note: Bars are horizontally offset for visual clarity" )
# # plt.errorbar(x + 0.2,counts_db_binned,yerr=y_error,fmt='o')
# plt.ylabel("Structure Counts")
# plt.yscale('log')
# plt.ylim(1e0,1e5)
# plt.legend(fontsize=14)
# plt.show()

# %%

# Choose an angle threshold for the magnetic and velocity deflections 
# (We can choose different ones)

B_angle_threshold = 20
V_angle_threshold = 20

# Get Angle Arrays
angle_array_db = get_angles_t(t_windows,Bvecs)
angle_array_vi = get_angles_t(t_windows,vivecs)

# Get Counts of each above their defined thresholds
counts_angles_db = get_counts(angle_array_db,B_angle_threshold) 
counts_angles_vi = get_counts(angle_array_vi,V_angle_threshold) 

# Binning & Plotting
counts_angles_binned_db = aggregate_counts_by_bounds(x_midpoints, counts_angles_db,
							 target_edges)

counts_angles_binned_vi = aggregate_counts_by_bounds(x_midpoints, counts_angles_vi,
							 target_edges)


# Create x positions that match the number of bins
# label = np.arange(len(counts_angles_binned_db))

# labels = ['1–10', '10–100', '100–1000']
# labels = labels #Using the same definition from above
# label = np.arange(len(labels))
# print(len(label))
# print(len(counts_angles_binned_db))
# print(len(counts_angles_binned_vi))
# print(len(labels))
# d_windows=np.array([4,10,100,150,200,250,500,600,700,750,1000])

d_windows = make_d_windows()
labels = make_labels(d_windows) 
labels = np.array(labels)
x = np.arange(len(labels))
plt.figure(figsize=(14,5))
plt.bar(x + 0.3 ,counts_angles_binned_db, width=0.3, 
color="skyblue", edgecolor="black", label=fr'$\theta_B > {B_angle_threshold}^o$ Structures')
plt.bar( x + 0.1 ,counts_angles_binned_vi, width=0.3, 
color="pink", edgecolor="black", label=fr'$\theta_V > {V_angle_threshold}^o$ Structures')
# plt.bar(label + 0.2 ,counts_angles_binned_ni, width=0.6, 
# color="skyblue", edgecolor="black", label=fr'$\theta > {angle_threshold}^o$ Structures')
plt.xticks(x, labels, rotation=45, ha='right')
plt.xlabel(r"Scale Size ($R_E$)",fontsize=16)
plt.xticks(x, labels)
plt.ylabel("Counts" ,fontsize=16)
plt.yscale('log')
plt.ylim(1e0,1e6)
plt.legend(fontsize=14)
plt.show()



# # print(len(label))
# print(len(counts_angles_binned_db))
# print(len(counts_angles_binned_vi))

# %%
threshold = 0.40


Bx = Bvecs[:,0]
By = Bvecs[:,1]
Bz = Bvecs[:,2]

vx = vivecs[:,0]
vy = vivecs[:,1]
vz = vivecs[:,2]


dBx = get_scaldeltas(t_windows, Bx)
dBy = get_scaldeltas(t_windows, By)
dBz = get_scaldeltas(t_windows, Bz)

dvx = get_scaldeltas(t_windows, vx)
dvy = get_scaldeltas(t_windows, vy)
dvz = get_scaldeltas(t_windows, vz)


counts_Bx = get_counts(dBx, threshold)
counts_By = get_counts(dBy, threshold)
counts_Bz = get_counts(dBz, threshold)

counts_vx = get_counts(dvx, threshold)
counts_vy = get_counts(dvy, threshold)
counts_vz = get_counts(dvz, threshold)

counts_dbx_binned = aggregate_counts_by_bounds(x_midpoints, counts_Bx, target_edges)
counts_dby_binned = aggregate_counts_by_bounds(x_midpoints, counts_By, target_edges)
counts_dbz_binned = aggregate_counts_by_bounds(x_midpoints, counts_Bz, target_edges)
counts_dvx_binned = aggregate_counts_by_bounds(x_midpoints, counts_vx, target_edges)
counts_dvy_binned = aggregate_counts_by_bounds(x_midpoints, counts_vy, target_edges)
counts_dvz_binned = aggregate_counts_by_bounds(x_midpoints, counts_vz, target_edges)

# d_windows = make_d_windows()
# d_windows=np.array([4,10,100,200,300,400,500,600,700,750,1000])
d_windows = make_d_windows()
labels = make_labels(d_windows) 
labels = np.array(labels)
x = np.arange(len(labels))


#Plotting the dB/B components
plt.figure(figsize=(14, 5))
# x - 0.2 shifting bar to the left
width = 0.16

plt.bar(x - 0.1, counts_dbx_binned, width=width, label=r"Sw $\delta B_x/B_x$")
plt.bar(x + 0.00, counts_dby_binned, width=width, label=r"Sw $\delta B_y/B_y$")
plt.bar(x + 0.1, counts_dbz_binned, width=width, label=r"Sw $\delta B_z/B_z$")

plt.xticks(x, labels)
plt.title("Solar Wind Component Fluctuations")
plt.xlabel(fr"Scale Size ($R_E$) Choosen Fliter: {threshold * 100}" "%" "\n" 
		   "Note: Bars are horizontally offset for visual clarity" )
plt.ylabel("Structure Counts")
plt.yscale('log')
plt.ylim(1e3,1e6)
plt.legend(fontsize= 15 , loc = "upper center" ,bbox_to_anchor=(0.5,1.02),ncol= 6)
plt.show()



#%%
#Plotting the dv/v components

plt.figure(figsize=(8, 4))
# x - 0.2 shifting bar to the left
width = 0.16

plt.bar(x - 0.1, counts_dvx_binned, width=width, label=r"Sw $\delta V_x/V_x$")
plt.bar(x + 0.00, counts_dvy_binned, width=width, label=r"Sw $\delta V_y/V_y$")
plt.bar(x + 0.1, counts_dvz_binned, width=width, label=r"Sw $\delta V_z/V_z$")
plt.xticks(x, labels)
plt.title("$\delta v/v$ Components")
plt.xlabel(fr"Scale Size ($R_E$) Choosen Fliter: {threshold * 100}" "%" "\n" 
		   "Note: Bars are horizontally offset for visual clarity" )
plt.ylabel("Structure Counts")
plt.yscale('log')
plt.ylim(1e0,1e6)
plt.legend(fontsize= 15 , loc = "upper center" ,bbox_to_anchor=(0.5,1.02),ncol= 6)
plt.show()
#%%


# t_windows_test = [4,50,100,150,200,250,300,350,400,450,500]
# angle_threshold = 0.05

# theta_Bx = Bvecs[:,0]
# theta_By = Bvecs[:,1]
# theta_Bz = Bvecs[:,2]

# theta_vx = vivecs[:,0]
# theta_vy = vivecs[:,1]
# theta_vz = vivecs[:,2]


# theta_dBx = get_angles_t(t_windows_test, theta_Bx)
# theta_dBy = get_angles_t(t_windows_test, theta_By)
# theta_dBz = get_angles_t(t_windows_test, theta_Bz)

# theta_dvx = get_angles_t(t_windows_test, theta_vx)
# theta_dvy = get_angles_t(t_windows_test, theta_vy)
# theta_dvz = get_angles_t(t_windows_test, theta_vz)


# counts_Bx = get_counts(theta_dBx, angle_threshold)
# counts_By = get_counts(theta_dBy, angle_threshold)
# counts_Bz = get_counts(theta_dBz, angle_threshold)

# counts_vx = get_counts(theta_dvx, angle_threshold)
# counts_vy = get_counts(theta_dvy, angle_threshold)
# counts_vz = get_counts(theta_dvz, angle_threshold)


# target_edges = np.array([1, 10, 100, 1000])
# Re_per_min = sw_km_to_Re(windspeed)
# t_windows = np.floor(d_windows / Re_per_min)
# x_midpoints = (d_windows[:-1] + d_windows[1:]) / 2


# counts_dbx_binned = aggregate_counts_by_bounds(x_midpoints, counts_Bx, target_edges)
# counts_dby_binned = aggregate_counts_by_bounds(x_midpoints, counts_By, target_edges)
# counts_dbz_binned = aggregate_counts_by_bounds(x_midpoints, counts_Bz, target_edges)
# counts_dvx_binned = aggregate_counts_by_bounds(x_midpoints, counts_vx, target_edges)
# counts_dvy_binned = aggregate_counts_by_bounds(x_midpoints, counts_vy, target_edges)
# counts_dvz_binned = aggregate_counts_by_bounds(x_midpoints, counts_vz, target_edges)


# labels = ["1-10", "10-100", "100-1000"]
# x = np.arange(len(labels))

# plt.figure(figsize=(10, 6))
# # x - 0.2 shifting bar to the left
# width = 0.16

# plt.bar(x - 0.25, counts_dbx_binned, width=width, label=r"Sw $\delta B_x$")
# plt.bar(x - 0.15, counts_dby_binned, width=width, label=r"Sw $\delta B_y$")
# plt.bar(x - 0.05, counts_dbz_binned, width=width, label=r"Sw $\delta B_z$")

# plt.bar(x + 0.05, counts_dvx_binned, width=width, label=r"Sw $\delta V_x$")
# plt.bar(x + 0.15, counts_dvy_binned, width=width, label=r"Sw $\delta V_y$")
# plt.bar(x + 0.25, counts_dvz_binned, width=width, label=r"Sw $\delta V_z$")
# plt.xticks(x, labels)
# plt.title("Solar Wind Structures Above A Set Threshold")
# plt.xlabel(fr"Scale Size ($R_E$) Choosen Fliter: {threshold}" "%" "\n" 
# 		   "Note: Bars are horizontally offset for visual clarity" )
# plt.ylabel("Structure Counts")
# plt.yscale('log')
# plt.ylim(0,1e6)
# plt.legend(fontsize= 10 , loc = "upper center" ,bbox_to_anchor=(0.5,1.02),ncol= 6)
# plt.show()

# %%
# Cell for plotting Variable vs Counts
#

# %%
