 #%%
import pyspedas
import numpy as np
import matplotlib.pyplot as plt
import sys
import math 

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
def get_angles(t_windows,vecvar):
	angle_array = []
	for i in range(len(t_windows) - 1):
		angle = np.zeros(len(timeax))

		T_small = t_windows[i]
		T_large = t_windows[i+1]

		vec_small = get_vecmean(vecvar, T_small*minuteint)
		vec_large = get_vecmean(vecvar, T_large*minuteint)

		for j in range(len(vec_small)):
			angle[j] = get_angle(vec_small[j],vec_large[j])

		angle_array.append(angle)
	return angle_array



#%%

trange=['2019-01-01/00:00', '2019-01-20/00:00']
# trange=['2019-01-20/00:00', '2019-01-30/00:00']
swe_vars = pyspedas.projects.wind.swe(trange=trange)
# denstiy 
mfi_vars = pyspedas.projects.wind.mfi(trange=trange)

# Define variable names for convenience
B_name='BGSE'
ni_name='N_elec'
vi_name = 'U_eGSE'
pos_name = 'DIST'

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

#%%
windspeed = np.nanmean(1e-3*vmag) # Get mean wind speed mag in km/s
v = sw_km_to_Re(windspeed) # Convert wind speed from km/s to Re/min

# d_windows=np.array([4,10,100,1000]) # Re
d_windows=np.array([4,50,100,150,200,250,300,350,400,450,500]) # Re
t_windows = []
for i in d_windows: t_windows.append(int(np.round(i/v)))
t_windows = np.array(t_windows)


t_windows_test = t_windows

dBmag = get_scaldeltas(t_windows_test, Bmag)
dni = get_scaldeltas(t_windows_test, ni)
dvmag = get_scaldeltas(t_windows_test, vmag)
dBvecs = get_vecdeltas(t_windows_test, Bvecs)
dvvecs = get_vecdeltas(t_windows_test, vivecs)
dB_angles = get_angles(t_windows_test, Bvecs)
dv_angles = get_angles(t_windows_test, vivecs) 

#%%
# Example of a mask with multiple conditions
threshold = 0.05
angle_threshold = 25 # degrees
mask0 = (dvmag[0] >= threshold) & (dBmag[0]>=threshold) & (dB_angles[0]>=angle_threshold)

# Plot Fluctuations
plt.plot(dBmag[0],label='$\delta B/B$')
plt.plot(dvmag[0],label='$\delta v/v$')
plt.xlabel('Time (unitless)')
plt.axhline(y=threshold, color='r', linestyle='-',label = r'Fluctuation threshold')
plt.legend()
plt.ylim(0,1)

print("Number of points above threshold: ", len(dBmag[0][mask0]))

# %%
# Plot fluctuation angles
plt.plot(dB_angles[0],label=r'$\theta_{B}$')
plt.axhline(y=angle_threshold, color='r', linestyle='-')
plt.legend()
plt.ylim(0,180)
# %%

# Next steps: Histograms...
