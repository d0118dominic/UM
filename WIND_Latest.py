
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
	for i in range(len(t_windows) - 1):
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

#%%

trange=['2019-01-01/00:00', '2019-01-20/00:00']
# trange=['2019-01-20/00:00', '2019-01-30/00:00']
swe_vars = pyspedas.projects.wind.swe(trange=trange)
# denstiy 
mfi_vars = pyspedas.projects.wind.mfi(trange=trange)

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


#%%

#%%
windspeed = np.nanmean(1e-3*vmag) # Get mean wind speed mag in km/s
v = sw_km_to_Re(windspeed) # Convert wind speed from km/s to Re/min

# d_windows=np.array([4,10,100,1000]) # Re
d_windows = np.array([4,10,100,1000,10000])
t_windows = []
for i in d_windows: t_windows.append(int(np.round(i/v)))
t_windows = np.array(t_windows)

t_windows_test = t_windows
#%%
dBmag = get_scaldeltas(t_windows_test, Bmag)
dni = get_scaldeltas(t_windows_test, ni)
dvmag = get_scaldeltas(t_windows_test, vmag)
dBvecs = get_vecdeltas(t_windows_test, Bvecs)
dvvecs = get_vecdeltas(t_windows_test, vivecs)
dB_angles = get_angles_t(t_windows_test, Bvecs)
dv_angles = get_angles_t(t_windows_test, vivecs) 


#%%
# threshold = 0.05
# target_edges = np.array([1,10,100,1000,10000])
# x_midpoints = (d_windows[:-1] + d_windows[1:]) / 2
# counts_db = get_counts(dBmag, threshold)

# counts_dv      = get_counts(dvmag, threshold)
# counts_dni     = get_counts(dni, threshold)
# counts_dBvecs  = get_counts(dBvecs, threshold)
# counts_dvvecs  = get_counts(dvvecs, threshold)

# counts_db_binned = aggregate_counts_by_bounds(x_midpoints, counts_db, target_edges)
# counts_dv_binned = aggregate_counts_by_bounds(x_midpoints, counts_dv, target_edges)
# counts_dni_binned = aggregate_counts_by_bounds(x_midpoints, counts_dni, target_edges)
# counts_dBvecs_binned = aggregate_counts_by_bounds(x_midpoints, counts_dBvecs, target_edges)
# counts_dBvecs_binned = aggregate_counts_by_bounds(x_midpoints, counts_dvvecs, target_edges)
# # Example of a mask with multiple conditions
# labels = ["1-10", "10-100", "100-1000","1000-10000"]
# # "1000-10000"
# x = np.arange(len(d_windows))
# counts_dni_mean = np.mean(counts_dni)

# Counts_db_plt = counts_db_binned[0] + counts_db_binned[1] + counts_db_binned[2] +  counts_db_binned[3]
# Counts_dv_plt = counts_dv_binned[0] + counts_dv_binned[1] + counts_dv_binned[2] + counts_dv_binned[3]
# Counts_dni_plt = counts_dni_binned[0] + counts_dni_binned[1] + counts_dni_binned[2] + counts_dni_binned[3]


# # plt.figure(figsize=(9, 5))
# # x - 0.2 shifting bar to the left
# plt.bar(x + 0.2 ,counts_db_binned, width=0.6, 
# color="skyblue", edgecolor="black", label=fr"Sw |B| Structures counts:{Counts_db_plt}")
# plt.bar(x + 0.3, counts_dv_binned, width=0.6,color="navy",
# 		edgecolor="black",label=fr"Sw |V| Structures : {Counts_dv_plt}")
# plt.bar(x, counts_dni_binned, width=0.6,color="pink",
# 		edgecolor="black",label=fr"Sw ni Structures counts: {Counts_dni_plt}")
# plt.xticks(x, labels)
# plt.title("Solar Wind Structures Above A Set Threshold")
# plt.xlabel(fr"Scale Size ($R_E$) Choosen Fliter: {threshold * 100}" "%" "\n" 
# 		   "Note: Bars are horizontally offset for visual clarity" )
# # plt.errorbar(x + 0.2,counts_db_binned,yerr=y_error,fmt='o')
# plt.ylabel("Structure Counts")
# plt.yscale('log')
# plt.ylim(0,1e5)
# plt.legend(fontsize=12)
# plt.show()

# labels = ['4-10 $R_E$', '10-100 $R_E$', '100-1000 $R_E$', '1000-10000 $R_E$']
# plt.figure(figsize=(8,5))
# plt.plot(counts_db_binned,color='blue',marker='o', label=fr"$\delta B/B$ counts : {Counts_db_plt}")
# plt.plot(counts_dv_binned,color='red',marker='o',label=fr"$\delta v/v$ counts {Counts_dv_plt} ")
# plt.plot(counts_dni_binned,color='purple',marker='o',label=fr"$\delta n/n$ counts {Counts_dni_plt}")
# plt.xticks([0, 1, 2, 3], labels)
# plt.yscale('log')
# plt.ylim(0,1e5)
# plt.xlabel('Scale Size ($R_E$)')
# plt.ylabel('Number of Fluctuations Above Threshold')
# plt.title(fr'Fluctuation Counts (Threshold = {threshold})')
# plt.legend(fontsize=16)
# plt.show()
#%%

# counts_db = get_counts(dBmag, threshold)
# counts_dv      = get_counts(dvmag, threshold)
# counts_dni     = get_counts(dni, threshold)
# counts_dBvecs  = get_counts(dBvecs, threshold)
# counts_dvvecs  = get_counts(dvvecs, threshold)



# plt.figure(figsize=(8,5))
# plt.plot(counts_db, marker='o', label=r'$\delta B/B$')
# plt.plot(counts_dv, marker='o', label=r'$\delta v/v$')
# plt.plot(counts_dni, marker='o', label=r'$\delta n/n$')
# # plt.hist(counts_dBvecs, marker='o', label=r'$\Delta \mathbf{B}$ rotation')
# # plt.hist(counts_dvvecs, marker='o', label=r'$\Delta \mathbf{v}$ rotation')
# plt.xlabel('Time Window Index')
# plt.ylabel('Number of Fluctuations Above Threshold')
# plt.title(fr'Fluctuation Counts (Threshold = {threshold})')
# plt.legend()
# plt.show()
#%%



# angle_threshold = 25 # degrees
# mask0 = (dvmag[0] >= threshold) & (dBmag[0]>=threshold) & (dB_angles[0]>=angle_threshold)

# # Plot Fluctuations
# plt.plot(dBmag[0], label=r'$\delta B/B$ (Magnetic fluctuation)')
# plt.plot(dvmag[0], label=r'$\delta v/v$ (Velocity fluctuation)')
# plt.xlabel('Time Index')
# plt.title(fr'Magnetic and Velocity Fluctuations ($\theta > {angle_threshold}^\circ$)')
# plt.axhline(y=threshold, color='r', linestyle='-',label = r'Fluctuation threshold')
# plt.legend()
# plt.ylim(0,1)

# print("Number of points above threshold: ", len(dBmag[0][mask0]))

# %%
# Plot fluctuation angles
# plt.plot(dB_angles[0],label=r'$\theta_{B}$')
# plt.axhline(y=angle_threshold, color='r', linestyle='-')
# plt.legend()
# plt.ylim(0,180)
# %%
#Condition



#cond for multi thresholds 

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

def make_d_windows(start = 0,stop=100,increment= 1000):
	n = 11
	a = np.zeros(n)
	for i in range(len(a)):
		a[i] = 100*i
	a[0] = 4
	return a

def make_labels(d_windows):
    return [f"{d_windows[i]}-{d_windows[i+1]}" for i in range(len(d_windows) - 1)]

d_windows = make_d_windows()

# Get t_windows (list of integers)
t_windows_low = []
for i in d_windows: t_windows_low.append(int(np.round(i/v_low)))
t_windows = np.array(t_windows_low)


t_windows_high = []
for i in d_windows: t_windows_high.append(int(np.round(i/v_high)))
t_windows = np.array(t_windows_high)


t_windows_mean = []
for i in d_windows: t_windows_mean.append(int(np.round(i/v_mean)))
t_windows = np.array(t_windows_mean)
print("High: ",v_high)
print("Low: ", v_low)
print("Mean: ", sw_mean)


print("Mean: " ,len(t_windows_mean))
print("High:" ,len(t_windows_high))
print("Low: " ,len(t_windows_low))
#%%



# Get the Delta arrays with the t_windows vals
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
# Set a percent threshold 
threshold = 0.05
TH = threshold




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

# print("counts_db_low: " , counts_db_low)
# print("counts_dv_low: " , counts_dv_low)
# print("counts_dni_low: " , counts_dni_low)
# print("counts_dbvec_low: " , counts_dBvecs_low)
# print("counts_dvvec_low: " , counts_dvvecs_low)







#%%
# All the binning and plotting methods (idk the details)

# print(r"Bin db :", counts_db_binned)
# labels = ["1-10", "10-100", "100-1k","1k-10k"]
# labels = ["1-100", "100-200", "200-300","300-400"]
# labels = ["1-100", "100-200", "200-300","300-400","400-500","500-600","600-700","700-800","800-900","900-1000", 
# "1k-1.1k","1.1-1.2k","1.2-1.3k","1.3-1.4k","1.4-1.5k","1.5-1.6k","1.6-1.7k","1.7-1.8k","1.8-1.9k","1.9-2k"]
# labels = ["1-100", "100-200", "200-300","300-400","400-500","500-600","600-700","700-800","800-900","900-1000"] 
# labels = ["1-2k", "2-3k", "3-4k","4-5k","5-6k","6-7k","7-8k","8-9k","9-10k"]
# 

# def y_error_bars(var,counts,p):
# 	y_error = []

# 	for var in counts:
		
# 		result = var * p

# 		y_error.append(result)
# 	return y_error

#==========================================================================#

target_edges = d_windows 
x_midpoints = (d_windows[:-1] + d_windows[1:]) / 2

labels = labels = make_labels(d_windows) if 'd_windows' in locals() else [f"W{i+1}" for i in range(10)]

x = np.arange(10)
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


#%%
#just for db 
high_db =  np.array(counts_db_high)
low_db =  np.array(counts_db_low)
mean_db =  np.array(counts_db)


lower_error = np.abs(mean_db - low_db)
upper_error =  np.abs(high_db - mean_db)
asy_error_db = [lower_error,upper_error]
#%%

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
plt.yscale('log')
plt.xlim()
plt.ylim(1,100000)
plt.show()
#%%


angle_threshold = 25 # degrees
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

# %%
#%%

p = 0.25

dni_error_bars = y_error_bars(dni,counts_dni,p)
bmag_error_bars = y_error_bars(dBmag,counts_dv,p)
vmag_error_bars = y_error_bars(dvmag,counts_db,p)

# xy = [dni * 0.25 for dni in counts_dni ]
# yz = [dvmag * 0.25 for dvmag in counts_dv]
# zx_error = [dBmag * 0.25 for dBmag in counts_db]

labels = make_labels(d_windows) 

x = np.arange(len(labels))

# plt.bar(x-0.25,counts_dni_high,label = 'Density Structures')
	    # ,width=0.5, yerr=50,error_kw ={'capsize': 5, 'ecolor': 'black'} )
# counts_dni_mean = np.mean(counts_dni)
# print(r"counts_dni: " ,counts_dni_mean)

# I'm not including the labels here because they crowd the axis 
# plt.figure(figsize=(14, 10))
# plt.bar(x-0.25,counts_dni,label = 'Density Structures' ,width=0.5, yerr=dni_error_bars ,error_kw ={'capsize': 5, 'ecolor': 'black'} )
# plt.bar(x,counts_db,label = 'Magnetic Structures',width=0.5 , yerr =bmag_error_bars ,error_kw ={'capsize': 5, 'ecolor': 'black'} )
# plt.bar(x+0.25,counts_dv,label = 'Velocity Structures',width=0.5,yerr = vmag_error_bars , error_kw={'capsize' : 5,'ecolor': 'black'})
# plt.xticks(x,labels)
# plt.xlabel(fr"Mesoscale Occurance Rate  ($R_E$) Thershold: {TH * 100}%",fontsize=20)
# plt.xlim(-1,10)
# plt.legend(fontsize=20)
# plt.yscale('log')




#%%

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

counts_db_mask = get_counts(dBmag_mean, threshold)
counts_dv_mask  = get_counts(dvmag_mean, threshold)
counts_dni_mask   = get_counts(dni_mean, threshold)
counts_dBvecs_mask = get_counts(dBvecs_mean, threshold)
counts_dvvecs_mask  = get_counts(dvvecs_mean, threshold)

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



width1 = 0.7
Move_Bar_Left = 0.05
Move_Bar_Right = 0.05
plt.figure(figsize=(12, 5))
plt.bar(arranged_clean_Bmag_A - Move_Bar_Left , clean_Bscal_A,
width=width1, alpha=0.7, color="blue", label=r"$\Delta B_{scal}$")
plt.bar(arranged_clean_dni_A + Move_Bar_Right, clean_dni_A,alpha=0.7, width=width1,color = "pink", label=r"$\Delta DNI$")
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
plt.tight_layout()
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
plt.title(fr"Scale Size ($R_E$) vmag >= {threshold_correlation * 100} AND  dni >= {threshold_correlation * 100} ", fontsize=12)
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
# Multiple histograms for multiple scale ranges
# This method should probably altered if there are a lot of scale bins
# Example here uses 4 

# Adjust Thresholds on dB/B and dv/v.  
# If you don't want any threshold on one, set it to zero
Bthreshold = 0.05
vthreshold = 0.02

# Set number of bins (50 seems fine)
bins = 50

# Set your mask to filter the data  
# Right now its set up to exclude all data that doesn't meet both thresholds
mask = [(dBmag[i] >= Bthreshold) & (dvmag[i] >= vthreshold) for i in range(len(dBmag))]

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
# %%
