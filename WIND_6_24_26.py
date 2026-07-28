
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
def get_angles_t(t_windows,vecvar):
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

# A version of get counts which can take multiple thresholds on normalized B,v, and n
# They default to zero so if you don't define one there is automatically no threshold on that variable
def get_counts(data_arr,threshold):
	counts_list = []
	#threshold = 0.05
	for i in range(len(t_windows) - 1):
		mask = (data_arr[i]>=threshold)
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
d_windows = np.array([4,10,100,1000,10000])
t_windows = []
for i in d_windows: t_windows.append(int(np.round(i/v)))
t_windows = np.array(t_windows)

t_windows_test = t_windows

dBmag = get_scaldeltas(t_windows_test, Bmag)
dni = get_scaldeltas(t_windows_test, ni)
dvmag = get_scaldeltas(t_windows_test, vmag)
dBvecs = get_vecdeltas(t_windows_test, Bvecs)
dvvecs = get_vecdeltas(t_windows_test, vivecs)
dB_angles = get_angles_t(t_windows_test, Bvecs)
dv_angles = get_angles_t(t_windows_test, vivecs) 

#%%
# Example of a mask with multiple conditions
# threshold = 0.05
# angle_threshold = 25 # degrees
# mask0 = (dvmag[0] >= threshold) & (dBmag[0]>=threshold) & (dB_angles[0]>=angle_threshold)

# # Plot Fluctuations
# plt.plot(dBmag[0],label='$\delta B/B$')
# plt.plot(dvmag[0],label='$\delta v/v$')
# plt.xlabel('Time (unitless)')
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

#Define d_windows in Re
d_windows = np.array([4,10,100,1000,10000])

# Define d_windows in Re/min
# d_windows = np.array([4,100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000])
# d_windows = np.array([4,100,200,300,400,500,600,700,800,900,1000])#,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000])
# d_windows = np.array([1000,2000,3000,4000,5000,6000,7000,8000,9000,10000])#,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000])

# def make_d_windows(start = 0,stop=10,increment = 100):
# 	n = 101
# 	a = np.zeros(n)
# 	for i in range(len(a)):
# 		a[i] = 10*i
# 	a[0] = 4
# 	return a

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

# Set v to either sw_mean, sw_low, or sw_high
v =  sw_mean


def make_labels(d_windows):
    return [f"{d_windows[i]}-{d_windows[i+1]}" for i in range(len(d_windows) - 1)]

# d_windows = make_d_windows()

# Get t_windows (list of integers)
t_windows = []
for i in d_windows: t_windows.append(int(np.round(i/v)))
t_windows = np.array(t_windows)


#%%

# Get the Delta arrays with the t_windows vals
dBmag = get_scaldeltas(t_windows , Bmag)
dni = get_scaldeltas(t_windows , ni)
dvmag = get_scaldeltas(t_windows , vmag)
dBvecs = get_vecdeltas(t_windows, Bvecs)
dvvecs = get_vecdeltas(t_windows, vivecs)
dB_angles = get_angles_t(t_windows, Bvecs)
dv_angles = get_angles_t(t_windows, vivecs) 
#%%

# Cell for plotting Variable vs Counts
# Multiple histograms for multiple scale ranges
# This method should probably altered if there are a lot of scale bins
# Example here uses 4 

# Adjust Thresholds on dB/B and dv/v.  
# If you don't want any threshold on one, set it to zero
Bthreshold = 0.0
vthreshold = 0.0

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

#%%

#%%



#Errorbar method example (seems to work, so long as high>mean>low)

# x = np.arange(10)


# mean = [4, 5,3, 0, 0, 0, 0, 0, 0, 0]

# high = [7, 6, 3.6, 0, 0, 0, 0, 0, 0, 0]
# low = [3, 1, 2.8, 0, 0, 0, 0, 0, 0, 0]


# lower_error = (np.abs(np.array(mean) - np.array(low)))
# upper_error = (np.abs(high) - np.array(mean))
# asy_error = [lower_error,upper_error]


# plt.figure(figsize=(15,8))
# plt.bar(x,mean,color ="royalblue",label="sw:mean")
# plt.errorbar(x, mean, yerr=np.abs(asy_error), fmt='none', 
#              ecolor='navy', elinewidth=2, capsize=5, label="solar wind : range")


# plt.legend()


# plt.yscale('linear')


#%%
# Set a percent threshold 
# threshold = 0.1

bthreshold = 0.2
vthreshold = 0.05
nthreshold = 0.01


# How can we change this to handle more aribitrary masks?

def get_counts(data_arr,mask):
	mask = [(dBmag[i] >= Bthreshold) & (dvmag[i] >= vthreshold) & (dni[i]>=nthreshold) for i in range(len(dBmag))]
    counts_list = []
    #threshold = 0.05
    for i in range(len(t_windows) - 1):
        number1 = len(data_arr[i][mask[i]])
        counts_list.append(number1)
    return counts_list



# # Get counts arrays
counts_db      = get_counts(dBmag, mask)
counts_dv      = get_counts(dvmag, mask)
counts_dni     = get_counts(dni, mask)
counts_dBvecs  = get_counts(dBvecs, mask)
counts_dvvecs  = get_counts(dvvecs, mask)





# #%%





# #%%
# # All the binning and plotting methods (idk the details)
# target_edges = d_windows 
# x_midpoints = (d_windows[:-1] + d_windows[1:]) / 2

# # counts_db_binned = aggregate_counts_by_bounds(x_midpoints, counts_db, target_edges)
# # counts_dv_binned = aggregate_counts_by_bounds(x_midpoints, counts_dv, target_edges)
# # counts_dni_binned = aggregate_counts_by_bounds(x_midpoints, counts_dni, target_edges)
# # counts_dBvecs_binned = aggregate_counts_by_bounds(x_midpoints, counts_dBvecs, target_edges)
# # counts_dBvecs_binned = aggregate_counts_by_bounds(x_midpoints, counts_dvvecs, target_edges)
# #counts_angles_binned = aggregate_counts_by_bounds(x_midpoints, counts_angles, target_edges)
# lower_error = [0,0,0]
# upper_error = [5,5,5]
# y_error = 50
# # print(r"Bin db :", counts_db_binned)
# # labels = ["1-10", "10-100", "100-1k","1k-10k"]
# # labels = ["1-100", "100-200", "200-300","300-400"]
# # labels = ["1-100", "100-200", "200-300","300-400","400-500","500-600","600-700","700-800","800-900","900-1000", 
# # "1k-1.1k","1.1-1.2k","1.2-1.3k","1.3-1.4k","1.4-1.5k","1.5-1.6k","1.6-1.7k","1.7-1.8k","1.8-1.9k","1.9-2k"]
# # labels = ["1-100", "100-200", "200-300","300-400","400-500","500-600","600-700","700-800","800-900","900-1000"] 
# # labels = ["1-2k", "2-3k", "3-4k","4-5k","5-6k","6-7k","7-8k","8-9k","9-10k"]
# # 
# labels = make_labels(d_windows) 

# x = np.arange(len(labels))

# # counts_dni_mean = np.mean(counts_dni)
# # print(r"counts_dni: " ,counts_dni_mean)

# # I'm not including the labels here because they crowd the axis 
# plt.bar(x,counts_dni,label = 'Density Structures')
# plt.bar(x,counts_db,label = 'Magnetic Structures')
# plt.bar(x,counts_dv,label = 'Velocity Structures')
# plt.xlim(-1,20)
# plt.ylim(1,1e4)
# plt.legend()
# plt.yscale('log')



















# #%%

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

# # %%

# # Choose an angle threshold for the magnetic and velocity deflections 
# # (We can choose different ones)
# B_angle_threshold = 45
# V_angle_threshold = 5

# # Get Angle Arrays
# angle_array_db = get_angles_t(t_windows,Bvecs)
# angle_array_vi = get_angles_t(t_windows,vivecs)

# # Get Counts of each above their defined thresholds
# counts_angles_db = get_counts(angle_array_db,B_angle_threshold) 
# counts_angles_vi = get_counts(angle_array_vi,V_angle_threshold) 

# # Binning & Plotting
# counts_angles_binned_db = aggregate_counts_by_bounds(x_midpoints, counts_angles_db,
# 							 target_edges)

# counts_angles_binned_vi = aggregate_counts_by_bounds(x_midpoints, counts_angles_vi,
# 							 target_edges)



# labels = labels #Using the same definition from above
# label = np.arange(len(labels))

# plt.bar(label + 0.1 ,counts_angles_binned_db, width=0.3, 
# color="skyblue", edgecolor="black", label=fr'$\theta_B > {B_angle_threshold}^o$ Structures')
# plt.bar(label - 0.1 ,counts_angles_binned_vi, width=0.3, 
# color="pink", edgecolor="black", label=fr'$\theta_V > {V_angle_threshold}^o$ Structures')
# # plt.bar(label + 0.2 ,counts_angles_binned_ni, width=0.6, 
# # color="skyblue", edgecolor="black", label=fr'$\theta > {angle_threshold}^o$ Structures')
# plt.xticks(label, labels)
# plt.xlabel(r"Scale Size ($R_E$)",fontsize=16)
# plt.ylabel("Counts" ,fontsize=16)
# plt.yscale('log')
# plt.ylim(1e0,1e6)
# plt.legend(fontsize=14)
# plt.show()

# # %%
# threshold = 0.4


# Bx = Bvecs[:,0]
# By = Bvecs[:,1]
# Bz = Bvecs[:,2]

# vx = vivecs[:,0]
# vy = vivecs[:,1]
# vz = vivecs[:,2]


# dBx = get_scaldeltas(t_windows, Bx)
# dBy = get_scaldeltas(t_windows, By)
# dBz = get_scaldeltas(t_windows, Bz)

# dvx = get_scaldeltas(t_windows, vx)
# dvy = get_scaldeltas(t_windows, vy)
# dvz = get_scaldeltas(t_windows, vz)


# counts_Bx = get_counts(dBx, threshold)
# counts_By = get_counts(dBy, threshold)
# counts_Bz = get_counts(dBz, threshold)

# counts_vx = get_counts(dvx, threshold)
# counts_vy = get_counts(dvy, threshold)
# counts_vz = get_counts(dvz, threshold)

# counts_dbx_binned = aggregate_counts_by_bounds(x_midpoints, counts_Bx, target_edges)
# counts_dby_binned = aggregate_counts_by_bounds(x_midpoints, counts_By, target_edges)
# counts_dbz_binned = aggregate_counts_by_bounds(x_midpoints, counts_Bz, target_edges)
# counts_dvx_binned = aggregate_counts_by_bounds(x_midpoints, counts_vx, target_edges)
# counts_dvy_binned = aggregate_counts_by_bounds(x_midpoints, counts_vy, target_edges)
# counts_dvz_binned = aggregate_counts_by_bounds(x_midpoints, counts_vz, target_edges)



# labels = labels
# x = np.arange(len(labels))


# #Plotting the dB/B components
# plt.figure(figsize=(8, 4))
# # x - 0.2 shifting bar to the left
# width = 0.16

# plt.bar(x - 0.1, counts_dbx_binned, width=width, label=r"Sw $\delta B_x/B_x$")
# plt.bar(x + 0.00, counts_dby_binned, width=width, label=r"Sw $\delta B_y/B_y$")
# plt.bar(x + 0.1, counts_dbz_binned, width=width, label=r"Sw $\delta B_z/B_z$")

# plt.xticks(x, labels)
# plt.title("Solar Wind Component Fluctuations")
# plt.xlabel(fr"Scale Size ($R_E$) Choosen Fliter: {threshold * 100}" "%" "\n" 
# 		   "Note: Bars are horizontally offset for visual clarity" )
# plt.ylabel("Structure Counts")
# plt.yscale('log')
# plt.ylim(1e3,1e6)
# plt.legend(fontsize= 15 , loc = "upper center" ,bbox_to_anchor=(0.5,1.02),ncol= 6)
# plt.show()



# #%%
# #Plotting the dv/v components

# plt.figure(figsize=(8, 4))
# # x - 0.2 shifting bar to the left
# width = 0.16


# plt.bar(x - 0.1, counts_dvx_binned, width=width, label=r"Sw $\delta V_x/V_x$")
# plt.bar(x + 0.00, counts_dvy_binned, width=width, label=r"Sw $\delta V_y/V_y$")
# plt.bar(x + 0.1, counts_dvz_binned, width=width, label=r"Sw $\delta V_z/V_z$")
# plt.xticks(x, labels)
# plt.title("$\delta v/v$ Components")
# plt.xlabel(fr"Scale Size ($R_E$) Choosen Fliter: {threshold * 100}" "%" "\n" 
# 		   "Note: Bars are horizontally offset for visual clarity" )
# plt.ylabel("Structure Counts")
# plt.yscale('log')
# plt.ylim(1e0,1e6)
# plt.legend(fontsize= 15 , loc = "upper center" ,bbox_to_anchor=(0.5,1.02),ncol= 6)
# plt.show()






# #%%


# # t_windows_test = [4,50,100,150,200,250,300,350,400,450,500]
# # angle_threshold = 0.05

# # theta_Bx = Bvecs[:,0]
# # theta_By = Bvecs[:,1]
# # theta_Bz = Bvecs[:,2]

# # theta_vx = vivecs[:,0]
# # theta_vy = vivecs[:,1]
# # theta_vz = vivecs[:,2]


# # theta_dBx = get_angles_t(t_windows_test, theta_Bx)
# # theta_dBy = get_angles_t(t_windows_test, theta_By)
# # theta_dBz = get_angles_t(t_windows_test, theta_Bz)

# # theta_dvx = get_angles_t(t_windows_test, theta_vx)
# # theta_dvy = get_angles_t(t_windows_test, theta_vy)
# # theta_dvz = get_angles_t(t_windows_test, theta_vz)


# # counts_Bx = get_counts(theta_dBx, angle_threshold)
# # counts_By = get_counts(theta_dBy, angle_threshold)
# # counts_Bz = get_counts(theta_dBz, angle_threshold)

# # counts_vx = get_counts(theta_dvx, angle_threshold)
# # counts_vy = get_counts(theta_dvy, angle_threshold)
# # counts_vz = get_counts(theta_dvz, angle_threshold)


# # target_edges = np.array([1, 10, 100, 1000])
# # Re_per_min = sw_km_to_Re(windspeed)
# # t_windows = np.floor(d_windows / Re_per_min)
# # x_midpoints = (d_windows[:-1] + d_windows[1:]) / 2


# # counts_dbx_binned = aggregate_counts_by_bounds(x_midpoints, counts_Bx, target_edges)
# # counts_dby_binned = aggregate_counts_by_bounds(x_midpoints, counts_By, target_edges)
# # counts_dbz_binned = aggregate_counts_by_bounds(x_midpoints, counts_Bz, target_edges)
# # counts_dvx_binned = aggregate_counts_by_bounds(x_midpoints, counts_vx, target_edges)
# # counts_dvy_binned = aggregate_counts_by_bounds(x_midpoints, counts_vy, target_edges)
# # counts_dvz_binned = aggregate_counts_by_bounds(x_midpoints, counts_vz, target_edges)


# # labels = ["1-10", "10-100", "100-1000"]
# # x = np.arange(len(labels))

# # plt.figure(figsize=(10, 6))
# # # x - 0.2 shifting bar to the left
# # width = 0.16

# # plt.bar(x - 0.25, counts_dbx_binned, width=width, label=r"Sw $\delta B_x$")
# # plt.bar(x - 0.15, counts_dby_binned, width=width, label=r"Sw $\delta B_y$")
# # plt.bar(x - 0.05, counts_dbz_binned, width=width, label=r"Sw $\delta B_z$")

# # plt.bar(x + 0.05, counts_dvx_binned, width=width, label=r"Sw $\delta V_x$")
# # plt.bar(x + 0.15, counts_dvy_binned, width=width, label=r"Sw $\delta V_y$")
# # plt.bar(x + 0.25, counts_dvz_binned, width=width, label=r"Sw $\delta V_z$")
# # plt.xticks(x, labels)
# # plt.title("Solar Wind Structures Above A Set Threshold")
# # plt.xlabel(fr"Scale Size ($R_E$) Choosen Fliter: {threshold}" "%" "\n" 
# # 		   "Note: Bars are horizontally offset for visual clarity" )
# # plt.ylabel("Structure Counts")
# # plt.yscale('log')
# # plt.ylim(0,1e6)
# # plt.legend(fontsize= 10 , loc = "upper center" ,bbox_to_anchor=(0.5,1.02),ncol= 6)
# # plt.show()

# # %%
