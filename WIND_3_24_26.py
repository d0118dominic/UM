 #%%
import pyspedas
import numpy as np
import matplotlib.pyplot as plt
import sys
import math 

from pyspedas import get_data,store_data
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



# List of Functions 


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
#funciton testing 
# Vx = vivecs[:, 0]
# Km_sw_test = V[:, 0]
# sw_km_to_Re(Km_sw_test)
# print(sw_km_to_Re(Km_sw_test)) 

# #############

t_windows_test = [1,2,3,5,6] # mins

print(Vsw_Tw(t_windows_test))


#%%
# WIND data (if you want to use WIND data, run this cell and skip the next cell)

# highest time resolutions for swe instrument: 92s (ions) & 6-12 s (electrons) 
# Assuming solar wind velocity of 400 km/s, this corresponds to structure sizes of 6 & 0.4-0.8 Re

trange=['2019-01-01/00:00', '2019-01-30/00:00']
# trange=['2019-01-20/00:00', '2019-01-30/00:00']
swe_vars = pyspedas.projects.wind.swe(trange=trange)
# denstiy 
mfi_vars = pyspedas.projects.wind.mfi(trange=trange)


# tdp_vars = pyspedas.projects.wind.threedp(trange=trange)

# tplot(['N_elec', 'T_elec','U_eGSE','BGSE'])

B_name='BGSE'
ni_name='N_elec'
vi_name = 'U_eGSE'
pos_name = 'DIST' # Distance from earth, in units of earth radii.  1 Re ~ 4.3e-5 AU





#%%

# ACE data

# highest time resolutions for swe instrument: 64s (ions) & 128 s (electrons) 
# Assuming solar wind velocity of 400 km/s, this corresponds to structure sizes of 4 & 8 Re

# trange = ['2018-11-9/00:00', '2018-11-20/00:00']
# swe_vars = pyspedas.projects.ace.swe(trange=trange)
# mfi_vars = pyspedas.projects.ace.mfi(trange=trange)

# B_name = 'BGSEc'
# ni_name = 'Np'
# vi_name = 'Vp'
# pos_name = 'SC_pos_GSE' # Disance from earth in units of km.  1 km ~ 6.7e-9 AU


# tplot([B_name,ni_name,vi_name,pos_name])





#%%
# Tplot variable basics

# To do more than plot the timeseries data directly, we need to extract the data from the tplot variable format

# We can use the get_data routine from pyspedas to retrieve the actual data structure
Bdata = get_data(B_name)
ndata = get_data(ni_name)
vdata = get_data(vi_name)

# If you look at the length of each part of the Bdata array, you'll see it has a weird format of [big number, same big number, 3]
print([len(Bdata[0]),len(Bdata[1]),len(Bdata[2])])

# The first (0th) element of Bdata is just a sequence of evenly spaced numbers corresponding to an array that represents time 
# (you can try plt.plot(Bdata[0]) and see its a line of constant slope)
# We can actually also call this time array with Bdata.times, whch is equivalent to Bdata[0]

# The 2nd element of Bdata is the actual data of interest
plt.plot(Bdata[1]) # This plot will show real data, which in this case is 3 components of the vector magnetic field (in nanoTesla units by default)

# I honestly don't recall what the purpose of the 3rd element is.  I don't know why it is necessary
# ndata doesn't have a 3rd element, so it has something to do with how many elements of the variable are measured at a given timestep (one density, but 3 components of the mag field)

# To turn a variable into a tplot array, we do the opposite of get_data, store_data.
store_data('Magnetic_Fields', data = {'x':Bdata.times,'y':Bdata[1]})


#%%

# Different instruments may have a different time resolution, some make more measurements per second than others
# Comparing the Magnetic field array and Density array (both covering the same time interval), you can see the magnetic field array is much longer (because it has higher time resolution)
print('Length of Magnetic Field Array: ' + str(len(Bdata[1])))
print('Length of Density Array: ' + str(len(ndata[1])))

# If we are interested in calculating quantities derived from multiple different instruments, we need all variables to have the same time resolution

# So let's choose which resolution we'll covert (interpolate) everything to.  Typically, it's best to interpolate to whichever variable has the lowest resolution (interpolating down)

# Here we choose the density tplot variable, assigning it to be the interpolation variable
# switching from ndata 
interpvar_name = B_name

# I like to define a timeax variable as well, representing the time variable
interpdata = get_data(interpvar_name)
timeax = interpdata.times


# To interpolate a tplot variable to another tplot variable, we use pyspedas's tinterpol and assign new names to the interpolated versions
tinterpol(B_name,interpvar_name,newname='B')
tinterpol(ni_name,interpvar_name,newname='ni') # (Obviously, interpolating a variable to itself isn't necessary, but I keep it here because sometimes I want to change the interpvar)
tinterpol(vi_name,interpvar_name,newname='vi') # 

# Now they should have the same length
B_name_new = get_data('B')
ni_name_new = get_data('ni')
vi_name_new = get_data('vi')
print('Length of New Magnetic Field Array: ' + str(len(B_name_new[1])))
print('Length of New Density Array: ' + str(len(ni_name_new[1])))


# The last thing to do before we start calculating things is reforming the tplot variables into simpler numpy arrays to work with.  This is why I made the reform() function.
# At the same time, I use factors to convert the variables from their default units to their SI units.  This just makes things easier in my experience 
Bvecs = 1e-9*reform(B_name_new) # Converting from nanotesla to Tesla
vivecs = 1e3*reform(vi_name_new) #converting km/s to m/s
ni = 1e6*reform(ni_name_new) # Converting from cm^-3 to m^-3

##%% 


# Let's start by calculating some things: the magnetic field magnitude and the alfven velocity

# define the arrays we'll overwrite
Bmag = np.zeros(len(timeax)) # this just creates the array of the correct length
vmag = np.zeros(len(timeax)) 
va = np.zeros(len(timeax))

# overwrite arrays with built-in and pre-defined functions
for i in range(len(timeax)):
	Bmag[i] = np.linalg.norm(Bvecs[i]) # np.linalg.norm automatically calculates the magnitude of a vector
	vmag[i] = np.linalg.norm(vivecs[i]) # np.linalg.norm automatically calculates the magnitude of a vector
	va[i] = get_va(Bvecs[i],ni[i],mi)  # using the get_va() functionVa = |B|/sqrt(mu0*n*m)  


# Now store new variables as tplot variables and plot them.  This is typically the step where I'd convert the units if I want to express something other than SI units.
store_data('Bvecs', data = {'x':timeax,'y':1e9*Bvecs}) # with a factor to convert Tesla to nanotesla 
store_data('Bmag', data = {'x':timeax,'y':1e9*Bmag}) # with a factor to convert Tesla to nanotesla 
store_data('vmag', data = {'x':timeax,'y':1e-3*vmag}) 
store_data('va', data = {'x':timeax,'y':1e-3*va}) # with a factor to convert m/s to km/s
tplot(['Bvecs','Bmag','va'])


#%%

# define minuteinterval, which represents the number of array elements equivalent to 1 minute (rounding up to the nearest integer)
minuteinterval = minute_int(timeax,trange)

# To change the averaging window, I like to just say how long (in minutes) i want the averaging window to be.  
# Try adjusting the value of n_minutes to see how the plots change (it has to be an integer). 
n_minutes_lg = 100
 # Larger mean interval
n_minutes_sml = 50 

# Now we can calculate the averaged values & store the tplot variables

# Long Average
Bvecs_lgmean = get_vecmean(Bvecs,n_minutes_lg*minuteinterval)
vivecs_lgmean = get_vecmean(vivecs,n_minutes_lg*minuteinterval)
Bmag_lgmean = get_mean(Bmag,n_minutes_lg*minuteinterval)
vmag_lgmean = get_mean(vmag,n_minutes_lg*minuteinterval)
n_lgmean = get_mean(ni,n_minutes_lg*minuteinterval)
va_lgmean = get_mean(va,n_minutes_lg*minuteinterval)
store_data('Bvecs_mean', data = {'x':timeax,'y':1e9*Bvecs_lgmean}) # converted to nT
store_data('vivecs_mean', data = {'x':timeax,'y':1e-3*vivecs_lgmean}) # converted to km/s
store_data('Bmag_mean', data = {'x':timeax,'y':1e9*Bmag_lgmean}) # converted to nT
store_data('vmag_mean', data = {'x':timeax,'y':1e-3*vmag_lgmean}) # converted to km/s
store_data('n_mean', data = {'x':timeax,'y':1e-6*n_lgmean}) # with a factor to convert m/s to km/s
store_data('va_mean', data = {'x':timeax,'y':1e-3*va_lgmean}) # with a factor to convert m/s to km/s

# Short Average
Bvecs_smlmean = get_vecmean(Bvecs,n_minutes_sml*minuteinterval)
vivecs_smlmean = get_vecmean(vivecs,n_minutes_sml*minuteinterval)
Bmag_smlmean = get_mean(Bmag,n_minutes_sml*minuteinterval)
vmag_smlmean = get_mean(vmag,n_minutes_sml*minuteinterval)
n_smlmean = get_mean(ni,n_minutes_sml*minuteinterval)
va_smlmean = get_mean(va,n_minutes_sml*minuteinterval)
store_data('Bvecs', data = {'x':timeax,'y':1e9*Bvecs_smlmean}) # converted to nT
store_data('vivecs', data = {'x':timeax,'y':1e-3*vivecs_smlmean}) # converted to km/s
store_data('Bmag', data = {'x':timeax,'y':1e9*Bmag_smlmean}) # converted to nT
store_data('vmag', data = {'x':timeax,'y':1e-3*vmag_smlmean}) # converted to km/s
store_data('n', data = {'x':timeax,'y':1e-6*n_smlmean}) # with a factor to convert m/s to km/s
store_data('va', data = {'x':timeax,'y':1e-3*va_smlmean}) # with a factor to convert m/s to km/s


# Storing non-averaged versions
store_data('Bvecs_inst', data = {'x':timeax,'y':1e9*Bvecs}) # converted to nT
store_data('Bmag_inst', data = {'x':timeax,'y':1e9*Bmag}) # converted to nT
#store_data('n_inst', data = {'x':timeax,'y':1e-6*n}) # with a factor to convert m/s to km/s
store_data('va_inst', data = {'x':timeax,'y':1e-3*va}) # with a factor to convert m/s to km/s



# Define empty arrays for the delta variables
dB,dB_norm = np.zeros_like(Bvecs),np.zeros_like(Bvecs)
dv,dv_norm = np.zeros_like(Bvecs),np.zeros_like(Bvecs)
dBmag,dBmag_norm = np.zeros_like(Bmag),np.zeros_like(Bmag)
dvmag,dvmag_norm = np.zeros_like(Bmag),np.zeros_like(Bmag)
dn,dn_norm = np.zeros_like(ni),np.zeros_like(ni)

# Get the deltas between the large and small mean
for i in range(len(timeax)):
	dB[i], dB_norm[i]= get_delta(Bvecs_smlmean[i],Bvecs_lgmean[i])  # dB & dB/|B| (vectors)
	dv[i], dv_norm[i]= get_delta(vivecs_smlmean[i],vivecs_lgmean[i])  # dB & dB/|B| (vectors)
	dBmag[i],dBmag_norm[i] = get_deltascalar(Bmag_smlmean[i],Bmag_lgmean[i]) # 
	dvmag[i],dvmag_norm[i] = get_deltascalar(vmag_smlmean[i],vmag_lgmean[i]) # 
	dn[i],dn_norm[i] = get_deltascalar(n_smlmean[i],n_lgmean[i])

# Store delta variables
store_data('dB', data = {'x':timeax,'y':1e9*dB}) # converted to nT
store_data('dB_norm', data = {'x':timeax,'y':dB_norm}) 

store_data('dBmag', data = {'x':timeax,'y':1e9*dBmag}) # converted to nT
store_data('dBmag_norm', data = {'x':timeax,'y':dBmag_norm}) 

store_data('dv', data = {'x':timeax,'y':1e-3*dv}) 
store_data('dv_norm', data = {'x':timeax,'y':dv_norm}) 

store_data('dvmag', data = {'x':timeax,'y':1e-3*dvmag}) 
store_data('dvmag_norm', data = {'x':timeax,'y':dvmag_norm}) 

store_data('dn', data = {'x':timeax,'y':1e-6*dn}) # converted to cm^-3
store_data('dn_norm', data = {'x':timeax,'y':dn_norm}) 

# Example plots showing variable, its average, and the normalized delta values
pyspedas.options('dB_norm','y_range',[-0.5,0.5])
tplot(['Bmag_inst','Bmag','Bmag_mean','dBmag_norm'])
# tplot(['Bvecs_inst','Bvecs','Bvecs_mean','dB_norm'])
# tplot(['ni','n_mean','dn_norm'])


# Here I am defining a threshold for the normalized delta variables
filter_threshold = 0.2

# Any data below the threshold is erased, leaving behind a new array with only the large spikes
# The reason for this is we want to 'count' the data that exceed the threshold
dBmag_norm_filtered = filter_deflections(dBmag_norm,filter_threshold)
dn_norm_filtered = filter_deflections(dn_norm,filter_threshold)

# We are left with a shorter array that contains only the normalized spikes above threshold
print('Length of Full Array: ' + str(len(dBmag_norm)))
print('Length of Filtered Array: ' + str(len(dBmag_norm_filtered)))



#%%
#################### dBmag_count =[] #################
# Here I am setting up a loop to conduct the previous cell's analysis for a series of averaging windows

## BYZ 
# I make a t_windows array, which lists 11 different averaging durations (in minutes)
# 11 averaging durations means there are 10 gaps between them
solarwind = 400

t_windows=np.array([1,27,54,81,108,135,162,189,216,243,270]) # minutes

# d_windows array isn't really meant to be used right now. Just for reference 
# It corresponds to the approx structure size (in Re) corresponding to the same element in the t_windows array
# example: if our large and small windows are 81 and 54 minutes, then we are (roughly) selecting for structures between 200 and 300 Re
d_windows=np.array([4,100,200,300,400,500,600,700,800,900,1000]) # Re

# Count arrays below are generated for source bins centered at these Re values
# source_bins = np.array([50,150,250,350,450,550,650,750,850,950])
source_bins = np.array([4,100,200,300,400,500,600,700,800,900])
# Explicit structure-size bounds to display on bar charts
# span_edges = np.array([0,200,400,600,800,1000,1200])
span_edges = np.array([4,100,200,300,400,500,600,700,800,900,1000])

# First define an empty list which will eventually have 10 elements representing the count between each pair of averaging windows
dBmag_count =[]
dvmag_count =[]
dn_count =[]
# bx y z count ^  
#dBx_count =[]
#dBy_count =[]
#dbz_count  =[]
# Now the loop.  for each element of t_windows, we assign t_windows[i+1] to be the large window and t_windows[i] to be the small one.  
for i in range(len(t_windows)-1):

	n_minutes_lg = t_windows[i+1]
	n_minutes_sml = t_windows[i] 

	# Long Average
	Bvecs_lgmean = get_vecmean(Bvecs,n_minutes_lg*minuteinterval)
	vivecs_lgmean = get_vecmean(vivecs,n_minutes_lg*minuteinterval)
	Bmag_lgmean = get_mean(Bmag,n_minutes_lg*minuteinterval)
	vmag_lgmean = get_mean(vmag,n_minutes_lg*minuteinterval)
	n_lgmean = get_mean(ni,n_minutes_lg*minuteinterval)
	va_lgmean = get_mean(va,n_minutes_lg*minuteinterval)

	# Short Average
	Bvecs_smlmean = get_vecmean(Bvecs,n_minutes_sml*minuteinterval)
	vivecs_smlmean = get_vecmean(vivecs,n_minutes_sml*minuteinterval)
	Bmag_smlmean = get_mean(Bmag,n_minutes_sml*minuteinterval)
	vmag_smlmean = get_mean(vmag,n_minutes_sml*minuteinterval)
	n_smlmean = get_mean(ni,n_minutes_sml*minuteinterval)
	va_smlmean = get_mean(va,n_minutes_sml*minuteinterval)

	# vector
	dB,dB_norm = np.zeros_like(Bvecs),np.zeros_like(Bvecs)
	dBmag,dBmag_norm = np.zeros_like(Bmag),np.zeros_like(Bmag)
	dn,dn_norm = np.zeros_like(ni),np.zeros_like(ni)

	# Get the deltas
	for j in range(len(timeax)):
		dB[j], dB_norm[j]= get_delta(Bvecs_smlmean[j],Bvecs_lgmean[j])  # dB & dB/|B| (vectors)
		dv[j], dv_norm[j]= get_delta(vivecs_smlmean[j],vivecs_lgmean[j])  # dB & dB/|B| (vectors)
		dBmag[j],dBmag_norm[j] = get_deltascalar(Bmag_smlmean[j],Bmag_lgmean[j]) # 
		dvmag[j],dvmag_norm[j] = get_deltascalar(vmag_smlmean[j],vmag_lgmean[j]) # 
		dn[j],dn_norm[j] = get_deltascalar(n_smlmean[j],n_lgmean[j])


	filter_threshold = 0.1

	dBmag_norm_filtered = filter_deflections(dBmag_norm,filter_threshold)
	dvmag_norm_filtered = filter_deflections(dvmag_norm,filter_threshold)
	dbx_norm_filtered  = filter_deflections(dB_norm[:,0],filter_threshold)
	
	dn_norm_filtered = filter_deflections(dn_norm,filter_threshold)
	#dBx_count.append(len(dbx_norm_filtered))
	dBmag_count.append(len(dBmag_norm_filtered))
	dvmag_count.append(len(dvmag_norm_filtered))
	dn_count.append(len(dn_norm_filtered))

	# source_centers = [50,150,250]
	# dBmag_count = [83,42,17]
	# at 50 Re → 83 events
dBmag_span = aggregate_counts_by_bounds(source_bins, dBmag_count, span_edges)
dn_span = aggregate_counts_by_bounds(source_bins, dn_count, span_edges)
dVmag_span = aggregate_counts_by_bounds(source_bins, dvmag_count, span_edges)

span_widths = np.diff(span_edges)
span_labels = get_span_labels(span_edges)
bar_centers = 0.5*(span_edges[:-1] + span_edges[1:])
#span_edges = np.array([0,200,400,600,800,1000])
plt.bar(span_edges[:-1], dBmag_span, width=span_widths, edgecolor='black', align='edge',alpha=0.1,label='Magnetic Magnitude structures',color='b')
plt.bar(span_edges[:-1], dn_span, width=span_widths, edgecolor='black', align='edge',alpha=0.3,label='Density structures',color='k')
plt.bar(span_edges[:-1], dVmag_span, width=span_widths, edgecolor='black', align='edge',alpha=0.4, label = 'Velocity Magnitude structures',color='c')
plt.xticks(bar_centers, span_labels, rotation=45, ha='center')
plt.xlabel(r'Scale Size $(R_e)$ | Filter Threshold = 10%')
plt.ylabel('Count')
plt.yscale('log')
plt.tight_layout()
plt.legend()
plt.show()


#Regular Line plot version of the above histograms

span_centers = 0.5*(span_edges[:-1] + span_edges[1:])
plt.plot(span_centers, dBmag_span,label='Magnetic Magnitude structures' ,alpha=0.1 ,color='b')
plt.plot(span_centers,dn_span,label='Density structures'  ,alpha=0.4 , color='k' )
plt.plot(span_centers,dVmag_span,label='Velocity Magnitude structures' , alpha=0.3, color='c')
plt.xlabel(r'Scale Size $(R_e)$ | Filter Threshold = 10%')
plt.legend()
plt.yscale('log')
# %%

# %%'

# Bvecs = 1e-9 * reform(B_name_new)
# print(Bvecs[0]) ~ 0 , 1 , 2 

# each ind comp of <B> MAG
#%%

## BYZ 
# I make a t_windows array, which lists 11 different averaging durations (in minutes)
# 11 averaging durations means there are 10 gaps between them
solarwind = 400
t_windows=np.array([1,27,54,81,108,135,162,189,216,243,270]) # minutes

# d_windows array isn't really meant to be used right now. Just for reference 
# It corresponds to the approx structure size (in Re) corresponding to the same element in the t_windows array
# example: if our large and small windows are 81 and 54 minutes, then we are (roughly) selecting for structures between 200 and 300 Re
d_windows=np.array([4,100,200,300,400,500,600,700,800,900,1000]) # Re

# Count arrays below are generated for source bins centered at these Re values
# source_bins = np.array([50,150,250,350,450,550,650,750,850,950])
source_bins = np.array([4,100,200,300,400,500,600,700,800,900])
# Explicit structure-size bounds to display on bar charts
# span_edges = np.array([0,200,400,600,800,1000,1200])
span_edges = np.array([4,100,200,300,400,500,600,700,800,900,1000])

dBx_count =[]
dBy_count =[]
dBz_count =[]

for i in range(len(t_windows)-1):

	n_minutes_lg = t_windows[i+1]
	n_minutes_sml = t_windows[i] 

	# Long Average
	Bvecs_lgmean = get_vecmean(Bvecs,n_minutes_lg*minuteinterval)
	vivecs_lgmean = get_vecmean(vivecs,n_minutes_lg*minuteinterval)
	Bmag_lgmean = get_mean(Bmag,n_minutes_lg*minuteinterval)
	vmag_lgmean = get_mean(vmag,n_minutes_lg*minuteinterval)
	n_lgmean = get_mean(ni,n_minutes_lg*minuteinterval)
	va_lgmean = get_mean(va,n_minutes_lg*minuteinterval)

	# Short Average
	Bvecs_smlmean = get_vecmean(Bvecs,n_minutes_sml*minuteinterval)
	vivecs_smlmean = get_vecmean(vivecs,n_minutes_sml*minuteinterval)
	Bmag_smlmean = get_mean(Bmag,n_minutes_sml*minuteinterval)
	vmag_smlmean = get_mean(vmag,n_minutes_sml*minuteinterval)
	n_smlmean = get_mean(ni,n_minutes_sml*minuteinterval)
	va_smlmean = get_mean(va,n_minutes_sml*minuteinterval)

	# vector
	dB,dB_norm = np.zeros_like(Bvecs),np.zeros_like(Bvecs)
	dBmag,dBmag_norm = np.zeros_like(Bmag),np.zeros_like(Bmag)
	dn,dn_norm = np.zeros_like(ni),np.zeros_like(ni)

	# Get the deltas
	for j in range(len(timeax)):
		dB[j], dB_norm[j]= get_delta(Bvecs_smlmean[j],Bvecs_lgmean[j])  # dB & dB/|B| (vectors)
		dv[j], dv_norm[j]= get_delta(vivecs_smlmean[j],vivecs_lgmean[j])  # dB & dB/|B| (vectors)
		dBmag[j],dBmag_norm[j] = get_deltascalar(Bmag_smlmean[j],Bmag_lgmean[j]) # 
		dvmag[j],dvmag_norm[j] = get_deltascalar(vmag_smlmean[j],vmag_lgmean[j]) # 
		dn[j],dn_norm[j] = get_deltascalar(n_smlmean[j],n_lgmean[j])

	filter_threshold = 0.4

	dBz_filtered = filter_deflections(dB_norm[:,2],filter_threshold)
	dby_norm_filtered = filter_deflections(dB_norm[:,1],filter_threshold)
	dbx_norm_filtered  = filter_deflections(dB_norm[:,0],filter_threshold)
	
	dn_norm_filtered = filter_deflections(dn_norm,filter_threshold)
	dBx_count.append(len(dbx_norm_filtered))
	dBy_count.append(len(dby_norm_filtered))
	dBz_count.append(len(dBz_filtered))


dBx_span = aggregate_counts_by_bounds(source_bins, dBx_count, span_edges)
dBy_span = aggregate_counts_by_bounds(source_bins, dBy_count, span_edges)
dBz_span = aggregate_counts_by_bounds(source_bins, dBz_count, span_edges)

#span_edges = np.array([0,200,400,600,800,1000])
span_widths = np.diff(span_edges)
span_labels = get_span_labels(span_edges)
bar_centers = 0.5*(span_edges[:-1] + span_edges[1:])


plt.figure(figsize=(10, 6))
plt.bar(span_edges[:-1], dBx_span, width=span_widths, edgecolor='black', align='edge',
alpha=0.3, label='dBx structures', color='b')
plt.bar(span_edges[:-1], dBy_span, width=span_widths, edgecolor='black',
align='edge', alpha=0.2, label='dBy structures', color='r')
plt.bar(span_edges[:-1], dBz_span, width=span_widths, edgecolor='black',
align='edge', alpha=0.1, label='dBz structures', color='g')
plt.xticks(bar_centers, span_labels, rotation=45, ha='center')
plt.xlabel(r'Scale Size $(R_e)$ |   Filter Threshold = 40%')
plt.ylabel('Count')
plt.yscale('log')
plt.tight_layout()
plt.legend()
plt.show()

## line plt  dBx , x , z counts
span_centers = 0.5*(span_edges[:-1] + span_edges[1:])
plt.plot(span_centers, dBx_span,label='dBx structures' ,color='b' , alpha=0.4)
plt.plot(span_centers,dBy_span,label='dBy structures', color='r' , alpha=0.4)
plt.plot(span_centers,dBz_span,label='dBz structures'  , color='g' , alpha=0.4)
plt.xlabel(r'Scale Size $(R_e)$ |   Filter Threshold = 40%')
plt.legend()
plt.yscale('log')
# %%

################ this for the vX Y Z ################################
# I make a t_windows array, which lists 11 different averaging durations (in minutes)
# 11 averaging durations means there are 10 gaps between them
t_windows=np.array([1,27,54,81,108,135,162,189,216,243,270]) # minutes

# d_windows array isn't really meant to be used right now. Just for reference 
# It corresponds to the approx structure size (in Re) corresponding to the same element in the t_windows array
# example: if our large and small windows are 81 and 54 minutes, then we are (roughly) selecting for structures between 200 and 300 Re
d_windows=np.array([4,100,200,300,400,500,600,700,800,900,1000]) # Re

# Count arrays below are generated for source bins centered at these Re values
source_bins = np.array([50,150,250,350,450,550,650,750,850,950])

# Explicit structure-size bounds to display on bar charts
span_edges = np.array([4,100,200,300,400,500,600,700,800,900,1000])


# Vx = vivecs[:, 0]
# Vy = vivecs[:, 1]
# Vz = vivecs[:, 2]
dVx_count =[]
dVy_count =[]
dVz_count =[]

for i in range(len(t_windows)-1):

	n_minutes_lg = t_windows[i+1]
	n_minutes_sml = t_windows[i] 

	# Long Average
	Bvecs_lgmean = get_vecmean(Bvecs,n_minutes_lg*minuteinterval)
	vivecs_lgmean = get_vecmean(vivecs,n_minutes_lg*minuteinterval)
	Bmag_lgmean = get_mean(Bmag,n_minutes_lg*minuteinterval)
	vmag_lgmean = get_mean(vmag,n_minutes_lg*minuteinterval)
	n_lgmean = get_mean(ni,n_minutes_lg*minuteinterval)
	va_lgmean = get_mean(va,n_minutes_lg*minuteinterval)

	# Short Average
	Bvecs_smlmean = get_vecmean(Bvecs,n_minutes_sml*minuteinterval)
	vivecs_smlmean = get_vecmean(vivecs,n_minutes_sml*minuteinterval)
	Bmag_smlmean = get_mean(Bmag,n_minutes_sml*minuteinterval)
	vmag_smlmean = get_mean(vmag,n_minutes_sml*minuteinterval)
	n_smlmean = get_mean(ni,n_minutes_sml*minuteinterval)
	va_smlmean = get_mean(va,n_minutes_sml*minuteinterval)

	# vector
	dB,dB_norm = np.zeros_like(Bvecs),np.zeros_like(Bvecs)
	dBmag,dBmag_norm = np.zeros_like(Bmag),np.zeros_like(Bmag)
	dn,dn_norm = np.zeros_like(ni),np.zeros_like(ni)

	# Get the deltas
	for j in range(len(timeax)):
		dB[j], dB_norm[j]= get_delta(Bvecs_smlmean[j],Bvecs_lgmean[j])  # dB & dB/|B| (vectors)
		dv[j], dv_norm[j]= get_delta(vivecs_smlmean[j],vivecs_lgmean[j])  # dB & dB/|B| (vectors)
		dBmag[j],dBmag_norm[j] = get_deltascalar(Bmag_smlmean[j],Bmag_lgmean[j]) # 
		dvmag[j],dvmag_norm[j] = get_deltascalar(vmag_smlmean[j],vmag_lgmean[j]) # 
		dn[j],dn_norm[j] = get_deltascalar(n_smlmean[j],n_lgmean[j])

	filter_threshold = 0.05

	dvz_filtered = filter_deflections(dv_norm[:,2],filter_threshold)
	dvy_norm_filtered = filter_deflections(dv_norm[:,1],filter_threshold)
	dvx_norm_filtered  = filter_deflections(dv_norm[:,0],filter_threshold)
	
	dn_norm_filtered = filter_deflections(dn_norm,filter_threshold)
	dVx_count.append(len(dvx_norm_filtered))
	dVy_count.append(len(dvy_norm_filtered))
	dVz_count.append(len(dvz_filtered))


dVx_span = aggregate_counts_by_bounds(source_bins, dVx_count, span_edges)
dVy_span = aggregate_counts_by_bounds(source_bins, dVy_count, span_edges)
dVz_span = aggregate_counts_by_bounds(source_bins, dVz_count, span_edges)

span_widths = np.diff(span_edges)
span_labels = get_span_labels(span_edges)
bar_centers = 0.5*(span_edges[:-1] + span_edges[1:])

plt.figure(figsize=(10, 6))
plt.bar(span_edges[:-1], dVx_span, width=span_widths, edgecolor='black', align='edge',
alpha=0.5, label='dVx structures', color='b')
plt.bar(span_edges[:-1], dVy_span, width=span_widths, edgecolor='black', align='edge',
alpha=0.3, label='dVy structures', color='r')
plt.bar(span_edges[:-1], dVz_span, width=span_widths, edgecolor='black',
align='edge', alpha=0.1, label='dVz structures', color='g')
plt.xticks(bar_centers, span_labels, rotation=45, ha='center')
# hard coded in swich thresholds each run 
plt.xlabel(r'Scale Size $(R_e)$ |   Filter Threshold = 5%')
plt.ylabel('Count')
plt.yscale('log')
plt.tight_layout()
plt.legend()
plt.show()

span_centers = 0.5*(span_edges[:-1] + span_edges[1:])
plt.plot(span_centers, dVx_span,label='dVx structures' ,color='b')
plt.plot(span_centers,dVy_span,label='dVy structures', color='k')
plt.plot(span_centers,dVz_span,label='dVz structures'  , color='c' )
plt.xlabel(r'Scale Size $(R_e)$ |   Filter Threshold = 5%')
plt.ylabel('Count')
plt.legend()
plt.yscale('log')

# %%


# %%
Vx = vivecs[:, 0]
Vy = vivecs[:, 1]
Vz = vivecs[:, 2]
import matplotlib.pyplot as plt
plt.hist(Vx, bins=10, alpha=0.4,color = "red" , label="Vx" )
plt.hist(Vy, bins=10, alpha=0.4, color = "blue" ,label="Vy")
plt.hist(Vz, bins=10, alpha=0.4, color = "green" ,label="Vz")
plt.legend()
plt.show()
# %%

# b = Bdata    # b.times, b.y
# v = vdata

# B_all = b.y[:]         # full data array
# t_all = b.times[:]   

# plt.hist(B_all, bins=10,alpha=10,color="blue" , label= "Practice")
# plt.hist(B_all, bins=10,alpha=10,color="green" , label= "Practice") 

#%%

# %%
import numpy as np


# Example variable names
B_name = 'BGSE'
ni_name = 'N_elec'
vi_name = 'U_eGSE'
pos_name = 'DIST'

# Put all variables in a list
var_list = [B_name, ni_name, vi_name, pos_name]

print("Sampling rate / temporal resolution for each variable:\n")

for var_name in var_list:
    test = get_data(var_name)

    time = test.times
    total_measurements = len(time)

    # Need at least 2 time points
    if total_measurements < 2:
        print(f"{var_name} -> not enough data points")
        continue

				# last and first item in time
    total_time_seconds = time[-1] - time[0]

    # Avoid division by zero
    if total_time_seconds == 0:
        print(f"{var_name} -> total time is zero")
        continue

    # Measurements per second
    sensitivity = total_measurements / total_time_seconds

    # Seconds per measurement
    delta_t = 1 / sensitivity

    # Measurements per minute
    measurements_per_min = sensitivity * 60

    print(f"Variable: {var_name}")
    print(f"  Total measurements: {total_measurements}")
    print(f"  Total time: {total_time_seconds:.2f} s")
    print(f"  Sensitivity: {sensitivity:.6f} measurements/s")
    print(f"  Time between measurements: {delta_t:.2f} s")
    print(f"  Measurements per minute: {measurements_per_min:.2f}")
    print("-" * 50)
#https://numpy.org/doc/stable/user/basics.indexing

# %%
# %%

# 
 # %%
# 

# Fixed comparison windows (do not reuse loop vars)
cmp_short_min = 5
cmp_long_min  = 240

Bmag_short_cmp = get_mean(Bmag, cmp_short_min * minuteinterval)
Bmag_long_cmp  = get_mean(Bmag, cmp_long_min  * minuteinterval)

short_Bmag_nt = 1e9 * Bmag_short_cmp
long_Bmag_nt  = 1e9 * Bmag_long_cmp
t_days = (timeax - timeax[0]) / 86400.0

# normalized difference: (short - long) / long
eps = 1e-30
dBmag_norm_cmp = (Bmag_short_cmp - Bmag_long_cmp) / np.where(np.abs(Bmag_long_cmp) > eps, Bmag_long_cmp, np.nan)

# same y-scale for top two panels
ymin = float(np.nanmin([np.nanmin(short_Bmag_nt), np.nanmin(long_Bmag_nt)]))
ymax = float(np.nanmax([np.nanmax(short_Bmag_nt), np.nanmax(long_Bmag_nt)]))

fig, ax = plt.subplots(3, 1, figsize=(12, 8), sharex=True)

ax[0].plot(t_days, short_Bmag_nt, 'k-', lw=0.8)
ax[0].set_title(f'Short avg ({cmp_short_min} min)')
ax[0].set_ylabel('Bmag (nT)')
ax[0].set_ylim(ymin, ymax)
ax[0].grid(alpha=0.25)

ax[1].plot(t_days, long_Bmag_nt, 'k-', lw=0.8)
ax[1].set_title(f'Long avg ({cmp_long_min} min)')
ax[1].set_ylabel('Bmag (nT)')
ax[1].set_ylim(ymin, ymax)
ax[1].grid(alpha=0.25)

ax[2].plot(t_days, dBmag_norm_cmp, color='tab:blue', lw=0.8)
ax[2].axhline(0.0, color='k', lw=0.8, alpha=0.7)
ax[2].set_title('Normalized difference: (short - long) / long')
ax[2].set_ylabel('dBmag_norm')
ax[2].set_xlabel('Days from start')
ax[2].grid(alpha=0.25)

plt.tight_layout()
plt.show()

#existing code
# %%
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
		var_small = get_mean(scalvar, T_small)
		var_large = get_mean(scalvar, T_large)

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
		var_small = get_vecmean(vecvar, T_small)
		var_large = get_vecmean(vecvar, T_large)

		delta_var = (var_small - var_large) / var_large
		delta_var_array.append(delta_var)

	return np.abs(delta_var_array)	



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

# Km_sw_test = Vx[:, 0]

Km_sw_test = 400

sw_km_to_Re(Km_sw_test)
#sw_km_to_Re(Vx)

# print( "Result hard coded 400 R_E :" , sw_km_to_Re(Km_sw_test))
# print( "Result Vx ~ in the x direction:" , np.abs(sw_km_to_Re(Vx))) 

#############
t_windows_test =np.array([])

t_windows_test = t_windows=np.array([1,27,54,81,108,135,162,189,216,243,270])

print("Result of Vsw test: " ,Vsw_Tw(t_windows_test))


#print("Result of the Vector helper function : " , get_vecdeltas(t_windows_test,vivecs))

#print("Result of the ni helper function : " , get_nideltas(t_windows_test,ni))

a = get_scaldeltas(t_windows_test, Bmag)
b = get_nideltas(t_windows_test, ni)
Bscal = a
dni = b
#creates a bool cond Bscal is greater than (random int)  then its in the mask 
mask = Bscal > 0.1
#puts true into the mask if cond is meet 
Bscal[mask]
# Show dni values at the SAME locations where Bscal > 5
dni[mask]
mask = Bscal >= 0.1
Bscal[mask]
dni[mask]
mask = dni < 0.1
Bscal[mask]
print("Mask:", mask)
print("Filtered a:", Bscal[mask])
print("Filtered b:", dni[mask])

new_Flitered_Bscal = Bscal[mask]
new_Flitered_dni = dni[mask]

#-----------------------------------------#



#flatten()
# plt.hist(np.ravel(Bscal), bins=10, alpha=0.4, color='red', label='Bmag')
# plt.hist(np.ravel(dni), bins=10, alpha=0.4, color='blue', label='ni')
# plt.legend()
# plt.show()

# range is = 0 ≤ x ≤ 0.1
# plt.hist(np.ravel(Flitered_Bscal), bins=100,
#          alpha=0.4, color='red', label='Bmag')
 
# plt.hist(np.ravel(Flitered_dni), bins=100, 
#          alpha=0.4, color='blue', label='ni')
# plt.ylabel("Counts of ni and Bmag")
# plt.xlabel("Normlized Bmag and ni at 10%")
# plt.legend()
# plt.show()
# for i in range(len(Bscal)):
# 	print(i, Bscal[i], dni[i], mask[i])
#flatten()
# .shape 
# np.ravel converts a 
# multi-dimensional array 
# into a one-dimensional (flattened) array
# arr = np.array([[1, 2, 3], [4, 5, 6]]) -> Output: [1, 2, 3, 4, 5, 6]
# Bvecs = 1e-9*reform(B_name_new) # Converting from nanotesla to Tesla
# vivecs = 1e3*reform(vi_name_new) #converting km/s to m/s
# ni = 1e6*reform(ni_name_new)

#%%
# scatter plt of Bscal and ni
x = np.ravel(new_Flitered_Bscal)

y = np.ravel(new_Flitered_dni)

plt.scatter(y,x, color = "#FEB88D" , label = "$\Delta n_i$" , alpha= 1 )
plt.scatter(x,y, color = "#6650F8" , label = "$\Delta scal(vx)$" ,alpha=0.7 )

plt.xlabel("$\Delta n_i / n_i$ Proton Density Relative Fluctuation")
plt.ylabel("$\Delta Bscal$ Mag field fluxations")
plt.title("Relationship between dni and Bscal at 10% threshold")
plt.legend(loc = "upper right")

plt.show()




# %%


# vivecs = reform(get_data('vi')) 

# # 2. Pass the full 3D array into the vector function
# velocity_deltas = get_vecdeltas(t_windows_test, vivecs)

# print("Result of the Vector helper function: ", velocity_deltas)

# # 3. Extract your individual components AFTER the calculation is done
# delta_Vx = velocity_deltas[:, 0]
# # delta_Vy = velocity_deltas[:, 1]
# # delta_Vz = velocity_deltas[:, 2]

# # #d_Vector = np.ravel(velocity_deltas)
# # d_Vector = velocity_deltas
# # #creates a bool cond Bscal is greater than (random int)  then its in the mask 
# # mask_vec = d_Vector > 0.1
# #puts true into the mask if cond is meet 
# d_Vector[mask_vec]
# # Show dni values at the SAME locations where Bscal > 5

# mask_vec = d_Vector >= 0.1
# d_Vector[mask_vec]

# print("Filtered Components:", d_Vector[mask_vec])

# new_Filtered_vector = d_Vector[mask_vec]

# x = np.ravel(new_Filtered_vector)
# y = np.ravel(new_Flitered_dni)
# plt.scatter(x,y, color = "#FEB88D" , label = "$\Delta n_i$" , alpha= 1 )
# plt.scatter(x,y, color = "#6650F8" , label = "$\Delta vector$" ,alpha=0.7 )

# plt.xlabel("$\Delta n_i / n_i$ Proton Density Relative Fluctuation")
# plt.ylabel("$\Delta vector$ V flux")
# plt.title("Relationship between dni and vector at 10% threshold")
# plt.legend(loc = "upper right")

# # plt.show()
# %%
# # 1. FIX: Uncomment this line to actually load the 3D data array
# vivecs = reform(get_data('vi')) 

# # 2. Pass the full 3D array into the vector function
# velocity_deltas = get_vecdeltas(t_windows_test, vivecs)

# print("Result of the Vector helper function: ", velocity_deltas)

# # 3. Extract your individual components AFTER the calculation is done
# delta_Vx = velocity_deltas[:, 0]
# delta_Vy = velocity_deltas[:, 1]
# delta_Vz = velocity_deltas[:, 2]
# # %%

# %%
