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

def duration_sec(trange):
	from datetime import datetime as dt
	start = dt.strptime(trange[0], '%Y-%m-%d/%H:%M')
	stop = dt.strptime(trange[1], '%Y-%m-%d/%H:%M')
	duration = stop-start
	duration_s = duration.total_seconds()
	duration_m = duration_s/60 #minutes
	duration_h = duration_m/60 #hours
	return duration_s


# takes an array of some length and a timerange in [start,stop] format and determines how many elements of the array corresponds to 1 minute 
# (so an array of length 120 representing an hour of data would return a value of 2, because 2 data points of that array represent a minute)
def minute_int(timeax,trange):
	steps = len(timeax)
	minutes = duration(trange)
	minute = np.ceil(steps/minutes)
	return int(minute)

def second_int(timeax,trange):
	steps = len(timeax)
	seconds = duration_sec(trange)
	second = np.ceil(steps/seconds)
	return int(second)

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
		var_small = get_mean(scalvar, T_small)
		var_large = get_mean(scalvar, T_large)
		#var_small = get_mean(scalvar, T_small*minuteinterval)
		#var_large = get_mean(scalvar,T_large*minuteinterval)

		delta_var = (var_small - var_large) / var_large
		delta_var_array.append(delta_var)

	return np.abs(np.array(delta_var_array))

def get_scaldeltas (t_windows,scalvar):
	#vx = V[:, 0]
	delta_var_array = []
	# protect against going out of bounds from the array 
	minuteinterval = minute_int(timeax,trange)
	secondinterval = second_int(timeax,trange)
	for i in range(len(t_windows) - 1):

		# These are the time windows <B>_ts - <B>_tL / <B>_tL 
		T_small = t_windows[i]
		T_large = t_windows[i+1]
		# Compute moving-average magnetic field values using
		# the small and large time windows
		var_small = get_mean(scalvar, T_small*minuteinterval)
		var_large = get_mean(scalvar,T_large*minuteinterval)

		delta_var = (var_small - var_large) / var_large
		delta_var_array.append(delta_var)

	return np.abs(np.array(delta_var_array))


def get_vecdeltas (t_windows, vecvar):
	#vx = V[:, 0]
	delta_var_array = []
	# protect against going out of bounds from the array
	minuteinterval = minute_int(timeax,trange)
	secondinterval = second_int(timeax,trange)

	for i in range(len(t_windows) - 1):
		# These are the time windows <B>_ts - <B>_tL / <B>_tL 
		T_small = t_windows[i]
		T_large = t_windows[i+1]
		# Compute moving-average magnetic field values using
		# the small and large time windows
		var_small = get_vecmean(vecvar, T_small*minuteinterval)  
		var_large = get_vecmean(vecvar, T_large*minuteinterval) 

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

#%%

trange=['2019-01-01/00:00', '2019-01-30/00:00']
# trange=['2019-01-20/00:00', '2019-01-30/00:00']
swe_vars = pyspedas.projects.wind.swe(trange=trange)
# denstiy 
mfi_vars = pyspedas.projects.wind.mfi(trange=trange)
#%%
B_name='BGSE'
ni_name='N_elec'
vi_name = 'U_eGSE'
pos_name = 'DIST'
 # Distance from earth, in units of earth radii.  1 Re ~ 4.3e-5 AU
#B_name_new = get_data('B')
#ni_name_new = get_data('ni')
#vi_name_new = get_data('vi')

B_name_new = get_data(B_name)
ni_name_new = get_data(ni_name)
vi_name_new = get_data(vi_name)
#%%
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
interpvar_name = vi_name

# I like to define a timeax variable as well, representing the time variable
interpdata = get_data(interpvar_name)
timeax = interpdata.times


# To interpolate a tplot variable to another tplot variable, we use pyspedas's tinterpol and assign new names to the interpolated versions
tinterpol(B_name,interpvar_name,newname='B')
tinterpol(ni_name,interpvar_name,newname='ni') # (Obviously, interpolating a variable to itself isn't necessary, but I keep it here because sometimes I want to change the interpvar)
tinterpol(vi_name,interpvar_name,newname='vi') # 

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
Bscal = get_scaldeltas(t_windows_test, Bmag)
dni = get_nideltas(t_windows_test, ni)
Velocity_mag = get_scaldeltas(t_windows_test, vmag)

#creates a bool cond Bscal is greater than (random int)  then its in the mask 
mask = Bscal > 0.4
mask_vmag = Velocity_mag > 0.01
#puts true into the mask if cond is meet 
Bscal[mask]
Velocity_mag[mask_vmag]
# Show dni values at the SAME locations where Bscal > 5
dni[mask]
mask = Bscal >= 0.4
mask_vmag = Velocity_mag >= 0.4
Bscal[mask]
dni[mask]
mask = dni < 0.4
Bscal[mask]

print("Mask for Bscal:", mask)
print("Flitered c" , Velocity_mag[mask_vmag])
print("Filtered a:", Bscal[mask])
print("Filtered b:", dni[mask])




#Choose which timescale range to look at
Bscal_0 = Bscal[0]  # Changed to 0
dni_0 = dni[0]      # Changed to 0
vmag_0 = Velocity_mag[0] # Changed to 0

Bscal_1 = Bscal[1]
dni_1 = dni[1]
vmag_1 = Velocity_mag[1]

# 3. Make Clean, Consistent Masks
mask0 = (Bscal_0 >= 0.05) & (dni_0 >= 0.05)
mask1 = (Bscal_1 >= 0.05) & (dni_1 >= 0.05)

# Fixed typo: changed dni_1 to dni_0 to match timescale 0
mask_vmag_0 = (vmag_0 >= 0.4) & (dni_0 >= 0.4) 
mask_vmag_1 = (vmag_1 >= 0.4) & (dni_1 >= 0.4)

# 4. Filter Arrays
new_Flitered_Bscal0 = Bscal_0[mask0]
new_Flitered_dni0 = dni_0[mask0]

new_Flitered_Bscal1 = Bscal_1[mask1]
new_Flitered_dni1 = dni_1[mask1]

new_Filtered_vector0 = vmag_0[mask_vmag_0]
new_Filtered_vector1 = vmag_1[mask_vmag_1]

#Plotting 
plt.scatter(new_Flitered_Bscal0,new_Flitered_dni0,s=10,color='b' )
plt.scatter(new_Flitered_Bscal1,new_Flitered_dni1,s=10,color='k')
plt.plot(Bscal_0,Bscal_0,color='r', label= 'Ref line')
plt.xlabel(r'$\delta B/B$')
plt.ylabel(r'$\delta n/n$')
#plt.legend(loc = "upper center", frontsize = "small") 
plt.ylim(0,1)
plt.xlim(0,1)

#-----------------------------------------#


# into a one-dimensional (flattened) array
# arr = np.array([[1, 2, 3], [4, 5, 6]]) -> Output: [1, 2, 3, 4, 5, 6]

# Bvecs = 1e-9*reform(B_name_new) # Converting from nanotesla to Tesla
# vivecs = 1e3*reform(vi_name_new) #converting km/s to m/s
# ni = 1e6*reform(ni_name_new)

#%%
#Get data

Bscal = get_scaldeltas(t_windows_test, Bmag)
dni = get_nideltas(t_windows_test, ni)
Velocity_mag = get_scaldeltas(t_windows_test, vmag)

#isolate separate masks
mask_B = Bscal >= 0.05
mask_vmag = Velocity_mag >= 0.05
mask_dni = dni >= 0.05

#Testing the masking 
# print("Filtered Bscal :", Bscal[mask_B])
# print("Filtered vmag:", Velocity_mag[mask_vmag])
# print("Filtered dni :", dni[mask_dni])
Bscal_0 = Bscal[0]
Bscal_1 = Bscal[1]

dni_0 = dni[0]
dni_1 = dni[1]
dni_2 = dni[2]

vmag_0 = Velocity_mag[0]
vmag_1 = Velocity_mag[1]

#Boolean Masks 
mask0 = (Bscal_0 >= 0.1) & (dni_0 >= 0.1)
mask1 = (Bscal_1 >= 0.1) & (dni_1 >= 0.1)
#Vmag is the magnitude (length) of the entire velocity vector.
Vmag_mask0 = (vmag_0 >= 0.1) & (dni_0 >= 0.1)
Vmag_mask1 = (vmag_1 >= 0.1) & (dni_1 >= 0.1)

# Apply masks to filter data
new_Filtered_Bscal0 = Bscal_0[mask0]
new_Filtered_dni0 = dni_0[mask0]

new_Filtered_Bscal1 = Bscal_1[mask1]
new_Filtered_dni1 = dni_1[mask1]

new_Filtered_vmag0 = vmag_0[Vmag_mask0]
new_Filtered_dni_vmag0 = dni_0[Vmag_mask0] 

new_Filtered_vmag1 = vmag_1[Vmag_mask1]
new_Filtered_dni_vmag1 = dni_1[Vmag_mask1] 

#Plotting
plt.scatter(new_Filtered_Bscal0, new_Filtered_dni0, s=10, color='b', label='Timescale 0')
plt.scatter(new_Filtered_Bscal1, new_Filtered_dni1, s=10, color='k', label='Timescale 1')
plt.plot(Bscal_0, Bscal_0, color='r', label='Reference Line')

plt.xlabel(r'$\delta B/B$')
plt.ylabel(r'$\delta n/n$')
plt.ylim(0, 1)
plt.xlim(0, 1)

# Fixed typo from 'frontsize' to 'fontsize'
plt.legend(loc="upper center", fontsize="small") 
plt.show()

plt.scatter(new_Filtered_vmag0, new_Filtered_dni_vmag0, s=10, color='b', label='Timescale 0')
plt.scatter(new_Filtered_vmag1, new_Filtered_dni_vmag1, s=10, color='k', label='Timescale 1')
plt.plot(Bscal_0, Bscal_0, color='r', label='Reference Line')

plt.xlabel(r'$\delta Bvmag/Bvmag$')
plt.ylabel(r'$\delta n/n$')
plt.ylim(0, 1)
plt.xlim(0, 1)

# Fixed typo from 'frontsize' to 'fontsize'
plt.legend(loc="upper center", fontsize="small") 
plt.show()
# %%


vivecs = reform(get_data('vi')) 

# 2. Pass the full 3D array into the vector function
velocity_deltas = get_vecdeltas(t_windows_test, vivecs)
density_deltas = get_nideltas(t_windows_test, ni)

print("Result of the Vector helper function: ", velocity_deltas)

dvx_0 = velocity_deltas[0,:,0]  #Timescale 0, time_samples  , X-component
dvy_0 = velocity_deltas[0,:,1]  
dvz_0 = velocity_deltas[0,:,2]  

dni_0 = density_deltas[0,:]


mask0x = (dvx_0 >= 0.1) & (dni_0 >= 0.1)


mask = mask0x
# mask_new0y = mask0y
# mask_new0z = mask0z

new_Filtered_dvx_0 = dvx_0[mask]
new_Filtered_dvy_0 = dvy_0[mask]

new_Filtered_dvz_0 = dvz_0[mask]
new_Filtered_dni_0 = dni_0[mask]

total_points_0 =  new_Filtered_dvx_0.shape[:]

# velocity comp in timestep 1 
dvx_1 = velocity_deltas[1,:,0]  #Timescale 1, X-component
dvy_1 = velocity_deltas[1,:,1]  
dvz_1 = velocity_deltas[1,:,2] 

dni_1 = density_deltas[1,:]

mask1x = (dvx_1 >= 0.1) & (dni_1 >= 0.1)

mask1 = mask1x

new_Filtered_dvx_1 = dvx_1[mask1]
new_Filtered_dvy_1 = dvy_1[mask1]
new_Filtered_dvz_1 = dvz_1[mask1]
new_Filtered_dni_1 = dni_1[mask1]
total_points_1 =  new_Filtered_dvx_1.shape[:]


# axis ~(x,y)



plt.scatter(new_Filtered_dvx_0, new_Filtered_dni_0, 
            s=10, color='g', label=f'Timescale 0 ~ Num={total_points_0}')
plt.scatter(new_Filtered_dvx_1, new_Filtered_dni_1, 
			s=10, color='b', label=f'Timescale 1 ~ Num={total_points_1}')
plt.plot(dvx_0 ,dvx_0,  color='r', label='Reference Line')
plt.xlabel(r'$\delta v_x / v$     Filter: 10%')
plt.ylabel(r'$\delta ni / ni$')
# linear
plt.yscale('log')
plt.xscale('log')
plt.ylim(0.1, 3)
plt.xlim(0.1, 3)
plt.legend()




# scales sizes in the v
#np.sqrt(dvx**2 + dvy**2 + dvz**2)
#%%
plt.scatter(new_Filtered_dvy_0, new_Filtered_dni_0, 
            s=10, color='g', label=f'Timescale 0 ~ Num={total_points_0}')
plt.scatter(new_Filtered_dvy_1, new_Filtered_dni_1, 
			s=10, color='b', label=f'Timescale 1 ~ Num={total_points_1}')
plt.plot(dvy_0 ,dvy_0,  color='r', label='Reference Line')
plt.xlabel(r'$\delta v_y / v$     Filter: 10%')
plt.ylabel(r'$\delta ni / ni$')

# linear
plt.yscale('log')
plt.xscale('log')
plt.ylim(0.1, 3)
plt.xlim(0.1, 3)
plt.legend()

#%%
plt.scatter(new_Filtered_dvz_0, new_Filtered_dni_0, 
            s=10, color='g', label=f'Timescale 0 ~ Num={total_points_0}')
plt.scatter(new_Filtered_dvz_1, new_Filtered_dni_1, 
			s=10, color='b', label=f'Timescale 1 ~ Num={total_points_1}')
plt.plot(dvz_0 ,dvz_0,  color='r', label='Reference Line')
plt.xlabel(r'$\delta v_z / v$     Filter: 10%')
plt.ylabel(r'$\delta ni / ni$')
# linear
plt.yscale('log')
plt.xscale('log')
plt.ylim(0.1, 3)
plt.xlim(0.1, 3)
plt.legend()


# %%

TimeArr = t_windows=np.array([1,27,54,81,108,135,162,189,216,243,270])
#

Vx = vivecs[:, 0]
Vy = vivecs[:, 1]
Vz = vivecs[:, 2]

plt.hist(Vx, bins=10, alpha=0.4,color = "red" , label="Vx" )
plt.hist(Vy, bins=10, alpha=0.4, color = "blue" ,label="Vy")
plt.hist(Vz, bins=10, alpha=0.4, color = "green" ,label="Vz")
plt.legend()
plt.show()

SolarWind = abs(Vx) 
ThosandRE_Km = 6000000 # 1,000 Re = 6mill 
min = 60

RE_sec = ThosandRE_Km / SolarWind 
print(RE_sec)

Min_to_hours = abs(RE_sec) / min

print(Min_to_hours)


# %%
#Einstein Summation np.einsum
# def get_vector_angle(vecA, vecB):
#     """
#     Calculates the angle (in degrees) between two vector time series.
#     vecA and vecB should be NumPy arrays of shape (N, 3).
#     """
#     #  Compute dot product row by row
#     dot_product = np.einsum('ij,ij->i', vecA, vecB)
    
#    # Compute magnitudes
#     norm_A = np.linalg.norm(vecA, axis=1)
#     # 3. Compute the angle in degrees
#     norm_B = np.linalg.norm(vecB, axis=1)
    
#     # Avoid division by zero
#     with np.errstate(invalid='ignore', divide='ignore'):
#         cos_theta = dot_product / (norm_A * norm_B)
#         cos_theta = np.clip(cos_theta, -1.0, 1.0)
        
#     angles_deg = np.degrees(np.arccos(cos_theta))
#     return angles_deg

# #	psp code for ref 

def get_angle(vec,meanvec):
	term1 = np.dot(vec,meanvec)
	term2 = np.dot(np.linalg.norm(vec),np.linalg.norm(meanvec))
	term3 = term1/term2
	term4 = np.arccos(term3)*180/np.pi
	return term4

TimeArr = t_windows=np.array([1,27,54,81,108,135,162,189,216,243,270])

Bscal = get_scaldeltas(TimeArr, Bmag)
dni = get_nideltas(TimeArr, ni)
Velocity_mag = get_scaldeltas(TimeArr, vmag)

dni_mask = dni >= 0.05
B_mask = Bscal>= 0.05
bulk_mask_v = Velocity_mag >= 0.05

dni[dni_mask]
Velocity_mag[bulk_mask_v]
Bscal[B_mask]
print("dni: ", dni[dni_mask])

print("Filtered Bscal :", Bscal[B_mask])

fl_B = Bscal[B_mask]
fl_dni = dni[dni_mask]
fl_Bulk_v = Velocity_mag[bulk_mask_v]


#use the dot product formula
# custom_bins = [0.05,1, 10, 100, 1000]
plt.figure(figsize=(10, 6))
plt.hist(fl_B, bins=10, alpha=0.4,color = "dodgerblue" , label="Mag" )
plt.hist(fl_dni, bins=10, alpha=0.4,color = "coral" , label="dni" )
plt.hist(fl_Bulk_v, bins=10, alpha=0.4,color = "pink" , label="Velocity" )
#plt.hist(theta_B_zero, bins=10, alpha=0.4, color="purple", label="Theta (B-V Angle)")
plt.yscale('log')
plt.xscale('log')
#matplotlib.colors.CSS4_COLORS
#bounds using powers of 10
plt.xlabel("Fractional Change ($\delta$ = |$\Delta$Var| / Var$_{Large}$)\n "
'(Filtered ~ |$\Delta$Desnity| , |$\Delta$Mag|, |$\Delta$Velocity| > 5%) ', fontsize=16)
plt.ylabel('Data Point Count (Log Scale)', fontsize=16)
bin_labels = ["1-10", "10-100", "100-1000"]
plt.title('Distribution of Mesoscale Fluctuations at 1 AU', fontsize=16, fontweight='bold')
plt.xlim(0.05, 1000.0)
plt.ylim(0.1, 2000.0)
plt.xticks([1, 10, 100, 1000], ["1", "10", "100", "1000"])
plt.legend()
plt.show()
# %%
#trange=['2019-01-01/00:00', '2019-01-30/00:00']

# B_name_new = get_data(B_name)
# ni_name_new = get_data(ni_name)
# vi_name_new = get_data(vi_name)


# for i in range len(vi_name_new):

# 	theta_vb = get_angle(vi_name_new[i],B_name_new[i])
	



# TimeArr = t_windows=np.array([1,27,54,81,108,135,162,189,216,243,270])
# get_mean(vi_name_new)
# theta_V = get_angle(vi_name_new,)
# theta_V_zero = theta_V[0] 

# theta_B = get_vector_angle(t_windows,B_name_new)
# theta_B_zero = theta_B[0] 

# print(" Theta V after finding theta:", theta_V)
# print("zero :", theta_V_zero)
# print(" Theta B after finding theta:", theta_B)
# print("zero :", theta_B_zero)


# %%
# 1. Define your exact bin edges


# 2. Draw the histograms using these custom bins

# %%
def aggregate_counts_by_bounds(source_centers, source_counts, target_edges):


	
    aggregated = []
    for i in range(len(target_edges) - 1):
        left = target_edges[i]
        right = target_edges[i + 1]

        mask = (source_centers >= left) & (source_centers < right)

       
        total = np.sum(np.array(source_counts)[mask])

        # 2. Check if the sum is NaN or empty; if it is, make it 0
        if np.isnan(total) or mask.sum() == 0:
            aggregated.append(0)
        else:
            aggregated.append(int(total))  

    return np.array(aggregated)
	
source_bins = np.array([0.05,1,10,100,1000])
span_edges = np.array([0.05,1,10,100,1000])
source_centers = np.sqrt(source_bins[:-1] * source_bins[1:])



TimeArr = t_windows=np.array([1,27,54,81,108,135,162,189,216,243,270])

Bscal = get_scaldeltas(TimeArr, Bmag)
dni = get_nideltas(TimeArr, ni)
Velocity_mag = get_scaldeltas(TimeArr, vmag)

dni_mask = dni >= 0.05
B_mask = Bscal >= 0.05
bulk_mask_v = Velocity_mag >= 0.05

dni[dni_mask]
Velocity_mag[bulk_mask_v]
Bscal[B_mask]
print("dni: ", dni[dni_mask])

print("Filtered Bscal :", Bscal[B_mask])

fl_B = Bscal[B_mask]
fl_dni = dni[dni_mask]
fl_Bulk_v = Velocity_mag[bulk_mask_v]



dBmag_span = aggregate_counts_by_bounds(source_bins, Bscal, span_edges)
dn_span = aggregate_counts_by_bounds(source_bins, dni, span_edges)
dVmag_span = aggregate_counts_by_bounds(source_bins, Velocity_mag, span_edges)

# span_widths = np.diff(span_edges)
# span_labels = get_span_labels(span_edges)
# bar_centers = 0.5*(span_edges[:-1] + span_edges[1:])
# #span_edges = np.array([0,200,400,600,800,1000])
# plt.bar(span_edges[:-1], dBmag_span, width=span_widths, edgecolor='black', align='edge',alpha=0.1,label='Magnetic Magnitude structures',color='blue')
# plt.bar(span_edges[:-1], dn_span, width=span_widths, edgecolor='black', align='edge',alpha=0.3,label='Density structures',color='green')
# plt.bar(span_edges[:-1], dVmag_span, width=span_widths, edgecolor='black', align='edge',alpha=0.4, label = 'Velocity Magnitude structures',color='pink')
# plt.xticks(bar_centers, span_labels, rotation=45, ha='center')
# plt.xlabel(r'Scale Size $(R_e)$ | Filter Threshold = 10%')
# plt.ylabel('Count')
# plt.yscale('log')
# plt.tight_layout()
# plt.legend()
# plt.show()

test_bins = [1,10,100,100]
plt.hist(dBmag_span, bins=test_bins, alpha=0.4, color="k", label="Mag_test")
plt.xscale("log")
plt.yscale("log")
plt.legend()
plt.show()
print(dBmag_span)
print(dn_span)
print(dVmag_span)



# %%
