#%%
import pyspedas
import numpy as np
import matplotlib.pyplot as plt
import pytplot
from pytplot import get_data, store_data,timespan,tplot_options
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


# %%

# Define the [start,stop] time for the interval you are interested in 
#trange = ['2021-08-12/00:30', '2021-08-12/02:30']

# Now use the pyspedas load routines to get the data from that interval.  In this case, we'll get data relevant to magnetic fields and ions
# Bfld_vars = pyspedas.projects.psp.fields(trange=trange, level='l2', time_clip=True)
# spi_vars = pyspedas.projects.psp.spi(trange=trange,level='l2',time_clip=True)

# # Assign strings from pyspedas to more convenient variables
# B_name = 'psp_fld_l2_mag_RTN'
# vi_name = 'psp_spi_VEL_RTN_SUN'
# ni_name = 'psp_spi_DENS'

#%%

# Tplot basics

# # Tplot is useful for quickly plotting a timeseries of the data you just retrieved, for example...
# tplot(B_name)

# # Tplot can also accept a list of pyspedas variables so long as they are bracketed, for example...
# tplot([B_name,vi_name,ni_name])

#%%
import pyspedas
from pyspedas import tplot

# highest time resolutions for swe instrument: 92s (ions) & 6-12 s (electrons) 

trange=['2019-01-09/00:00', '2019-01-20/00:00']
swe_vars = pyspedas.projects.wind.swe(trange=trange)
mfi_vars = pyspedas.projects.wind.mfi(trange=trange)

tplot(['N_elec', 'T_elec','U_eGSE','BGSE'])

B_name='BGSE'
ni_name='N_elec'
vi_name = 'U_eGSE'
pos_name = 'DIST' # Distance from earth, in units of earth radii.  1 Re ~ 4.3e-5 AU





#%%

import pyspedas
from pyspedas import tplot


# highest time resolutions for swe instrument: 64s (ions) & 128 s (electrons) 

trange = ['2013-11-1/00:00', '2013-11-6/00:00']
swe_vars = pyspedas.projects.ace.swe(trange=trange)
mfi_vars = pyspedas.projects.ace.mfi(trange=trange)

B_name = 'BGSEc'
ni_name = 'Np'
vi_name = 'Vp'
pos_name = 'SC_pos_GSE' # Disance from earth in units of km.  1 km ~ 6.7e-9 AU


tplot([B_name,ni_name,vi_name,pos_name])





#%%
# Tplot variable basics

# To do more than plot the timeseries data directly, we need to extract the data from the tplot variable format

# We can use the get_data routine from pyspedas to retrieve the actual data structure
Bdata = get_data(B_name)
ndata = get_data(ni_name)

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
interpvar_name = ni_name

# I like to define a timeax variable as well, representing the time variable
interpdata = get_data(interpvar_name)
timeax = interpdata.times


# To interpolate a tplot variable to another tplot variable, we use pyspedas's tinterpol and assign new names to the interpolated versions
tinterpol(B_name,interpvar_name,newname='B')
tinterpol(ni_name,interpvar_name,newname='ni') # (Obviously, interpolating a variable to itself isn't necessary, but I keep it here because sometimes I want to change the interpvar)

# Now they should have the same length
B_name_new = get_data('B')
ni_name_new = get_data('ni')
print('Length of New Magnetic Field Array: ' + str(len(B_name_new[1])))
print('Length of New Density Array: ' + str(len(ni_name_new[1])))


# The last thing to do before we start calculating things is reforming the tplot variables into simpler numpy arrays to work with.  This is why I made the reform() function.
# At the same time, I use factors to convert the variables from their default units to their SI units.  This just makes things easier in my experience 
Bvecs = 1e-9*reform(B_name_new) # Converting from nanotesla to Tesla
ni = 1e6*reform(ni_name_new) # Converting from cm^-3 to m^-3

##%%

# Let's start by calculating some things: the magnetic field magnitude and the alfven velocity

# define the arrays we'll overwrite
Bmag = np.zeros(len(timeax)) # this just creates the array of the correct length
va = np.zeros(len(timeax))

# overwrite arrays with built-in and pre-defined functions
for i in range(len(timeax)):
	Bmag[i] = np.linalg.norm(Bvecs[i]) # np.linalg.norm automatically calculates the magnitude of a vector
	va[i] = get_va(Bvecs[i],ni[i],mi)  # using the get_va() functionVa = |B|/sqrt(mu0*n*m)  


# Now store new variables as tplot variables and plot them.  This is typically the step where I'd convert the units if I want to express something other than SI units.
store_data('Bvecs', data = {'x':timeax,'y':1e9*Bvecs}) # with a factor to convert Tesla to nanotesla 
store_data('Bmag', data = {'x':timeax,'y':1e9*Bmag}) # with a factor to convert Tesla to nanotesla 
store_data('va', data = {'x':timeax,'y':1e-3*va}) # with a factor to convert m/s to km/s
tplot(['Bvecs','Bmag','va'])


#%%

# define minuteinterval, which represents the number of array elements equivalent to 1 minute (rounding up to the nearest integer)
minuteinterval = minute_int(timeax,trange)

# To change the averaging window, I like to just say how long (in minutes) i want the averaging window to be.  
# Try adjusting the value of n_minutes to see how the plots change (it has to be an integer). 
n_minutes_lg = 120 # Larger mean interval
n_minutes_sml = 1 #

# Now we can calculate the averaged values & store the tplot variables

# Long Average
Bvecs_mean = get_vecmean(Bvecs,n_minutes_lg*minuteinterval)
Bmag_mean = get_mean(Bmag,n_minutes_lg*minuteinterval)
n_mean = get_mean(ni,n_minutes_lg*minuteinterval)
va_mean = get_mean(va,n_minutes_lg*minuteinterval)

store_data('Bvecs_mean', data = {'x':timeax,'y':1e9*Bvecs_mean}) # converted to nT
store_data('Bmag_mean', data = {'x':timeax,'y':1e9*Bmag_mean}) # converted to nT
store_data('n_mean', data = {'x':timeax,'y':1e-6*n_mean}) # with a factor to convert m/s to km/s
store_data('va_mean', data = {'x':timeax,'y':1e-3*va_mean}) # with a factor to convert m/s to km/s

# Short Average
Bvecs = get_vecmean(Bvecs,n_minutes_sml*minuteinterval)
Bmag = get_mean(Bmag,n_minutes_sml*minuteinterval)
n = get_mean(ni,n_minutes_sml*minuteinterval)
va = get_mean(va,n_minutes_sml*minuteinterval)

store_data('Bvecs', data = {'x':timeax,'y':1e9*Bvecs}) # converted to nT
store_data('Bmag', data = {'x':timeax,'y':1e9*Bmag}) # converted to nT
store_data('n', data = {'x':timeax,'y':1e-6*n}) # with a factor to convert m/s to km/s
store_data('va', data = {'x':timeax,'y':1e-3*va}) # with a factor to convert m/s to km/s




# Here's a plot of the vector magnetic field vs its averaged version as an example.
tplot(['Bvecs','Bvecs_mean'])
tplot(['n','n_mean'])












# %%

# Now lets define some delta variables.
# in this case, we'll use the difference between the measured value and the 'averaged' one
# I defined the actual method of this in the get_delta and get_deltascalar functions, which we will now use

# Define empty arrays first
dB,dB_norm = np.zeros_like(Bvecs),np.zeros_like(Bvecs)
dBmag,dBmag_norm = np.zeros_like(Bmag),np.zeros_like(Bmag)
dn,dn_norm = np.zeros_like(ni),np.zeros_like(ni)

# Get the deltas
for i in range(len(timeax)):
	dB[i], dB_norm[i]= get_delta(Bvecs[i],Bvecs_mean[i])  # dB & dB/|B| (vectors)
	dBmag[i],dBmag_norm[i] = get_deltascalar(Bmag[i],Bmag_mean[i]) # d|B| & d|B|/|B|
	dn[i],dn_norm[i] = get_deltascalar(n[i],n_mean[i])

# Store variables
store_data('dB', data = {'x':timeax,'y':1e9*dB}) # converted to nT
store_data('dB_norm', data = {'x':timeax,'y':dB_norm}) 

store_data('dBmag', data = {'x':timeax,'y':1e9*dBmag}) # converted to nT
store_data('dBmag_norm', data = {'x':timeax,'y':dBmag_norm}) 

store_data('dn', data = {'x':timeax,'y':1e-6*dn}) # converted to cm^-3
store_data('dn_norm', data = {'x':timeax,'y':dn_norm}) 

# Example plots showing variable, its average, and the normalized delta values
tplot(['Bvecs','Bvecs_mean','dB_norm'])
tplot(['ni','n_mean','dn_norm'])

# %%
