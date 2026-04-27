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


def filter_deflections(var,threshold):
	# newvar = np.array([])
	newvar = []
	for i in range(len(var)):
		if (abs(var[i])>=threshold):
			newvar.append(var[i]) 
		else: 
			pass
	return np.array(newvar)


#%%

# WIND data (if you want to use WIND data, run this cell and skip the next cell)

# highest time resolutions for swe instrument: 92s (ions) & 6-12 s (electrons) 
# Assuming solar wind velocity of 400 km/s, this corresponds to structure sizes of 6 & 0.4-0.8 Re

trange=['2019-01-01/00:00', '2019-01-10/00:00']
# trange=['2019-01-20/00:00', '2019-01-30/00:00']
swe_vars = pyspedas.projects.wind.swe(trange=trange)
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
interpvar_name = ni_name

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
n_minutes_lg = 54
 # Larger mean interval
n_minutes_sml = 27 #

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
# store_data('n_inst', data = {'x':timeax,'y':1e-6*n}) # with a factor to convert m/s to km/s
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
tplot(['Bmag','Bmag_mean','dBmag_norm'])
# tplot(['Bvecs_inst','Bvecs','Bvecs_mean','dB_norm'])
# tplot(['ni','n_mean','dn_norm'])


# Here I am defining a threshold for the normalized delta variables
filter_threshold = 0.05

# Any data below the threshold is erased, leaving behind a new array with only the large spikes
# The reason for this is we want to 'count' the data that exceed the threshold
dBmag_norm_filtered = filter_deflections(dBmag_norm,filter_threshold)
dBx_norm_filtered = filter_deflections(dB_norm[:,0],filter_threshold)
dBy_norm_filtered = filter_deflections(dB_norm[:,1],filter_threshold)
dBz_norm_filtered = filter_deflections(dB_norm[:,2],filter_threshold)


dvmag_norm_filtered = filter_deflections(dvmag_norm,filter_threshold)
dvx_norm_filtered = filter_deflections(dv_norm[:,0],filter_threshold)
dvy_norm_filtered = filter_deflections(dv_norm[:,1],filter_threshold)
dvz_norm_filtered = filter_deflections(dv_norm[:,2],filter_threshold)

dn_norm_filtered = filter_deflections(dn_norm,filter_threshold)

# We are left with a shorter array that contains only the normalized spikes above threshold
print('Length of Full Array: ' + str(len(dBmag_norm)))
print('Length of Filtered Array: ' + str(len(dBmag_norm_filtered)))



#%%

# Here I am setting up a loop to conduct the previous cell's analysis for a series of averaging windows


# I make a t_windows array, which lists 11 different averaging durations (in minutes)
# 11 averaging durations means there are 10 gaps between them
t_windows=np.array([1,27,54,81,108,135,162,189,216,243,270]) # minutes

# d_windows array isn't really meant to be used right now. Just for reference 
# It corresponds to the approx structure size (in Re) corresponding to the same element in the t_windows array
# example: if our large and small windows are 81 and 54 minutes, then we are (roughly) selecting for structures between 200 and 300 Re
d_windows=np.array([4,100,200,300,400,500,600,700,800,900,1000]) # Re

# I'm sure there's a better way to do this, but for now I made a bins array that makes the plotting of counts and size more convenient
bins = [50,150,250,350,450,550,650,750,850,950]


# First define an empty list which will eventually have 10 elements representing the count between each pair of averaging windows
dBmag_count =[]
dvmag_count =[]
dvx_count =[]
dvy_count =[]
dvz_count =[]
dBx_count =[]
dBy_count =[]
dBz_count =[]
dn_count =[]

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


	filter_threshold = 0.01

	dBmag_norm_filtered = filter_deflections(dBmag_norm,filter_threshold)
	dBx_norm_filtered = filter_deflections(dB_norm[:,0],filter_threshold)
	dBy_norm_filtered = filter_deflections(dB_norm[:,1],filter_threshold)
	dBz_norm_filtered = filter_deflections(dB_norm[:,2],filter_threshold)


	dvmag_norm_filtered = filter_deflections(dvmag_norm,filter_threshold)
	dvx_norm_filtered = filter_deflections(dv_norm[:,0],filter_threshold)
	dvy_norm_filtered = filter_deflections(dv_norm[:,1],filter_threshold)
	dvz_norm_filtered = filter_deflections(dv_norm[:,2],filter_threshold)
	dn_norm_filtered = filter_deflections(dn_norm,filter_threshold)

	dBmag_count.append(len(dBmag_norm_filtered))
	dBx_count.append(len(dvx_norm_filtered))
	dBy_count.append(len(dvy_norm_filtered))
	dBz_count.append(len(dvz_norm_filtered))
	dvmag_count.append(len(dvmag_norm_filtered))
	dn_count.append(len(dn_norm_filtered))

    # dvy_count.append(len(dvx_norm_filtered)
    # dvz_count.append(len(dvz_norm_filtered))
    # dvx_count.append(len(dvx_norm_filtered))

bin_width = bins[1] - bins[0] if len(bins) > 1 else 25
plt.bar(bins, dBx_count, width=bin_width, edgecolor='blue', align='edge',alpha=0.2,label='Bx structures',color='blue')
plt.bar(bins, dBy_count, width=bin_width, edgecolor='green', align='edge',alpha=0.2,label='By',color='green')
plt.bar(bins, dBz_count, width=bin_width, edgecolor='red', align='edge',alpha=0.2, label = 'Bz structures',color='red')
plt.xlabel(r'Scale Size $(R_e)$')
plt.ylabel('Count')
plt.yscale('log')
plt.legend()


plt.show()



# %%

# UNDER CONSTRUCTION

#Regular Line plot version of the above histograms
# Here I am setting up a loop to conduct the previous cell's analysis for a series of averaging windows


Re = 6378.14 #km
v_sw = 400 #km/s
v = v_sw*(60/Re)

d_windows=np.array([4,100,200,300,400,500,600,700,800,900,1000]) # Re
d_windows=np.array([4,50,100,150,200,250,300,350,400,450,500]) # Re
t_windows = []
for i in d_windows: t_windows.append(int(np.round(i/v)))
t_windows = np.array(t_windows)

bins = np.zeros(len(d_windows)-1)

for i in range(len(d_windows)-1): bins[i] = d_windows[i] + (d_windows[i+1]-d_windows[i])/2

#%%


# First define an empty list which will eventually have 10 elements representing the count between each pair of averaging windows
dBmag_count =[]
dvmag_count =[]
dn_count =[]

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

	dBmag_norm_filtered = filter_deflections(dBmag_norm,filter_threshold)
	dvmag_norm_filtered = filter_deflections(dvmag_norm,filter_threshold)
	dn_norm_filtered = filter_deflections(dn_norm,filter_threshold)

	dBmag_count.append(len(dBmag_norm_filtered))
	dvmag_count.append(len(dvmag_norm_filtered))
	dn_count.append(len(dn_norm_filtered))


bin_width = bins[1] - bins[0] if len(bins) > 1 else 25
plt.bar(bins, dBmag_count, width=bin_width, edgecolor='black', align='edge',alpha=0.1,label='Magnetic Magnitude structures',color='b')
plt.bar(bins, dn_count, width=bin_width, edgecolor='black', align='edge',alpha=0.3,label='Density structures',color='k')
plt.bar(bins, dvmag_count, width=bin_width, edgecolor='black', align='edge',alpha=0.3, label = 'Velocity Magnitude structures',color='c')
plt.xlabel(r'Scale Size $(R_e)$')
plt.ylabel('Count')
plt.yscale('log')
plt.legend()


plt.show()

# %%
