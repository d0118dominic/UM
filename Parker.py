#%%
import pyspedas
import numpy as np
import matplotlib.pyplot as plt
import pytplot
from pytplot import get_data, store_data,timespan
from pyspedas import tplot
from pyspedas import tinterpol

me = 9.1094e-31 #kg
mi = 1837*me
mu0 = 1.2566370e-06  #;m kg / C^2
eps0 = 8.85e-12   # C^2/Nm^2
e = 1.602e-19 #C
Z = 1 # 1 for H+, 2 for He2+
gamma = 5/3
kb = 1.380649e-23




# Functions to get various useful parameters
#%%

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

def duration(trange):
	from datetime import datetime as dt
	start = dt.strptime(trange[0], '%Y-%m-%d/%H:%M')
	stop = dt.strptime(trange[1], '%Y-%m-%d/%H:%M')
	duration = stop-start
	duration_s = duration.total_seconds()
	duration_m = duration_s/60 #minutes
	duration_h = duration_m/60 #hours
	return duration_m

def mean_int(timeax,trange):
	steps = len(timeax)
	minutes = duration(trange)
	mean = np.ceil(steps/minutes)
	return int(mean)

def filter_angle(angle,floor,ceiling):
	angle_reduced = angle
	for i in range(len(timeax)): 
		if (angle[i] < floor or angle[i]>ceiling): angle_reduced[i] = np.nan
	return angle_reduced

# Velocities
def get_va(B,n,m):
	va = np.linalg.norm(B)/np.sqrt(mu0*n*m)
	return va
def get_vs(T,m):
	vs = np.sqrt(gamma*Z*kb*T/m)
	return vs
def get_vth(T,m):
	vth = np.sqrt(kb*T/m)
	return vth
def get_vpar(v,B):
	vpar = np.dot(v,B)/np.linalg.norm(B)
# Scalar Pressures
def get_pm(B):
	pm = 0.5*(mu0**-1)*np.linalg.norm(B)**2
	return pm
def get_pth(n,T):
	pth = gamma*n*T
	return pth

# Normalized Parameters
def get_brnorm(B):
	brnorm = B[0]/np.linalg.norm(B)
	return brnorm

def get_deflection(B):
	brnorm = get_brnorm(B)
	theta = np.arccos(brnorm)*180/np.pi
	return theta

def get_angle(vec,meanvec):
	term1 = np.dot(vec,meanvec)
	term2 = np.dot(np.linalg.norm(vec),np.linalg.norm(meanvec))
	term3 = term1/term2
	term4 = np.arccos(term3)*180/np.pi
	return term4


# NEED TESTING.  Having pyspedas path issues...##
def get_parperps(n,T,B):  #B ant T tensor coord systems need to match for this
	trace = T[0] + T[1] + T[2]
	term1 = (T[0]*B[0]**2 + T[1]*B[1]**2 + T[2]*B[2]**2)/(np.linalg.norm(B)**2)
	term2 = 2*(T[3]*B[0]*B[1] + T[4]*B[0]*B[2] + T[5]*B[1]*B[2])/(np.linalg.norm(B)**2)
	Tpar = term1+term2
	Tperp=0.5*(trace-Tpar)
	Ppar = n*kb*Tpar
	Pperp = n*kb*Tperp
	return Tpar,Tperp,Ppar,Pperp

def get_vxb(v,B):
	vxb = np.cross(v,B)
	return vxb

def get_K(m,n,v):
	K = 0.5*m*n*v**3
	return K

def get_H(v,P):
	H = 0.5*v*np.trace(P) + np.dot(v,P)
	return H

def get_S(E,B):
	S = np.cross(E,B)/mu0
	return S



def get_mean(var,int):   #Take a timeseries and compute the mean (basically smooth outsmall flucs over some interval)
	box = np.ones(int)/int
	smoothed_var = np.convolve(var,box,mode='same')
	return smoothed_var

def get_vecmean(vec,int):   # Vector mean
	vec1_mean = get_mean(vec[:,0],minutes*meaninterval)
	vec2_mean = get_mean(vec[:,1],minutes*meaninterval)
	vec3_mean = get_mean(vec[:,2],minutes*meaninterval)
	vec_mean = np.zeros_like(vec)
	for i in range(len(vec_mean)): vec_mean[i] = np.array([vec1_mean[i],vec2_mean[i],vec3_mean[i]])
	return vec_mean
# %%
# Get Mag Data

# trange = ['2024-10-29/00:00', '2024-10-31/00:00']
# trange = ['2018-11-1', '2018-11-10'] #Dudok de wit 2020 Full interval
#trange = ['2021-04-28', '2021-04-30'] # Encounter 8 (some sub-Alfvenic)
#trange = ['2021-08-09/12:00', '2021-08-10/00:00'] # Encounter 9 (some sub-Alfvenic)
#trange = ['2022-02-25', '2022-02-28'] #Dudok de wit 2020 Full interval

# trange = ['2018-11-05/00:00', '2018-11-05/03:00'] # Bale 2019 event (includes Sr)
trange = ['2021-08-11/09:00', '2021-08-12/09:00'] # Soni 2024 Parker interval
#trange = ['2024-09-30/00:00', '2024-09-30/23:59'] # E21
# trange = ['2024-12-24/00:00', '2024-12-25/00:00'] # E22
# trange = ['2025-03-22/00:00', '2025-03-23/00:00'] # E23
# trange = ['2025-06-19/00:00', '2025-06-20/00:00'] # E24

Bfld_vars = pyspedas.projects.psp.fields(trange=trange, level='l2', time_clip=True)
ACfld_vars = pyspedas.projects.psp.fields(trange=trange, datatype='dfb_ac_spec', level='l2',time_clip=True)
DCfld_vars = pyspedas.projects.psp.fields(trange=trange, datatype='dfb_dc_spec', level='l2',time_clip=True)

# spc_vars = pyspedas.projects.psp.spc(trange=trange, datatype='l2', level='l2')
swp_vars = pyspedas.projects.psp.spi(trange=trange,level='l3',time_clip=True)

#%%1G
# Reform all data to simple arrays and convert to SI units 
B_name = 'psp_fld_l2_mag_RTN'
vi_name = 'psp_spi_VEL_RTN_SUN'
Bxyz_name = 'psp_spi_MAGF_INST'
vxyz_name = 'psp_spi_VEL_INST'
TiTensor_name = 'psp_spi_T_TENSOR_INST'
Ti_name = 'psp_spi_TEMP'
ni_name = 'psp_spi_DENS'

interpvar_name = vi_name
timeax = pytplot.get_data(interpvar_name).times
meaninterval = mean_int(timeax,trange)

tinterpol(B_name,interpvar_name,newname='B')
tinterpol(Bxyz_name,interpvar_name,newname='Bxyz')
tinterpol(vxyz_name,interpvar_name,newname='vxyz')
tinterpol(vi_name,interpvar_name,newname='vi')
tinterpol(Ti_name,interpvar_name,newname='Ti')
tinterpol(ni_name,interpvar_name,newname='ni')
tinterpol(TiTensor_name,interpvar_name,newname='TiTensor')

Bvecs = 1e-9*reform(pytplot.get_data('B'))
Bxyz = 1e-9*reform(pytplot.get_data('Bxyz'))
vxyz = 1e3*reform(pytplot.get_data('vxyz'))
vivecs = 1e3*reform(pytplot.get_data('vi'))
ni = 1e6*reform(pytplot.get_data('ni'))
Ti = 1.602e-19*reform(pytplot.get_data('Ti'))
TiTensor = 1.602e-19*reform(pytplot.get_data('TiTensor'))

# %%
# Calculate B Magnitude & create Normalized Br/|B|
Br,Bt,Bn = np.zeros_like(Bvecs[:,0]),np.zeros_like(Bvecs[:,0]),np.zeros_like(Bvecs[:,0])
vr,vt,vn = np.zeros_like(Bvecs[:,0]),np.zeros_like(Bvecs[:,0]),np.zeros_like(Bvecs[:,0])

B_mag = np.zeros_like(Bvecs[:,0])
v_mag = np.zeros_like(Bvecs[:,0])

Br_norm = np.zeros_like(Bvecs[:,0])
Br_mean = np.zeros_like(Bvecs[:,0])
vr_norm = np.zeros_like(Bvecs[:,0])
E_conv = np.zeros_like(Bvecs)

S = np.zeros_like(Bvecs)
H = np.zeros_like(Bvecs)
K = np.zeros_like(Bvecs)
ExB = np.zeros_like(Bvecs)

P_mag = np.zeros_like(Br_norm)
P_th = np.zeros_like(ni)

va = np.zeros_like(Br_norm)
vth = np.zeros_like(Br_norm)
vs = np.zeros_like(Br_norm)
Tpar = np.zeros_like(Br_norm)
Tperp = np.zeros_like(Br_norm)
Ppar = np.zeros_like(Br_norm)
Pperp = np.zeros_like(Br_norm)

beta = np.zeros([len(timeax),2])
Br_norm = np.zeros([len(timeax),2])
viandva = np.zeros([len(timeax),2])
v_ratio = np.zeros([len(timeax),2])
pressures = np.zeros([len(timeax),3])
parperps = np.zeros([len(timeax),2])
line = np.ones_like(ni)
theta = np.zeros([len(timeax),2])
theta_v = np.zeros([len(timeax),2])

for i in range(len(Bvecs)):
	Br[i], Bt[i], Bn[i] = Bvecs[i,0], Bvecs[i,1], Bvecs[i,2]
	vr[i], vt[i], vn[i] = vivecs[i,0], vivecs[i,1], vivecs[i,2]
	B_mag[i] = np.linalg.norm(Bvecs[i,:])
	v_mag[i] = np.linalg.norm(vivecs[i,:])
	
	E_conv[i] = -get_vxb(vivecs[i],Bvecs[i])
	S[i] = get_S(E_conv[i],Bvecs[i])
	K[i] = get_K(mi,ni[i],vivecs[i])

	Tpar[i],Tperp[i],Ppar[i],Pperp[i] = get_parperps(ni[i],TiTensor[i],Bxyz[i])
	
	ExB[i] = get_vxb(E_conv[i],Bvecs[i]) # Will probably need to revisit under better assumptions
	Br_norm[i] = [get_brnorm(Bvecs[i]),0]
	vr_norm[i] = get_brnorm(vivecs[i])
	P_mag[i] = get_pm(Bvecs[i])
	P_th[i] = get_pth(ni[i],Ti[i])
	va[i] = get_va(Bvecs[i],ni[i],mi)
	vs[i] = get_vs(Ti[i],mi)
	vth[i]  = get_vth(Ti[i],mi)
	beta[i] = [P_th[i]/P_mag[i],1]
	v_ratio[i] = [np.linalg.norm(vivecs[i])/va[i],1]
	viandva[i] = [np.linalg.norm(vivecs[i]),va[i]]
	pressures[i] = [P_mag[i] + P_th[i],P_mag[i],P_th[i]]
	parperps[i] = [Ppar[i], Pperp[i]]
	theta[i] = [get_deflection(Bvecs[i]),90]
	theta_v[i] = [get_deflection(vivecs[i]),90]

#%%
# Get minute averaged quantities.  meaninterval = 1 min ##
minutes = 1

Bvecs_mean = get_vecmean(Bvecs,minutes*meaninterval)
vivecs_mean = get_vecmean(vivecs,minutes*meaninterval)
Bmag_mean = get_mean(B_mag,minutes*meaninterval)
vmag_mean = get_mean(v_mag,minutes*meaninterval)
n_mean = get_mean(ni,minutes*meaninterval)
va_mean = Bmag_mean/np.sqrt(mu0*mi*n_mean)
ma_mean = vmag_mean/va_mean
beta_mean = get_mean(beta[:,0], minutes*meaninterval)

plt.plot(ma_mean)

#%%
# Get deflection angle & related proxies
angle = np.zeros_like(timeax)
for i in range(len(timeax)): angle[i] = get_angle(Bvecs[i],Bvecs_mean[i])

angle_reduced = filter_angle(angle,10,20)
plt.plot(angle_reduced)
plt.ylim(0,180)

#%%
# Br_short = get_mean(Br,meaninterval)
# Vr_short = get_mean(Vr,meaninterval)

#deflection_param = Vr/Vr_mean 
#plt.plot((Br_short - Br_mean)/Br_mean)
#plt.plot((Vr_short - Vr_mean)/Vr_mean)
plt.xlabel('Time')
plt.ylabel('')
plt.plot(ma_mean,label='Mach Number')
plt.plot(beta_mean,label='Pth/Pm')
plt.axhline(y=1,color='k',linestyle = '--')
plt.ylim(0,4)
plt.legend()

# %%
# Make Tplot variables

#--------------------------------------------------------------------------------------
# Energy Fluxes.  This is still shaky, as it assumes E = -vxB then uses that E for ExB
# In principle, B x (-vxB) = 0 by definition, so idk why it becomes a non-zero signal
# Better versions of this should use direct Efield data or maybe use background parker spiral for the ExB
# but not for the calculation of -vxB

store_data('S', data = {'x':timeax,'y':S})
store_data('Sr', data = {'x':timeax,'y':S[:,0]})
store_data('St', data = {'x':timeax,'y':S[:,1]})
store_data('Sn', data = {'x':timeax,'y':S[:,2]})
store_data('Sr_norm', data = {'x':timeax,'y':S[:,0]/(P_mag)}) 
store_data('K', data = {'x':timeax,'y':K})
store_data('Kr', data = {'x':timeax,'y':K[:,0]})
store_data('Kt', data = {'x':timeax,'y':K[:,1]})
store_data('Kn', data = {'x':timeax,'y':K[:,2]})
store_data('Kr_norm', data = {'x':timeax,'y':K[:,0]/(P_th)})
store_data('E_conv', data = {'x':timeax,'y':E_conv})
#---------------------------------------------------------------------------------------



# ------Angles---------------------------------
#store_data('theta', data = {'x':timeax,'y':theta})
#store_data('theta_v', data = {'x':timeax,'y':theta_v})
#store_data('deflection', data = {'x':timeax,'y':deflection_param})
#----------------------------------------------



# ----Pressures&Temps-------------------------------
store_data('Pt', data = {'x':timeax,'y':P_th})
store_data('Pm', data = {'x':timeax,'y':P_mag})
store_data('Ppar', data = {'x':timeax,'y':Ppar})
store_data('Pperp', data = {'x':timeax,'y':Pperp})
store_data('pressures', data = {'x':timeax,'y':pressures})
store_data('parperps', data = {'x':timeax,'y':parperps})

#---------------------------------------------



## -----Various Dimensionless Parameters -------
# Alfvenic Mach Number (Ma)
store_data('Ma', data = {'x':timeax,'y':v_ratio})

# Plasma (proton) beta
store_data('beta', data = {'x':timeax,'y':beta})

# Br/|B|
store_data('Br_norm', data = {'x':timeax,'y':Br_norm})

# vr/|v|
store_data('vr_norm', data = {'x':timeax,'y':vr_norm})



store_data('vs_ratio', data = {'x':timeax,'y':vs/va})

#-----------------------------------------------



# ----Ion moments-----
store_data('vi', data = {'x':timeax,'y':vivecs})
store_data('n', data = {'x':timeax,'y':ni})
store_data('T', data = {'x':timeax,'y':Ti})
#----------------------


# Both Velocities on one plot
store_data('viandva', data = {'x':timeax,'y':viandva})


# %%
# Construct Timeseries plots for eventual paper
# Plot 1: 

pyspedas.options('Ma', 'ytitle', 'Vp/Va')
pyspedas.options('Ma', 'ylog', 1)
pyspedas.options('beta', 'ylog', 1)

pyspedas.ylim('Ma',0.1,10)
pyspedas.ylim('beta', 0.1,10)
pyspedas.ylim('theta',0,180)

pyspedas.tsmooth('Ma',0) # creates a 'Ma-s' variable
pyspedas.options('Ma', 'linestyle', ['-','--'])
pyspedas.options('Ma', 'color', 'k')
pyspedas.options('theta', 'linestyle', ['-','--'])
pyspedas.options('theta', 'color', 'k')
pyspedas.options('theta', 'ytitle', 'Deflection Angle (degrees)')
pyspedas.options('beta','ytitle','Plasma Beta')
pyspedas.options('beta', 'linestyle', ['-','--'])
pyspedas.options('beta', 'color', 'k')


pyspedas.options('Br_norm','color', 'k')
pyspedas.options('Br_norm','linestyle', ['-','--'])
tplot(['theta', 'beta',])

# %%

# Scatter Plot w/ colorbar
cm = plt.cm.get_cmap('seismic')
sc=plt.scatter(ma_mean,theta[:,0],c=S[:,0],s=3,cmap=cm,vmin=-0.01,vmax=0.01)
#sc=plt.scatter(ma_mean,Br_brmean,c=S[:,0],s=3,cmap=cm,vmin=-0.3,vmax=0.3)

plt.colorbar(sc,label="sr")
plt.xlim(0,2)
plt.ylim(0,180)


plt.axhline(y=90,c='k')
plt.axvline(x=1,c='k')
plt.xlabel('Alfven Mach Number Ma')
plt.ylabel('Deflection Angle')
plt.show()



# %%


# %%
# %%
