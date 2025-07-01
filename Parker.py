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
k = 1.380649e-23

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



# Functions to get various useful parameters
#%%
# Velocities
def get_va(B,n,m):
	va = np.linalg.norm(B)/np.sqrt(mu0*n*m)
	return va
def get_vs(T,m):
	vs = np.sqrt(gamma*Z*k*T/m)
	return vs
def get_vth(T,m):
	vth = np.sqrt(k*T/m)
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

# NEED TESTING ##
def get_ppar(n,T,B):  #B ant T tensor coord systems need to match for this
	term1 = (T[0]*B[0]**2 + T[1]*B[1]**2 + T[2]*B[2]**2)/(np.linalg.norm(B)**2)
	term2 = 2*(T[3]*B[0]*B[1] + T[4]*B[0]*B[2] + T[5]*B[1]*B[2])/(np.linalg.norm(B)**2)
	T_par = term1+term2
	P_par = n*T_par
	return P_par

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

def get_Tcomps(Ttensor,B):
	pass
	return

def get_mean(var,int):   #Take a timeseries and compute the mean (basically smooth outsmall flucs over some interval)
	box = np.ones(int)/int
	smoothed_var = np.convolve(var,box,mode='same')
	return smoothed_var
# %%
# Get Mag Data

# trange = ['2024-10-29/00:00', '2024-10-31/00:00']
# trange = ['2018-11-1', '2018-11-10'] #Dudok de wit 2020 Full interval
#trange = ['2021-04-28', '2021-04-30'] # Encounter 8 (some sub-Alfvenic)
#trange = ['2021-08-09/12:00', '2021-08-10/00:00'] # Encounter 9 (some sub-Alfvenic)
#trange = ['2022-02-25', '2022-02-28'] #Dudok de wit 2020 Full interval
trange = ['2018-11-05/00:00', '2018-11-05/03:00'] # Bale 2019 event (includes Sr)
#trange = ['2021-08-11/09:00', '2021-08-12/09:00'] # Soni 2024 Parker interval

Bfld_vars = pyspedas.projects.psp.fields(trange=trange, level='l2', time_clip=True)
ACfld_vars = pyspedas.projects.psp.fields(trange=trange, datatype='dfb_ac_spec', level='l2',time_clip=True)
DCfld_vars = pyspedas.projects.psp.fields(trange=trange, datatype='dfb_dc_spec', level='l2',time_clip=True)

# spc_vars = pyspedas.projects.psp.spc(trange=trange, datatype='l2', level='l2')
swp_vars = pyspedas.projects.psp.spi(trange=trange,level='l3',time_clip=True)

#%%1G
# Reform all data to simple arrays and convert to SI units 
B_name = 'psp_fld_l2_mag_RTN'
vi_name = 'psp_spi_VEL_RTN_SUN'
vsc_name = 'psp_spi_VEL_SC'
TiTensor_name = 'psp_spi_T_TENSOR_INST'
Ti_name = 'psp_spi_TEMP'
ni_name = 'psp_spi_DENS'

interpvar_name = vi_name
timeax = pytplot.get_data(interpvar_name).times

tinterpol(B_name,interpvar_name,newname='B')
tinterpol(vi_name,interpvar_name,newname='vi')
tinterpol(vsc_name,interpvar_name,newname='vsc')
tinterpol(Ti_name,interpvar_name,newname='Ti')
tinterpol(ni_name,interpvar_name,newname='ni')
tinterpol(TiTensor_name,interpvar_name,newname='TiTensor')

Bvecs = 1e-9*reform(pytplot.get_data('B'))
vivecs = 1e3*reform(pytplot.get_data('vi'))
vscvecs = 1e3*reform(pytplot.get_data('vsc'))
ni = 1e6*reform(pytplot.get_data('ni'))
Ti = 1.602e-19*reform(pytplot.get_data('Ti'))
TiTensor = 1.602e-19*reform(pytplot.get_data('TiTensor'))

# %%
# Calculate B Magnitude & create Normalized Br/|B|

Br_norm = np.zeros_like(Bvecs[:,0])
Br = np.zeros_like(Bvecs[:,0])
Vr = np.zeros_like(Bvecs[:,0])
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

PiTensor = np.zeros_like(TiTensor)


beta = np.zeros([len(timeax),2])
Br_norm = np.zeros([len(timeax),2])
viandva = np.zeros([len(timeax),2])
v_ratio = np.zeros([len(timeax),2])
pressures = np.zeros([len(timeax),3])
line = np.ones_like(ni)
theta = np.zeros([len(timeax),2])
theta_v = np.zeros([len(timeax),2])

for i in range(len(Bvecs)):
	Br[i] = Bvecs[i,0]
	Vr[i] = vivecs[i,0]
	E_conv[i] = -get_vxb(vivecs[i],Bvecs[i])
	S[i] = get_S(E_conv[i],Bvecs[i])
	K[i] = get_K(mi,ni[i],vivecs[i])
	
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
	PiTensor[i] = ni[i]*TiTensor[i]
	viandva[i] = [np.linalg.norm(vivecs[i]),va[i]]
	pressures[i] = [P_mag[i] + P_th[i],P_mag[i],P_th[i]]
	theta[i] = [get_deflection(Bvecs[i]),90]
	theta_v[i] = [get_deflection(vivecs[i]),90]

Br_mean = get_mean(Br,5)
Vr_mean = get_mean(Vr,5)
plt.plot((Br - Br_mean)/Br_mean)
plt.plot((Vr - Vr_mean)/Vr_mean)
plt.xlabel('Time')
plt.ylabel('dBr/Br , dVr/Vr')
plt.ylim(-10,10)

# %%
# Make Tplot variables

#--------------------------------------------------------------------------------------
# Energy Fluxes.  This is still shaky, as it assumes E = -vxB then uses that E for ExB
# In principle, E x (-vxB) = 0 by definition, so idk why it becomes a non-zero signal
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
store_data('theta', data = {'x':timeax,'y':theta})
store_data('theta_v', data = {'x':timeax,'y':theta_v})
#----------------------------------------------



# ----Pressures-------------------------------
store_data('Pt', data = {'x':timeax,'y':P_th})
store_data('Pm', data = {'x':timeax,'y':P_mag})
store_data('pressures', data = {'x':timeax,'y':pressures})
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
tplot(['Ma', 'beta','theta','Br_norm'])

# %%

# Scatter Plot w/ colorbar
cm = plt.cm.get_cmap('seismic')
sc=plt.scatter(v_ratio[:,0],theta[:,0],c=beta[:,0],s=1,vmin=0,vmax=2,cmap=cm)

plt.colorbar(sc,label="Radial Poynting Flux Sr")
plt.xlim(0,2)
plt.ylim(0,180)

plt.axhline(y=90,c='k')
plt.axvline(x=1,c='k')
plt.xlabel('Alfven Mach Number Ma')
plt.ylabel('Deflection Angle (degrees)')
plt.show()



# %%


# %%


# %%
