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

# NEED TESTING ##
def get_ppar(n,T,B):  #B ant T tensor coord systems need to match for this
	term1 = (T[0]*B[0]**2 + T[1]*B[1]**2 + T[2]*B[2]**2)/(np.linalg.norm(B)**2)
	term2 = 2*(T[3]*B[0]*B[1] + T[4]*B[0]*B[2] + T[5]*B[1]*B[2])/(np.linalg.norm(B)**2)
	T_par = term1+term2
	P_par = n*T_par
	return P_par



def get_Tcomps(Ttensor,B):
	pass
	return

# %%
# Get Mag Data

trange = ['2024-10-29/00:00', '2024-10-31/00:00']
# trange = ['2018-11-1', '2018-11-10'] #Dudok de wit 2020 Full interval
trange = ['2022-02-25', '2022-02-28'] #Dudok de wit 2020 Full interval
trange = ['2021-08-11/09:00', '2021-08-12/09:00'] # Soni 2024 Parker interval

Bfld_vars = pyspedas.projects.psp.fields(trange=trange, level='l2', time_clip=True)
ACfld_vars = pyspedas.projects.psp.fields(trange=trange, datatype='dfb_ac_spec', level='l2',time_clip=True)
DCfld_vars = pyspedas.projects.psp.fields(trange=trange, datatype='dfb_dc_spec', level='l2',time_clip=True)

spc_vars = pyspedas.projects.psp.spc(trange=trange, datatype='l2', level='l2')
swp_vars = pyspedas.projects.psp.spi(trange=trange,level='l3',time_clip=True)

#%%1G
# Reform all data to simple arrays and convert to SI units 
B_name = 'psp_fld_l2_mag_RTN'
vi_name = 'psp_spi_VEL_RTN_SUN'
vsc_name = 'psp_spi_VEL_SC'
TiTensor_name = 'psp_spi_T_TENSOR_INST'
Ti_name = 'psp_spi_TEMP'
ni_name = 'psp_spi_DENS'

interpvar_name = B_name
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
vr_norm = np.zeros_like(Bvecs[:,0])
P_mag = np.zeros_like(Br_norm)
P_th = np.zeros_like(ni)
va = np.zeros_like(Br_norm)
vth = np.zeros_like(Br_norm)
vs = np.zeros_like(Br_norm)
beta = np.zeros_like(ni)
v_ratio = np.zeros_like(ni)
PiTensor = np.zeros_like(TiTensor)
viandva = np.zeros([len(timeax),2])

for i in range(len(Bvecs)):
	Br_norm[i] = np.abs(get_brnorm(Bvecs[i]))
	vr_norm[i] = np.abs(get_brnorm(vivecs[i]))
	P_mag[i] = get_pm(Bvecs[i])
	P_th[i] = get_pth(ni[i],Ti[i])
	va[i] = get_va(Bvecs[i],ni[i],mi)
	vs[i] = get_vs(Ti[i],mi)
	vth[i]  = get_vth(Ti[i],mi)
	beta[i] = P_th[i]/P_mag[i]
	v_ratio[i] = np.linalg.norm(vivecs[i])/va[i]
	PiTensor[i] = ni[i]*TiTensor[i]
	viandva[i] = [np.linalg.norm(vivecs[i]),va[i]]
#plt.plot(beta,color='b')
# plt.plot(v_ratio,color='r',label='Va/Vi')
# plt.plot(beta,color='b',label='Pth/Pm')
plt.plot(v_ratio)
# plt.ylim(0,5)
plt.legend()
# plt.axhline(y=1,linestyle='--',color='k')
# %%
# Make it a Tplot variable & plot it
store_data('viandva', data = {'x':timeax,'y':viandva})
store_data('Br_norm', data = {'x':timeax,'y':Br_norm})
store_data('vr_norm', data = {'x':timeax,'y':vr_norm})
tplot(['Br_norm','vr_norm','viandva'])

# %%
