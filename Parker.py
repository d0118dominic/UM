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

# Convenient Plot Functions
def quickplot():
	pyspedas.options(vi_name,'legend_names',['vR','vT','vN'])
	pyspedas.tplot([B_name,vi_name,Ti_name,ni_name])
	return
def machplot():
	mach = np.zeros([len(timeax),2])
	for i in range(len(timeax)):
		mach[i] = [ma_mean[i],1]
	store_data('Mach', data = {'x':timeax,'y':mach})
	pyspedas.options('Mach', 'ytitle', 'Alfven Mach Number')
	pyspedas.options('Mach', 'ylog', 1)
	pyspedas.options('Mach', 'linestyle', ['-','--'])
	pyspedas.options('Mach', 'color', ['k','r'])
	pyspedas.ylim('Mach',0.1,10)
	pyspedas.options('velocities','ytitle','Ion & Alfven Velocities')
	pyspedas.options('velocities','color',['r','g'])
	pyspedas.options('velocities','legend_names',['|vi|','|va|'])
	pyspedas.tplot(['Mach','velocities'])
	return
def fluxplot():
	pyspedas.options(['S','K'], 'legend_names',['R','T','N'])
	pyspedas.tplot(['S','K'])
	return
def presplot():
	pyspedas.options('pressures', 'legend_names',['Total','Pm','Pth'])
	pyspedas.options('pressures','color',['k','g','r'])
	pyspedas.options('beta', 'color',['k','k'])
	pyspedas.options('beta', 'linestyle',['-','--'])
	pyspedas.tplot(['pressures','beta'])
	return
def angleplot():
	mach = np.zeros([len(timeax),2])
	for i in range(len(timeax)):
		mach[i] = [ma[i],1]
	store_data('Mach', data = {'x':timeax,'y':mach})
	pyspedas.options('Mach', 'ytitle', 'Alfven Mach Number')
	pyspedas.options('Mach', 'ylog', 1)
	pyspedas.options('Mach', 'linestyle', ['-','--'])
	pyspedas.options('Mach', 'color', ['k','r'])
	pyspedas.ylim('Mach',0.1,10)
	pyspedas.options('angle','ytitle','Deflection Angle (degrees)')
	pyspedas.tplot(['Mach','angle','dB_norm'])
def betaplot():
	pyspedas.options('beta', 'color',['k','k'])
	pyspedas.options('beta', 'linestyle',['-','--'])
	pyspedas.tplot('beta')
def duration(trange):
	from datetime import datetime as dt
	start = dt.strptime(trange[0], '%Y-%m-%d/%H:%M')
	stop = dt.strptime(trange[1], '%Y-%m-%d/%H:%M')
	duration = stop-start
	duration_s = duration.total_seconds()
	duration_m = duration_s/60 #minutes
	duration_h = duration_m/60 #hours
	return duration_m


# Various mean,filter,binning functions
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

def bin_swvel(vel):
	vel_binned = vel
	for i in range(len(timeax)):
		if vel[i,0] <= 450e3: vel[i,0]=1
		elif vel[i,0] > 450e3: vel[i,0]=2
	return vel_binned

def bin_angle(angle):
	angle_binned = angle
	for i in range(len(angle)):
		if (angle[i]>=0 and angle[i]<=15): angle_binned[i] = 8
		elif (angle[i]>=15 and angle[i]<30): angle_binned[i] = 23 
		elif (angle[i]>=30 and angle[i]<45): angle_binned[i] = 38 
		elif (angle[i]>=45 and angle[i]<60): angle_binned[i] = 53 
		elif (angle[i]>=60 and angle[i]<75): angle_binned[i] = 68 
		elif (angle[i]>=75 and angle[i]<90): angle_binned[i] = 83 
		elif (angle[i]>=90 and angle[i]<105): angle_binned[i] = 98 
		elif (angle[i]>=105 and angle[i]<120): angle_binned[i] = 113 
		elif (angle[i]>=120 and angle[i]<135): angle_binned[i] = 128
		elif (angle[i]>=135 and angle[i]<150): angle_binned[i] = 143 
		elif (angle[i]>=150 and angle[i]<165): angle_binned[i] = 158 
		elif (angle[i]>=165 and angle[i]<180): angle_binned[i] = 173
		else: pass
	return angle_binned


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
	pth = kb*n*T
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

def get_parperps(n,T,B):  #B ant T tensor coord systems need to match for this
	trace = T[0] + T[1] + T[2]
	term1 = (T[0]*B[0]**2 + T[1]*B[1]**2 + T[2]*B[2]**2)/(np.linalg.norm(B)**2)
	term2 = 2*(T[3]*B[0]*B[1] + T[4]*B[0]*B[2] + T[5]*B[1]*B[2])/(np.linalg.norm(B)**2)
	Tpar = term1+term2
	Tperp=0.5*(trace-Tpar)
	Ppar = n*kb*Tpar
	Pperp = n*kb*Tperp
	return Tpar,Tperp,Ppar,Pperp


def get_par(v1,v2):
	par = np.dot(v1,v2)/np.linalg.norm(v2)
	return par
def get_voltagepairs():
	l_eff = 3.5
	pairxy = reform(get_data('psp_fld_l2_dfb_wf_dVdc_sc'))
	return pairxy/l_eff
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

# mean, vector mean, and delta vector (a-a_mean)
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

def get_delta(vec,vec_mean):
	dvec = vec - vec_mean
	dvec_norm = dvec/np.linalg.norm(vec)
	return dvec,dvec_norm

def store_alldata():
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
	store_data('angle', data = {'x':timeax,'y':angle})
	#----------------------------------------------

	# ----Pressures&Temps--(something is wrong)-----------------------------
	store_data('Pt', data = {'x':timeax,'y':P_th})
	store_data('Pm', data = {'x':timeax,'y':P_mag})
	store_data('Ppar', data = {'x':timeax,'y':Ppar})
	store_data('Pperp', data = {'x':timeax,'y':Pperp})
	store_data('pressures', data = {'x':timeax,'y':pressures})
	store_data('parperps', data = {'x':timeax,'y':parperps})
	#---------------------------------------------

	## -----Various Dimensionless Parameters -------
	# Alfvenic Mach Number (Ma)
	store_data('Ma', data = {'x':timeax,'y':ma})

	# Plasma (proton) beta
	store_data('beta', data = {'x':timeax,'y':beta})

	# vB alignment (alfvenicity?)
	store_data('vB_alignment', data = {'x':timeax,'y':vB_alignment})

	# dBr/|B| & dvr/|v| 
	store_data('dBr_norm', data = {'x':timeax,'y':dB_norm[:,0]})
	store_data('dvr_norm', data = {'x':timeax,'y':dv_norm[:,0]})

	# dB/|B| & dv/|v| 
	store_data('dB_norm', data = {'x':timeax,'y':dB_norm})
	store_data('dv_norm', data = {'x':timeax,'y':dv_norm})
	#-----------------------------------------------

	# ----Ion moments-----
	store_data('vi', data = {'x':timeax,'y':vivecs})
	store_data('n', data = {'x':timeax,'y':ni})
	store_data('T', data = {'x':timeax,'y':Ti})
	#----------------------

	store_data('velocities', data = {'x':timeax,'y':viandva})
	store_data('|vi|', data = {'x':timeax,'y':va})
	store_data('|va|', data = {'x':timeax,'y':v_mag})
	return

# %%
# Get Mag Data
# Encounters below 0.05 AU
# Encounters below 0.08 AU 
# Encounters below 0.05 AU
# Encounters below 0.06 AU
# Encounters below 0.07 AU
#trange = ['2024-10-29/00:00', '2024-10-31/00:00']
#trange = ['2018-11-3/00:00', '2018-11-6/00:00'] #Dudok de wit 2020 Full interval
#trange = ['2021-04-28/00:00', '2021-04-30/00:00'] # Encounter 8 (some sub-Alfvenic)
#trange = ['2021-08-09/12:00', '2021-08-10/00:00'] # Encounter 9 (some sub-Alfvenic)
trange = ['2021-08-09/20:00', '2021-08-10/05:00'] # Encounter 9 (some sub-Alfvenic)
#trange = ['2022-02-25', '2022-02-28'] #Dudok de wit 2020 Full interval

#

sub_alfs = 	[['2021-08-09/21:30','2021-08-10/00:00'], # 2.5 hrs Encounter 9
			 ['2021-11-21/21:00','2021-11-22/01:00'], # 4 hrs Encounter 10
			 ['2021-11-22/03:30','2021-11-22/10:00'], # 6.5 hrs Encounter 10
			 ['2022-02-25/20:00','2022-02-25/23:30'], # 3.5 hrs Encounter 11
			 ['2022-06-01/18:00','2022-06-02/08:00'], # 14 hrs Encounter 12 
			 ['2022-09-06/06:00','2022-09-06/16:00'], # 10 hrs Encounter 13
			 ['2022-09-06/18:00','2022-09-07/12:00'], # 18 hrs Encounter 13
			 ['2022-12-11/00:00','2022-12-11/18:00'], # 18 hrs Encounter 14
			 ['2023-03-16/12:00','2023-03-17/06:00'], # 18 hrs Encounter 15
			 ]


sup_alfs = 	[['2021-08-10/00:30','2021-08-10/06:00'], # 5.5 hrs Encounter 9
			 ['2021-11-21/16:00','2021-11-21/21:00'], # 5 hrs Encounter 10
			 ['2021-11-22/01:00','2021-11-22/02:30'], # 1.5 hrs Encounter 10
			 ['2022-02-26/07:30','2022-02-26/08:30'], # 1 hr Encounter 11
			 ['2022-06-02/13:00','2022-06-02/20:00'], # 7 hrs Encounter 12
			 ['2022-09-07/18:00','2022-09-08/18:00'], # 24 hrs Encounter 13
			 ['2022-12-11/01:00','2022-12-10/09:00'], # 8 hrs Encounter 14
			 ]

#sub_alfs =  ['2024-09-30/03:00','2024-09-30/11:00']
trange=sup_alfs[6]


# Alfven crossings (<2 hr)
#trange = ['2021-11-21/21:00','2021-11-21/22:00']
#trange = ['2021-08-10/00:15', '2021-08-10/00:45'] # Encounter 9 (some sub-Alfvenic)
# trange=['2023-12-29/01:30','2023-12-29/03:00'] # E18
#trange=['2024-03-29/22:00','2024-03-29/23:30']  # E19 
# trange = ['2024-06-29/11:00', '2024-06-29/13:00'] # E20
# trange = ['2024-09-28/00:00', '2024-10-05/12:00'] # E21 


#trange = ['2021-08-11/09:00', '2021-08-12/09:00'] # Soni 2024 Parker interval
#trange = ['2024-09-28/00:00', '2024-10-05/12:00'] # E21 
# trange = ['2024-06-29/00:00', '2024-07-01/12:00'] # E20
# trange=['2024-03-27/00:00','2024-04-02/00:00']  # E19 
# trange=['2023-12-28/00:00','2023-12-30/00:00'] # E18


#trange = ['2018-11-05/00:00', '2018-11-05/03:00'] # Bale 2019 event (includes Sr)
##trange = ['2024-12-24/00:00', '2024-12-25/00:00'] # E22 
# trange = ['2025-03-22/00:00', '2025-03-23/00:00'] # E23
# trange = ['2025-06-19/00:00', '2025-06-20/00:00'] # E24
# 


#Need to make a list of chosen Alfvenic & sub-alfvenic intervals
#Most people choose a handful of intervals by eye

Bfld_vars = pyspedas.projects.psp.fields(trange=trange, level='l2', time_clip=True)
spi_vars = pyspedas.projects.psp.spi(trange=trange,level='l3',time_clip=True)
# spe_vars = pyspedas.projects.psp.spi(trange=trange,level='l3',time_clip=True)
# spe_vars = pyspedas.projects.psp.spe(trange=trange,level='l2',time_clip=True)

#voltages_vars = pyspedas.projects.psp.fields(trange=trange, datatype='dfb_wf_dvdc', level='l2',time_clip=True)
#On DC datatype: 'sqn_rfs_V1V2 has some kind of electron density & core temp, but looks weird
# spc_vars = pyspedas.projects.psp.spc(trange=trange, datatype='l2', level='l2')
#%%1G
# Reform all data to simple arrays and convert to SI units 
B_name = 'psp_fld_l2_mag_RTN'
vi_name = 'psp_spi_VEL_RTN_SUN'
Bxyz_name = 'psp_spi_MAGF_INST'
vxyz_name = 'psp_spi_VEL_INST'
TiTensor_name = 'psp_spi_T_TENSOR_INST'
Ti_name = 'psp_spi_TEMP'
ni_name = 'psp_spi_DENS'
# voltages_name = 'psp_fld_l2_dfb_wf_dVdc_sc'
position_name = 'psp_spi_SUN_DIST'

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
# tinterpol(voltages_name,interpvar_name,newname='voltages')
tinterpol(position_name,interpvar_name,newname='position')
#%%
Bvecs = 1e-9*reform(pytplot.get_data('B'))
Bxyz = 1e-9*reform(pytplot.get_data('Bxyz'))
vxyz = 1e3*reform(pytplot.get_data('vxyz'))
vivecs = 1e3*reform(pytplot.get_data('vi'))
ni = 1e6*reform(pytplot.get_data('ni'))
Ti = 1.602e-19*reform(pytplot.get_data('Ti'))
TiTensor = 1.602e-19*reform(pytplot.get_data('TiTensor'))
# voltages = reform(get_data('voltages'))
position = reform(get_data('position'))/695700 #Solar radii

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
ma = np.zeros(len(timeax))
pressures = np.zeros([len(timeax),3])
parperps = np.zeros([len(timeax),2])
line = np.ones_like(ni)
theta = np.zeros([len(timeax),2])
theta_v = np.zeros([len(timeax),2])

# E_antennas = np.zeros_like(voltages)

for i in range(len(Bvecs)):
	Br[i], Bt[i], Bn[i] = Bvecs[i,0], Bvecs[i,1], Bvecs[i,2]
	vr[i], vt[i], vn[i] = vivecs[i,0], vivecs[i,1], vivecs[i,2]
	B_mag[i] = np.linalg.norm(Bvecs[i,:])
	v_mag[i] = np.linalg.norm(vivecs[i,:])
	
	E_conv[i] = -get_vxb(vivecs[i],Bvecs[i])
	# E_antennas[i] = voltages[i]/3.5
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
	ma[i] = v_mag[i]/va[i]
	beta[i] = [P_th[i]/P_mag[i],1]
	viandva[i] = [np.linalg.norm(vivecs[i]),va[i]]
	pressures[i] = [P_mag[i] + P_th[i],P_mag[i],P_th[i]]
	parperps[i] = [Ppar[i], Pperp[i]]

#%%
# Get averaged quantities.  meaninterval = 1 min ##
minutes = 10

Bvecs_mean = get_vecmean(Bvecs,minutes*meaninterval)
vivecs_mean = get_vecmean(vivecs,minutes*meaninterval)
Bmag_mean = get_mean(B_mag,minutes*meaninterval)
vmag_mean = get_mean(v_mag,minutes*meaninterval)
n_mean = get_mean(ni,minutes*meaninterval)
va_mean = Bmag_mean/np.sqrt(mu0*mi*n_mean)
ma_mean = vmag_mean/va_mean
beta_mean = get_mean(beta[:,0], minutes*meaninterval)

# Get dB,dv + components, etc. and vB alignment variable (proxy for alfvenicity)
dB,dB_norm = np.zeros_like(Bvecs),np.zeros_like(Bvecs)
dv,dv_norm = np.zeros_like(vivecs),np.zeros_like(vivecs)
dB_par,dB_perp = np.zeros_like(ni),np.zeros_like(ni)
dv_par,dv_perp = np.zeros_like(ni),np.zeros_like(ni)
B_par = np.zeros_like(ni)
v_par = np.zeros_like(ni)
dB_par_norm = np.zeros_like(ni)
dB_perp_norm = np.zeros_like(ni)
dv_par = np.zeros_like(ni)
dv_par_norm = np.zeros_like(ni)
dv_perp_norm = np.zeros_like(ni)
dB_norm_mag,dv_norm_mag = np.zeros_like(ni),np.zeros_like(ni)
vB_alignment = np.zeros_like(ni)
dvdB_alignment = np.zeros_like(ni)
angle = np.zeros_like(ni)

for i in range(len(timeax)):
	dvdB_alignment[i] = np.abs(np.dot(dB[i],dv[i]))/(np.linalg.norm(dB[i])*np.linalg.norm(dv[i]))
	vB_alignment[i] = np.abs(np.dot(Bvecs[i],vivecs[i]))/(B_mag[i]*v_mag[i]) #this one seems most useful

	# Magnetic deflection angle
	angle[i] = get_angle(Bvecs[i],Bvecs_mean[i])

	# B,v components	
	B_par[i] = get_par(Bvecs[i],Bvecs_mean[i]) # All par to mean B 
	v_par[i] = get_par(vivecs[i],Bvecs_mean[i])
	
	# dB,dv & components
	dB[i], dB_norm[i]= get_delta(Bvecs[i],Bvecs_mean[i])  # dB & dB/|B| (vectors)
	dv[i], dv_norm[i]= get_delta(vivecs[i],vivecs_mean[i]) # dv & dv/|v| (vectors)
	dB_norm_mag[i] = np.linalg.norm(dB_norm[i]) # |dB|/|B| (scalar)
	dv_norm_mag[i] = np.linalg.norm(dv_norm[i]) # |dv|/|v| (scalar)
	dB_par[i] = get_par(dB[i],Bvecs_mean[i])
	dv_par[i] = get_par(dv[i],Bvecs_mean[i])
	dB_perp[i] = np.sqrt(np.linalg.norm(dB[i])**2 - dB_par[i]**2)
	dv_perp[i] = np.sqrt(np.linalg.norm(dv[i])**2 - dv_par[i]**2)
	dB_par_norm[i] = dB_par[i]/Bmag_mean[i]
	dv_par_norm[i] = dv_par[i]/vmag_mean[i]
	dB_perp_norm[i] = dB_perp[i]/Bmag_mean[i]
	dv_perp_norm[i] = dv_perp[i]/vmag_mean[i]


low = 30
high = 180
angle_reduced = np.zeros_like(angle)
for i in range(len(timeax)): angle[i] = get_angle(Bvecs[i],Bvecs_mean[i])
angle_reduced = filter_angle(angle,low,high) #If you want all angles, do 0,180
angle=angle_reduced

# Make Tplot variables
store_alldata()
machplot()
quickplot()
#%%

bins=100
fig,ax = plt.subplots(1,2,figsize=(10,5))
ax[0].hist(np.log10(ma_mean),bins=bins)
ax[0].set_xlim(-1,1)
ax[0].set_xlabel('Log10(Ma)')
ax[0].set_ylabel('Counts')
ax[1].hist(angle,bins=bins)
ax[1].set_xlim(0,180)
ax[1].set_xlabel('Deflection Angle (degrees)')
ax[1].set_ylabel('Counts')






#%%
fig, ax = plt.subplots(1,3,figsize=(12,4))
ax[0].hist2d(Br/Bmag_mean,vr/vmag_mean,range=[[-2,2],[-0.5,1.5]],bins=100)
ax[0].set_xlabel('Br/|B|')
ax[0].set_ylabel('vr/|v|')

ax[1].hist2d(Bt/Bmag_mean,vt/vmag_mean,bins=100)
ax[1].set_xlabel('Bt/|B|')
ax[1].set_ylabel('vt/|v|')


ax[2].hist2d(Bn/Bmag_mean,vn/vmag_mean,bins=100)
ax[2].set_xlabel('Bn/|B|')
ax[2].set_ylabel('vn/|v|')

for i in range(3):
	ax[i].axvline(x=1)
	ax[i].axhline(y=1)

#%%

fig, ax = plt.subplots(1,2,figsize=(10,5))
ax[0].hist2d(abs(dB_par_norm),abs(dB_perp_norm),range=[[0,0.1],[0,0.1]],bins=70)
ax[0].set_xlabel('|dB_par|/|B|')
ax[0].set_ylabel('|dB_perp|/|B|')

ax[1].hist2d(abs(dv_par_norm),abs(dv_perp_norm),range=[[0,0.2],[0,0.2]],bins=70)
ax[1].set_xlabel('|dv_par|/|v|')
ax[1].set_ylabel('|dv_perp|/|v|')

for i in range(2):
	ax[i].axvline(x=1)
	ax[i].axhline(y=1)

#%%

fig, ax = plt.subplots(1,1,figsize=(8,8))
ax.hist2d(np.log10(ma_mean),angle,range=[[-2,2],[0,180]],bins=500,cmin=1,cmap='viridis')
ax.set_xlabel('Log10 Alfven Mach Number')
ax.set_ylabel('Deflection Angle (degrees)')
ax.axhline(y=90,color='w')
ax.axvline(x=0,color='w')
ax.set_xlim(-0.5,0.5)



#%%

plt.scatter(np.log10(ma_mean),angle,s=3)
plt.xlabel('Log10(Ma)')
plt.axhline(y=90,color='k')
plt.axvline(x=0,color='k')
plt.xlim(-1,1)


#%%
cm = plt.cm.get_cmap('seismic')
sc=plt.scatter(ma_mean,angle,c=S[:,0],s=0.11,cmap=cm,vmin=0,vmax=2)

plt.colorbar(sc,label="beta")
plt.xscale('log')
plt.xlim(0.1,10)
plt.ylim(0,180)
plt.axhline(y=90,c='k')
# plt.axhline(y=low,c='b',linestyle='dashed',linewidth=0.3)
plt.axvline(x=1,c='k')
plt.xlabel('Alfven Mach Number Ma')
plt.ylabel('Deflection Angle')
plt.show()
# np.random.seed(42)  # for reproducibility
# x = position
# y = ma_mean

# # Create the 2D histogram
# plt.hist2d(x, y, bins=(0, 50),cmap='seismic')

# # Add a color bar to show the density scale
# plt.colorbar(label='Count')

# # Add labels and title
# plt.xlabel('Sr')
# plt.ylabel('Kr')
# plt.xlim()
# plt.title('2D Histogram of X and Y')

# # Show the plot
# plt.show()
# # %%

# %%
