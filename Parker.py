#%%
import pyspedas
import numpy as np
import matplotlib.pyplot as plt
import pytplot
from pytplot import get_data, store_data,timespan,tplot_options
from pyspedas import tplot,tplot_options
from pyspedas import tinterpol
import cdflib

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
	pyspedas.tplot(['Mach','position','velocities'])
	return

def Tempplot():
	mach = np.zeros([len(timeax),2])
	for i in range(len(timeax)):
		mach[i] = [ma_mean[i],1]
	store_data('Mach', data = {'x':timeax,'y':mach})
	pyspedas.options('Tpar', 'ytitle', r'$T_\parallel$ (KeV)')
	pyspedas.options('Tperp', 'ytitle', r'$T_\perp$ (KeV)')
	pyspedas.options('Tparperp', 'ytitle', r'$T_\parallel/T_\perp$')
	# pyspedas.ylim('Tpar',0,2)
	pyspedas.ylim('Tparperp',0,2)
	pyspedas.tplot([B_name,'Tparperp','Tpar','Tperp'])
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
	char_size = 13
	mach = np.zeros([len(timeax),2])
	pos = np.zeros([len(timeax),2])
	ang = np.zeros([len(timeax),2])
	vels = np.zeros([len(timeax),2])

	for i in range(len(timeax)):
		mach[i] = [ma_mean[i],1]
		pos[i] = [position[i],10]
		ang[i] = [angle[i],90]
		vels[i] = [1e-3*vmag_mean[i],1e-3*va_mean[i]]
	position[-1] = position[-2]
	store_data('Mach', data = {'x':timeax,'y':mach})
	store_data('Position',data = {'x':timeax,'y':position})
	store_data('Angle',data = {'x':timeax,'y':ang})
	store_data('Bvecs_mean', data = {'x':timeax,'y':1e9*Bvecs_mean})
	store_data('vvecs_mean', data = {'x':timeax,'y':1e-3*vivecs_mean})
	store_data('vmag_mean', data = {'x':timeax,'y':1e-3*vmag_mean})
	store_data('vmag', data = {'x':timeax,'y':1e-3*v_mag})
	store_data('vels', data = {'x':timeax,'y':vels})
	pyspedas.options('vels','ytitle', 'Velocities' + '\n [km/s]')
	pyspedas.options('vels','color', ['k','r'])
	pyspedas.options('vels','legend_names', [r'$\langle|v|\rangle$',r'$\langle v_a \rangle$'])
	pyspedas.options('Bvecs_mean','legend_names', [r'$R$',r'$T$',r'$N$'])
	pyspedas.options('vvecs_mean','legend_names', [r'$R$',r'$T$',r'$N$'])
	pyspedas.options(B_name,'ytitle', r'$\vec{B}$' + '\n [nT]')
	pyspedas.options('Bvecs_mean', 'ytitle', r'$\vec{\langle B\rangle}$' + '\n [nT]')
	pyspedas.options('vvecs_mean', 'ytitle', r'$\vec{\langle v\rangle}$' + '\n [km/s]')
	pyspedas.options('vmag_mean', 'ytitle', r'$\langle |v| \rangle$' + '\n [km/s]')
	pyspedas.options('Mach', 'ytitle', r'$M_a$')
	pyspedas.options('Mach', 'ylog', 1)
	pyspedas.options('Mach', 'linestyle', ['-','--'])
	pyspedas.options('Mach', 'color', ['k','r'])
	pyspedas.options('Position', 'linestyle', ['-','--'])
	pyspedas.options('Position', 'color', ['k','r'])
	pyspedas.options('Angle', 'linestyle', ['-','--'])
	pyspedas.options('Angle', 'color', ['k','r'])
	pyspedas.ylim('Mach',0.3,3)
	pyspedas.ylim('Bvecs_mean',-1000,2500)
	pyspedas.ylim('vvecs_mean',-200,500)
	pyspedas.ylim('Angle',0,120)
	pyspedas.ylim('vels',0,1000)
	pyspedas.ylim('Position',30,50)
	pyspedas.options('Angle','ytitle',r'$\theta$'+'\n [degrees]')
	pyspedas.options('Position','ytitle',r'$R/R_{sun}$')
	tplot_options('axis_font_size',17)
	tplot_options('charsize',13)
	tplot_options('aspect',13)
	pyspedas.options(['Position','Mach','Angle','Bvecs_mean','vvecs_mean','vels'],'char_size',17)
	pyspedas.tplot(['Position','Bvecs_mean','vvecs_mean','Mach','Angle'])
def betaplot():
	bpar = np.zeros([len(timeax),2])
	ang = np.zeros([len(timeax),2])
	Tratio = np.zeros([len(timeax),2])
	for i in range(len(timeax)):
		bpar[i] = [beta_par[i],1]
		ang[i] = [angle[i],90]
		Tratio[i] = [abs(Tperp[i]/Tpar[i]),1]
	store_data('Angle',data = {'x':timeax,'y':ang})
	store_data('Position',data = {'x':timeax,'y':position})
	store_data('bpar',data={'x':timeax,'y':bpar})
	store_data('Tratio',data={'x':timeax,'y':Tratio})
	pyspedas.options('bpar', 'color',['k','k'])
	pyspedas.ylim('bpar',0,12)
	pyspedas.options('beta', 'linestyle',['-','--'])
	pyspedas.tplot([B_name,'Position','Tratio','bpar','Angle'])
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


def FOV_filter():
	upper = phis[0,1]
	lower = phis[0,-2]

	peak_phi = phis[np.arange(len(timeax)),np.argmax(ephi,axis=0)]
	mask = (peak_phi>=upper) | (peak_phi<=lower)
	Ti[mask],Tpar[mask],Tperp[mask] = np.nan,np.nan,np.nan
	return 

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
	# pth = kb*n*T
	pth = n*T
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
	Ppar = n*Tpar
	Pperp = n*Tperp
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

# Needed to change get_H.  the Pressure tensor data product is expressed as a vector w 6 elements (not a tensor structure)
# So Manually get the trace of the tensor.  
# Simple version for now:  No v.P for now
def get_H(v,P):
	traceP = P[0]+P[1]+P[2]
	vdPx = v[0]*(P[0]+P[3]+P[4])
	vdPy = v[1]*(P[1]+P[3]+P[5])
	vdPz = v[2]*(P[2]+P[4]+P[5])
	vdotP = np.array([vdPx,vdPy,vdPz])
	H = 0.5*v*traceP + vdotP
	return H
def get_S(E,B):
	S = np.cross(E,B)/mu0
	return S

# mean, vector mean, and delta vector (a-a_mean)
def get_mean(var,int):   #Take a timeseries and compute the mean (basically smooth outsmall flucs over some interval)
	box = np.ones(int)/int
	smoothed_var = np.convolve(var,box,mode='same')
	return smoothed_var


def get_mean(var, int):
    """Take a timeseries and compute the mean (basically smooth out small flucs over some interval)
    Handles NaN values by computing local averages only from valid data points.
    """
    import numpy as np
    from scipy.ndimage import uniform_filter1d
    
    # Create a copy to avoid modifying original
    var_copy = np.array(var, dtype=float)
    
    # Create mask for valid (non-NaN) values
    valid_mask = ~np.isnan(var_copy)
    
    # Replace NaNs with 0 for convolution
    var_filled = np.where(valid_mask, var_copy, 0)
    
    # Convolve the data and the mask
    box = np.ones(int) / int
    smoothed_sum = np.convolve(var_filled, box, mode='same')
    smoothed_count = np.convolve(valid_mask.astype(float), box, mode='same')
    
    # Divide by actual count of valid points in each window
    # Avoid division by zero
    smoothed_var = np.where(smoothed_count > 0, smoothed_sum / smoothed_count, np.nan)
    
    return smoothed_var


def get_med(var,int):
	from scipy import ndimage
	return ndimage.median_filter(var, size=int, mode='nearest')

def get_vecmean(vec,interval):   # Vector mean
	vec1_mean = get_mean(vec[:,0],interval)
	vec2_mean = get_mean(vec[:,1],interval)
	vec3_mean = get_mean(vec[:,2],interval)
	vec_mean = np.zeros_like(vec)
	for i in range(len(vec_mean)): vec_mean[i] = np.array([vec1_mean[i],vec2_mean[i],vec3_mean[i]])
	return vec_mean


def get_vecmed(vec,int):   # Vector mean
	vec1_mean = get_med(vec[:,0],minutes*meaninterval)
	vec2_mean = get_med(vec[:,1],minutes*meaninterval)
	vec3_mean = get_med(vec[:,2],minutes*meaninterval)
	vec_mean = np.zeros_like(vec)
	for i in range(len(vec_mean)): vec_mean[i] = np.array([vec1_mean[i],vec2_mean[i],vec3_mean[i]])
	return vec_mean


def get_delta(vec,vec_mean):
	dvec = vec - vec_mean
	dvec_norm = dvec/np.linalg.norm(vec_mean)
	return dvec,dvec_norm


def get_deltascalar(var,var_mean):
	dvar = var - var_mean
	dvar_norm = dvar/var_mean
	return dvar,dvar_norm

# Why are there so many different versions of cross helicity?  
def get_crosshelicity(v,B,n,m): #vector v & B
	z_plus = v + B/np.sqrt(n*m*mu0)
	z_minus = v - B/np.sqrt(n*m*mu0)
	term1 = np.linalg.norm(z_plus)**2 - np.linalg.norm(z_minus)**2
	term2 = np.linalg.norm(z_plus)**2 + np.linalg.norm(z_minus)**2
	sigma_c = term1/term2
	return sigma_c

#Somehow always = 1 (need to resolve)
def get_residenergy(v,B,n,m): # vector dv & dB (Alfven units??)
	# term1 = np.linalg.norm(dv)**2 - np.linalg.norm(dB)**2
	# term2 = np.linalg.norm(dv)**2 + np.linalg.norm(dB)**2
	term1 = np.linalg.norm(dv)**2 - np.linalg.norm(dB/np.sqrt(n*m*mu0))**2
	term2 = np.linalg.norm(dv)**2 + np.linalg.norm(dB/np.sqrt(n*m*mu0))**2
	sigma_r = term1/term2
	return sigma_r 


def store_alldata():

	store_data('Bvecs_mean', data = {'x':timeax,'y':1e9*Bvecs_mean})
	store_data('vmag_mean', data = {'x':timeax,'y':1e-3*vmag_mean})
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
	store_data('T', data = {'x':timeax,'y':Ti/1.602e-19})
	store_data('Tpar', data = {'x':timeax,'y':Tpar/1.602e-19})
	store_data('Tperp', data = {'x':timeax,'y':Tperp/1.602e-19})
	store_data('Tparperp', data = {'x':timeax,'y':Tpar/Tperp})

	#----------------------

	store_data('velocities', data = {'x':timeax,'y':viandva})
	store_data('|vi|', data = {'x':timeax,'y':va})
	store_data('|va|', data = {'x':timeax,'y':v_mag})
	return

# %%

# %%

# %%

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
# trange = ['2021-08-09/20:00', '2021-08-10/05:00'] # Encounter 9 (some sub-Alfvenic)
#trange = ['2022-02-25', '2022-02-28'] #Dudok de wit 2020 Full interval

#

# sub_alfs = 	[['2021-08-09/21:30','2021-08-10/00:00'], # 2.5 hrs Encounter 9   
# 			 ['2021-11-21/21:00','2021-11-22/01:00'], # 4 hrs Encounter 10    
# 			 ['2021-11-22/03:30','2021-11-22/10:00'], # 6.5 hrs Encounter 10  
# 			 ['2022-02-25/20:00','2022-02-25/23:30'], # 3.5 hrs Encounter 11  
# 			 ['2022-06-01/18:00','2022-06-02/08:00'], # 14 hrs Encounter 12 

# sup_alfs = 	[['2021-08-10/00:30','2021-08-10/06:00'], # 5.5 hrs Encounter 9
# 			 ['2021-11-21/16:00','2021-11-21/21:00'], # 5 hr:w
# Encounter 10
# 			 ['2021-11-22/01:00','2021-11-22/02:30'], # 1.5 hrs Encounter 10
# 			 ['2022-02-26/07:30','2022-02-26/08:30'], # 1 hr Encounter 11
# 			 ['2022-06-02/13:00','2022-06-02/20:00'], # 7 hrs Encounter 12



sub_alfs =  [['2022-09-06/06:00','2022-09-06/16:00'], # 10 hrs Encounter 13
			 ['2022-09-06/18:00','2022-09-07/12:00'], # 18 hrs Encounter 13
			 ['2022-12-11/04:00','2022-12-11/16:00'], # 12 hrs Encounter 14
			 ['2023-03-16/12:00','2023-03-17/06:00'], # 18 hrs Encounter 15
			 ['2023-06-20/01:00','2023-06-21/01:00'], # 24 hrs Encounter 16
			 ['2023-09-27/06:00','2023-09-27/15:00'], # 9 hrs Encounter 17
			 ['2023-12-29/04:00','2023-12-29/14:00'], # 10 hrs Encounter 18
			 ['2024-03-29/06:00','2024-03-29/21:00'], # 15 hrs Encounter 19
			 ['2024-06-30/03:00','2024-06-30/18:00'], # 15 hrs Encounter 20
			 ['2024-09-28/11:00','2024-09-28/18:00'], # 7 hrs Encounter 21
 		     ['2024-12-24/00:00', '2024-12-25/00:00'], # 24 hours Encounter 22 
	         ['2025-03-22/06:00', '2025-03-23/12:00'] # 30 hrs Encounter 23
			 ]


sup_alfs=	[['2022-09-07/18:00','2022-09-08/18:00'], # 24 hrs Encounter 13
			 ['2022-12-10/01:00','2022-12-10/09:00'], # 8 hrs Encounter 14
			 ['2023-03-16/00:00','2023-03-16/11:00'], # 11 hrs Encounter 15
			 ['2023-06-24/00:00','2023-06-25/00:00'], # 24  hrs Encounter 16
			 ['2023-09-28/20:00','2023-09-30/09:00'], # 37 hrs Encounter 17
			 ['2023-12-25/14:00','2023-12-25/23:00'], # 9 hrs Encounter 18
			 ['2024-04-01/10:00','2024-04-02/01:00'], # 15 hrs Encounter 19
			 ['2024-07-01/20:00','2024-07-03/00:00'], # 28 hrs Encounter 20
			 ['2024-10-03/00:00','2024-10-03/12:00'], # 12 hrs Encounter 21
			 ['2025-03-21/00:00', '2025-03-22/03:00'], # 27 hrs Encounter 23
			 ]


#May need to re-think the definition of near-alfs.  May artificially produce many >90 deflections
near_alfs_old = [['2022-09-05/11:00','2022-09-05/17:00'], # 6 hrs Encounter 13
			 ['2023-06-23/00:00','2023-06-23/18:00'], # 18 hrs Encounter 16
			 ['2023-09-28/09:00','2023-09-28/18:00'], # 9 hrs Encounter 17
			 ['2023-12-28/00:00','2023-12-29/00:00'], # 24 hrs Encounter 18
			 ['2024-03-28/22:00','2024-03-29/05:00'], # 7 hrs Encounter 19
			 ['2024-06-29/00:00','2024-06-29/09:00'], # 9 hrs Encounter 20 
			 ['2024-10-01/15:00','2024-10-02/00:00'], # 9 hrs Encounter 21 
]


near_alfs = [['2022-09-05/11:00','2022-09-05/17:00'], # 6 hrs Encounter 13
			 ['2023-09-28/09:00','2023-09-28/18:00'], # 9 hrs Encounter 17
			 ['2023-12-27/06:00','2023-12-28/00:00'], # 18 hrs Encounter 18 (done)
			 ['2024-03-28/22:00','2024-03-29/03:00'], # 7 hrs Encounter 19
			 ['2024-06-29/00:00','2024-06-29/09:00'], # 9 hrs Encounter 20 
			 ['2024-09-30/13:00','2024-10-01/00:00'], # 9 hrs Encounter 21 
]

new_near_sups = [
]



# 
more_near = []
# Collecting near Alfvenic intervals in the |Log(Ma)| ~ 0-0.15 range (max 0.2)


near_sups = [['2024-09-27/04:30','2024-09-27/06:30'], # 2 hr Encounter 21
			 ['2024-10-03/23:00','2024-10-04/00:00'], # 1 hr Encounter 21
			 ['2024-10-04/09:00','2024-10-04/10:00'], # 1 hr Encounter 21
			 ['2024-06-29/04:45','2024-06-29/06:15'], # 1.5 hr Encounter 20
			 ['2024-07-02/03:00','2024-07-02/08:00'], # 5 hr Encounter 20
			 ['2024-03-28/03:30','2024-03-28/06:00'], # 2.5 hr Encounter 19
			 ['2023-12-29/22:00','2023-12-30/03:00'], # 5 hrs Encounter 18
			 ['2023-09-29/00:00','2023-09-29/06:00'], # 6 hrs Encounter 17
			 ['2023-09-28/18:00','2023-09-28/21:00'], # 3 hrs Encounter 17
			 ]

near_subs = [['2024-09-27/07:30','2024-09-27/8:30'], # 1 hr Encounter 21
			 ['2024-09-30/23:00','2024-10-01/01:00'],# 2 hr Encounter 21
			 ['2024-06-29/16:00','2024-06-29/21:00'],# 5 hr Encounter 20
			 ['2024-03-30/07:00','2024-03-30/11:00'], # 4 hr Encounter 19
			 ['2023-12-28/12:00','2023-12-28/15:00'], # 3 hrs Encounter 18
			 ['2023-09-27/09:00','2023-09-29/18:00'], # ? hrs Encounter 17
			 ]


sup_alfs_alt=[['2022-09-07/18:00','2022-09-08/18:00'], # 24 hrs Encounter 13
			 ['2022-12-10/01:00','2022-12-10/09:00'], # 8 hrs Encounter 14
			 ['2023-03-16/06:00','2023-03-16/11:00'], # 5 hrs Encounter 15
			 ['2023-06-24/09:00','2023-06-25/00:00'], # 15  hrs Encounter 16
			 ['2023-09-28/20:00','2023-09-30/09:00'], # 37 hrs Encounter 17
			 ['2023-12-25/14:00','2023-12-25/23:00'], # 9 hrs Encounter 18
			 ['2024-04-01/10:00','2024-04-02/01:00'], # 15 hrs Encounter 19
			 ['2024-07-02/03:00','2024-07-03/00:00'], # 21 hrs Encounter 20
			 ['2024-10-03/00:00','2024-10-03/12:00'], # 12 hrs Encounter 21
			 ['2025-03-21/00:00', '2025-03-22/00:00'], # 27 hrs Encounter 23
			 ]

#sub_alfs =  ['2024-09-30/03:00','2024-09-30/11:00']
# trange=sub_alfs[]

# Alfven crossings (<2 hr)
#trange = ['2021-11-21/21:00','2021-11-21/22:00']
# trange = ['2021-08-10/00:15', '2021-08-10/00:45'] # Encounter 9 (some sub-Alfvenic)
# trange=['2023-12-29/01:30','2023-12-29/03:00'] # E18
#trange=['2024-03-29/22:00','2024-03-29/23:30']  # E19 
#trange = ['2024-06-29/11:00', '2024-06-29/13:00'] # E20
# trange = ['2024-09-28/00:00', '2024-10-05/12:00'] # E21 


trange = ['2021-08-12/00:30', '2021-08-12/02:30'] # Soni 2024 Parker interval
#trange = ['2022-09-04/00:00','2022-09-04/12:00'] 
#trange = ['2024-09-28/00:00', '2024-10-05/12:00'] # E21 
# trange = ['2024-06-29/00:00', '2024-07-01/12:00'] # E20
# trange=['2024-03-27/00:00','2024-04-02/00:00']  # E19 
# trange=['2023-12-28/00:00','2023-12-30/00:00'] # E18


#trange = ['2018-11-05/00:00', '2018-11-05/03:00'] # Bale 2019 event (includes Sr)
# trange = near_alfs[0]
#trange = ['2024-12-22/00:00', '2024-12-22/12:00'] # E22 
#trange = ['2025-03-21/10:00', '2025-03-21/14:00'] # E23
# trange = ['2025-06-19/00:00', '2025-06-20/00:00'] # E24
#trange = ['2021-04-26/00:00', '2021-05-02/00:00'] 
#t#range = ['2023-06-20/00:00','2023-06-23/00:00']


# Intervals with some beta_par>1
trange = ['2024-03-25/04:00','2024-03-25/09:00'] #Enc 19
trange = ['2024-07-04/22:00','2024-07-05/00:00'] #Enc 20 #HCS just before 
trange = ['2024-09-25/06:00','2024-09-26/06:00'] #Enc 21
trange = ['2024-12-19/00:00', '2024-12-20/00:00'] #Enc 22
trange = ['2025-03-27/00:00', '2025-03-28/00:00'] # Maybe
trange = ['2023-06-24/00:00','2023-06-24/12:00'] # Maybe

betaparlist = [['2024-03-25/04:00','2024-03-25/09:00'],
			   ['2024-07-04/22:00','2024-07-05/00:00'],
			   ['2024-09-25/06:00','2024-09-26/06:00'],
			#    ['2024-09-25/06:00','2024-09-26/12:00'],
			   ['2024-12-19/00:00', '2024-12-20/00:00'],
			   ['2025-03-18/00:00', '2025-03-19/00:00'],
			   ['2025-03-17/00:00', '2025-03-17/12:00'],
			   ['2023-12-23/12:00','2023-12-24/12:00']]

# trange=sup_alfs_alt[3]
# trange = ['2024-07-04/21:00','2024-07-05/00:00'] #Enc20 #Maybe brief HCS
#Need to make a list of chosen Alfvenic & sub-alfvenic intervals
#Most people choose a handful of intervals by eye
# trange=sup_alfs[9]
trange = ['2022-09-03/00:00','2022-09-10/00:00']
trange = ['2023-09-28/09:00','2023-09-28/18:00']
trange = betaparlist[2]
trange=['2022-09-05/00:00','2022-09-05/12:00']

# longints = [['2025-06-15/00:00','2025-06-23/00:00'],['2025-03-18/00:00','2025-03-27/00:00']]
Bfld_vars = pyspedas.projects.psp.fields(trange=trange, level='l2', time_clip=True)
spi_vars = pyspedas.projects.psp.spi(trange=trange,level='l3',time_clip=True,get_support_data=True)
# spe_vars = pyspedas.projects.psp.spe(trange=trange,level='l2',time_clip=True)

# qtn_vars = pyspedas.projects.psp.fields(trange=trange,level='l3',datatype='sqtn_rfs_V1V2',time_clip=True)

# spe_vars = pyspedas.projects.psp.spi(trange=trange,level='l3',time_clip=True)
#spe_vars = pyspedas.projects.psp.spe(trange=trange,level='l2',datatype='sqtn_rfs_V1V2',time_clip=True)

#voltages_vars = pyspedas.projects.psp.fields(trange=trange, datatype='dfb_wf_dvdc', level='l2',time_clip=True)
#On DC datatype: 'sqn_rfs_V1V2 has some kind of electron density & core temp, but looks weird
# spc_vars = pyspedas.projects.psp.spc(trange=trange, datatype='l2', level='l2')
##%%1G



# Reform all data to simple arrays and convert to SI units 
B_name = 'psp_fld_l2_mag_RTN'
vi_name = 'psp_spi_VEL_RTN_SUN'
Bxyz_name = 'psp_spi_MAGF_INST'
vxyz_name = 'psp_spi_VEL_INST'
TiTensor_name = 'psp_spi_T_TENSOR_INST'
Ti_name = 'psp_spi_TEMP'
phivals_name = 'psp_spi_PHI_VALS'
ephi_name = 'psp_spi_EFLUX_VS_PHI'
ni_name = 'psp_spi_DENS'
# ni_name = 'electron_density'
# voltages_name = 'psp_fld_l2_dfb_wf_dVdc_sc'
position_name = 'psp_spi_SUN_DIST'

interpvar_name = vi_name
timeax = pytplot.get_data(interpvar_name).times
meaninterval = mean_int(timeax,trange)


for name in [B_name, Bxyz_name, vxyz_name, vi_name, Ti_name, ni_name, TiTensor_name, position_name]:
    data = pytplot.get_data(name)
    if data is None:
        print(f"WARNING: no data found for {name}")
        continue
    times = data.times
    _, unique_idx = np.unique(times, return_index=True)
    if len(unique_idx) < len(times):
        print(f"Removing {len(times)-len(unique_idx)} duplicate timestamps from {name}")
        pytplot.store_data(name, data={'x': times[unique_idx], 'y': data.y[unique_idx]})


tinterpol(B_name,interpvar_name,newname='B')
tinterpol(Bxyz_name,interpvar_name,newname='Bxyz')
tinterpol(vxyz_name,interpvar_name,newname='vxyz')
tinterpol(vi_name,interpvar_name,newname='vi')
tinterpol(Ti_name,interpvar_name,newname='Ti')
tinterpol(ni_name,interpvar_name,newname='ni')
tinterpol(TiTensor_name,interpvar_name,newname='TiTensor')
tinterpol(phivals_name,interpvar_name,newname='phivals')
tinterpol(ephi_name,interpvar_name,newname='ephi')
# tinterpol(voltages_name,interpvar_name,newname='voltages')
tinterpol(position_name,interpvar_name,newname='position')
##%%
Bvecs = 1e-9*reform(pytplot.get_data('B'))
Bxyz = 1e-9*reform(pytplot.get_data('Bxyz'))
vxyz = 1e3*reform(pytplot.get_data('vxyz'))
vivecs = 1e3*reform(pytplot.get_data('vi'))
ni = 1e6*reform(pytplot.get_data('ni'))
Ti = 1.602e-19*reform(pytplot.get_data('Ti'))
TiTensor = 1.602e-19*reform(pytplot.get_data('TiTensor'))
phis = reform(pytplot.get_data('phivals'))
ephi = reform(pytplot.get_data('ephi')).T
#Ti = reform(pytplot.get_data('Ti'))/kb # Kelvin Units
#TiTensor = reform(pytplot.get_data('TiTensor'))/kb # Kelvin Units
# voltages = reform(get_data('voltages'))
position = reform(get_data('position'))/695700 #Solar radii

##%%

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
beta_par = np.zeros_like(ni)

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
	# H[i] = get_H(vivecs[i],ni[i]*TiTensor[i])
	H[i] = get_H(vxyz[i],ni[i]*TiTensor[i])

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
	beta_par[i] = Ppar[i]/P_mag[i]


# FOV filtering
FOV_filter()



##%%
# Get averaged quantities.  meaninterval = 1 min ##
minutes = 60

# Means
Bvecs_mean = get_vecmean(Bvecs,minutes*meaninterval)
vivecs_mean = get_vecmean(vivecs,minutes*meaninterval)
Bmag_mean = get_mean(B_mag,minutes*meaninterval)
vmag_mean = get_mean(v_mag,minutes*meaninterval)
n_mean = get_mean(ni,minutes*meaninterval)
va_mean = Bmag_mean/np.sqrt(mu0*mi*n_mean)
# va_mean = Bmag_mean/np.sqrt(mu0*mi*ni)
ma_mean = vmag_mean/va_mean
beta_mean = get_mean(beta[:,0], minutes*meaninterval)


# Medians
# Bvecs_mean = get_vecmed(Bvecs,minutes*meaninterval)
# vivecs_mean = get_vecmed(vivecs,minutes*meaninterval)
# Bmag_mean = get_med(B_mag,minutes*meaninterval)
# vmag_mean = get_med(v_mag,minutes*meaninterval)
# n_mean = get_med(ni,minutes*meaninterval)
# va_mean = Bmag_mean/np.sqrt(mu0*mi*n_mean)
# ma_mean = vmag_mean/va_mean
# beta_mean = get_med(beta[:,0], minutes*meaninterval)

# Get dB,dv + components, etc. and vB alignment variable (proxy for alfvenicity)
dB,dB_norm = np.zeros_like(Bvecs),np.zeros_like(Bvecs)
dv,dv_norm = np.zeros_like(vivecs),np.zeros_like(vivecs)
dB_par,dB_perp = np.zeros_like(ni),np.zeros_like(ni)
dv_par,dv_perp = np.zeros_like(ni),np.zeros_like(ni)
B_par = np.zeros_like(ni)
v_par = np.zeros_like(ni)
dn,dn_norm = np.zeros_like(ni),np.zeros_like(ni)
dB_par_norm = np.zeros_like(ni)
dB_perp_norm = np.zeros_like(ni)
dv_par = np.zeros_like(ni)
dv_par_norm = np.zeros_like(ni)
dv_perp_norm = np.zeros_like(ni)
dB_norm_mag,dv_norm_mag = np.zeros_like(ni),np.zeros_like(ni)
vB_alignment = np.zeros_like(ni)
dvdB_alignment = np.zeros_like(ni)
angle = np.zeros_like(ni)
sigma_c = np.zeros_like(ni)
sigma_r = np.zeros_like(ni)

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
	

	# dn
	dn[i],dn_norm[i] = get_deltascalar(ni[i],n_mean[i])
	sigma_r[i] = get_residenergy(vivecs[i],Bvecs[i],ni[i],mi)
	sigma_c[i] = get_crosshelicity(vivecs[i],Bvecs[i],ni[i],mi)
	
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


low = 0
high = 180
angle_reduced = np.zeros_like(angle)
for i in range(len(timeax)): angle[i] = get_angle(Bvecs[i],Bvecs_mean[i])
angle_reduced = filter_angle(angle,low,high) #If you want all angles, do 0,180
angle=angle_reduced





spi_vars

# Make Tplot variables
store_alldata()
# betaplot()
# angleplot()
# angleplot()
# quickplot()
#betaplot()
#%%
fig,ax = plt.subplots(2,1,figsize=(10,5))
ax[0].plot(np.log10(ma_mean))
ax[0].axhline(y=0,color='k')
ax[1].plot(angle)
# %%
fluxplot()
# %%
fig, ax = plt.subplots(2,1,figsize=(10,10))
ax[0].plot(ma_mean)
ax[0].plot(ma)
ax[0].axhline(y=1,color='k',linestyle='dashed')
ax[1].plot(angle)
ax[1].axhline(y=90,color='k',linestyle='dashed')
# %%

 #%%

bins=30
fig,ax = plt.subplots(1,2,figsize=(10,5))
ax[0].hist(np.log10(ma_mean),bins=bins)
ax[0].set_xlim(-1,1)
ax[0].set_xlabel('Log10(Ma)')
ax[0].set_ylabel('Counts')
ax[0].axvline(x=0,color='k')
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




# %%
