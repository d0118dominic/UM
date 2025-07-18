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
def quickplot():
	pyspedas.tplot([B_name,vi_name,Ti_name,ni_name])
	return
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
# Encounters below 0.05 AU



# Encounters below 0.08 AU 


# Encounters below 0.05 AU
# Encounters below 0.06 AU
# Encounters below 0.07 AU



# trange = ['2024-10-29/00:00', '2024-10-31/00:00']
#trange = ['2018-11-1', '2018-11-10'] #Dudok de wit 2020 Full interval
#trange = ['2021-04-28/00:00', '2021-04-30/00:00'] # Encounter 8 (some sub-Alfvenic)
#trange = ['2021-08-09/12:00', '2021-08-10/00:00'] # Encounter 9 (some sub-Alfvenic)
#trange = ['2021-08-09/20:00', '2021-08-10/05:00'] # Encounter 9 (some sub-Alfvenic)
#trange = ['2022-02-25', '2022-02-28'] #Dudok de wit 2020 Full interval

#trange = ['2018-11-05/00:00', '2018-11-05/03:00'] # Bale 2019 event (includes Sr)
#


sub_alfs =  [['2022-09-06/06:00','2022-09-06/16:00'], # 10 hrs Encounter 13
			 ['2022-09-06/18:00','2022-09-07/12:00'], # 18 hrs Encounter 13
			 ['2022-12-11/00:00','2022-12-11/18:00'], # 18 hrs Encounter 14
			 ['2023-03-16/12:00','2023-03-17/06:00'], # 18 hrs Encounter 15
			 ['2023-06-20/01:00','2023-06-21/01:00'], # 24 hrs Encounter 16
			 ['2023-09-27/06:00','2023-09-27/15:00'], # 9 hrs Encounter 17
			 ['2023-12-29/04:00','2023-12-29/14:00'], # 10 hrs Encounter 18
			 ['2024-03-29/06:00','2024-03-29/21:00'], # 15 hrs Encounter 19
			 ]


sup_alfs=	[['2022-09-07/18:00','2022-09-08/18:00'], # 24 hrs Encounter 13
			 ['2022-12-10/01:00','2022-12-10/09:00'], # 8 hrs Encounter 14
			 ['2023-03-16/00:00','2023-03-16/11:00'], # 11 hrs Encounter 15
			 ['2023-06-24/00:00','2023-06-25/00:00'], # 24  hrs Encounter 16
			 ['2023-09-30/00:00','2023-09-30/09:00'], # 9 hrs Encounter 17
			 ['2023-12-25/14:00','2023-12-25/23:00'], # 9 hrs Encounter 18
			 ['2024-04-01/10:00','2024-04-02/01:00'], # 15 hrs Encounter 19
			 ]

#sub_alfs =  ['2024-09-30/03:00','2024-09-30/11:00']

recent_perihelia = [['2024-09-29/00:00', '2024-10-01/12:00'], #E21
					['2024-06-29/00:00', '2024-07-01/12:00'], #E20
					['2024-03-29/00:00','2024-03-31/00:00'], #E19
					['2023-12-28/00:00','2023-12-30/00:00']] #E18

# trange = ['2021-08-11/09:00', '2021-08-12/09:00'] # Soni 2024 Parker interval
# trange = ['2024-09-30/00:00', '2024-09-30/23:59'] # E21 
#trange = ['2024-12-24/00:00', '2024-12-25/00:00'] # E22 
# trange = ['2025-03-22/00:00', '2025-03-23/00:00'] # E23
# trange = ['2025-06-19/00:00', '2025-06-20/00:00'] # E24

#Need to make a list of chosen Alfvenic & sub-alfvenic intervals
#Most people choose a handful of intervals by eye


# eventlist = subs

# Alfven crossings (<2 hr)
c1= ['2021-11-21/21:00','2021-11-21/22:00']
c2= ['2021-08-10/00:15', '2021-08-10/00:45'] # Encounter 9 (some sub-Alfvenic)
c3=['2023-12-29/01:30','2023-12-29/03:00'] # E18
c4=['2024-03-29/22:00','2024-03-29/23:30']  # E19 
c5= ['2024-06-29/11:00', '2024-06-29/13:00'] # E20
eventlist=[c1,c2,c3,c4,c5]
eventlist=sup_alfs+sub_alfs



#%%
allangles=np.array([])
allpositions=np.array([])
allmachs=np.array([])
allSr=np.array([])
allKr=np.array([])
allmachs=np.array([])
allbetas=np.array([])
alldB = np.array([])
alldv = np.array([])
alldvpar = np.array([])
alldvperp = np.array([])
alldBpar = np.array([])
alldBperp = np.array([])

for i in range(len(eventlist)):
	trange=eventlist[i]
	# trange=trange
	Bfld_vars = pyspedas.projects.psp.fields(trange=trange, level='l2', time_clip=True)
	swp_vars = pyspedas.projects.psp.spi(trange=trange,level='l3',time_clip=True)

	#voltages_vars = pyspedas.projects.psp.fields(trange=trange, datatype='dfb_wf_dvdc', level='l2',time_clip=True)
	#On DC datatype: 'sqn_rfs_V1V2 has some kind of electron density & core temp, but looks weird
	# spc_vars = pyspedas.projects.psp.spc(trange=trange, datatype='l2', level='l2')

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
	
	Bvecs = 1e-9*reform(pytplot.get_data('B'))
	Bxyz = 1e-9*reform(pytplot.get_data('Bxyz'))
	vxyz = 1e3*reform(pytplot.get_data('vxyz'))
	vivecs = 1e3*reform(pytplot.get_data('vi'))
	ni = 1e6*reform(pytplot.get_data('ni'))
	Ti = 1.602e-19*reform(pytplot.get_data('Ti'))
	TiTensor = 1.602e-19*reform(pytplot.get_data('TiTensor'))
	# voltages = reform(get_data('voltages'))
	position = reform(get_data('position'))/695700 #Solar radii



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


	low = 45
	high = 135
	angle_reduced = np.zeros_like(angle)
	for i in range(len(timeax)): angle[i] = get_angle(Bvecs[i],Bvecs_mean[i])
	angle_reduced = filter_angle(angle,low,high) #If you want all angles, do 0,180

	allangles = np.concat([allangles,angle_reduced])
	allmachs = np.concat([allmachs,ma_mean])
	allpositions = np.concat([allpositions,position])
	allSr = np.concat([allSr,S[:,0]])
	allKr = np.concat([allKr,K[:,0]])
	allbetas = np.concat([allbetas,beta_mean])
	

# %%
plt.plot(allpositions)
# plt.plot(allangles)
# %%


#%%
bins=100
fig,ax = plt.subplots(1,3,figsize=(18,6))
ax[0].hist(np.log10(allmachs),bins=bins,range=[-1.5,1.5],log=False)
ax[0].set_xlim(-1.5,1.5)
ax[0].set_xlabel('Log10(Ma)')
ax[0].set_ylabel('Counts')
ax[0].axvline(x=0,color='k')
ax[1].hist(allpositions,bins=bins,range=[0,50],log=False)
ax[1].set_xlim(10,45)
ax[1].set_xlabel('R (solar Radii)')
ax[1].set_ylabel('Counts')

ax[2].hist(allangles,bins=bins,range=[0,180],log=True)
ax[2].set_xlim(0,180)
ax[2].set_xlabel('Deflection Angle (degrees)')
ax[2].set_ylabel('Counts')


#%%

plt.scatter(np.log10(allmachs),allangles,s=1)
plt.xlabel('Log10(Ma)')
plt.ylabel('Deflection Angle (degrees)')
plt.axhline(y=90,color='k')
plt.axvline(x=0,color='k')
plt.xlim(-0.8,0.8)
#%%

fig, ax = plt.subplots(1,1,figsize=(8,8))
ax.hist2d(np.log10(allmachs),allangles,range=[[-0.5,0.5],[0,180]],bins=10,cmin=1,cmap='viridis')
ax.set_xlabel('Log10 Alfven Mach Number')
ax.set_ylabel('Deflection Angle (degrees)')
ax.set_facecolor('k')
ax.axhline(y=90,color='w')
ax.axvline(x=0,color='w')
ax.set_xlim(-0.5,0.5)


#%%
# cm = plt.cm.get_cmap('seismic')
# sc=plt.scatter(allmachs,allangles,s=1)
cm = plt.cm.get_cmap('seismic')
sc=plt.scatter(allmachs,allangles,c=allbetas,s=1,cmap=cm,vmin=0,vmax=2)

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
# %%

cm = plt.cm.get_cmap('seismic')
sc=plt.scatter(allpositions,allangles,s=1)


plt.ylim(0,180)
plt.xlim(9,20)
plt.axhline(y=90,c='k')
plt.axhline(y=low,c='b',linestyle='dashed',linewidth=0.3)
plt.axvline(x=1,c='k')
plt.xlabel('Position')
plt.ylabel('Deflection Angle')
plt.show()
# %%
