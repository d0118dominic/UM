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
# From E21,E9,E10,E11,E12,13,14,15
sub_alfs =  [['2024-09-30/03:00','2024-09-30/11:00'],
			 ['2021-08-09/21:30','2021-08-10/00:00'],
			 ['2021-11-21/21:00','2021-11-22/01:00'],
			 ['2021-11-22/03:30','2021-11-22/10:00'],
			 ['2022-02-25/20:00','2022-02-25/23:30'],
			 ['2022-06-01/18:00','2022-06-02/08:00'],
			 ['2022-09-06/06:00','2022-09-06/16:00'],
			 ['2022-09-06/18:00','2022-09-07/12:00'],
			 ['2022-12-11/00:00','2022-12-11/18:00'],
			 ['2023-03-16/12:00','2023-03-17/06:00'],
			 ]

sub_alfs=sub_alfs

allangles=np.array([])
allpositions=np.array([])
allmachs=np.array([])
allSr=np.array([])
allKr=np.array([])
allmachs=np.array([])
for i in range(len(sub_alfs)):

	trange = sub_alfs[i]
	Bfld_vars = pyspedas.projects.psp.fields(trange=trange, level='l2', time_clip=True)
	swp_vars = pyspedas.projects.psp.spi(trange=trange,level='l3',time_clip=True)

	B_name = 'psp_fld_l2_mag_RTN'
	vi_name = 'psp_spi_VEL_RTN_SUN'
	Bxyz_name = 'psp_spi_MAGF_INST'
	vxyz_name = 'psp_spi_VEL_INST'
	TiTensor_name = 'psp_spi_T_TENSOR_INST'
	Ti_name = 'psp_spi_TEMP'
	ni_name = 'psp_spi_DENS'
	voltages_name = 'psp_fld_l2_dfb_wf_dVdc_sc'
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
	tinterpol(voltages_name,interpvar_name,newname='voltages')
	tinterpol(position_name,interpvar_name,newname='position')

	Bvecs = 1e-9*reform(pytplot.get_data('B'))
	Bxyz = 1e-9*reform(pytplot.get_data('Bxyz'))
	vxyz = 1e3*reform(pytplot.get_data('vxyz'))
	vivecs = 1e3*reform(pytplot.get_data('vi'))
	ni = 1e6*reform(pytplot.get_data('ni'))
	Ti = 1.602e-19*reform(pytplot.get_data('Ti'))
	TiTensor = 1.602e-19*reform(pytplot.get_data('TiTensor'))
	voltages = reform(get_data('voltages'))
	position = reform(get_data('position'))

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

	E_antennas = np.zeros_like(voltages)

	for i in range(len(Bvecs)):
		Br[i], Bt[i], Bn[i] = Bvecs[i,0], Bvecs[i,1], Bvecs[i,2]
		vr[i], vt[i], vn[i] = vivecs[i,0], vivecs[i,1], vivecs[i,2]
		B_mag[i] = np.linalg.norm(Bvecs[i,:])
		v_mag[i] = np.linalg.norm(vivecs[i,:])
		
		E_conv[i] = -get_vxb(vivecs[i],Bvecs[i])
		E_antennas[i] = voltages[i]/3.5
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

	minutes = 10
	Bvecs_mean = get_vecmean(Bvecs,minutes*meaninterval)
	vivecs_mean = get_vecmean(vivecs,minutes*meaninterval)
	Bmag_mean = get_mean(B_mag,minutes*meaninterval)
	vmag_mean = get_mean(v_mag,minutes*meaninterval)
	n_mean = get_mean(ni,minutes*meaninterval)
	va_mean = Bmag_mean/np.sqrt(mu0*mi*n_mean)
	ma_mean = vmag_mean/va_mean
	beta_mean = get_mean(beta[:,0], minutes*meaninterval)



	ma = np.zeros_like(theta)
	S_alt = np.zeros_like(S)
	for i in range(len(Bvecs)):
		E_conv[i] = -get_vxb(vxyz[i],Bxyz[i])
		ma[i] = [ma_mean[i],1]

	angle = np.zeros_like(timeax)
	for i in range(len(timeax)): angle[i] = get_angle(Bvecs[i],Bvecs_mean[i])

	low = 25
	high = 180

	angle = np.zeros_like(timeax)
	for i in range(len(timeax)): angle[i] = get_angle(Bvecs[i],Bvecs_mean[i])
	angle_reduced = filter_angle(angle,low,high) #If you want all angles, do 0,180
	allangles = np.concat([allangles,angle_reduced])
	allmachs = np.concat([allmachs,ma_mean])
	allpositions = np.concat([allpositions,position])
	allSr = np.concat([allSr,S[:,0]])
	allKr = np.concat([allKr,K[:,0]])



#%%





plt.plot(allangles)
#%%
plt.plot(allpositions)
#%%
plt.plot(allSr)
#%%
plt.plot(allmachs)
#%%


cm = plt.cm.get_cmap('seismic')
sc=plt.scatter(allmachs,allSr,s=1)

plt.xlim(0,3)
plt.xlabel('Alfven Mach Number')
plt.ylabel('Sr')
plt.show()



#%%

# Scatter Plot w/ colorbar

#angle_reduced = bin_angle(angle)
cm = plt.cm.get_cmap('seismic')
sc=plt.scatter(allmachs,allangles,c=allSr,s=3,cmap=cm,vmin=-0.1,vmax=0.1)
# sc=plt.scatter(allpositions,allmachs,s=1)
#sc=plt.scatter(ma_mean,Br_brmean,c=S[:,0],s=3,cmap=cm,vmin=-0.3,vmax=0.3)

plt.colorbar(sc,label="Radial Poynting Flux")
plt.xlim(0,1)
plt.ylim(0,90)


plt.axhline(y=90,c='k')
plt.axhline(y=low,c='b',linestyle='dashed',linewidth=0.3)
plt.axvline(x=1,c='k')
plt.xlabel('Alfven Mach Number Ma')
plt.ylabel('Deflection Angle')
plt.show()

# %%
