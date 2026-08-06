#%%
import pyspedas
import numpy as np
import matplotlib.pyplot as plt
import pytplot
import matplotlib.cm as cm
from pytplot import get_data, store_data,timespan
from pyspedas import tplot
from pyspedas import tinterpol
from scipy.stats import binned_statistic_2d
from matplotlib.colors import ListedColormap
import cdflib
import matplotlib.colors as colors

me = 9.1094e-31 #kg
mi = 1837*me
mu0 = 1.2566370e-06  #;m kg / C^2
eps0 = 8.85e-12   # C^2/Nm^2
e = 1.602e-19 #C
Z = 1 # 1 for H+, 2 for He2+
gamma = 5/3
kb = 1.380649e-23


# Functions to get various useful parameters
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


def stability_condition(betapar, a, b, beta0):
    denom = (betapar - beta0)**b
    Tperppar = 1 + a/denom
    return Tperppar

def quickplot():
	pyspedas.tplot([B_name,vi_name,Ti_name,ni_name])
	return
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
def get_vecmean(vec,int):   # Vector mean
	vec1_mean = get_mean(vec[:,0],minutes*meaninterval)
	vec2_mean = get_mean(vec[:,1],minutes*meaninterval)
	vec3_mean = get_mean(vec[:,2],minutes*meaninterval)
	vec_mean = np.zeros_like(vec)
	for i in range(len(vec_mean)): vec_mean[i] = np.array([vec1_mean[i],vec2_mean[i],vec3_mean[i]])
	return vec_mean
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
def get_delta(vec,vec_mean):
	dvec = vec - vec_mean
	dvec_norm = dvec/np.linalg.norm(vec_mean)
	return dvec,dvec_norm
def get_deltascalar(var,var_mean):
	dvar = var - var_mean
	dvar_norm = dvar/var_mean
	return dvar,dvar_norm

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
def get_pm(B):
	pm = 0.5*(mu0**-1)*np.linalg.norm(B)**2
	return pm
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


def get_angle(vec,meanvec):
	term1 = np.dot(vec,meanvec)
	term2 = np.dot(np.linalg.norm(vec),np.linalg.norm(meanvec))
	term3 = term1/term2
	term4 = np.arccos(term3)*180/np.pi
	return term4


def get_vecmean(vec,int):   # Vector mean
	vec1_mean = get_mean(vec[:,0],minutes*meaninterval)
	vec2_mean = get_mean(vec[:,1],minutes*meaninterval)
	vec3_mean = get_mean(vec[:,2],minutes*meaninterval)
	vec_mean = np.zeros_like(vec)
	for i in range(len(vec_mean)): vec_mean[i] = np.array([vec1_mean[i],vec2_mean[i],vec3_mean[i]])
	return vec_mean


# Specific to WIND
def get_parperps(n,T,B):  #B ant T tensor coord systems need to match for this
	trace = T[0] + T[3] + T[5]
	term1 = (T[0]*B[0]**2 + T[3]*B[1]**2 + T[5]*B[2]**2)/(np.linalg.norm(B)**2)
	term2 = 2*(T[1]*B[0]*B[1] + T[2]*B[0]*B[2] + T[4]*B[1]*B[2])/(np.linalg.norm(B)**2)
	Tpar = term1+term2
	Tperp=0.5*(trace-Tpar)
	Ppar = n*Tpar
	Pperp = n*Tperp
	return Tpar,Tperp,Ppar,Pperp
def get_crosshelicity(v,B,n,m): #vector v & B
	z_plus = v + B/np.sqrt(n*m*mu0)
	z_minus = v - B/np.sqrt(n*m*mu0)
	term1 = np.linalg.norm(z_plus)**2 - np.linalg.norm(z_minus)**2
	term2 = np.linalg.norm(z_plus)**2 + np.linalg.norm(z_minus)**2
	sigma_c = term1/term2
	return sigma_c

#Somehow always = 1 (need to resolve)
def get_residenergy(v,B,n,m): # vector dv & dB (Alfven units??)
	term1 = np.linalg.norm(v)**2 - np.linalg.norm(B/np.sqrt(n*m*mu0))**2
	term2 = np.linalg.norm(v)**2 + np.linalg.norm(B/np.sqrt(n*m*mu0))**2
	sigma_r = term1/term2
	return sigma_r 











#%%


ev1 = ['2025-03-26/00:30', '2025-03-30/13:40'] #Full interval
ev1_reduced = ['2025-03-28/00:00', '2025-03-31/00:00'] #Full interval

# enc23coronalhole_fast = [['2025-03-26/14:50','2025-03-29/04:20']] # Encounter 21
# enc23coronalhole_slow = [['2025-03-26/02:20','2025-03-26/11:20']] # Encounter 21

enc23coronalhole_fast = [['2025-03-27/12:00','2025-03-28/00:00']] # Encounter 21
enc23coronalhole_slow = [['2025-03-30/12:00','2025-03-31/12:00']] # Encounter 21
enc23coronalhole_bound = [['2025-03-26/06:00','2025-03-27/12:00']] # Encounter 21
enc23coronalhole_full = [['2025-03-26/00:00','2025-03-27/12:00']] # Encounter 21


enc23coronalhole_fast = [['2025-03-27/12:00','2025-03-28/00:00']] # Encounter 21
enc23coronalhole_shoulder = [['2025-03-26/14:00','2025-03-27/00:00']] # Encounter 21
enc23coronalhole_trailing = [['2025-03-30/12:00','2025-03-31/12:00']] # Encounter 21
enc23coronalhole_preceding = [['2025-03-25/22:00','2025-03-26/13:00']] # Encounter 21

eventlist=enc23coronalhole_shoulder


allB_list = []
allv_list = []
alln_list = []
allT_list = []
allpositions_list = []
allbeta_par_list = []
allangles_list = []
allTpar_list = []
allTperp_list = []
allTparperp_list = []
allangles_list = []
alldB_list = []
alldv_list = []
alldn_list = []
alldT_list = []
alldvsqrdmag_list = []
alldBsqrdmag_list = []

# allB = np.array([])
# allv = np.array([])
# alln = np.array([])
# allT = np.array([])
# allpositions = np.array([])
# allbeta_par = np.array([])
# allangles = np.array([])
# allTpar = np.array([])
# allTperp = np.array([])
# allTparperp = np.array([])

for i in range(len(eventlist)):
	trange = eventlist[i]

	ion_vars = pyspedas.projects.wind.swe(trange=trange,datatype='h1',time_clip=True)
	electron_vars = pyspedas.projects.wind.swe(trange=trange,datatype='h5',time_clip=True)
	# sms_vars = pyspedas.projects.wind.sms(trange=trange, datatype='k0')
	mfi_vars = pyspedas.projects.wind.mfi(trange=trange,time_clip=True)
	B_name='BGSE'
	# ni_name='N_elec'
	ni_name='Proton_Np_moment'
	# vi_name = 'U_eGSE'
	vi_name = 'Proton_V_moment'
	pos_name = 'DIST' # Distance from earth, in units of earth radii.  1 Re ~ 4.3e-5 AU
	Ptensor_name = 'P_eGSE'
	Tscalar_name = 'T_elec'
	W_name = 'Proton_W_moment' 
	Wpar_name = 'Proton_Wpar_moment' 
	Wperp_name = 'Proton_Wperp_moment' 

	interpvar_name = ni_name
	timeax = pytplot.get_data(interpvar_name).times
	meaninterval = mean_int(timeax,trange)
	tinterpol(B_name, interpvar_name, newname='B')
	tinterpol(Ptensor_name, interpvar_name, newname='PTensor')
	tinterpol(ni_name, interpvar_name, newname='ni')
	tinterpol(vi_name, interpvar_name, newname='vi')
	tinterpol(W_name, interpvar_name, newname='W')
	tinterpol(Wpar_name, interpvar_name, newname='Wpar')
	tinterpol(Wperp_name, interpvar_name, newname='Wperp')

	Bvecs = 1e-9 * reform(pytplot.get_data('B'))
	PiTensor = 1e-1 * reform(pytplot.get_data('PTensor'))
	ni = 1e6 * reform(pytplot.get_data('ni'))
	vi = 1e3 * reform(pytplot.get_data('vi'))
	# T = 1.602e-19*(8.617e-5)*reform(pytplot.get_data('T'))
	W = 1e3*reform(pytplot.get_data('W'))
	Wpar = 1e3*reform(pytplot.get_data('Wpar'))
	Wperp = 1e3*reform(pytplot.get_data('Wperp'))
	T = (mi*W**2)/(2*kb)
	Tpar = (mi*Wpar**2)/(2*kb)
	Tperp = (mi*Wperp**2)/(2*kb)

	# Tpar = np.zeros_like(Bvecs[:, 0])
	# Tperp = np.zeros_like(Tpar)
	Pperp = np.zeros_like(Tpar)
	Ppar = np.zeros_like(Tpar)
	P_mag = np.zeros_like(Tpar)
	beta_par = np.zeros_like(Tpar)
	angle = np.zeros_like(Tpar)
	dB,dB_norm = np.zeros_like(Bvecs), np.zeros_like(Bvecs)
	dv,dv_norm = np.zeros_like(vi), np.zeros_like(vi)
	dvsqrdmag,dvsqrdmag_norm = np.zeros_like(Tpar), np.zeros_like(Tpar)
	dBsqrdmag,dBsqrdmag_norm = np.zeros_like(Tpar), np.zeros_like(Tpar)
	dn,dn_norm = np.zeros_like(ni), np.zeros_like(ni)
	dT,dT_norm = np.zeros_like(T), np.zeros_like(T)



	minutes = 60
	# Means
	Bvecs_mean = get_vecmean(Bvecs,minutes*meaninterval)
	# vivecs_mean = get_vecmean(vi,minutes*meaninterval)
	n_mean = get_mean(ni,minutes*meaninterval)
	v_mean = get_mean(vi,minutes*meaninterval)
	T_mean = get_mean(T,minutes*meaninterval)

	for j in range(len(Bvecs)):  # Fixed: was shadowing outer loop variable i
		# Tpar[j], Tperp[j], Ppar[j], Pperp[j] = get_parperps(ni[j], PiTensor[j] / ni[j], Bvecs[j])

		P_mag[j] = get_pm(Bvecs[j])
		# beta_par[j] = Ppar[j]/P_mag[j]
		beta_par[j] = (3/2)*ni[j]*kb*Tpar[j]/P_mag[j]
		angle[j] = get_angle(Bvecs[j], Bvecs_mean[j])
		dB[j],dB_norm[j] = get_delta(Bvecs[j], Bvecs_mean[j])
		# dv[j],dv_norm[j] = get_delta(vi[j], vivecs_mean[j])
		dn[j],dn_norm[j] = get_deltascalar(ni[j], n_mean[j])
		dT[j],dT_norm[j] = get_deltascalar(T[j], T_mean[j])
		dvsqrdmag[j],dvsqrdmag_norm[j] = get_deltascalar(np.linalg.norm(vi[j])**2, np.linalg.norm(vivecs_mean[j])**2)
		dBsqrdmag[j],dBsqrdmag_norm[j] = get_deltascalar(np.linalg.norm(Bvecs[j])**2, np.linalg.norm(Bvecs_mean[j])**2)

	
	allB_list.append(Bvecs)
	alln_list.append(ni)
	allT_list.append(T)
	allv_list.append(vi)
	allbeta_par_list.append(beta_par)
	allangles_list.append(angle)	
	allTpar_list.append(Tpar)
	allTperp_list.append(Tperp)
	allTparperp_list.append(Tpar / Tperp)
	allangles_list.append(angle)
	alldB_list.append(dB)
	# alldv_list.append(dv)
	alldn_list.append(dn)
	alldT_list.append(dT)
	alldvsqrdmag_list.append(dvsqrdmag)
	alldBsqrdmag_list.append(dBsqrdmag)

allB = np.concatenate(allB_list, axis=0)
allv = np.concatenate(allv_list, axis=0)
alln = np.concatenate(alln_list, axis=0)
allT = np.concatenate(allT_list, axis=0)
# allpositions = np.concatenate(allpositions_list, axis=0)
allbeta_par = np.concatenate(allbeta_par_list, axis=0)
# allangles = np.concatenate(allangles_list, axis=0)
allTpar = np.concatenate(allTpar_list, axis=0)
allTperp = np.concatenate(allTperp_list, axis=0)
allTparperp = np.concatenate(allTparperp_list, axis=0)
allangles = np.concatenate(allangles_list, axis=0)
alldB = np.concatenate(alldB_list, axis=0)
# alldv = np.concatenate(alldv_list, axis=0)
alldn = np.concatenate(alldn_list, axis=0)
alldT = np.concatenate(alldT_list, axis=0)
alldvsqrdmag = np.concatenate(alldvsqrdmag_list, axis=0)
alldBsqrdmag = np.concatenate(alldBsqrdmag_list, axis=0)

allvmags = allv
allBmags = np.zeros_like(allB[:,0])
for i in range(len(allv)):
	# allvmags[i] = np.linalg.norm(allv[i])
	allBmags[i] = np.linalg.norm(allB[i])


def get_windvals():
	n = 1e-6*np.nanmean(alln)	
	nstd = 1e-6*np.nanstd(alln)
	T = np.nanmean(allT)/11604
	Tstd = np.nanstd(allT)/11604
	v = 1e-3*np.nanmean(allv)
	vstd = 1e-3*np.nanstd(allv)
	B = 1e9*np.nanmean(allBmags)
	Bstd = 1e9*np.nanstd(allBmags)
	dB = 1e9*np.sqrt(np.nanmean(alldBsqrdmag))
	dBstd = 1e9*np.sqrt(np.nanstd(alldBsqrdmag))

	return n,nstd,T,Tstd,v,vstd,B,Bstd,dB,dBstd

def get_wind_anisos():
	betapar = allbeta_par
	Tparperp = allTparperp
	return betapar, Tparperp

# if (eventlist == enc23coronalhole_fast):
# 	nwind_fast,nstdwind_fast,Twind_fast,Tstdwind_fast,vwind_fast,vstdwind_fast,Bwind_fast,Bstdwind_fast,dBwind_fast,dBstdwind_fast = get_windvals()
# 	betaparswind_fast, Tparperpwind_fast = get_wind_anisos()
# 	print("Fast Variables")
# elif (eventlist == enc23coronalhole_slow):
# 	nwind_slow,nstdwind_slow,Twind_slow,Tstdwind_slow,vwind_slow,vstdwind_slow,Bwind_slow,Bstdwind_slow,dBwind_slow,dBstdwind_slow = get_windvals()
# 	betaparswind_slow, Tparperpwind_slow = get_wind_anisos()
# 	print("Slow Variables")
# else:
# 	nwind,nstdwind,Twind,Tstdwind,vwind,vstdwind,Bwind,Bstdwind,dBwind,dBstdwind = get_windvals()
# 	betaparswind, Tparperpwind = get_wind_anisos()
# 	print("Variables")
if (eventlist == enc23coronalhole_fast):
	nwind_fast,nstdwind_fast,Twind_fast,Tstdwind_fast,vwind_fast,vstdwind_fast,Bwind_fast,Bstdwind_fast,dBwind_fast,dBstdwind_fast = get_windvals()
	betaparswind_fast, Tparperpwind_fast = get_wind_anisos()
	print("Fast Variables")
elif (eventlist == enc23coronalhole_preceding):
	nwind_preceding,nstdwind_preceding,Twind_preceding,Tstdwind_preceding,vwind_preceding,vstdwind_preceding,Bwind_preceding,Bstdwind_preceding,dBwind_preceding,dBstdwind_preceding = get_windvals()
	betaparswind_preceding, Tparperpwind_preceding = get_wind_anisos()
	print("Preceding Variables")
elif (eventlist == enc23coronalhole_trailing):
	nwind_trailing,nstdwind_trailing,Twind_trailing,Tstdwind_trailing,vwind_trailing,vstdwind_trailing,Bwind_trailing,Bstdwind_trailing,dBwind_trailing,dBstdwind_trailing = get_windvals()
	betaparswind_trailing, Tparperpwind_trailing = get_wind_anisos()
	print("Trailing Variables")
elif (eventlist == enc23coronalhole_shoulder):
	nwind_shoulder,nstdwind_shoulder,Twind_shoulder,Tstdwind_shoulder,vwind_shoulder,vstdwind_shoulder,Bwind_shoulder,Bstdwind_shoulder,dBwind_shoulder,dBstdwind_shoulder = get_windvals()
	betaparswind_shoulder, Tparperpwind_shoulder = get_wind_anisos()
	print("Shoulder Variables")
else:
	nwind,nstdwind,Twind,Tstdwind,vwind,vstdwind,Bwind,Bstdwind,dBwind,dBstdwind = get_windvals()
	betaparswind, Tparperpwind = get_wind_anisos()
	print("Variables")
# %%






plt.hist(allangles,bins=80)
plt.xlim(0,180)
plt.yscale('log')
plt.axvline(np.nanmean(allangles),color='k',label = r'$\langle\theta\rangle$')
plt.axvline(np.nanquantile(allangles,0.90),color='r',linestyle='dashed',label = r'$90^{th}$ Percentile')
plt.xlabel(r'$\theta$', fontsize=15)
plt.legend()

print(np.nanmean(allangles))
print(np.nanstd(allangles))
print(np.nanquantile(allangles,0.90))
# mask = allpositions >0
# brazilplot(1e3*allbeta_par, r'$\beta_\parallel$', r'$T_{\perp}/T_{\parallel}$', r'$v_x$')
# %%

# mask = 1e-3*allvmags >= 400

#%%

# Conjunction Stats

# # Normalized values for conjunctions
# n_psp = 48.7
# n_err_psp = 36.8
# n_solo = 25.4
# n_err_solo = 10.4
# n_wind = 23.16
# n_err_wind = 17.7


# T_psp = 5.47
# T_err_psp = 3
# T_solo = 24.2
# T_err_solo = 9.7
# T_wind = 45.2
# T_err_wind = 13.3

# v_psp = 353.8
# v_err_psp = 56.2
# v_solo = 602.3
# v_err_solo = 118.5
# v_wind = 459.7
# v_err_wind = 105.3	

# nvals = np.array([n_psp,n_solo,n_wind])
# n_errs = np.array([n_err_psp,n_err_solo,n_err_wind])
# Tvals = np.array([T_psp,T_solo,T_wind])
# T_errs = np.array([T_err_psp,T_err_solo,T_err_wind])
# vvals = np.array([v_psp,v_solo,v_wind])
# v_errs = np.array([v_err_psp,v_err_solo,v_err_wind])

# x = np.array([1,2,3])



# plt.errorbar(x,nvals,yerr=n_errs,fmt='o',label=r'$n$',capsize=5)
# plt.xticks(x,['PSP','SOLO','WIND'])
# plt.xlabel('Spacecraft')
# plt.legend()
# plt.errorbar(x,Tvals,yerr=T_errs,fmt='o',label=r'$T$',capsize=5)
# plt.xticks(x,['PSP','SOLO','WIND'])
# plt.xlabel('Spacecraft')
# plt.legend()

# plt.errorbar(x,vvals,yerr=v_errs,fmt='o',label=r'$v$',capsize=5)
# plt.xticks(x,['PSP','SOLO','WIND'])
# plt.xlabel('Spacecraft')
# plt.legend()
# plt.xticks(x,['PSP','SOLO','WIND'])
# plt.xlabel('Spacecraft')
# plt.legend()

# plt.errorbar(x,nvals/nvals[0],yerr=n_errs/nvals[0],fmt='o-',label=r'$R^2 n/n_{psp}$',capsize=5)
# plt.xticks(x,['PSP','SOLO','WIND'])
# plt.xlabel('Spacecraft')
# plt.legend()
# plt.errorbar(x,Tvals/Tvals[0],yerr=T_errs/Tvals[0],fmt='o-',label=r'$R^{4/3} T/T_{psp}$',capsize=5)
# plt.xticks(x,['PSP','SOLO','WIND'])
# plt.xlabel('Spacecraft')
# plt.legend()

# plt.errorbar(x,vvals**2/(vvals[0]**2),yerr=v_errs**2/vvals[0]**2,fmt='o-',label=r'$v^2/v^2_{psp}$',capsize=5)
# plt.xticks(x,['PSP','SOLO','WIND'])
# plt.xlabel('Spacecraft')
# plt.legend()
# plt.xticks(x,['PSP','SOLO','WIND'])
# plt.xlabel('Spacecraft')
# plt.axhline(y=1,linestyle='dashed',color='k')
# plt.ylim(0,11)
# # plt.legend()



# %%

#Need to update with WIND ion results

#And double-check those dB terms

#New Conjunction Stats - Fast Wind
x = np.array([1,2,3])

def get_fastvars():
	narray = 1e6*np.array([3051,48,5])
	n_stdarray = 1e6*np.array([451,11,4])
	Tarray = 1.602e-19*np.array([148,51,16])
	T_stdarray = 1.602e-19*np.array([61,13,5])
	varray = 1e3*np.array([452,683,479])
	v_stdarray = 1e3*np.array([47, 83, 107])
	Barray = 1e-9*np.array([1966,56,9])
	B_stdarray = 1e-9*np.array([288, 4, 5])
	dB_array = 1e-9*np.array([367, 26, 5])
	dB_stdarray = 1e-9*np.array([253, 5, 4])
	rarray = np.array([0.05/3,0.3/3,1/3])
	# rarray = np.array([1,1,1])  # For now, just use the actual values for each spacecraft.  Can normalize later if needed.

	narray_norm = (rarray**2)*narray
	n_stdarray_norm = (rarray**2)*n_stdarray
	Tarray_norm = (rarray**(4/3))*Tarray
	T_stdarray_norm = (rarray**(4/3))*T_stdarray
	varray_norm = varray
	v_stdarray_norm = v_stdarray
	Barray_norm = (rarray**2)*Barray
	B_stdarray_norm = (rarray**2)*B_stdarray
	dB_array_norm = (rarray**2)*dB_array
	dB_stdarray_norm = (rarray**2)*dB_stdarray

	var_n = narray_norm/narray_norm[0]
	std_n = n_stdarray_norm/narray_norm[0]
	var_T = Tarray_norm/Tarray_norm[0]
	std_T = T_stdarray_norm/Tarray_norm[0]
	var_v = varray_norm/varray_norm[0]
	std_v = v_stdarray_norm/varray_norm[0]
	var_B = Barray_norm/Barray_norm[0]
	std_B = B_stdarray_norm/Barray_norm[0]
	var_dB = dB_array_norm/dB_array_norm[0]
	std_dB = dB_stdarray_norm/dB_array_norm[0]

	x = np.array([1,2,3])
	return narray_norm,Tarray_norm,varray_norm,Barray_norm,dB_array_norm, n_stdarray_norm, T_stdarray_norm, v_stdarray_norm, B_stdarray_norm, dB_stdarray_norm


def get_slowvars():
	# # narray = 1e6*np.array([4346,216,10])
	# narray = 1e6*np.array([4346,216,25])
	# # n_stdarray = 1e6*np.array([781,48,5])
	# n_stdarray = 1e6*np.array([781,48,8])
	# # Tarray = 1.602e-19*np.array([69,8,18])
	# Tarray = 1.602e-19*np.array([69,8,15])
	# # T_stdarray = 1.602e-19*np.array([20,3,3])
	# T_stdarray = 1.602e-19*np.array([20,3,3])
	# # varray = 1e3*np.array([228,253,421])
	# varray = 1e3*np.array([228,253,371])
	# # v_stdarray = 1e3*np.array([20,30,85])
	# v_stdarray = 1e3*np.array([20,30,29])
	# Barray = 1e-9*np.array([1172,35,16])
	# B_stdarray = 1e-9*np.array([187, 4, 5])
	# dB_array = 1e-9*np.array([263, 17, 10])
	# dB_stdarray = 1e-9*np.array([162, 9, 5])
	# rarray = np.array([0.05/3,0.3/3,1/3])
	# # rarray = np.array([1,1,1])  # For now, just use the actual values for each spacecraft.  Can normalize later if needed.





	narray_norm = (rarray**2)*narray
	n_stdarray_norm = (rarray**2)*n_stdarray
	Tarray_norm = (rarray**(4/3))*Tarray
	T_stdarray_norm = (rarray**(4/3))*T_stdarray
	varray_norm = varray
	v_stdarray_norm = v_stdarray
	Barray_norm = (rarray**2)*Barray
	B_stdarray_norm = (rarray**2)*B_stdarray
	dB_array_norm = (rarray**2)*dB_array
	dB_stdarray_norm = (rarray**2)*dB_stdarray

	# var_n = narray_norm/narray_norm[0]
	# std_n = n_stdarray_norm/narray_norm[0]
	# var_T = Tarray_norm/Tarray_norm[0]
	# std_T = T_stdarray_norm/Tarray_norm[0]
	# var_v = varray_norm/varray_norm[0]
	# std_v = v_stdarray_norm/varray_norm[0]
	# var_B = Barray_norm/Barray_norm[0]
	# std_B = B_stdarray_norm/Barray_norm[0]
	# var_dB = dB_array_norm/dB_array_norm[0]
	# std_dB = dB_stdarray_norm/dB_array_norm[0]

	return narray_norm,Tarray_norm,varray_norm,Barray_norm,dB_array_norm, n_stdarray_norm, T_stdarray_norm, v_stdarray_norm, B_stdarray_norm, dB_stdarray_norm


nfast,Tfast,vfast,Bfast,dBfast,nfast_std,Tfast_std,vfast_std,Bfast_std,dBfast_std = get_fastvars()
nslow,Tslow,vslow,Bslow,dBslow,nslow_std,Tslow_std,vslow_std,Bslow_std,dBslow_std = get_slowvars()


n = nfast + nslow
T = (Tfast + Tslow)/2
bulkfast = 0.5*nfast*mi*vfast**2
bulkslow = 0.5*nslow*mi*vslow**2
thermfast = (3/2)*kb*nfast*Tfast
thermslow = (3/2)*kb*nslow*Tslow
magfast = (0.5*Bfast**2)/(mu0)
magslow = (0.5*Bslow**2)/(mu0)
magmean = (magfast + magslow)/2
deltamagfast = (0.5*dBfast**2)/(mu0)
deltamagslow = (0.5*dBslow**2)/(mu0)
plasmafast = thermfast + bulkfast
plasmaslow = thermslow + bulkslow

plasma = plasmafast + plasmaslow
fieldfast = magfast + deltamagfast
fieldslow = magslow + deltamagslow
field = fieldfast + fieldslow
#%%
# Number Density
plt.errorbar(x,nfast,yerr=nfast_std,fmt='o-',color = 'r',label = r'$n^{\prime}_{fast}$',capsize=5)
plt.errorbar(x,nslow,yerr=nslow_std,fmt='o-',color='b',label = r'$n^{\prime}_{slow}$',capsize=5)
plt.plot(x,(nslow+nfast)/2,color='g',marker='s',linestyle='dashed',label = r'$n^{\prime}_{avg}$')
plt.title('Adjusted Number Density',fontsize=15)
plt.axvline(x=1,linestyle='dashed',color='k',linewidth=0.2)
plt.axvline(x=2,linestyle='dashed',color='k',linewidth=0.2)
plt.axvline(x=3,linestyle='dashed',color='k',linewidth=0.2)
plt.xticks(x,['PSP \n $\sim0.05 \  AU$','SOLO \n $\sim0.3 \  AU$','WIND \n $\sim1 \ AU$	'],fontsize=12)
plt.legend(fontsize=15)

#%%
# Temperature
plt.errorbar(x,Tfast,yerr=Tfast_std,fmt='o-',color = 'r',label = r'$T^{\prime}_{fast}$',capsize=5)
plt.errorbar(x,Tslow,yerr=Tslow_std,fmt='o-',color='b',label = r'$T^{\prime}_{slow}$',capsize=5)
plt.plot(x,(Tslow+Tfast)/2,color='g',marker='s',linestyle='dashed',label = r'$T^{\prime}_{avg}$')
plt.title('Adjusted Temperature',fontsize=15)
plt.axvline(x=1,linestyle='dashed',color='k',linewidth=0.2)
plt.axvline(x=2,linestyle='dashed',color='k',linewidth=0.2)
plt.axvline(x=3,linestyle='dashed',color='k',linewidth=0.2)
plt.xticks(x,['PSP \n $\sim0.05 \  AU$','SOLO \n $\sim0.3 \  AU$','WIND \n $\sim1 \ AU$	'],fontsize=12)
plt.legend(fontsize=15)



#%%
# Velocity
plt.errorbar(x,1e-3*vfast,yerr=1e-3*vfast_std,fmt='o-',color = 'r',label = r'$v_{fast}$',capsize=5)
plt.errorbar(x,1e-3*vslow,yerr=1e-3*vslow_std,fmt='o-',color='b',label = r'$v_{slow}$',capsize=5)
plt.errorbar(x,1e-3*(vslow+vfast)/2,color = 'g',marker='s',linestyle='dashed',label = r'$v_{avg}$')

# plt.ylim(200,800)
plt.title('Velocity [km/s]',fontsize=15)
plt.axvline(x=1,linestyle='dashed',color='k',linewidth=0.2)
plt.axvline(x=2,linestyle='dashed',color='k',linewidth=0.2)
plt.axvline(x=3,linestyle='dashed',color='k',linewidth=0.2)
plt.ylim(200,800)
plt.yticks(fontsize=12)
plt.xticks(x,['PSP \n $\sim0.05 \  AU$','SOLO \n $\sim0.3 \  AU$','WIND \n $\sim1 \ AU$	'],fontsize=12)
plt.legend(fontsize=15)

#%%
# Magnetic Field Strength
plt.errorbar(x,Bfast,yerr=Bfast_std,fmt='o-',color = 'r',label = r'$|B|_{fast}$',capsize=5)
plt.errorbar(x,Bslow,yerr=Bslow_std,fmt='o-',color='b',label = r'$|B|_{slow}$',capsize=5)
plt.errorbar(x,(Bslow+Bfast)/2,color = 'g',marker='s',linestyle='dashed',label = r'$|B|_{avg}$')
plt.title('Adjusted Magnetic Field Strength [T]')
plt.ylim(0,2.5e-9)
plt.axvline(x=1,linestyle='dashed',color='k',linewidth=0.2)
plt.axvline(x=2,linestyle='dashed',color='k',linewidth=0.2)
plt.axvline(x=3,linestyle='dashed',color='k',linewidth=0.2)
plt.yticks(fontsize=12)
plt.xticks(x,['PSP \n $\sim0.05 \  AU$','SOLO \n $\sim0.3 \  AU$','WIND \n $\sim1 \ AU$	'],fontsize=12)
plt.legend(fontsize=15)

#%%
# Magnetic Fluctuation Strength
plt.errorbar(x,dBfast,yerr=dBfast_std,fmt='o-',color = 'r',label = r'$\delta B_{fast}$',capsize=5)
plt.errorbar(x,dBslow,yerr=dBslow_std,fmt='o-',color='b',label = r'$\delta B_{slow}$',capsize=5)
plt.errorbar(x,(dBslow+dBfast)/2,color = 'g',marker='s',linestyle='dashed',label = r'$\delta B_{avg}$')
plt.ylim(0,2.5e-9)
plt.title('Adjusted Magnetic Fluctuation Magnitude [T]')
plt.axvline(x=1,linestyle='dashed',color='k',linewidth=0.2)
plt.axvline(x=2,linestyle='dashed',color='k',linewidth=0.2)
plt.axvline(x=3,linestyle='dashed',color='k',linewidth=0.2)
plt.yticks(fontsize=12)
plt.xticks(x,['PSP \n $\sim0.05 \  AU$','SOLO \n $\sim0.3 \  AU$','WIND \n $\sim1 \ AU$	'],fontsize=12)
plt.legend(fontsize=15)



#%%
# Magnetic Fluctuation Ratio
fastratio_err = (dBfast/Bfast)*np.sqrt((dBfast_std/dBfast)**2 + (Bfast_std/Bfast)**2)
slowsratio_err = (dBslow/Bslow)*np.sqrt((dBslow_std/dBslow)**2 + (Bslow_std/Bslow)**2)

plt.errorbar(x,dBfast/Bfast,yerr=fastratio_err,fmt='o-',color = 'r',label = r'$|\delta B/B|_{fast}$',capsize=5)
plt.errorbar(x,dBslow/Bslow,yerr=slowsratio_err,fmt='o-',color='b',label = r'$|\delta B/B|_{slow}$',capsize=5)
plt.plot(x,(dBslow/Bslow+dBfast/Bfast)/2,color='g',marker='s',linestyle='dashed',label = r'$|\delta B/B|_{avg}$')
plt.title('$\delta B/B$',fontsize=15)
plt.axvline(x=1,linestyle='dashed',color='k',linewidth=0.2)
plt.axvline(x=2,linestyle='dashed',color='k',linewidth=0.2)
plt.axvline(x=3,linestyle='dashed',color='k',linewidth=0.2)
plt.yticks(fontsize=12)
plt.xticks(x,['PSP \n $\sim0.05 \  AU$','SOLO \n $\sim0.3 \  AU$','WIND \n $\sim1 \ AU$	'],fontsize=12)
plt.legend(fontsize=15)



#%%
# Bulk Energy Density
bulkfast_err = bulkfast*np.sqrt((nfast_std/nfast)**2 + 2*(vfast_std/vfast)**2)
bulkslow_err = bulkslow*np.sqrt((nslow_std/nslow)**2 + 2*(vslow_std/vslow)**2)

plt.errorbar(x,bulkfast,yerr=bulkfast_err,fmt='o-',color = 'r',label = r'Fast Wind',capsize=5)
plt.errorbar(x,bulkslow,yerr=bulkslow_err,fmt='o-',color='b',label = r'Slow Wind',capsize=5)
plt.plot(x,(bulkfast+bulkslow)/2,color='g',marker='s',linestyle='dashed',label = r'Mean')
plt.title('Adjusted Bulk Energy Density [J/m^3]',fontsize=15)
plt.axvline(x=1,linestyle='dashed',color='k',linewidth=0.2)
plt.axvline(x=2,linestyle='dashed',color='k',linewidth=0.2)
plt.axvline(x=3,linestyle='dashed',color='k',linewidth=0.2)
plt.xticks(x,['PSP','SOLO','WIND'])
plt.legend(fontsize=15)

#%%
# Thermal Energy Density
thermfast_err = thermfast*np.sqrt((nfast_std/nfast)**2 + (Tfast_std/Tfast)**2)
thermslow_err = thermslow*np.sqrt((nslow_std/nslow)**2 + (Tslow_std/Tslow)**2)

plt.errorbar(x,thermfast,yerr=thermfast_err,fmt='o-',color = 'r',label = r'Fast Wind',capsize=5)
plt.errorbar(x,thermslow,yerr=thermslow_err,fmt='o-',color='b',label = r'Slow Wind',capsize=5)
plt.plot(x,(thermfast+thermslow)/2,color='g',marker='s',linestyle='dashed',label = r'Mean')
plt.title('Adjusted Thermal Energy Density',fontsize=15)
plt.axvline(x=1,linestyle='dashed',color='k',linewidth=0.2)
plt.axvline(x=2,linestyle='dashed',color='k',linewidth=0.2)
plt.axvline(x=3,linestyle='dashed',color='k',linewidth=0.2)
plt.xticks(x,['PSP','SOLO','WIND'])
plt.legend(fontsize=15)

# %%
# Magnetic Energy Density
magfast_err = magfast*np.sqrt(2*(Bfast_std/Bfast)**2)
magslow_err = magslow*np.sqrt(2*(Bslow_std/Bslow)**2)
plt.errorbar(x,magfast,yerr=magfast_err,fmt='o-',color = 'r',label = r'Fast Wind',capsize=5)
plt.errorbar(x,magslow,yerr=magslow_err,fmt='o-',color='b',label = r'Slow Wind',capsize=5)
plt.plot(x,(magfast+magslow)/2,color='g',marker='s',linestyle='dashed',label = r'Mean')
plt.title('Adjusted Magnetic Energy Density')
plt.ylim(0,2e-12)
plt.axvline(x=1,linestyle='dashed',color='k',linewidth=0.2)
plt.axvline(x=2,linestyle='dashed',color='k',linewidth=0.2)
plt.axvline(x=3,linestyle='dashed',color='k',linewidth=0.2)
plt.xticks(x,['PSP','SOLO','WIND'])
plt.xlabel('Spacecraft')
plt.legend(fontsize=15)

#%%
#Magnetic Fluctuation Energy Density
deltamagfast_err = deltamagfast*np.sqrt(2*(dBfast_std/dBfast)**2)
deltamagslow_err = deltamagslow*np.sqrt(2*(dBslow_std/dBslow)**2)
plt.errorbar(x,deltamagfast,yerr=deltamagfast_err,fmt='o-',color = 'r',label = r'Fast Wind',capsize=5)
plt.errorbar(x,deltamagslow,yerr=deltamagslow_err,fmt='o-',color='b',label = r'Slow Wind',capsize=5)
plt.plot(x,(deltamagfast+deltamagslow)/2,color='g',marker='s',linestyle='dashed',label = r'Mean')
plt.title('Magnetic Fluctuation Energy Density [J/m^3]')
plt.axvline(x=1,linestyle='dashed',color='k',linewidth=0.2)
plt.axvline(x=2,linestyle='dashed',color='k',linewidth=0.2)
plt.axvline(x=3,linestyle='dashed',color='k',linewidth=0.2)
plt.xticks(x,['PSP','SOLO','WIND'])
plt.xlabel('Spacecraft')
plt.legend(fontsize=15)
#%%
#%%


# %%
#%%

plt.plot(x,thermfast/nfast,color='r',marker='o',label = r'Fast')
plt.plot(x,thermslow/nslow,color = 'b',marker = 'o',label = r'Slow')
plt.title('Thermal Energy per Particle [J/part.]')
plt.axvline(x=1,linestyle='dashed',color='k',linewidth=0.2)
plt.axvline(x=2,linestyle='dashed',color='k',linewidth=0.2)
plt.axvline(x=3,linestyle='dashed',color='k',linewidth=0.2)
plt.xticks(x,['PSP','SOLO','WIND'])
plt.legend(fontsize=15)
# %%
plt.plot(x,deltamagfast/nfast,color='r',marker='o',label = r'Fast Wind')
plt.plot(x,deltamagslow/nslow,color = 'b',marker = 'o',label = r'Slow Wind')
plt.title('Magnetic Fluctuation Energy per Particle [J/part.]')
plt.axvline(x=1,linestyle='dashed',color='k',linewidth=0.2)
plt.axvline(x=2,linestyle='dashed',color='k',linewidth=0.2)
plt.axvline(x=3,linestyle='dashed',color='k',linewidth=0.2)
plt.xticks(x,['PSP','SOLO','WIND'])
plt.legend(fontsize=15)

#%%

plt.plot(x,magfast/nfast,color='r',marker='o',label = r'Fast Wind')
plt.plot(x,magslow/nslow,color = 'b',marker = 'o',label = r'Slow Wind')
plt.title('Magnetic Energy per Particle [J/part.]')
plt.axvline(x=1,linestyle='dashed',color='k',linewidth=0.2)
plt.axvline(x=2,linestyle='dashed',color='k',linewidth=0.2)
plt.axvline(x=3,linestyle='dashed',color='k',linewidth=0.2)
plt.xticks(x,['PSP','SOLO','WIND'])
plt.legend(fontsize=15)

#%%


#%%