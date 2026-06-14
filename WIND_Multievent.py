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

eventlist=[ev1]
#%%



# %%
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

	swe_vars = pyspedas.projects.wind.swe(trange=trange)
	mfi_vars = pyspedas.projects.wind.mfi(trange=trange)
	B_name='BGSE'
	ni_name='N_elec'
	vi_name = 'U_eGSE'
	pos_name = 'DIST' # Distance from earth, in units of earth radii.  1 Re ~ 4.3e-5 AU
	Ptensor_name = 'P_eGSE'
	Tscalar_name = 'T_elec' 

	interpvar_name = ni_name
	timeax = pytplot.get_data(interpvar_name).times
	meaninterval = mean_int(timeax,trange)
	tinterpol(B_name, interpvar_name, newname='B')
	tinterpol(Ptensor_name, interpvar_name, newname='PTensor')
	tinterpol(ni_name, interpvar_name, newname='ni')
	tinterpol(vi_name, interpvar_name, newname='vi')
	tinterpol(Tscalar_name, interpvar_name, newname='T')

	Bvecs = 1e-9 * reform(pytplot.get_data('B'))
	PiTensor = 1e-1 * reform(pytplot.get_data('PTensor'))
	ni = 1e6 * reform(pytplot.get_data('ni'))
	vi = 1e-3 * reform(pytplot.get_data('vi'))
	T = 1.602e-19*(8.617e-5)*reform(pytplot.get_data('T'))

	Tpar = np.zeros_like(Bvecs[:, 0])
	Tperp = np.zeros_like(Tpar)
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
	vivecs_mean = get_vecmean(vi,minutes*meaninterval)
	n_mean = get_mean(ni,minutes*meaninterval)
	T_mean = get_mean(T,minutes*meaninterval)

	for j in range(len(Bvecs)):  # Fixed: was shadowing outer loop variable i
		Tpar[j], Tperp[j], Ppar[j], Pperp[j] = get_parperps(ni[j], PiTensor[j] / ni[j], Bvecs[j])

		P_mag[j] = get_pm(Bvecs[j])
		beta_par[j] = Ppar[j]/P_mag[j]
		angle[j] = get_angle(Bvecs[j], Bvecs_mean[j])
		dB[j],dB_norm[j] = get_delta(Bvecs[j], Bvecs_mean[j])
		dv[j],dv_norm[j] = get_delta(vi[j], vivecs_mean[j])
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
	alldv_list.append(dv)
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
alldv = np.concatenate(alldv_list, axis=0)
alldn = np.concatenate(alldn_list, axis=0)
alldT = np.concatenate(alldT_list, axis=0)
alldvsqrdmag = np.concatenate(alldvsqrdmag_list, axis=0)
alldBsqrdmag = np.concatenate(alldBsqrdmag_list, axis=0)

allvmags = np.zeros_like(allv[:,0])
allBmags = np.zeros_like(allB[:,0])
for i in range(len(allv)):
	allvmags[i] = np.linalg.norm(allv[i])
	allBmags[i] = np.linalg.norm(allB[i])



#%%
# Normalized arrays
Tvals = np.array([9,38,72])
nvals = np.array([49.6,25.1,23.4])
vvals = np.array([354,602,459])
Pvals = nvals*Tvals
nvsqrdvals = nvals*(vvals**2)
plt.plot(Pvals/Pvals[0],marker='o',linestyle='-',label = r'$P/P_0$')
plt.plot(Tvals/Tvals[0],marker='o',linestyle='-',label = r'$T/T_0$')
plt.plot(nvals/nvals[0],marker='o',linestyle='-',label = r'$n/n_0$')
plt.plot(nvsqrdvals/nvsqrdvals[0],marker='o',linestyle='-',label = r'$nv^2/nv^2_0$')
plt.axhline(y=1,linestyle='dashed',color='k')
plt.legend()

allmagperpart = 6.242e18*(allBmags**2)/(2*mu0*alln)
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
brazilhist(1e3*abs(allv[:,0]), z_label = r'$V \ (km/s)$',nbins=80,vmin=0,vmax=800,count=False)

#%%

# Conjunction Stats

# Normalized values for conjunctions
n_psp = 48.7
n_err_psp = 36.8
n_solo = 25.4
n_err_solo = 10.4
n_wind = 23.16
n_err_wind = 17.7


T_psp = 5.47
T_err_psp = 3
T_solo = 24.2
T_err_solo = 9.7
T_wind = 45.2
T_err_wind = 13.3

v_psp = 353.8
v_err_psp = 56.2
v_solo = 602.3
v_err_solo = 118.5
v_wind = 459.7
v_err_wind = 105.3	

nvals = np.array([n_psp,n_solo,n_wind])
n_errs = np.array([n_err_psp,n_err_solo,n_err_wind])
Tvals = np.array([T_psp,T_solo,T_wind])
T_errs = np.array([T_err_psp,T_err_solo,T_err_wind])
vvals = np.array([v_psp,v_solo,v_wind])
v_errs = np.array([v_err_psp,v_err_solo,v_err_wind])

x = np.array([1,2,3])



plt.errorbar(x,nvals,yerr=n_errs,fmt='o',label=r'$n$',capsize=5)
plt.xticks(x,['PSP','SOLO','WIND'])
plt.xlabel('Spacecraft')
plt.legend()
plt.errorbar(x,Tvals,yerr=T_errs,fmt='o',label=r'$T$',capsize=5)
plt.xticks(x,['PSP','SOLO','WIND'])
plt.xlabel('Spacecraft')
plt.legend()

plt.errorbar(x,vvals,yerr=v_errs,fmt='o',label=r'$v$',capsize=5)
plt.xticks(x,['PSP','SOLO','WIND'])
plt.xlabel('Spacecraft')
plt.legend()
plt.xticks(x,['PSP','SOLO','WIND'])
plt.xlabel('Spacecraft')
plt.legend()
# %%
plt.errorbar(x,nvals/nvals[0],yerr=n_errs/nvals[0],fmt='o-',label=r'$R^2 n/n_{psp}$',capsize=5)
plt.xticks(x,['PSP','SOLO','WIND'])
plt.xlabel('Spacecraft')
plt.legend()
plt.errorbar(x,Tvals/Tvals[0],yerr=T_errs/Tvals[0],fmt='o-',label=r'$R^{4/3} T/T_{psp}$',capsize=5)
plt.xticks(x,['PSP','SOLO','WIND'])
plt.xlabel('Spacecraft')
plt.legend()

plt.errorbar(x,vvals**2/(vvals[0]**2),yerr=v_errs**2/vvals[0]**2,fmt='o-',label=r'$v^2/v^2_{psp}$',capsize=5)
plt.xticks(x,['PSP','SOLO','WIND'])
plt.xlabel('Spacecraft')
plt.legend()
plt.xticks(x,['PSP','SOLO','WIND'])
plt.xlabel('Spacecraft')
plt.axhline(y=1,linestyle='dashed',color='k')
plt.ylim(0,11)
plt.legend()



# %%
