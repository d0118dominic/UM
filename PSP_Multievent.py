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
def get_vecmean(vec,int):   # Vector mean
	vec1_mean = get_mean(vec[:,0],minutes*meaninterval)
	vec2_mean = get_mean(vec[:,1],minutes*meaninterval)
	vec3_mean = get_mean(vec[:,2],minutes*meaninterval)
	vec_mean = np.zeros_like(vec)
	for i in range(len(vec_mean)): vec_mean[i] = np.array([vec1_mean[i],vec2_mean[i],vec3_mean[i]])
	return vec_mean
def get_mean(var,int):   #Take a timeseries and compute the mean (basically smooth outsmall flucs over some interval)
	box = np.ones(int)/int
	smoothed_var = np.convolve(var,box,mode='same')
	return smoothed_var
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



def get_parperps(n,T,B):  #B ant T tensor coord systems need to match for this
	trace = T[0] + T[1] + T[2]
	term1 = (T[0]*B[0]**2 + T[1]*B[1]**2 + T[2]*B[2]**2)/(np.linalg.norm(B)**2)
	term2 = 2*(T[3]*B[0]*B[1] + T[4]*B[0]*B[2] + T[5]*B[1]*B[2])/(np.linalg.norm(B)**2)
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


ev1 = [['2025-03-22/20:10','2025-03-23/09:20']] 
ev1_reduced = ['2025-03-23/00:00','2025-03-23/09:20']

eventlist=[ev1_reduced]
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

for i in range(len(eventlist)):
	trange = eventlist[i]
	Bfld_vars = pyspedas.projects.psp.fields(trange=trange, level='l2', time_clip=True)
	swp_vars = pyspedas.projects.psp.spi(trange=trange,level='l3',get_support_data=True,time_clip=True)
	qtn_vars = pyspedas.projects.psp.fields(trange=trange,level='l3',datatype='sqtn_rfs_V1V2',time_clip=True)
	
	Bvec_name = 'psp_fld_l2_mag_RTN'
	B_name = Bvec_name
	PiTensor_name = 'P_RTN'
	ni_name = 'electron_density'
	vi_name = 'psp_spi_VEL_RTN_SUN'
	Tscalar_name = 'T'
	Ti_name = Tscalar_name
	phivals_name = 'psp_spi_PHI_VALS'
	ephi_name = 'psp_spi_EFLUX_VS_PHI'

	interpvar_name = ni_name
	timeax = pytplot.get_data(interpvar_name).times
	meaninterval = mean_int(timeax,trange)
	tinterpol(Bvec_name, interpvar_name, newname='B')
	tinterpol(PiTensor_name, interpvar_name, newname='PiTensor')
	tinterpol(ni_name, interpvar_name, newname='ni')
	tinterpol(vi_name, interpvar_name, newname='vi')
	tinterpol(Tscalar_name, interpvar_name, newname='T')

	Bvecs = 1e-9 * reform(pytplot.get_data('B'))
	PiTensor = 1e6 * reform(pytplot.get_data('PiTensor'))
	ni = 1e6 * reform(pytplot.get_data('ni'))
	vi = 1e3 * reform(pytplot.get_data('vi'))
	T = 1.602e-19* reform(pytplot.get_data('T'))

	Tpar = np.zeros_like(Bvecs[:, 0])
	Tperp = np.zeros_like(Tpar)
	Pperp = np.zeros_like(Tpar)
	Ppar = np.zeros_like(Tpar)
	P_mag = np.zeros_like(Tpar)
	beta_par = np.zeros_like(Tpar)
	angle = np.zeros_like(Tpar)



	minutes = 60
	# Means
	Bvecs_mean = get_vecmean(Bvecs,minutes*meaninterval)
	vivecs_mean = get_vecmean(vi,minutes*meaninterval)

	for j in range(len(Bvecs)):  # Fixed: was shadowing outer loop variable i
		Tpar[j], Tperp[j], Ppar[j], Pperp[j] = get_parperps(ni[j], PiTensor[j] / ni[j], Bvecs[j])

		P_mag[j] = get_pm(Bvecs[j])
		beta_par[j] = Ppar[j]/P_mag[j]
		angle[j] = get_angle(Bvecs[j], Bvecs_mean[j])

	
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

allB = np.concatenate(allB_list, axis=0)
allv = np.concatenate(allv_list, axis=0)
alln = np.concatenate(alln_list, axis=0)
allT = np.concatenate(allT_list, axis=0)
# allpositions = np.concatenate(allpositions_list, axis=0)
allbeta_par = np.concatenate(allbeta_par_list, axis=0)
allangles = np.concatenate(allangles_list, axis=0)
allTpar = np.concatenate(allTpar_list, axis=0)
allTperp = np.concatenate(allTperp_list, axis=0)
allTparperp = np.concatenate(allTparperp_list, axis=0)
allangles = np.concatenate(allangles_list, axis=0)
#%%
allvmags = np.zeros_like(allv[:,0])
allBmags = np.zeros_like(allB[:,0])
for i in range(len(allv)):
	allvmags[i] = np.linalg.norm(allv[i])
	allBmags[i] = np.linalg.norm(allB[i])
# %%

allmagperpart = 6.242e18*(allBmags**2)/(2*mu0*alln)

# %%


# mask = allpositions >0
brazilplot(1e3*allbeta_par, r'$\beta_\parallel$', r'$T_{\perp}/T_{\parallel}$', r'$v_x$')
# %%

# mask = 1e-3*allvmags >= 400
brazilhist(1e-3*allvmags, z_label = r'$V \ (km/s)$',nbins=80,vmin=0,vmax=1000,count=False)
# %%
