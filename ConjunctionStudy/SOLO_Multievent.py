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

def get_delta(vec,vec_mean):
	dvec = vec - vec_mean
	dvec_norm = dvec/np.linalg.norm(vec_mean)
	return dvec,dvec_norm
def get_deltascalar(var,var_mean):
	dvar = var - var_mean
	dvar_norm = dvar/var_mean
	return dvar,dvar_norm

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


ev1 = ['2025-03-28/19:30', '2025-04-03/14:50'] #Full interval
ev1_reduced = ['2025-03-28/19:30', '2025-04-01/00:00'] #Full interval


# enc23coronalhole_fast = [['2025-03-30/00:00','2025-04-01/00:00']] # Encounter 21
# enc23coronalhole_slow = [['2025-03-23/00:40','2025-03-24/05:20']] # Encounter 21


enc23coronalhole_fast = [['2025-03-29/01:20','2025-04-01/16:20']] # Encounter 21
enc23coronalhole_slow = [['2025-04-04/00:00','2025-04-06/00:00']] # Encounter 21
enc23coronalhole_bound = [['2025-03-26/00:00','2025-03-29/00:00']] # Encounter 21
enc23coronalhole_full = [['2025-03-23/12:00','2025-04-07/06:00']] # Encounter 21

enc23coronalhole_fast = [['2025-03-29/01:20','2025-04-01/16:20']] # Encounter 21
enc23coronalhole_shoulder = [['2025-03-26/20:00','2025-03-28/12:00']] # Encounter 21
enc23coronalhole_preceding = [['2025-03-23/18:00','2025-03-24/12:00']] # Encounter 21
enc23coronalhole_trailing = [['2025-04-05/00:00','2025-04-07/00:00']] # Encounter 21
enc23coronalhole_rarefaction = [['2025-04-02/00:00','2025-04-04/00:00']] # Encounter 21

# enc23coronalhole_preceding = [['2025-03-22/18:0
# 0','2025-03-25/12:00']] # Encounter 21

eventlist=enc23coronalhole_trailing




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
allcrosshelicity_list = []
allresidenergy_list = []

for i in range(len(eventlist)):
	trange = eventlist[i]
	mag_vars = pyspedas.projects.solo.mag(trange=trange, datatype='rtn-normal',get_support_data=False, time_clip=True,no_update=True)
	swa_vars = pyspedas.projects.solo.swa(trange=trange, datatype='pas-grnd-mom',get_support_data=False, time_clip=True,no_update=True)    

	Bvec_name = 'B_RTN'
	B_name = Bvec_name
	PiTensor_name = 'P_RTN'
	ni_name = 'N'
	vi_name = 'V_RTN'
	Tscalar_name = 'T'
	Ti_name = Tscalar_name

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
	dB,dB_norm = np.zeros_like(Bvecs), np.zeros_like(Bvecs)
	dv,dv_norm = np.zeros_like(vi), np.zeros_like(vi)
	dvsqrdmag,dvsqrdmag_norm = np.zeros_like(Tpar), np.zeros_like(Tpar)
	dBsqrdmag,dBsqrdmag_norm = np.zeros_like(Tpar), np.zeros_like(Tpar)
	dn,dn_norm = np.zeros_like(ni), np.zeros_like(ni)
	dT,dT_norm = np.zeros_like(T), np.zeros_like(T)
	sigma_c,sigma_r = np.zeros_like(ni),np.zeros_like(ni)



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
		sigma_r[j] = get_residenergy(dv[j],dB[j],ni[j],mi)
		sigma_c[j] = get_crosshelicity(dv[j],dB[j],ni[j],mi)

	
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
	allcrosshelicity_list.append(abs(sigma_c))
	allresidenergy_list.append(abs(sigma_r))

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
alldB = np.concatenate(alldB_list, axis=0)
alldv = np.concatenate(alldv_list, axis=0)
alldn = np.concatenate(alldn_list, axis=0)
alldT = np.concatenate(alldT_list, axis=0)
alldvsqrdmag = np.concatenate(alldvsqrdmag_list, axis=0)
alldBsqrdmag = np.concatenate(alldBsqrdmag_list, axis=0)
allcrosshelicity = np.concatenate(allcrosshelicity_list,axis=0)
allresidenergy = np.concatenate(allresidenergy_list,axis=0)
allCpar = (alln*allTpar*allBmags**2)/(alln**3)
allCperp = (alln*allTperp)/(alln*allBmags)

allvmags = np.zeros_like(allv[:,0])
allBmags = np.zeros_like(allB[:,0])
for i in range(len(allv)):
	allvmags[i] = np.linalg.norm(allv[i])
	allBmags[i] = np.linalg.norm(allB[i])


# def make_errors(var, upper,lower):
# 	low = var-lower
# 	upp = upper-var
# 	return [low,upp]


def get_solopcts():
	nupper = 1e-6*np.nanpercentile(alln,75)
	nlower = 1e-6*np.nanpercentile(alln,25)
	npct = np.array([nlower,nupper])

	Tupper = np.nanpercentile(allT,75)/(1.602e-19)
	Tlower = np.nanpercentile(allT,25)/(1.602e-19)
	Tpct = np.array([Tlower,Tupper])

	vupper = 1e-3*np.nanpercentile(allvmags,75)
	vlower = 1e-3*np.nanpercentile(allvmags,25)
	vpct = np.array([vlower,vupper])
	
	Bupper = 1e9*np.nanpercentile(allBmags,75)
	Blower = 1e9*np.nanpercentile(allBmags,25)
	Bpct = np.array([Blower,Bupper])

	
	dBupper = 1e9*np.nanpercentile(np.sqrt(alldBsqrdmag),75)
	dBlower = 1e9*np.nanpercentile(np.sqrt(alldBsqrdmag),25)
	dBpct = np.array([dBlower,dBupper])

	chupper = np.nanpercentile(allcrosshelicity,75)
	chlower = np.nanpercentile(allcrosshelicity,25)
	chpct = np.array([chlower,chupper])
	
	reupper = np.nanpercentile(allresidenergy,75)
	relower = np.nanpercentile(allresidenergy,25)
	repct = np.array([relower,reupper])
	return npct,Tpct,vpct,Bpct,dBpct,chpct,repct



def get_solovals():
	n = 1e-6*np.nanmean(alln)	
	nstd = 1e-6*np.nanstd(alln)
	T = np.nanmean(allT)/(1.602e-19)
	Tstd = np.nanstd(allT)/(1.602e-19)
	v = 1e-3*np.nanmean(allvmags)
	vstd = 1e-3*np.nanstd(allvmags)
	B = 1e9*np.nanmean(allBmags)
	Bstd = 1e9*np.nanstd(allBmags)
	dB = 1e9*np.sqrt(np.nanmean(alldBsqrdmag))
	dBstd = 1e9*np.sqrt(np.nanstd(alldBsqrdmag))
	ch = np.nanmean(allcrosshelicity)
	chstd = np.nanstd(allcrosshelicity)
	re = np.nanmean(allresidenergy)
	restd = np.nanstd(allresidenergy)

	return n,nstd,T,Tstd,v,vstd,B,Bstd,dB,dBstd,ch,chstd,re,restd

def get_solo_anisos():
	betapar = allbeta_par
	Tparperp = allTparperp
	return betapar, Tparperp


def get_CGLvars():
	Cpar = np.nanmean(allCpar)
	Cparpct = np.array([np.nanpercentile(Cpar,25),np.nanpercentile(Cpar,75)])
	Cperppct = np.array([np.nanpercentile(Cperp,25),np.nanpercentile(Cperp,75)])
	Cperp = np.nanmean(allCperp)
	Cparstd = np.nanstd(allCpar)
	Cperpstd = np.nanstd(allCperp)
	return Cpar,Cparstd,Cperp,Cperpstd
# if (eventlist == enc23coronalhole_fast):
# 	nsolo_fast,nstdsolo_fast,Tsolo_fast,Tstdsolo_fast,vsolo_fast,vstdsolo_fast,Bsolo_fast,Bstdsolo_fast,dBsolo_fast,dBstdsolo_fast = get_solovals()
# 	betaparssolo_fast, Tparperpsolo_fast = get_solo_anisos()
# 	print("Fast Variables")
# elif (eventlist == enc23coronalhole_slow):
# 	nsolo_slow,nstdsolo_slow,Tsolo_slow,Tstdsolo_slow,vsolo_slow,vstdsolo_slow,Bsolo_slow,Bstdsolo_slow,dBsolo_slow,dBstdsolo_slow = get_solovals()
# 	betaparssolo_slow, Tparperpsolo_slow = get_solo_anisos()
# 	print("Slow Variables")
# else:
# 	nsolo,nstdsolo,Tsolo,Tstdsolo,vsolo,vstdsolo,Bsolo,Bstdsolo,dBsolo,dBstdsolo = get_solovals()
# 	betaparssolo, Tparperpsolo = get_solo_anisos()
# 	print("Variables")
if (eventlist == enc23coronalhole_fast):
	nsolo_fast,nstdsolo_fast,Tsolo_fast,Tstdsolo_fast,vsolo_fast,vstdsolo_fast,Bsolo_fast,Bstdsolo_fast,dBsolo_fast,dBstdsolo_fast,chsolo_fast,chstdsolo_fast,resolo_fast,restdsolo_fast = get_solovals()
	betaparssolo_fast, Tparperpsolo_fast = get_solo_anisos()
	Cparsolo_fast,Cparstdsolo_fast,Cperpsolo_fast,Cperpstdsolo_fast = get_CGLvars()
	print("Fast Variables")
elif (eventlist == enc23coronalhole_preceding):
	nsolo_preceding,nstdsolo_preceding,Tsolo_preceding,Tstdsolo_preceding,vsolo_preceding,vstdsolo_preceding,Bsolo_preceding,Bstdsolo_preceding,dBsolo_preceding,dBstdsolo_preceding,chsolo_preceding,chstdsolo_preceding,resolo_preceding,restdsolo_preceding = get_solovals()
	betaparssolo_preceding, Tparperpsolo_preceding = get_solo_anisos()
	Cparsolo_preceding,Cparstdsolo_preceding,Cperpsolo_preceding,Cperpstdsolo_preceding = get_CGLvars()
	print("Preceding Variables")
elif (eventlist == enc23coronalhole_trailing):
	nsolo_trailing,nstdsolo_trailing,Tsolo_trailing,Tstdsolo_trailing,vsolo_trailing,vstdsolo_trailing,Bsolo_trailing,Bstdsolo_trailing,dBsolo_trailing,dBstdsolo_trailing,chsolo_trailing,chstdsolo_trailing,resolo_trailing,restdsolo_trailing = get_solovals()
	betaparssolo_trailing, Tparperpsolo_trailing = get_solo_anisos()
	Cparsolo_trailing,Cparstdsolo_trailing,Cperpsolo_trailing,Cperpstdsolo_trailing = get_CGLvars()
	print("Trailing Variables")
elif (eventlist == enc23coronalhole_shoulder):
	nsolo_shoulder,nstdsolo_shoulder,Tsolo_shoulder,Tstdsolo_shoulder,vsolo_shoulder,vstdsolo_shoulder,Bsolo_shoulder,Bstdsolo_shoulder,dBsolo_shoulder,dBstdsolo_shoulder,chsolo_shoulder,chstdsolo_shoulder,resolo_shoulder,restdsolo_shoulder = get_solovals()
	betaparssolo_shoulder, Tparperpsolo_shoulder = get_solo_anisos()
	Cparsolo_shoulder,Cparstdsolo_shoulder,Cperpsolo_shoulder,Cperpstdsolo_shoulder = get_CGLvars()
	print("Shoulder Variables")
elif (eventlist == enc23coronalhole_rarefaction):
	nsolo_rarefaction,nstdsolo_rarefaction,Tsolo_rarefaction,Tstdsolo_rarefaction,vsolo_rarefaction,vstdsolo_rarefaction,Bsolo_rarefaction,Bstdsolo_rarefaction,dBsolo_rarefaction,dBstdsolo_rarefaction,chsolo_rarefaction,chstdsolo_rarefaction,resolo_rarefaction,restdsolo_rarefaction = get_solovals()
	betaparssolo_rarefaction, Tparperpsolo_rarefaction = get_solo_anisos()
	Cparsolo_rarefaction,Cparstdsolo_rarefaction,Cperpsolo_rarefaction,Cperpstdsolo_rarefaction = get_CGLvars()
	print("Rarefaction Variables")
else:
	nsolo,nstdsolo,Tsolo,Tstdsolo,vsolo,vstdsolo,Bsolo,Bstdsolo,dBsolo,dBstdsolo,chsolo,chstdsolo,resolo,restdsolo = get_solovals()
	betaparssolo, Tparperpsolo = get_solo_anisos()
	Cparsolo,Cparstdsolo,Cperpsolo,Cperpstdsolo = get_CGLvars()
	print("Variables")
# allmagperpart = 6.242e18*(allBmags**2)/(2*mu0*alln)

# %%
