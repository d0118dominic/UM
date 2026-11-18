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


# enc23coronalhole_fast = [['2025-03-29/01:20','2025-04-01/16:20']] # Encounter 21
# enc23coronalhole_slow = [['2025-04-04/00:00','2025-04-06/00:00']] # Encounter 21
# enc23coronalhole_bound = [['2025-03-26/00:00','2025-03-29/00:00']] # Encounter 21
enc23coronalhole_full = [['2025-03-23/12:00','2025-04-07/06:00']] # Encounter 21

enc23coronalhole_fast = [['2025-03-29/01:20','2025-04-01/16:20']] # Encounter 21
enc23coronalhole_shoulder = [['2025-03-26/20:00','2025-03-28/12:00']] # Encounter 21
enc23coronalhole_preceding = [['2025-03-23/18:00','2025-03-24/12:00']] # Encounter 21
enc23coronalhole_trailing = [['2025-04-05/00:00','2025-04-07/00:00']] # Encounter 21
enc23coronalhole_rarefaction = [['2025-04-02/00:00','2025-04-04/00:00']] # Encounter 21

# enc23coronalhole_preceding = [['2025-03-23/00:00','2025-03-26/00:00']] # Encounter 21
# 0','2025-03-25/12:00']] # Encounter 21

# enc23coronalhole_rarefaction = [['2025-04-03/00:00','2025-04-05/00:00']] # Encounter 21
# enc23coronalhole_shoulder = [['2025-03-25/20:00','2025-03-27/12:00']] # Encounter 21
eventlist=enc23coronalhole_full




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
	allcrosshelicity_list.append(sigma_c)
	allresidenergy_list.append(sigma_r)

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

allvmags = np.zeros_like(allv[:,0])


allBmags = np.zeros_like(allB[:,0])


for i in range(len(allv)):
	allvmags[i] = np.linalg.norm(allv[i])
	# allBmags[i] = np.linalg.norm(allB[i])
	allBmags[i] = abs(allB[i,0])

allCpar = (alln*allTpar*allBmags**2)/(alln**3)
allCperp = (alln*allTperp)/(alln*allBmags)

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
	n = 1e-6*np.nanmedian(alln)	
	nstd = 1e-6*np.nanstd(alln)
	T = np.nanmedian(allT)/(1.602e-19)
	Tstd = np.nanstd(allT)/(1.602e-19)
	v = 1e-3*np.nanmedian(allvmags)
	vstd = 1e-3*np.nanstd(allvmags)
	B = 1e9*np.nanmedian(allBmags)
	Bstd = 1e9*np.nanstd(allBmags)
	dB = 1e9*np.sqrt(np.nanmedian(alldBsqrdmag))
	dBstd = 1e9*np.sqrt(np.nanstd(alldBsqrdmag))
	ch = np.nanmedian(allcrosshelicity)
	chstd = np.nanstd(allcrosshelicity)
	re = np.nanmedian(allresidenergy)
	restd = np.nanstd(allresidenergy)

	return n,nstd,T,Tstd,v,vstd,B,Bstd,dB,dBstd,ch,chstd,re,restd

def get_solo_anisos():
	betapar = allbeta_par
	Tparperp = allTparperp
	return betapar, Tparperp


def get_solo_CGLvars():
	allCpar = (alln*allTpar*allBmags**2)/(alln**3)
	allCperp = (alln*allTperp)/(alln*allBmags)

	Cpar = np.nanmedian(allCpar)
	Cperp = np.nanmedian(allCperp)
	Cparstd = np.nanstd(allCpar)
	Cperpstd = np.nanstd(allCperp)
	Cparpct = np.array([np.nanpercentile(allCpar,25),np.nanpercentile(allCpar,75)])
	Cperppct = np.array([np.nanpercentile(allCperp,25),np.nanpercentile(allCperp,75)])
	return Cpar,Cparstd,Cperp,Cperpstd,Cparpct,Cperppct


if (eventlist == enc23coronalhole_fast):
	nsolo_fast,nstdsolo_fast,Tsolo_fast,Tstdsolo_fast,vsolo_fast,vstdsolo_fast,Bsolo_fast,Bstdsolo_fast,dBsolo_fast,dBstdsolo_fast,chsolo_fast,chstdsolo_fast,resolo_fast,restdsolo_fast = get_solovals()
	betaparssolo_fast, Tparperpsolo_fast = get_solo_anisos()
	Cparsolo_fast,Cparstdsolo_fast,Cperpsolo_fast,Cperpstdsolo_fast,Cparpctsolo_fast,Cperppctsolo_fast = get_solo_CGLvars()
	npctsolo_fast,Tpctsolo_fast,vpctsolo_fast,Bpctsolo_fast,dBpctsolo_fast,chpctsolo_fast,repctsolo_fast = get_solopcts()
	print("Fast Variables")
elif (eventlist == enc23coronalhole_preceding):
	nsolo_preceding,nstdsolo_preceding,Tsolo_preceding,Tstdsolo_preceding,vsolo_preceding,vstdsolo_preceding,Bsolo_preceding,Bstdsolo_preceding,dBsolo_preceding,dBstdsolo_preceding,chsolo_preceding,chstdsolo_preceding,resolo_preceding,restdsolo_preceding = get_solovals()
	betaparssolo_preceding, Tparperpsolo_preceding = get_solo_anisos()
	Cparsolo_preceding,Cparstdsolo_preceding,Cperpsolo_preceding,Cperpstdsolo_preceding,Cparpctsolo_preceding,Cperppctsolo_preceding = get_solo_CGLvars()
	npctsolo_preceding,Tpctsolo_preceding,vpctsolo_preceding,Bpctsolo_preceding,dBpctsolo_preceding,chpctsolo_preceding,repctsolo_preceding = get_solopcts()
	print("Preceding Variables")
elif (eventlist == enc23coronalhole_trailing):
	nsolo_trailing,nstdsolo_trailing,Tsolo_trailing,Tstdsolo_trailing,vsolo_trailing,vstdsolo_trailing,Bsolo_trailing,Bstdsolo_trailing,dBsolo_trailing,dBstdsolo_trailing,chsolo_trailing,chstdsolo_trailing,resolo_trailing,restdsolo_trailing = get_solovals()
	betaparssolo_trailing, Tparperpsolo_trailing = get_solo_anisos()
	Cparsolo_trailing,Cparstdsolo_trailing,Cperpsolo_trailing,Cperpstdsolo_trailing,Cparpctsolo_trailing,Cperppctsolo_trailing = get_solo_CGLvars()
	npctsolo_trailing,Tpctsolo_trailing,vpctsolo_trailing,Bpctsolo_trailing,dBpctsolo_trailing,chpctsolo_trailing,repctsolo_trailing = get_solopcts()
	print("Trailing Variables")
elif (eventlist == enc23coronalhole_shoulder):
	nsolo_shoulder,nstdsolo_shoulder,Tsolo_shoulder,Tstdsolo_shoulder,vsolo_shoulder,vstdsolo_shoulder,Bsolo_shoulder,Bstdsolo_shoulder,dBsolo_shoulder,dBstdsolo_shoulder,chsolo_shoulder,chstdsolo_shoulder,resolo_shoulder,restdsolo_shoulder = get_solovals()
	betaparssolo_shoulder, Tparperpsolo_shoulder = get_solo_anisos()
	Cparsolo_shoulder,Cparstdsolo_shoulder,Cperpsolo_shoulder,Cperpstdsolo_shoulder,Cparpctsolo_shoulder,Cperppctsolo_shoulder = get_solo_CGLvars()
	npctsolo_shoulder,Tpctsolo_shoulder,vpctsolo_shoulder,Bpctsolo_shoulder,dBpctsolo_shoulder,chpctsolo_shoulder,repctsolo_shoulder = get_solopcts()
	print("Shoulder Variables")
elif (eventlist == enc23coronalhole_rarefaction):
	nsolo_rarefaction,nstdsolo_rarefaction,Tsolo_rarefaction,Tstdsolo_rarefaction,vsolo_rarefaction,vstdsolo_rarefaction,Bsolo_rarefaction,Bstdsolo_rarefaction,dBsolo_rarefaction,dBstdsolo_rarefaction,chsolo_rarefaction,chstdsolo_rarefaction,resolo_rarefaction,restdsolo_rarefaction = get_solovals()
	betaparssolo_rarefaction, Tparperpsolo_rarefaction = get_solo_anisos()
	Cparsolo_rarefaction,Cparstdsolo_rarefaction,Cperpsolo_rarefaction,Cperpstdsolo_rarefaction,Cparpctsolo_rarefaction,Cperppctsolo_rarefaction = get_solo_CGLvars()
	npctsolo_rarefaction,Tpctsolo_rarefaction,vpctsolo_rarefaction,Bpctsolo_rarefaction,dBpctsolo_rarefaction,chpctsolo_rarefaction,repctsolo_rarefaction = get_solopcts()
	print("Rarefaction Variables")
else:
	nsolo,nstdsolo,Tsolo,Tstdsolo,vsolo,vstdsolo,Bsolo,Bstdsolo,dBsolo,dBstdsolo,chsolo,chstdsolo,resolo,restdsolo = get_solovals()
	betaparssolo, Tparperpsolo = get_solo_anisos()
	Cparsolo,Cparstdsolo,Cperpsolo,Cperpstdsolo,Cparpctsolo,Cperppctsolo = get_solo_CGLvars()
	npctsolo,Tpctsolo,vpctsolo,Bpctsolo,dBpctsolo,chpctsolo,repctsolo = get_solopcts()
	print("Variables")

# %%
def brazil_hist2d(z, x_label=r'$\beta_\parallel$', y_label=r'$T_{\perp}/T_{\parallel}$',
     z_label='', title='', mask=None, nbins=50, mincount=10,
     xrange=[-3,3], yrange=[-1,1], cmap='coolwarm', vmin=None, vmax=None):

    mirror_params = [0.77, 0.76, -0.016]  # [a, b, beta0]
    firehose_params = [-1.4, 1.0, -0.11]  # [a, b, beta0]
    cyclotron_params = [0.45, 0.4, -0.0004]  # [a, b, beta0]
    parfirehose_params = [-0.47, 0.53, 0.59]  # [a, b, beta0]

    beta_par_range = np.logspace(xrange[0], xrange[1], 100)
    mirror_threshold = stability_condition(beta_par_range, *mirror_params)
    firehose_threshold = stability_condition(beta_par_range, *firehose_params)
    cyclotron_threshold = stability_condition(beta_par_range, *cyclotron_params)
    parfirehose_threshold = stability_condition(beta_par_range, *parfirehose_params)

    if mask is None:
        mask = np.ones(allTparperp.shape, dtype=bool)

    x = allbeta_par[mask]
    y = 1/allTparperp[mask]
    z = z[mask]

    valmask = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    x_valid = x[valmask]
    y_valid = y[valmask]
    z_valid = z[valmask]

    xbins = np.logspace(xrange[0], xrange[1], nbins)
    ybins = np.logspace(yrange[0], yrange[1], nbins)

    hist, xedges, yedges, _ = binned_statistic_2d(x_valid, y_valid, z_valid, statistic='mean', bins=[xbins, ybins])
    counts, _, _, _ = binned_statistic_2d(x_valid, y_valid, z_valid, statistic='count', bins=[xbins, ybins])
    hist_masked = np.ma.masked_where(counts <= mincount, hist)

    if vmin is None:
        vmin = np.nanpercentile(hist_masked.compressed(), 2) if hist_masked.compressed().size else np.nanpercentile(z_valid, 2)
    if vmax is None:
        vmax = np.nanpercentile(hist_masked.compressed(), 98) if hist_masked.compressed().size else np.nanpercentile(z_valid, 98)

    fig, ax = plt.subplots(figsize=(10, 8))

    im = ax.pcolormesh(xedges, yedges, hist_masked.T, cmap=cmap, vmin=vmin, vmax=vmax)

    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label(z_label, fontdict={'size': 15})

    ax.set_facecolor('lightgray')
    ax.set_xlabel(x_label, fontdict={'size': 15})
    ax.set_ylabel(y_label, fontdict={'size': 15})
    ax.set_title(title, fontdict={'size': 15})
    ax.set_xlim(10**(xrange[0]), 10**(xrange[1]))
    ax.set_ylim(10**(yrange[0]), 10**(yrange[1]))
    ax.axhline(y=1, color='k', linestyle='-')
    ax.axvline(x=1, color='k', linestyle='-')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.plot(beta_par_range, mirror_threshold, label='Mirror', color='k', linestyle='dotted', linewidth=2)
    ax.plot(beta_par_range, firehose_threshold, label='Oblique FH', color='k', linestyle='--', linewidth=2)
    ax.plot(beta_par_range, cyclotron_threshold, label='IC', color='k', linestyle='dotted', linewidth=1)
    ax.plot(beta_par_range, parfirehose_threshold, label='Parallel FH', color='k', linestyle='--', linewidth=1)
    ax.legend()
    plt.show()
    return fig, ax
# %%
# brazil_hist2d(allresidenergy, z_label=r'$\sigma_r$', cmap='seismic',vmin=-1,vmax=1,nbins=50)

# brazil_hist2d(allCpar/np.nanmean(allCpar), z_label=r'$C_\parallel/ \langle C_\parallel \rangle$',vmin=0,vmax=2, cmap='plasma',nbins=60)
# brazil_hist2d(allCperp/np.nanmean(allCperp), z_label=r'$C_\perp/ \langle C_\perp\rangle$', vmin=0,vmax=1.4,cmap='plasma',nbins=60)
#%%
allSp = allT/(alln**(2/3))

brazil_hist2d(allv[:,0],cmap = 'jet',z_label = '$|vR|$')
# brazil_hist2d(alln,cmap = 'jet',z_label = '$n$')
brazil_hist2d(allSp,cmap = 'jet',z_label = '$Sp$',xrange = [-2,1],yrange = [-0.5,1],vmin=4e-23,vmax=8e-23)
# brazil_hist2d(allBmags,cmap = 'jet',z_label = '$B$')

# %%

#%%
# brazil_hist2d(abs(alldn),cmap = 'jet',z_label = '$\delta n$')
# brazil_hist2d(abs(alldT),cmap = 'jet',z_label = '$\delta T$')
brazil_hist2d(abs(alldBsqrdmag),cmap = 'jet',z_label = '$\delta B^2$')
# brazil_hist2d(allcrosshelicity,cmap = 'plasma',z_label = '$|\sigma_c|$')
# brazil_hist2d(abs(alldvsqrdmag),cmap = 'jet',z_label = '$\delta v^2$')
# %%
