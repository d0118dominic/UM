
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
	if not isinstance(var[1][0], np.ndarray):
		newvar = np.zeros(len(var[0]))
	elif isinstance(var[1][0][0], np.ndarray):
		newvar = np.zeros([len(var[0]), len(var[1][0]), len(var[1][0][0])])
	elif isinstance(var[1][0], np.ndarray):
		newvar = np.zeros([len(var[0]), len(var[1][0])])
	for i in range(len(var[0]) - 1):
		newvar[i] = var[1][i]
	return newvar
 
def get_delta(vec, vec_mean):
	dvec = vec - vec_mean
	dvec_norm = dvec / np.linalg.norm(vec_mean)
	return dvec, dvec_norm
 
def get_deltascalar(var, var_mean):
	dvar = var - var_mean
	dvar_norm = dvar / var_mean
	return dvar, dvar_norm
 
def stability_condition(betapar, a, b, beta0):
	denom = (betapar - beta0)**b
	Tperppar = 1 + a / denom
	return Tperppar
 
def quickplot():
	pyspedas.tplot([B_name, vi_name, Ti_name, ni_name])
	return
 
def duration(trange):
	from datetime import datetime as dt
	start = dt.strptime(trange[0], '%Y-%m-%d/%H:%M')
	stop = dt.strptime(trange[1], '%Y-%m-%d/%H:%M')
	duration = stop - start
	duration_s = duration.total_seconds()
	duration_m = duration_s / 60  # minutes
	duration_h = duration_m / 60  # hours
	return duration_m
 
def mean_int(timeax, trange):
	steps = len(timeax)
	minutes = duration(trange)
	mean = np.ceil(steps / minutes)
	return int(mean)
 
def get_vecmean(vec, int):   # Vector mean
	vec1_mean = get_mean(vec[:, 0], minutes * meaninterval)
	vec2_mean = get_mean(vec[:, 1], minutes * meaninterval)
	vec3_mean = get_mean(vec[:, 2], minutes * meaninterval)
	vec_mean = np.zeros_like(vec)
	for i in range(len(vec_mean)):
		vec_mean[i] = np.array([vec1_mean[i], vec2_mean[i], vec3_mean[i]])
	return vec_mean
 
def get_mean(var, int):   # Take a timeseries and compute the mean
	box = np.ones(int) / int
	smoothed_var = np.convolve(var, box, mode='same')
	return smoothed_var
 
def get_delta(vec, vec_mean):
	dvec = vec - vec_mean
	dvec_norm = dvec / np.linalg.norm(vec_mean)
	return dvec, dvec_norm
 
def get_deltascalar(var, var_mean):
	dvar = var - var_mean
	dvar_norm = dvar / var_mean
	return dvar, dvar_norm
 
def get_brnorm(B):
	brnorm = B[0] / np.linalg.norm(B)
	return brnorm
 
def get_deflection(B):
	brnorm = get_brnorm(B)
	theta = np.arccos(brnorm) * 180 / np.pi
	return theta
 
def get_angle(vec, meanvec):
	term1 = np.dot(vec, meanvec)
	term2 = np.dot(np.linalg.norm(vec), np.linalg.norm(meanvec))
	term3 = term1 / term2
	term4 = np.arccos(term3) * 180 / np.pi
	return term4
 
def get_pm(B):
	pm = 0.5 * (mu0**-1) * np.linalg.norm(B)**2
	return pm
 
def get_va(B, n, m):
	va = np.linalg.norm(B) / np.sqrt(mu0 * n * m)
	return va
 
def get_vs(T, m):
	vs = np.sqrt(gamma * Z * kb * T / m)
	return vs
 
def get_vth(T, m):
	vth = np.sqrt(kb * T / m)
	return vth
 
def get_vpar(v, B):
	vpar = np.dot(v, B) / np.linalg.norm(B)
 
def get_angle(vec, meanvec):
	term1 = np.dot(vec, meanvec)
	term2 = np.dot(np.linalg.norm(vec), np.linalg.norm(meanvec))
	term3 = term1 / term2
	term4 = np.arccos(term3) * 180 / np.pi
	return term4
 
def get_vecmean(vec, int):   # Vector mean
	vec1_mean = get_mean(vec[:, 0], minutes * meaninterval)
	vec2_mean = get_mean(vec[:, 1], minutes * meaninterval)
	vec3_mean = get_mean(vec[:, 2], minutes * meaninterval)
	vec_mean = np.zeros_like(vec)
	for i in range(len(vec_mean)):
		vec_mean[i] = np.array([vec1_mean[i], vec2_mean[i], vec3_mean[i]])
	return vec_mean
 
def get_parperps(n, T, B):
	# NOTE FOR MAVEN: SWIA onboardsvymom only provides a scalar temperature,
	# not a full pressure tensor. This function expects a 6-component tensor
	# [Txx, Tyy, Tzz, Txy, Txz, Tyz].
	# To use this properly, load SWIA computed moments (swia, datatype='onboardsvymom')
	# and compute a pseudo-tensor, or use STATIC L2 data for the full ion tensor.
	# Here we treat T as a diagonal isotropic tensor: [T, T, T, 0, 0, 0].
	trace = T[0] + T[1] + T[2]
	term1 = (T[0] * B[0]**2 + T[1] * B[1]**2 + T[2] * B[2]**2) / (np.linalg.norm(B)**2)
	term2 = 2 * (T[3] * B[0] * B[1] + T[4] * B[0] * B[2] + T[5] * B[1] * B[2]) / (np.linalg.norm(B)**2)
	Tpar = term1 + term2
	Tperp = 0.5 * (trace - Tpar)
	Ppar = n * Tpar
	Pperp = n * Tperp
	return Tpar, Tperp, Ppar, Pperp
 
def get_crosshelicity(v, B, n, m):  # vector v & B
	z_plus = v + B / np.sqrt(n * m * mu0)
	z_minus = v - B / np.sqrt(n * m * mu0)
	term1 = np.linalg.norm(z_plus)**2 - np.linalg.norm(z_minus)**2
	term2 = np.linalg.norm(z_plus)**2 + np.linalg.norm(z_minus)**2
	sigma_c = term1 / term2
	return sigma_c
 
def get_residenergy(v, B, n, m):  # vector dv & dB (Alfven units??)
	term1 = np.linalg.norm(v)**2 - np.linalg.norm(B / np.sqrt(n * m * mu0))**2
	term2 = np.linalg.norm(v)**2 + np.linalg.norm(B / np.sqrt(n * m * mu0))**2
	sigma_r = term1 / term2
	return sigma_r
 
 
#%%
 
# -----------------------------------------------------------------
# MAVEN event time range
# MAVEN operates in Mars orbit; data is available from 2014 onward.
# Replace with your event of interest.
# -----------------------------------------------------------------
ev1 = ['2025-03-28/19:30', '2025-04-03/14:50'] #Full interval
# ev1 = ['2025-03-29/12:00', '2025-03-30/00:00'] #Full interval

eventlist = [ev1]
 
#%%
 
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
 
# for i in range(len(eventlist)):
trange = eventlist[0]
trange=['2025-03-01', '2025-03-02']

mag_vars = pyspedas.projects.maven.mag(trange=trange)
swia_vars = pyspedas.projects.maven.swia(trange=trange, datatype='onboardsvymom')

# MAVEN variable names
Bvec_name  = 'OB_B'            # MAG magnetic field vector, payload frame [nT]
B_name = Bvec_name
ni_name    = 'density'         # SWIA onboard moment: ion density [cm^-3]
vi_name    = 'velocity_mso'    # SWIA onboard moment: ion velocity, MSO frame [km/s]
Tscalar_name = 'temperature'   # SWIA onboard moment: scalar ion temperature [eV]
Ti_name    = Tscalar_name

# MAVEN SWIA does not provide a pressure tensor in onboardsvymom.
# We will build an isotropic pseudo-tensor [Txx,Tyy,Tzz,Txy,Txz,Tyz]
# from the scalar temperature after interpolation (see below).
# If you have STATIC data or SWIA computed 3D moments with a tensor,
# replace 'PiTensor_name' with the appropriate variable name.
PiTensor_name = None  # No native tensor from SWIA onboardsvymom

interpvar_name = ni_name
timeax = pytplot.get_data(interpvar_name).times
meaninterval = mean_int(timeax, trange)

tinterpol(Bvec_name, interpvar_name, newname='B')
tinterpol(ni_name, interpvar_name, newname='ni')
tinterpol(vi_name, interpvar_name, newname='vi')
tinterpol(Tscalar_name, interpvar_name, newname='T')

# Unit conversions
Bvecs = 1e-9 * reform(pytplot.get_data('B'))          # T
ni    = 1e6  * reform(pytplot.get_data('ni'))          # m^-3  (from cm^-3)
vi    = 1e3  * reform(pytplot.get_data('vi'))          # m/s   (from km/s)
T     = 1.602e-19 * reform(pytplot.get_data('T'))      # J     (from eV)

# Build isotropic pseudo pressure tensor from scalar T: [Txx,Tyy,Tzz,Txy,Txz,Tyz]
# This treats the plasma as isotropic — replace with real tensor if available.
npts = len(ni)
PiTensor = np.zeros((npts, 6))
PiTensor[:, 0] = T  # Txx = T
PiTensor[:, 1] = T  # Tyy = T
PiTensor[:, 2] = T  # Tzz = T
# Off-diagonal elements remain zero (isotropic assumption)

Tpar  = np.zeros_like(Bvecs[:, 0])
Tperp = np.zeros_like(Tpar)
Pperp = np.zeros_like(Tpar)
Ppar  = np.zeros_like(Tpar)
P_mag = np.zeros_like(Tpar)
beta_par = np.zeros_like(Tpar)
angle = np.zeros_like(Tpar)
dB, dB_norm = np.zeros_like(Bvecs), np.zeros_like(Bvecs)
dv, dv_norm = np.zeros_like(vi), np.zeros_like(vi)
dvsqrdmag, dvsqrdmag_norm = np.zeros_like(Tpar), np.zeros_like(Tpar)
dBsqrdmag, dBsqrdmag_norm = np.zeros_like(Tpar), np.zeros_like(Tpar)
dn, dn_norm = np.zeros_like(ni), np.zeros_like(ni)
dT, dT_norm = np.zeros_like(T), np.zeros_like(T)

minutes = 60
# Means
Bvecs_mean  = get_vecmean(Bvecs, minutes * meaninterval)
vivecs_mean = get_vecmean(vi, minutes * meaninterval)
n_mean = get_mean(ni, minutes * meaninterval)
T_mean = get_mean(T, minutes * meaninterval)

for j in range(len(Bvecs)):
    Tpar[j], Tperp[j], Ppar[j], Pperp[j] = get_parperps(ni[j], PiTensor[j] / ni[j], Bvecs[j])

    P_mag[j]   = get_pm(Bvecs[j])
    beta_par[j] = Ppar[j] / P_mag[j]
    angle[j]   = get_angle(Bvecs[j], Bvecs_mean[j])
    dB[j], dB_norm[j] = get_delta(Bvecs[j], Bvecs_mean[j])
    dv[j], dv_norm[j] = get_delta(vi[j], vivecs_mean[j])
    dn[j], dn_norm[j] = get_deltascalar(ni[j], n_mean[j])
    dT[j], dT_norm[j] = get_deltascalar(T[j], T_mean[j])
    dvsqrdmag[j], dvsqrdmag_norm[j] = get_deltascalar(np.linalg.norm(vi[j])**2, np.linalg.norm(vivecs_mean[j])**2)
    dBsqrdmag[j], dBsqrdmag_norm[j] = get_deltascalar(np.linalg.norm(Bvecs[j])**2, np.linalg.norm(Bvecs_mean[j])**2)

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
 
allB         = np.concatenate(allB_list, axis=0)
allv         = np.concatenate(allv_list, axis=0)
alln         = np.concatenate(alln_list, axis=0)
allT         = np.concatenate(allT_list, axis=0)
allbeta_par  = np.concatenate(allbeta_par_list, axis=0)
allangles    = np.concatenate(allangles_list, axis=0)
allTpar      = np.concatenate(allTpar_list, axis=0)
allTperp     = np.concatenate(allTperp_list, axis=0)
allTparperp  = np.concatenate(allTparperp_list, axis=0)
allangles    = np.concatenate(allangles_list, axis=0)
alldB        = np.concatenate(alldB_list, axis=0)
alldv        = np.concatenate(alldv_list, axis=0)
alldn        = np.concatenate(alldn_list, axis=0)
alldT        = np.concatenate(alldT_list, axis=0)
alldvsqrdmag = np.concatenate(alldvsqrdmag_list, axis=0)
alldBsqrdmag = np.concatenate(alldBsqrdmag_list, axis=0)
 
#%%
allvmags = np.zeros_like(allv[:, 0])
allBmags = np.zeros_like(allB[:, 0])
for i in range(len(allv)):
	allvmags[i] = np.linalg.norm(allv[i])
	allBmags[i] = np.linalg.norm(allB[i])
 
#%%
allmagperpart = 6.242e18 * (allBmags**2) / (2 * mu0 * alln)