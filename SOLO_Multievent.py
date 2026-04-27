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

def get_delta(vec,vec_mean):
	dvec = vec - vec_mean
	dvec_norm = dvec/np.linalg.norm(vec_mean)
	return dvec,dvec_norm
def get_deltascalar(var,var_mean):
	dvar = var - var_mean
	dvar_norm = dvar/var_mean
	return dvar,dvar_norm

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
#%%


ev1 = ['2022-08-10/00:00', '2022-09-02/00:00'] #Full interval
ev2 = ['2022-02-23/00:00', '2022-02-28/00:00'] #Full interval

eventlist=[ev1,ev2]
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

for i in range(len(eventlist)):
    trange = eventlist[i]
    mag_vars = pyspedas.projects.solo.mag(trange=trange, datatype='rtn-normal',get_support_data=True, time_clip=False)
    swa_vars = pyspedas.projects.solo.swa(trange=trange, datatype='pas-grnd-mom',get_support_data=True, time_clip=False)    

    Bvec_name = 'B_RTN'
    PiTensor_name = 'P_RTN'
    ni_name = 'N'
    vi_name = 'V_RTN'
    Tscalar_name = 'T'

    interpvar_name = ni_name
    timeax = pytplot.get_data(interpvar_name).times
    tinterpol(Bvec_name, interpvar_name, newname='B')
    tinterpol(PiTensor_name, interpvar_name, newname='PiTensor')
    tinterpol(ni_name, interpvar_name, newname='ni')
    tinterpol(vi_name, interpvar_name, newname='vi')
    tinterpol(Tscalar_name, interpvar_name, newname='T')

    Bvecs = 1e-9 * reform(pytplot.get_data('B'))
    PiTensor = 1e6 * reform(pytplot.get_data('PiTensor'))
    ni = 1e6 * reform(pytplot.get_data('ni'))
    vi = 1e-3 * reform(pytplot.get_data('vi'))
    T = 1e6 * reform(pytplot.get_data('T'))

    Tpar = np.zeros_like(Bvecs[:, 0])
    Tperp = np.zeros_like(Tpar)
    Pperp = np.zeros_like(Tpar)
    Ppar = np.zeros_like(Tpar)
    P_mag = np.zeros_like(Tpar)
    beta_par = np.zeros_like(Tpar)

    for j in range(len(Bvecs)):  # Fixed: was shadowing outer loop variable i
        Tpar[j], Tperp[j], Ppar[j], Pperp[j] = get_parperps(ni[j], PiTensor[j] / ni[j], Bvecs[j])

        P_mag[j] = get_pm(Bvecs[j])
        beta_par[j] = Ppar[j]/P_mag[j]

    allB_list.append(Bvecs)
    alln_list.append(ni)
    allT_list.append(T)
    allv_list.append(vi)
    # allpositions_chunks.append(positions)
    # allbeta_par_chunks.append(beta_par)
    # allangles_chunks.append(angles)
    allbeta_par_list.append(beta_par)
    allTpar_list.append(Tpar)
    allTperp_list.append(Tperp)
    allTparperp_list.append(Tpar / Tperp)

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
# %%

# %%

def brazilplot(z, x_label='', y_label='', z_label='', title=''):    
    mirror_params = [0.77, 0.76, -0.016]  # [a, b, beta0]
    firehose_params = [-1.4, 1.0, -0.11]  # [a, b, beta0]
    cyclotron_params = [0.45, 0.4, -0.0004]  # [a, b, beta0]
    parfirehose_params = [-0.47, 0.53, 0.59]  # [a, b, beta0]

    beta_par_range = np.logspace(-3, 3, 100)  # Adjust range and resolution as needed
    mirror_threshold = stability_condition(beta_par_range, *mirror_params)
    firehose_threshold = stability_condition(beta_par_range, *firehose_params)
    cyclotron_threshold = stability_condition(beta_par_range, *cyclotron_params)
    parfirehose_threshold = stability_condition(beta_par_range, *parfirehose_params)

    plt.figure(figsize=(10, 8))
    sc = plt.scatter(allbeta_par[mask], 1/allTparperp[mask], c=z[mask], cmap='plasma', s=0.01)
    plt.colorbar(sc, label=z_label)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(title)
    plt.xlim(0.001, 1000)  # Adjust as needed
    plt.ylim(0.1, 10)  # Adjust as needed
    plt.axhline(y=1, color='k', linestyle='-')
    plt.axvline(x=1, color='k', linestyle='-')
    plt.xscale('log')
    plt.yscale('log')
    plt.plot(beta_par_range, mirror_threshold, label='Mirror', color='k',linestyle='dotted', linewidth=2)
    plt.plot(beta_par_range, firehose_threshold, label='Oblique FH', color='k', linestyle='--', linewidth=2)
    plt.plot(beta_par_range, cyclotron_threshold, label='IC', color='k',linestyle='dotted', linewidth=1)
    plt.plot(beta_par_range, parfirehose_threshold, label='Parallel FH', color='k', linestyle='--', linewidth=1)
    plt.legend()
    plt.show()

mask = 1e3*allv[:,0] < 900
brazilplot(1e3*allv[:,0], r'$\beta_\parallel$', r'$T_{\parallel}/T_{\perp}$', r'$v_x$')
# %%
