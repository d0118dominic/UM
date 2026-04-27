#%%

import pyspedas
import numpy as np
import matplotlib.pyplot as plt
import pytplot
from pytplot import get_data, store_data,timespan,tplot_options
from pyspedas import tplot,tplot_options
from pyspedas import tinterpol,options
from pyspedas.mms import mec,fgm,fpi,edp,scm
 
me = 9.1094e-31 #kg
mi = 1837*me
mu0 = 1.2566370e-06  #;m kg / C^2
mu_0 = mu0
eps0 = 8.85e-12   # C^2/Nm^2
e = 1.602e-19 #C
Z = 1 # 1 for H+, 2 for He2+
gamma = 5/3
kb = 1.380649e-23
 
# %%
 
# Get Data
trange = ['2017-07-11/22:34:00', '2017-07-11/22:34:06'] #Classic EDR
trange = ['2017-08-10/12:18:00', '2017-08-10/12:19:00']
# trange = ['2017-06-17/20:24:00', '2017-06-17/20:24:15'] #Generator Event
# trange = ['2016-12-09/09:03:54', '2016-12-09/09:03:55'] # E-only Event
# trange = ['2016-12-09/09:02:00', '2016-12-09/09:05:00'] # 
# trange = ['2017-07-11/22:33:45', '2017-07-11/22:34:15'] #
trange = ['2015-10-16/13:07:00', '2015-10-16/13:07:05'] # #Oct 16th EDR
 
probe  = 2
trange = trange
 
fgm_vars = fgm(probe = probe, data_rate = 'brst', trange=trange,time_clip=True)
edp_vars = edp(probe = probe,data_rate = 'brst',trange=trange,time_clip=True) 
fpi_vars = fpi(probe = probe,data_rate = 'brst',trange=trange,time_clip=True)
scm_vars = scm(probe = probe, data_rate = 'brst', trange=trange,time_clip=True)
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
 
 
def get_critscales(E,B,m,n,q,T):
    tau_KS = (1/q)*((2*B/(n*mu_0))**(1/3))*((m/E)**(2/3))     # tau_KS (kinetic ETD vs Poynting)
    l_KS   = (1/q)*((B/(n*mu_0))**(2/3))*((m/(2*E))**(1/3))   # l_KS
    tau_HS = (1/q)*(2*m*B/(5*mu_0*n*T))                        # tau_HS (thermal ETD vs Poynting)
    l_HS   = (1/q)*(2*m*E)*(B/(5*mu_0*n*T))**2                 # l_HS (T already in Joules, no kb)
    tau_D  = m/(q*B)                                            # tau_D (inverse cyclotron frequency)
    l_D    = m*E/(2*q*B**2)                                     # l_D
    tau_KH = (1/q)*np.sqrt(5*m*T/E**2)                         # tau_KH (kinetic vs thermal ETD)
    l_KH   = (1/q)*(5*T/(2*E))                                  # l_KH (T already in Joules, no kb)
    return tau_KS, l_KS, tau_HS, l_HS, tau_D, l_D, tau_KH, l_KH
 
 
# Field names variables
B_name = 'mms' + str(probe) + '_fgm_b_gse_brst_l2'
E_name = 'mms' + str(probe) + '_edp_dce_gse_brst_l2'
vi_name = 'mms' + str(probe) + '_' + 'dis' + '_bulkv_gse_brst'
ve_name = 'mms' + str(probe) + '_' + 'des' + '_bulkv_gse_brst'
B_name = 'mms' + str(probe) + '_fgm_b_gse_brst_l2'
E_name = 'mms' + str(probe) + '_edp_dce_gse_brst_l2'
ne_name = 'mms' + str(probe) + '_' + 'des' + '_numberdensity_brst'
ni_name = 'mms' + str(probe) + '_' + 'dis' + '_numberdensity_brst'
scm_name = 'mms' + str(probe) + '_scm_acb_gse_scb_brst_l2'
Pi_name = 'mms' + str(probe) + '_dis_prestensor_gse_brst'
Pe_name = 'mms' + str(probe) + '_des_prestensor_gse_brst'
Ti_name = 'mms' + str(probe) + '_dis_temptensor_gse_brst'
Te_name = 'mms' + str(probe) + '_des_temptensor_gse_brst'
 
##%%
var_name = B_name
timeax = pytplot.get_data(var_name).times
 
tinterpol(B_name,var_name, newname='B')
tinterpol(E_name,var_name, newname='E')
tinterpol(vi_name,var_name, newname='vi')
tinterpol(ve_name,var_name, newname='ve')
tinterpol(scm_name,var_name, newname='B_scm')
tinterpol(Pi_name,var_name, newname='Pi')
tinterpol(Pe_name,var_name, newname='Pe')
tinterpol(Ti_name,var_name, newname='Ti')
tinterpol(Te_name,var_name, newname='Te')
tinterpol(ni_name,var_name, newname='ni')
tinterpol(ne_name,var_name, newname='ne')
 
B = 1e-9*reform(get_data('B'))
E = 1e-3*reform(get_data('E'))
vi = 1e3*reform(get_data('vi'))
ve = 1e3*reform(get_data('ve'))
B_scm = 1e-9*reform(get_data('B_scm'))
ni = 1e6*reform(get_data('ni'))
ne = 1e6*reform(get_data('ne'))
Pi = 1e-9*reform(get_data('Pi'))
Pe = 1e-9*reform(get_data('Pe'))
Ti = 1.602e-19*reform(get_data('Ti'))
Te = 1.602e-19*reform(get_data('Te'))
 
#%%
# Initialize arrays — electrons
tau_KS_e = np.zeros(len(B))
l_KS_e   = np.zeros(len(B))
tau_HS_e = np.zeros(len(B))
l_HS_e   = np.zeros(len(B))
tau_D_e  = np.zeros(len(B))
l_D_e    = np.zeros(len(B))
tau_KH_e = np.zeros(len(B))
l_KH_e   = np.zeros(len(B))
 
# Initialize arrays — ions
tau_KS_i = np.zeros(len(B))
l_KS_i   = np.zeros(len(B))
tau_HS_i = np.zeros(len(B))
l_HS_i   = np.zeros(len(B))
tau_D_i  = np.zeros(len(B))
l_D_i    = np.zeros(len(B))
tau_KH_i = np.zeros(len(B))
l_KH_i   = np.zeros(len(B))
 
for i in range(len(B)):
    tau_KS_e[i], l_KS_e[i], tau_HS_e[i], l_HS_e[i], tau_D_e[i], l_D_e[i], tau_KH_e[i], l_KH_e[i] = get_critscales(np.linalg.norm(E[i]),np.linalg.norm(B[i,:-1]),me,ne[i],e,np.trace(Te[i])/3)
    tau_KS_i[i], l_KS_i[i], tau_HS_i[i], l_HS_i[i], tau_D_i[i], l_D_i[i], tau_KH_i[i], l_KH_i[i] = get_critscales(np.linalg.norm(E[i]),np.linalg.norm(B[i,:-1]),mi,ni[i],e,np.trace(Ti[i])/3)
 
# #%%
fig, ax = plt.subplots(2,1,figsize=(10,8))
ax[0].plot(tau_KS_e, label=r'$\tau_{KSe}$', color='blue')
ax[0].plot(tau_D_e,  label=r'$\tau_{De}$',  color='cyan')
ax[0].plot(tau_KS_i, label=r'$\tau_{KSi}$', color='red')
ax[0].plot(tau_D_i,  label=r'$\tau_{Di}$',  color='orange')
ax[0].set_yscale('log')
ax[0].set_ylabel('Time Scale (s)')
ax[0].legend()
ax[1].plot(1e-3*l_KS_e, label=r'$\ell_{KSe}$', color='blue')
ax[1].plot(1e-3*l_D_e,  label=r'$\ell_{De}$',  color='cyan')
ax[1].plot(1e-3*l_KS_i, label=r'$\ell_{KSi}$', color='red')
ax[1].plot(1e-3*l_D_i,  label=r'$\ell_{Di}$',  color='orange')
ax[1].set_yscale('log')
ax[1].set_ylabel('Length Scale (km)')
ax[1].legend()
ax[1].set_xlabel('Time Index')
 
# Store all scales
store_data('tau_KS_e', data={'x': timeax, 'y': tau_KS_e})
store_data('tau_KS_i', data={'x': timeax, 'y': tau_KS_i})
store_data('l_KS_e',   data={'x': timeax, 'y': l_KS_e})
store_data('l_KS_i',   data={'x': timeax, 'y': l_KS_i})
 
store_data('tau_HS_e', data={'x': timeax, 'y': tau_HS_e})
store_data('tau_HS_i', data={'x': timeax, 'y': tau_HS_i})
store_data('l_HS_e',   data={'x': timeax, 'y': l_HS_e})
store_data('l_HS_i',   data={'x': timeax, 'y': l_HS_i})
 
store_data('tau_D_e',  data={'x': timeax, 'y': tau_D_e})
store_data('tau_D_i',  data={'x': timeax, 'y': tau_D_i})
store_data('l_D_e',    data={'x': timeax, 'y': l_D_e})
store_data('l_D_i',    data={'x': timeax, 'y': l_D_i})
 
store_data('tau_KH_e', data={'x': timeax, 'y': tau_KH_e})
store_data('tau_KH_i', data={'x': timeax, 'y': tau_KH_i})
store_data('l_KH_e',   data={'x': timeax, 'y': l_KH_e})
store_data('l_KH_i',   data={'x': timeax, 'y': l_KH_i})
 
# #%%
fig, ax = plt.subplots(1,1,figsize=(10,8))
ax.plot(1/tau_KS_e, label=r'$\omega_{KSe}$', color='blue')
ax.plot(1/tau_D_e,  label=r'$\omega_{De}$',  color='cyan')
ax.plot(1/tau_KS_i, label=r'$\omega_{KSi}$', color='red')
ax.plot(1/tau_D_i,  label=r'$\omega_{Di}$',  color='orange')
ax.legend()
ax.set_yscale('log')
ax.set_ylabel('Frequency ($s^{-1}$)')
 
# %%
store_data('omega_KS_e', data={'x': timeax, 'y': 1/tau_KS_e})
store_data('omega_KS_i', data={'x': timeax, 'y': 1/tau_KS_i})
store_data('omega_D_e',  data={'x': timeax, 'y': 1/tau_D_e})
store_data('omega_D_i',  data={'x': timeax, 'y': 1/tau_D_i})
 
# Combined electron lengthscales
electron_lengthscales = np.array([l_KS_e, l_D_e, l_HS_e, l_KH_e])
ion_lengthscales      = np.array([l_KS_i, l_D_i, l_HS_i, l_KH_i])
store_data('le', data={'x': timeax, 'y': electron_lengthscales.T})
store_data('li', data={'x': timeax, 'y': ion_lengthscales.T})
options('le', 'legend_names', [r'$\ell_{KSe}$', r'$\ell_{De}$', r'$\ell_{HSe}$', r'$\ell_{KHe}$'])
options('le', 'ytitle', 'Length Scale (m)')
options('le', 'thick', [1,1,1,1])
options('le', 'color', ['b','r','g','m'])
options('le', 'ylog', True)
 
options('li', 'legend_names', [r'$\ell_{KSi}$', r'$\ell_{Di}$', r'$\ell_{HSi}$', r'$\ell_{KHi}$'])
options('li', 'ytitle', 'Length Scale (m)')
options('li', 'thick', [1,1,1,1])
options('li', 'color', ['b','r','g','m'])
options('li', 'ylog', True)
 
# Combined electron timescales
electron_timescales = np.array([tau_KS_e, tau_D_e, tau_HS_e, tau_KH_e])
ion_timescales      = np.array([tau_KS_i, tau_D_i, tau_HS_i, tau_KH_i])
 
store_data('te', data={'x': timeax, 'y': electron_timescales.T})
store_data('ti', data={'x': timeax, 'y': ion_timescales.T})
options('te', 'legend_names', [r'$\tau_{KSe}$', r'$\tau_{De}$', r'$\tau_{HSe}$', r'$\tau_{KHe}$'])
options('te', 'ytitle', 'Time Scale (s)')
options('te', 'thick', [1,1,1,1])
options('te', 'color', ['b','r','g','m'])
options('te', 'ylog', True)
 
options('ti', 'legend_names', [r'$\tau_{KSi}$', r'$\tau_{Di}$', r'$\tau_{HSi}$', r'$\tau_{KHi}$'])
options('ti', 'ytitle', 'Time Scale (s)')
options('ti', 'thick', [1,1,1,1])
options('ti', 'color', ['b','r','g','m'])
options('ti', 'ylog', True)
 
options(['ti','te','li','le'], 'legend_size', 22)
 
tplot('le')
tplot('li')
tplot('te')
tplot('ti')
 
# %%

electron_lengthscales = np.array([l_KS_e, l_D_e])
ion_lengthscales      = np.array([l_KS_i, l_D_i])
store_data('le_simple', data={'x': timeax, 'y': electron_lengthscales.T})
store_data('li_simple', data={'x': timeax, 'y': ion_lengthscales.T})
options('le_simple', 'legend_names', [r'$\ell_{KSe}$', r'$\ell_{De}$'])
options('le_simple', 'ytitle', 'Length Scale (m)')
options('le_simple', 'thick', [1,1])
options('le_simple', 'color', ['b','r'])                    
options('le_simple', 'ylog', True)
options('li_simple', 'legend_names', [r'$\ell_{KSi}$', r'$\ell_{Di}$'])
options('li_simple', 'ytitle', 'Length Scale (m)')
options('li_simple', 'thick', [1,1])
options('li_simple', 'color', ['b','r'])                    
options('li_simple', 'ylog', True)

options(['li_simple','le_simple'], 'legend_size', 22)

tplot([B_name,'le_simple','li_simple'])
# %%

# %%
