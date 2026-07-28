#%%
import numpy as np
import matplotlib.pyplot as plt

#Need to update with WIND ion results

#And double-check those dB terms

#New Conjunction Stats - Fast Wind
x = np.array([1,2,3])

def get_fastvars():
	narray = 1e6*np.array([npsp_fast,nsolo_fast,nwind_fast])
	n_stdarray = 1e6*np.array([nstdpsp_fast,nstdsolo_fast,nstdwind_fast])
	Tarray = 1.602e-19*np.array([Tpsp_fast,Tsolo_fast,Twind_fast])
	T_stdarray = 1.602e-19*np.array([Tstdpsp_fast,Tstdsolo_fast,Tstdwind_fast])
	varray = 1e3*np.array([vpsp_fast,vsolo_fast,vwind_fast])
	v_stdarray = 1e3*np.array([vstdpsp_fast,vstdsolo_fast,vstdwind_fast])
	Barray = 1e-9*np.array([Bpsp_fast,Bsolo_fast,Bwind_fast])
	B_stdarray = 1e-9*np.array([Bstdpsp_fast,Bstdsolo_fast,Bstdwind_fast])
	dB_array = 1e-9*np.array([dBpsp_fast,dBsolo_fast,dBwind_fast])
	dB_stdarray = 1e-9*np.array([dBstdpsp_fast,dBstdsolo_fast,dBstdwind_fast])
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
	narray = 1e6*np.array([npsp_slow,nsolo_slow,nwind_slow])
	n_stdarray = 1e6*np.array([nstdpsp_slow,nstdsolo_slow,nstdwind_slow])
	Tarray = 1.602e-19*np.array([Tpsp_slow,Tsolo_slow,Twind_slow])
	T_stdarray = 1.602e-19*np.array([Tstdpsp_slow,Tstdsolo_slow,Tstdwind_slow])
	varray = 1e3*np.array([vpsp_slow,vsolo_slow,vwind_slow])
	v_stdarray = 1e3*np.array([vstdpsp_slow,vstdsolo_slow,vstdwind_slow])
	Barray = 1e-9*np.array([Bpsp_slow,Bsolo_slow,Bwind_slow])
	B_stdarray = 1e-9*np.array([Bstdpsp_slow,Bstdsolo_slow,Bstdwind_slow])
	dB_array = 1e-9*np.array([dBpsp_slow,dBsolo_slow,dBwind_slow])
	dB_stdarray = 1e-9*np.array([dBstdpsp_slow,dBstdsolo_slow,dBstdwind_slow])
	rarray = np.array([0.05/3,0.3/3,1/3])


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
# # # Magnetic Fluctuation Ratio
# fastratio_err = (dBfast/Bfast)*np.sqrt((dBfast_std/dBfast)**2 + (Bfast_std/Bfast)**2)
# slowsratio_err = (dBslow/Bslow)*np.sqrt((dBslow_std/dBslow)**2 + (Bslow_std/Bslow)**2)

# plt.errorbar(x,dBfast/Bfast,yerr=fastratio_err,fmt='o-',color = 'r',label = r'$|\delta B/B|_{fast}$',capsize=5)
# plt.errorbar(x,dBslow/Bslow,yerr=slowsratio_err,fmt='o-',color='b',label = r'$|\delta B/B|_{slow}$',capsize=5)
# plt.plot(x,(dBslow/Bslow+dBfast/Bfast)/2,color='g',marker='s',linestyle='dashed',label = r'$|\delta B/B|_{avg}$')
# plt.title('$\delta B/B$',fontsize=15)
# plt.axvline(x=1,linestyle='dashed',color='k',linewidth=0.2)
# plt.axvline(x=2,linestyle='dashed',color='k',linewidth=0.2)
# plt.axvline(x=3,linestyle='dashed',color='k',linewidth=0.2)
# plt.yticks(fontsize=12)
# plt.xticks(x,['PSP \n $\sim0.05 \  AU$','SOLO \n $\sim0.3 \  AU$','WIND \n $\sim1 \ AU$	'],fontsize=12)
# plt.legend(fontsize=15)



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

