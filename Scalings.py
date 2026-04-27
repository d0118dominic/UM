#%%
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


emframe = np.zeros([3,3,3])
emnframe = np.zeros_like(emframe)

# EM,E,M values for every run (mean of last 10% of values in timeseries)
emframe[0,0] = [5.0026188,2.9880660,2.0145528]
emframe[0,1] = [3.3744791,2.1916420,1.1828374]
emframe[0,2] = [2.3615732,1.6898993,0.67167377]

emframe[1,0] = [8.5524759,5.1672764,3.3851984]
emframe[1,1] = [5.7791672,3.7635241,2.0156434]
emframe[1,2] = [4.0311985,2.9133432,1.1178551]

emframe[2,0] = [14.772674,9.0728302,5.6998425]
emframe[2,1] = [9.9524193,6.6284862,3.3239348]
emframe[2,2] = [7.0424581,5.1953168,1.871413]

# Create version of em array normalized to the reference run (N0,T0)
emnframe = emframe/emframe[0,0]

# Make len(3) arrays to idicate the rows (Dens) and Columns (Temps) in normalized units
Dens = np.array([1,2,4])
Temps = np.array([1,0.5,0.25])

# Make separate arrays for ue,um variables in a given row (samedens series) or column (sametemp series)
# Example: um_samedens0 array has all the um values in the 0th (N0) row 
um_samedens0 = np.array([emnframe[0,0,2],emnframe[0,1,2],emnframe[0,2,2]])
ue_samedens0 = np.array([emnframe[0,0,1],emnframe[0,1,1],emnframe[0,2,1]])
um_samedens1 = np.array([emnframe[1,0,2],emnframe[1,1,2],emnframe[1,2,2]])
ue_samedens1 = np.array([emnframe[1,0,1],emnframe[1,1,1],emnframe[1,2,1]])
um_samedens2 = np.array([emnframe[2,0,2],emnframe[2,1,2],emnframe[2,2,2]])
ue_samedens2 = np.array([emnframe[2,0,1],emnframe[2,1,1],emnframe[2,2,1]])

um_sametemp0 = np.array([emnframe[0,0,2],emnframe[1,0,2],emnframe[2,0,2]])
ue_sametemp0 = np.array([emnframe[0,0,1],emnframe[1,0,1],emnframe[2,0,1]])
um_sametemp1 = np.array([emnframe[0,1,2],emnframe[1,1,2],emnframe[2,1,2]])
ue_sametemp1 = np.array([emnframe[0,1,1],emnframe[1,1,1],emnframe[2,1,1]])
um_sametemp2 = np.array([emnframe[0,2,2],emnframe[1,2,2],emnframe[2,2,2]])
ue_sametemp2 = np.array([emnframe[0,2,1],emnframe[1,2,1],emnframe[2,2,1]])

# Put each set in a list, so averaging of values can be done later (probably convoluted, but it works)
uelist_samedens = [ue_samedens0,ue_samedens1,ue_samedens2]
umlist_samedens = [um_samedens0,um_samedens1,um_samedens2]
uelist_sametemp = [ue_sametemp0,ue_sametemp1,ue_sametemp2]
umlist_sametemp = [um_sametemp0,um_sametemp1,um_sametemp2]

#%%

# Power law function
def power_law(x, a, b):
    return a * x**b

# Function to fit to a power law
def get_fits(x,y):
    # Fit  data
    params, covariance = curve_fit(power_law, x, y)
    a_fit, exponent = params
    #print(f"Rounded Exponent: {round(exponent,1):.4f}")

    # Generate smooth curve for plotting
    x_smooth = np.linspace(min(x), max(x), 100)
    y_fitted = power_law(x_smooth, a_fit, exponent)
    
    return x_smooth,y_fitted,exponent,a_fit

# Function to return exponents for all sets, and their averages
def get_exps(x,uelist,umlist):
    expelist = np.zeros(3)
    expmlist = np.zeros(3)

    for i in range(3):
        a,b,expe,c = get_fits(x,uelist[i])
        a,b,expm,c = get_fits(x,umlist[i])
        expelist[i] = expe
        expmlist[i] = expm

    print(f"UE Exponents: {expelist}")
    print(f"Mean: {np.mean(expelist)}")
    print(f"Rounded: {round(np.mean(expelist),2)}")
    print(" ")
    print(f"UM Exponents: {expmlist}")
    print(f"Mean: {np.mean(expmlist)}")
    print(f"Rounded: {round(np.mean(expmlist),2)}")
    return round(np.mean(expelist),2), round(np.mean(expmlist),2)

# Summarise the exponents 
print("------Same Density, Temperature Variation----------- ")
Texp_E, Texp_M = get_exps(Temps,uelist_samedens,umlist_samedens)
print("----------------------------------------------------")
print("   ")
print("------Same Temp, Density Variation----------- ")
Nexp_E, Nexp_M = get_exps(Dens,uelist_sametemp,umlist_sametemp)
print("----------------------------------------------------")
#%%

# Get fit parameters for ue,um,ue/um for a chosen row and a chosen column (for simplicity, using (0,0) )
xsme_samedens, ysme_samedens, expe_samedens,fite_samedens = get_fits(Temps,ue_samedens0)
xsmm_samedens, ysmm_samedens, expm_samedens,fitm_samedens = get_fits(Temps,um_samedens0)
xsmr_samedens, ysmr_samedens, expr_samedens,fitr_samedens = get_fits(Temps,ue_samedens0/um_samedens0)

xsme_sametemp, ysme_sametemp, expe_sametemp,fite_sametemp = get_fits(Dens,ue_sametemp0)
xsmm_sametemp, ysmm_sametemp, expm_sametemp,fitm_sametemp = get_fits(Dens,um_sametemp0)
xsmr_sametemp, ysmr_sametemp, expr_sametemp,fitr_sametemp = get_fits(Dens,ue_sametemp0/um_sametemp0)

#%%


# Plot variation of um,ue,ue/um with Temperature (density invariant)
# Choosing row 0 in this case

font=15
ticksize=15
lw = 2

# Plot actual data 
plt.figure(figsize=[8,8])
plt.scatter(Temps,um_samedens0, marker='o',color='b',label="$U_M$")
plt.scatter(Temps,ue_samedens0,marker='o',color='r',label = "$U_E$")
plt.scatter(Temps,ue_samedens0/um_samedens0,marker='o',color='g',label="$U_E/U_M$")
plt.plot(Temps,Temps,color='k',label='$U_{th}$',linewidth=lw)
plt.grid(axis='both')

# Plot fits
plt.plot(xsme_samedens, ysme_samedens, '--',color ='r', label=rf"Fit: $T^{{{expe_samedens:.2f}}}$",linewidth=lw)
plt.plot(xsmm_samedens, ysmm_samedens, '--',color='b', label=rf"Fit: $T^{{{expm_samedens:.2f}}}$",linewidth=lw)
plt.plot(xsmr_samedens, ysmr_samedens, '--',color='g', label=rf"Fit: $T^{{{expr_samedens:.2f}}}$",linewidth=lw)
plt.title("Runs of Same Density",fontsize=font+5)
plt.xlabel("$T/T_0$",fontsize=font)
plt.ylabel("Normalized Value",fontsize=font)
plt.legend(ncol=2,fontsize=font)
ticks = [0,0.25,0.5,0.75,1,1.25]
plt.xticks(ticks,size=ticksize)
plt.yticks(size=ticksize)
plt.xlim(0,1.25)
plt.ylim(0,2)

#%%

# Plot variation of um,ue,ue/um with Density (temperature invariant)
# Choosing column 0 in this case

font=15
ticksize=15
lw=2

# Plot actual data 
plt.figure(figsize=[8,8])
plt.scatter(Dens,um_sametemp0, color='b',marker='o', label="$U_M$")
plt.scatter(Dens,ue_sametemp0,color='r',marker='o',label = "$U_E$")
plt.scatter(Dens,ue_sametemp0/um_sametemp0,color='g',marker='o',label="$U_E /U_M $")
plt.plot(Dens,Dens,color='k',label = '$U_{th}$',linewidth=lw)
plt.grid(axis='both')

# Plot fits
plt.plot(xsme_sametemp, ysme_sametemp, '--',color='r', label=rf"Fit: $N^{{{expe_sametemp:.2f}}}$",linewidth=lw)
plt.plot(xsmm_sametemp, ysmm_sametemp, '--',color='b', label=rf"Fit: $N^{{{expm_sametemp:.2f}}}$",linewidth=lw)
plt.plot(xsmr_sametemp, ysmr_sametemp, '--',color='g', label=rf"Fit: $N^{{{expr_sametemp:.2f}}}$",linewidth=lw)
plt.title("Runs of Same Temperature",fontsize=font+5)
# plt.grid()
plt.xlabel("$N/N_0$",fontsize=font)
plt.ylabel("Normalized Value",fontsize=font)
plt.legend(ncol=2,fontsize=font)
ticks = [1,2,3,4,5]
plt.xticks(ticks,size=ticksize)
plt.yticks(size=ticksize)
plt.xlim(0,5)
plt.ylim(0,5)

#%%


# For manual adjustment of exponents if wanted
# Nexp_E = 0.8
# Texp_E = 0.4
# Nexp_M = 0.8
# Texp_M = 0.8

def get_energies(V,T):

    Nexp_EM = (Nexp_E+Nexp_M)/2
    Texp_EM = (Texp_E+Texp_M)/2
    uem = (V**(-Nexp_EM))*(T**(Texp_EM)) 
    um = (V**(-Nexp_M))*(T**(Texp_M)) 
    ue = (V**(-Nexp_E))*(T**(Texp_E))  
    up = T/V

    uem_tot = uem*V 
    um_tot = um*V 
    ue_tot = ue*V 
    up_tot = up*V 
    return uem,um,ue,up,uem_tot,um_tot,ue_tot,up_tot



def solve_energies(V):

    Nexp_EM = (Nexp_E+Nexp_M)/2
    Texp_EM = (Texp_E+Texp_M)/2


    #V_chor, T_itherm, T_ibar
    V_chor = 1
    T_itherm = 1
    T_ibar = V

    #T_ad
    T_ad = V**(-2/3)

    #T_imag, T_ielec, T_iem
    T_imag = V**(Nexp_M/Texp_M)
    T_ielec = V**(Nexp_E/Texp_E)
    T_iem= V**(Nexp_EM/Texp_EM)

    return V_chor, T_ad, T_itherm, T_ibar, T_imag, T_ielec, T_iem





#Need to change these conditions based on experimental values

V = np.linspace(0,10,100000)

V_chor,T_ad,T_itherm,T_ibar,T_imag,T_ielec,T_iem = solve_energies(V)

# %%

import numpy as np
import matplotlib.pyplot as plt


uem_ad,um_ad,ue_ad,up_ad,uem_ad_tot,um_ad_tot,ue_ad_tot,up_ad_tot = get_energies(V,T_ad)
uem_itherm,um_itherm,ue_itherm,up_itherm,uem_itherm_tot,um_itherm_tot,ue_itherm_tot,up_itherm_tot = get_energies(V,T_itherm)
uem_ibar,um_ibar,ue_ibar,up_ibar,uem_ibar_tot,um_ibar_tot,ue_ibar_tot,up_ibar_tot = get_energies(V,T_ibar)
uem_imag,um_imag,ue_imag,up_imag,uem_imag_tot,um_imag_tot,ue_imag_tot,up_imag_tot = get_energies(V,T_imag)
uem_ielec,um_ielec,ue_ielec,up_ielec,uem_ielec_tot,um_ielec_tot,ue_ielec_tot,up_ielec_tot = get_energies(V,T_ielec)
uem_iem,um_iem,ue_iem,up_iem,uem_iem_tot,um_iem_tot,ue_iem_tot,up_iem_tot = get_energies(V,T_iem)

#%%

font = 10
fig, axs = plt.subplots(2,3, figsize=(20,10))
# axs[0,0].plot(V, uem_itherm, label = r"$P_EM$")
axs[0,0].plot(V, ue_itherm, label = r"$P_E$",color='r')
axs[0,0].plot(V, um_itherm, label = r"$P_M$",color='b')
axs[0,0].plot(V,up_itherm,label=r"$P_{th}$",color='k')
axs[0,0].title.set_text('IsoThermal')

# axs[0,1].plot(V, uem_ibar, label = r"$P_{EM}$")
axs[0,1].plot(V, ue_ibar, label = r"$P_E$",color='r')
axs[0,1].plot(V, um_ibar, label = r"$P_M$",color='b')
axs[0,1].plot(V,up_ibar,label=r"$P_{th}$",color='k')
axs[0,1].title.set_text('IsoBaric')

# axs[0,2].plot(V, uem_ad, label = "PEM")
axs[0,2].plot(V, ue_ad, label = r"$P_E$",color='r')
axs[0,2].plot(V, um_ad, label = r"$P_M$",color='b')
axs[0,2].plot(V, up_ad, label=r"$P_{th}$",color='k')
axs[0,2].title.set_text('Adiabatic')

# axs[1,0].plot(V, uem_ielec, label = "PEM")
axs[1,0].plot(V, ue_ielec, label = r"$P_E$",color='r')
axs[1,0].plot(V, um_ielec, label = r"$P_M$",color='b')
axs[1,0].plot(V,up_ielec,label=r"$P_{th}$",color='k')
axs[1,0].title.set_text('IsoElectric')


# axs[1,1].plot(V, uem_imag, label = "PEM")
axs[1,1].plot(V, ue_imag, label = r"$P_E$",color='r')
axs[1,1].plot(V, um_imag, label = r"$P_M$",color='b')
axs[1,1].plot(V, up_imag, label=r"$P_{th}$",color='k')
axs[1,1].title.set_text('IsoMagnetic')


# axs[1,2].plot(V, uem_iem, label = "PEM")
axs[1,2].plot(V, ue_iem, label = r"$P_E$",color='r')
axs[1,2].plot(V, um_iem, label = r"$P_M$",color='b')
axs[1,2].plot(V, up_iem, label=r"$P_{th}$",color='k')
axs[1,2].title.set_text('IsoElectroMagnetic')

for i in range(2):
    for j in range(3):
        axs[i,j].grid()
        axs[i,j].legend(loc='best', prop={'size': 14})
        axs[i,j].set_xlabel(r'$V/V_0$',fontsize=15)  
        axs[i,j].set_ylabel(r'$P/P_0$',fontsize=15)  
        axs[i,j].set_xlim(0,2)
        axs[i,j].set_ylim(0,2)


plt.subplots_adjust(hspace=0.3)  # Increase this value for more space

#%%
fig, axs = plt.subplots(2,2, figsize=(15,15))
axs[0,0].plot(V, uem_ad, label = "Adiabatic",color='c')
axs[0,0].plot(V, uem_ibar, label = "Isobaric",color='k')
axs[0,0].plot(V, uem_itherm, label = "Isothermal",color='g')
axs[0,0].plot(V,uem_ielec,label="Isoelectric",color='r')
axs[0,0].plot(V,uem_imag,label="Isomagnetic",color='b')
axs[0,0].plot(V,uem_iem,label="Isoelectromagnetic",color='m')
axs[0,0].set_title(r'$P_{EM}$ Curves',fontsize=15)

axs[0,1].plot(V, up_ad, label = "Adiabatic",color='c')
axs[0,1].plot(V, up_ibar, label = "Isobaric",color='k')
axs[0,1].plot(V, up_itherm, label = "Isothermal",color='g')
axs[0,1].plot(V,up_ielec,label="Isoelectric",color='r')
axs[0,1].plot(V,up_imag,label="Isomagnetic",color='b')
axs[0,1].plot(V,up_iem,label="Isoelectromagnetic",color='m')
axs[0,1].set_title(r'$P_{th}$ Curves',fontsize=15)

axs[1,0].plot(V, ue_ad/um_ad, label = "Adiabatic",color='c')
axs[1,0].plot(V, ue_ibar/um_ibar, label = "Isobaric",color='k')
axs[1,0].plot(V, ue_itherm/um_itherm, label = "Isothermal",color='g')
axs[1,0].plot(V,ue_ielec/um_ielec,label="Isoelectric",color='r')
axs[1,0].plot(V,ue_imag/um_imag,label="Isomagnetic",color='b')
axs[1,0].plot(V,ue_iem/um_iem,label="Isoelectromagnetic",color='m')
axs[1,0].set_title(r'$P_E/P_M$ Curves',fontsize=15)

axs[1,1].plot(V, up_ad/uem_ad, label = "Adiabatic",color='c')
axs[1,1].plot(V, up_ibar/uem_ibar, label = "Isobaric",color='k')
axs[1,1].plot(V, up_itherm/uem_itherm, label = "Isothermal",color='g')
axs[1,1].plot(V,up_ielec/uem_ielec,label="Isoelectric",color='r')
axs[1,1].plot(V,up_imag/uem_imag,label="Isomagnetic",color='b')
axs[1,1].plot(V,up_iem/uem_iem,label="Isoelectromagnetic",color='m')
axs[1,1].set_title(r'$P_{th}/P_{EM}$ Curves',fontsize=15)


# axs[1,0].plot(V, uem_ad_tot, label = "adiabatic")
# axs[1,0].plot(V, uem_ibar_tot, label = "isobaric")
# axs[1,0].plot(V, uem_itherm_tot, label = "isothermal")
# axs[1,0].plot(V,uem_ielec_tot,label="isoelectric")
# axs[1,0].plot(V,uem_imag_tot,label="isomagnetic")
# # axs[1,0].plot(V,uem_iem_tot,label="isoelectromagnetic")
# axs[1,0].title.set_text('Uem Curves')

# axs[1,1].plot(V, up_ad_tot, label = "adiabatic")
# axs[1,1].plot(V, up_ibar_tot, label = "isobaric")
# axs[1,1].plot(V, up_itherm_tot, label = "isothermal")
# axs[1,1].plot(V,up_ielec_tot,label="isoelectric")
# axs[1,1].plot(V,up_imag_tot,label="isomagnetic")
# # axs[1,1].plot(V,up_iem_tot,label="isoelectromagnetic")
# axs[1,1].title.set_text('Up Curves')

# axs[1,2].plot(V, ue_ad_tot/um_ad_tot, label = "adiabatic")

# axs[1,2].plot(V, ue_ibar_tot/um_ibar_tot, label = "isobaric")
# axs[1,2].plot(V, ue_itherm_tot/um_itherm_tot, label = "isothermal")
# axs[1,2].plot(V,ue_ielec_tot/um_ielec_tot,label="isoelectric")
# axs[1,2].plot(V,ue_imag_tot/um_imag_tot,label="isomagnetic")
# # axs[1,2].plot(V,ue_iem_tot/um_iem_tot,label="isoelectromagnetic")
# axs[1,2].title.set_text('Ue/Um Curves')

# axs[1,3].plot(V, up_ad_tot/uem_ad_tot, label = "adiabatic")
# axs[1,3].plot(V, up_ibar_tot/uem_ibar_tot, label = "isobaric")
# axs[1,3].plot(V, up_itherm_tot/uem_itherm_tot, label = "isothermal")
# axs[1,3].plot(V,up_ielec_tot/uem_ielec_tot,label="isoelectric")
# axs[1,3].plot(V,up_imag_tot/uem_imag_tot,label="isomagnetic")

# axs[1,3].plot(V,up_iem_tot/uem_iem_tot,label="isoelectromagnetic")
# axs[1,3].title.set_text('U_p/Uem Curves')

for i in range(2):
    for j in range(2):
        axs[i,j].grid()
        axs[i,j].legend(loc='best', prop={'size': 10})
        axs[i,j].set_xlabel(r'$V/V_0$',fontsize=15)  
        axs[i,j].set_xlim(0,2)
        axs[i,j].set_ylim(0,2)

#axs[0,0].set_ylabel(r'$P_{EM}/P_{EM_0}$',fontsize=15)
#axs[0,1].set_ylabel(r'$P_{th}/P_{th_0}$',fontsize=15)
#axs[1,0].set_ylabel(r'$P/P_0$',fontsize=15)


plt.subplots_adjust(hspace=0.3)  # Increase this value for more space

# %%
# PV Diagrams - Carnot Cycle
fig, axs = plt.subplots(2,1, figsize=(10,10))


axs[0].plot(V, up_itherm+0.3, label = "Itherm",color='r')
axs[0].plot(V, up_itherm+1, label = "Itherm",color='r')
axs[0].plot(V+0.2, up_ad, label = "Itherm",color='g')
axs[0].plot(V+0.5, up_ad, label = "Itherm",color='g')



axs[1].plot(V, ue_iem/um_iem,color='r')
axs[1].plot(V, uem_imag/uem_ielec,color='b')


axs[0].set_ylim(0,3)
axs[1].set_ylim(0,3)


axs[0].set_xlim(0,2)
axs[1].set_xlim(0,2)
axs[1].axhline(1,color='k',linestyle='--')
# plt.plot(V, um_itherm_tot, label = "UM")
# plt.plot(V,up_itherm_tot,label="Ut")




# %%



