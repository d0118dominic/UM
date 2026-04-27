
#%%
import matplotlib


import matplotlib.pyplot as plt
import numpy as np


matplotlib.rcParams['figure.figsize'] = (8, 8)

import scipy.constants as constants
mu0 = constants.mu_0
eps0 = constants.epsilon_0
mi = constants.proton_mass
me = constants.electron_mass

c = constants.c
e = constants.e
kB = constants.k
#%%

def get_critscales(E,B,m,n,q):
    tm = (1/q)*((2*B/(n*mu_0))**(1/3))*((m/E)**(2/3))        #Tanh
    lm = (1/q)*((B/(n*mu_0))**(2/3))*((m/(2*E))**(1/3))        #Tanh
    td = m/(q*B)
    ld = m*E/(2*q*B**2)
    return tm,td,lm,ld

def get_crittimescales(E,B,m,n,q):
    tm = (1/q)*((2*B/(n*mu_0))**(1/3))*((m/E)**(2/3))        #Tanh
    td = m/(q*B)
    return tm,td,lm,ld

def get_maxenflux(r, t, E, B, m, n, q, T):
    
    # existing terms
    Kdr = np.sqrt(2*(n**2)/m) * (q*E*r)**(3/2)
    Kdt = n*(q*E*t)**3 / (2*m**2)
    Hdr = (5/2) * n * kb * T * np.sqrt(2*q*E*r/m)   
    Hdt = (5/2) * n * kb * T * (q*E*t/m)            
    KD  = (n*m/2) * (E/B)**3
    HD  = (5/2) * n * kb * T * (E/B)        
    S   = E*B / mu_0

    # enthalpy flux terms

    return Kdr, Kdt, KD, S, HD, Hdr, Hdt


def colorplot(x,y,z):
    xgrid,ygrid=np.meshgrid(x,y)
    zgrid = z(xgrid,ygrid)

    fig,ax=plt.subplots()
    im = ax.pcolormesh(xgrid,ygrid,zgrid,cmap=cmap,shading='auto')

    cbar = fig.colorbar(im,ax=ax)
    cbar.set_label(str(z))

    ax.set_xlabel(str(x))
    ax.set_ylabel(str(y))
    return fig, ax

q = e
mu_0 = mu0
n = 1e6*10 #m^-3
m = mi 

E = np.logspace(-5,5,100) #V/m
B = np.logspace(-10,1,100) # T

# fig, ax = colorplot(E,B,get_critscales(E,B,m,n,q))
plt.show()



#%%

def get_critscales(E,B,m,n,q):
    mu0 = constants.mu_0
    q = constants.e
    m=mi
    tm = (1/q)*((2*B/(n*mu0))**(1/3))*((m/E)**(2/3))        #Tanh
    lm = (1/q)*((B/(n*mu0))**(2/3))*((m/(2*E))**(1/3))        #Tanh
    td = m/(q*B)
    ld = m*E/(2*q*B**2)
    return tm,td,lm,ld


envs = {
    "EDR":         (50,   20),
    "Magnetotail": (1,    10),
    "Solar Wind":  (5,    5),
    "Inner MS":    (10,   200),
    "Bow Shock":   (0.1,  20),
    # "Solar Corona":   (1e3*0.1,  100),
}
n=1e7
import matplotlib.colors as colors
E = np.logspace(-5,1,1000) #V/m
B = np.logspace(-10,-5,1000) # T

EE,BB = np.meshgrid(E,B)
tm,td,lm,ld = get_critscales(EE,BB,mi,n,q)
norm = colors.LogNorm(vmin=1e-2,vmax=1e2)
fig,ax=plt.subplots()
im=ax.pcolormesh(1e3*EE,1e9*BB,tm/td,cmap = 'seismic',norm=norm)



# for label, (E_pt, B_pt) in envs.items():
#     ax.scatter(E_pt, B_pt, s=80, zorder=5,
#                edgecolors='white', linewidths=0.8, color='k')
#     ax.annotate(label, xy=(E_pt, B_pt), xytext=(6, 4),
#                 textcoords='offset points', fontsize=8, color='white')


# E_path = np.array([1e-4, 1e-3, 1e-2, 1e-1,1e-1])  # V/m
# B_path = np.array([1e-9, 1e-8, 1e-7, 1e-6,1e-7])  # T
# ax.plot(1e3 * E_path, 1e9 * B_path, 
#         color='lime', linewidth=2, marker='o', markersize=5, zorder=5)
t = np.linspace(0, 1, 200)
E_path = 1e-4 * (1e2)*t   # exponential sweep in E
B_path = 1e-9 * (1e3)**t   # exponential sweep in B
ax.plot(1e3 * E_path, 1e9 * B_path, color='k', lw=2, zorder=5)


ax.set_xlabel('$|E| \ (mV/m)$',fontsize=15)
ax.set_ylabel('$|B| \  (nT)$',fontsize=15)
ax.set_xscale('log')
ax.set_yscale('log')




cbar = plt.colorbar(im)
cbar.set_label(r'$t_M/t_d$',fontsize=15)
plt.show()
#%%
# Current sheet half-width
lam = 1.0  # normalized units

# 1D coordinate across the sheet
y = np.linspace(-5*lam, 5*lam, 1000)

# Harris sheet profiles
B0 = 20e-9      # T (20 nT asymptotic lobe field)
n0 = 1e6        # m^-3 (1 cm^-3 sheet density)
nb = 0.1e6      # m^-3 background/lobe density
E0 = 1e-3       # V/m (peak reconnection E field)

B = B0 * np.tanh(y / lam)
n = n0 * (1/np.cosh(y/lam))**2 + nb
E = E0 * (1/np.cosh(y/lam))**2  # peaks at center, falls off like density

plt.plot(B)
plt.plot(1e-5*E)

#%%
tm_path, td_path, lm_path, ld_path = get_critscales(abs(E), abs(B), mi, n, q)

plt.plot(lm_path/ld_path)
# plt.plot(ld_path)
plt.ylim(0,2)

#%%



plt.plot(tm_path)
plt.plot(td_path)


#%%


import numpy as np
import matplotlib.pyplot as plt

# tau/tauD axis - log spaced
tau_ratio = np.logspace(-5, 5, 1000)

# theta in each limit
theta_kinetic = np.arctan((tau_ratio)**3)
theta_thermal = np.arctan(tau_ratio)

# convert to degrees
theta_kinetic_deg = np.degrees(theta_kinetic)
theta_thermal_deg = np.degrees(theta_thermal)

fig, ax = plt.subplots(figsize=(6, 5))

ax.plot(tau_ratio, theta_kinetic_deg, color='royalblue', lw=2.5, label=r'Kinetic limit')
ax.plot(tau_ratio, theta_thermal_deg, color='firebrick', lw=2.5, linestyle='-', label=r'Thermal limit')

# critical scale line
ax.axvline(x=1, color='k', lw=1.5, linestyle=':', alpha=0.7, label=r'$\tau = \tau_D$')
ax.axhline(y=45, color='gray', lw=1, linestyle=':', alpha=0.5)



# pi/2 and 0 labels on y axis
ax.set_yticks([0, 15, 30, 45, 60, 75, 90])
ax.set_yticklabels(['0', '15°', '30°', '45°', '60°', '75°', '90°'])

ax.set_xscale('log')
ax.set_xlim(1e-3, 1e3)
ax.set_ylim(0, 90)

ax.set_xlabel(r'$\tau / \tau_D$', fontsize=13)
ax.set_ylabel(r'$\theta$', fontsize=13)
ax.set_title(r'Angle Between Total Plasma Energy Transport and MTD transport', fontsize=10)

ax.legend(fontsize=10, loc='upper left')


plt.tight_layout()
plt.savefig('theta_vs_tau.png', dpi=150, bbox_inches='tight')
plt.show()



# ell/ellD axis - log spaced
ell_ratio = np.logspace(-5, 5, 1000)

# theta in each limit
theta_kinetic = np.arctan((ell_ratio)**(3/2))
theta_thermal = np.arctan((ell_ratio)**(1/2))

# convert to degrees
theta_kinetic_deg = np.degrees(theta_kinetic)
theta_thermal_deg = np.degrees(theta_thermal)

fig, ax = plt.subplots(figsize=(6, 5))

ax.plot(ell_ratio, theta_kinetic_deg, color='royalblue', lw=2.5, label=r'Kinetic limit')
ax.plot(ell_ratio, theta_thermal_deg, color='firebrick', lw=2.5, linestyle='-', label=r'Thermal limit')

# critical scale line
ax.axvline(x=1, color='k', lw=1.5, linestyle=':', alpha=0.7, label=r'$\mathscr{l} = \mathscr{l}_D$')
ax.axhline(y=45, color='gray', lw=1, linestyle=':', alpha=0.5)



# pi/2 and 0 labels on y axis
ax.set_yticks([0, 15, 30, 45, 60, 75, 90])
ax.set_yticklabels(['0', '15°', '30°', '45°', '60°', '75°', '90°'])

ax.set_xscale('log')
ax.set_xlim(1e-5, 1e5)
ax.set_ylim(0, 90)

ax.set_xlabel(r'$\mathscr{l} / \mathscr{l}_D$', fontsize=13)
ax.set_ylabel(r'$\theta$', fontsize=13)
ax.set_title(r'Angle between Total Plasma Energy Transport and MTD transport', fontsize=10)

ax.legend(fontsize=10, loc='upper left')


plt.tight_layout()
plt.savefig('theta_vs_tau.png', dpi=150, bbox_inches='tight')
plt.show()
# %%


# ell/ellD axis - log spaced
ell_ratio = np.logspace(-5, 5, 1000)

jdeods = (ell_ratio)**(3/2)

# convert to degrees
fig, ax = plt.subplots(figsize=(6, 5))

ax.plot(ell_ratio, jdeods, color='royalblue', lw=2.5, label = r'$\frac{\delta J \delta E}{\delta S/\ell}$')
ax.plot(ell_ratio, ell_ratio/ell_ratio, color='k', lw=2.5,label = r'$\frac{\delta J \delta E}{\delta K_E/\ell}$')

# critical scale line
ax.axvline(x=1, color='k', lw=1.5, linestyle=':', alpha=0.7)
ax.axhline(y=1, color='k', lw=1.5, linestyle=':', alpha=0.5)



# pi/2 and 0 labels on y axis
ax.set_yticks([0, 0.5, 1, 1.5, 2])
ax.set_yticklabels(['0', '0.5', '1', '1.5', '2'])

ax.set_xscale('log')
ax.set_xlim(1e-2, 1e2)
ax.set_ylim(0, 2)

ax.set_xlabel(r'$\mathscr{l} / \mathscr{l}_{KS}$', fontsize=13)
# ax.set_ylabel(r'$\frac{\delta J \delta E}{\delta S/\ell}$', fontsize=17)
# ax.set_title(r'Angle between Total Plasma Energy Transport and MTD transport', fontsize=10)

ax.legend(fontsize=17, loc='upper left')


plt.tight_layout()
# plt.savefig('theta_vs_tau.png', dpi=150, bbox_inches='tight')
plt.show()


# ell/ellD axis - log spaced
ell_ratio = np.logspace(-5, 5, 1000)

jdeodhe = (ell_ratio)
jdeodpdv = (3/2)*(ell_ratio)

# convert to degrees
fig, ax = plt.subplots(figsize=(6, 5))

ax.plot(ell_ratio, jdeodhe, color='royalblue', lw=2.5,label = r'$\frac{\delta J \delta E}{\delta H_E/\ell}$')
ax.plot(ell_ratio, jdeodpdv, color='firebrick', lw=2.5,label = r'$\frac{\delta J \delta E}{p\delta v/\ell}$')
ax.plot(ell_ratio, ell_ratio/ell_ratio, color='k', lw=2.5,label = r'$\frac{\delta J \delta E}{\delta K_E/\ell}$')

# critical scale line
ax.axvline(x=1, color='k', lw=1.5, linestyle=':', alpha=0.5)
ax.axhline(y=1, color='k', lw=1.5, linestyle=':', alpha=0.5)



# pi/2 and 0 labels on y axis
ax.set_yticks([0, 0.5, 1, 1.5, 2])
ax.set_yticklabels(['0', '0.5', '1', '1.5', '2'])

ax.set_xscale('log')
ax.set_xlim(1e-2, 1e2)
ax.set_ylim(0, 2)

ax.set_xlabel(r'$\mathscr{l} / \mathscr{l}_{KH}$', fontsize=13)
# ax.set_ylabel(r'$\frac{\delta J \delta E}{\delta H_E/\ell}$', fontsize=17)
# ax.set_title(r'Angle between Total Plasma Energy Transport and MTD transport', fontsize=10)

ax.legend(fontsize=17, loc='upper left')


plt.tight_layout()
# plt.savefig('theta_vs_tau.png', dpi=150, bbox_inches='tight')
plt.show()
# %%

# tau/tauD axis - log spaced
tau_ratio = np.logspace(-5, 5, 1000)
# tau_ratio = np.linspace(0.1, 10, 1000)

# theta in each limit
jdeoduk = (2*tau_ratio/tau_ratio)
jdeodut = (2*(tau_ratio**2))
jdeodup = (tau_ratio**2)



fig, ax = plt.subplots(figsize=(6, 5))

ax.plot(tau_ratio, jdeoduk, color='royalblue', lw=2.5, label=r'$\frac{\delta J \delta E}{\delta u_K/\tau}$')
ax.plot(tau_ratio, jdeodut, color='firebrick', lw=2.5, label=r'$\frac{\delta J \delta E}{\delta u_T/\tau}$')
ax.plot(tau_ratio, jdeodup, color='orange', lw=2.5, label=r'$\frac{\delta J \delta E}{\delta u_P/\tau}$')

# critical scale line
ax.axvline(x=1, color='k', lw=1.5, linestyle=':', alpha=0.7)
ax.axhline(y=45, color='gray', lw=1, linestyle=':', alpha=0.5)



# pi/2 and 0 labels on y axis
ax.set_yticks([0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4])
ax.set_yticklabels(['0', '0.5', '1', '1.5', '2', '2.5', '3', '3.5', '4'])

ax.set_xscale('log')
ax.set_xlim(1e-1, 1e1)
ax.set_ylim(0, 4)

ax.set_xlabel(r'$\tau / \tau_{KH}$', fontsize=13)
ax.axhline(y=1, color='k', lw=1, linestyle=':', alpha=0.5)
# ax.set_ylabel(r'$\theta$', fontsize=13)
# ax.set_title(r'Angle Between Total Plasma Energy Transport and MTD transport', fontsize=10)

ax.legend(fontsize=17, loc='upper left')


plt.tight_layout()
plt.savefig('theta_vs_tau.png', dpi=150, bbox_inches='tight')
plt.show()

# %%

#__________________________

#%%


# %%
