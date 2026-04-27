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



def dataplot2():
	bins=100
	fig,ax = plt.subplots(1,3,figsize=(18,6))
	ax[0].hist(np.log10(allmachs),bins=bins,range=[-1.5,1.5],log=False)	
	ax[0].set_xlim(-0.5,0.5)
	ax[0].set_xlabel('Log10(Ma)')
	ax[0].set_ylabel('Counts')
	ax[0].axvline(x=0,color='k')
	ax[1].hist(allpositions,bins=bins,range=[0,50],log=False)
	ax[1].set_xlim(10,55)
	ax[1].set_xlabel('R (solar Radii)')
	ax[1].set_ylabel('Counts')

	ax[2].hist(allangles,bins=bins,range=[0,180],log=True)
	ax[2].set_xlim(0,180)
	ax[2].set_xlabel('Deflection Angle (degrees)')
	ax[2].set_ylabel('Counts')
	return()

def dataplot():
	bins=50
	font=12
	fig,ax = plt.subplots(1,2,figsize=(12,5))
	plt.subplots_adjust(wspace=0.3)  

	ax[0].hist(np.log10(allmachs),bins=bins,range=[-1.5,1.5],log=False)
	ax[0].set_xlim(-0.6,0.6)
	ax[0].set_xlabel(r'$log_{10}(Ma)$',fontsize=font)
	ax[0].set_ylabel('Counts',fontsize=font)
	ax[0].axvline(x=0,color='k')
	ax[1].hist(allpositions,bins=bins,range=[0,50],log=True)
	ax[1].set_xlim(10,45)
	ax[1].set_ylim(1e3,1e5)
	ax[1].set_xlabel(r'$R/R_{sun}$',fontsize=font)
	ax[1].set_ylabel('Counts',fontsize=font)

	# ax[2].hist(allangles,bins=bins,range=[0,180],log=True)
	# ax[2].set_xlim(0,180)
	# ax[2].set_xlabel('Deflection Angle (degrees)')
	# ax[2].set_ylabel('Counts')
	return()

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
def filter_angle(angle,floor,ceiling):
	angle_reduced = angle
	for i in range(len(timeax)): 
		if (angle[i] < floor or angle[i]>ceiling): angle_reduced[i] = np.nan
	return angle_reduced

def filter_deflections(var,floor):
	filtered_var = np.zeros_like(var)
	for i in range(len(var)):
		if (alldv_norm_mag[i] < floor or alldB_norm_mag[i] < floor):
			filtered_var[i] = np.nan
		else:
			filtered_var[i] = var[i]
	return filtered_var

def FOV_filter():
	upper = phis[0,1]
	lower = phis[0,-2]

	peak_phi = phis[np.arange(len(timeax)),np.argmax(ephi,axis=0)]
	mask = (peak_phi>=upper) | (peak_phi<=lower)
	Ti[mask],Tpar[mask],Tperp[mask] = np.nan,np.nan,np.nan
	return 

##### Untested ####
def filter_winds(var,sigma_threshold,v_threshold):
	slowAlf = np.zeros_like(var)
	fastAlf = np.zeros_like(var)
	nonAlf = np.zeros_like(var)
	for i in range(len(var)):
		if (abs(allcrosshelicity[i]) > sigma_threshold and allvmags[i] > v_threshold):
			fastAlf[i] = var[i] 
			slowAlf[i] = np.nan 
			nonAlf[i] = np.nan
		elif (abs(allcrosshelicity[i]) > sigma_threshold and allvmags[i] < v_threshold):
			fastAlf[i] = np.nan
			slowAlf[i] = var[i] 
			nonAlf[i] = np.nan
		elif (abs(allcrosshelicity[i]) < sigma_threshold and allvmags[i] < v_threshold):
			nonAlf[i] = var[i]
			fastAlf[i] = np.nan
			slowAlf[i] = np.nan 
		else:
			nonAlf[i] = np.nan
			fastAlf[i] = np.nan
			slowAlf[i] = np.nan
	return slowAlf,fastAlf,nonAlf

def filter_machs(var):
	sub = np.zeros_like(var)
	sup = np.zeros_like(var)
	for i in range(len(var)):
		if (allmachs[i]<1):
			sub[i] = var[i] 
			sup[i] = np.nan 
		elif (allmachs[i]>1):
			sub[i] = np.nan 
			sup[i] = var[i] 
		else:
			sub[i] = np.nan 
			sup[i] = np.nan 
	return sub,sup
# Velocities

def get_delta(vec,vec_mean):
	dvec = vec - vec_mean
	dvec_norm = dvec/np.linalg.norm(vec_mean)
	return dvec,dvec_norm
def get_deltascalar(var,var_mean):
	dvar = var - var_mean
	dvar_norm = dvar/var_mean
	return dvar,dvar_norm

def get_entropy(T,n):
	S = T / (n**(2/3))
	return S
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



#### IN PROGRESS
def get_vperps(dv,v,B):
	v_hat = v/np.linalg.norm(v) 
	B_hat = B/np.linalg.norm(B) # coord 1 (along B)
	vxb_hat = np.cross(v_hat,B_hat) # coord 2 (along vxB)
	lastcoord = np.cross(B_hat,vxb_hat) # coord 3 (perp to B and vxB)


	dv_par = np.dot(dv,B_hat) # dv component parallel to B
	dv_perp1 = np.dot(dv,vxb_hat) # dv component  parallel to vxB, perp to B
	dv_perp2 = np.dot(dv,lastcoord) #dv perp to B and vxB


	dv_perp = (np.linalg.norm(dv)**2 - dv_par**2)**0.5
	dv_perp1 = dv_perp

	return dv_par,dv_perp1,dv_perp2
####################################






# Scalar Pressures
def get_pm(B):
	pm = 0.5*(mu0**-1)*np.linalg.norm(B)**2
	return pm
def get_pth(n,T):
	pth = n*T
	return pth
# Normalized Parameters
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
# NEED TESTING.  Having pyspedas path issues...##
def get_parperps(n,T,B):  #B ant T tensor coord systems need to match for this
	trace = T[0] + T[1] + T[2]
	term1 = (T[0]*B[0]**2 + T[1]*B[1]**2 + T[2]*B[2]**2)/(np.linalg.norm(B)**2)
	term2 = 2*(T[3]*B[0]*B[1] + T[4]*B[0]*B[2] + T[5]*B[1]*B[2])/(np.linalg.norm(B)**2)
	Tpar = term1+term2
	Tperp=0.5*(trace-Tpar)
	Ppar = n*Tpar
	Pperp = n*Tperp
	return Tpar,Tperp,Ppar,Pperp
def get_par(v1,v2):   # Get magnitude of parallel component 
	par = np.dot(v1,v2)/np.linalg.norm(v2)
	return abs(par)
def get_voltagepairs():
	l_eff = 3.5
	pairxy = reform(get_data('psp_fld_l2_dfb_wf_dVdc_sc'))
	return pairxy/l_eff
def get_vxb(v,B):
	vxb = np.cross(v,B)
	return vxb
def get_K(m,n,v):
	K = 0.5*m*n*v**3
	return K
def get_H(v,P):
	traceP = P[0]+P[1]+P[2]
	vdPx = v[0]*(P[0]+P[3]+P[4])
	vdPy = v[1]*(P[1]+P[3]+P[5])
	vdPz = v[2]*(P[2]+P[4]+P[5])
	vdotP = np.array([vdPx,vdPy,vdPz])
	H = 0.5*v*traceP + vdotP
	return H
def get_S(E,B):
	S = np.cross(E,B)/mu0
	return S


def get_med(var,int):
	from scipy import ndimage
	return ndimage.median_filter(var, size=int, mode='nearest')
def get_mean(var,int):   #Take a timeseries and compute the mean (basically smooth outsmall flucs over some interval)
	box = np.ones(int)/int
	smoothed_var = np.convolve(var,box,mode='same')
	return smoothed_var

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

def get_vecmed(vec,int):   # Vector median
	vec1_mean = get_med(vec[:,0],minutes*meaninterval)
	vec2_mean = get_med(vec[:,1],minutes*meaninterval)
	vec3_mean = get_med(vec[:,2],minutes*meaninterval)
	vec_mean = np.zeros_like(vec)
	for i in range(len(vec_mean)): vec_mean[i] = np.array([vec1_mean[i],vec2_mean[i],vec3_mean[i]])
	return vec_mean

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


# %%
# Get Mag Data

sub_alfs =  [['2022-09-06/06:00','2022-09-06/16:00'], # 10 hrs Encounter 13
			 ['2022-09-06/18:00','2022-09-07/12:00'], # 18 hrs Encounter 13
			 ['2022-12-11/04:00','2022-12-11/16:00'], # 12 hrs Encounter 14
			 ['2023-03-16/12:00','2023-03-17/06:00'], # 18 hrs Encounter 15
			 ['2023-06-20/01:00','2023-06-21/01:00'], # 24 hrs Encounter 16
			 ['2023-09-27/06:00','2023-09-27/15:00'], # 9 hrs Encounter 17
			 ['2023-12-29/04:00','2023-12-29/14:00'], # 10 hrs Encounter 18
			 ['2024-03-29/06:00','2024-03-29/21:00'], # 15 hrs Encounter 19
			 ['2024-06-30/03:00','2024-06-30/18:00'], # 15 hrs Encounter 20
			 ['2024-09-28/11:00','2024-09-28/18:00'], # 7 hrs Encounter 21
 		     ['2024-12-24/00:00', '2024-12-25/00:00'], # 24 hours Encounter 22 
	         ['2025-03-22/06:00', '2025-03-23/12:00'] # 30 hrs Encounter 23
			 ]


sup_alfs=	[['2022-09-07/18:00','2022-09-08/18:00'], # 24 hrs Encounter 13
			 ['2022-12-10/01:00','2022-12-10/09:00'], # 8 hrs Encounter 14
			 ['2023-03-16/00:00','2023-03-16/11:00'], # 11 hrs Encounter 15
			 ['2023-06-24/00:00','2023-06-25/00:00'], # 24  hrs Encounter 16
			 ['2023-09-28/20:00','2023-09-30/09:00'], # 37 hrs Encounter 17
			 ['2023-12-25/14:00','2023-12-25/23:00'], # 9 hrs Encounter 18
			 ['2024-04-01/10:00','2024-04-02/01:00'], # 15 hrs Encounter 19
			 ['2024-07-01/20:00','2024-07-03/00:00'], # 28 hrs Encounter 20
			 ['2024-10-03/00:00','2024-10-03/12:00'], # 12 hrs Encounter 21
			 ['2025-03-21/00:00', '2025-03-22/03:00'], # 27 hrs Encounter 23
			 ]


near_alfs = [['2022-09-05/11:00','2022-09-05/17:00'], # 6 hrs Encounter 13
			 ['2023-06-23/00:00','2023-06-23/18:00'], # 18 hrs Encounter 16
			 ['2023-09-28/09:00','2023-09-28/18:00'], # 9 hrs Encounter 17
			 ['2023-12-28/00:00','2023-12-29/00:00'], # 24 hrs Encounter 18
			 ['2024-03-28/22:00','2024-03-29/05:00'], # 7 hrs Encounter 19
			 ['2024-06-29/00:00','2024-06-29/09:00'], # 9 hrs Encounter 20 
			 ['2024-10-01/15:00','2024-10-02/00:00'] # 9 hrs Encounter 21 
			]



#No CS crossings
sup_alfs_alt=[['2022-09-07/18:00','2022-09-08/18:00'], # 24 hrs Encounter 13
			 ['2022-12-10/01:00','2022-12-10/09:00'], # 8 hrs Encounter 14
			 ['2023-03-16/06:00','2023-03-16/11:00'], # 5 hrs Encounter 15
			 ['2023-06-24/09:00','2023-06-25/00:00'], # 15  hrs Encounter 16 ####
			 ['2023-09-28/20:00','2023-09-30/09:00'], # 37 hrs Encounter 17
			 ['2023-12-25/14:00','2023-12-25/23:00'], # 9 hrs Encounter 18
			 ['2024-04-01/10:00','2024-04-02/01:00'], # 15 hrs Encounter 19
			 ['2024-07-02/03:00','2024-07-03/00:00'], # 21 hrs Encounter 20
			 ['2024-10-03/00:00','2024-10-03/12:00'], # 12 hrs Encounter 21
			 ['2025-03-21/00:00', '2025-03-22/00:00'], # 27 hrs Encounter 23
			 ]


betaparlist = [['2024-03-25/04:00','2024-03-25/09:00'],
			   ['2024-07-04/22:00','2024-07-05/00:00'],
			   ['2024-09-25/05:00','2024-09-26/05:00'],
			   ['2024-12-19/00:00', '2024-12-20/00:00'],
			   ['2025-03-18/02:00', '2025-03-19/00:00'],
			   ['2025-03-17/00:00', '2025-03-17/12:00'],
			   ['2023-12-23/12:00','2023-12-24/00:00']]

encounter13 = [['2022-09-03/00:00','2022-09-03/01:00'],['2022-09-03/02:00','2022-09-03/03:00']]
# encounter13 = [['2022-09-01/00:00','2022-09-11/00:00'],['2022-09-03/00:00','2022-09-03/01:00']]

encounter24 = [['2025-06-17/00:00', '2025-06-21/00:00'], ['2025-06-21/00:00', '2025-06-21/01:00']]


# near_sups = [['2024-09-27/04:30','2024-09-27/06:30'], # 2 hr Encounter 21
# 			 ['2024-10-03/23:00','2024-10-04/00:00'], # 1 hr Encounter 21
# 			 ['2024-10-04/09:00','2024-10-04/10:00'], # 1 hr Encounter 21
# 			 ['2024-06-29/04:45','2024-06-29/06:15'], # 1.5 hr Encounter 20
# 			 ['2024-07-02/03:00','2024-07-02/08:00'], # 5 hr Encounter 20
# 			 ['2024-03-28/03:30','2024-03-28/06:00'], # 2.5 hr Encounter 19
# 			 ['2023-12-29/22:00','2023-12-30/03:00'], # 5 hrs Encounter 18
# 			 ['2023-09-29/00:00','2023-09-29/06:00'], # 6 hrs Encounter 17
# 			 ['2023-09-28/18:00','2023-09-28/21:00'], # 3 hrs Encounter 17
# 			 ['2025-03-21/10:00', '2025-03-21/14:00'], # 4 hrs Encounter 23
# 			 ]

# near_subs = [['2024-09-27/07:30','2024-09-27/8:30'], # 1 hr Encounter 21
# 			 ['2024-09-30/23:00','2024-10-01/01:00'],# 2 hr Encounter 21
# 			 ['2024-06-29/16:00','2024-06-29/21:00'],# 5 hr Encounter 20
# 			 ['2024-03-30/07:00','2024-03-30/11:00'], # 4 hr Encounter 19
# 			 ['2023-12-28/12:00','2023-12-28/15:00'], # 3 hrs Encounter 18
# 			 ]




# longints = [['2025-06-15/00:00','2025-06-23/00:00'],['2025-03-18/00:00','2025-03-24/00:00']]
# nears = near_subs+near_subs
# recent_perihelia = [['2024-09-29/00:00', '2024-10-01/12:00'], #E21
# 					['2024-06-29/00:00', '2024-07-01/12:00'], #E20
# 					['2024-03-29/00:00','2024-03-31/00:00'], #E19
# 					['2023-12-28/00:00','2023-12-30/00:00']] #E18


#Need to make a list of chosen Alfvenic & sub-alfvenic intervals
#Most people choose a handful of intervals by eye


# eventlist = subs

# Alfven crossings (<2 hr)
# short_crossings = [['2021-11-21/21:00','2021-11-21/22:00'],
# 				   ['2021-08-10/00:15', '2021-08-10/00:45'], # Encounter 9 (some sub-Alfvenic)
#                    ['2023-12-29/01:30','2023-12-29/03:00'], # E18
#                    ['2024-03-29/22:00','2024-03-29/23:30'],  # E19 
#                    ['2024-06-29/11:00', '2024-06-29/13:00']] # E20


#eventlist = sub_alfs+sup_alfs+near_alfs
# eventlist=near_alfs
# eventlist = sup_alfs+sub_alfs+near_alfs
# eventlist = sup_alfs_alt#+sub_alfs
eventlist = sup_alfs_alt+sub_alfs+near_alfs+betaparlist
# eventlist = sub_alfs
eventlist = encounter24
# eventlist = sup_alfs_alt[0:2]

# Sup-alfs with the large beta_par instability
# 7 (July 1-2 2024, 21:00-03:00) (some large-scale Br reversal)
# 9 (March 22 2025, 0-3:00) (might be the HCS)
# 2 ( March 16 2023 0-4:00) (more big Br reversalss)
# 3 (kinda) (June 24 2023 )

#%%

allBmags = np.array([])
allvmags = np.array([])
allSmags = np.array([])
allKmags = np.array([])
allHmags = np.array([])

allSpar = np.array([])
allKpar = np.array([])
allHpar = np.array([])

allKmags = np.array([])
allHmags = np.array([])
alldBmag = np.array([])
alldvmag = np.array([])
alldB = np.array([])
alldv = np.array([])
alldv_par_norm = np.array([])
alldv_perp_norm = np.array([])
alldv_par1_norm = np.array([])
alldv_perp1_norm = np.array([])
alldv_perp2_norm = np.array([])

alldB_par_norm = np.array([])
alldB_perp_norm = np.array([])
alldB_norm_mag = np.array([])
alldv_norm_mag = np.array([])
alldn_norm = np.array([])


allvr = np.array([])
allvt = np.array([])
allvn = np.array([])
allBr = np.array([])
allBt = np.array([])
allBn = np.array([])

allva = np.array([])
allva_inst = np.array([])
allbeta_par = np.array([])
allmagcomp = np.array([])

allangles=np.array([])
allpositions=np.array([])
allmachs=np.array([])
allSr=np.array([])
allKr=np.array([])
allHr=np.array([])
alldSmag=np.array([])
alldKmag=np.array([])
alldS_norm=np.array([])
alldK_norm=np.array([])
allmachs=np.array([])
allbetas=np.array([])
allvbalignment = np.array([])
allcrosshelicity = np.array([])
allresidenergy = np.array([])
allentropy = np.array([])

#Newvars
alln = np.array([])
allT = np.array([])
allTpar = np.array([])
allTperp = np.array([])
allTparperp = np.array([])

alldP_norm = np.array([])
alldPth_norm = np.array([])
alldPm_norm = np.array([])
alldSp_norm = np.array([])

for i in range(len(eventlist)):
	trange=eventlist[i]
	Bfld_vars = pyspedas.projects.psp.fields(trange=trange, level='l2', time_clip=True)
	swp_vars = pyspedas.projects.psp.spi(trange=trange,level='l3',get_support_data=True,time_clip=True)
	qtn_vars = pyspedas.projects.psp.fields(trange=trange,level='l3',datatype='sqtn_rfs_V1V2',time_clip=True)
	# alph_vars = pyspedas.projects.psp.spi(trange=trange,level='l3',datatype='sf0a_l3_mom',time_clip=True)
	# voltages_vars = pyspedas.projects.psp.fields(trange=trange, datatype='dfb_wf_dvdc', level='l2',time_clip=True)
	#On DC datatype: 'sqn_rfs_V1V2 has some kind of electron density & core temp, but looks weird
	# spc_vars = pyspedas.projects.psp.spc(trange=trange, datatype='l2', level='l2')
	
	# Reform all data to simple arrays and convert to SI units 
	B_name = 'psp_fld_l2_mag_RTN'
	vi_name = 'psp_spi_VEL_RTN_SUN'
	Bxyz_name = 'psp_spi_MAGF_INST'
	vxyz_name = 'psp_spi_VEL_INST'
	TiTensor_name = 'psp_spi_T_TENSOR_INST'
	Ti_name = 'psp_spi_TEMP'
	ni_name = 'psp_spi_DENS'
	ni_name = 'electron_density'
	phivals_name = 'psp_spi_PHI_VALS'
	ephi_name = 'psp_spi_EFLUX_VS_PHI'
	# voltages_name = 'psp_fld_l2_dfb_wf_dVdc_sc'
	position_name = 'psp_spi_SUN_DIST'

	interpvar_name = vi_name
	timeax = pytplot.get_data(interpvar_name).times
	meaninterval = mean_int(timeax,trange)

    # Handling troublesome qtn indices
	for name in [B_name, Bxyz_name, vxyz_name, vi_name, Ti_name, ni_name, TiTensor_name, position_name]:
		data = pytplot.get_data(name)
		if data is None:
			print(f"WARNING: no data found for {name}")
			continue
		times = data.times
		_, unique_idx = np.unique(times, return_index=True)
		if len(unique_idx) < len(times):
			print(f"Removing {len(times)-len(unique_idx)} duplicate timestamps from {name}")
			pytplot.store_data(name, data={'x': times[unique_idx], 'y': data.y[unique_idx]})

##%%
	tinterpol(B_name,interpvar_name,newname='B')
	tinterpol(Bxyz_name,interpvar_name,newname='Bxyz')
	tinterpol(vxyz_name,interpvar_name,newname='vxyz')
	tinterpol(vi_name,interpvar_name,newname='vi')
	tinterpol(Ti_name,interpvar_name,newname='Ti')
	tinterpol(ni_name,interpvar_name,newname='ni')
	tinterpol(TiTensor_name,interpvar_name,newname='TiTensor')
	tinterpol(phivals_name,interpvar_name,newname='phivals')
	tinterpol(ephi_name,interpvar_name,newname='ephi')
	# tinterpol(voltages_name,interpvar_name,newname='voltages')
	tinterpol(position_name,interpvar_name,newname='position')
	
	Bvecs = 1e-9*reform(pytplot.get_data('B'))
	Bxyz = 1e-9*reform(pytplot.get_data('Bxyz'))
	vxyz = 1e3*reform(pytplot.get_data('vxyz'))
	vivecs = 1e3*reform(pytplot.get_data('vi'))
	ni = 1e6*reform(pytplot.get_data('ni'))
	Ti = 1.602e-19*reform(pytplot.get_data('Ti'))
	TiTensor = 1.602e-19*reform(pytplot.get_data('TiTensor')) #Comes in xyz
	phis = reform(pytplot.get_data('phivals'))
	ephi = reform(pytplot.get_data('ephi')).T
	# PiTensor = ni*TiTensor
	# voltages = reform(get_data('voltages'))
	position = reform(get_data('position'))/695700 #Solar radii



	# Calculate B Magnitude & create Normalized Br/|B|
	Br,Bt,Bn = np.zeros_like(Bvecs[:,0]),np.zeros_like(Bvecs[:,0]),np.zeros_like(Bvecs[:,0])
	vr,vt,vn = np.zeros_like(Bvecs[:,0]),np.zeros_like(Bvecs[:,0]),np.zeros_like(Bvecs[:,0])

	B_mag = np.zeros_like(Bvecs[:,0])
	v_mag = np.zeros_like(Bvecs[:,0])

	Br_norm = np.zeros_like(Bvecs[:,0])
	Br_mean = np.zeros_like(Bvecs[:,0])
	vr_norm = np.zeros_like(Bvecs[:,0])
	E_conv = np.zeros_like(Bvecs)

	S = np.zeros_like(Bvecs)
	H = np.zeros_like(Bvecs)
	K = np.zeros_like(Bvecs)
	ExB = np.zeros_like(Bvecs)

	P_mag = np.zeros_like(ni)
	P_th = np.zeros_like(ni)
	entropy = np.zeros_like(ni)
	beta_par = np.zeros_like(ni)

	va = np.zeros_like(Br_norm)
	vth = np.zeros_like(Br_norm)
	vs = np.zeros_like(Br_norm)
	Tpar = np.zeros_like(Br_norm)
	Tperp = np.zeros_like(Br_norm)
	Ppar = np.zeros_like(Br_norm)
	Pperp = np.zeros_like(Br_norm)

	beta = np.zeros([len(timeax),2])
	Br_norm = np.zeros([len(timeax),2])
	viandva = np.zeros([len(timeax),2])
	ma = np.zeros(len(timeax))
	pressures = np.zeros([len(timeax),3])
	parperps = np.zeros([len(timeax),2])
	line = np.ones_like(ni)
	theta = np.zeros([len(timeax),2])
	theta_v = np.zeros([len(timeax),2])

	# E_antennas = np.zeros_like(voltages)

	for i in range(len(Bvecs)):
		Br[i], Bt[i], Bn[i] = Bvecs[i,0], Bvecs[i,1], Bvecs[i,2]
		vr[i], vt[i], vn[i] = vivecs[i,0], vivecs[i,1], vivecs[i,2]
		B_mag[i] = np.linalg.norm(Bvecs[i,:])
		v_mag[i] = np.linalg.norm(vivecs[i,:])
		
		E_conv[i] = -get_vxb(vivecs[i],Bvecs[i])
		# E_antennas[i] = voltages[i]/3.5
		S[i] = get_S(E_conv[i],Bvecs[i]) #RTN
		K[i] = get_K(mi,ni[i],vivecs[i]) #RTN
		H[i] = get_H(vxyz[i],ni[i]*TiTensor[i]) #XYZ

		Tpar[i],Tperp[i],Ppar[i],Pperp[i] = get_parperps(ni[i],TiTensor[i],Bxyz[i])

		entropy[i] = get_entropy(Ti[i],ni[i])
		
		ExB[i] = get_vxb(E_conv[i],Bvecs[i]) # Will probably need to revisit under better assumptions
		Br_norm[i] = [get_brnorm(Bvecs[i]),0]
		vr_norm[i] = get_brnorm(vivecs[i])
		P_mag[i] = get_pm(Bvecs[i])
		P_th[i] = get_pth(ni[i],Ti[i])
		va[i] = get_va(Bvecs[i],ni[i],mi)
		vs[i] = get_vs(Ti[i],mi)
		vth[i]  = get_vth(Ti[i],mi)
		ma[i] = v_mag[i]/va[i]
		beta[i] = [P_th[i]/P_mag[i],1]
		viandva[i] = [np.linalg.norm(vivecs[i]),va[i]]
		pressures[i] = [P_mag[i] + P_th[i],P_mag[i],P_th[i]]
		parperps[i] = [Ppar[i], Pperp[i]]
		beta_par[i] = Ppar[i]/P_mag[i]




	# FOV filtering
	FOV_filter()



 #------------------------------------------------------------------------
	# Get averaged (or median) quantities.  meaninterval = 1 min ##
	minutes = 60

	# Means
	Bvecs_mean = get_vecmean(Bvecs,minutes*meaninterval)
	Bvecs_bigmean = get_vecmean(Bvecs,2*minutes*meaninterval)
	Bxyz_mean = get_vecmean(Bxyz,minutes*meaninterval)
	vivecs_mean = get_vecmean(vivecs,minutes*meaninterval)
	S_mean = get_vecmean(S,minutes*meaninterval)
	K_mean = get_vecmean(K,minutes*meaninterval)
	H_mean = get_vecmean(H,minutes*meaninterval)
	Bmag_mean = get_mean(B_mag,minutes*meaninterval)
	Bmag_bigmean = get_mean(B_mag,2*minutes*meaninterval)
	vmag_mean = get_mean(v_mag,minutes*meaninterval)
	n_mean = get_mean(ni,minutes*meaninterval)
	# ma_mean = get_mean(v_mag/va,minutes*meaninterval)


	# # Medians
	#Bvecs_med = get_vecmed(Bvecs,minutes*meaninterval)
	# Bxyz_mean = get_vecmed(Bxyz,minutes*meaninterval)
	# vivecs_mean = get_vecmed(vivecs,minutes*meaninterval)
	#Bmag_mean = get_med(B_mag,minutes*meaninterval)
	#vmag_mean = get_med(v_mag,minutes*meaninterval)
	# n_mean = get_med(ni,minutes*meaninterval)
	# beta_mean = get_med(beta[:,0], minutes*meaninterval)
	# S_mean = get_vecmed(S,minutes*meaninterval)
	# K_mean = get_vecmed(K,minutes*meaninterval)
	# H_mean = get_vecmed(H,minutes*meaninterval)
	
	
	
    # Get va and Ma using the chosen n, B	
	# va_mean = Bmag_mean/np.sqrt(mu0*mi*n_mean)
	va_mean = Bmag_mean/np.sqrt(mu0*mi*n_mean)
	# va_mean = abs(Bvecs_mean[:,0])/np.sqrt(mu0*mi*ni) #Uses Br instead of B_mag
	ma_mean = vmag_mean/va_mean 
	beta_mean = get_mean(beta[:,0], minutes*meaninterval)
	
	P_mean = get_mean(P_th + P_mag,minutes*meaninterval) 
	Pth_mean = get_mean(P_th,minutes*meaninterval) 
	Pmag_mean = get_mean(P_mag,minutes*meaninterval) 

#-----------------------------------------------------------------------



	for i in range(len(Bvecs)):
		E_conv[i] = -get_vxb(vivecs[i],Bvecs[i])
		S[i] = get_S(E_conv[i],Bvecs[i]) #RTN
		K[i] = get_K(mi,ni[i],vivecs[i]) #RTN



	# Get dB,dv + components, etc. and vB alignment variable (proxy for alfvenicity)
	dB,dB_norm = np.zeros_like(Bvecs),np.zeros_like(Bvecs)
	dv,dv_norm = np.zeros_like(vivecs),np.zeros_like(vivecs)
	dS,dS_norm = np.zeros_like(vivecs),np.zeros_like(vivecs)
	dK,dK_norm = np.zeros_like(vivecs),np.zeros_like(vivecs)
	dH,dH_norm = np.zeros_like(vivecs),np.zeros_like(vivecs)
	dBmag = np.zeros_like(ni)
	dvmag = np.zeros_like(ni)
	dSmag = np.zeros_like(ni)
	dKmag = np.zeros_like(ni)
	dHmag = np.zeros_like(ni)
	Smag = np.zeros_like(ni)
	Kmag = np.zeros_like(ni)
	Hmag = np.zeros_like(ni)
	
	S_par = np.zeros_like(ni)
	K_par = np.zeros_like(ni)
	H_par = np.zeros_like(ni)



	dB_par,dB_perp = np.zeros_like(ni),np.zeros_like(ni)
	dv_par,dv_perp = np.zeros_like(ni),np.zeros_like(ni)

	dv_par1, dv_perp1,dv_perp2 = np.zeros_like(ni),np.zeros_like(ni),np.zeros_like(ni)
	dv_par1_norm,dv_perp1_norm,dv_perp2_norm = np.zeros_like(ni),np.zeros_like(ni),np.zeros_like(ni)

	dn,dn_norm = np.zeros_like(ni),np.zeros_like(ni)

	dP,dP_norm = np.zeros_like(ni),np.zeros_like(ni)
	dPth,dPth_norm = np.zeros_like(ni),np.zeros_like(ni)
	dPm,dPm_norm = np.zeros_like(ni),np.zeros_like(ni)


	B_par = np.zeros_like(ni)
	v_par = np.zeros_like(ni)
	dB_par_norm = np.zeros_like(ni)
	dB_perp_norm = np.zeros_like(ni)
	dv_par = np.zeros_like(ni)
	dv_par_norm = np.zeros_like(ni)
	dv_perp_norm = np.zeros_like(ni)
	dB_norm_mag,dv_norm_mag = np.zeros_like(ni),np.zeros_like(ni)
	vB_alignment = np.zeros_like(ni)
	dvdB_alignment = np.zeros_like(ni)
	angle = np.zeros_like(ni)
	sigma_c = np.zeros_like(ni)
	sigma_r = np.zeros_like(ni)
	magcomp = np.zeros_like(ni)

	for i in range(len(timeax)):
		dvdB_alignment[i] = np.abs(np.dot(dB[i],dv[i]))/(np.linalg.norm(dB[i])*np.linalg.norm(dv[i]))
		vB_alignment[i] = np.abs(np.dot(Bvecs[i],vivecs[i]))/(B_mag[i]*v_mag[i]) #this one seems most useful

		# Magnetic deflection angle
		angle[i] = get_angle(Bvecs[i],Bvecs_mean[i])


		B_par[i] = get_par(Bvecs[i],Bvecs_mean[i]) # All par to mean B 
		v_par[i] = get_par(vivecs[i],Bvecs_mean[i]) 



		
		# dB,dv & components
		dB[i], dB_norm[i]= get_delta(Bvecs[i],Bvecs_mean[i])  # dB & dB/|B| (vectors)
		dv[i], dv_norm[i]= get_delta(vivecs[i],vivecs_mean[i]) # dv & dv/|v| (vectors)
		
		dH[i], dH_norm[i]= get_delta(H[i],H_mean[i]) # dH & dH/|H| (vectors)
		dS[i], dS_norm[i]= get_delta(S[i],S_mean[i]) # dS & dS/|S| (vectors)
		dK[i], dK_norm[i]= get_delta(K[i],K_mean[i]) # dK & dK/|K| (vectors)
		
		dBmag[i] = np.linalg.norm(dB[i])
		dvmag[i] = np.linalg.norm(dv[i])
		dSmag[i] = np.linalg.norm(dS[i])
		dKmag[i] = np.linalg.norm(dK[i])
		dHmag[i] = np.linalg.norm(dH[i])

		Smag[i] = np.linalg.norm(S[i])
		Kmag[i] = np.linalg.norm(K[i])
		Hmag[i] = np.linalg.norm(H[i])

		S_par[i] = get_par(S[i], Bvecs_mean[i])
		K_par[i] = get_par(K[i], Bvecs_mean[i])
		H_par[i] = get_par(H[i], Bxyz_mean[i])
		
		dB_norm_mag[i] = np.linalg.norm(dB_norm[i]) # |dB|/|B| (scalar)
		dv_norm_mag[i] = np.linalg.norm(dv_norm[i]) # |dv|/|v| (scalar)
		dB_par[i] = get_par(dB[i],Bvecs_mean[i])
		dv_par[i] = get_par(dv[i],Bvecs_mean[i])



		dB_perp[i] = np.sqrt(np.linalg.norm(dB[i])**2 - dB_par[i]**2)

		magcomp[i] = (dB_par[i]**2) / (dB_par[i]**2 + dB_perp[i]**2)

		#########
		dv_perp[i] = np.sqrt(np.linalg.norm(dv[i])**2 - dv_par[i]**2)
		#dv_perp[i] = dv[i,1]  #Using dvn as proxy

		dv_par1[i], dv_perp1[i],dv_perp2[i] = get_vperps(dv[i],vivecs_mean[i],Bvecs_mean[i])
		##########

		dB_par_norm[i] = dB_par[i]/Bmag_mean[i]
		dv_par_norm[i] = dv_par[i]/vmag_mean[i]
		dB_perp_norm[i] = dB_perp[i]/Bmag_mean[i]
		dv_perp_norm[i] = dv_perp[i]/vmag_mean[i]
		
		dv_par1_norm[i] = dv_par1[i]/vmag_mean[i]
		dv_perp1_norm[i] = dv_perp1[i]/vmag_mean[i]
		dv_perp2_norm[i] = dv_perp2[i]/vmag_mean[i]

		dn[i],dn_norm[i] = get_deltascalar(ni[i],n_mean[i])
		
		dP[i],dP_norm[i] = get_deltascalar(P_th[i] + P_mag[i],P_mean[i])
		dPth[i],dPth_norm[i] = get_deltascalar(P_th[i],Pth_mean[i])
		dPm[i],dPm_norm[i] = get_deltascalar(P_mag[i],Pmag_mean[i])

		# sigma_r[i] = get_residenergy(vivecs[i],Bvecs[i],ni[i],mi)
		# sigma_c[i] = get_crosshelicity(vivecs[i],Bvecs[i],ni[i],mi)

		sigma_r[i] = get_residenergy(dv[i],dB[i],ni[i],mi)
		sigma_c[i] = get_crosshelicity(dv[i],dB[i],ni[i],mi)




	low =0
	high = 180
	angle_reduced = np.zeros_like(angle)
	for i in range(len(timeax)): angle[i] = get_angle(Bvecs[i],Bvecs_mean[i])
	angle_reduced = filter_angle(angle,low,high) #If you want all angles, do 0,180




	allangles = np.concat([allangles,angle_reduced])
	allmachs = np.concat([allmachs,ma_mean])
	allpositions = np.concat([allpositions,position])
	
	allvr = np.concat([allvr,vr])
	allvt = np.concat([allvr,vt])
	allvn = np.concat([allvr,vn])
	allBr = np.concat([allvr,Br])
	allBt = np.concat([allvr,Bt])
	allBn = np.concat([allvr,Bn])

	allSr = np.concat([allSr,S[:,0]])
	allKr = np.concat([allKr,K[:,0]])
	allHr = np.concat([allHr,H[:,0]])
	allbetas = np.concat([allbetas,beta_mean])
	alldBmag = np.concat([alldBmag,dBmag])
	alldvmag = np.concat([alldvmag,dvmag])

	alldSmag = np.concat([alldSmag,dSmag])
	alldKmag = np.concat([alldKmag,dKmag])
	alldHmag = np.concat([alldKmag,dKmag])
	allSmags = np.concat([allSmags,Smag])
	allKmags = np.concat([allKmags,Kmag])
	allHmags = np.concat([allHmags,Hmag])

	allSpar = np.concat([allSpar,S_par])
	allKpar = np.concat([allKpar,K_par])
	allHpar = np.concat([allHpar,H_par])
	
	alldB_norm_mag = np.concat([alldB_norm_mag,dB_norm_mag])
	alldv_norm_mag = np.concat([alldv_norm_mag,dv_norm_mag])
	alldB_par_norm = np.concat([alldB_par_norm,dB_par_norm])
	alldv_par_norm = np.concat([alldv_par_norm,dv_par_norm])
	alldB_perp_norm = np.concat([alldB_perp_norm,dB_perp_norm])
	alldv_perp_norm = np.concat([alldv_perp_norm,dv_perp_norm])
	alldv_par1_norm = np.concat([alldv_par1_norm,dv_par1_norm])
	alldv_perp1_norm = np.concat([alldv_perp1_norm,dv_perp1_norm])
	alldv_perp2_norm = np.concat([alldv_perp2_norm,dv_perp2_norm])

	allmagcomp = np.concat([allmagcomp, magcomp])

	

	alldn_norm = np.concat([alldn_norm,dn_norm])
	alldP_norm = np.concat([alldP_norm,dP_norm])
	alldPth_norm = np.concat([alldPth_norm,dPth_norm])
	alldPm_norm = np.concat([alldPm_norm,dPm_norm])

	allT = np.concat([allT,Ti])
	alln = np.concat([alln,ni])
	allTpar = np.concat([allTpar,Tpar])
	allTperp = np.concat([allTperp,Tperp])
	allTparperp = np.concat([allTparperp,Tpar/Tperp])
	allbeta_par = np.concat([allbeta_par, beta_par])
	
	
	allBmags = np.concat([allBmags,Bmag_mean])
	allvmags = np.concat([allvmags,vmag_mean])
	allva = np.concat([allva,va_mean])
	allva_inst = np.concat([allva_inst,va])
	# allvbalignment = np.concat([allvbalignment,vB_alignment])
	allcrosshelicity = np.concat([allcrosshelicity,sigma_c])
	allresidenergy = np.concat([allresidenergy,sigma_r])
	allentropy = np.concat([allentropy,entropy])
	


# %%



#%%

r_avg = np.mean(allpositions)
r_norm = allpositions/r_avg

norm_factor = (1/r_norm)#**2 # Normalized 1/r^2 = (<r>/r)^2 


#Quick normalized energy flux fluctuation params 
dK_norm_mag = alldKmag/allKmags
dS_norm_mag = alldSmag/allSmags





# def filter_combined(var,def_floor,sigma_threshold,v_threshold):
# 	floor = def_floor
# 	thr1 = sigma_threshold
# 	thr2 = v_threshold

# 	filtered_var = filter_deflections(var,floor)
# 	slowvar,fastvar,nonvar = filter_winds(filtered_var,thr1,thr2)

# 	return slowvar,fastvar,nonvar

# floor,th1,th2 = 0.05,0.8,3.5e5
# slowent,fastent,nonent = filter_combined(allentropy,floor,th1,th2)
# slowSr,fastSr,nonSr = filter_combined(allSr,floor,th1,th2)
# slowKr,fastKr,nonKr = filter_combined(allKr,floor,th1,th2)


floor = 0.05 # Choose the minimum deflection magnitude 

# Filter all vars of interest (probably only need a few)
filtered_angle = filter_deflections(allangles,floor)
filtered_mach = filter_deflections(allmachs,floor)
# filtered_ent = filter_deflections(allentropy,floor)



# slowmachs,fastmachs,nonmachs = filter_winds(filtered_mach,0.8,3.5e5)
# slowangs,fastangs,nonangs = filter_winds(filtered_angle,0.8,3.5e5)
# slowent,fastent,nonent = filter_winds(filtered_ent,0.8,3.5e5)

#%%

#%%


#%%

x = np.log10(filtered_mach) # x axis
y = allangles # y axis
z = allSmags # variable for colorbar
#z = allvr*(allBt_i**2 + allBn_i**2)/mu0

# z = alldv_perp_norm # variable for colorbar

# Define bins
bins=30
x_bins = np.linspace(-1.2, 1.2, bins)
y_bins = np.linspace(0, 180, bins)

# Calculate the mean of 'z' for each bin
mean_z, x_edge, y_edge, bin_number = binned_statistic_2d(x, y, z, statistic='mean', bins=[x_bins, y_bins])

# Plot the result
plt.imshow(mean_z.T, origin='lower', extent=[x_edge[0], x_edge[-1], y_edge[0], y_edge[-1]], cmap='plasma',aspect='auto',vmin=0.,vmax=0.2)
plt.colorbar(label='|S|')
plt.xlim(-1.2,1.2)
plt.ylim(0,160)
plt.xlabel('Log10(Ma)')
plt.ylabel('Deflection Angle (degrees)')
plt.axhline(y=90,color='k')
plt.axvline(x=0,color='k')
plt.show()


z = allKmags # variable for colorbar
# z = alldv_perp_norm # variable for colorbar

# Define bins

# Calculate the mean of 'z' for each bin
mean_z, x_edge, y_edge, bin_number = binned_statistic_2d(x, y, z, statistic='mean', bins=[x_bins, y_bins])

# Plot the result
plt.imshow(mean_z.T, origin='lower', extent=[x_edge[0], x_edge[-1], y_edge[0], y_edge[-1]], cmap='plasma',aspect='auto',vmin=0.,vmax=0.2)
plt.colorbar(label='|K|')
plt.xlim(-1.2,1.2)
plt.ylim(0,160)
plt.xlabel('Log10(Ma)')
plt.ylabel('Deflection Angle (degrees)')
plt.axhline(y=90,color='k')
plt.axvline(x=0,color='k')
plt.show()

z = allHmags # variable for colorbar
# z = alldv_perp_norm # variable for colorbar

# Define bins

# Calculate the mean of 'z' for each bin
mean_z, x_edge, y_edge, bin_number = binned_statistic_2d(x, y, z, statistic='mean', bins=[x_bins, y_bins])

# Plot the result
plt.imshow(mean_z.T, origin='lower', extent=[x_edge[0], x_edge[-1], y_edge[0], y_edge[-1]], cmap='plasma',aspect='auto',vmin=0.,vmax=0.2)
plt.colorbar(label='|H|')
plt.xlim(-1.2,1.2)
plt.ylim(0,160)
plt.xlabel('Log10(Ma)')
plt.ylabel('Deflection Angle (degrees)')
plt.axhline(y=90,color='k')
plt.axvline(x=0,color='k')
plt.show()

# %%
#%%

#%%



#%%

fig,ax = plt.subplots(1,2,figsize=(10,5))
ax[0].scatter(np.log10(allmachs),allresidenergy,s=0.01)
ax[1].scatter(np.log10(allmachs),allSr/allKr,s=0.01)
ax[0].set_xlim(-0.8,0.8)
ax[1].set_xlim(-0.8,0.8)
ax[0].set_ylim(0,10)
ax[1].set_ylim(0,10)
# plt.axhline(y=0,color='k')
ax[0].axvline(x=0,color='k')
ax[1].axvline(x=0,color='k')
ax[0].axhline(y=1,color='k')
ax[1].axhline(y=1,color='k')
ax[0].set_ylabel("dS/dK")
ax[1].set_ylabel("Sr/Kr")
ax[0].set_xlabel("Log10(Ma)")
ax[1].set_xlabel("Log10(Ma)")




#%%
# Scatterplot - dv/va and dB/B
fig,ax = plt.subplots(1,3,figsize=(14,3))
ax[0].scatter(np.log10(allmachs),alldB_norm_mag,s=0.01)
ax[0].set_ylabel("dB/B")
ax[0].set_xlim(-0.8,0.8)
ax[0].set_ylim(0.,2)
ax[0].axvline(x=0,color='k')
ax[0].set_xlabel("Log10(Ma)")
ax[0].axhline(y=1,linestyle='dotted',color='k')


ax[2].scatter(np.log10(allmachs),alldv_norm_mag,s=0.01)
ax[2].set_xlim(-0.8,0.8)
ax[2].set_ylim(0.,2)
ax[2].axvline(x=0,color='k')
ax[2].set_ylabel("dv/v")
ax[2].set_xlabel("Log10(Ma)")
ax[2].axhline(y=1,linestyle='dotted',color='k')


ax[1].scatter(np.log10(allmachs),alldvmag/allva,s=0.01)
ax[1].set_xlim(-0.8,0.8)
ax[1].set_ylim(0.,2)
ax[1].axvline(x=0,color='k')
ax[1].set_ylabel("dv/va")
ax[1].set_xlabel("Log10(Ma)")
ax[1].axhline(y=1,linestyle='dotted',color='k')






# %%





# %%

#%%

# beta dependence on Ma and deflection angle

x = np.log10(allmachs) # x axis
y = allangles # y axis
z = alldSmag # variable for colorbar

# Define bins
bins=50
x_bins = np.linspace(-1, 1, bins)
y_bins = np.linspace(0, 180, bins)

# Calculate the mean of 'z' for each bin
mean_z, x_edge, y_edge, bin_number = binned_statistic_2d(x, y, z, statistic='mean', bins=[x_bins, y_bins])

# Plot the result
plt.imshow(mean_z.T, origin='lower', extent=[x_edge[0], x_edge[-1], y_edge[0], y_edge[-1]], cmap='plasma',aspect='auto',vmin=0.,vmax=0.1)
plt.colorbar(label='dS')
plt.xlim(-0.6,1)
plt.ylim(20,160)
plt.xlabel('Log10(Ma)')
plt.ylabel('Deflection Angle (degrees)')
plt.axhline(y=90,color='k')
plt.axvline(x=0,color='k')
plt.show()
# %%

# %%
#%%
# Plot 2 
var = alldvmag/allva
var2 = alldv_norm_mag
range=[0,2]
bins=40
plt.hist(np.log10(alldB_norm_mag),bins=bins,range=range,log=True,histtype="step")
plt.hist(np.log10(var),bins=bins,log=True,range=range,histtype="step")
plt.hist(np.log10(var2),bins=bins,log=True,range=range,histtype="step")
plt.ylim(0.1,1000000)
#%%

#%%
dataplot()

#%%

fig, ax = plt.subplots(1,1,figsize=(8,8))
ax.hist2d(np.log10(allmachs),allSmags,range=[[-0.5,0.5],[0,0.6]],bins=100,cmin=1,cmap='viridis')
ax.set_xlabel('Log10 Alfven Mach Number')
ax.set_ylabel('Radial Poynting Flux (SI units W/m^2)')
ax.set_facecolor('k')
ax.axhline(y=90,color='w')
ax.axvline(x=0,color='w')
ax.set_xlim(-0.5,0.5)
#%%


fig, ax = plt.subplots(1,1,figsize=(8,8))
ax.hist2d(np.log10(allmachs),allKmags,range=[[-0.5,0.5],[0,0.6]],bins=100,cmin=1,cmap='viridis')
ax.set_xlabel('Log10 Alfven Mach Number')
ax.set_ylabel('Radial Kinetic Energy Flux (SI units W/m^2)')
ax.set_facecolor('k')
ax.axhline(y=90,color='w')
ax.axvline(x=0,color='w')
ax.set_xlim(-0.5,0.5)


fig, ax = plt.subplots(1,1,figsize=(8,8))
ax.hist2d(np.log10(allmachs),allHmags,range=[[-0.5,0.5],[0,0.6]],bins=100,cmin=1,cmap='viridis')
ax.set_xlabel('Log10 Alfven Mach Number')
ax.set_ylabel('Radial Kinetic Energy Flux (SI units W/m^2)')
ax.set_facecolor('k')
ax.axhline(y=90,color='w')
ax.axvline(x=0,color='w')
ax.set_xlim(-0.5,0.5)



#%%
# %%


# %%
def interp(arr1,arr2):
	# Original array
	original_array = arr1
	original_length = len(arr1)

	# Desired new length
	new_length = len(arr2)

	# Create x-coordinates for the original array
	x_original = np.linspace(0, original_length - 1, original_length)

	# Create x-coordinates for the new, interpolated array
	x_new = np.linspace(0, original_length - 1, new_length)

	# Perform interpolation
	interpolated_array = np.interp(x_new, x_original, original_array)

	return interpolated_array

#print("Original array:", original_array)
#
# print("Interpolated array:", interpolated_array)




# %%

allBr_i = interp(allBr,allvr)
allBt_i = interp(allBt,allvr)
allBn_i = interp(allBn,allvr)