# Script to obtain Energy fluxes for all spacecraft and get their Div/Curl

# Produces all quantities without error, but this is only a frst pass.  Might be some hidden issues
#%%

#%%
import pyspedas
import pandas as pd
from pytplot import tplot
import numpy as np
from pytplot import options
from pyspedas import tinterpol
from pytplot import tplot_options
from pyspedas.mms import mec,fgm,fpi,edp,scm
from pytplot import get_data, store_data,timespan
# from pyspedas.analysis.tsmooth import tsmooth
from matplotlib.pyplot import plot
from matplotlib.pyplot import scatter

me = 9.1094e-31 #kg
mi = 1837*me
mu0 = 1.2566370e-06  #;m kg / C^2
eps0 = 8.85e-12   # C^2/Nm^2
e = 1.602e-19 #C

lmn_1016 = np.array([
[0.3665, -0.1201, 0.9226],
[0.5694, -0.7553, -0.3245],
[0.7358, 0.6443, -0.2084]])
lmn_0711 = np.array([
[0.94822, -0.25506, -0.18926],
[0.18182, 0.92451, -0.334996],
[0.26042, 0.28324, 0.92301]])
lmn_0617 = np.array([
[0.93, 0.3, -0.2],
[-0.27, 0.2, -0.94],
[-0.24, 0.93, 0.27]])
lmn_1209 = np.array([
[-0.091,0.87,0.49],
[-0.25,-0.49,0.83],
[0.96,-0.05,0.2700]])
# lmn_0128 = np.array([
# [-0.315,-0.6516,-0.6896],
# [0.001,-0.7269,0.6863],
# [-0.949,0.2156,0.2298]])


lmn_0128 = np.array([    # from wilder (2018)
[0.29,0.91,0.51],
[-0.35,0.59,-0.72],
[-0.89,0.03,0.45]])

lmn_0810 = np.array([
[0.985, -0.141, 0.097],
[0.152, 0.982, -0.109],
[-0.080, 0.122, 0.989]])

I = np.identity(3)


frame = lmn_1209  # Choose the LMN coordinate system you want to use.  Choose I to keep GSE

# Get Data
# trange = ['2017-08-10/12:18:00', '2017-08-10/12:19:00']
# timespan('2017-08-10 12:18:20', 20, keyword='seconds')


# trange = ['2017-07-11/22:33:30', '2017-07-11/22:34:30']
# timespan('2017-07-11 22:34:00', 6, keyword='seconds')


# Hubbert Events (aside from June 17th)
# trange = ['2017-07-20/09:59:00 2017-07-20/09:10:00']
# trange = ['2017-06-19/09:43:00', '2017-06-19/09:45:00']


# trange = ['2016-11-09/13:39:00', '2016-11-09/13:40:00']  # Shock.  Sep ~ 16 de


trange = ['2016-12-09/09:03:54', '2016-12-09/09:03:54.6']  # Electron. Sep ~ 5 de|
# trange = ['2016-12-09/09:03:30', '2016-12-09/09:03:54.6']  # Electron. Sep ~ 5 de
timespan('2016-12-09 09:03:54',0.6, keyword='seconds')

# trange = ['2017-01-28/09:08:55', '2017-01-28/09:09:05']
# timespan('2017-01-28 09:08:55', 10, keyword='seconds')

# trange = ['2017-01-28/05:32:30', '2017-01-28/05:32:34']
# timespan('2017-01-28 05:32:31', 2.5, keyword='seconds')
#trange = ['2016-11-09/13:38:00', '2016-11-09/13:40:00']
# trange = ['2017-06-17/20:23:30', '2017-06-17/20:24:30']

# trange = ['2017-06-17/20:24:05', '2017-06-17/20:24:10']
# timespan('2017-06-17 20:24:05', 5, keyword='seconds')
# trange = ['2015-10-16/13:06:30', '2015-10-16/13:07:20']
# timespan('2015-10-16 13:06:50', 20, keyword='seconds')

probe  = [1,2,3,4] 
trange = trange
#%%
#fpi_vars = fpi(probe = probe,data_rate = 'brst',trange=trange,time_clip=False)
# Weird FPI error.  It loads the data, but throws an error that interrupts.  Just wait for it to interuppt then 
# then run the remaining blocks.  

#%%
fgm_vars = fgm(probe = probe, data_rate = 'brst', trange=trange,time_clip=True)
edp_vars = edp(probe = probe,data_rate = 'brst',trange=trange,time_clip=True) 
fpi_vars = fpi(probe = probe,data_rate = 'brst',trange=trange,time_clip=True)
mec_vars = mec(probe = probe,trange=trange,data_rate='brst',time_clip=True)

#%%

# interpfld = 'mms1_fgm_b_gse_brst_l2'
interpfld = 'mms1_edp_dce_gse_brst_l2'
#interpfld = 'mms1_dis_numberdensity_brst'
#interpfld = 'mms1_des_numberdensity_brst'
# Function to change shape of fields to np arrays (will change them back eventually)
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

# LMN Conversion function
def convert_lmn(gse, lmn):
    vec = np.zeros_like(gse)
    for i in range(3):
        vec[i] = np.dot(lmn[i],gse)
    return vec

def Eprime(E,v,B):
	Ep = E + np.cross(v,B)
	return Ep


def get_j(n,vi,ve):
	q = 1.6e-19 # charge unit in SI
	j = q*n*(vi - ve)
	return j
# Functions for energy densites and fluxes #
# Still unsure, but should be correct within ord of mag
# Poynting S = ExB
# Kinetic = 0.5*(n*mv^2)*v
# Enthalpy = 0.5*v*Tr(P)  + v dot P
# Q = ?
def therm_dens(P):
	ut = 0.5*np.trace(P)
	return ut

def kinetic_dens(m,n,v):
	uk = 0.5*m*n*np.linalg.norm(v)**2
	return uk

def B_dens(B):
	# um = 0.5*(np.linalg.norm(B)**2)*mu0**-1
	#um = 0.5*(np.linalg.norm(B)**2)/(mu0)
	um = 0.5*( np.linalg.norm(B[0])**2 + np.linalg.norm(B[1])**2 + np.linalg.norm(B[2])**2)/(mu0)
	return um

def E_dens(E):
	ue = 0.5*eps0*np.linalg.norm(E)**2
	return ue


# REMEMBER TO FIX!
def kin_flux(m,n,v):
	K = 0.5*m*n*v**3
	#K = 0.5*m*n*v*np.linalg.norm(v)**2# this uses |v|*v
	# K = v
	return K

# This works now even though its literally identical to the old version that didnt work..
def enth_flux(v,P):
	H = 0.5*v*np.trace(P) + np.dot(v,P)
	return H

def Poynt_flux(E,B):
	S = np.cross(E,B)/mu0
	return S

# Energy Flux Transport
# Thermal
def Therm_Transport(v,P):
	H = enth_flux(v,P)
	u = therm_dens(P)
	tft = H/u
	return tft

def Kin_Transport(m,n,v):
	K = kin_flux(m,n,v)
	u = kinetic_dens(m,n,v)
	kft = K/u
	return kft

def EM_Transport(E,B):
	S = Poynt_flux(E,B)
	uB = B_dens(B)
	uE = E_dens(E)
	emft = S/(uB+uE)
	return emft


# Add MFT 
def MFT(E,B):
	pass
	return None
#


# Under Construction
#///////////////////////////////////////////////////////////////////////
def PStrain(Plist,vlist,klist):
	i = 0
	vmlist = []
	pmlist=[]
	var = np.zeros(1)
	P = sum(Plist)/4
	for i in range(4):
		vmlist.append(np.array([vlist[i]]))
		pmlist.append(np.array([Plist[i]]))
		var = var + np.dot(pmlist[i],klist[i]*vmlist[i].T)
	
	Pst = np.dot(P,grad(vlist,klist)) # version 1, take average pressure, dot with the gradient
	Pst = var   # version 2, pressure taken from each sc, dot product for each SC

	return Pst
#///////////////////////////////////////////////////////////////////////////



def curv(veclist, klist):      # Modified version of the gradient calculation to get the curvature (1/Rc)
    i = 0
    mlist = []  # We need to turn these regular arrays into matrix objects
    vec = np.zeros(3)
    B_avg = sum(veclist)/4
    for i in range(4):
        mlist.append(np.array([veclist[i]]))
        vec = vec + np.dot(B_avg/np.linalg.norm(B_avg),klist[i]*mlist[i].T/np.linalg.norm(mlist[i]))
    crv = np.linalg.norm(vec)    
    return crv





# Functions for Power densities (J.E, J.E', -PdotGrad(u), du/dt, Div(S))
def get_j(n,vi,ve):
	q = 1.602e-19 # charge unit in SI
	j = q*n*(vi - ve)
	return j

# Put a J.E function here!!
	
def sepvec(a,b):
    sepvec = b-a
    return sepvec

# Gets a list of the reciprocal vectors [k1,k2,k3,k4] in LMN (give it GSE)
def recip_vecs(pos1,pos2,pos3,pos4):
    
    # Define every possible separation vector
    # There has to be a cleaner way to do this part.  Keep for now
    pos12 = sepvec(pos1,pos2)
    pos13 = sepvec(pos1,pos3)
    pos14 = sepvec(pos1,pos4)
    pos23 = sepvec(pos2,pos3)
    pos24 = sepvec(pos2,pos4)
    pos34 = sepvec(pos3,pos4)
    pos21 = -pos12
    pos31 = -pos13
    pos41 = -pos14
    pos32 = -pos23
    pos42 = -pos24
    pos43 = -pos34

    # Function to calculate k vector given separation vectors a,b,c
    # where the form is k = axb/(c dot axb)
    def kvec(a,b,c):
        crossp = np.cross(a,b)
        denom = np.dot(c,crossp)
        k = crossp/denom
        return k 

    # Get each k vector and put in a list  (Also convert Kvec to LMN coord)
    k1 = convert_lmn(kvec(pos23,pos24,pos21),frame)
    k2 = convert_lmn(kvec(pos34,pos31,pos32),frame)
    k3 = convert_lmn(kvec(pos41,pos42,pos43),frame)
    k4 = convert_lmn(kvec(pos12,pos13,pos14),frame)

    klist = np.array([k1,k2,k3,k4])

    return klist


# Define divergence given a veclist and klist
# where veclist is a list of some vector quantity measured at [MMS1,MMS2,MMS3,MMS4]
# and klist is the list of reciprocal vectors [k1,k2,k3,k4]
def grad(veclist, klist):
    i = 0
    mlist = []  # We need to turn these regular arrays into matrix objects
    grd = np.zeros([3,3])
    for i in range(4):
        mlist.append(np.array([veclist[i]]))
        grd = grd + klist[i]*mlist[i].T #Fix this too!!
    return grd

def grad_scalar(scalist, klist):
    i = 0
    grd = np.zeros(3)
    for i in range(4):
        grd = grd + klist[i]*scalist[i] #Fix this too!!
    return grd

def div(veclist, klist):
	i = 0
	div,divx,divy,divz = 0,0,0,0
	for i in range(4):
		divx = divx + klist[i,0]*veclist[i][0]
		divy = divy + klist[i,1]*veclist[i][1]
		divz = divz + klist[i,2]*veclist[i][2]
		div = div + np.dot(klist[i],veclist[i])
	# div = divx+divy+divz
	return np.array([divx,divy,divz,div])  

#////////////////////////////////////////////////////
def tensdiv(tenslist, klist):
	i = 0
	div = np.zeros(3)
	dxxx,dxyx,dxzx,dyxy,dyyy,dyzy,dzxz,dzyz,dzzz = 0,0,0,0,0,0,0,0,0
	for i in range(4):
		dxxx += klist[i,0]*tenslist[i][0,0]
		dxyx += klist[i,0]*tenslist[i][1,0]
		dxzx += klist[i,0]*tenslist[i][2,0]
		
		dyxy += klist[i,1]*tenslist[i][0,1]
		dyyy += klist[i,1]*tenslist[i][1,1]
		dyzy += klist[i,1]*tenslist[i][2,1]

		dzxz += klist[i,2]*tenslist[i][0,2]
		dzyz += klist[i,2]*tenslist[i][1,2]
		dzzz += klist[i,2]*tenslist[i][2,2]
		div = div + np.dot(klist[i],tenslist[i])
	return div 
#/////////////////////////////////////////////////////



def curl(veclist, klist):
    i = 0
    crl = np.array([0,0,0])
    for i in range(4):
        crl = crl + np.cross(klist[i],veclist[i])
    return crl




def curv(veclist, klist):      # Modified version of the gradient calculation to get the curvature (1/Rc)
    i = 0
    mlist = []  # We need to turn these regular arrays into matrix objects
    vec = np.zeros(3)
    B_avg = sum(veclist)/4
    for i in range(4):
        mlist.append(np.array([veclist[i]]))
        vec = vec + np.dot(B_avg/np.linalg.norm(B_avg),klist[i]*mlist[i].T/np.linalg.norm(mlist[i]))
    crv = np.linalg.norm(vec)    
    return crv





##%%


def AllDensities(probe):
	# Field names variables
	B_name = 'mms' + str(probe) + '_fgm_b_gse_brst_l2'
	E_name = 'mms' + str(probe) + '_edp_dce_gse_brst_l2'
	vi_name = 'mms' + str(probe) + '_' + 'dis' + '_bulkv_gse_brst'
	ve_name = 'mms' + str(probe) + '_' + 'des' + '_bulkv_gse_brst'
	B_name = 'mms' + str(probe) + '_fgm_b_gse_brst_l2'
	E_name = 'mms' + str(probe) + '_edp_dce_gse_brst_l2'
	ne_name = 'mms' + str(probe) + '_' + 'des' + '_numberdensity_brst'
	ni_name = 'mms' + str(probe) + '_' + 'dis' + '_numberdensity_brst'
	Pi_name = 'mms' + str(probe) + '_dis_prestensor_gse_brst'
	Pe_name = 'mms' + str(probe) + '_des_prestensor_gse_brst'
	pos_name = 'mms'+ str(probe) + '_mec_r_gse'
	var_name = interpfld
	tinterpol(B_name,var_name, newname='B')
	tinterpol(E_name,var_name, newname='E')
	tinterpol(vi_name,var_name, newname='vi')
	tinterpol(ve_name,var_name, newname='ve')
	tinterpol(Pi_name,var_name, newname='Pi')
	tinterpol(Pe_name,var_name, newname='Pe')
	tinterpol(ni_name,var_name, newname='ni')
	tinterpol(ne_name,var_name, newname='ne')
	tinterpol(pos_name, var_name,newname='pos')
	B = 1e-9*reform(get_data('B'))
	E = 1e-3*reform(get_data('E'))
	vi = 1e3*reform(get_data('vi'))
	ve = 1e3*reform(get_data('ve'))
	ni = 1e6*reform(get_data('ni'))
	ne = 1e6*reform(get_data('ne'))
	Pi = 1e-9*reform(get_data('Pi'))
	Pe = 1e-9*reform(get_data('Pe'))
	ndata = len(vi)
	
	# EM Energy Density
	uE = np.zeros_like(ni)
	uB = np.zeros_like(ni)
	uem = np.zeros_like(ni)
	for i in range(ndata-1):
		uE[i] = E_dens(E[i])  #Array element sequence error?
		uB[i] = B_dens(B[i])
		uem[i] = uE[i] + uB[i]
	
	# Kinetic Energy Density
	uke = np.zeros_like(ni)
	uki = np.zeros_like(ni)
	uk = np.zeros_like(ni)
	for i in range(ndata-1):
		uke[i] = kinetic_dens(me,ne[i],ve[i])
		uki[i] = kinetic_dens(mi,ni[i],vi[i])
		uk[i] = uke[i] + uki[i]

	# Thermal Energy Density
	ute = np.zeros_like(ni)
	uti = np.zeros_like(ni)
	ut = np.zeros_like(ni)
	for i in range(ndata-1):
		ute[i] = therm_dens(Pe[i])
		uti[i] = therm_dens(Pi[i])
		ut[i] = ute[i] + uti[i]

	# Plasma Energy Density
	ue = np.zeros_like(ni)
	ui = np.zeros_like(ni)
	up = np.zeros_like(ni)
	for i in range(ndata-1):
		ue[i] = ute[i] + uke[i]
		ui[i] = uti[i] + uki[i]
		up[i] = ue[i] + ui[i]
	return uem,uE,uB,up,ue,ui,ut,ute,uti,uk,uke,uki

# Function to read data, interpolate, convert to SI units, and come out in the form of np arrays
# Cadence order (high -> low): edp & scm (same) -> fgm -> fpi-des -> fpi-dis
def AllFluxes(probe):
	# frame = lmn_1209  # convert to LMN at the end if desired
	# Field names variables
	B_name = 'mms' + str(probe) + '_fgm_b_gse_brst_l2'
	E_name = 'mms' + str(probe) + '_edp_dce_gse_brst_l2'
	vi_name = 'mms' + str(probe) + '_' + 'dis' + '_bulkv_gse_brst'
	ve_name = 'mms' + str(probe) + '_' + 'des' + '_bulkv_gse_brst'
	B_name = 'mms' + str(probe) + '_fgm_b_gse_brst_l2'
	E_name = 'mms' + str(probe) + '_edp_dce_gse_brst_l2'
	ne_name = 'mms' + str(probe) + '_' + 'des' + '_numberdensity_brst'
	ni_name = 'mms' + str(probe) + '_' + 'dis' + '_numberdensity_brst'
	Pi_name = 'mms' + str(probe) + '_dis_prestensor_gse_brst'
	Pe_name = 'mms' + str(probe) + '_des_prestensor_gse_brst'
	pos_name = 'mms'+ str(probe) + '_mec_r_gse'
	qi_name = 'mms' + str(probe) + '_dis_heatq_gse_brst'
	qe_name = 'mms' + str(probe) + '_des_heatq_gse_brst'
	var_name = interpfld
	tinterpol(B_name,var_name, newname='B')
	tinterpol(E_name,var_name, newname='E')
	tinterpol(vi_name,var_name, newname='vi')
	tinterpol(ve_name,var_name, newname='ve')
	tinterpol(Pi_name,var_name, newname='Pi')
	tinterpol(Pe_name,var_name, newname='Pe')
	tinterpol(ni_name,var_name, newname='ni')
	tinterpol(ne_name,var_name, newname='ne')
	tinterpol(pos_name, var_name,newname='pos')
	tinterpol(qi_name, var_name,newname='qi')
	tinterpol(qe_name, var_name,newname='qe')
	B = 1e-9*reform(get_data('B'))
	E = 1e-3*reform(get_data('E'))
	vi = 1e3*reform(get_data('vi'))
	ve = 1e3*reform(get_data('ve'))
	ni = 1e6*reform(get_data('ni'))
	ne = 1e6*reform(get_data('ne'))
	Pi = 1e-9*reform(get_data('Pi'))
	Pe = 1e-9*reform(get_data('Pe'))
	qe = 1e-3*reform(get_data('qe'))
	qi = 1e-3*reform(get_data('qi'))
	ndata = len(vi)

	
	
	# X-line Velocity
	vxl = np.zeros(3) # Default to zero
	
	if str(frame) == str(lmn_1209):
		vxl = np.array([-42.98e3,246.57e3,-25.56e3]) #ion center


    ###############################
	# Use x-line frame v
	vep = np.zeros_like(ve)
	vip = np.zeros_like(vi)
	for i in range(ndata):
		vep[i] = ve[i] - vxl
		vip[i] = vi[i] - vxl
	ve = vep
	vi = vip
    #################################
    #################################
	# Use x-line frame E
	Ep = np.zeros_like(E)
	for i in range(ndata-1):
		Ep[i] = Eprime(E[i],vxl,B[i][:-1])
	E = Ep
    ####################################


	# Poynting Flux
	S = np.zeros_like(E)
	for i in range(ndata-1):
		S[i] = convert_lmn(Poynt_flux(E[i],B[i,:-1]),frame)

	# Electron Energy Flux
	Ke = np.zeros_like(E)
	He = np.zeros_like(E)
	for i in range(ndata-1):
		He[i] = convert_lmn(enth_flux(ve[i],Pe[i]),frame) #for some reason this works now
		#He[i] = convert_lmn(0.5*ve[i]*np.trace(Pe[i]) + np.dot(ve[i],Pe[i]),frame)
		Ke[i] = convert_lmn(kin_flux(me,ne[i],ve[i]),frame)  # this uses v^3  #unclear which is more correct. Eastwood uses v^3

	# Ion Energy Flux
	Ki = np.zeros_like(E)
	Hi = np.zeros_like(E)
	for i in range(ndata-1):
		Hi[i] = convert_lmn(enth_flux(vi[i],Pi[i]),frame)
		#Hi[i] = convert_lmn(0.5*vi[i]*np.trace(Pi[i]) + np.dot(vi[i],Pi[i]),frame)
		Ki[i] = convert_lmn(kin_flux(mi,ni[i],vi[i]),frame)


	# Heat Flux 
	for i in range(ndata-1):
		qi[i] = convert_lmn(qi[i],frame)
		qe[i] = convert_lmn(qe[i],frame)

	##############

	return S,He,Hi,Ke,Ki,qe,qi

#/////////////////////////////////////////////////////////////////////////////////////////////////////////
def AllTransports(probe):
	# frame = lmn_1209  # convert to LMN at the end if desired
	# Field names variables
	B_name = 'mms' + str(probe) + '_fgm_b_gse_brst_l2'
	E_name = 'mms' + str(probe) + '_edp_dce_gse_brst_l2'
	vi_name = 'mms' + str(probe) + '_' + 'dis' + '_bulkv_gse_brst'
	ve_name = 'mms' + str(probe) + '_' + 'des' + '_bulkv_gse_brst'
	B_name = 'mms' + str(probe) + '_fgm_b_gse_brst_l2'
	E_name = 'mms' + str(probe) + '_edp_dce_gse_brst_l2'
	ne_name = 'mms' + str(probe) + '_' + 'des' + '_numberdensity_brst'
	ni_name = 'mms' + str(probe) + '_' + 'dis' + '_numberdensity_brst'
	Pi_name = 'mms' + str(probe) + '_dis_prestensor_gse_brst'
	Pe_name = 'mms' + str(probe) + '_des_prestensor_gse_brst'
	pos_name = 'mms'+ str(probe) + '_mec_r_gsm'
	var_name = interpfld
	tinterpol(B_name,var_name, newname='B')
	tinterpol(E_name,var_name, newname='E')
	tinterpol(vi_name,var_name, newname='vi')
	tinterpol(ve_name,var_name, newname='ve')
	tinterpol(Pi_name,var_name, newname='Pi')
	tinterpol(Pe_name,var_name, newname='Pe')
	tinterpol(ni_name,var_name, newname='ni')
	tinterpol(ne_name,var_name, newname='ne')
	tinterpol(pos_name, var_name,newname='pos')
	B = 1e-9*reform(get_data('B'))
	E = 1e-3*reform(get_data('E'))
	vi = 1e3*reform(get_data('vi'))
	ve = 1e3*reform(get_data('ve'))
	ni = 1e6*reform(get_data('ni'))
	ne = 1e6*reform(get_data('ne'))
	Pi = 1e-9*reform(get_data('Pi'))
	Pe = 1e-9*reform(get_data('Pe'))
	ndata = len(vi)
	
	# Get Grads
	keft,heft,kift,hift,emft = np.zeros([ndata,3]),np.zeros([ndata,3]),np.zeros([ndata,3]),np.zeros([ndata,3]),np.zeros([ndata,3])

	# "Flux Transports"  (still unclear whether physically meaningful)
	for i in range(ndata-1):
		keft[i] = convert_lmn(Kin_Transport(me,ne[i],ve[i]),frame)
		kift[i] = convert_lmn(Kin_Transport(mi,ni[i],vi[i]),frame)
		heft[i] = convert_lmn(Therm_Transport(ve[i],Pe[i]),frame)
		hift[i] = convert_lmn(Therm_Transport(vi[i],Pi[i]),frame)
		# emft[i] = convert_lmn(EM_Transport(E[i],B[i:-1]),frame)  #Dont know why this part isnt working
	# emft = S1/uB1

	return 	keft,heft,kift,hift,emft

#////////////////////////////////////////////////////////////////////////////////////////////////////	
# All Divs for energy flux divergences
def AllDivs():
	pos_names = [] # Names of positions 
	posits = [] # get_data for each position
	# Get field and mec (position) data
	for i in range(4):
		pos_names.append('mms'+ str(probe[i]) + '_mec_r_gse')
		tinterpol(pos_names[i],interpfld,newname = 'pos' + str(i+1)) #interpolate
		posits.append(get_data('pos' + str(i+1)))
	timeax = posits[0].times
	ndata = len(timeax)
	
	# Reform data into np arrays in SI units (flds and posits)
	for i in range(4):
		posits[i] = 1e3*reform(posits[i]) # [pos1,pos2,pos3,pos4]
	# pos1,pos2,pos3,pos4 = posits[0],posits[1],posits[2],posits[3]
	pos1,pos2,pos3,pos4 = posits[0],posits[1],posits[2],posits[3]
	
	# Get Divs
	div_S,div_Hi,div_He,div_Ki,div_Ke = np.zeros([ndata,4]),np.zeros([ndata,4]),np.zeros([ndata,4]),np.zeros([ndata,4]),np.zeros([ndata,4])
	div_Qi,div_Qe = np.zeros_like(div_Ki),np.zeros_like(div_Ke)
	for i in range(ndata-1):
		Slist = [S1[i],S2[i],S3[i],S4[i]]
		Hilist = [Hi1[i],Hi2[i],Hi3[i],Hi4[i]]
		Helist = [He1[i],He2[i],He3[i],He4[i]]
		Kilist = [Ki1[i],Ki2[i],Ki3[i],Ki4[i]]
		Kelist = [Ke1[i],Ke2[i],Ke3[i],Ke4[i]]
		Qelist = [qe1[i],qe2[i],qe3[i],qe4[i]]
		Qilist = [qi1[i],qi2[i],qi3[i],qi4[i]]
		klist = recip_vecs(pos1[i],pos2[i],pos3[i],pos4[i])

		div_S[i] = div(Slist,klist)
		div_Hi[i] = div(Hilist,klist)
		div_He[i] = div(Helist,klist)
		div_Ki[i] = div(Kilist,klist)
		div_Ke[i] = div(Kelist,klist)
		div_Qi[i] = div(Qilist,klist)
		div_Qe[i] = div(Qelist,klist)

	return div_S,div_He,div_Hi,div_Ke,div_Ki,div_Qe,div_Qi


def AllGrads():
	pos_names = [] # Names of positions 
	posits = [] # get_data for each position
	# Get field and mec (position) data
	for i in range(4):
		pos_names.append('mms'+ str(probe[i]) + '_mec_r_gse')
		tinterpol(pos_names[i],interpfld,newname = 'pos' + str(i+1)) #interpolate
		posits.append(get_data('pos' + str(i+1)))
	timeax = posits[0].times
	ndata = len(timeax)
	
	# Reform data into np arrays in SI units (flds and posits)
	for i in range(4):
		posits[i] = 1e3*reform(posits[i]) # [pos1,pos2,pos3,pos4]
	pos1,pos2,pos3,pos4 = posits[0],posits[1],posits[2],posits[3]

	
	# Get Grads
	grd_uem,grd_up,grd_ue,grd_ui,grd_ute = np.zeros([ndata,3]),np.zeros([ndata,3]),np.zeros([ndata,3]),np.zeros([ndata,3]),np.zeros([ndata,3])
	for i in range(ndata-1):
		uemlist = [uem1[i],uem2[i],uem3[i],uem4[i]]
		uplist = [up1[i],up2[i],up3[i],up4[i]]
		uelist = [ue1[i],ue2[i],ue3[i],ue4[i]]
		uilist = [ui1[i],ui2[i],ui3[i],ui4[i]]
		utelist = [ute1[i],ute2[i],ute3[i],ute4[i]]
		klist = recip_vecs(pos1[i],pos2[i],pos3[i],pos4[i])

		grd_uem[i] =  grad_scalar(uemlist,klist)
		grd_up[i] = grad_scalar(uplist,klist)
		grd_ue[i] = grad_scalar(uelist,klist)
		grd_ui[i] = grad_scalar(uilist,klist)
		grd_ute[i] = grad_scalar(utelist,klist)

	return grd_uem,grd_up,grd_ue,grd_ui,grd_ute



#//////////////////////////////////////////////////////////////////////////////////////////
	# Pressure Strains
def PressureStrains():
	pos_names = [] # Names of positions 
	posits = [] # get_data for each position
	# Get field and mec (position) data
	for i in range(4):
		pos_names.append('mms'+ str(probe[i]) + '_mec_r_gse')
		tinterpol(pos_names[i],interpfld,newname = 'pos' + str(i+1)) #interpolate
		posits.append(get_data('pos' + str(i+1)))
	timeax = posits[0].times
	ndata = len(timeax)
	# Reform data into np arrays in SI units (flds and posits)
	for i in range(4):
		posits[i] = 1e3*reform(posits[i]) # [pos1,pos2,pos3,pos4]
	pos1,pos2,pos3,pos4 = posits[0],posits[1],posits[2],posits[3]

	# Get kvectors (automatically in LMN)
	for i in range(ndata-1):
		klist = recip_vecs(pos1[i],pos2[i],pos3[i],pos4[i])
	
	
	# Field names variables
	velist,vilist,Pelist,Pilist = [],[],[],[]
	for i in range(4):
		vi_name = 'mms' + str(probe[i]) + '_' + 'dis' + '_bulkv_gse_brst'
		ve_name = 'mms' + str(probe[i]) + '_' + 'des' + '_bulkv_gse_brst'
		Pi_name = 'mms' + str(probe[i]) + '_dis_prestensor_gse_brst'
		Pe_name = 'mms' + str(probe[i]) + '_des_prestensor_gse_brst'
		var_name = interpfld
		tinterpol(vi_name,var_name, newname='vi')
		tinterpol(ve_name,var_name, newname='ve')
		tinterpol(Pi_name,var_name, newname='Pi')
		tinterpol(Pe_name,var_name, newname='Pe')
		vi = 1e3*reform(get_data('vi'))
		ve = 1e3*reform(get_data('ve'))
		Pi = 1e-9*reform(get_data('Pi'))
		Pe = 1e-9*reform(get_data('Pe'))
# This needs to be fixed.  The "convert lmn" needs to operate on a single time, not the whole thing at once
		for i in range(ndata-1):
			vi[i] = convert_lmn(vi[i],frame)
			ve[i] = convert_lmn(ve[i],frame)
			Pi[i] = convert_lmn(Pi[i],frame)
			Pe[i] = convert_lmn(Pe[i],frame)
		
		vilist.append(vi)
		velist.append(ve)
		Pilist.append(Pi)
		Pelist.append(Pe)


	e_strain = PStrain(Pilist,vilist,klist) 
	i_strain = PStrain(Pelist,velist,klist) 

	return i_strain, e_strain

#////////////////////////////////////////////////////////////////////////////////////

#%%
S1,He1,Hi1,Ke1,Ki1,qe1,qi1 = AllFluxes(1)  
S2,He2,Hi2,Ke2,Ki2,qe2,qi2 = AllFluxes(2)  
S3,He3,Hi3,Ke3,Ki3,qe3,qi3 = AllFluxes(3)  
S4,He4,Hi4,Ke4,Ki4,qe4,qi4 = AllFluxes(4) 
uem1,uE1,uB1,up1,ue1,ui1,ut1,ute1,uti1,uk1,uke1,uki1 = AllDensities(1)
uem2,uE2,uB2,up2,ue2,ui2,ut2,ute2,uti2,uk2,uke2,uki2 = AllDensities(2)
uem3,uE3,uB3,up3,ue3,ui3,ut3,ute3,uti3,uk3,uke3,uki3 = AllDensities(3)
uem4,uE4,uB4,up4,ue4,ui4,ut4,ute4,uti4,uk4,uke4,uki4 = AllDensities(4)
div_S,div_He,div_Hi,div_Ke,div_Ki,div_Qe,div_Qi = AllDivs()
grd_uem,grd_up,grd_ue,grd_ui,grd_ute = AllGrads()
#%%

ndata = len(S1)

pos_names = [] # Names of positions 
fld_interp_names = [] # interpolated version of field names 
posits = [] # get_data for each position

# Get field and mec (position) data
for i in range(4):
	pos_names.append('mms'+ str(probe[i]) + '_mec_r_gse')
	tinterpol(pos_names[i],interpfld,newname = 'pos' + str(i+1)) #interpolate
	posits.append(get_data('pos' + str(i+1)))

# N data points and time axis (setting to different vars here bc flds will change form)
# Also define shape of curl vs divergence (vector vs scalar)
timeax = posits[0].times
ndata = len(timeax)
crl = np.zeros([ndata,3])
divr = np.zeros([ndata])

Slist = [S1,S2,S3,S4]
Hilist = [Hi1,Hi2,Hi3,Hi4]
Helist = [He1,He2,He3,He4]
Kilist = [Ki1,Ki2,Ki3,Ki4]
Kelist = [Ke1,Ke2,Ke3,Ke4]
qelist = [qe1,qe2,qe3,qe4]
qilist = [qi1,qi2,qi3,qi4]
names = []


# Conversion factor outside everything
cnvrt = 1e9


for i in range(4):
	store_data('S'+str(i+1), data = {'x':timeax, 'y': cnvrt*Slist[i]})
	store_data('Ke'+str(i+1), data = {'x':timeax, 'y': cnvrt*Kelist[i]})
	store_data('He'+str(i+1), data = {'x':timeax, 'y': cnvrt*Helist[i]})
	store_data('Ki'+str(i+1), data = {'x':timeax, 'y': cnvrt*Kilist[i]})
	store_data('Hi'+str(i+1), data = {'x':timeax, 'y': cnvrt*Hilist[i]})
	store_data('qi'+str(i+1), data = {'x':timeax, 'y': cnvrt*qilist[i]})
	store_data('qe'+str(i+1), data = {'x':timeax, 'y': cnvrt*qelist[i]})
	names = names + ['S'+str(i+1), 'Ke'+str(i+1),'He'+str(i+1),'Ki'+str(i+1),'Hi'+str(i+1), 'qi'+str(i+1),'qe'+str(i+1)]

uemlist = [uem1,uem2,uem3,uem4]
uplist = [up1,up2,up3,up4]
uelist = [ue1,ue2,ue3,ue4]
uilist = [ui1,ui2,ui3,ui4]
utelist = [ute1,ute2,ute3,ute4]
utilist = [uti1,uti2,uti3,uti4]
store_data('uem', data = {'x':timeax, 'y': 0.25*cnvrt*(uemlist[0]+uemlist[1]+uemlist[2]+uemlist[3])})
store_data('up', data = {'x':timeax, 'y': 0.25*cnvrt*(uplist[0]+uplist[1]+uplist[2]+uplist[3])})
store_data('ue', data = {'x':timeax, 'y': 0.25*cnvrt*(uelist[0]+uelist[1]+uelist[2]+uelist[3])})
store_data('ui', data = {'x':timeax, 'y': 0.25*cnvrt*(uilist[0]+uilist[1]+uilist[2]+uilist[3])})
store_data('ute', data = {'x':timeax, 'y': 0.25*cnvrt*(utelist[0]+utelist[1]+utelist[2]+utelist[3])})
store_data('uti', data = {'x':timeax, 'y': 0.25*cnvrt*(utilist[0]+utilist[1]+utilist[2]+utilist[3])})
store_data('uke', data = {'x':timeax, 'y': 0.25*cnvrt*(uke1+uke2+uke3+uke4)})



store_data('uke4', data = {'x':timeax, 'y': cnvrt*uke4,'dy': 0.1*cnvrt*uke4})
store_data('ute4', data = {'x':timeax, 'y': cnvrt*ute4})
store_data('uem4', data = {'x':timeax, 'y': cnvrt*uem4})
store_data('uboth4',data = {'x':timeax, 'y': cnvrt*(uti4+ute4+uem4)})



# Specific variables for MMS2 (crossed 12/09 event center) ############################
store_data('uke2', data = {'x':timeax, 'y': cnvrt*uke2})
store_data('ute2', data = {'x':timeax, 'y': cnvrt*ute2})
store_data('uem2', data = {'x':timeax, 'y': cnvrt*uem2})
store_data('uboth2',data = {'x':timeax, 'y': cnvrt*(ute2+uem2)})
unames = ['uem','up','ue','ui','ute','uke','uti']

u_comp = np.zeros([ndata,4])

for i in range(len(u_comp)):
	u_comp[i] = np.array([uem2[i]+2*ute2[i]+2*uti2[i],uem2[i]+2*ute2[i],uem2[i],2*ute2[i]])
	#u_comp[i] = np.array([uem2[i]+ute2[i],uem2[i]+2*ute2[i],uem2[i]+2*ute2[i]/3,uem2[i],ute2[i],2*ute2[i],2*ute2[i]/3])
store_data('u_compare', data = {'x':timeax,'y':cnvrt*u_comp})
options('u_compare','thick',1.5)
options('u_compare','yrange',[0,10])
options('u_compare','color',['k','k','b','g'])


Pe_name = 'mms2_des_prestensor_gse_brst'
Te_name = 'mms2_des_temptensor_gse_brst'
Ti_name = 'mms2_dis_temptensor_gse_brst'
var_name = interpfld
tinterpol(Pe_name,var_name, newname='Pe')
tinterpol(Te_name,var_name, newname='Te')
tinterpol(Ti_name,var_name, newname='Ti')
Pe2 = reform(get_data('Pe'))
Te2 = reform(get_data('Te'))
Ti2 = reform(get_data('Ti'))
Pe2_ = np.zeros([ndata,3])
Te2_ = np.zeros([ndata,3])
Ti2_ = np.zeros([ndata,3])
Pe2_mag = np.zeros(ndata)
Te2_mag = np.zeros(ndata)
Ti2_mag = np.zeros(ndata)

# for i in range(ndata-1):
# 	for j in range(3):
# 		Pe2_[i,j] = Pe2[i,j,j]
# 	Pe2[i] = convert_lmn(Pe2_[i],frame) 

for i in range(ndata-1):
	for j in range(3):
		Pe2_[i,j] = Pe2[i,j,j]
		Te2_[i,j] = Te2[i,j,j]
		Ti2_[i,j] = Ti2[i,j,j]
	Pe2_mag[i] = Pe2_[i,0]+Pe2_[i,1]+Pe2_[i,2]
	Te2_mag[i] = Te2_[i,0]+Te2_[i,1]+Te2_[i,2]
	Ti2_mag[i] = Ti2_[i,0]+Ti2_[i,1]+Ti2_[i,2]
	Pe2_[i] = convert_lmn(Pe2_[i],frame) 
	Te2_[i] = convert_lmn(Te2_[i],frame)
	Ti2_[i] = convert_lmn(Ti2_[i],frame)



store_data('Pe2', data = {'x':timeax, 'y':Pe2_})
store_data('Te2', data = {'x':timeax, 'y':Te2_})
store_data('Ti2', data = {'x':timeax, 'y':Ti2_})

store_data('Pe2_mag', data = {'x':timeax, 'y':Pe2_mag/3})
store_data('Te2_mag', data = {'x':timeax, 'y':Te2_mag/3})
options(['Pe2','Te2','Pe2_mag','Te2_mag','mms2_des_numberdensity_brst','Ti2'],'thick',1.5)
options('Pe2','yrange',[0,0.5])
options('Te2','yrange',[0,160])
options('Ti2','yrange',[-1000,1000])
options('mms2_des_numberdensity_brst','yrange',[0,30])
options(['Pe2','Te2'],'color',['b','g','r'])
# tplot(['Te2','mms2_des_numberdensity_brst'])
#tplot(['mms2_des_numberdensity_brst','Te2_mag', 'Pe2_mag'])
###################################
# tplot(['Te2','Ti2'],'yrange',[0,100])
# tplot(['Pe2','Te2'])
#tplot(['Te2'])
#tplot(['mms2_des_numberdensity_brst'])
#%%
cnvrt = 1e9
store_data('divS', data = {'x':timeax, 'y': cnvrt*div_S})     # ALL CONVERTED TO nW/m^3
store_data('divHi', data = {'x':timeax, 'y': cnvrt*div_Hi})
store_data('divHe', data = {'x':timeax, 'y': cnvrt*div_He})
store_data('divKi', data = {'x':timeax, 'y': cnvrt*div_Ki})
store_data('divKe', data = {'x':timeax, 'y': cnvrt*div_Ke})
store_data('divQi', data = {'x':timeax, 'y': cnvrt*div_Qi})
store_data('divQe', data = {'x':timeax, 'y': cnvrt*div_Qe})
divnames = ['divS','divHi','divHe','divKi','divKe','divQi','divQe']

store_data('grd_uem', data = {'x':timeax,'y': cnvrt*grd_uem})
store_data('grd_up', data = {'x':timeax,'y': cnvrt*grd_up})
store_data('grd_ue', data = {'x':timeax,'y': cnvrt*grd_ue})
store_data('grd_ui', data = {'x':timeax,'y': cnvrt*grd_ui})
store_data('grd_ute', data = {'x':timeax,'y': cnvrt*grd_ute})
grdnames = ['grd_uem','grd_up','grd_ue','grd_ui','grd_ute']







#%%
# Getting  Curl of B 
# Get all Bs and positions in np form and interpolted together
# [0,1,2,3] = MMS[1,2,3,4]
B_names,E_names,Pe_names,Pi_names,Te_names,Ti_names,ne_names, ni_names, vi_names,ve_names,ni_names,ne_names,pos_names,\
Binterp_names, viinterp_names,veinterp_names,niinterp_names,neinterp_names,\
Bflds, Eflds, viflds,veflds, neflds, niflds, eTemps,iTemps,ePres,iPres,posits\
= [],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]
#%%
# Factors to convert fld & position to SI units (change fldfact depending on fld, keep posfact as is)
Bfact = 1e-9  
Pfact = 1e-9  
vfact = 1e3
Efact = 1e-3 
nfact = 1e6  
Tfact = 1.602e-19
posfact = 1e3 

probes = probe
# Get field and mec (position) data
for i in range(4):
	B_names.append('mms' + str(probes[i]) + '_fgm_b_gse_brst_l2') #Change this if you want to calc for another vec quantity  
	pos_names.append('mms'+ str(probes[i]) + '_mec_r_gsm')
	ve_names.append('mms' + str(probes[i]) + '_des_bulkv_gse_brst')
	vi_names.append('mms' + str(probes[i]) + '_dis_bulkv_gse_brst')
	E_names.append('mms' + str(probes[i]) + '_edp_dce_gse_brst_l2')
	ne_names.append('mms' + str(probes[i]) + '_des_numberdensity_brst')
	ni_names.append('mms' + str(probes[i]) + '_dis_numberdensity_brst')
	Te_names.append('mms' + str(probes[i]) + '_des_temptensor_gse_brst')
	Ti_names.append('mms' + str(probes[i]) + '_dis_temptensor_gse_brst')
	Pe_names.append('mms' + str(probes[i]) + '_des_prestensor_gse_brst')
	Pi_names.append('mms' + str(probes[i]) + '_dis_prestensor_gse_brst')
	tinterpol(B_names[i],interpfld,newname = 'B' + str(i+1))
	tinterpol(E_names[i],interpfld,newname = 'E' + str(i+1))
	tinterpol(ve_names[i],interpfld,newname = 've' + str(i+1))
	tinterpol(vi_names[i],interpfld,newname = 'vi' + str(i+1))
	tinterpol(ne_names[i],interpfld,newname = 'ne' + str(i+1))
	tinterpol(ni_names[i],interpfld,newname = 'ni' + str(i+1))
	tinterpol(pos_names[i],interpfld,newname = 'pos' + str(i+1))
	tinterpol(Te_names[i],interpfld,newname = 'Te' + str(i+1))
	tinterpol(Ti_names[i],interpfld,newname = 'Ti' + str(i+1))
	tinterpol(Pe_names[i],interpfld,newname = 'Pe' + str(i+1))
	tinterpol(Pi_names[i],interpfld,newname = 'Pi' + str(i+1))
	Bflds.append(get_data('B' + str(i+1)))
	Eflds.append(get_data('E' + str(i+1)))
	veflds.append(get_data('ve' + str(i+1)))
	viflds.append(get_data('vi' + str(i+1)))
	niflds.append(get_data('ni' + str(i+1)))
	neflds.append(get_data('ne' + str(i+1)))
	posits.append(get_data('pos' + str(i+1)))
	eTemps.append(get_data('Te' + str(i+1)))
	iTemps.append(get_data('Ti' + str(i+1)))
	ePres.append(get_data('Pe' + str(i+1)))
	iPres.append(get_data('Pi' + str(i+1)))
#%%
# N data points and time axis (setting to different vars here bc flds will change form)
# Also define shape of curl vs divergence (vector vs scalar)
# a = get_data(interpfld)
# timeax = a.times
# ndata = len(timeax)
# crl = np.zeros([ndata,3])
timeax = Bflds[0].times
ndata = len(timeax)
crl = np.zeros([ndata,3])
divJ = np.zeros([ndata,4])
divE = np.zeros([ndata,4])
div_nve = np.zeros([ndata,4])
ndiv_ve = np.zeros_like(div_nve)
v_dot_gradne = np.zeros([ndata])
v_dot_gradni = np.zeros([ndata])
v_dot_gradTe = np.zeros([ndata])
v_dot_gradTi = np.zeros([ndata])
v_dot_graduem = np.zeros([ndata])
v_dot_gradute = np.zeros([ndata])
v_dot_gradBLsqd = np.zeros([ndata])
v_dot_gradBMsqd = np.zeros([ndata])
v_dot_gradBNsqd = np.zeros([ndata])
#%%
# Reform data into np arrays in SI units (flds and posits)
for i in range(4):
    Bflds[i] = Bfact*reform(Bflds[i]) # [fld1,fld2,fld3,fld4]
    Eflds[i] = Efact*reform(Eflds[i]) # [fld1,fld2,fld3,fld4]
    veflds[i] = vfact*reform(veflds[i]) # [fld1,fld2,fld3,fld4]
    viflds[i] = vfact*reform(viflds[i]) # [fld1,fld2,fld3,fld4]
    niflds[i] = nfact*reform(niflds[i]) # [fld1,fld2,fld3,fld4]
    neflds[i] = nfact*reform(neflds[i]) # [fld1,fld2,fld3,fld4]
    posits[i] = posfact*reform(posits[i]) # [pos1,pos2,pos3,pos4]
    eTemps[i] = Tfact*reform(eTemps[i])
    iTemps[i] = Tfact*reform(iTemps[i])
    ePres[i] = Pfact*reform(ePres[i])
    iPres[i] = Pfact*reform(iPres[i])

#%%


# Make a set of E variables in LMN
E1 = np.zeros_like(Eflds[0])
E2 = np.zeros_like(Eflds[0])
E3 = np.zeros_like(Eflds[0])
E4 = np.zeros_like(Eflds[0])
for i in range(len(Eflds[0])):
	E1[i] = convert_lmn(Eflds[0][i],frame)
	E2[i] = convert_lmn(Eflds[1][i],frame)
	E3[i] = convert_lmn(Eflds[2][i],frame)
	E4[i] = convert_lmn(Eflds[3][i],frame)
################################################

ndif1 = np.zeros_like(niflds[0])
ndif2 = np.zeros_like(niflds[0])
ndif3 = np.zeros_like(niflds[0])
ndif4 = np.zeros_like(niflds[0])
ndif_avg = np.zeros_like(niflds[0])
ne_avg = np.zeros_like(niflds[0])
ni_avg = np.zeros_like(niflds[0])
for i in range(len(neflds[0])):
	ndif1[i] = niflds[0][i] - neflds[0][i]
	ndif2[i] = niflds[1][i] - neflds[1][i]
	ndif3[i] = niflds[2][i] - neflds[2][i]
	ndif4[i] = niflds[3][i] - neflds[3][i]
	ndif_avg[i] = (ndif1[i] + ndif2[i] + ndif3[i] + ndif4[i])/4
	ne_avg[i] = (neflds[0][i] + neflds[1][i] + neflds[2][i] + neflds[3][i])/4 
	ni_avg[i] = (niflds[0][i] + niflds[1][i] + niflds[2][i] + niflds[3][i])/4 
#%%
Te1 = np.zeros_like(ndif1)
Te2 = np.zeros_like(ndif1)
Te3 = np.zeros_like(ndif1)
Te4 = np.zeros_like(ndif1)
Ti1 = np.zeros_like(ndif1)
Ti2 = np.zeros_like(ndif1)
Ti3 = np.zeros_like(ndif1)
Ti4 = np.zeros_like(ndif1)
#%%
for i in range(ndata-1):
	Te1[i] = np.trace(eTemps[0][i])
	Te2[i] = np.trace(eTemps[1][i])
	Te3[i] = np.trace(eTemps[2][i])
	Te4[i] = np.trace(eTemps[3][i])

	Ti1[i] = np.trace(iTemps[0][i])
	Ti2[i] = np.trace(iTemps[1][i])
	Ti3[i] = np.trace(iTemps[2][i])
	Ti4[i] = np.trace(iTemps[3][i])
# difflist = [ndif1,ndif2,ndif3,ndif4]

#%%
# If fld is 4D (usually because total is included), chop off the 4th term
if len(Bflds[0][0]) == 4:   # Just put this here so its convenient to minimize in vscode
    Bfld1 = np.zeros([ndata,3])
    Bfld2 = np.zeros([ndata,3])
    Bfld3 = np.zeros([ndata,3])
    Bfld4 = np.zeros([ndata,3])
    for i in range(ndata-1):
        Bfld1[i] = Bflds[0][i][:-1]
        Bfld2[i] = Bflds[1][i][:-1] 
        Bfld3[i] = Bflds[2][i][:-1]
        Bfld4[i] = Bflds[3][i][:-1]
    Bflds = [Bfld1,Bfld2,Bfld3,Bfld4]

# Make a set of B variables in LMN
B1 = np.zeros_like(Eflds[0])
B2 = np.zeros_like(Bflds[0])
B3 = np.zeros_like(Bflds[0])
B4 = np.zeros_like(Bflds[0])
for i in range(len(Bflds[0])):
	B1[i] = convert_lmn(Bflds[0][i],frame)
	B2[i] = convert_lmn(Bflds[1][i],frame)
	B3[i] = convert_lmn(Bflds[2][i],frame)
	B4[i] = convert_lmn(Bflds[3][i],frame)
################################################

# Calculate Grad(n)


#Calculate Div(V)



# "moments" j, assuming ne~ni  (Converted to LMN)
j_moms = []
for sc in range(4):
	jm = np.zeros_like(Eflds[0])
	for i in range(ndata-1):
		# jm[i] = get_j(neflds[sc][i],viflds[sc][i],veflds[sc][i])
		jm[i] = get_j(neflds[sc][i],viflds[sc][i],veflds[sc][i])
	j_moms.append(jm)

## Get Curlometer J ##
for i in range(ndata-1):
    veclist = [B1[i],B2[i],B3[i],B4[i]]
    klist = recip_vecs(posits[0][i],posits[1][i],posits[2][i],posits[3][i])
    crl[i] = curl(veclist,klist)
j_crl = crl/mu0


###########

















vxl = np.zeros(3) # Default to zero
if str(frame) == str(lmn_1209):
	vxl = np.array([-42.98e3,246.57e3,-25.56e3])
	#vxl = np.array([0,0,0]) # version to play around with


# Use x-line frame E
Epflds = []
Edrifts = []
for sc in range(4):
	Ep = np.zeros_like(Eflds[0])
	for i in range(ndata-1):
		Ep[i] = Eprime(Eflds[sc][i],vxl,Bflds[sc][i])
		# Edrift[i] = Eprime(Eflds[sc][i],vdrift[i],Bflds[sc][i]])
	Epflds.append(Ep)
	# Edrifts.append(Edrift)



   
# Make a set of Ep variables in LMN
Ep1 = np.zeros_like(Eflds[0])
Ep2 = np.zeros_like(Eflds[0])
Ep3 = np.zeros_like(Eflds[0])
Ep4 = np.zeros_like(Eflds[0])
for i in range(len(Eflds[0])):
	Ep1[i] = convert_lmn(Epflds[0][i],frame)
	Ep2[i] = convert_lmn(Epflds[1][i],frame)
	Ep3[i] = convert_lmn(Epflds[2][i],frame)
	Ep4[i] = convert_lmn(Epflds[3][i],frame)

Epflds = [Ep1,Ep2,Ep3,Ep4] # Now Ep fields are in LMN
################################################


# Make a set of J variables in LMN
j1 = np.zeros_like(j_moms[0])
j2 = np.zeros_like(j_moms[0])
j3 = np.zeros_like(j_moms[0])
j4 = np.zeros_like(j_moms[0])
for i in range(len(j_moms[0])):
	j1[i] = convert_lmn(j_moms[0][i],frame)
	j2[i] = convert_lmn(j_moms[1][i],frame)
	j3[i] = convert_lmn(j_moms[2][i],frame)
	j4[i] = convert_lmn(j_moms[3][i],frame)
################################################



nvs,vs = [],[]
for sc in range(4):
	nv = np.zeros_like(Eflds[0])
	v = np.zeros_like(Eflds[0])
	for i in range(ndata-1):
		# jm[i] = get_j(neflds[sc][i],viflds[sc][i],veflds[sc][i])
		v[i] = veflds[sc][i]
		nv[i] = neflds[sc][i]*v[i]
	nvs.append(nv)
	vs.append(v)


nv1 = np.zeros_like(j_moms[0])
nv2 = np.zeros_like(j_moms[0])
nv3 = np.zeros_like(j_moms[0])
nv4 = np.zeros_like(j_moms[0])

v1 = np.zeros_like(j_moms[0])
v2 = np.zeros_like(j_moms[0])
v3 = np.zeros_like(j_moms[0])
v4 = np.zeros_like(j_moms[0])
v_avg = np.zeros_like(j_moms[0])
for i in range(len(j_moms[0])):
	nv1[i] = convert_lmn(nvs[0][i],frame)
	nv2[i] = convert_lmn(nvs[1][i],frame)
	nv3[i] = convert_lmn(nvs[2][i],frame)
	nv4[i] = convert_lmn(nvs[3][i],frame)
	v1[i] = convert_lmn(vs[0][i],frame)
	v2[i] = convert_lmn(vs[1][i],frame)
	v3[i] = convert_lmn(vs[2][i],frame)
	v4[i] = convert_lmn(vs[3][i],frame)
	v_avg[i] = (v1[i] + v2[i] + v3[i] + v4[i])/4
################################################

# DivJ & div(nv)#
for i in range(ndata-1):
	jlist = [j1[i],j2[i],j3[i],j4[i]] 
	nvlist = [nv1[i],nv2[i],nv3[i],nv4[i]]
	vlist = [v1[i],v2[i],v3[i],v4[i]]
	nelist = [neflds[0][i],neflds[1][i],neflds[2][i],neflds[3][i]]
	nilist = [niflds[0][i],niflds[1][i],niflds[2][i],niflds[3][i]]
	Telist = [Te1[i],Te2[i],Te3[i],Te4[i]]
	BLsqdlist = [B1[i][0]**2,B2[i][0]**2,B3[i][0]**2,B4[i][0]**2]
	BMsqdlist = [B1[i][1]**2,B2[i][1]**2,B3[i][1]**2,B4[i][1]**2]
	BNsqdlist = [B1[i][2]**2,B2[i][2]**2,B3[i][2]**2,B4[i][2]**2]
	klist = recip_vecs(posits[0][i],posits[1][i],posits[2][i],posits[3][i])
	divJ[i] = div(jlist,klist) 
	div_nve[i] = div(nvlist,klist)
	ndiv_ve[i] = ne_avg[i]*div(vlist,klist)
	# v_dot_gradn[i] = np.dot(v_avg[i],grad_scalar(nlist,klist))
	v_dot_gradne[i] = np.dot(vxl,grad_scalar(nelist,klist))
	v_dot_gradni[i] = np.dot(vxl,grad_scalar(nilist,klist))
	v_dot_gradTe[i] = np.dot(vxl,grad_scalar(Telist,klist))
	v_dot_graduem[i] = np.dot(vxl,grd_uem[i])
	v_dot_gradute[i] = np.dot(vxl,grd_ute[i])

	v_dot_gradBLsqd[i] = np.dot(vxl,grad_scalar(BLsqdlist,klist))/(2*mu0)
	v_dot_gradBMsqd[i] = np.dot(vxl,grad_scalar(BMsqdlist,klist))/(2*mu0)
	v_dot_gradBNsqd[i] = np.dot(vxl,grad_scalar(BNsqdlist,klist))/(2*mu0)

Hdrift1 = np.zeros_like(He1)
Hdrift2 = np.zeros_like(He2)
Hdrift3 = np.zeros_like(He3)
Hdrift4 = np.zeros_like(He4)

for i in range(ndata-1):
	Hdrift1[i] = enth_flux(np.cross(Ep1[i],B1[i])/(np.linalg.norm(B1[i])**2),ePres[0][i])
	Hdrift2[i] = enth_flux(np.cross(Ep2[i],B2[i])/(np.linalg.norm(B2[i])**2),ePres[1][i])
	Hdrift3[i] = enth_flux(np.cross(Ep3[i],B3[i])/(np.linalg.norm(B3[i])**2),ePres[2][i])
	Hdrift4[i] = enth_flux(np.cross(Ep4[i],B4[i])/(np.linalg.norm(B4[i])**2),ePres[3][i])

store_data('Hd1',data = {'x':timeax,'y':1e3*Hdrift1})
store_data('Hd2',data = {'x':timeax,'y':1e3*Hdrift2})
store_data('Hd3',data = {'x':timeax,'y':1e3*Hdrift3})
store_data('Hd4',data = {'x':timeax,'y':1e3*Hdrift4})

store_data('Hdiff', data = {'x':timeax,'y':1e3*(He2 - Hdrift2)})

# store_data('Te2',data = {'x':timeax,'y':Te2})
# store_data('ne2',data = {'x':timeax,'y':nelist[1]})

# tplot(['Te2','ne2'])


#%%
###### This is specific to the 0.6 sec interval #######
# Make a manual d/dt here I can see that it will basically balance out to zero 
ptspersec = ndata/0.6

uem_avg = np.zeros_like(uem1)
ute_avg = np.zeros_like(uem1)
Te_avg = np.zeros_like(uem1)
Ti_avg = np.zeros_like(uem1)
BLsqd_avg = np.zeros_like(uem1)
BMsqd_avg = np.zeros_like(uem1)
BNsqd_avg = np.zeros_like(uem1)
#///////////////////////////////////////////////////////////////////////
for i in range(ndata-1):
	uem_avg[i] = 0.25*(uem1[i] + uem2[i] + uem3[i] + uem4[i]) 
	ute_avg[i] = 0.25*(ute1[i] + ute2[i] + ute3[i] + ute4[i])
	Te_avg[i] = 0.25*(Te1[i] + Te2[i] + Te3[i] + Te4[i]) 
	Ti_avg[i] = 0.25*(Ti1[i] + Ti2[i] + Ti3[i] + Ti4[i])
	BLsqd_avg[i] = 0.25*(B1[i][0]**2 + B2[i][0]**2 +B3[i][0]**2 +B4[i][0]**2)/(2*mu0) 
	BMsqd_avg[i] = 0.25*(B1[i][1]**2 + B2[i][1]**2 +B3[i][1]**2 +B4[i][1]**2)/(2*mu0)
	BNsqd_avg[i] = 0.25*(B1[i][2]**2 + B2[i][2]**2 +B3[i][2]**2 +B4[i][2]**2)/(2*mu0) 

	# BLsqd_avg[i] = 0.25*(B_dens(B1[i,0])+B_dens(B2[i,0]) +B_dens(B3[i,0]) +B_dens(B4[i,0])) 
	# BMsqd_avg[i] = 0.25*(B_dens(B1[i,1])+B_dens(B2[i,1]) +B_dens(B3[i,1]) +B_dens(B4[i,1]))
	# BNsqd_avg[i] = 0.25*(B_dens(B1[i,2])+B_dens(B2[i,2]) +B_dens(B3[i,2]) +B_dens(B4[i,2])) 

derivt_uem = np.zeros(ndata)
derivt_ute = np.zeros(ndata)
derivt_ne = np.zeros(ndata)
derivt_Te = np.zeros(ndata)
derivt_ni = np.zeros(ndata)
derivt_Ti = np.zeros(ndata)
derivt_BLsqd = np.zeros(ndata)
derivt_BMsqd = np.zeros(ndata)
derivt_BNsqd = np.zeros(ndata)
for i in range(1,ndata-1):
	derivt_uem[i] = 0.5*ptspersec*(uem_avg[i+1] - uem_avg[i-1])
	derivt_ute[i] = 0.5*ptspersec*(ute_avg[i+1] - ute_avg[i-1])
	derivt_ne[i] = 0.5*ptspersec*(ne_avg[i+1] - ne_avg[i-1])
	derivt_Te[i] = 0.5*ptspersec*(Te_avg[i+1] - Te_avg[i-1])
	derivt_ni[i] = 0.5*ptspersec*(ni_avg[i+1] - ni_avg[i-1])
	derivt_Ti[i] = 0.5*ptspersec*(Ti_avg[i+1] - Ti_avg[i-1])
	derivt_BLsqd[i] = 0.5*ptspersec*(BLsqd_avg[i+1] - BLsqd_avg[i-1])
	derivt_BMsqd[i] = 0.5*ptspersec*(BMsqd_avg[i+1] - BMsqd_avg[i-1])
	derivt_BNsqd[i] = 0.5*ptspersec*(BNsqd_avg[i+1] - BNsqd_avg[i-1])

store_data('Convec_duemdt',data = {'x':timeax,'y':1e9*(derivt_uem - v_dot_graduem)})
store_data('Convec_dutedt',data = {'x':timeax,'y':1e9*(derivt_ute - v_dot_gradute)})

convect = np.zeros([ndata,3])
for i in range(ndata):
	convect[i][0] = 0.5e9*Te_avg[i]*(0*derivt_ne[i] - v_dot_gradne[i])
	convect[i][1] = 0.5e9*ne_avg[i]*(0*derivt_Te[i] - v_dot_gradTe[i])
	convect[i][2] = 1e9*(0*derivt_ute[i] - v_dot_gradute[i])
	# convect[i][2] = convect[i][0] + convect[i][1]


Bsqdterms = np.zeros([ndata,4])
for i in range(ndata):
	Bsqdterms[i][0] = (derivt_BLsqd[i] - v_dot_gradBLsqd[i])
	Bsqdterms[i][1] = (derivt_BMsqd[i] - v_dot_gradBMsqd[i])
	# Bsqdterms[i][1] = 0
	Bsqdterms[i][2] = (derivt_BNsqd[i] - v_dot_gradBNsqd[i])
	Bsqdterms[i][3] = (derivt_uem[i] - v_dot_graduem[i])
	# Bsqdterms[i][3] = 0 
	# Bsqdterms[i][3] = Bsqdterms[i][0]+Bsqdterms[i][1]+Bsqdterms[i][2]


convectboth = np.zeros([ndata,3])
for i in range(ndata):
	convectboth[i][0] = 0.5*Ti_avg[i]*(derivt_ni[i] - v_dot_gradni[i]) + 0.5*ni_avg[i]*(derivt_Ti[i] - v_dot_gradTi[i])
	convectboth[i][1] = 0.5*Te_avg[i]*(derivt_ne[i] - v_dot_gradne[i]) + 0.5*ne_avg[i]*(derivt_Te[i] - v_dot_gradTe[i])
	convectboth[i][2] = convectboth[i][0] + convectboth[i][1]
	# convectboth[i][3] = -(derivt_uem[i] - v_dot_graduem[i])
	# convect[i][2] = 1e9*(derivt_ute[i] - v_dot_gradute[i])

dnboth = np.zeros([ndata,2])
for i in range(ndata):
	dnboth[i][0] = derivt_ni[i] - v_dot_gradni[i]
	dnboth[i][1] = derivt_ne[i] - v_dot_gradne[i]
# store_data('Convec_dnedt',data = {'x':timeax,'y':0.5e9*Ti_avg*(derivt_ni - v_dot_gradni)})

store_data('Convec_dndt',data = {'x':timeax,'y':1e-6*dnboth})
store_data('Convec_dnidt',data = {'x':timeax,'y':1e-6*(derivt_ni - v_dot_gradni)})
store_data('Convec_dnedt',data = {'x':timeax,'y':1e-6*(derivt_ne - v_dot_gradne)})
store_data('Convec_diff',data = {'x':timeax,'y':1e-6*(derivt_ne - v_dot_gradne - derivt_ni + v_dot_gradni)})

options('Convec_dndt', 'thick',1.5)
options('Convec_dndt', 'color',['r','b'])
options('Convec_dndt', 'yrange',[-20,100])
# store_data('Convec_dTedt',data = {'x':timeax,'y':0.5e9*ni_avg*(derivt_Ti - v_dot_gradTi)})
# store_data('JdE4',data = {'x':timeax,'y':1e9*(JdotE4)})


store_data('convect',data = {'x':timeax,'y':convect})
store_data('convec_both', data = {'x':timeax,'y':1e9*(convectboth)})
store_data('Convec_uemconts',data = {'x':timeax,'y':1e9*Bsqdterms} )
options('Convec_uemconts','thick',1.5)
options('Convec_uemconts','color',['b','g','r','k'])
options('Convec_uemconts','yrange',[-8,5])

options('convec_both', 'thick',1.5)
options('convec_both', 'color',['r','b','k'])
options('convect','thick',1.5)
options('convect','color',['b','r','k'])
options('convect','yrange',[-0.4,2])

# options(['Convec_duemdt','Convec_dutedt','Convec_dnedt','Convec_dTedt','Convec_dnidt','Convec_dTidt'],'thick',1.5)
# options(['Convec_duemdt','Convec_dutedt','Convec_dnedt','Convec_dTedt'],'yrange',[-2,2])
options(['JdotE_LMN','Convec_duemdt','Convec_dutedt'],'thick',1.5)
options('Convec_duemdt','yrange',[-8,5])
options('Convec_dutedt','yrange',[-0.4,2])
tplot(['convect'])
#/////////////////////////////////////////////////////////////////
#%%






#%%
# DivE #
for i in range(ndata-1):
    Elist = [E1[i],E2[i],E3[i],E4[i]]
    klist = recip_vecs(posits[0][i],posits[1][i],posits[2][i],posits[3][i])
    divE[i] = div(Elist,klist) 

# "moments" versions of J.E and J.E' for all SC 
jdem = np.zeros_like(neflds)
jdepm = np.zeros_like(neflds)
for sc in range(4):
	for i in range(ndata-1):
		jdem[sc][i] = np.dot(j_moms[sc][i],Eflds[sc][i]) #something is going wrong here
		jdepm[sc][i] = np.dot(j_moms[sc][i],Epflds[sc][i])   

JdotE_mom = sum(jdem)/4
# JdotEp_mom = sum(jdepm)/4
JdotEp_mom = sum(jdepm)/4 # For some reason, the other sc measure zero somewhere, so the average is diminished
# Keeping the only SC with nonzero results for now

# Make a set of J.E variables in LMN (with individual contributions in vector form)
je1 = np.zeros_like(j_moms[0])
je2 = np.zeros_like(j_moms[0])
je3 = np.zeros_like(j_moms[0])
je4 = np.zeros_like(j_moms[0])
for i in range(len(j_moms[0])):
	je1[i] = j1[i]*Ep1[i]
	je2[i] = j2[i]*Ep2[i]
	je3[i] = j3[i]*Ep3[i]
	je4[i] = j4[i]*Ep4[i]

############
# temporary fix for 8/10 event:
# np.concatenate(Eflds[1],Eflds[1][-2])
##########################

#E_avg = sum(Eflds)/4
Ep_avg = sum(Epflds)/4

#JdotE = np.zeros(len(E_avg))
JdotEp = np.zeros(len(Ep_avg))
JdotE_LMN = np.zeros([ndata,4])
JdE2 = np.zeros([ndata,4])
JdE4 = np.zeros([ndata,4])
ddt = np.zeros([ndata,3])
jdeavg = np.zeros(ndata)
for i in range(len(Ep_avg)):
	#JdotE[i] = np.dot(j_crl[i],Ep_avg[i])
	# JdotE_LMN[i] = np.array([j_crl[i]*Ep_avg[i],j_crl[i],Ep_avg[i]])
	JdotE_LMN[i] = np.array([j_crl[i][0]*Ep_avg[i][0],j_crl[i][1]*Ep_avg[i][1],j_crl[i][2]*Ep_avg[i][2],np.dot(j_crl[i],Ep_avg[i])])
	JdE2[i] = np.array([je2[i,0],je2[i,1],je2[i,2],je2[i,0]+je2[i,1]+je2[i,2]])
	JdE4[i] = np.array([je4[i,0],je4[i,1],je4[i,2],je4[i,0]+je4[i,1]+je4[i,2]])
	ddt[i] = np.array([-JdotE_LMN[i,3]-div_S[i,3],JdotE_LMN[i,3]-div_He[i,3]-div_Ke[i,3]-div_Qe[i,3],JdotE_LMN[i,3]-div_Hi[i,3]-div_Ki[i,3]-div_Qi[i,3]])
	#JdotEp[i] = np.dot(j_crl[i],Ep_avg[i]) # This version uses J_curl dot average Eprime
# timespan('2016-12-09 09:03:54', 1, keyword='seconds')

# store_data('Poynt', data = {'x':timeax,'y':1e9*np.array([JdotE,div_S,-JdotE-div_S])})
store_data('jdote2', data = {'x':timeax,'y':1e9*JdE2[:,3]})  # CONVERTED TO nW/m^3
store_data('jdote4', data = {'x':timeax,'y':1e9*JdE4[:,3]})  # CONVERTED TO nW/m^3
store_data('jdotep_mom', data = {'x':timeax,'y':cnvrt*JdotEp_mom})  # CONVERTED TO nW/m^3
store_data('jdote_mom', data = {'x':timeax,'y':cnvrt*JdotE_mom})  # CONVERTED TO nW/m^3


store_data('du/dt', data = {'x':timeax,'y': -cnvrt*(JdotE_LMN[:,3]+div_S[:,3])})
store_data('dup/dt', data = {'x':timeax,'y': -cnvrt*(-JdotE_LMN[:,3]+div_He[:,3]+div_Ke[:,3]+div_Hi[:,3] + div_Ki[:,3])}) #still need to include Qe

store_data('d/dt', data = {'x':timeax,'y': cnvrt*ddt}) #still need to include Qe

##STill playing around with these variables ##############
store_data('dne/dt', data = {'x':timeax,'y': -(1e-6)*div_nve[:,3]}) # Units 1/scm^3 s
store_data('divEeps0', data = {'x':timeax,'y': 1e-6*eps0*divE}) # Units 1/cm^3 s
options('divEeps0', 'colors',['b','g','r','k'])
# store_data('ne_fromdiv', data = {'x':timeax,'y': 22-1e-6*eps0*divE[:,0]/(e)}) # Units 1/cm^3 s
store_data('ndif_avg',data = {'x':timeax,'y': 1e-6*ndif_avg})
store_data('rhodif_avg',data = {'x':timeax,'y': 1e-6*ndif_avg*e}) #C/cm^3
store_data('ndivV',data = {'x':timeax,'y': 1e-6*ndiv_ve[:,3]}) #C/cm^3
store_data('vgradn',data = {'x':timeax,'y': 1e-6*v_dot_gradne}) #C/cm^3
store_data('vgraduem',data = {'x':timeax,'y': 1e9*v_dot_graduem}) #nJ/m^3
store_data('vgradute',data = {'x':timeax,'y': 1e9*v_dot_gradute}) #nJ/m^3


# tplot(['vgradute','vgraduem'])
#####################################################################################


store_data('je1', data = {'x':timeax,'y':cnvrt*(je1)})
store_data('je2', data = {'x':timeax,'y':cnvrt*(je2)})
store_data('je3', data = {'x':timeax,'y':cnvrt*(je3)})
store_data('je4', data = {'x':timeax,'y':cnvrt*(je4)})
store_data('JdotE_LMN', data = {'x':timeax,'y':cnvrt*(JdotE_LMN)})

store_data('B2_LMN',data = {'x':timeax,'y':cnvrt*B2})
store_data('B1_LMN',data = {'x':timeax,'y':cnvrt*B1})
store_data('B3_LMN',data = {'x':timeax,'y':cnvrt*B3})
store_data('B4_LMN',data = {'x':timeax,'y':cnvrt*B4})
Bterms = ['B1_LMN','B2_LMN','B3_LMN','B4_LMN']

store_data('E2_LMN',data = {'x':timeax,'y':cnvrt*Ep2})
store_data('E1_LMN',data = {'x':timeax,'y':cnvrt*Ep1})
store_data('E3_LMN',data = {'x':timeax,'y':cnvrt*Ep3})
store_data('E4_LMN',data = {'x':timeax,'y':cnvrt*Ep4})
Eterms = ['E1_LMN','E2_LMN','E3_LMN','E4_LMN']

ENBM = Ep2[:,2]*B2[:,1]
EMBN = Ep2[:,1]*B2[:,2]
ELBM = Ep2[:,0]*B2[:,1]
EMBL = Ep2[:,1]*B2[:,0]

SL = np.zeros([ndata,2])
SN = np.zeros([ndata,2])

for i in range(ndata):
	SL[i,0] = EMBN[i]
	SL[i,1] = -ENBM[i]
	SN[i,0] = ELBM[i]
	SN[i,1] = -EMBL[i]

store_data('SL2_comps', data = {'x':timeax,'y': 1e3*SL/mu0 })
store_data('SN2_comps', data = {'x':timeax,'y': 1e3*SN/mu0 })
options('SL2_comps', 'color', ['r','g'])
options('SN2_comps', 'color', ['g','r'])

options(['SL2_comps','SN2_comps'], 'thick', 1.5)
options(['SL2_comps','SN2_comps'], 'yrange', [-1,1])

tplot(['SL2_comps','SN2_comps'])
options(Eterms,'thick',1.5)
options(Bterms,'thick',1.5)
# tplot(['jdote','jdotep','jdotep_mom'])
#%%
jeterms = ['je1','je2','je3','je4']
Sterms = ['S1','S2','S3','S4']
Heterms = ['He1','He2','He3','He4']
Hiterms = ['Hi1','Hi2','Hi3','Hi4']
Keterms = ['Ke1','Ke2','Ke3','Ke4']
Kiterms = ['Ki1','Ki2','Ki3','Ki4']
qiterms = ['qi1','qi2','qi3','qi4']
qeterms = ['qe1','qe2','qe3','qe4']
options(qiterms,'thick',1.5)
options(qeterms,'thick',1.5)
options(jeterms, 'thick',1.5)
options(Sterms, 'thick',1.5)
options(Heterms, 'thick',1.5)
options(Hiterms, 'thick',1.5)
options(Keterms, 'thick',1.5)
options(Kiterms, 'thick',1.5)
pterms = ['divS','jdote','du/dt','Poynt Terms']
terms = ['divS','jdote','du/dt','divHe','divKe']
options(pterms,'thick', 1.5)

divterms = ['divS','divHe','divKe','divHi','divKi','divQe','JdotE_LMN']
#options(divterms,'Color',['b','g','r','k'])
#options(divterms,'legend_names', ['L','M','N','Total'])
options(divterms,'thick',1.5)
options('du/dt','thick',1.5)
options('dup/dt','thick',1.5)
options('dn/dt','thick',1.5)
options('d/dt','thick',1.5)
options('d/dt','Color',['k','b','r'])
# options('d/dt','yrange',[-200,200])


# options(terms,'yrange',[-8,8])
# options(terms,'yrange',[-180,180])
# options('divKe', 'yrange',[-1,1])
# options('divHe', 'yrange',[-10,10])
# options('divS', 'yrange',[-10,10])
# options('du/dt', 'yrange',[-10,10])
# options('jdote', 'yrange',[-10,10])
#options(pterms,'yrange',[-1,1])
width = 10
tsmooth('divS', width=width,new_names = 'div_S', preserve_nans=0)
tsmooth('JdotE_LMN', width=width,new_names = 'j_dot_e', preserve_nans=0)
tsmooth('du/dt', width=width,new_names = 'du_dt', preserve_nans=0)
tplot_options('axis_font_size', 15)
tplot_options('vertical_spacing', 0.15)
options(['JdotE_LMN','divS','divHe','divQe','divKe'],'thick,1.5')
tplot(['JdotE_LMN','divS','divHe','divQe','divKe'])

options(['dn/dt','divE'], 'thick',1.5)
# tplot(['dn/dt','divE'])
#tplot(['JdotE_LMN','divHe','d/dt'])
#tplot('d/dt')
#tplot(['j_dot_e','div_S','du_dt','divHe'])
# options(['jdote','jdotep','jdote_mom','jdotep_mom'],'thick',1.5)
# tplot(['jdote','jdotep','jdote_mom','jdotep_mom'])
#tplot(Heterms)
# Clean up interp schemes bc j.e looks a little larger now for the 12/09 evbent
# %%

options(['JdotE_LMN','divS','divHe','divQe','divKe','d/dt'], 'thick', 1.5)
options('divHe','yrange',[-60,60])
options('divQe','yrange',[-60,60])
options('divKe','yrange',[-0.4,0.4])
tplot(['divS','divHe','divQe','divKe'])
# tplot(['JdotE_LMN','Convec_uemconts','convect'])
#tplot('Convec_uemconts')
#tplot('convect')
#tplot('Convec_dndt')
# %%


