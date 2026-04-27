#%%


import pyspedas
import numpy as np
import matplotlib.pyplot as plt
import pytplot
from pytplot import get_data, store_data,timespan,tplot_options
from pyspedas import tplot,tplot_options
from pyspedas import tinterpol
from datetime import datetime
from scipy.interpolate import interp1d

me = 9.1094e-31 #kg
mi = 1837*me
mu0 = 1.2566370e-06  #;m kg / C^2
eps0 = 8.85e-12   # C^2/Nm^2
e = 1.602e-19 #C
Z = 1 # 1 for H+, 2 for He2+
gamma = 5/3
kb = 1.380649e-23
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


# trange = ['2022-08-30/00:00', '2022-08-30/00:40'] 
trange = ['2022-08-30/00:00', '2022-09-21/00:00'] #Full interval
trange = ['2022-02-27/00:00', '2022-02-28/00:00'] #Full interval

import pyspedas
from pyspedas import tplot



mag_vars = pyspedas.projects.solo.mag(trange=trange, datatype='rtn-normal',time_clip=False)
swa_vars = pyspedas.projects.solo.swa(trange=trange,datatype = 'pas-grnd-mom',time_clip=False)

swa_vars
# %%

Bvec_name = 'B_RTN'
PiTensor_name = 'P_RTN'
ni_name = 'N'


interpvar_name = ni_name
timeax = pytplot.get_data(interpvar_name).times
tinterpol(Bvec_name,interpvar_name,newname='B')
tinterpol(PiTensor_name,interpvar_name,newname='PiTensor')
tinterpol(ni_name,interpvar_name,newname='ni')

Bvecs = 1e-9*reform(pytplot.get_data('B'))
PiTensor = 1e6*reform(pytplot.get_data('PiTensor')) #Convert to J/m^3
ni = 1e6*reform(pytplot.get_data('ni'))

Tpar = np.zeros_like(Bvecs[:,0])
Tperp = np.zeros_like(Tpar)
Pperp = np.zeros_like(Tpar)
Ppar = np.zeros_like(Tpar)

for i in range(len(Bvecs)):
	Tpar[i],Tperp[i],Ppar[i],Pperp[i] = get_parperps(ni[i],PiTensor[i]/ni[i],Bvecs[i])

store_data('Tpar', data = {'x':timeax,'y':6.242e18*Tpar}) #Converted from J to eV
store_data('Tperp', data = {'x':timeax,'y':6.242e18*Tperp})
#%%
time_perp, Tperp = get_data('Tperp')  # adjust variable name as needed
time_par, Tpar = get_data('Tpar')      # adjust variable name as needed


start_time = time_perp[0]
end_time = time_perp[-1]
time_1min = np.arange(start_time, end_time, 60)  # 60 seconds



f_perp = interp1d(time_perp, Tperp, bounds_error=False, fill_value=np.nan)
f_par = interp1d(time_par, Tpar, bounds_error=False, fill_value=np.nan)

Tperp_1min = f_perp(time_1min)
Tpar_1min = f_par(time_1min)

# Convert unix time to datetime components
from datetime import datetime, timezone

output_data = []
for i, t in enumerate(time_1min):
    dt = datetime.fromtimestamp(t, tz=timezone.utc)
    
    # Format: year, month, day, hour, min, sec, Tperp, Tpar
    line = [
        dt.year,
        dt.month,
        dt.day,
        dt.hour,
        dt.minute,
        dt.second + dt.microsecond/1e6,  # seconds with fractional part
        Tperp_1min[i],
        Tpar_1min[i]
    ]
    output_data.append(line)

output_data = np.array(output_data)

# Save to ascii
np.savetxt('Solo_temperatures.txt', output_data, 
           fmt='%4d %02d %02d %02d %02d %09.4f   %.4e   %.4e',
           header='Y M D H Min Sec Tperp(eV) Tpar(eV)',
           comments='# ')

print("Data saved to solo_temperatures.txt")
# %%
