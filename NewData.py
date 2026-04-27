#%%


import pyspedas
import numpy as np
import matplotlib.pyplot as plt
import pytplot
from pytplot import get_data, store_data,timespan
from pyspedas import tplot
from pyspedas import tinterpol
from spacepy import pycdf
import datetime
import cdflib


me = 9.1094e-31 #kg
mi = 1837*me
mu0 = 1.2566370e-06  #;m kg / C^2
eps0 = 8.85e-12   # C^2/Nm^2
e = 1.602e-19 #C
Z = 1 # 1 for H+, 2 for He2+
gamma = 5/3
kb = 1.380649e-23



# %%
# import spacepy.pycdf as pycdf
# cdf_file1 = pycdf.CDF('~/psp_data/sweap/spi/l3/spi_sf00_l3_mom/2024/psp_swp_spi_sf00_l3_mom_20241224_v00.cdf')
# cdf_file2 = pycdf.CDF('~/psp_data/sweap/spi/l3/spi_sf00_l3_mom/2024/psp_swp_spi_sf00_l3_mom_20241030_v00.cdf')



#%%


#trange=['2023-12-28/00:00','2023-12-29/00:00'] # E18
# trange = ['2024-12-24/00:00', '2024-12-24/12:00'] # E22 
trange = ['2025-06-15/00:00', '2025-06-16/00:00'] # E24 


spi_vars = pyspedas.projects.psp.spi(trange=trange,level='l3',time_clip=True,no_update=True)


tplot('psp_spi_TEMP')
# %%
import cdflib
from pytplot import store_data,options,tplot


cdf_file_path = os.path.expanduser('~/psp_data/sweap/spi/l3/spi_sf00_l3_mom/2024/psp_swp_spi_sf00_l3_mom_20241224_v00.cdf')
cdf_file = cdflib.CDF(cdf_file_path)

print(cdf_file.cdf_info())


# %%
epoch_data = cdf_file.varget('Epoch')

temp_data = cdf_file.varget('TEMP')



store_data('TempVar', data = {'x':epoch_data, 'y':temp_data})


# Set Y-axis label

# Set plot title (optional)
options('TempVar', 'ytitle', 'Parker Solar Probe TEMP Data')



tplot('TempVar')
#%%

from datetime import datetime
epoch_dat = epoch_data[0]

real_date = datetime.utcfromtimestamp(epoch_dat)




# info = cdf_file.cdf_info()
# print("CDF Info:")
# print(info)

#%%

import pandas as pd

# Your example nanosecond timestamp
epoch_data_nanoseconds = 788356867282672128

# Convert to a datetime object
# The 'unit' parameter tells pandas what unit the integer represents
real_date = pd.to_datetime(epoch_data_nanoseconds, unit='ns')

# Print the resulting datetime object
print(real_date)



