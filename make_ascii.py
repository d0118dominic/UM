#%%
import pyspedas
from pytplot import get_data
import numpy as np
from datetime import datetime

# Load PSP SWEAP data (proton temperatures)
# Adjust trange to your needs
# trange = ['2022-09-06/00:00:00', '2022-09-06/23:59:59']

# Load the data
# psp_data = pyspedas.psp.spc(trange=trange, datatype='l3i', level='l3', time_clip=True)

# Get the temperature data
# Variable names might be 'psp_spc_Tperp' and 'psp_spc_Tpar' or similar
# Check what's loaded with: from pytplot import tplot_names; print
# (tplot_names())

time_perp, Tperp = get_data('Tperp')  # adjust variable name as needed
time_par, Tpar = get_data('Tpar')      # adjust variable name as needed

# Resample to 1-minute cadence if needed
# (if data is already 1-min cadence, skip this)
from scipy.interpolate import interp1d

# Create 1-minute time grid
start_time = time_perp[0]
end_time = time_perp[-1]
time_1min = np.arange(start_time, end_time, 60)  # 60 seconds

# Interpolate to 1-min cadence
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
np.savetxt('psp_temperatures.txt', output_data, 
           fmt='%4d %02d %02d %02d %02d %09.4f   %.4e   %.4e',
           header='Y M D H Min Sec Tperp(eV) Tpar(eV)',
           comments='# ')

print("Data saved to psp_temps.txt")



# %%
