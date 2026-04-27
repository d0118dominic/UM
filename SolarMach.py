#%%
from solarmach import SolarMACH, print_body_list, get_sw_speed, get_gong_map, calculate_pfss_solution
import datetime

# Set the date of interest
date = datetime.datetime(2025, 9, 21) # Close PSP SO approach
# date = datetime.datetime(2021, 8, 11,17,0,0) # Soni Conjunction
date = datetime.datetime(2022, 2, 27,0,0,0) 

# List of bodies/spacecraft to plot
bodies = ["PSP", "SOLO","WIND","STEREO-A"]
# bodies = ["PSP"]

# Create the SolarMACH object with the required arguments
sm = SolarMACH(date, body_list=bodies)


# Plot the configuration
sm.plot(plot_sun_body_line=True,markers='numbers',long_offset=0,long_sector=[140,200],long_sector_vsw=[200,400],long_sector_color='grey')
# sm.plot(plot_sun_body_line=False,markers='numbers',long_offset=0,background_spirals=[6,318])

# %%

sm4a = SolarMACH(date=date, body_list=bodies, coord_sys='Stonyhurst')
df = sm4a.coord_table
display(df)

# %%
display(gong_map)
