#%%
from solarmach import SolarMACH, print_body_list, get_sw_speed, get_gong_map, calculate_pfss_solution
import datetime

# Set the date of interest
date = datetime.datetime(2025, 9, 21) # Close PSP SO approach
# date = datetime.datetime(2021, 8, 11,17,0,0) # Soni Conjunction
# date = datetime.datetime(2022, 2, 27,0,0,0) 
# date = datetime.datetime(2025, 3, 24,12,0,0) 

date_psp = datetime.datetime(2025, 3, 23, 00, 10) 
date_solo = datetime.datetime(2025, 3, 31, 0, 0)
date_wind = datetime.datetime(2025, 3, 28, 0, 0)
date_maven = datetime.datetime(2025, 3, 29, 12, 0) #Ballpark estimate of when MAVEN will be in stream conjunction

# List of bodies/spacecraft to plot
# bodies = ["PSP", "SOLO","WIND", "BepiColombo"]
psp_body = ["PSP"]
solo_body = ["SOLO"]
wind_body = ["WIND"]
maven_body = ["MAVEN"]
allbods = ["PSP", "SOLO","WIND", "MAVEN"]
# bodies = ["PSP"]

# Create the SolarMACH object with the required arguments
sm_psp = SolarMACH(date_psp, body_list=psp_body)
sm_solo = SolarMACH(date_solo, body_list=solo_body)
sm_wind = SolarMACH(date_wind, body_list=wind_body)
sm_maven = SolarMACH(date_maven, body_list=maven_body)

# Plot the configuration
# sm.plot(plot_sun_body_line=True,markers='numbers',long_offset=0,long_sector=[140,200],long_sector_vsw=[200,600],long_sector_color='grey')
sm_psp.plot(plot_sun_body_line=False,markers=False,long_offset=0,background_spirals=[6,318])
sm_solo.plot(plot_sun_body_line=False,markers=False,long_offset=-105,background_spirals=[6,318])
sm_wind.plot(plot_sun_body_line=False,markers=False,long_offset=-65,background_spirals=[6,318])
sm_maven.plot(plot_sun_body_line=False,markers=False,long_offset=-85,background_spirals=[6,318])

sm_maven.coord_table
# %%

date_psp1 = datetime.datetime(2025, 3, 22, 20, 10) 
date_psp2 = datetime.datetime(2025, 3, 23, 9, 20) 

date_solo1 = datetime.datetime(2025, 3, 28, 19, 30)
date_solo2 = datetime.datetime(2025, 4, 3, 14, 50)


date_wind1 = datetime.datetime(2025, 3, 26, 0, 30)
date_wind2 = datetime.datetime(2025, 3, 30, 13, 40)

psp1 = SolarMACH(date=date_psp1, body_list=psp_body, coord_sys='Stonyhurst')
psp2 = SolarMACH(date=date_psp2, body_list=psp_body, coord_sys='Stonyhurst')
solo1 = SolarMACH(date=date_solo1, body_list=solo_body, coord_sys='Stonyhurst')
solo2 = SolarMACH(date=date_solo2, body_list=solo_body, coord_sys='Stonyhurst')
wind1 = SolarMACH(date=date_wind1, body_list=wind_body, coord_sys='Stonyhurst')
wind2 = SolarMACH(date=date_wind2, body_list=wind_body, coord_sys='Stonyhurst')
df_psp1 = psp1.coord_table
df_psp2 = psp2.coord_table
df_solo1 = solo1.coord_table
df_solo2 = solo2.coord_table
df_wind1 = wind1.coord_table
df_wind2 = wind2.coord_table
#
# %%

display(df_psp1.coord_table)

# %%
