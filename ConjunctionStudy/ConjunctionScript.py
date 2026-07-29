#%%

#%%

PSPinterval = ['2025-03-22/12:00','2025-03-25/12:00']
SOLOinterval = ['2025-03-22/12:00','2025-03-25/12:00']
WINDinterval = ['2025-03-22/12:00','2025-03-25/12:00']

def get_PSP(trange):
	Bfld_vars = pyspedas.projects.psp.fields(trange=trange, level='l2', time_clip=True)
	swp_vars = pyspedas.projects.psp.spi(trange=trange,level='l3',get_support_data=True,time_clip=True)
	qtn_vars = pyspedas.projects.psp.fields(trange=trange,level='l3',datatype='sqtn_rfs_V1V2',time_clip=True)

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

    return PSP_data 

def get_SOLO(trange):
    return SOLO_data

def get_WIND(trange):
    return WIND_data