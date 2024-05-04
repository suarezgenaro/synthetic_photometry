# compute synthetic magnitude and fluxes from spectra using several filters
from sys import exit
from astropy.io import ascii
import numpy as np

def synthetic_phot(wl, flux, filters, flux_unit, eflux=None): 
	'''
	wl: wavelength in um
	flux and eflux (optional) in units specified by flux_unit
	flux_unit: flux and error units
				'erg/s/cm2/A'
				'Jy'
	filters: filters to derive synthetic photometry following SVO filter names 

	Example
	>>> import synthetic_phot as synthetic_phot
	>>> # assume we have a spectrum wavelength (wl in um), flux (in erg/s/cm2/A), and flux error (eflux) and we want synthetic photometry for several filters
	>>> filters = (['Spitzer/IRAC.I1', 'WISE/WISE.W1']) # filters of interest
	>>>	# run the code
	>>> out_synthetic_phot = synthetic_phot.synthetic_phot(wl=wl, flux=flux, eflux=eflux, flux_unit='erg/s/cm2/A', filters=filters)
	>>> # output
	>>> eff_wl = out_synthetic_phot['lambda_eff(um)'] # effective wavelength (um) for each filter
	>>> eff_width = out_synthetic_phot['width_eff(um)'] # effective width (um) for each filter
	>>> flux_syn = out_synthetic_phot['syn_flux(erg/s/cm2/A)'] # synthetic flux (erg/s/cm2/A) for each filter
	>>> eflux_syn = out_synthetic_phot['esyn_flux(erg/s/cm2/A)'] # synthetic flux errors (erg/s/cm2/A) for each filter
	>>> flux_Jy_syn = out_synthetic_phot['syn_flux(Jy)'] # synthetic flux (Jy) for each filter
	>>> eflux_Jy_syn = out_synthetic_phot['esyn_flux(Jy)'] # synthetic flux errors (Jy) for each filter
	>>> mag_syn = out_synthetic_phot['syn_mag'] # synthetic magnitude for each filter
	>>> emag_syn = out_synthetic_phot['esyn_mag'] # synthetic magnitude error for each filter

	'''

	mask_nonan = ~np.isnan(flux)
	wl = wl[mask_nonan]
	flux = flux[mask_nonan]
	if (eflux is not None): eflux = eflux[mask_nonan]

	# when filters parameter is given as a string, convert the variable to a list so len(filters) returns 1 rather than the string length
	if (type(filters) is str): filters = ([filters])

	# read filters response curves and zero points
	filter_parms_file = '/home/gsuarez/TRABAJO/PDFs/Papers/data/Filter_Transmissions/filters_properties_SVO'
	filter_parms = ascii.read(filter_parms_file)
	name_filt = filter_parms['filter']
	f0_filt = filter_parms['zero(Jy)'] # in Jy

	# arrays to store relevant information
	syn_flux_Jy = np.zeros(len(filters))
	esyn_flux_Jy = np.zeros(len(filters))
	syn_flux_erg = np.zeros(len(filters))
	esyn_flux_erg = np.zeros(len(filters))
	syn_mag = np.zeros(len(filters))
	esyn_mag = np.zeros(len(filters))
	lambda_eff = np.zeros(len(filters))
	width_eff = np.zeros(len(filters))
	for k in range(len(filters)): # iterate for each filter
		# read filter response
		path_filters = '/home/gsuarez/TRABAJO/PDFs/Papers/data/Filter_Transmissions/'
		filter_response = ascii.read(path_filters+filters[k].replace('/', '_')+'.dat')
	
		try: filter_response
		except NameError: print(f'ERROR: NO FILTER RESPONSE FILE FOR FILTER {filters[k]}'), exit()

		filter_wl = filter_response['col1'] / 1e4 # in um
		filter_flux = filter_response['col2'] # filter response (named filter_flux just for ease)
	
		# verify the spectrum fully covers the filter response
		if ((filter_wl.min()<wl.min()) | (filter_wl.max()>wl.max())): print(f'CAVEAT: NO FULL SPECTRAL COVERAGE FOR FILTER {filters[k]}')
		if ((filter_wl.min()<wl.max()) | (filter_wl.max()>wl.min())): # filter covered partially or fully
			# wavelength dispersion of the spectrum in the filter wavelength range
			mask_wl = (wl>=filter_wl.min()) & (wl<=filter_wl.max())
	
			wl_disp = wl[mask_wl][1:] - wl[mask_wl][:-1] # dispersion of spectrum spectra (um)
			wl_disp = np.append(wl_disp, wl_disp[-1]) # add an element equal to the last row to keep the same shape as the wl array
	
			# synthetic photometry
			# resample filter response to the spectrum wavelength
			filter_flux_resam = np.interp(wl[mask_wl], filter_wl, filter_flux) # dimensionless
			
			# normalize the response curve (it was dimensionless but now it has 1/um units)
			filter_flux_resam_norm = filter_flux_resam / sum(filter_flux_resam*wl_disp) # 1/um
			
			# synthetic flux density
			syn_flux = sum(flux[mask_wl]*filter_flux_resam_norm*wl_disp) # in input flux units (erg/s/cm2/A or Jy)
			if (eflux is not None): esyn_flux = np.median(eflux[mask_wl]/flux[mask_wl]) * syn_flux # synthetic flux error as the median fractional flux uncertainties in the filter passband
			
			# compute the effective wavelength and effective width for each filter
			lambda_eff[k] = sum(wl[mask_wl]*filter_flux_resam*flux[mask_wl]*wl_disp) / sum(filter_flux_resam*flux[mask_wl]*wl_disp) # um
			width_eff[k] = sum(filter_flux_resam*wl_disp) / filter_flux_resam.max() # um
			
			# convert flux to magnitudes
			# first from erg/s/cm2/A to Jy (if needed) and then from Jy to mag
			if (flux_unit=='erg/s/cm2/A'):
				#syn_flux_Jy = syn_flux*(lambda_eff*1e-4)**2 / 3e-21 # Jy
				syn_flux_Jy[k] = 3.33564095e+04*syn_flux*(lambda_eff[k]*1e4)**2 # Jy (the above option gives the same result)
				if (eflux is not None): esyn_flux_Jy[k] = esyn_flux/syn_flux * syn_flux_Jy[k] # Jy
	
				syn_flux_erg[k] = syn_flux # erg/s/cm2/A
				if (eflux is not None): esyn_flux_erg[k] = esyn_flux # erg/s/cm2/A
	
			if (flux_unit=='Jy'): # convert Jy to erg/s/cm2/A to be an output
				syn_flux_erg[k] = syn_flux/(3.33564095e+04*(lambda_eff[k]*1e4)**2) # erg/s/cm2/A
				if (eflux is not None): esyn_flux_erg[k] = esyn_flux/syn_flux * syn_flux_erg[k] # erg/s/cm2/A
	
				syn_flux_Jy[k] = syn_flux # Jy
				if (eflux is not None): esyn_flux_Jy[k] = esyn_flux # Jy
	
			# from Jy to mag
			mask = name_filt==filters[k]
			if any(mask) is False: print(f'\nERROR: NO ZERO POINT FOR FILTER {filters[k]}'), exit()
			syn_mag[k] = -2.5*np.log10(syn_flux_Jy[k]/f0_filt[mask]) # in mag
			if (eflux is not None): esyn_mag[k] = (2.5/np.log(10))*np.sqrt((esyn_flux_Jy[k]/syn_flux_Jy[k])**2)#+(ephot_F0/phot_F0)**2) # in mag
		
			del filter_response # remove variable with filter response so it won't exit if an input filter name doesn't match an existing one

			## 2MASS photometry in erg/s/cm2/A
			## first from mag to Jy
			#if (fil_syn_flux=='J'):
			#	twoMASS_Jy = phot_F0 * 10**(J_2MASS/(-2.5))
			#	etwoMASS_Jy = 10**(-J_2MASS//2.5)*np.sqrt(ephot_F0**2+((-1*np.log(10)/2.5)*phot_F0*eJ_2MASS)**2)
			#if (fil_syn_flux=='H'):
			#	twoMASS_Jy = phot_F0 * 10**(H_2MASS/(-2.5))
			#	etwoMASS_Jy = 10**(-H_2MASS//2.5)*np.sqrt(ephot_F0**2+((-1*np.log(10)/2.5)*phot_F0*eH_2MASS)**2)
			#if (fil_syn_flux=='K'):
			#	twoMASS_Jy = phot_F0 * 10**(K_2MASS/(-2.5))
			#	etwoMASS_Jy = 10**(-K_2MASS//2.5)*np.sqrt(ephot_F0**2+((-1*np.log(10)/2.5)*phot_F0*eK_2MASS)**2)
		
			## from Jy to erg/s/cm2/A
			#twoMASS_erg[i] = twoMASS_Jy / (3.33564095e+04*(lambda_eff*1e4)**2)
			#etwoMASS_erg[i] = twoMASS_erg[i] * etwoMASS_Jy/twoMASS_Jy # in erg/s/cm2/A

		out = {'syn_flux(Jy)': syn_flux_Jy, 'syn_flux(erg/s/cm2/A)': syn_flux_erg, 'syn_mag': syn_mag, 'lambda_eff(um)': lambda_eff, 'width_eff(um)': width_eff}
		if (eflux is not None): out['esyn_flux(Jy)'] = esyn_flux_Jy
		if (eflux is not None): out['esyn_flux(erg/s/cm2/A)'] = esyn_flux_erg
		if (eflux is not None): out['esyn_mag'] = esyn_mag

	return out
