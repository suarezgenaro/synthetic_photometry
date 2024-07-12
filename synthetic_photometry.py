from astropy.io import ascii
from astropy.table import Table
import numpy as np
import os
from sys import exit

def synthetic_photometry(wl, flux, filters, flux_unit, eflux=None): 
	'''
	Compute synthetic magnitudes and fluxes from spectra for different filters

	Input parameters:
	wl: float array
		wavelength in um
	flux: float array
		fluxes in units specified by flux_unit
	eflux: flux array (optional) 
		flux errors in units specified by flux_unit
	flux_unit: str
		flux and flux error units
			'erg/s/cm2/A'
			'Jy'
	filters: str
		filters (following SVO filter IDs) to derive synthetic photometry

	Returns
	------
	out: python dictionary with the following parameters for each filter
		out['lambda_eff(um)']: effective wavelength in micron
		out['width_eff(um)']: effective width in micron
		out['syn_flux(erg/s/cm2/A)']: synthetic flux in erg/s/cm2/A
		out['esyn_flux(erg/s/cm2/A)']: synthetic flux errors in erg/s/cm2/A (if input flux errors are provided)
		out['syn_flux(Jy)']: synthetic flux in Jy
		out['esyn_flux(Jy)']: synthetic flux errors in Jy (if input flux errors are provided)
		out['syn_mag']: synthetic magnitude
		out['esyn_mag']: synthetic magnitude error (if input flux errors are provided)

	Example
	>>> import synthetic_photometry
	>>> # assume we have a spectrum wavelength (wl in um), flux (in erg/s/cm2/A), and flux error (eflux) 
	>>>	# and we want synthetic photometry for several filters
	>>> filters = (['Spitzer/IRAC.I1', 'WISE/WISE.W1']) # filters of interest
	>>>	# run the code
	>>> out = synthetic_photometry.synthetic_photometry(wl=wl, flux=flux, eflux=eflux, flux_unit='erg/s/cm2/A', filters=filters)
	>>> # output
	>>> eff_wl = out['lambda_eff(um)'] # effective wavelength (um) for all filters
	>>> eff_width = out['width_eff(um)'] # effective width (um) for all filters
	>>> flux_syn = out['syn_flux(erg/s/cm2/A)'] # synthetic flux (erg/s/cm2/A) for all filters
	>>> eflux_syn = out['esyn_flux(erg/s/cm2/A)'] # synthetic flux errors (erg/s/cm2/A) for all filters
	>>> flux_Jy_syn = out['syn_flux(Jy)'] # synthetic flux (Jy) for all filters
	>>> eflux_Jy_syn = out['esyn_flux(Jy)'] # synthetic flux errors (Jy) for all filters
	>>> mag_syn = out['syn_mag'] # synthetic magnitude for all filters
	>>> emag_syn = out['esyn_mag'] # synthetic magnitude error for all filters

	Author: Genaro Suárez

	Modification history
		2024/05/24	SVO VOTable table is read from a link with the most updated filter parameters
		2024/05/09	included VOTable table with all SVO filter zero points
		2024/05/07	filter transmissions are read and downloaded from SVO, if not already stored locally
		2024/04/20	function created
		2021      	functional code not as a function

	'''

	path_synthetic_photometry = os.path.dirname(__file__)+'/'

	mask_nonan = ~np.isnan(flux)
	wl = wl[mask_nonan]
	flux = flux[mask_nonan]
	if (eflux is not None): eflux = eflux[mask_nonan]

	# when filters parameter is given as a string, convert the variable to a list so len(filters) returns 1 rather than the string length
	if (type(filters) is str): filters = ([filters])

	# read filters transmission curves and zero points
	svo_table = f'{path_synthetic_photometry}/FPS_info.xml'
	if os.path.exists(svo_table): 
		svo_data = Table.read(svo_table, format='votable') # open downloaded table with filters' info
	else:
		svo_data = Table.read('https://svo.cab.inta-csic.es/files/svo/Public/HowTo/FPS/FPS_info.xml', format='votable') # this link will be updated as soon as new filters are added to FPS. 
		svo_data.write(svo_table, format='votable') # save table to avoid reading from the web the table each time the code is run, which can take a few seconds
	filterID = svo_data['filterID'] # VSO ID
	ZeroPoint = svo_data['ZeroPoint'] # Jy

	# arrays to store relevant information
	syn_flux_Jy = np.zeros(len(filters))
	esyn_flux_Jy = np.zeros(len(filters))
	syn_flux_erg = np.zeros(len(filters))
	esyn_flux_erg = np.zeros(len(filters))
	syn_mag = np.zeros(len(filters))
	esyn_mag = np.zeros(len(filters))
	lambda_eff = np.zeros(len(filters))
	width_eff = np.zeros(len(filters))
	zero_point = np.zeros(len(filters))
	# assign NaN values (this will be the output for an input filter name not recognized by the SVO)
	syn_flux_Jy[:] = np.nan
	esyn_flux_Jy[:] = np.nan
	syn_flux_erg[:] = np.nan
	esyn_flux_erg[:] = np.nan
	syn_mag[:] = np.nan
	esyn_mag[:] = np.nan
	lambda_eff[:] = np.nan
	width_eff[:] = np.nan
	zero_point[:] = np.nan
	for k in range(len(filters)): # iterate for each filter
		# check first if the filter name is on the VSO
		if not filters[k] in filterID: print(f'Caveat for {filters[k]}: FILTER NOT RECOGNIZED BY THE SVO, so will be ignored')
		else:
			# read filter transmission
			# check if the filter transmission exits locally already
			path_filter_transmissions = f'{path_synthetic_photometry}/filter_transmissions/'
			if not os.path.exists(path_filter_transmissions): os.makedirs(path_filter_transmissions) # make directory (if not existing) to store filter transmissions
			filter_transmission_name = filters[k].replace('/', '_')+'.dat' # when filter name includes '/' replace it by '_'
			if not os.path.exists(path_filter_transmissions+filter_transmission_name): # filter transmission does not exits yet
				print(f'\nreading and storing filter {filters[k]} directly from the SVO')
				# read filter transmission directly from SVO
				page = f'http://svo2.cab.inta-csic.es/theory/fps/fps.php?ID={filters[k]}'
				filter_transmission = Table.read(page, format='votable')
				# save filter transmission if it doesn't exist already
				ascii.write(filter_transmission, path_filter_transmissions+filter_transmission_name, format='no_header', formats={'Wavelength': '%.1f', 'Transmission': '%.10f'})
	
			filter_transmission = ascii.read(path_filter_transmissions+filter_transmission_name) # read locally stored filter transmission
	
			filter_wl = filter_transmission['col1'] / 1e4 # in um
			filter_flux = filter_transmission['col2'] # filter transmission (named filter_flux just for ease)
		
			# verify the spectrum fully covers the filter transmission
			if ((filter_wl.max()<wl.min()) | (filter_wl.min()>wl.max())): # filter out of the spectrum coverage
				print(f'Caveat for {filters[k]}: NO WAVELENGTH COVERAGE, so synthetic photometry won\'t be obtained')
			else: # filter fully or partially covered
				if ((filter_wl.min()<wl.min()) & (filter_wl.max()>wl.min())): # blue-end of the filter partially covered
					print(f'Caveat for {filters[k]}: NO FULL SPECTRAL COVERAGE, so the synthetic value is a lower limit')
				if ((filter_wl.max()>wl.max()) & (filter_wl.min()<wl.max())): # red-end of the filter partially covered
					print(f'Caveat for {filters[k]}: NO FULL SPECTRAL COVERAGE, so the synthetic value is a lower limit')

				# wavelength dispersion of the spectrum in the filter wavelength range
				mask_wl = (wl>=filter_wl.min()) & (wl<=filter_wl.max())
		
				wl_disp = wl[mask_wl][1:] - wl[mask_wl][:-1] # dispersion of spectrum spectra (um)
				wl_disp = np.append(wl_disp, wl_disp[-1]) # add an element equal to the last row to keep the same shape as the wl array

				# synthetic photometry
				# resample filter transmission to the spectrum wavelength
				filter_flux_resam = np.interp(wl[mask_wl], filter_wl, filter_flux) # dimensionless

				# normalize the transmission curve (it was dimensionless but now it has 1/um units)
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
				mask = filterID==filters[k]
				if any(mask) is False: print(f'\nERROR: NO ZERO POINT FOR FILTER {filters[k]}'), exit()
				syn_mag[k] = -2.5*np.log10(syn_flux_Jy[k]/ZeroPoint[mask]) # in mag
				if (eflux is not None): esyn_mag[k] = (2.5/np.log(10))*np.sqrt((esyn_flux_Jy[k]/syn_flux_Jy[k])**2)#+(ephot_F0/phot_F0)**2) # in mag

				zero_point[k] = ZeroPoint[mask]

				del filter_transmission # remove variable with filter transmission so it won't exit if an input filter name doesn't match an existing one

	out_synthetic_photometry = {'syn_flux(Jy)': syn_flux_Jy, 'syn_flux(erg/s/cm2/A)': syn_flux_erg, 'syn_mag': syn_mag, 'lambda_eff(um)': lambda_eff, 'width_eff(um)': width_eff, 'zero_point(Jy)': zero_point}
	if (eflux is not None): 
		out_synthetic_photometry['esyn_flux(Jy)'] = esyn_flux_Jy
		out_synthetic_photometry['esyn_flux(erg/s/cm2/A)'] = esyn_flux_erg
		out_synthetic_photometry['esyn_mag'] = esyn_mag

	return out_synthetic_photometry 
