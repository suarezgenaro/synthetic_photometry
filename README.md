# Synthetic Photometry
Compute synthetic photometry (magnitudes and fluxes) from spectra for different filters

## Input parameters
* wl: wavelength in um
* flux and eflux (optional) in units specified by flux_unit
* flux_unit: flux and error units
	* 'erg/s/cm2/A'
	* 'Jy'
* filters: filters to derive synthetic photometry following SVO filter names 

## Returns
------
* out_synthetic_photometry: python dictionary with the following parameters for each filter
	* out_synthetic_photometry['lambda_eff(um)']: effective wavelength in micron
	* out_synthetic_photometry['width_eff(um)']: effective width in micron
	* out_synthetic_photometry['syn_flux(erg/s/cm2/A)']: synthetic flux in erg/s/cm2/A
	* out_synthetic_photometry['esyn_flux(erg/s/cm2/A)']: synthetic flux errors in erg/s/cm2/A (if input flux errors are provided)
	* out_synthetic_photometry['syn_flux(Jy)']: synthetic flux in Jy
	* out_synthetic_photometry['esyn_flux(Jy)']: synthetic flux errors in Jy (if input flux errors are provided)
	* out_synthetic_photometry['syn_mag']: synthetic magnitude
	* out_synthetic_photometry['esyn_mag']: synthetic magnitude error (if input flux errors are provided)

## Example
```
import synthetic_photometry as synthetic_photometry
```
assume we have a spectrum wavelength (wl in um), flux (in erg/s/cm2/A), and flux error (eflux) and we want synthetic photometry for several filters
```
filters = (['Spitzer/IRAC.I1', 'WISE/WISE.W1']) # filters of interest
```
run the code
```
out_synthetic_photometry = synthetic_photometry.synthetic_photometry(wl=wl, flux=flux, eflux=eflux, flux_unit='erg/s/cm2/A', filters=filters)
```
# output
* eff_wl = out_synthetic_photometry['lambda_eff(um)'] # effective wavelength (um) for each filter
* eff_width = out_synthetic_photometry['width_eff(um)'] # effective width (um) for each filter
* flux_syn = out_synthetic_photometry['syn_flux(erg/s/cm2/A)'] # synthetic flux (erg/s/cm2/A) for each filter
* eflux_syn = out_synthetic_photometry['esyn_flux(erg/s/cm2/A)'] # synthetic flux errors (erg/s/cm2/A) for each filter
* flux_Jy_syn = out_synthetic_photometry['syn_flux(Jy)'] # synthetic flux (Jy) for each filter
* eflux_Jy_syn = out_synthetic_photometry['esyn_flux(Jy)'] # synthetic flux errors (Jy) for each filter
* mag_syn = out_synthetic_photometry['syn_mag'] # synthetic magnitude for each filter
* emag_syn = out_synthetic_photometry['esyn_mag'] # synthetic magnitude error for each filter
