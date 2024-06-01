# Synthetic Photometry
Compute synthetic photometry (magnitudes and fluxes) from spectra for different filters

## Input parameters
* wl: wavelength in um
* flux and eflux (optional) in units specified by flux_unit
* flux_unit: flux and error units
	* 'erg/s/cm2/A'
	* 'Jy'
* filters: filters to derive synthetic photometry following [SVO filter names](http://svo2.cab.inta-csic.es/theory/fps/)

## Returns
------
* out: python dictionary with the following parameters for each filter
	* out['lambda_eff(um)']: effective wavelength in micron
	* out['width_eff(um)']: effective width in micron
	* out['syn_flux(erg/s/cm2/A)']: synthetic flux in erg/s/cm2/A
	* out['esyn_flux(erg/s/cm2/A)']: synthetic flux errors in erg/s/cm2/A (if input flux errors are provided)
	* out['syn_flux(Jy)']: synthetic flux in Jy
	* out['esyn_flux(Jy)']: synthetic flux errors in Jy (if input flux errors are provided)
	* out['syn_mag']: synthetic magnitude
	* out['esyn_mag']: synthetic magnitude error (if input flux errors are provided)

## Example
```
from synthetic_photometry import synthetic_photometry
```
assume we have a spectrum wavelength (wl in um), flux (in erg/s/cm2/A), and flux error (eflux) and we want synthetic photometry for several filters
```
filters = (['Spitzer/IRAC.I1', 'WISE/WISE.W1']) # filters of interest
```
run the code
```
out = synthetic_photometry(wl=wl, flux=flux, eflux=eflux, flux_unit='erg/s/cm2/A', filters=filters)
```
output
```
eff_wl = out['lambda_eff(um)'] # effective wavelength (um) for all filters
eff_width = out['width_eff(um)'] # effective width (um) for all filters
flux_syn = out['syn_flux(erg/s/cm2/A)'] # synthetic flux (erg/s/cm2/A) for all filters
eflux_syn = out['esyn_flux(erg/s/cm2/A)'] # synthetic flux errors (erg/s/cm2/A) for all filters
flux_Jy_syn = out['syn_flux(Jy)'] # synthetic flux (Jy) for all filters
eflux_Jy_syn = out['esyn_flux(Jy)'] # synthetic flux errors (Jy) for all filters
mag_syn = out['syn_mag'] # synthetic magnitude for all filters
emag_syn = out['esyn_mag'] # synthetic magnitude error for all filters
```
