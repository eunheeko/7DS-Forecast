import numpy as np
from scipy import interpolate 
from astropy.table import Table
from astropy import units as u

# import os
# from svnds import PATH_DATA

###### Constants
H = 6.626e-27 # erg/Hz
C = 3e10      # cm/s 

###### Telescope
D = 50.8             # Aperture size [cm]
D_OBSCURATION = 29.8   # Central Obscuration (diameter) [cm]
EFL = 1537.3           # Effective FL [mm]

D_EFF = np.sqrt(D**2 - D_OBSCURATION**2)
D_EFF

###### Detector
ARRAY = 'CMOS'       # detector array type - verified, websete
DQ_RN = 3.           # [e], readout noise - estimated, websete
I_DARK = 0.01        # [e/s], dark current - seems to be very low ???
PIXEL_SIZE = 3.76    # [um], "pitch" - verified, websete
THETA_PIXEL = PIXEL_SIZE / EFL * 206.265  # [arcsec], pixel scale- estimated, websete 
NXPIX, NYPIX = 9576, 6388  # [pixels], detector format, approx. 9k x 6k - verified, websete

FWHM_SEEING = 1.5     # [arcsec]
FWHM0 = FWHM_SEEING   # analysis aperture size

# How many pixels does a point source occupy?
# Effective number of pixels for a Gaussian PSF with FWHM0
NPIX_PTSRC = np.pi*(FWHM0/THETA_PIXEL)**2

###### Lens
eff_lambdas = np.arange(400, 1100, 100) / 1000 #um
eff_lens1 = np.array([0.997, 0.997, 0.998, 0.998, 0.998, 0.998, 0.998])
eff_lens2 = np.array([0.974, 0.995, 0.997, 0.997, 0.998, 0.998, 0.998])
eff_lens3 = np.array([0.996, 0.998, 0.998, 0.998, 0.998, 0.998, 0.998])
eff_vis = np.array([0.988, 0.998, 0.999, 0.996, 0.974, 0.936, 0.900])
eff_al = np.array([0.86, 0.95, 0.95, 0.92, 0.85, 0.83, 0.85])

eff_optics = eff_lens1 * eff_lens2 * eff_lens3 * eff_vis**5 * eff_al**2
ref_lambdas = np.linspace(0.366, 1.0, 1001) #um
smooth_opt = interpolate.interp1d(eff_lambdas, eff_optics, kind = 'quadratic')
eff_optics_smooth = np.interp(ref_lambdas, np.linspace(0.4, 1, 1001), smooth_opt(np.linspace(0.4, 1, 1001)))

qe_lam = [366, 405, 450, 492, 550, 589, 694, 768, 853, 905, 1060]
qe_lam = np.array(qe_lam) / 1000
qe_tot = [28.5, 68.4, 82.5, 84.3, 75.6, 67.8, 42.9, 24.7, 13.9, 10.1, 0.6]
qe_tot = np.array(qe_tot) / 100

T_qe = Table([qe_lam, qe_tot], names = ['wavelength', 'QE'])
T_qe['wavelength'].unit = u.um
T_qe['wavelength'].format = '8.4f'
smooth_qe = interpolate.interp1d(T_qe['wavelength'], T_qe['QE'], kind = 'quadratic')

eff_fpa = np.interp(eff_lambdas, T_qe['wavelength'], T_qe['QE'])
eff_fpa_smooth = smooth_qe(ref_lambdas)

#  Total Efficeincy
eff_total = eff_optics * eff_fpa
# Smooth total efficiency
eff_total_smooth = eff_optics_smooth * eff_fpa_smooth
