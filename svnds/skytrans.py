import os

from astropy.table import Table
from scipy.ndimage import gaussian_filter
from svnds import PATH_DATA

##### Sky
path_skytable = os.path.join(PATH_DATA, "systematics/", "skytable.fits")
# path_skytabel = ("./data/sys")
sky_tbl = Table.read(path_skytable)
sky_lambdas = sky_tbl['lam']/1e3 #um
trans_smooth = gaussian_filter(sky_tbl['trans'], 10)

wl_nm_sky = sky_tbl['lam']          # nm
wl_um_sky = wl_nm_sky / 1e3       # micron
wl_cm_sky = wl_um_sky / 1e4       # cm
wl_am_sky = wl_angstrom_sky = wl_nm_sky * 10  # angstrom
nu_sky = 3e18 / wl_angstrom_sky   # Hz

I_lambda_sky = sky_tbl['flux']      # [ph/s/m2/micron/arcsec2] photon reate
f_lambda_sky = I_lambda_sky * (H*C/wl_cm_sky) / (1e2**2) / (1e4)  # erg/s/cm2/A
f_nu_sky = f_lambda_sky * wl_angstrom_sky * (wl_cm_sky/C) / (1e-23 * 1e-6)  # micro Jansky