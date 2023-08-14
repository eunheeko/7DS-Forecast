import os
import numpy as np
from scipy.integrate import trapezoid


#Load EL-COSMOS catalog
def load_catalog(path_elcosmos, zone):
    """load_catalog returns the spec (EL-COSMOS identifier) for a given zone.

    Args:
        path_elcosmos (string): path of the EL-COSMOS data
        zone (int): EL-COSMOS zone

    Returns:
        int: all identifiers in an input zone
    """
    path = os.path.join(path_elcosmos, f"/zone_{zone:d}")

    specs = os.listdir(path)
    specs = [x for x in specs if '.fits' in x]
    specs.sort()

    return specs

#Synthetic photometry
def synth_phot(wave, flux, wave_lvf, resp_lvf, tol=1e-3, return_photonrate = False):
    """
    Quick synthetic photometry routine.

    Parameters
    ----------
    wave : `numpy.ndarray`
        wavelength of input spectrum.
    flux : `numpy.ndarray`
        flux density of input spectrum in f_nu unit
        if `return_countrate` = True, erg/s/cm2/Hz is assumed
    wave_lvf : `numpy.ndarray`
        wavelength of the response function
    resp_lvf : `numpy.ndarray`
        response function. assume that this is a QE.
    tol : float, optional
        Consider only wavelength range above this tolerence (peak * tol).
        The default is 1e-3.

    Returns
    -------
    synthethic flux density in the input unit
        if return_photonrate = True, photon rates [ph/s/cm2]

    """
    index_filt, = np.where(resp_lvf > resp_lvf.max()*tol)

    index_flux, = np.where(np.logical_and( wave > wave_lvf[index_filt].min(), 
                                           wave < wave_lvf[index_filt].max() ))

    wave_resamp = np.concatenate( (wave[index_flux], wave_lvf[index_filt]) )
    wave_resamp.sort()
    wave_resamp = np.unique(wave_resamp)

    flux_resamp = np.interp(wave_resamp, wave, flux)
    resp_resamp = np.interp(wave_resamp, wave_lvf, resp_lvf)

    if return_photonrate:
        h_planck = 6.626e-27 # erg/Hz
        return trapezoid(resp_resamp / wave_resamp * flux_resamp, wave_resamp) / h_planck

    return trapezoid(resp_resamp / wave_resamp * flux_resamp, wave_resamp) \
         / trapezoid(resp_resamp / wave_resamp, wave_resamp)

from . import info_7ds as const_7ds
f_nu_sky = const_7ds.f_nu_sky
wl_um_sky = const_7ds.wl_um_sky
Deff = const_7ds.D_EFF
theta_pixel = const_7ds.THETA_PIXEL
I_dark = const_7ds.I_DARK
Npix_ptsrc = const_7ds.NPIX_PTSRC
dQ_RN = const_7ds.DQ_RN

def pointsrc_current_sds(mag_src, wave_sys, resp_sys):
    """
    Calculate SN for a point source

    Input
        mag_src: AB mag of the source, scalar
        Tsamp: individual exposure time [sec], can be scalar or array

    WARNING: !!! ALL VARIABLES ARE GLOBALLY DECLARED !!!
    """

    f_nu_src = f_nu_sky*0 + 10**(-0.4*(mag_src + 48.6))  # erg/s/cm2/Hz
#     f_nu_src = 10**(-0.4*(mag_src + 48.6))  # erg/s/cm2/Hz

    f_nu_sky_erg = f_nu_sky*(1e-23*1e-6)                     # erg/s/cm2/Hz/arcsec2

    photon_rate_src = synth_phot(wl_um_sky, f_nu_src, wave_sys, resp_sys, return_photonrate=True)  # ph/s/cm2
    photon_rate_sky = synth_phot(wl_um_sky, f_nu_sky_erg, wave_sys, resp_sys, return_photonrate=True)  # ph/s/cm2/arcsec2

    I_photo_src = photon_rate_src * (np.pi/4*Deff**2)                     # [e/s] per aperture (no aperture loss)
    I_photo_sky = photon_rate_sky * (np.pi/4*Deff**2) * (theta_pixel**2)  # [e/s] per pixel 

    return I_photo_src, I_photo_sky, I_dark

def pointsrc_sn_sds(I_photo_src, I_photo_sky, I_dark, Tsamp):

    Naper = Npix_ptsrc 

    Q_photo_src = I_photo_src * Tsamp
    Q_photo_sky = I_photo_sky * Tsamp
    Q_photo_dark = I_dark * Tsamp

    SN = Q_photo_src / np.sqrt(Q_photo_src + Naper*Q_photo_sky + Naper*Q_photo_dark + Naper*dQ_RN**2)

    return SN
