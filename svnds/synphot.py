import os
import numpy as np
import tqdm
import pickle
from astropy.table import Table, vstack

from scipy.integrate import trapezoid
from .skytrans import f_nu_sky, wl_um_sky

from svnds import PATH_DATA

elcat_light = Table.read(os.path.join(PATH_DATA, "elcosmos/elcosmos_light.csv"))

# Load EL-COSMOS catalog
def load_catalog(path_elcosmos, zone):
    """load_catalog returns the spec (EL-COSMOS identifier) for a given zone.

    Args:
        path_elcosmos (string): path of the EL-COSMOS data
        zone (int): EL-COSMOS zone

    Returns:
        int: all identifiers in an input zone
    """
    path = os.path.join(path_elcosmos, f"zone_{zone:d}/")

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

def merge_tbl(path_save, name_merged_tbl):
    synphots = os.listdir(path_save)
    synphots = [x for x in synphots if f'synphot_' in x and 'zone' in x]
    synphots.sort()

    if len(synphots) > 9:
        raise Exception("Number of tables for mergin can not exceed the nomber of zones (9)")

    synphot_tot = None

    for ii, f in enumerate(synphots):
        
        tbl = Table.read(path_save + f)
        
        if ii == 0:
            synphot_tot = tbl
        else:
            synphot_tot = vstack([synphot_tot, tbl])   

    synphot_tot.write(path_save + name_merged_tbl, overwrite = True)
    print("Final table is created")
    return
    
Tsamps = np.array([180, 180 * 364/14, 180 * 364/14 * 3, 180 * 364/14 * 5, 180 * 364/14 * 7, 180 * 365, 180 * 365 * 5])
Tsamps = (Tsamps * 0.5).astype(int)
class Syn7DS():

    def __init__(self, path_elcosmos, path_save, Tsamps = Tsamps):
        """__init__ _summary_

        _extended_summary_

        Args:
            path_elcosmos (string): path to the el-cosmos files
            path_save (string): path to save output files
        """
        self.survey = '7DS'
        self.path_elcosmos = path_elcosmos
        self.path_save = path_save
        self.lambda_7ds = np.arange(4000., 9000., 125)

        self.Tsamps = Tsamps
        
        from .info_7DS import D_EFF, THETA_PIXEL, I_DARK, NPIX_PTSRC, DQ_RN
        self.Deff = D_EFF
        self.theta_pixel = THETA_PIXEL
        self.I_dark = I_DARK
        self.Npix_ptsrc = NPIX_PTSRC
        self.dQ_RN = DQ_RN


        with open(os.path.join(PATH_DATA, "filters/SDS/filters_corrected"), 'rb') as fr:
            self.filters_corrected = pickle.load(fr)  
        return

    def synphot(self, zones = range(1, 10), magscale = False, mag_ref = 18):
        
        if magscale:
            self.mag_ref = mag_ref
            self.flux_ref = 10**((self.mag_ref + 48.6)/(-2.5))

        for zone in zones:
            self.synphot_zone(zone, magscale)

            print(f'zone {zone:d} completed!')

        return

    def synphot_zone(self, zone, magscale):

        specs = load_catalog(self.path_elcosmos, zone)
        
        spec_total = []
        flux_total = []
        sn_Tsamps_total = []

        for sp in tqdm.tqdm(specs):

            #ID
            spec_total.append(sp[4:-5])
            spec_path = self.path_elcosmos + f"/zone_{zone:d}/" + sp

            spec = Table.read(spec_path)
        
            #############
            if magscale:
                idx_elcat = np.where(elcat_light['ID'] == int(sp[4:-5]))[0][0]
                
                mag_rband = elcat_light['r_HSC'][idx_elcat]
                flux_rband = 10**((mag_rband + 48.6)/(-2.5))
                fl_factor = self.flux_ref / flux_rband
            #############

            
            wl = spec['wavelength'] # Anstrom
            f_lambda = spec['flux'] # erg/s/cm2/A

            # raw data 
            f_nu = f_lambda * wl * (wl / 2.99792e18) / (1e-23 * 1e-6)  # micro Jansky
            wl = wl / 10000      # micron

            #synthetic photometry: 7DS
            flux_survey = np.zeros_like(self.lambda_7ds, dtype = float)
            sn_Tsamps_survey = np.zeros(shape = ( len(self.lambda_7ds) * len(Tsamps), ), dtype = float)

            for ii, wl_cen in enumerate(self.lambda_7ds):

                wave_cen = f'{int(wl_cen):d}'

                wave_lvf = self.filters_corrected['wave_' + wave_cen]
                resp_lvf = self.filters_corrected['resp_' + wave_cen]

                fl = synth_phot(wl, f_nu, wave_lvf, resp_lvf) #micro Jy

                fl_erg = fl * 1e-6 * 1e-23
                
                #############
                if magscale:
                    fl_erg *= fl_factor
                #############
                
                flux_survey[ii] = fl_erg #erg

                I_photo_src, I_photo_sky, I_dark = self.pointsrc_current_sds(-2.5 * np.log10(fl_erg) - 48.6, wave_lvf, resp_lvf)

                for it, tt in enumerate(Tsamps):
                    sn = self.pointsrc_sn_sds(I_photo_src, I_photo_sky, I_dark, tt)
                    sn_Tsamps_survey[ii*len(Tsamps) + it] = sn

            flux_total.append(flux_survey)
            sn_Tsamps_total.append(sn_Tsamps_survey)

        spec_total = np.array(spec_total)
        flux_total = np.array(flux_total)
        sn_Tsamps_total = np.array(sn_Tsamps_total)

        # create table
        # columns
        cols = []
        cols.append(spec_total) #plate

        names = ['spec']

        for ii, wl_cen in enumerate(self.lambda_7ds):

            wave_cen = f'{int(wl_cen):d}'

            cols.append(flux_total[:, ii])
            names.append('flux_' + wave_cen)

            for it, tt in enumerate(Tsamps):

                cols.append(flux_total[:, ii] / sn_Tsamps_total[:, ii*len(Tsamps)+it])
                names.append('flux_' + wave_cen + f'_err_{tt:d}')

        filename = self.path_save + 'synphot_' + self.survey 
        if magscale:
            filename += f'_scale_{self.mag_ref:d}'
        filename += f'_zone_{zone:d}.csv'

        tbl_synphot = Table(cols, names = names)
        tbl_synphot.write(filename, overwrite = True) #'synphot_uband_zone9.csv'

        return

    def pointsrc_current_sds(self, mag_src, wave_sys, resp_sys):
        """
        Calculate SN for a point source

        Input
            mag_src: AB mag of the source, scalar
            Tsamp: individual exposure time [sec], can be scalar or array

        """

        f_nu_src = f_nu_sky*0 + 10**(-0.4*(mag_src + 48.6))  # erg/s/cm2/Hz
    #     f_nu_src = 10**(-0.4*(mag_src + 48.6))  # erg/s/cm2/Hz

        f_nu_sky_erg = f_nu_sky*(1e-23*1e-6)                     # erg/s/cm2/Hz/arcsec2

        photon_rate_src = synth_phot(wl_um_sky, f_nu_src, wave_sys, resp_sys, return_photonrate=True)  # ph/s/cm2
        photon_rate_sky = synth_phot(wl_um_sky, f_nu_sky_erg, wave_sys, resp_sys, return_photonrate=True)  # ph/s/cm2/arcsec2

        I_photo_src = photon_rate_src * (np.pi/4*self.Deff**2)                     # [e/s] per aperture (no aperture loss)
        I_photo_sky = photon_rate_sky * (np.pi/4*self.Deff**2) * (self.theta_pixel**2)  # [e/s] per pixel 

        return I_photo_src, I_photo_sky, self.I_dark
    
    def pointsrc_sn_sds(self, I_photo_src, I_photo_sky, I_dark, Tsamp):

        Naper = self.Npix_ptsrc 

        Q_photo_src = I_photo_src * Tsamp
        Q_photo_sky = I_photo_sky * Tsamp
        Q_photo_dark = I_dark * Tsamp

        SN = Q_photo_src / np.sqrt(Q_photo_src + Naper*Q_photo_sky + Naper*Q_photo_dark + Naper*self.dQ_RN**2)

        return SN

# class SynSPx():

#     def __init__():
#         return
    
class SynSvy():

    def __init__(self, survey, path_elcosmos, path_save):

        self.survey = survey

        if self.survey == 'LSST':
            from .info_LSST import FILTERS, MAG5, savepath, COL_WAVE, COL_RESP, FACTOR_WAVE, NAMES
            self.filters = FILTERS
            self.mag5 = MAG5
            self.savpath = savepath
            self.col_wave = COL_WAVE
            self.col_resp = COL_RESP
            self.factor_wave = FACTOR_WAVE
            self.names = NAMES
        elif self.survey == 'Euclid':
            from .info_EUCLID import FILTERS, MAG5, savepath, COL_WAVE, COL_RESP, FACTOR_WAVE, NAMES
            self.filters = FILTERS
            self.mag5 = MAG5
            self.savpath = savepath
            self.col_wave = COL_WAVE
            self.col_resp = COL_RESP
            self.factor_wave = FACTOR_WAVE
            self.names = NAMES
        elif self.survey == 'VIKING':
            from .info_VIKING import FILTERS, MAG5, savepath, COL_WAVE, COL_RESP, FACTOR_WAVE, NAMES
            self.filters = FILTERS
            self.mag5 = MAG5
            self.savpath = savepath
            self.col_wave = COL_WAVE
            self.col_resp = COL_RESP
            self.factor_wave = FACTOR_WAVE
            self.names_filter = NAMES
        else:
            raise Exception("Options: LSST, Euclid, and VIKING")
        
        self.path_elcomos = path_elcosmos
        self.path_save = path_save

        return
    

    def synphot(self, zones = range(1, 10), magscale = False, mag_ref = 18):
        
        if magscale:
            self.mag_ref = mag_ref
            self.flux_ref = 10**((self.mag_ref + 48.6)/(-2.5))

        for zone in zones:
            self.synphot_zone(zone, magscale)

            print(f'zone {zone:d} completed!')

        return

    def synphot_zone(self, zone, magscale):

        specs = load_catalog(self.path_elcosmos, zone)
        
        spec_total = []
        flux_total = []
        flux_err_total = []

        for sp in tqdm.tqdm(specs):
            #ID
            spec_total.append(sp[4:-5])
            spec_path = self.path_elcosmos + f"/zone_{zone:d}/" + sp

            spec = Table.read(spec_path)
        
            #############
            if magscale:
                idx_elcat = np.where(elcat_light['ID'] == int(sp[4:-5]))[0][0]
                
                mag_rband = elcat_light['r_HSC'][idx_elcat]
                flux_rband = 10**((mag_rband + 48.6)/(-2.5))
                fl_factor = self.flux_ref / flux_rband
            #############

            wl = spec['wavelength'] # Anstrom
            f_lambda = spec['flux'] # erg/s/cm2/A

            # raw data 
            f_nu = f_lambda * wl * (wl / 2.99792e18) / (1e-23 * 1e-6)  # micro Jansky
            wl = wl / 10000      # micron

            #synthetic photometry:
            flux_survey = np.zeros(len(self.names_filter), dtype = float)
            flux_err_survey = np.zeros(len(self.names_filter), dtype = float)

            for ii, wl_cen in enumerate(self.names_filter):
                
                wave_cen = wl_cen

                wave_lvf = self.filters[ii][self.col_wave] * self.factor_wave #nm to micron
                resp_lvf = self.filters[ii][self.col_resp]

                fl = synth_phot(wl, f_nu, wave_lvf, resp_lvf) #micro Jy

                fl_erg = fl * 1e-6 * 1e-23 #erg
                
                #############
                if magscale:
                    fl_erg *= fl_factor
                #############

                flux_survey[ii] = fl_erg #erg

                xx = 10**(0.4 * (-2.5 * np.log10(fl_erg) - 48.6 - self.mag5[ii]))
                gamma = 0.04
                sigma = np.sqrt( gamma * xx**2) #mag
                sigma_flux = np.abs( np.log(10) / 2.5 * sigma * fl_erg ) #erg
                
                flux_err_survey[ii] = sigma_flux #erg

            flux_total.append(flux_survey)
            flux_err_total.append(flux_err_survey)

        spec_total = np.array(spec_total)
        flux_total = np.array(flux_total)
        flux_err_total = np.array(flux_err_total)

        ##save table
        cols = []
        cols.append(spec_total) #plate

        names = ['spec']

        for ii, wl_cen in enumerate(self.names_filter):

            wave_cen = wl_cen

            cols.append(flux_total[:, ii])
            cols.append(flux_err_total[:, ii])
            names.append('flux_' + wave_cen)
            names.append('flux_' + wave_cen + '_err')

        filename = self.path_save + 'synphot_' + self.survey 
        if magscale:
            filename += f'_scale_{self.mag_ref:d}'
        filename += f'_zone_{zone:d}.csv'

        tbl_synphot = Table(cols, names = names)
        tbl_synphot.write(filename, overwrite = True)