import svnds.info_7DS as const_7ds
import svnds.info_SPHEREx as const_spx

import numpy as np
import matplotlib.pyplot as plt

from astropy.modeling.models import Gaussian1D


class Gen7DS():
    
    def __init__(self):
        '''
        Initializes a Construction pipeline. 
        INPUT: 
        RETURNS:
            None
        '''
        self.ref_lambdas = const_7ds.ref_lambdas
        self.eff_total_smooth = const_7ds.eff_total_smooth
        self.sky_lambdas = const_7ds.sky_lambdas
        self.trans_smooth = const_7ds.trans_smooth
        
    def tophat_trans(self, x, center=0, fwhm=1, smoothness=0.2):

        from scipy.special import erf, erfc

        t_left  = erfc(+((2*(x-center)/fwhm)-1)/smoothness)/2 
        t_right = erfc(-((2*(x-center)/fwhm)+1)/smoothness)/2

        return (t_left*t_right)


    def generate_7dsfilters(self, lambda_7ds, fwhm):
        
        # Ideal top-hat filter
        filters_original = {}

        for ii, wl_cen in enumerate(lambda_7ds):
            fwhm = fwhm #anstrom
            wave_lvf = np.linspace(0.1, 1.0, 1001)
            resp_lvf = self.tophat_trans(wave_lvf, center=wl_cen/1e4, fwhm=fwhm/1e4)

            wave_cen = f'{int(wl_cen):d}'

            filters_original['wave_' + wave_cen] = wave_lvf
            filters_original['resp_' + wave_cen] = resp_lvf

        # Corrected filter
        filters_corrected = {}
        
        for ii, wl_cen in enumerate(lambda_7ds):
            wave_cen = f'{int(wl_cen):d}'

            wave_lvf = filters_original['wave_' + wave_cen]
            resp_lvf = filters_original['resp_' + wave_cen]

            resp_sys = resp_lvf.copy()

            ### consider efficiencies
            intp_qe = np.interp(wave_lvf, self.ref_lambdas, self.eff_total_smooth)
            intp_trans = np.interp(wave_lvf, self.sky_lambdas, self.trans_smooth)

            #sky & CMOS QE
            resp_sys = resp_sys * intp_qe * intp_trans 

            filters_corrected['wave_' + wave_cen] = wave_lvf
            filters_corrected['resp_' + wave_cen] = resp_sys
        
        return filters_original, filters_corrected
        
    def plot_filters(self, lambda_7ds, filters_original, filters_corrected):
        
        colors = plt.cm.Spectral(np.linspace(1, 0, len(filters_original)//2))

        fig, ax = plt.subplots(figsize = (8, 5))

        for ii, wl_cen in enumerate(lambda_7ds):
            wave_cen = f'{int(wl_cen):d}'

            wave_ori = filters_original['wave_' + wave_cen]
            resp_ori = filters_original['resp_' + wave_cen]

            wave_cor = filters_corrected['wave_' + wave_cen]
            resp_cor = filters_corrected['resp_' + wave_cen]

            ax.plot(wave_ori, resp_ori, color = colors[ii], linewidth = 2, linestyle = '--')
            ax.plot(wave_cor, resp_cor, color = colors[ii], linewidth = 2)

        ax.set_xlim(0.35, 0.95)

        ax.set_xlabel('wavelength [$\mu m$]', fontsize = 15)
        ax.set_ylabel('Response', fontsize = 15)

class GenSPx():
    
    def __init__(self):
        '''
        Initializes a Construction pipeline. 
        INPUT: 
        RETURNS:
            None
        '''
        self.SPHEREx_lambda_min = const_spx.SPHEREx_lambda_min
        self.SPHEREx_lambda_max = const_spx.SPHEREx_lambda_max
        self.SPHEREx_R = const_spx.SPHEREx_R
        self.Nchan = const_spx.NCHAN
        self.eff_total = const_spx.eff_total

    def generate_spxfilters(self):
        
        # Ideal Gaussian filter
        filter_spherex = {}

        for iband, (l1, l2, R) in enumerate(zip(self.SPHEREx_lambda_min, self.SPHEREx_lambda_max, self.SPHEREx_R)):
            
            i_steps  = np.arange(self.Nchan, dtype=float)
            lambda_i = l1 * (((2*R+1)/(2*R-1))**i_steps)

            for ii, wl_cen in enumerate(lambda_i):
                stddev = wl_cen / R / 2.35
                lvf_profile = Gaussian1D(amplitude=1.0, mean=wl_cen, stddev=stddev)
                wave_lvf = np.linspace(0.7, 5.1, 5001)

                filter_spherex[f'wave_band{iband+1:d}_{ii+1}'] = wave_lvf
                filter_spherex[f'resp_band{iband+1:d}_{ii+1}'] = lvf_profile(wave_lvf)

        # Corrected filter
        filter_spherex_corrected = {}

        for iband, (l1, l2, R) in enumerate(zip(self.SPHEREx_lambda_min, self.SPHEREx_lambda_max, self.SPHEREx_R)):
            i_steps  = np.arange(self.Nchan, dtype=float)
            lambda_i = l1 * (((2*R+1)/(2*R-1))**i_steps)

            for ii, wl_cen in enumerate(lambda_i):

                wave_lvf = filter_spherex[f'wave_band{iband+1:d}_{ii+1}']
                resp_lvf = filter_spherex[f'resp_band{iband+1:d}_{ii+1}']
                resp_sys = resp_lvf.copy()

                resp_sys *= self.eff_total[iband]

                filter_spherex_corrected[f'wave_band{iband+1:d}_{ii+1}'] = wave_lvf
                filter_spherex_corrected[f'resp_band{iband+1:d}_{ii+1}'] = resp_sys
         
        return filter_spherex, filter_spherex_corrected
        
    def plot_filters(self, filters_spherex, filters_spherex_corrected):
        
        colors = plt.cm.Spectral(np.linspace(1, 0, len(filters_spherex)//2))

        fig, ax = plt.subplots(figsize = (8, 5))
        
        for iband, (l1, l2, R) in enumerate(zip(self.SPHEREx_lambda_min, self.SPHEREx_lambda_max, self.SPHEREx_R)):
            i_steps  = np.arange(self.Nchan, dtype=float)
            lambda_i = l1 * (((2*R+1)/(2*R-1))**i_steps)


            for ii, wl_cen in enumerate(lambda_i):
                wave_ori = filters_spherex[f'wave_band{iband+1:d}_{ii+1}']
                resp_ori = filters_spherex[f'resp_band{iband+1:d}_{ii+1}']
                
                wave_cor = filters_spherex_corrected[f'wave_band{iband+1:d}_{ii+1}']
                resp_cor = filters_spherex_corrected[f'resp_band{iband+1:d}_{ii+1}']
                
                ax.plot(wave_cor, resp_cor, color = colors[iband * 6 + ii], linewidth = 2, linestyle = '--')
                ax.plot(wave_cor, resp_cor, color = colors[iband * 6 + ii], linewidth = 2)


        ax.set_xlabel('wavelength [$\mu m$]', fontsize = 15)
        ax.set_ylabel('Response', fontsize = 15)
