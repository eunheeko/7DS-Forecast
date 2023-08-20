import eazy
import pickle
import os
from svnds import PATH_DATA
import numpy as np

def create_res(path_save, surveys = ['7DS', 'SPHEREx']):
    
    fsets = ''

    if '7DS' in surveys:
        lambda_7ds = np.arange(4000., 9000., 125)

        path_filter = os.path.join(PATH_DATA, "filters/SDS/filters_corrected")

        with open(path_filter, 'rb') as fr:
            filters_corrected = pickle.load(fr)

        for ii, wl_cen in enumerate(lambda_7ds):
            wave_cen = f'{int(wl_cen):d}'
            
            wave = filters_corrected['wave_' + wave_cen] * 1e4
            resp = filters_corrected['resp_' + wave_cen]
                
            fset = eazy.filters.FilterDefinition(wave = wave, throughput = resp, name = f'F{int(wl_cen):d}W250')
            fsets = fsets + fset.for_filter_file() + '\n'

    if 'SPHEREx' in surveys:

        from svnds.info_SPHEREx import SPHEREx_lambda_min, SPHEREx_lambda_max, SPHEREx_R

        lmin = SPHEREx_lambda_min
        lmax = SPHEREx_lambda_max
        resolving_power = SPHEREx_R

        path_filter = os.path.join(PATH_DATA, "filters/SPHEREx/filters_SPHEREx_corrected")

        with open(path_filter, 'rb') as fr:
            filters_spherex_corrected = pickle.load(fr)

        for iband, (l1, l2, R) in enumerate(zip(lmin, lmax, resolving_power)):
            i_steps  = np.arange(16, dtype=float)
            lambda_i = l1 * (((2*R+1)/(2*R-1))**i_steps)
            
            for ii, wl_cen in enumerate(lambda_i):

                wave = filters_spherex_corrected[f'wave_band{iband+1:d}_{ii+1}'] * 1e4
                resp = filters_spherex_corrected[f'resp_band{iband+1:d}_{ii+1}']

                fset = eazy.filters.FilterDefinition(wave = wave, throughput = resp, name = f'SPx_band{iband+1:d}_{ii+1}')
                fsets = fsets + fset.for_filter_file() + '\n'
    
    if 'Euclid' in surveys:
        from svnds.info_EUCLID import FILTERS, COL_WAVE, COL_RESP, NAMES
        filters = FILTERS
        col_wave = COL_WAVE
        col_resp = COL_RESP
        names_filter = NAMES

        fsets = _add_res(fsets, filters, col_wave, col_resp, names_filter)
    
    if 'LSST' in surveys:
        from svnds.info_LSST import FILTERS, COL_WAVE, COL_RESP, NAMES
        filters = FILTERS
        col_wave = COL_WAVE
        col_resp = COL_RESP
        names_filter = NAMES

        fsets = _add_res(fsets, filters, col_wave, col_resp, names_filter)

    if 'VIKING' in surveys:
        from svnds.info_VIKING import FILTERS, COL_WAVE, COL_RESP, NAMES
        filters = FILTERS
        col_wave = COL_WAVE
        col_resp = COL_RESP
        names_filter = NAMES

        fsets = _add_res(fsets, filters, col_wave, col_resp, names_filter)

    with open(path_save, 'w') as f:
        f.write(fsets)

    return

def _add_res(fsets, filters, col_wave, col_resp, names_filter):
    
    for ii, wl_cen in enumerate(names_filter):
    
        wave = filters[ii][col_wave] * 10 #nm to AA
        resp = filters[ii][col_resp]
            
        fset = eazy.filters.FilterDefinition(wave = wave, throughput = resp, name = 'F_' + wl_cen)
        fsets = fsets + fset.for_filter_file() + '\n'
    
    return fsets


# def create_translate(path_save):
#     return

def create_cat(path_save, surveys = ['7DS', 'SPHEREx']):

    header = _add_header(surveys = surveys)
    values = _add_rows()

    with open(path_save, 'w') as f:
        f.write(header + values)

    return

def _add_header(surveys = ['7DS', 'SPHEREx']):

    header = '#'
    header += '\t' + 'id'

    if '7DS' in surveys:
        lambda_7ds = np.arange(4000., 9000., 125)
        
        for ii, wl_cen in enumerate(lambda_7ds):

            wave_cen = int(wl_cen)
            
            name = f'F_F{wave_cen:d}W250'
            header += '\t' + name
            
            name = f'E_F{int(wl_cen):d}W250'
            header += '\t' + name

    if 'SPHEREx' in surveys:
        from svnds.info_SPHEREx import SPHEREx_lambda_min, SPHEREx_lambda_max, SPHEREx_R

        lmin = SPHEREx_lambda_min
        lmax = SPHEREx_lambda_max
        resolving_power = SPHEREx_R

        for iband, (l1, l2, R) in enumerate(zip(lmin, lmax, resolving_power)):
            i_steps  = np.arange(16, dtype=float)
            lambda_i = l1 * (((2*R+1)/(2*R-1))**i_steps)
            
            for ii, wl_cen in enumerate(lambda_i):
                name = f'F_band{iband+1:d}_{ii+1}'
                header += '\t' + name

                name = f'E_band{iband+1:d}_{ii+1}'
                header += '\t' + name

    if 'Euclid' in surveys:
        from svnds.info_EUCLID import NAMES
        names_filter = NAMES
        for ii, wl_cen in enumerate(names_filter):    
            wave_cen = wl_cen
            
            name = 'F_' + wave_cen 
            header += '\t' + name
            
            name = 'E_' + wave_cen
            header += '\t' + name

    if 'LSST' in surveys:
        from svnds.info_LSST import NAMES
        names_filter = NAMES
        for ii, wl_cen in enumerate(names_filter):    
            wave_cen = wl_cen
            
            name = 'F_' + wave_cen 
            header += '\t' + name
            
            name = 'E_' + wave_cen
            header += '\t' + name

    if 'VIKING' in surveys:
        from svnds.info_VIKING import NAMES
        names_filter = NAMES
        for ii, wl_cen in enumerate(names_filter):    
            wave_cen = wl_cen
            
            name = 'F_' + wave_cen 
            header += '\t' + name
            
            name = 'E_' + wave_cen
            header += '\t' + name


    header += '\n'

    return header

def _add_rows():

    values = ''

    return values