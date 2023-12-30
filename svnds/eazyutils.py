import numpy as np
from astropy.table import Table
from astropy.io import ascii

import os
from svnds import PATH_DATA_ELCOSMOS

def find_ylims(xmin, xmax, xs, ys):
    
    y_mins = []
    y_maxs = []
    
    for ix, x in enumerate(xs):
        
        mask = (x >= xmin) & (x <= xmax)
        y_new = ys[ix][mask]
        
        y_mins.append(min(y_new))
        y_maxs.append(max(y_new))
    
    ymin = min(y_mins)
    ymax = max(y_maxs)
    return ymin / 3, ymax * 3


def load_tempfilt(root):
    #tempfilt
    with open(root + '.tempfilt','rb') as f:
        s = np.fromfile(file=f,dtype=np.int32, count=4)
        NFILT=s[0]
        NTEMP=s[1]
        NZ=s[2]
        NOBJ=s[3]
        tempfilt = np.fromfile(file=f,dtype=np.double,count=NFILT*NTEMP*NZ).reshape((NZ,NTEMP,NFILT)).transpose()
        lc = np.fromfile(file=f,dtype=np.double,count=NFILT)
        zgrid = np.fromfile(file=f,dtype=np.double,count=NZ)
        fnu = np.fromfile(file=f,dtype=np.double,count=NFILT*NOBJ).reshape((NOBJ,NFILT)).transpose()
        efnu = np.fromfile(file=f,dtype=np.double,count=NFILT*NOBJ).reshape((NOBJ,NFILT)).transpose()

    tempfilts  = {'NFILT':NFILT,'NTEMP':NTEMP,'NZ':NZ,'NOBJ':NOBJ,\
                 'tempfilt':tempfilt,'lc':lc,'zgrid':zgrid,'fnu':fnu,'efnu':efnu}
    
    return tempfilts

def load_coeff(root):
    # coeff file
    with open(root+'.coeff','rb') as f:
        s = np.fromfile(file=f,dtype=np.int32, count=4)
        NFILT=s[0]
        NTEMP=s[1]
        NZ=s[2]
        NOBJ=s[3]
        coeffs = np.fromfile(file=f,dtype=np.double,count=NTEMP*NOBJ).reshape((NOBJ,NTEMP)).transpose()
        izbest = np.fromfile(file=f,dtype=np.int32,count=NOBJ)
        tnorm = np.fromfile(file=f,dtype=np.double,count=NTEMP)

    coeffs = {'NFILT':NFILT,'NTEMP':NTEMP,'NZ':NZ,'NOBJ':NOBJ,\
              'coeffs':coeffs,'izbest':izbest,'tnorm':tnorm}
    
    return coeffs

def load_tempsed(root):
    # temp_sed file
    with open(root+'.temp_sed','rb') as f:
        s = np.fromfile(file=f,dtype=np.int32, count=3)
        NTEMP=s[0]
        NTEMPL=s[1]
        NZ=s[2]
        templam = np.fromfile(file=f,dtype=np.double,count=NTEMPL)
        temp_seds = np.fromfile(file=f,dtype=np.double,count=NTEMPL*NTEMP).reshape((NTEMP,NTEMPL)).transpose()
        da = np.fromfile(file=f,dtype=np.double,count=NZ)
        db = np.fromfile(file=f,dtype=np.double,count=NZ)

    temp_seds = {'NTEMP':NTEMP,'NTEMPL':NTEMPL,'NZ':NZ,\
              'templam':templam,'temp_seds':temp_seds,'da':da,'db':db}
    
    return temp_seds


def load_pz(root):
    ## PDF of redshift
    with open(root + '.pz','rb') as f:
        s = np.fromfile(file=f,dtype=np.int32, count=2)

        NZ=s[0]
        NOBJ=s[1]

        chi2fit = np.fromfile(file = f, dtype = np.double, count = NZ * NOBJ).reshape((NOBJ, NZ)).transpose()


    pzs = {'NZ':NZ,'NOBJ':NOBJ,\
              'chi2fit':chi2fit}
    
    return pzs

## calculate the pdf of photometric redshift
def load_zpdf(idx, pzs, zgrid):
    
    chi2_idx = pzs['chi2fit'][:, idx]

    prior = np.ones(pzs['NZ'])
    zpdf_idx = np.exp(-0.5 * (chi2_idx - min(chi2_idx))) * prior
    zpdf_idx /= np.trapz(zpdf_idx, zgrid)
    zpdf_idx /= sum(zpdf_idx)
    
    return zpdf_idx



def load_obs(idx, temp_restframe, z_best, templam, das, dbs):
    
    temp_idx = temp_restframe[:, idx]
    z_idx = z_best[idx]
    
    temp_igm = temp_idx.copy() #* (1 + z_idx)**2
    templam_idx = templam.copy() * (1 + z_idx)
    
    lim1 = np.where(templam_idx < 912)
    temp_igm[lim1] *= 0

    lim2 = np.where((templam_idx  >= 912) & (templam_idx  < 1026))
    db = dbs[idx]
    temp_igm[lim2] *= db

    lim3 = np.where((templam_idx >= 1026) & (templam_idx < 1216))
    da = das[idx]
    temp_igm[lim3] *= da
    
    
    return templam_idx, temp_igm / (1e-23 * 1e-6)


def load_inputs(idx, tempfilts):
    
    phot_lam = tempfilts['lc']
    phot_flux = tempfilts['fnu'][:,idx] / (1e-23 * 1e-6)
    phot_fluxerr = tempfilts['efnu'][:, idx] / (1e-23 * 1e-6)
    
    return phot_lam, phot_flux, phot_fluxerr



elcat_light = Table.read(os.path.join(PATH_DATA_ELCOSMOS, "elcosmos_light.csv"))

# synphot = Table.read("/data8/EL_COSMOS/synphots/SDS/original/synphot_7ds_yrs_all_tot.csv") #synthetic photometry
# elcat = ascii.read('/data8/EL_COSMOS/ELCOSMOS_v1.cat')

# sed_ids = [x for x in synphot['spec']]
# sed_ids_int = [int(x) for x in sed_ids] # synphot[np.argsort(sed_ids_int)]

# args = np.argsort(np.argsort(sed_ids_int))
# elcat = elcat[args]

def load_elcosmos_spectra(path_elcosmos, spec_id):
    
    ## EL-COSMOS
    # spec_id = catz[idx]['id']
    zone = elcat_light[elcat_light['ID'] == spec_id]['zone'][0]
    spec_path = path_elcosmos + f"zone_{zone}/sed_" + str(spec_id) + '.fits'
    
    spec = Table.read(spec_path)

    wl = spec['wavelength'].value # Anstrom
    f_lambda = spec['flux'].value # erg/s/cm2/A
    f_nu = f_lambda * wl * (wl / 2.99792e18) / (1e-23 * 1e-6)  # micro Jansky
    
    return wl, f_nu

def load_results(root_cat):
    
    tempfilts = load_tempfilt(root_cat)
    coeffs = load_coeff(root_cat)
    temp_seds = load_tempsed(root_cat)
    pzs = load_pz(root_cat)
    
    temp_restframe = (np.dot(temp_seds['temp_seds'],
                             coeffs['coeffs']).T * (temp_seds['templam']/5500.)**2
                     ).T
    
    zgrid = tempfilts['zgrid']
    z_best = zgrid[coeffs['izbest']]
    templam = temp_seds['templam']
    
    dbs = 1. - temp_seds['db'][coeffs['izbest']]
    das = 1. - temp_seds['da'][coeffs['izbest']]
    
    return tempfilts, coeffs, temp_seds, pzs, temp_restframe, zgrid, z_best, templam, dbs, das

