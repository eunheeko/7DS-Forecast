import numpy as np
# from scipy import interpolate 
# from scipy.ndimage import gaussian_filter

from astropy.table import Table
from astropy import units as u


###### Constants
C_UMS = 3e14               # c in um/s
C = 3e8                    # m/s
H = 6.626e-34              # Planck constant   [J/Hz]
K = 1.38e-23               # Boltzman constant [J/K]
rad2arcsec = (180/np.pi*3600) # 206265 arcsec
arcsec2rad = 1/rad2arcsec

###### Telescope
D = 20.             # effetive diameter [cm]
F = 3.0             # F-number
EFL = D*F           # effective focal length [cm]
WFE = 0.25          # wave front error [um]
RMS_POINTING = 1.0  # pointing accuray [arcsec]

###### Detector
ARRAY = 'HgCdTe'   # detector array type
NPIX = 2048        # [pixels], detector format
dQ_CDS = 12.5      # [e]
I_DARK = 0.01      # [e/s], dark current
PIXEL_SIZE = 18.0  # [um], "pitch"
T_SAMP = 1.5        # sampling time of IR detectors [sec]

# Pixel scale, FOV in spatial direction
rad2arcsec = (180/np.pi*3600) # 206265 arcsec
arcsec2rad = 1/rad2arcsec
THETA_PIXEL = rad2arcsec * PIXEL_SIZE*0.0001/EFL  # [arcsec], pixel scale
THETA_X = THETA_PIXEL * NPIX / 3600               # [deg], FOV in spatial direction (x)

SPHEREx_R = np.array([41, 41, 41, 35, 110, 130])                     # Resolving power
SPHEREx_eff_LVF = np.array([0.97, 0.97, 0.88, 0.86, 0.78, 0.72])     # LVF efficiency

NCHAN = 16         # number of channels per band (or steps)

SPHEREx_lambda_min = np.array([0.75, 1.11, 1.64, 2.42, 3.82, 4.42])  # starting wavelength
SPHEREx_lambda_max = SPHEREx_lambda_min * ((2*SPHEREx_R+1)/(2*SPHEREx_R-1))**NCHAN
SPHEREx_lambda_mid = (SPHEREx_lambda_min + SPHEREx_lambda_max)/2


THETA_SPEC = NPIX / NCHAN * THETA_PIXEL / 60  # [arcmin]

FOV_SPEC = (THETA_SPEC/60) * THETA_X          # [deg^2]

# Length of spectrum per channel
NPIX_PER_CHANNEL = (NPIX) / (NCHAN)

# diffraction-limited PSF size = 1.22 (lambda/D)
theta_diffraction = 1.22 * (SPHEREx_lambda_mid*1e-4) / D * rad2arcsec

# Area of a pixel in steradian
PIXEL_SR = (THETA_PIXEL*arcsec2rad)**2

# (Area) x (solid angle) per pixel [m^2 sr]
AOMEGA = np.pi * (D/2/100)**2 * PIXEL_SR

##### Efficiency
eff_mirrors_Au = (0.965)**3  # Gold coating, 3 mirrors
eff_dichroic = 0.98          # splitter
eff_LVF = SPHEREx_eff_LVF
eff_fpa = 0.75               # Detector Quantum Efficiency (QE)

T_scope = 80.    # temperature of the telescope [K]
T_FPA   = 50.    # temperature of the focal plane array (FPA) [K]

eff_opt   = eff_mirrors_Au * eff_dichroic * eff_LVF
eff_total = eff_opt * eff_fpa

##### Mission cycle
T_mission = 2  # year
resolutionElement_per_survey = 1  # Nyquist = 0.5
Sh_Redun  = 2 / resolutionElement_per_survey # per year [2=visit twice]

# Survey inefficiency margin
# How much will we lose the observing time due to unexpected circumstances?
Sh_Inefficiency = 1.2    # 1.0 = perpect, 1.2 = 20% is wasted

# All-sky steps per year = (4pi / FOV_spec)
Area_allsky = 4*np.pi*(180/np.pi)**2  # [deg^2] = 4pi steradian
# all-sky survey를 위해 1년 동안 필요한 pointing/step 개수
Nsteps_per_year = (Area_allsky/FOV_SPEC) * Sh_Redun * Sh_Inefficiency


Tmin_orbit = 98.6                 # [min], time per orbit 
TM_downlink = 60.                 # [sec/orbit] Downlink 시간
SAA_time = 415                    # [sec/orbit] South Atlantic Anomaly 시간

N_orbits_per_year = 365.25*24*60/Tmin_orbit

# Large steps (여러 개의 small step들로 이루어진)
lg_steps_per_orbit = 8                          # [/orbit] 한 바퀴당 all-sky에 쓸 large step 개수
lg_step_time = Tmin_orbit*60/lg_steps_per_orbit # [sec] large step당 시간

lg_SS_time = 90           # [sec]      large slew 당 필요한 시간 (spacecraft 성능)
sm_SS_time =  8           # [sec]      small slew 당 필요한 시간 (spacecraft 성능)

# Fraction of time to be used fro all-sky
# 전체 시간 중 얼마나 전천탐사에 시간을 쓸 것인가?
frac_allsky = 0.8   

# small steps per one large step to cover all-sky
# 전천을 완전히 커버하기 위해, 큰 step 당 필요한 small step 개수
sm_steps = Nsteps_per_year / N_orbits_per_year / lg_steps_per_orbit  # [/lg step]

T_usable_per_orbit = ( Tmin_orbit*60                                   # Total orbit
                     - lg_steps_per_orbit * lg_SS_time                 # 큰 스텝 사이의 이동시간
                     - lg_steps_per_orbit * (sm_steps-1) * sm_SS_time  # 작은 스텝 사이의 이동시간 (all-sky)
                     - TM_downlink                                     # downlink
                     - SAA_time )                                      # SAA

Tint = T_usable_per_orbit / (lg_steps_per_orbit * sm_steps)
Tint *= frac_allsky

######################################## Noise sources
################## Readout noise
dQ_RN_sh = dQ_CDS*np.sqrt(6*T_SAMP/Tint) # [e]

################## Zodiac light => functional form


################## Thermal background from telescope & FPA
# Surface brightness in (nu Inu) [nW/m2/sr]
# Telescope
hc_kTlambda = H * (C_UMS/SPHEREx_lambda_max) / (K*T_scope)
nuInu_scope = (2*H/C**2) * (C_UMS / SPHEREx_lambda_max)**4 / (np.exp(hc_kTlambda) - 1) / 1e-9

# FPA
hc_kTlambda = H * (C_UMS/SPHEREx_lambda_max) / (K*T_FPA)
nuInu_FPA   = (2*H/C**2) * (C_UMS / SPHEREx_lambda_max)**4 / (np.exp(hc_kTlambda) - 1) / 1e-9
# Iphoto

# Count rates [e/s]
I_scope = 1e-9 * nuInu_scope * np.pi*(PIXEL_SIZE*1e-6)**2/SPHEREx_R*eff_LVF*eff_fpa/(H*C_UMS/SPHEREx_lambda_max)
I_FPA   = 1e-9 * nuInu_FPA   * np.pi*(PIXEL_SIZE*1e-6)**2          *eff_fpa/(H*C_UMS/SPHEREx_lambda_max)


######################################################################################################
################# Nchan = 16 ==> iband
I_STEPS  = np.arange(NCHAN, dtype=float)

###########  ibands 16 x Nchanc 6 array
LAMBDA_IS = np.array(
    [list( SPHEREx_lambda_min * ( ( (2*SPHEREx_R+1) / (2*SPHEREx_R-1) )**x ) )
    for x in I_STEPS]
    , dtype = object)

LAMBDA_IS = LAMBDA_IS.astype(float)
def nuInu_ZL(lambda_um, f_ZL=1.7):
    # very rough approximation for ZL
    # nuInu(sky): fit for zodiacal light [nW/m2/sr]
    # f_ZL = a fudge factor for margin
    A_scat = 3800
    T_scat = 5500
    b_scat = 0.4
    A_therm = 5000
    T_therm = 270
    nuInu = f_ZL *(A_scat*(lambda_um**b_scat)*((lambda_um)**(-4))/(np.exp(H*C_UMS/(K*T_scat*lambda_um))-1)
                +A_therm*1000000          *((lambda_um)**(-4))/(np.exp(H*C_UMS/(K*T_therm*lambda_um))-1)                      )
    return nuInu

nuInu_skys = nuInu_ZL(LAMBDA_IS) #[ nW/m2/sr]
##### sky (zodi here) & instrument & dark + readout noses [e/s]
I_skys = 1e-9 * nuInu_skys*AOMEGA*eff_opt*eff_fpa/(SPHEREx_R * H * C_UMS/LAMBDA_IS)
I_noises = I_skys + I_scope + I_FPA

dQ_noises = np.sqrt( 1.2*(I_noises + I_DARK)*Tint )
dI_noises = np.sqrt( dQ_noises**2 + dQ_RN_sh**2 ) / Tint

#extended souece, [nW/m2/sr]
dnuInu_noises = (dI_noises/I_noises)*(nuInu_skys + nuInu_scope + nuInu_FPA)

# point source
FWHM_diffraction = 1.22*LAMBDA_IS/(D*1e4) * rad2arcsec
FWHM_wfe = FWHM_diffraction * np.sqrt(np.exp((2*np.pi*WFE/LAMBDA_IS)**2))
FWHM_jitter = RMS_POINTING * 2.35 * np.ones_like(LAMBDA_IS)
FWHM = np.sqrt(FWHM_wfe**2 + FWHM_jitter**2)

# N(pixels) for a point-source
Npix_ptsrc = np.pi*(FWHM/THETA_PIXEL)**2

# Fluxes noises within apertures - dFnu [uJy]
dFnu_noises = np.sqrt(Npix_ptsrc) * 1e26 * 1e6 * PIXEL_SR * (dnuInu_noises*1e-9) * (LAMBDA_IS/C_UMS)