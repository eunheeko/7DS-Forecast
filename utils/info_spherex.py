import numpy as np
from scipy import interpolate 
from scipy.ndimage import gaussian_filter

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


