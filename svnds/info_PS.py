from astropy.table import Table
from svnds import PATH_DATA
import os

# import PS1 filters
path_filters = os.path.join(PATH_DATA, 'filters/PANSTARRS/')
g_PS = Table.read(path_filters + 'PS1_gp1.fit')
r_PS = Table.read(path_filters + 'PS1_rp1.fit')
i_PS = Table.read(path_filters + 'PS1_ip1.fit')
z_PS = Table.read(path_filters + 'PS1_zp1.fit')
y_PS = Table.read(path_filters + 'PS1_yp1.fit')
# w_PS = Table.read(path_filters + 'PS1_wp1.fit')

# Parameters
FILTERS = [g_PS, r_PS, i_PS, z_PS, y_PS]
NFILTER = len(FILTERS)
NAMES = ['g_PS', 'r_PS', 'i_PS', 'z_PS', 'y_PS']
MAG5 = [23.3, 23.2, 23.1, 22.3, 21.4] # S/R = 5 limiting magnitudes from VKING DR4

COL_WAVE = 'lambda'
FACTOR_WAVE = 1 / 1000 # nm to micron
COL_RESP = 'throughput'